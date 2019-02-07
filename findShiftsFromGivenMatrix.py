import numpy as np
from PCAVisual import PCAVisual
from Gene import Gene
import sys
import time
from Graphics import Graphics
import scipy.stats as stats


#run parameters
PERCENT_OF_SHIFT = 1.5
PERCENT_OF_READS_HIST = 0.5
RATIO_TRESHOLD = 0.2
#updated during the run:
THRESHOLD = 20
NOT_ANNOTATED = []
#
SAMPLES_PARTS = [8, 16]
PEAK_WINDOW = 1000


names = []


def readTheFile(path):
    """
    Reading the file from tha path and returning a list of Gene objects
     sorted by gene names
    :param path:
    :return:
    """
    global names
    file = open(path, 'r')
    names = file.readline().split()
    names = names[5:]
    line = file.readline()
    data = []
    while line != '':
        columns = line.split()
        reads = np.array([float(x) for x in columns[5:]])
        name = columns[0]
        chrm = columns[1]
        if chrm == "chrM":
            line = file.readline()
            continue
        start = int(columns[2])
        end = int(columns[3])
        strand = columns[4]
        data.append(Gene(name, reads, np.array([start, end]).astype(np.uint64), strand, chrm))
        line = file.readline()
    return list(sorted(data, key=lambda x: x.getName()))




def findAlternatives(sortedList):
    """
    Gets a sorted list of reads in format [geneName, reads] and returns only the items
    where more than one transcript was found. The reads of the transcripts will be stacked
    into one matrix so the format of the output is an array of item of the next format [geneName, M]
    where M stands for matrix
    :param sortedList:
    :return:
    """
    #zeroing the data below treshold
    global TRESHOLD
    # if THRESHOLD == 0:
    TRESHOLD = readsHistogram(sortedList)
    afterTresholdData = []
    for i in range(len(sortedList)):
        if np.mean(sortedList[i].getSamples()) >= TRESHOLD:
            afterTresholdData.append(sortedList[i])   #leaves only the reads only if the mean of the reads above TRESHOLD
    index = 0
    while index < (len(afterTresholdData) - 1):
        counter = 1
        while afterTresholdData[index].getName() == afterTresholdData[index + counter].getName():
            afterTresholdData[index].appendSamples(afterTresholdData[index + counter].getSamples())
            afterTresholdData[index].appendCoordinates(afterTresholdData[index + counter].getCoordinates())
            counter += 1
        index += counter
    alternatives = []
    for item in afterTresholdData:
        if len(item.getSamples().shape) > 1:
            alternatives.append(item)
    return alternatives


def readsHistogram(sortedList):
    """
    Using histogram and its cumulative histogram to find a treshold where PERCENT_OF_READS_HIST percent
    of the reads greater than the treshold
    :param alternatives:
    :return:
    """
    M = sortedList[0].getSamples()
    for item in sortedList:
        if max(item.getSamples()) < 5000:
            M = np.vstack((M, item.getSamples()))
    m = int(np.max(M))
    M = np.ndarray.flatten(M)
    # M = M[np.nonzero(M)]     #considering zeros or not?
    hist = np.histogram(M, bins= m)
    cumulative = np.cumsum(hist[0])
    mc = np.max(cumulative)
    treshold = 0
    indexesTreshold = np.where(cumulative <= mc * (1 - PERCENT_OF_READS_HIST))
    if indexesTreshold[0].size > 0:
        treshold = np.max(indexesTreshold)
    return treshold




def calculateFractions(alternatives):
    """
    Just a calculation of fractions: calculating the sum of each column and dividing each value in the
    column by its sum.
    :param alternatives: List with alternatively expressed genes
    :return: List with fractions of the reads from the original list
    """
    for item in alternatives:
        item.setMaxRead(np.max(item.getSamples()))
        item.setMeanRead(np.mean(item.getSamples()))
        Mt = np.transpose(item.getSamples())
        for row in Mt:
            s = np.sum(row)
            if s != 0:
                for i in range(len(row)):
                    row[i] /= s
        item.setSamples(np.transpose(Mt))
    return alternatives




def findShifts(alternatives):
    """
    After receiving the genes with alternative transcripts, for each row we will calculate
     the mean of fraction of each tissue. We will calculate the ratio between these means
     and in this way we will look for significant shifts.
    :param alternatives: List of genes and their alternative transcripts
    :return: Filtered list contains only the shifted genes
    """
    shifted = []
    maxShift = 0
    samples = {0: "Acute", 1:"Challenge", 2:"Chronic"}
    for item in alternatives:
        isShifted = False
        maxShift = 0.0
        counter = 0
        what = ""
        p_value = 0
        transcript = 0
        for row in item.getSamples():
            meanAcute = np.mean(row[:SAMPLES_PARTS[0]])
            meanChallenge = np.mean(row[SAMPLES_PARTS[0]:SAMPLES_PARTS[1]])
            meanChronic = np.mean(row[SAMPLES_PARTS[1]:])
            p = stats.kruskal(row[:SAMPLES_PARTS[0]], row[SAMPLES_PARTS[0]:SAMPLES_PARTS[1]],
                                    row[SAMPLES_PARTS[1]:])[1]
            # meanAmg = np.mean(row[:SAMPLES_PARTS[0]])
            # meanLH = np.mean(row[SAMPLES_PARTS[0]:SAMPLES_PARTS[1]])
            # meanNAC = np.mean(row[SAMPLES_PARTS[1]:SAMPLES_PARTS[2]])
            # meanPFC = np.mean(row[SAMPLES_PARTS[2]:SAMPLES_PARTS[3]])
            # meanSTR = np.mean(row[SAMPLES_PARTS[3]:])
            # means = [meanAmg, meanLH, meanNAC, meanPFC, meanSTR]
            means = [meanAcute, meanChallenge, meanChronic]
            for i in range(len(means)):
                mean1 = means[i]
                for j in range(len(means)):
                    mean2 = means[j]
                    if mean2 > 0 and (mean1 / mean2 >= PERCENT_OF_SHIFT) and \
                            (mean1 >= RATIO_TRESHOLD or mean2 >= RATIO_TRESHOLD)\
                            and counter not in item.getNonAnnotated() \
                            and (abs(int(item.getCoordinates()[counter][1]) - int(item.getCoordinates()[counter][0])) <= 500):
                        if (mean1 / mean2) > maxShift:
                            maxShift = (mean1 / mean2)
                            transcript = counter
                            p_value = p
                            what = samples[i] + "-" + samples[j]
                        isShifted = True
            counter += 1
        if isShifted:
            shifted.append(item)
            item.setMaxShift(maxShift)
            item.setNumTranscript(transcript)
            item.setWhatDiffers(what)
            item.setPValue(p_value)
    return shifted


def readAnnotation(path):
    """
    Returns a list of Gene objects sorted by names, the data is annotations
    :param path:
    :return:
    """
    file = open(path, 'r')
    file.readline()
    line = file.readline()
    data = []
    while line != '':
        columns = line.split()
        reads = np.array([])
        name = columns[-4]
        chrm = columns[2]
        start = columns[4]
        end = columns[5]
        strand = columns[3]
        data.append(Gene(name, reads, np.array([start, end]).astype(np.uint64), strand, chrm))
        line = file.readline()
    return list(sorted(data, key=lambda x: x.getName()))


def findAnnotatedShifts(shifted, annotation):
    global NOT_ANNOTATED
    for shifted_gene in shifted:
        counter = 1
        for coordinate in shifted_gene.getCoordinates():
            isAnnotated = False
            for annotated_gene in annotation:
                if shifted_gene.getName() == annotated_gene.getName():
                    if shifted_gene.getStrand() == '+' \
                        and annotated_gene.getCoordinates()[1] + PEAK_WINDOW > coordinate[1] \
                        and annotated_gene.getCoordinates()[1] - PEAK_WINDOW < coordinate[1]:
                        isAnnotated = True
                    elif shifted_gene.getStrand() == '-' \
                        and annotated_gene.getCoordinates()[0] + PEAK_WINDOW > coordinate[0] \
                        and annotated_gene.getCoordinates()[0] - PEAK_WINDOW < coordinate[0]:
                        isAnnotated = True
            if isAnnotated == False and shifted_gene.getName() not in NOT_ANNOTATED:
                NOT_ANNOTATED.append(shifted_gene.getName())
                shifted_gene.addNonAnnotated(counter)
            counter += 1


def writeShifted(shifted, path, name):
    counter = 0
    for item in shifted:
        if np.sum(item.getSamples()) > 0:
            counter += 1
    file = open(path[:len(path) - 3] + name, 'w')
    file.write("Number of genes shifted:" + str(counter))
    file.write("\nTreshold:" + str(TRESHOLD))
    file.write("\nShift:" + str(PERCENT_OF_SHIFT * 100) + "%")
    file.write ("\nNumber of non-annotated genes: " + str(len(NOT_ANNOTATED)) + '\n\n')
    file.write(names[0] + (20 - len(names[0])) * " ")
    for name in names[1:]:
        file.write(name + (12 - len(name)) * ' ')
    file.write('\n')
    for item in shifted:
        if np.sum(item.getSamples()) == 0:
            continue
        if item.getName() in NOT_ANNOTATED:
            file.write(str(item.getNonAnnotated()) + " " + item.getName())
            file.write(" " * (17 - len(item.getNonAnnotated()) - len(item.getName())))
        else:
            file.write(item.getName())
            file.write(" " * (20 - len(item.getName())))
        row = item.getSamples()[0]
        for read in row:
            file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
            file.write(7 * ' ')
        file.write('\n')
        for row in item.getSamples()[1:]:
            file.write(" " * 20)
            for read in row:
                file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
                file.write(7 * ' ')
            file.write('\n')
        file.write('\n')
    file.close()


    # def removeNotAnnotated(shifted):
    #     shifts = []
    #     for gene in shifted:
    #         if gene.getName() in NOT_ANNOTATED:




def main():
    global SAMPLES_PARTS
    SAMPLES_PARTS[0] = int(input("Number of samples of the first experiment: "))
    SAMPLES_PARTS[1] = SAMPLES_PARTS[0] + int(input("Number of samples of the second experiment: "))
    grph = Graphics()
    path = sys.argv[1]
    anotation_path = sys.argv[2]
    output_filename = sys.argv[3]
    print("Reading the file...")
    fromFile = readTheFile(path)

    data = []
    for item in fromFile:
        row = item.getSamples()
        data.append((row - np.mean(row)) / np.std(row))
    pca = PCAVisual(data, SAMPLES_PARTS)
    pca.show(path)


    alternatives = findAlternatives(fromFile)
    if len(alternatives) > 0:
        print("Found alternatives...")
    else:
        print("No alternatives, check the given arguments")
        raise SystemExit
    fracs = calculateFractions(alternatives)
    annotations = readAnnotation(anotation_path)
    print("Checks the annotations...")
    cur = time.time()
    findAnnotatedShifts(fracs, annotations)
    print("Time took to check the annotations: " + str((time.time() - cur)) + " seconds")
    # showFracsScatered(fracs)
    shifts = findShifts(fracs)
    # showMaxShifts(shifts, num_to_show=500, show_above=1.5, show_coordinates=True)
    # showFracsScatered(shifts)
    # findAnnotatedShifts(shifts, annotations)
    # shifts = removeNotAnnotated(shifts)
    print("Writing the output...")
    grph.dataToHeatMap(shifts, names)
    writeShifted(shifts, path, output_filename)
    data = []
    for item in shifts:
        for row in item.getSamples():
            data.append((row - np.mean(row)) / np.std(row))
    pca = PCAVisual(data, SAMPLES_PARTS)
    pca.show(path)



if __name__ == "__main__":
    main()




