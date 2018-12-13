import numpy as np
from PCAVisual import PCAVisual
from Anova import Anova
from Gene import Gene
import matplotlib.pyplot as plt


#run parameters
PERCENT_OF_SHIFT = 2.5
PERCENT_OF_READS_HIST = 0.5
RATIO_TRESHOLD = 0.1
#updated during the run:
THRESHOLD = 0
SAMPLES_PARTS = [2, 4, 8, 11]
PEAK_WINDOW = 2000


def readTheFile(path):
    """
    Reading the file from tha path and returning a list of Gene objects
     sorted by gene names
    :param path:
    :return:
    """
    file = open(path, 'r')
    file.readline()
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
        start = columns[2]
        end = columns[3]
        strand = columns[4]
        if abs(float(end) - float(start)) > 3000:    #check with Reut, it's not OK
            line = file.readline()
            continue
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
    # a1 = np.transpose(M)[-4]
    # a2 = np.transpose(M)[0]
    # plt.plot(a1, a2, 'ro')
    # plt.show()
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
    for item in alternatives:
        isShifted = False
        for row in item.getSamples():
            meanAmg = np.mean(row[:SAMPLES_PARTS[0]])
            meanLH = np.mean(row[SAMPLES_PARTS[0]:SAMPLES_PARTS[1]])
            meanNAC = np.mean(row[SAMPLES_PARTS[1]:SAMPLES_PARTS[2]])
            meanPFC = np.mean(row[SAMPLES_PARTS[2]:SAMPLES_PARTS[3]])
            meanSTR = np.mean(row[SAMPLES_PARTS[3]:])
            means = [meanAmg, meanLH, meanNAC, meanPFC, meanSTR]
            maxShift = 0
            for mean1 in means:
                for mean2 in means:
                    if mean1 / mean2 >= PERCENT_OF_SHIFT and (mean1 >= RATIO_TRESHOLD or mean2 >= RATIO_TRESHOLD):
                        if mean1 / mean2 > maxShift:
                            maxShift = mean1 / mean2
                        isShifted = True
        if isShifted:
            shifted.append(item)
            item.setMaxShift(maxShift)
    for item in shifted:
        print(item.getMaxShift())
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
        # if name == "Fntb":
        #     print("Fntb " + start + " " + end)
        strand = columns[3]
        data.append(Gene(name, reads, np.array([start, end]).astype(np.uint64), strand, chrm))
        line = file.readline()
    return list(sorted(data, key=lambda x: x.getName()))


def findAnnotatedShifts(shifted, annotation):
    notAnnotated = []
    for shifted_gene in shifted:
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
            if isAnnotated == False and shifted_gene.getName() not in notAnnotated:
                notAnnotated.append(shifted_gene.getName())
    print(len(notAnnotated))
    for gene in notAnnotated:
        print(gene)
    return notAnnotated
    #
    # for gene_s in shifted:
    #     for gene_a in annotation:
    #         if gene_s.getName() == gene_a.getName():
    #             coordinate = gene_s.getCoordinates()
    #             for i in range(len(coordinate)):
    #                 if gene_s.getStrand() == "+" and gene_a.getCoordinates()[1] + PEAK_WINDOW > coordinate[i][1]\
    #                     and gene_a.getCoordinates()[1] - PEAK_WINDOW < coordinate[i][1]:
    #                     annotated.append(Gene(gene_a.getName(), gene_s.getSamples()[i],
    #                                           coordinate[i], "+", gene_s.getChromosome()))
    #                 elif gene_s.getStrand() == "-" and gene_a.getCoordinates()[0] + PEAK_WINDOW > coordinate[i][0]\
    #                     and gene_a.getCoordinates()[0] - PEAK_WINDOW < coordinate[i][0]:
    #                     annotated.append(Gene(gene_a.getName(), gene_s.getSamples()[i],
    #                                           coordinate[i], "-", gene_s.getChromosome()))
    # temp = True
    # for gene_s in shifted:
    #     for coordinate in gene_s.getCoordinates():
    #         for gene_a in annotated:
    #             if coordinate[0] == gene_a.getCoordinates()[0] and coordinate[1] == gene_a.getCoordinates()[1]:
    #                 temp = False
    #         # if temp:
    #             # print(gene_s.getName(), coordinate)
    #         temp = True





def writeShifted(shifted, path, name):
    counter = 0
    for item in shifted:
        if np.sum(item.getSamples()) > 0:
            counter += 1
    file = open(path[:len(path) - 3] + name, 'w')
    file.write("Number of genes shifted:" + str(counter))
    file.write("\nTreshold:" + str(TRESHOLD))
    file.write("\nShift:" + str(PERCENT_OF_SHIFT * 100) + "%")
    file.write("\n\nGeneName            Amg1    Amg2    LH1     LH2     Nac1    Nac2    Nac3    Nac4	 PFC2	 PFC3	 PFC4	 STR1	 STR2	 STR3	 STR4\n")
    for item in shifted:
        if np.sum(item.getSamples()) == 0:
            continue
        file.write(item.getName())
        file.write(" " * (20 - len(item.getName())))
        row = item.getSamples()[0]
        for read in row:
            file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
            file.write('\t')
        file.write('\n')
        for row in item.getSamples()[1:]:
            file.write(" " * 20)
            for read in row:
                file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
                file.write('\t')
            file.write('\n')
        file.write('\n')
    file.close()




def main():
    path = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Project_AltPolyA\\data\\data_0h_Acute_by_tissues.window.txt"
    fromFile = readTheFile(path)
    alternatives = findAlternatives(fromFile)
    fracs = calculateFractions(alternatives)
    shifts = findShifts(fracs)
    annotations = readAnnotation("C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\data\\mm10_form_ucsc.txt")
    findAnnotatedShifts(shifts, annotations)
    # genes = []
    # for item in alternatives:
    #     for row in item.getSamples():
    #         anova = Anova(row, [2,4,8,11,15])
    #         p = anova.get_p_value()
    #         if p <= 0.004:
    #             genes.append(item.getName())
    # shiftWithP = []
    # for item in shifts:
    #     if item.getName() in genes:
    #         shiftWithP.append(item)
    writeShifted(shifts, path, "output.txt")
    # writeShifted(shifts, path)
    # pval = []
    # for item in shifts:
    #     for row in item[1]:
    #         anova = Anova(row, [2,4,8,11,15])
    #         pval.append(anova.get_p_value())
    # print(min(pval), max(pval))
    #from here PCA and ploting:
    # data = []
    # for item in shifts:
    #     for row in item.getSamples():
    #         data.append(row)
    # pca = PCAVisual(data)
    # pca.show(path)



if __name__ == "__main__":
    main()




