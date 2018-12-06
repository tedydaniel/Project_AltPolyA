import numpy as np
from PCAVisual import PCAVisual
from Anova import Anova
from Gene import Gene


#run parameters
PERCENT_OF_SHIFT = 2
PERCENT_OF_READS_HIST = 0.6
RATIO_TRESHOLD = 0.1
#updated during the run:
THRESHOLD = 0


def readTheFile(path):
    """
    Reading the file from tha path and returning a list sorted by gene names
    of the format [geneName, reads] where reads are an array of reads for each tissue.
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
        # name = columns[0]
        # chr = columns[1]
        # start = columns[2]
        # end = columns[3]
        # strand = columns[4]
        data.append([columns[0], np.array(reads)])
        line = file.readline()
    return list(sorted(data, key=lambda x: x[0]))


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
        if np.mean(sortedList[i][1]) >= TRESHOLD:
            afterTresholdData.append(sortedList[i])   #leaves only the reads only if the mean of the reads above TRESHOLD
    # print(afterTresholdData)
    # afterTresholdData = np.asarray(afterTresholdData)
    index = 0
    while index < (len(afterTresholdData) - 1):
        counter = 1
        while afterTresholdData[index][0] == afterTresholdData[index + counter][0] \
                and np.max(afterTresholdData[index][1]) > 0:
            if np.max(afterTresholdData[index + counter][1]) > 0:
                afterTresholdData[index][1] = np.vstack((afterTresholdData[index][1],
                                                         afterTresholdData[index + counter][1]))
            counter += 1
        index += counter
    alternatives = []
    for item in afterTresholdData:
        if len(item[1].shape) > 1:
            alternatives.append(item)
    return alternatives


def readsHistogram(sortedList):
    """
    Using histogram and its cumulative histogram to find a treshold where PERCENT_OF_READS_HIST percent
    of the reads greater than the treshold
    :param alternatives:
    :return:
    """
    M = sortedList[0][1]
    for item in sortedList:
        M = np.vstack((M, item[1]))
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
        Mt = np.transpose(item[1])
        for row in Mt:
            s = np.sum(row)
            if s != 0:
                for i in range(len(row)):
                    row[i] /= s
        item[1] = np.transpose(Mt)
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
    for item in alternatives:
        isShifted = False
        for row in item[1]:
            meanAmg = np.mean(row[:2])
            meanLH = np.mean(row[2:4])
            meanNAC = np.mean(row[4:8])
            meanPFC = np.mean(row[8:11])
            meanSTR = np.mean(row[11:])
            means = [meanAmg, meanLH, meanNAC, meanPFC, meanSTR]
            for mean1 in means:
                for mean2 in means:
                    if mean1 / mean2 >= PERCENT_OF_SHIFT and (mean1 > RATIO_TRESHOLD and mean2 > RATIO_TRESHOLD): #and?
                        isShifted = True
        if isShifted:
            shifted.append(item)
    return shifted


def writeShifted(shifted, path):
    counter = 0
    for item in shifted:
        if np.sum(item[1]) > 0:
            counter += 1
    file = open(path[:len(path) - 3] + "output.txt", 'w')
    file.write("Number of genes shifted:" + str(counter))
    file.write("\nTreshold:" + str(TRESHOLD))
    file.write("\nShift:" + str(PERCENT_OF_SHIFT * 100) + "%")
    file.write("\n\nGeneName            Amg1    Amg2    LH1     LH2    Nac1     Nac2     Nac3    Nac4	 PFC2	 PFC3	 PFC4	 STR1	 STR2	 STR3	 STR4\n")
    for item in shifted:
        if np.sum(item[1]) == 0:
            continue
        file.write(item[0])
        file.write(" " * (20 - len(item[0])))
        row = item[1][0]
        for read in row:
            file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
            file.write('\t')
        file.write('\n')
        for row in item[1][1:]:
            file.write(" " * 20)
            for read in row:
                file.write(str('{0:.3f}'.format(read)))  # number of digits after the dot
                file.write('\t')
            file.write('\n')
        file.write('\n')
    file.close()




def main():
    path = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Data\\data_0h_Acute_by_tissues.window.txt"
    fromFile = readTheFile(path)
    alternatives = findAlternatives(fromFile)
    shifts = findShifts(calculateFractions(alternatives))
    genes = []
    for item in alternatives:
        for row in item[1]:
            anova = Anova(row, [2,4,8,11,15])
            p = anova.get_p_value()
            if p <= 0.004:
                genes.append(item[0])
    shiftWithP = []
    for item in shifts:
        if item[0] in genes:
            shiftWithP.append(item)
    writeShifted(shiftWithP, path)
    # writeShifted(shifts, path)
    # pval = []
    # for item in shifts:
    #     for row in item[1]:
    #         anova = Anova(row, [2,4,8,11,15])
    #         pval.append(anova.get_p_value())
    # print(min(pval), max(pval))
    #from here PCA and ploting:
    data = []
    for item in shifts:
        for row in item[1]:
            data.append(row)
    pca = PCAVisual(data)
    pca.show(path)



if __name__ == "__main__":
    main()


