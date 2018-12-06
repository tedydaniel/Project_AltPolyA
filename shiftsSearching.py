import os
import numpy as np
GENE_NAME_COLUMN = 0
NUMBER_OF_READS_COLUMN = 5
START_OF_READ = 2
END_OF_READ = 3
RANGE_FOR_READS = 500



def readingTheFile(path):
    """
    Taking only the 'window' files, as they contain more information
    :param path: path to the directory
    :return: list of opened files containing 'window' in their name
    """
    tissues = []
    for filename in os.listdir(path):
        if "window" in filename:
            tissue = filename[: len(filename) - 11]
            tissues.append([open(path + '\\' + filename), tissue])
    return tissues


def readGenesAndReads(file, tissue):
    """
    :param file: the current file we want to read
    :return: list of tuples [geneName, numberOfReads]
    """
    file.readline() #the first line contains only the names of the columns
    line = file.readline()
    genesAndNumOfReads = []
    while line != '':
        columns = line.split()
        genesAndNumOfReads.append([columns[GENE_NAME_COLUMN], tissue, [int(float(columns[NUMBER_OF_READS_COLUMN])),
                                  int(columns[START_OF_READ]), int(columns[END_OF_READ])]])
        line = file.readline()
    return genesAndNumOfReads




def findingAlternatives(genesAndNumOfReads):
    """
    finds the alternative genes in O(n)   :)
    :param genesAndNumOfReads: List of tuples of gene and number
    of reads for some transcript of this gene
    :return: Sorted list of lists of genes which have an alternative polyA transcript
    """
    alternative = []
    sortedGenesAndNumOfReads = list(sorted(genesAndNumOfReads, key= lambda x : x[0]))
    index = 0

    while index <= (len(sortedGenesAndNumOfReads) - 2):
        counter = 1
        while sortedGenesAndNumOfReads[index][0] == sortedGenesAndNumOfReads[index + counter][0]\
            and ((index + counter) <= (len(sortedGenesAndNumOfReads) - 2)):
            sortedGenesAndNumOfReads[index].append(sortedGenesAndNumOfReads[index + counter][2])
            counter += 1
        if len(sortedGenesAndNumOfReads[index]) > 3:
            alternative.append(sortedGenesAndNumOfReads[index])
        index = index + counter
    return alternative



def calculateFractions (genesWithAlternatives):
    """
    :param genesWithAlternatives: list of lists of genes with alternative transcripts
    :return: same list but with relative fraction of reads instead of number of reads
    """
    fractions = []
    for item in genesWithAlternatives:
        temp = []
        temp.append(item[0])
        temp.append(item[1])
        sumOfReads = 0
        for toSum in item[2:]:
            sumOfReads += toSum[0]
        for i in item[2:]:
            temp.append([float(i[0]) / float(sumOfReads), i[1], i[2]])
        fractions.append(temp)
    return fractions

def outputTheResults(path, geneTissueMatrix):
    file = open(path + "\\filtered2.txt", 'w+')
    last = ''
    for item in geneTissueMatrix:
        file.write(item[0])
        file.write(' ' * (20 - len(item[0])))
        for i in range(len(item[1])):
            if i != 0:
                file.write(20 * ' ')
            file.write(item[1][i])
            file.write(' ' * (30 - len(item[1][i])))
            file.write('\t')
            for j in range(item[2].shape[1] - 1):
                file.write(str('{0:.3f}'.format(item[2][i][j])))  # number of digits after the dot
                file.write('\t')
            file.write('\n')
    file.close()


def findCloseReads(final):
    counter = 0
    geneTissueMatrix = []
    while counter < (len(final) - 2):
        temp = [final[counter][2:]]
        tissuesNames = [final[counter][1]]
        index = 1
        while (counter + index < (len(final))) and final[counter + index][0] == final[counter][0]:
            temp.append(final[counter + index][2:])
            tissuesNames.append(final[counter + index][1])
            index += 1
        maximum = 0
        maximumIndex = 0
        for i in range(len(temp) - 1):
            if len(temp[i]) > maximum:
                maximum = len(temp[i])
                maximumIndex = i
        matrix = np.zeros((len(temp), maximum))
        for i in range(maximum):
            matrix[maximumIndex][i] = temp[maximumIndex][i][0]
        for i in range(len(temp) - 1):
            for j in range(len(temp[i]) - 1):
                for k in range(maximum):
                    if abs(temp[i][j][1] - temp[maximumIndex][k][1]) <= RANGE_FOR_READS \
                            and abs(temp[i][j][2] - temp[maximumIndex][k][2]) <= RANGE_FOR_READS:
                        matrix[i][k] = temp[i][j][0]
        geneTissueMatrix.append([final[counter][0], tissuesNames, matrix])
        counter += index
    return geneTissueMatrix




def run(path):
    fractions = []
    sortedFractions = []
    for file in readingTheFile(path):
        print("Proccessing file " + file[1])
        fractions += calculateFractions(findingAlternatives(readGenesAndReads(file[0], file[1])))
        sortedFractions = list(sorted(fractions, key=lambda x: x[0]))
    toLeave = []
    for i in range(len(sortedFractions) - 2):
        if (sortedFractions[i][0] == sortedFractions[i + 1][0]) and (sortedFractions[i][0] not in toLeave):
            toLeave.append(sortedFractions[i][0])
    final = []
    for i in range(len(sortedFractions) - 1):
        if sortedFractions[i][0] in toLeave:
            final.append(sortedFractions[i])
    geneTissueMatrix = findCloseReads(final)

    outputTheResults(path, geneTissueMatrix)
    print("Done calculating fractions of alternative 3' ends")




path = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\data_0h_Acute_txt"
run(path)

"""
Advice for workflow:
1. Read the data for names of genes and number of reads with readGenesAndReads(file)
2. Find and filter only genes with alternative expressions findingAlternatives(genesAndReads)
3. Calculate the fractions of the alternative reads calculateFractions(alternatives)
4. Filter the genes to leave only the common ones.
"""