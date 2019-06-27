"""
This class will used as an object of Gene. Each gene
will have gene name, matrix of start and end coordinates for each transcript,
matrix of vectors of samples, strand and chromosome.
"""

import numpy as np


class Gene:


    def __init__(self, name, m_samples, m_coordinates, strand, chrm, m_cds = []):
        self.name = name
        self.m_samples = m_samples #numpy array
        self.m_coordinates = m_coordinates #numpy array
        self.m_cds = m_cds
        self.strand = strand
        self.chrm = chrm
        self.maxShift = -1.0
        self.nonAnnotatedTranscripts = []
        self.max_read = -1.0
        self.mean_read = -1.0
        self.num_transcript = -1 #number of the transcript that differs
        self.what_differs = "" #what two kinds of experiences differs
        self.p_value = -1
        self.relative_length = []
        self.percent_of_expression = -1
        self.relative_cds_length = []
        self.sequence = ""


    def calculate_lengths(self):
        """
        :return: calculates the difference between the first coordinate of each transcript
        and the first coordinate of the most expressed transcript
        """
        related = 0
        for transcript in range(self.m_samples.shape[0]):
            if np.mean(self.m_samples[transcript]) > np.mean(self.m_samples[related]):
                related = transcript
        for i in range(self.m_coordinates.shape[0]):
            if self.strand == '+':
                self.relative_length.append(self.m_coordinates[i][0] - self.m_coordinates[related][0])
                self.relative_cds_length.append(self.m_cds[related][1] - self.m_coordinates[i][1])
            else:
                self.relative_length.append(self.m_coordinates[related][1] - self.m_coordinates[i][1])
                self.relative_cds_length.append(self.m_coordinates[i][1] - self.m_cds[related][1])


    def getSequence(self):
        """
        Get the sequence of the added 3' UTR
        :return:
        """
        if self.sequence != "":
            return self.sequence

        file = open("reference\\Mus_musculus.GRCm38.dna.chromosome." + str(self.chrm[-1]) + ".fa", 'r')
        file.readline()
        line = ""
        start = np.min(self.m_coordinates) + 500
        end = np.max(self.m_coordinates)
        for i in range(int(start / 60)):
            file.readline()
        for i in range(int((end - start) / 60) + 1):
            line += file.readline()[:-1]
        file.close()
        return line




    def getLengths(self):
        if self.relative_length:
            return self.relative_length
        print("Not calculated, use calculate_lengths method")


    def getCDS(self):
        return self.m_cds


    def setCDS(self, cds):
        self.m_cds = cds


    def getPercentOfExpression(self):
        return self.percent_of_expression


    def setPercentOfExpression(self, new_percent):
        self.percent_of_expression = new_percent


    def addNonAnnotated(self, num):
        self.nonAnnotatedTranscripts.append(num)


    def getNonAnnotated(self):
        return self.nonAnnotatedTranscripts


    def getSamples(self):
        return self.m_samples


    def setSamples(self, samples):
        self.m_samples = samples


    def getCoordinates(self):
        return self.m_coordinates


    def getStrand(self):
        return self.strand


    def getName(self):
        return self.name


    def getChromosome(self):
        return self.chrm


    def appendSamples(self, samples):
        self.m_samples = np.vstack((self.m_samples, samples))


    def appendCoordinates(self, coords):
        self.m_coordinates = np.vstack((self.m_coordinates, coords))


    def setMaxShift(self, maxshift):
        self.maxShift = maxshift


    def getMaxShift(self):
        return self.maxShift

    def setMaxRead(self, maxr):
        self.max_read = maxr

    def getMaxRead(self):
        return self.max_read

    def setMeanRead(self, meanr):
        self.mean_read = meanr

    def getMeanRead(self):
        return self.mean_read

    def setNumTranscript(self, num):
        self.num_transcript = num

    def getNumTranscript(self):
        return self.num_transcript

    def setWhatDiffers(self, what):
        self.what_differs = what

    def getWhatDiffers(self):
        return self.what_differs

    def setPValue(self, p):
        self.p_value = p

    def getPValue(self):
        return self.p_value




