"""
This class will used as an object of Gene. Each gene
will have gene name, matrix of start and end coordinates for each transcript,
matrix of vectors of samples, strand and chromosome.
"""

import numpy as np


class Gene:


    def __init__(self, name, m_samples, m_coordinates, strand, chrm):
        self.name = name
        self.m_samples = m_samples
        self.m_coordinates = m_coordinates
        self.strand = strand
        self.chrm = chrm
        self.maxShift = -1.0
        self.nonAnnotatedTranscripts = []
        self.max_read = -1.0
        self.mean_read = -1.0
        self.num_transcript = -1
        self.what_differs = ""
        self.p_value = 0


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
        return  self.p_value




