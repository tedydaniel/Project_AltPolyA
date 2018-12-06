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


