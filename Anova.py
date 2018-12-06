"""
This is an ANOVA test for alternative polyA data.
Our factor is tissue. For now we looking at 5 different tissues.
So our levels are the different tissues and our response are the
fractials of the different polyA transcripts.
The data is going to be a vector of all the responses.
"""

import numpy as np


class Anova:


    def __init__(self, data, levels):
        self.data = data
        self.general_mean = np.mean(data)
        levels_temp = []
        self.means = np.zeros(len(levels) - 1)
        self.n_i = np.zeros(len(levels) - 1)
        counter = 0
        for i in range(1, len(levels)):
            levels_temp.append(data[levels[i - 1] : levels[i]])
            self.means[counter] = np.mean(data[levels[i - 1] : levels[i]])
            self.n_i[counter] = levels_temp[counter].shape[0]
            counter += 1
        self.levels = np.array(levels_temp)
        self.means = np.array(self.means)
        self.totSS = np.sum((self.data - self.general_mean)**2)
        self.totDF = data.shape[0] - 1
        self.betSS = float(np.sum(self.n_i * ((self.means - self.general_mean) ** 2)))
        self.betDF = float(self.n_i.shape[0] - 1)
        temp = np.copy(self.levels)
        for i in range(self.levels.shape[0]):
            temp[i] = np.sum(np.power(temp[i] - self.means[i], 2))
        self.withinSS = np.sum(temp)
        self.withinDF = float(np.sum(self.n_i - 1))


    def get_p_value(self):
        f = np.random.f(self.betDF, self.withinDF, 100000)
        betMS = self.betSS / self.betDF
        withinMS = self.withinSS / self.withinDF
        F = 0
        if withinMS !=0:
            F = betMS / withinMS
        p = np.where(f > F)[0].shape[0] / 100000
        return p


# anova = Anova(np.array([60.8, 67.0, 65.0, 68.6, 61.7, 68.7, 67.7, 75.0 , 73.3, 71.8, 69.6,
#                         77.1, 75.2, 71.5, 61.9, 64.2, 63.1, 66.7, 60.3]), [0,5,10,14,19])
# print(anova.get_p_value())







