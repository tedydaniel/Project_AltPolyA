import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os


class PCAVisual:


    def __init__(self, data, parts):
        self.data = data
        self.parts = parts


    def show(self, path):
        data = pd.DataFrame(data=np.transpose(self.data))

        x = StandardScaler().fit_transform(data)
        data = pd.DataFrame(data=x)
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(data)
        principalDf = pd.DataFrame(data=principalComponents
                                   , columns=['PC 1', 'PC 2'])
        t = np.transpose(principalDf.values)
        # t1 = np.array([t[0][:self.parts[0]], t[1][:self.parts[0]]])
        # t2 = np.array([t[0][self.parts[0]:self.parts[1]], t[1][self.parts[0]:self.parts[1]]])
        # t3 = np.array([t[0][self.parts[1]:], t[1][self.parts[1]:]])
        t1 = np.array([t[0][:2], t[1][:2]])
        t2 = np.array([t[0][2:4], t[1][2:4]])
        t3 = np.array([t[0][4:8], t[1][4:8]])
        t4 = np.array([t[0][8:11], t[1][8:11]])
        t5 = np.array([t[0][11:], t[1][11:]])
        plt.figure()
        plt.plot(t1[0], t1[1], 'ro', label='AMG')
        plt.plot(t2[0], t2[1], 'bo', label='LH')
        plt.plot(t3[0], t3[1], 'go', label='NAC')
        plt.plot(t4[0], t4[1], 'yo', label='PFC')
        plt.plot(t5[0], t5[1], 'mo', label='STR')
        # plt.plot(t1[0], t1[1], 'ro', label='Acute')
        # plt.plot(t2[0], t2[1], 'bo', label='Chronic')
        # plt.plot(t3[0], t3[1], 'go', label='Challenge')
        plt.legend(loc='best', numpoints=1)
        explained_var = pca.explained_variance_ratio_
        plt.xlabel("PC1(explained_variance = " + str('{0:.3f}'.format(explained_var[0])) + ")", weight='bold')
        plt.ylabel("PC2(explained_variance = " + str('{0:.3f}'.format(explained_var[1])) + ")", weight='bold')
        # print("Explained variance ratio: ", pca.explained_variance_ratio_)
        # plt.savefig(os.path.dirname(path) + '/PCAfig.jpg')
        plt.show()


