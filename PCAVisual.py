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
        t1 = np.array([t[0][:self.parts[0]], t[1][:self.parts[0]]])
        t2 = np.array([t[0][self.parts[0]:self.parts[1]], t[1][self.parts[0]:self.parts[1]]])
        # t3 = np.array([t[0][4:6], t[1][4:6]])
        # t4 = np.array([t[0][6:8], t[1][6:8]])
        t5 = np.array([t[0][self.parts[1]:], t[1][self.parts[1]:]])
        plt.figure()
        plt.plot(t1[0], t1[1], 'ro', label='Acute')
        plt.plot(t2[0], t2[1], 'bo', label='Challenge')
        # plt.plot(t3[0], t3[1], 'go', label='Nac')
        # plt.plot(t4[0], t4[1], 'yo', label='PFC')
        plt.plot(t5[0], t5[1], 'mo', label='Chronic')
        plt.legend(loc='best', numpoints=1)
        explained_var = pca.explained_variance_ratio_
        plt.xlabel("PC1(explained_variance = " + str('{0:.3f}'.format(explained_var[0])) + ")")
        plt.ylabel("PC2(explained_variance = " + str('{0:.3f}'.format(explained_var[1])) + ")")
        # print("Explained variance ratio: ", pca.explained_variance_ratio_)
        # plt.savefig(os.path.dirname(path) + '/PCAfig.jpg')
        plt.show()


