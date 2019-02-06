import numpy as np
from Gene import Gene
import matplotlib.pyplot as plt


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
    names = [names[0]] + names[5:]
    line = file.readline()
    data = []
    while line != '':
        columns = line.split()
        reads = np.array([float(x) for x in columns[5:]])
        if np.mean(reads) > 1000:
            line = file.readline()
            continue
        name = columns[0]
        chrm = columns[1]
        if chrm == "chrM":
            line = file.readline()
            continue
        start = columns[2]
        end = columns[3]
        strand = columns[4]
        # if abs(float(end) - float(start)) > 3000:    #check with Reut, it's not OK
        #     line = file.readline()
        #     continue
        data.append(Gene(name, reads, np.array([start, end]).astype(np.uint64), strand, chrm))
        line = file.readline()
    return list(sorted(data, key=lambda x: x.getName()))


def oldAndNewScatter(genes):
    x = []
    y = []
    for gene in genes:
        x.append(gene.getSamples()[0])
        y.append(gene.getSamples()[1])
    plt.scatter(x,y)
    plt.xlabel("new sequencing")
    plt.ylabel("old sequencing")
    plt.title("New vs Old sequencing reads")
    plt.show()



def acuteAndChronicScatter(genes):
    x = []
    y = []
    for gene in genes:
        x.append(gene.getSamples()[0])
        y.append(gene.getSamples()[1])
    plt.scatter(x,y, c='red')
    plt.xlabel("acute")
    plt.ylabel("chronic")
    plt.title("Acute vs Chronic 0h experiment")
    plt.show()


path1 = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Project_AltPolyA\\data\\str_old_new.window.txt"
path2 = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Project_AltPolyA\\data\\acute_chronic_str.window.txt"
genes1 = readTheFile(path1)
genes2 = readTheFile(path2)
oldAndNewScatter(genes1)
acuteAndChronicScatter(genes2)



