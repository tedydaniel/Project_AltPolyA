import matplotlib.pyplot as plt




x = [0.1333, 0.2666, 0.5333, 0.6666, 0.9333]

y = [0.28, 0.36, 0.4666, 0.68, 0.68]


plt.plot(x, y)
plt.show()




path1 = "data\\striatum_old.window.output.txt"
path2 = "data\\striatum.window.output.txt"

def readGeneNames(path):
    genes = []
    file = open(path, 'r')
    line = file.readline()
    for i in range(6):
        line = file.readline()
    while line != "":
        line = line.split(" ")
        if not line[0].isdigit():
            genes.append(line[0])
        line = file.readline()
        if line == "\n":
            line = file.readline()
    return genes


genes1 = readGeneNames(path1)
genes2 = readGeneNames(path2)
print(set(genes1).intersection(genes2))