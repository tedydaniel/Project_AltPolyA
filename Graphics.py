import numpy as np
from Gene import Gene
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF
import os


class Graphics:


    def showFracsScatered(self, fracs):
        M = fracs[0].getSamples()
        for i in range(1, len(fracs)):
            M = np.append(M, fracs[i].getSamples(), axis=0)
        M = np.transpose(M)
        MSE_M = np.zeros((M.shape[0], M.shape[0]))
        for i in range(M.shape[0]):
            for j in range(M.shape[0]):
                MSE_M[i][j] = np.sum(np.power(M[i] - M[j], 2))
        MSE_M = MSE_M * (255 / np.max(MSE_M))
        plt.imshow(MSE_M, interpolation='nearest', aspect='auto')
        plt.show()


    def dataToHeatMap(self, shifted, names):
        pdf = FPDF()
        counter = 1
        total = len(shifted)
        pdf.add_page()
        pdf.set_text_color(255, 255, 255)
        pdf.set_font('Arial', '', 12)
        shifted.sort(key=lambda x: x.getMaxShift())
        shifted = shifted[::-1]
        for gene in shifted:
            gene.calculate_lengths()
            ax = sns.heatmap(np.transpose(gene.getSamples()), vmin=0.0, vmax=1.0,
                             annot=True, fmt=".3g", linewidths=.5, yticklabels=names, xticklabels=gene.getLengths())
            plt.axes(ax)
            pdf.write(0, gene.getName() + "\n")
            plt.title(gene.getName() + "\n max_shift = " + str('{0:.3f}'.format(gene.getMaxShift())) + " , max_read = "
                      + str(gene.getMaxRead()) + " , mean_read = " + str('{0:.3f}'.format(gene.getMeanRead()))
                      + "\n transcript = " + str(gene.getNumTranscript()) + " , " + gene.getWhatDiffers()
                      + " , p_value = " + str('{0:.5f}'.format(gene.getPValue())))
            plt.yticks(rotation=0)
            plt.savefig("data\\heat_maps\\" + gene.getName() + ".png", dpi=65)

            pdf.image("data\\heat_maps\\" + gene.getName() + ".png")
            # pdf.
            # os.remove("data\\heat_maps\\" + gene.getName() + ".png")
            plt.close()
            print(str(counter) + " out of " + str(total))
            counter += 1
        pdf.output("all_hours_str_above20percent_130shift_annotated_with_lengths.pdf", "F")


    def showMaxShifts(self, shifts, num_to_show=100, show_coordinates=False, show_above=1.0, show_names=False):
        # global TRESHOLD
        max_shifts = []
        for shifted in shifts:
            if shifted.getMaxShift() >= show_above:
                max_shifts.append((shifted.getName(), shifted.getMaxShift(), shifted))
        max_shifts.sort(key=lambda x: x[1])
        max_shifts1 = []
        genes = []
        items = []
        for shifted in max_shifts:
            max_shifts1.append(shifted[1])
            genes.append(shifted[0])
            items.append(shifted[2])
        max_shifts1 = max_shifts1[-num_to_show:]
        genes = genes[-num_to_show:]
        file = open("maxShifts_" + str(THRESHOLD) + ".txt", 'w')
        for item in items:
            file.write(item.getName() + "   " + item.getChromosome() + "    " + str(item.getMaxShift()))
            print(item.getName(), item.getChromosome(), item.getMaxShift())
            if show_coordinates:
                file.write('\n')
                file.writelines(str(item.getCoordinates()))
                file.write('\n')
                print(item.getCoordinates())
        plt.tick_params(labelsize=7)
        plt.plot(max_shifts1)
        if show_names:
            x = range(len(max_shifts1))
            plt.xticks(x, genes, rotation=90)
        plt.show()


    def readTheFile(self, path):
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


    def oldAndNewScatter(self, genes):
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



    def acuteAndChronicScatter(self, genes):
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


# path1 = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Project_AltPolyA\\data\\str_old_new.window.txt"
# path2 = "C:\\Users\\Nudelman\\Desktop\\Files\\Project_CB\\Project_AltPolyA\\data\\acute_chronic_str.window.txt"
# genes1 = Graphics.readTheFile(path1)
# genes2 = Graphics.readTheFile(path2)
# Graphics.oldAndNewScatter(genes1)
# Graphics.acuteAndChronicScatter(genes2)



