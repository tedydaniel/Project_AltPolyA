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
        # file = open("all_hours_str_symbols_annotated.txt", 'w')
        for gene in shifted:
            # if gene.getNumTranscript() in gene.getNonAnnotated():
            #     continue
            # file.write(gene.getName() + "\n")
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
        # file.close()
        pdf.output("all_hours_str_above20percent_all_annotated_and_significant.pdf", "F")


    def histPValue(self, shifts):
        shifts.sort(key=lambda x: x.getPValue())
        forhist = []
        for gene in shifts:
            if not gene.getNonAnnotated():
                p = gene.getPValue()
                forhist.append(-np.log10(p))
        plt.hist(forhist, bins=100, cumulative=True)
        plt.xlabel('-log10(p-value)')
        plt.ylabel('number of genes')
        plt.title("Cumulative histogram of -log10(p-value) of the annotated genes(n=" + str(len(forhist)) + ")")
        plt.show()


    def scatterPvalFold(self, shifts):
        shifts.sort(key=lambda x: x.getPValue())
        forhist = []
        maxshift = []
        togoterm = []
        topdf = []
        for gene in shifts:
            if not gene.getNonAnnotated():
                togoterm.append(gene.getName())
                topdf.append(gene)
                p = gene.getPValue()
                forhist.append(-np.log10(p))
                maxshift.append(gene.getMaxShift())
        bluex = []
        bluey = []
        redx = []
        redy = []
        topdf2 = []
        file1 = open("interesting.txt", 'w')
        file2 = open("all.txt", 'w')
        for i in range(len(forhist)):
            file2.write(togoterm[i] + "\n")
            if maxshift[i] < np.exp(-15 * (forhist[i] - 2)) + 1.3:
                bluex.append(forhist[i])
                bluey.append(maxshift[i])
            else:
                topdf2.append(topdf[i])
                file1.write(togoterm[i] + '\n')
                redx.append(forhist[i])
                redy.append(maxshift[i])
        file1.close()
        file2.close()
        t = np.arange(0.0, 4.5, 0.1)
        s = [np.exp(-15 * (x - 2)) + 1.3 for x in t]
        plt.plot(t, s, linestyle='dashed')
        plt.ylim([0, 3])
        print(len(bluex), len(bluey), len(redx), len(redy))
        plt.scatter(bluex, bluey, color='b', label="Low fold change and high p-value")
        plt.scatter(redx, redy, color='r', label="High fold change and low p-value")
        plt.xlabel("-log(p-value)")
        plt.ylabel("fold change")
        plt.title("Scatter plot of p-value (Kruskal Wallis) vs Fold change")
        plt.legend(loc='upper left')
        plt.show()
        return topdf2


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
        
    
    def fold_change_and_percent(self, shifted):
        tograph = []
        for gene in shifted:
            if not gene.getNonAnnotated():
                tograph.append(gene)
        tograph.sort(key= lambda x: x.getMaxShift())
        tograph.reverse()
        fig, ax1 = plt.subplots()
        treshold = [x for x in tograph if x.getMaxShift() >= 1.3]
        ax1.plot(np.arange(len(treshold)), [x.getMaxShift() for x in treshold])
        ax1.set_ylabel("fold change", color='b')
        ax1.set_xlabel("genes ordered by fold change(n=" + str(len(treshold)) + ")")
        ax2 = ax1.twinx()
        ax2.plot(np.arange(len(treshold)), [x.getPValue() for x in treshold], 'ro', markersize=3)
        ax2.set_ylabel("p-value", color='r')
        ax2.set_ylim([0.0, 0.1])
        fig.tight_layout()
        plt.title("Relation between the fold change and the p-value\n"
                  " (fold change >= 1.3)")
        plt.savefig("data\\fold_vs_percent_expression.png")
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



