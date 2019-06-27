import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF


class Graphics:


    def show_fracs_scattered(self, fracs):
        """
        :param fracs:
        :return:
        """
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


    def data_to_heat_map(self, shifted, names,
                         filename = "all_hours_amg_above20percent_all_"
                                    "annotated_and_significant_with_01_and_15_fold.pdf"):
        """
        Creates a heat map of the samples for each gene within the 'shifted' list and outputs it to the pdf file
        with the name 'filename'. The labels for the samples are taken rom the 'names' list.
        :param shifted: list of 'Gene' objects.
        :param names: list with the names of the samples.
        :return:
        """
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
                             annot=True, fmt=".3g", linewidths=.5, yticklabels=names,
                             xticklabels=gene.getLengths(), cmap='coolwarm')
            plt.axes(ax)
            pdf.write(0, gene.getName() + "\n")
            plt.title(gene.getName() + "\n max_shift = " + str('{0:.3f}'.format(gene.getMaxShift())) + " , max_read = "
                      + str(gene.getMaxRead()) + " , mean_read = " + str('{0:.3f}'.format(gene.getMeanRead()))
                      + "\n transcript = " + str(gene.getNumTranscript()) + " , " + gene.getWhatDiffers()
                      + " , p_value = " + str('{0:.5f}'.format(gene.getPValue())) + "\ncds = " +
                      str(gene.relative_cds_length) )
            plt.yticks(rotation=0)
            plt.savefig("data\\heat_maps\\" + gene.getName() + ".png", dpi=65)
            pdf.image("data\\heat_maps\\" + gene.getName() + ".png")
            plt.close()
            print("Writing to PDF " + str(counter) + " out of " + str(total))
            counter += 1
        pdf.output(filename, "F")

    def hist_of_pvalue(self, shifts):
        """
        Shows the histogram and the cumulative histogram of teh given genes p-values.
        :param shifts: list of 'Gene' object.
        :return: shows two plots, the first one os the histogram of the p-values
        and the second one is the cumulative histogram of the same p-values.
        """
        forhist = []
        for gene in shifts:
            if not gene.getNonAnnotated():
                p = gene.getPValue()
                forhist.append(p)
        plt.hist(forhist, bins=150)
        plt.xlabel('p-value')
        plt.ylabel('number of genes')
        plt.title("Histogram of p-value of the annotated genes(n=" + str(len(forhist)) + ")")
        plt.show()
        forhist = []
        for gene in shifts:
            if not gene.getNonAnnotated():
                p = gene.getPValue()
                forhist.append(p)
        plt.hist(forhist, bins=150, cumulative=True)
        plt.xlabel('p-value')
        plt.ylabel('number of genes')
        plt.title("Cumulative histogram of p-value of the annotated genes(n=" + str(len(forhist)) + ")")
        plt.show()


    def length_by_regulation(self, genes):
        upreg = []
        downreg = []
        updown = []
        downup = []
        for gene in genes:
            gene.calculate_lengths()
            if abs(gene.getLengths()[gene.getNumTranscript() - 1]) > 10000:
                continue
            acute = np.mean(gene.getSamples()[gene.getNumTranscript() - 1][:SAMPLES_PARTS[0]])
            challenge = np.mean(gene.getSamples()[gene.getNumTranscript() - 1][SAMPLES_PARTS[0]:SAMPLES_PARTS[1]])
            chronic = np.mean(gene.getSamples()[gene.getNumTranscript() - 1][SAMPLES_PARTS[1]:])
            if acute <= chronic <= challenge:
                upreg.append(gene.getLengths()[gene.getNumTranscript() - 1])
            elif acute >= chronic >= challenge:
                downreg.append(gene.getLengths()[gene.getNumTranscript() - 1])
            elif acute >= chronic <= challenge:
                downup.append(gene.getLengths()[gene.getNumTranscript() - 1])
            else:
                updown.append(gene.getLengths()[gene.getNumTranscript() - 1])

        x_pos = np.arange(4)
        ctes = [np.mean(upreg), np.mean(downreg), np.mean(updown), np.mean(downup)]
        std_downreg = 0
        if downreg:
            std_downreg = np.std(downreg)
        fig, ax = plt.subplots()
        ax.bar(x_pos, ctes, align='center', alpha=0.5, ecolor='black', capsize=10)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(['Upreg', 'Downreg', 'Updown', 'Downup'])
        ax.set_ylabel("Tail relative length")
        plt.title("Variance of the tail lengths\ns.d. " + str(int(np.std(upreg))) +
                  " , " + str(int(std_downreg)) + " , " + str(int(np.std(updown))) +
                  " , " + str(int(np.std(downup))))
        plt.plot([0] * len(upreg), upreg, 'ko')
        plt.plot([1] * len(downreg), downreg, 'ko')
        plt.plot([2] * len(updown), updown, 'ko')
        plt.plot([3] * len(downup), downup, 'ko')
        plt.show()


    def show_change_in_box(self, gene, segments):
        """
        Changed 30.5
        :param gene:
        :param segments:
        :return:
        """
        acute = gene.getSamples()[gene.getNumTranscript() - 1][:segments[0]]
        challenge = gene.getSamples()[gene.getNumTranscript() - 1][segments[0]:segments[1]]
        chronic = gene.getSamples()[gene.getNumTranscript() - 1][segments[1]:]

        plt.boxplot([acute, chronic, challenge])
        # plt.title(gene.getName() + "\n" + "isoform " + str(gene.getNumTranscript()))
        plt.title(gene.getName() + "\n" + "The Short Isoform", weight='bold')
        plt.ylabel("fraction", weight='bold')
        plt.xticks(np.arange(0, 3) + 1, ["Acute", "Chronic", "Challenge"])
        plt.show()

        # sample1 = np.array([0.247, 0.159, 0.218])
        # sample2 = np.array([0.753, 0.841, 0.782])
        # plt.bar(np.arange(0, 6, 2) + 1, sample1, color='r', bottom=sample2, width=0.5, label="isoform1")
        # plt.bar(np.arange(0, 6, 2) + 1, sample2, width=0.5, label="isoform2")
        # plt.ylim([0, 1.1])
        # plt.xlim([-0.01, 7])
        # plt.ylabel("fraction")
        # plt.legend()
        # plt.title("Camk2a")
        # plt.xticks(np.arange(0, 6, 2) + 1.25, ["Acute", "Chronic", "Challenge"])
        # plt.show()

    def scatter_pval_to_fold(self, shifts, shift=1.5, logpval=1, out=True,
                        name1="sig_nac.txt", name2="all_nac.txt"):
        """
        :param shifts: list of the genes
        :param shift: the value of the shift considered significant
        :param logpval: log10 of pvalue considired significant
        :param out: if True then output the significant genes symbols to name1 file
        and all the genes symbols to name2 file.
        :param name1: file for significant symbols
        :param name2: file for all the symbols from 'shifts' list
        :return: list of the significant genes.
        """
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
        if out:
            file1 = open(name1, 'w')
            file2 = open(name2, 'w')
        for i in range(len(forhist)):
            # if maxshift[i] > 2 and forhist[i] > 1:
            #     print(togoterm[i])
            if out:
                file2.write(togoterm[i] + "\n")
            if maxshift[i] < shift or forhist[i] < logpval:
                bluex.append(forhist[i])
                bluey.append(maxshift[i])
            else:
                topdf2.append(topdf[i])
                if out:
                    file1.write(togoterm[i] + '\n')
                redx.append(forhist[i])
                redy.append(maxshift[i])
        if out:
            file1.close()
            file2.close()
        print("Non significant: " + str(len(bluex)) + "\nSignificant: " + str(len(redx)))
        plt.scatter(bluey, bluex, color='b', label="Low fold change and high p-value", s=40)
        plt.scatter(redy, redx, color='r', label="High fold change and low p-value", s=40)
        plt.figure(figsize=(8,3))
        plt.ylabel("-log10(p-value)", weight='bold')
        plt.xlabel("fold change", weight='bold')
        # plt.title("Scatter plot of p-value (Kruskal Wallis) vs Fold change (=1.5)")
        # plt.legend(loc='upper left')
        plt.show()
        return topdf2

    def show_max_shifts(self, shifts, num_to_show=100, show_coordinates=False, show_above=1.0, show_names=False):
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
        file = open("maxShifts.txt", 'w')
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
        
    
    def fold_change_and_pvalue(self, shifted, min_shift = 1.3):
        """
        Shows plot with 2 y-axes for annotated genes only. The left y axis represents
        the values of the maximal shift. The right y-axis represents the value of the p-value
        for the same gene. The genes on the x-axis sorted by their max shift.
        :param shifted: list of the genes
        :return:
        """
        tograph = []
        for gene in shifted:
            if not gene.getNonAnnotated():
                tograph.append(gene)
        tograph.sort(key= lambda x: x.getMaxShift())
        tograph.reverse()
        fig, ax1 = plt.subplots()
        treshold = [x for x in tograph if x.getMaxShift() >= min_shift]
        ax1.plot(np.arange(len(treshold)), [x.getMaxShift() for x in treshold])
        ax1.set_ylabel("fold change", color='b')
        ax1.set_xlabel("genes ordered by fold change(n=" + str(len(treshold)) + ")")
        ax2 = ax1.twinx()
        ax2.plot(np.arange(len(treshold)), [x.getPValue() for x in treshold], 'ro', markersize=3)
        ax2.set_ylabel("p-value", color='r')
        fig.tight_layout()
        plt.title("Relation between the fold change and the p-value\n"
                  " (fold change >= " + str(min_shift) + ")")
        plt.savefig("data\\fold_vs_pvalue.png")
        plt.show()

    def bars_plot(self, gene, segments):
        """
        :param gene: the gene of interest
        :param segments: list with the number of samples for each experience
        :return: shows the bar plot of the mean values with the values in form of black dots
        """
        acute = gene.getSamples()[gene.getNumTranscript() - 1][:segments[0]]
        challenge = gene.getSamples()[gene.getNumTranscript() - 1][segments[0]:segments[1]]
        chronic = gene.getSamples()[gene.getNumTranscript() - 1][segments[1]:]
        x_pos = np.arange(3)
        ctes = [np.mean(acute), np.mean(chronic), np.mean(challenge)]
        fig, ax = plt.subplots()
        ax.bar(x_pos, ctes, align='center', alpha=0.5, ecolor='black', capsize=10)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(['Acute', 'Chronic', 'Challenge'])
        ax.set_ylabel("Fraction (relative ratio) of the transcript")
        plt.title(gene.getName() + " transcript #" + str(gene.getNumTranscript()))
        plt.plot([0] * len(acute), acute, 'ko')
        plt.plot([1] * len(chronic), chronic, 'ko')
        plt.plot([2] * len(challenge), challenge, 'ko')
        plt.show()

    def calculate_distances_matrix(self, data, names):
        """
        Shows the heat map of the distances matrix. The
        :param data:
        :param names:
        :return:
        """
        data = np.transpose(data)
        distance = np.zeros((data.shape[0], data.shape[0]))
        ys = ["Acute", "Challenge", "Chronic"]
        plt.yticks([0, 1, 2], ys)
        for i in range(data.shape[0]):
            for j in range(data.shape[0]):
                distance[i, j] = np.linalg.norm(np.power((data[i] - data[j]), 2))
        distance += (1 - np.max(distance))
        ax = sns.heatmap(np.transpose(distance), vmin=0.0, vmax=1.0, yticklabels=names, xticklabels=names)
        plt.yticks(rotation=0)
        plt.xticks(rotation=90)
        plt.axes(ax)
        plt.show()




