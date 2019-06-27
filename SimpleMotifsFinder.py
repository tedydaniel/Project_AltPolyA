
"""
This finder using hashing to find motifs of any length with any number of mistakes.
"""

k = 18
delta = 1



class Family:

    def __init__(self):
        self.glob_hash = {}
        self.represent_counter = {}


    def find_family(self, word):
        """
        This function looking for a relative of the word by looking at all the found words
        and if the editing distance is at most 1 then the words are relative.
        The method is to hash all the possible variants of the word when each time one of the k
        letters replaced with 'N' and then hashed.
        :param word: the word to look its relatives
        :return:
        """
        rep_word = word
        for i in range(k):
            word_to_hash = word[:i] + "N" + word[i + 1:]
            if word_to_hash in self.glob_hash:
                self.represent_counter[self.glob_hash[word_to_hash]] += 1
                rep_word = self.glob_hash[word_to_hash]
            else:
                self.glob_hash[word_to_hash] = word
                self.represent_counter[word] = 1
        return (rep_word, self.represent_counter[rep_word])




start = 134000
end = 200000
file = open("reference\\Mus_musculus.GRCm38.dna.chromosome.1.fa", 'r')
file.readline()
line = ""
for i in range(int(start / 60)):
    file.readline()
for i in range(int(end - start / 60) + 1):
    line += file.readline()[:-1]





fm = Family()
while(line != ""):
    for i in range(k, len(line) + 1):
        word = line[i - k: i]
        if word == "N" * k:
            continue
        temp = fm.find_family(word)
        if temp[1] > 1:
            print(word + " " + str(temp))
    line = ""


