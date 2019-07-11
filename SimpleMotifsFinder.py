
"""
This finder using hashing to find motifs of any length with any number of mistakes.
"""

k = 20
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
        prev = []
        for i in range(k):
            word_to_hash = word[:i] + "N" + word[i + 1:]
            if word_to_hash in self.glob_hash:
                if self.glob_hash[word_to_hash] in prev:
                    continue
                prev.append(self.glob_hash[word_to_hash])
                self.represent_counter[self.glob_hash[word_to_hash]] += 1
                rep_word = self.glob_hash[word_to_hash]
            else:
                prev.append(word)
                self.glob_hash[word_to_hash] = word
                self.represent_counter[word] = 1
        return (rep_word, self.represent_counter[rep_word])


    def write_motifs(self, file_name = "motifs-" + str(k) + "-mers.txt"):
        file = open(file_name, "w")
        for motif in self.represent_counter:
            if self.represent_counter[motif] > 1:
                file.write(motif + " " + str(self.represent_counter[motif]) + "\n")
        file.close()



    def hash_sequence(self, line):
        ns = "N" * k
        for i in range(k, len(line) + 1):
            word = line[i - k: i]
            if word == ns:
                continue
            self.find_family(word)



