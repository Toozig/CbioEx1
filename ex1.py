



## import
#from Bio import SeqIO
import numpy as np

##constance
ORIGIN = 1
SCORE = 0
SCORE_AND_ORIGIN = 2
GAP = '-'
PRINT_LEN = 50
GET_SCORE = lambda x, y: scoring_dict[ord(x) + ord(y)]


class Alignment:
    """
     Holds the information for each cell in the alignment matrix
     """

    def __init__(self, mother=None, index=None, score=0):
        """
        Ctor for the Class
        """

        self.__mother = None
        self.__index = None
        self.__score = 0

    def set_mother(self, mother):
        """
        sets the origin to know the alignment for printing
        """
        self.__mother = mother

    def set_index(self, pair):
        """
        sets the base pairing of this aligment
        """
        self.__index = pair

    def set_score(self, score):
        """
        sets the score for this aligment
        """
        self.__score = score

    def get_mother(self):
        """
        gets the origin to know the alignment for printing
        """
        return self.__mother

    def get_index(self):
        """
        gets the base pairing of this aligment
        """
        return self.__index

    def get_score(self):
        """
        gets the score for this aligment
        """
        return self.__score

    def __repr__(self):
        return str(self.__score) +" " + str(self.__mother)


def readfasta(fastaFile):
    """
    Reads sequence from fasta file
    :param fastaFile: the sequnce file location & name  as a fasta file
    :return: a string of the given sequence
    """
    print(SeqIO.read(fastaFile, "fasta").seq)

def createScoringDict(file):
    matrix = np.genfromtxt(fname=file, delimiter="\t", skip_header=1, filling_values=1)
    matrix = np.delete(matrix,[0], axis=1)
    bases = ['A', 'C', 'G', 'T', '-']
    scoring_dict = {}
    for i in range(len(bases)):
        for j in range(len(bases)):
            val = ord(bases[i]) + ord(bases[j])
            scoring_dict[val] = int(matrix[i][j])
    return scoring_dict




def globalScore(seq_one, seq_two, scoring_dict):
    """
    This function evalute the score of two sequencestype=(int, SCORE_AND_ORIGIN)
    :param seq_one: string of sequence
    :param seq_two: string of sequence
    :param scoring_dict: dictionary with the given scores
    :return: the score of the best alligment
    """
    matrix = init_matrix(seq_one, seq_two, scoring_dict)
    for i in range(1, len(seq_one) + 1):
        for j in range(1, len(seq_two) + 1):
            matrix[i][j].set_index((i, j))
            alignment_tup = get_best_alligment(matrix, (i, j), scoring_dict, seq_one, seq_two)
            matrix[i][j].set_score(alignment_tup[SCORE])
            matrix[i][j].set_mother(alignment_tup[ORIGIN])

    return matrix



def init_matrix(seq_one, seq_two, scoring_dict):
    """
    creates the matrix for the global score algorithim
    """
    # todo put it in a tuple so we know its the begining
    matrix = [[Alignment() for i in range(len(seq_two) + 1)]
              for j in range((len(seq_one) + 1))]


    #initilaize matrix
    matrix[0][0].set_index((0, 0))
    for i in range(1, len(matrix[0])):
        gap_score = GET_SCORE(GAP, seq_two[i - 1])
        matrix[0][i].set_score(matrix[0][i - 1].get_score() + gap_score)
        matrix[0][i].set_index((0, i))


    for i in range(1, len(seq_one) + 1):
        gap_score = GET_SCORE(seq_one[i - 1], GAP)
        matrix[i][0].set_score(matrix[i - 1][0].get_score() + gap_score)
        matrix[i][0].set_index((i, 0))
    return matrix

def get_best_alligment(matrix, cur_idx, scoring_dict, seq_one, seq_two):
    """
    evaluates the value for this index best alligment
    """
    cur_score = GET_SCORE(seq_one[cur_idx[0] - 1], seq_two[cur_idx[1] - 1])
    options_dict = {matrix[cur_idx[0] - 1][cur_idx[1] - 1].get_score() + cur_score:
                        (cur_idx[0] - 1, cur_idx[1] - 1),
                    matrix[cur_idx[0] - 1][cur_idx[1]].get_score() + cur_score:
                        (cur_idx[0] - 1, cur_idx[1]),
                    matrix[cur_idx[0]][cur_idx[1] - 1].get_score() + cur_score:
                        (cur_idx[0], cur_idx[1] - 1)}
    max_score = max(options_dict)

    return max_score, options_dict[max_score]


def print_max_alignment(matrix, seq_one, seq_two):
    seq_one_alig = list()
    seq_two_alig = list()
    cur = matrix[len(seq_one)][len(seq_two)]
    while cur is not None:
        mother = matrix[cur.get_mother()[0]][cur.get_mother()[1]]
        mother_x = mother.get_index()[0]
        mother_y = mother.get_index()[1]
        if cur.get_index()[0] != mother_x and cur.get_index()[1] != mother_y:
            seq_one_alig.append(seq_one[cur.get_index()[0] - 1])
            seq_two_alig.append(seq_two[cur.get_index()[1] - 1])
        elif cur.get_index()[0] == mother_x and cur.get_index()[1] != mother_y:
            seq_one_alig.append(seq_one[cur.get_index()[0] - 1])
            seq_two_alig.append(seq_two[cur.get_index()[1] - 1])
            seq_one_alig.append(GAP)
            seq_two_alig.append(get_next_base(cur, mother, seq_two))
        else:
            seq_one_alig.append(seq_one[cur.get_index()[0] - 1])
            seq_two_alig.append(seq_two[cur.get_index()[1] - 1])
            seq_two_alig.append(GAP)
        cur = mother if cur.get_mother()[0] and cur.get_mother()[1] else None

    seq_one_alig.reverse()
    seq_two_alig.reverse()

    chunks_one = [seq_one_alig[x: x + PRINT_LEN] for x in range(0, len(seq_one_alig), PRINT_LEN)]
    chunks_two = [seq_two_alig[x: x + PRINT_LEN] for x in range(0, len(seq_two_alig), PRINT_LEN)]
    for i in range(len(chunks_one)):
        print(chunks_one[i])
        print(chunks_two[i])
        print ('\n')



if __name__ == '__main__':
    scoring_dict = createScoringDict("/cs/usr/toozig/year3/Cbio/cBIOEX1/score_matrix.tsv")
    # readfasta("/cs/usr/toozig/year3/Cbio/cBIOEX1/fastas/HelicoverpaArmigera-cMyc.fasta")
    a = 'AGCT'
    b = 'GCT'
    matrix = globalScore(a, b, scoring_dict)
    print_max_alignment(matrix, a, b)

