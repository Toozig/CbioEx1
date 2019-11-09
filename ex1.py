## imports
# from Bio import SeqIO
import numpy as np

##constants
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
        return str(self.__score) + " " +"("+ str(self.__mother)+')'


def readfasta(fastaFile):
    """
    Reads sequence from fasta file
    :param fastaFile: the sequnce file location & name  as a fasta file
    :return: a string of the given sequence
    """
    print(SeqIO.read(fastaFile, "fasta").seq)


def createScoringDict(file):
    """
    This function parses the scoring matrix=S into a dictionary that keeps the sums ins ASCII of two
    letters as keys, and their scores as values.
    :param file: the path of the file for the scoring matrix
    :return: the scoring dictionary
    """
    matrix = np.genfromtxt(fname=file, delimiter="\t", skip_header=1, filling_values=1)
    matrix = np.delete(matrix, [0], axis=1)
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
    for i in range(1, len(seq_two) + 1):
        for j in range(1, len(seq_one) + 1):
            matrix[i][j].set_index((i, j))
            alignment_tup = get_best_alligment(matrix, (i, j), scoring_dict, seq_one, seq_two)
            matrix[i][j].set_score(alignment_tup[SCORE])
            matrix[i][j].set_mother(matrix[alignment_tup[ORIGIN][0]][alignment_tup[ORIGIN][1]])

    return matrix


def printGlobalScore(matrix, seq_one, seq_two):
    """
    This function prints the score for the global alignment algorithm.
    :param matrix: the received matrix after performing the algorithm.
    :param seq_one: the first sequence
    :param seq_two: the second sequence
    :return: prints the score
    """
    print(matrix[len(seq_two) + 1][len(seq_one) + 1])


def init_matrix(seq_one, seq_two, scoring_dict):
    """
    inits the first row and first col of the matrix for the global score algorithm
    """
    matrix = [[Alignment() for i in range(len(seq_one) + 1)]
              for j in range((len(seq_two) + 1))]

    # initilaize matrix
    matrix[0][0].set_index((0, 0))
    for i in range(1, len(matrix[0])):
        gap_score = GET_SCORE(GAP, seq_one[i - 1])
        matrix[0][i].set_score(matrix[0][i - 1].get_score() + gap_score)
        matrix[0][i].set_index((0, i))
        matrix[0][i].set_mother(matrix[0][i - 1])

    for i in range(1, len(seq_two) + 1):
        gap_score = GET_SCORE(seq_two[i - 1], GAP)
        matrix[i][0].set_score(matrix[i - 1][0].get_score() + gap_score)
        matrix[i][0].set_index((i, 0))
        matrix[i][0].set_mother(matrix[i - 1][0])
    return matrix


def get_best_alligment(matrix, cur_idx, scoring_dict, seq_one, seq_two):
    """
    evaluates the value for this index best alligment
    """
    gap_reduction_left = GET_SCORE(seq_two[cur_idx[0] - 1], GAP)
    gap_reduction_up = GET_SCORE(GAP, seq_one[cur_idx[1] - 1])
    pair_score = GET_SCORE(seq_two[cur_idx[0] - 1], seq_one[cur_idx[1] - 1])
    left_score = matrix[cur_idx[0] - 1][cur_idx[1]].get_score() + gap_reduction_left
    up_score = matrix[cur_idx[0]][cur_idx[1] - 1].get_score() + gap_reduction_up
    diagonal_score = matrix[cur_idx[0] - 1][cur_idx[1] - 1].get_score() + pair_score

    options_dict = {diagonal_score: (cur_idx[0] - 1, cur_idx[1] - 1),
                    left_score: (cur_idx[0] - 1, cur_idx[1]),
                    up_score: (cur_idx[0], cur_idx[1] - 1)}
    max_score = max(options_dict)

    return max_score, options_dict[max_score]


def restore_alignment(alignment, seq1, seq2):
    """
    This function gets the last alignment in the matrix and restore the
    maximal score alignment
    :param alignment: the last rubic
    :param seq1: string of seq 1
    :param seq2: string of set two
    :return: tuple with the alignment (reversed!)
    """
    ali_seq1 = str()
    ali_seq2 = str()
    cur = alignment
    while cur.get_mother() is not None:
        mother = cur.get_mother()
        cur_idx = cur.get_index()
        pre_idx = mother.get_index()
        direction = (cur_idx[0] - pre_idx[0]) - (cur_idx[1] - pre_idx[1])
        if not direction:
            ali_seq1 += seq1[cur_idx[1] - 1]
            ali_seq2 += seq2[cur_idx[0] - 1]
        elif direction > 0:
            ali_seq1 += GAP
            ali_seq2 += seq2[cur_idx[0] - 1]
        else:
            ali_seq1 += seq1[cur_idx[1] - 1]
            ali_seq2 += GAP
        cur = mother

    return ali_seq1, ali_seq2


def print_global_alignment(ali_seq1, ali_seq2):
    is_seq1 = True
    counter = PRINT_LEN if len(ali_seq2) > PRINT_LEN else len(ali_seq2)
    saved_counter = counter
    idx = len(seq1)
    while idx:
        to_print = ali_seq1[idx-1] if is_seq1 else ali_seq2[idx-1]
        print(to_print, end="")
        if not counter - 1:
            print()
            if is_seq1:
                idx += saved_counter - 1
            else:
                print()
            is_seq1 = not is_seq1
            continue
        counter -= 1
        idx -= 1










    # is_seq1 = True
    # counter = PRINT_LEN
    # for i in range(len(ali_seq1)*2, 0, -1):
    #     if i == len(seq1):
    #         print()
    #         is_seq1 = not is_seq1
    #     if counter:
    #         to_print = ali_seq1[i-1 - len(ali_seq2)] if is_seq1 else ali_seq2[i-1]
    #         print(to_print, end="")
    #         counter -= 1
    #     else:
    #         print('\n')
    #         if not is_seq1:
    #             print('\n')
    #         is_seq1 = not is_seq1
    #         counter = PRINT_LEN
#


if __name__ == '__main__':
    scoring_dict = createScoringDict("/cs/usr/toozig/year3/Cbio/cBIOEX1/score_matrix.tsv")
    # readfasta("/cs/usr/toozig/year3/Cbio/cBIOEX1/fastas/HelicoverpaArmigera-cMyc.fasta")
    a = 'AGCT'
    b = 'GCT'
    matrix = globalScore(a, b, scoring_dict)
    for i in matrix:
        print(i)
    seq1, seq2 = restore_alignment(matrix[len(b)][len(a)], a, b)
    print_global_alignment(seq1, seq2)
