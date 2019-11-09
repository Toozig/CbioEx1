## imports
# from Bio import SeqIO
import numpy as np

##constants
INDEX_REDUCTION = 2
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
        return str(self.__score) + " " +"("+ str(self.__index)+')'


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


def globalScore(seq_one, seq_two, scoring_dict, alignment_type):
    """
    This function evalute the score of two sequencestype=(int, SCORE_AND_ORIGIN)
    :param seq_one: string of sequence
    :param seq_two: string of sequence
    :param scoring_dict: dictionary with the given scores
    :return: the score of the best alligment
    """
    if not len(seq_one) and not len(seq_two):
        return None
    max_score = None
    matrix = init_matrix(seq_one, seq_two, scoring_dict, alignment_type)
    for i in range(1, len(seq_two) + 1):
        for j in range(1, len(seq_one) + 1):
            matrix[i][j].set_index((i, j))
            score, mother = get_best_alligment(matrix, (i, j), scoring_dict, seq_one, seq_two)

            if alignment_type == "local":
                mother, score = check_local(i, j, mother, score, seq_one, seq_two)
                if not max_score or max_score.get_score() < score:
                    max_score = matrix[i][j]

            elif alignment_type == "overlap":
                if i == len(seq_two) or j == len(seq_one):
                    max_score = matrix[i][j] if not max_score \
                                                or max_score.get_score() < score else max_score

            matrix[i][j].set_score(score)
            matrix[i][j].set_mother(mother)
    # for i in matrix:
    #     print(i)
    result = matrix[len(seq_two) - 1][len(seq_one) - 1] if alignment_type == 'global' \
        else get_result(matrix, max_score, alignment_type)
    return result


def get_result(matrix, max_score, alignment_type):
    if alignment_type == 'local':
        return max_score
    i, j = max_score.get_index()
    if not i - (len(matrix) - 1):
        add_i, add_j = 0, 1

    else:
        add_i, add_j = 1, 0
    cur = max_score
    score = max_score.get_score()
    while sum(cur.get_index()) != len(matrix[0]) + len(matrix) - INDEX_REDUCTION:
        new_i, new_j = cur.get_index()
        new_i += add_i
        new_j += add_j
        temp = cur
        cur = matrix[new_i][new_j]
        cur.set_mother(temp)
        cur.set_score(score)
    return cur


def check_overlap_score(max_score, score):
    """
    supports the overlap scoring
    """
    if not max_score or max_score.get_score() < score:
        return True
    return False
def check_local(i, j, mother, score, seq_one, seq_two):
    """
    helper for the local alignment scoring
    """
    start_here_score = GET_SCORE(seq_one[j - 1], seq_two[i - 1]) #todo check what todo if negative
    start_here_score = start_here_score if start_here_score > 0 else 0
    if start_here_score > score:
        score = start_here_score
        mother = None
    return mother, score


def printGlobalScore(alignment):
    """
    This function prints the score for the global alignment algorithm.
    :param matrix: the received matrix after performing the algorithm.
    :param seq_one: the first sequence
    :param seq_two: the second sequence
    :return: prints the score
    """
    print(alignment.get_score())


def init_matrix(seq_one, seq_two, scoring_dict, alignment_type):
    """
    inits the first row and first col of the matrix for the global score algorithm
    """
    matrix = [[Alignment() for i in range(len(seq_one) + 1)]
              for j in range((len(seq_two) + 1))]

    # initilaize matrix
    matrix[0][0].set_index((0, 0))
    for i in range(1, len(matrix[0])):
        gap_score = GET_SCORE(GAP, seq_one[i - 1])
        matrix[0][i].set_index((0, i))
        mother = matrix[0][i - 1] if decide_mother(gap_score, alignment_type) else None
        matrix[0][i].set_mother(mother)
        score = decide_score(gap_score, mother, alignment_type)
        matrix[0][i].set_score(score)

    for i in range(1, len(seq_two) + 1):
        gap_score = GET_SCORE(seq_two[i - 1], GAP)
        matrix[i][0].set_index((i, 0))
        mother = matrix[i - 1][0] if decide_mother(gap_score, alignment_type) else None
        matrix[i][0].set_mother(mother)
        score = decide_score(gap_score, mother, alignment_type)
        matrix[i][0].set_score(score)
    return matrix

def decide_score(gap_score, mother, alignment_type):
    """
    decides the score for the initiating objects
    """
    if not mother or alignment_type == 'overlap':
        return 0
    return mother.get_score() + gap_score

def decide_mother(gap_score, alignment_type):
    """
    support function for the local alignment scoring
    """
    if alignment_type == 'global' or (alignment_type == 'local' and gap_score > 0):
        return True
    elif alignment_type == 'local' and gap_score < 0:
        return False
    else:
        return True



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
    origin = options_dict[max_score]
    return max_score, matrix[origin[0]][origin[1]]


def restore_alignment(alignment, seq_one, seq_two):
    """
    This function gets the last alignment in the matrix and restore the
    maximal score alignment
    :param alignment: the last rubic
    :param seq_one: string of seq 1
    :param seq_two: string of set two
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
            ali_seq1 += seq_one[cur_idx[1] - 1]
            ali_seq2 += seq_two[cur_idx[0] - 1]
        elif direction > 0:
            ali_seq1 += GAP
            ali_seq2 += seq_two[cur_idx[0] - 1]
        else:
            ali_seq1 += seq_one[cur_idx[1] - 1]
            ali_seq2 += GAP
        cur = mother

    return ali_seq1, ali_seq2


def print_global_alignment(alignment, seq_one, seq_two):
    ali_seq1, ali_seq2 = restore_alignment(alignment, seq_one, seq_two)
    is_seq1 = True
    counter = PRINT_LEN if len(ali_seq2) > PRINT_LEN else len(ali_seq2)
    saved_counter = counter
    idx = len(ali_seq2)
    while idx:
        to_print = ali_seq1[idx-1] if is_seq1 else ali_seq2[idx-1]
        print(to_print, end="")
        if not counter - 1:
            print()
            is_seq1 = not is_seq1
            if not is_seq1:
                idx += saved_counter - 1
                counter += saved_counter - 1
                continue
            else:
                if idx > PRINT_LEN:
                    saved_counter = PRINT_LEN - 1
                    counter = PRINT_LEN
                else:
                    saved_counter = idx - 1
                    counter = idx
                print()

        counter -= 1
        idx -= 1


#
# if __name__ == '__main__':
#     scoring_dict = createScoringDict("/cs/usr/toozig/year3/Cbio/cBIOEX1/score_matrix.tsv")
#     # readfasta("/cs/usr/toozig/year3/Cbio/cBIOEX1/fastas/HelicoverpaArmigera-cMyc.fasta")
#     a = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCTACGACGCAAACGCACGTCCGCACGGACTCGCTGCCGCGCGT'
#     b = 'CTACGACGCAAACGCACGTCCGCACGGACTCGCTGCCGCGCGTCTACGACGCAAACGCACGTCCGCACGGACTCGCTGCCGCGCGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
#     matrixx = globalScore(a, b, scoring_dict, 'overlap')
#     # print(matrixx)
#     print_global_alignment(matrixx, a,b)
#     # seq1111, seq1112 = restore_alignment(matrix[len(b)][len(a)], a, b)
#     # print_global_alignment(seq1, seq2)
