# imports
import argparse
import sys
import numpy as np
from itertools import groupby
from Bio import SeqIO

# constant s
UPPER = 1
LEFT = -1
DIAGONAL = 0
INDEX_REDUCTION = 2
ORIGIN = 1
SCORE = 0
SCORE_AND_ORIGIN = 2
GAP = '-'
PRINT_LEN = 50
GLOBAL = 0
LOCAL = 1
OVERLAP = -1

# lambda function
GET_SCORE = lambda x, y, scoring_dict: scoring_dict[ord(x) + ord(y)]


class Alignment:
    """
     Holds the information for each cell in the alignment matrix
     """

    def __init__(self, mother=None, index=None, score=0, origin =None):
        """
        Ctor for the Class
        """

        self.mother = mother
        self.index = index
        self.score = score
        self.origin = origin

    def set_mother(self, mother):
        """
        sets the origin to know the alignment for printing
        """
        self.mother = mother

    def set_index(self, pair):
        """
        sets the base pairing of this aligment
        """
        self.index = pair

    def set_score(self, score):
        """
        sets the score for this aligment
        """
        self.score = score

    def get_mother(self):
        """
        gets the origin to know the alignment for printing
        """
        return self.mother

    def get_index(self):
        """
        gets the base pairing of this aligment
        """
        return self.index

    def get_score(self):
        """
        gets the score for this aligment
        """
        return self.score

    def __repr__(self):
        idx = '^_^' if not self.mother else self.mother.get_index()
        return str(self.score) + " " +"("+ str(idx)+')'


def readfasta(fastaFile):
    """
    Reads sequence from fasta file
    :param fastaFile: the sequnce file location & name  as a fasta file
    :return: a string of the given sequence
    """
    return (SeqIO.read(fastaFile, "fasta").seq)


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


def calc_score(seq_one, seq_two, scoring_dict, alignment_type):
    """
    This function evaluates the score of two sequences type=(int, SCORE_AND_ORIGIN)
    :param seq_one: string of sequence
    :param seq_two: string of sequence
    :param scoring_dict: dictionary with the given scores
    :return: the score of the best alignment
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
                mother, score = check_local(i, j, mother, score, seq_one, seq_two, scoring_dict)
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
    result = matrix[len(seq_two)][len(seq_one)] if alignment_type == 'global' \
        else get_result(matrix, max_score, alignment_type)
    return result


def get_result(matrix, max_score, alignment_type):
    """
    This function is for the overlap scoring
    :param matrix:
    :param max_score:
    :param alignment_type:
    :return:
    """
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


def check_local(i, j, mother, score, seq_one, seq_two, scoring_dict):
    """
    helper for the local alignment scoring
    """
    start_here_score = GET_SCORE(seq_one[j - 1], seq_two[i - 1], scoring_dict) #todo check what todo if negative
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


#### new function from ido todo

def create_alignment(gap_score, mother_option, alignment_type, x_idx, y_idx, origin):
    """

    """
    mother = None if not LOCAL - alignment_type and gap_score < 0 else mother_option
    score = 0 if not mother or not alignment_type - OVERLAP else mother.get_score() + gap_score
    return Alignment(mother, (x_idx, y_idx), score, origin)

def init_array(seq_two, scoring_dict, alignment_type):
    alignment_list = [Alignment(None, (0, 0), 0)]

    for i in range(len(seq_two)):
        gap_score = GET_SCORE(GAP, seq_two[i], scoring_dict)
        alignment_list.append(create_alignment(gap_score, alignment_list[-1], alignment_type,0, i, LEFT))
    return alignment_list

def new_get_best_alligment(left_score, upper_score, diag_score, seq_one_base, seq_two_base, scoring_dict):
    """
    evaluates the value for this index best alligment
    """
    gap_reduction_left = GET_SCORE(seq_two_base, GAP, scoring_dict)
    gap_reduction_up = GET_SCORE(GAP,seq_one_base, scoring_dict)
    pair_score = GET_SCORE(seq_two_base, seq_one_base, scoring_dict)
    left_score = left_score + gap_reduction_left
    up_score = upper_score + gap_reduction_up
    diagonal_score = diag_score + pair_score

    options_dict = {diagonal_score: DIAGONAL,
                    left_score: LEFT,
                    up_score: UPPER}
    max_score = max(options_dict)
    origin = options_dict[max_score]
    return max_score, origin


def fill_sec_arr(cur_arr, cur_base, scoring_dict, alignment_type, cur_x, seq_two, len_seq_one):
    gap_score = GET_SCORE(GAP, cur_base, scoring_dict)
    alignment_list = [create_alignment(gap_score, cur_arr[0], alignment_type, 0, cur_x, UPPER)]
    max_alignment = None
    for i in range(len(seq_two)):
        max_score , origin = new_get_best_alligment(alignment_list[-1].get_score(), cur_arr[i + 1].get_score(),
                                                    cur_arr[i].get_score(), seq_two[i], cur_base, scoring_dict)
        if not origin:
            mother = cur_arr[i]
        elif not (UPPER - origin):
            mother = cur_arr[i + 1]
        else:
            mother = alignment_list[-1]
        alignment = Alignment(mother, (i, cur_x), max_score, origin)
        alignment_list.append(alignment)
        if alignment_type:
            is_bound = (not i + 1 - len(seq_two) or not cur_x + 1 - len_seq_one)
            if not alignment_type - OVERLAP and is_bound:
                max_alignment = alignment if not max_alignment or max_score > max_alignment.score else max_alignment
            elif not alignment_type - LOCAL:
                max_alignment = alignment if not max_alignment or max_score > max_alignment.score else max_alignment

    return alignment_list, max_alignment


def new_get_best_alignment(seq_one, seq_two, scoring_dict, alignment_type):

    cur_arr = init_array(seq_two, scoring_dict, alignment_type)
    max_alignment = None
    len_seq_one = len(seq_one)
    for i in range(len(seq_one)):
        cur_arr, temp_max_alignment = fill_sec_arr(cur_arr, seq_one[i], scoring_dict,
                                                   alignment_type, i, seq_two, len_seq_one)
        if  alignment_type:
            max_alignment = temp_max_alignment if not max_alignment or max_alignment.score < temp_max_alignment.score \
                else max_alignment
    return cur_arr[-1] if not max_alignment else new_get_result(seq_one, seq_two, max_alignment, alignment_type)



def new_get_result(seq_one, seq_two, max_score, alignment_type):
    """
    This function is for the overlap scoring
    :param matrix:
    :param max_score:
    :param alignment_type:
    :return:
    """
    if not alignment_type - LOCAL:
        return max_score
    i, j = max_score.get_index()
    if i - len(seq_two):
        add_i, add_j = 0, 1
        origin = UPPER
    else:
        add_i, add_j = 1, 0
        origin = LEFT
    cur = max_score
    score = max_score.score
    last_idx_sum = len(seq_one) + len(seq_two)
    while sum(cur.index) < last_idx_sum - INDEX_REDUCTION:
        (new_i, new_j)= cur.index
        new_i += add_i
        new_j += add_j
        temp = cur
        cur = Alignment(temp, (new_i, new_j), score, origin)
    return cur

def new_restore_alignment(alignment, seq_one, seq_two):
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
    while cur.mother is not None:
        mother = cur.mother
        cur_idx = cur.index
        pre_idx = mother.index
        direction = (cur_idx[0] - pre_idx[0]) - (cur_idx[1] - pre_idx[1])
        if not cur.origin:
            ali_seq1 += seq_one[cur_idx[1]]
            ali_seq2 += seq_two[cur_idx[0]]
        elif cur.origin < 0:
            ali_seq1 += GAP
            ali_seq2 += seq_two[cur_idx[0]]
        else:
            ali_seq1 += seq_one[cur_idx[1]]
            ali_seq2 += GAP
        cur = mother

    return ali_seq1, ali_seq2


def new_print_global_alignment(alignment, seq_one, seq_two):
    """
    This function prints the requested alignment, after calculating the score, by tracing back
    starting with the alignment object, until we reach the beginning
    :param alignment: the alignment object
    :param seq_one: first seq
    :param seq_two: second seq
    :return: nothing. this function prints the alignment
    """
    ali_seq1, ali_seq2 = new_restore_alignment(alignment, seq_one, seq_two)
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

################# end of ido func todo

def init_matrix(seq_one, seq_two, scoring_dict, alignment_type):
    """
    inits the first row and first col of the matrix for the global score algorithm
    """
    matrix = [[Alignment() for i in range(len(seq_one) + 1)]
              for j in range((len(seq_two) + 1))]

    # initilaize matrix
    matrix[0][0].set_index((0, 0))
    for i in range(1, len(matrix[0])):
        gap_score = GET_SCORE(GAP, seq_one[i - 1], scoring_dict)
        matrix[0][i].set_index((0, i))
        mother = matrix[0][i - 1] if decide_mother(gap_score, alignment_type) else None
        matrix[0][i].set_mother(mother)
        score = decide_score(gap_score, mother, alignment_type)
        matrix[0][i].set_score(score)

    for i in range(1, len(seq_two) + 1):
        gap_score = GET_SCORE(seq_two[i - 1], GAP, scoring_dict)
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
    gap_reduction_left = GET_SCORE(seq_two[cur_idx[0] - 1], GAP, scoring_dict)
    gap_reduction_up = GET_SCORE(GAP, seq_one[cur_idx[1] - 1], scoring_dict)
    pair_score = GET_SCORE(seq_two[cur_idx[0] - 1], seq_one[cur_idx[1] - 1], scoring_dict)
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
    """
    This function prints the requested alignment, after calculating the score, by tracing back
    starting with the alignment object, until we reach the beginning
    :param alignment: the alignment object
    :param seq_one: first seq
    :param seq_two: second seq
    :return: nothing. this function prints the alignment
    """
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



def main():
    """
    This main function parses the args from the command, and prints the score and alignment for the
    requested sequences, scoring matrix and the alignment type
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()
    scoring_dict = createScoringDict(command_args.score)
    seq_one = readfasta(command_args.seq_a)
    seq_two = readfasta(command_args.seq_b)
    # if command_args.align_type == 'global':
    #     raise NotImplementedError
    # elif command_args.align_type == 'local':
    #     raise NotImplementedError
    # elif command_args.align_type == 'overlap':
    #     raise NotImplementedError
    alignment = calc_score(seq_one, seq_two, scoring_dict, command_args.align_type)
#    print_global_alignment(alignment, seq_one, seq_two)
    printGlobalScore(alignment)
    # print the best alignment and score



if __name__ == '__main__':
    """
    runs the main function
    """
    # x = sys.argv[1]
    # main()
    scoring_dict = createScoringDict("score_matrix.tsv")
    a = 'CCTCGCTGCTGGTGTGC'
    b = 'CGCT'
    new_print_global_alignment(new_get_best_alignment(a, b, scoring_dict,OVERLAP),a,b)
