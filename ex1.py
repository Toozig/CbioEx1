



## import
from Bio import SeqIO
import numpy as np

##constance
GAP = '-'

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
            scoring_dict[val] = matrix[i][j]
    return scoring_dict




def globalScore(seq_one, seq_two, scoring_dict):
    """
    This function evalute the score of two sequences
    :param seq_one: string of sequence
    :param seq_two: string of sequence
    :param scoring_dict: dictionary with the given scores
    :return: the score of the best alligment
    """
    matrix = init_matrix(seq_one, seq_two)





def init_matrix(seq_one, seq_two):
    """
    creates the matrix for the global score algorithim
    """
    # todo put it in a tuple so we know its the begining
    matrix = np.empty((len(seq_one) + 1, len(seq_two) + 1))
    # initilaize matrix
    get_score = lambda x, y: scoring_dict[ord(x) + ord(y)]
    matrix[0][0] = 0
    for i in range(1, len(matrix[0])):
        gap_score = get_score(GAP, seq_one[i - 1])
        matrix[0][i] = matrix[0][i - 1] + gap_score
    for i in range(1, len(seq_two) + 1):
        gap_score = get_score(seq_two[i - 1], GAP)
        matrix[i][0] = matrix[i - 1][0] + gap_score
    return matrix


if __name__ == '__main__':
    scoring_dict = createScoringDict("C:\\Users\\133\\Documents\\school\\year3\\semester a\\CBIO\\ex1\\score_matrix.tsv")
    readfasta("/cs/usr/toozig/year3/Cbio/cBIOEX1/fastas/HelicoverpaArmigera-cMyc.fasta")

