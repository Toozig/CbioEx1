



## import
from biopython import SeqIO

##constance


def readfasta(fastaFile):
    """
    Reads sequence from fasta file
    :param fastaFile: the sequnce file location & name  as a fasta file
    :return: a string of the given sequence
    """
    sequence = SeqIO.parse(fastaFile, "fasta")
    a = 3
    return sequence





if __name__ == '__main__':
    readfasta("/cs/usr/toozig/year3/Cbio/cBIOEX1/fastas/HelicoverpaArmigera-cMyc.fasta")