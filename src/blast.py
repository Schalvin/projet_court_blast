import numpy as np
from Bio import SeqIO
import Bio.Align.substitution_matrices

# Provided arguments
WORD_LEN = 4
T_SCORE = 13

# Constants
MATRIX = Bio.Align.substitution_matrices.load('BLOSUM62')

def aa_sub_score(aa1, aa2):
    """Get the T score of substitution depending on the provided MATRIX

    Args:
        aa1 (str): first amino acid
        aa2 (str): second amino acid

    Returns:
        float: the substitution score
    """
    i = MATRIX.alphabet.index(aa1)
    j = MATRIX.alphabet.index(aa2)
    tscore = MATRIX[i][j]
    return tscore

#def close_words_gen(word):
#    word

def seq_to_words(fastaf):
    """_summary_

    Args:
        fastaf (str): path to fasta file containing the sequence to be searched

    Returns:
        dict: dictionnary containing all the words and close words obtained from the sequence
    """

    words = {}
    records = list(SeqIO.parse(fastaf, "fasta"))
    if len(records) > 2 :
        print("The fasta file provided contains more than one sequence. The first sequence of the file will be used for the blast search.")
    seq = str(records[0].seq)
    for i in range(0,len(seq)-WORD_LEN):
        word = seq[i:i+WORD_LEN]
        print(word)
        if word not in words:
            words[word] = []
        score = sum([aa_sub_score(aa, aa) for aa in word])
        if score > 13:
            words[word].append(i)
        else: 
            print('nope')
    return words

#def db_word_search(dbf, words):


if __name__ == "__main__" :
    seq_to_words("data/chloroplasts_pdb.fasta")