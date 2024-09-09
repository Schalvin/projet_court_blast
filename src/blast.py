import numpy as np
from Bio import SeqIO
import Bio.Align.substitution_matrices

# Provided arguments
WORD_LEN = 4
T_SCORE = 13
MAX_HITS_GAP_LEN = 20

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

# def close_words_gen(word):
#    word


def seq_to_words(fastaf):
    """Decomposition of a sequence into words of previously specified length

    Args:
        fastaf (str): path to fasta file containing the sequence to be searched

    Returns:
        dict: keys are the words obtained from the sequence, and value is the position in the sequence of the first letter of the word
    """

    words = {}
    records = list(SeqIO.parse(fastaf, "fasta"))
    if len(records) > 2:
        print("The fasta file provided contains more than one sequence. The first sequence of the file will be used for the blast search.")
    seq = str(records[0].seq)
    for i in range(0, len(seq)-WORD_LEN):
        word = seq[i:i+WORD_LEN]
        if word not in words:
            words[word] = []
        score = sum([aa_sub_score(aa, aa) for aa in word])
        if score > 13:
            words[word].append(i)
    return words


def double_hit_finder(words, seq):
    """Get all non over-lapping, double hits with maximum distance between first position of hit of MAX_HITS_GAP_LEN. 

    Args:
        words (dict): keys are words of length constant WORD_LEN, values are positions in the query sequence. generated using the seq_to_words() function
        seq (str): sequence to find hits in
    Returns:
        list: list of all double hits found in seq in the form of a tuple with 
            poss1 - position on sequence in first hit 
            posq1 - position on query in first hit
            poss2 - position on sequence in second hit
            posq2 - position on query in second hit
    """
    hits = {}
    double_hits = []
    for i in range(0, len(seq)-WORD_LEN):
        wordseq = seq[i:i+WORD_LEN]
        if wordseq in words:
            diags = [int(j)-i for j in words[wordseq]]
            for diag in diags:
                if diag not in hits:
                    # tuple containing the pos in the db seq followed by the pos in the query seq
                    hits[diag] = (i, diag+i)
                elif i-hits[diag][0] in range(WORD_LEN, MAX_HITS_GAP_LEN):
                    double_hits.append((hits[diag]+(i, diag+i)))
                    hits[diag] = (i, diag+i)
                elif i-hits[diag][0] < WORD_LEN:
                    continue
                else:
                    hits[diag] = (i, diag+i)
    return double_hits


def db_search(dbf, words):
    """Search of words in a sequence database in fasta format

    Args:
        dbf (str): path to fasta file containing the db
        words (dict): keys are words of length constant WORD_LEN, values are positions in the query sequence. generated using the seq_to_words() function

    Returns:
        dict: keys are the sequence IDs, 
    """
    results = {'index': ('poss1', 'posq1', 'poss2', 'posq2')}
    for dbrecord in SeqIO.parse(dbf, "fasta"):
        result = double_hit_finder(words, str(dbrecord))
        if len(result) > 1:
            results[dbrecord.id] = result
    return results


if __name__ == "__main__":
    words = seq_to_words("data/P06007.fasta")
    print(db_search("data/photosystem_viruses_archea.fasta", words))
