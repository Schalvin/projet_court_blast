import math
import numpy as np
from Bio import SeqIO
import Bio.Align.substitution_matrices

# Provided arguments
WORD_LEN = 3
T_SCORE = 13
MAX_HITS_GAP_LEN = 40
l = 0.3
K = 0.13
GAP_PENALTY = 11
GAP_EXT_PENALTY = 1

# Constants
MATRIX = Bio.Align.substitution_matrices.load('BLOSUM62')
similar_groups = [
    {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'},  # Hydrophobic
    {'S', 'T', 'N', 'Q'},                      # Polar
    {'D', 'E', 'K', 'R'},                      # Charged
    {'G', 'A'},                                # Small
    {'F', 'W', 'Y'},                           # Aromatic
    {'C', 'P', 'H'}                            # Special
]


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


def seq_to_words(seq):
    """Decomposition of a sequence into words of previously specified length

    Args:
        fastaf (str): path to fasta file containing the sequence to be searched

    Returns:
        dict: keys are the words obtained from the sequence, and value is the position in the sequence of the first letter of the word
    """

    words = {}
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
    double_hits = [seq]
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
        dict: keys are the sequence IDs, values are lists containing the target sequence, followed by double hits found in tuples with 
            poss1 - position on sequence in first hit 
            posq1 - position on query in first hit
            poss2 - position on sequence in second hit
            posq2 - position on query in second hit
    """
    results = {'index': ('poss1', 'posq1', 'poss2', 'posq2')}
    for dbrecord in SeqIO.parse(dbf, "fasta"):
        result = double_hit_finder(words, str(dbrecord.seq))
        results[dbrecord.id] = result
    return results


def normalized_score(score):
    """Normalizes alignment score depending on database

    Args:
        score (float): Non normalized score of alignement

    Returns:
        float: normalized alignment score
    """
    return (l*score - math.log(K))/math.log(2)


def e_value(norm_score, N):
    """Calculation of e_value for an alignment

    Args:
        norm_score (float): normalized alignment score with noramized_score() function
        N (int): result of the multiplication of the lengths of the aligned sequences

    Returns:
        float: e-value
    """
    return N/2**norm_score


def alignment(seqq, seqt):
    """_summary_

    Args:
        seqq (_type_): _description_
        seqt (_type_): _description_

    Returns:
        _type_: _description_
    """

    score_mat = np.zeros((len(seqq) + 1, len(seqt) + 1))
    direction_mat = np.zeros((len(seqq) + 1, len(seqt) + 1), dtype=int)

    STOP, LEFT, UP, DIAG = 0, 1, 2, 3

    max_score = 0
    max_pos = None

    # fill score and direction matrices
    for i in range(1, len(seqq) + 1):
        for j in range(1, len(seqt) + 1):

            # scores for substitution deletion and insertion
            sub = score_mat[i-1, j-1] + \
                aa_sub_score(seqq[i-1], seqt[j-1])
            dele = score_mat[i, j-1] + \
                (GAP_EXT_PENALTY if (
                    direction_mat[i, j-1] == 2) else GAP_PENALTY)
            ins = score_mat[i-1, j] + \
                (GAP_EXT_PENALTY if (
                    direction_mat[i-1, j] == 1) else GAP_PENALTY)

            score_mat[i][j] = max(0, sub, dele, ins)

            if score_mat[i][j] == 0:
                direction_mat[i][j] = STOP
            elif score_mat[i][j] == sub:
                direction_mat[i][j] = DIAG
            elif score_mat[i][j] == dele:
                direction_mat[i][j] = UP
            elif score_mat[i][j] == ins:
                direction_mat[i][j] = LEFT

            # highest score
            if score_mat[i][j] >= max_score:
                max_score = score_mat[i][j]
                max_pos = (i, j)

    # get the optimal alignment
    aseqq = []
    aseqt = []
    i, j = max_pos

    while direction_mat[i][j] != STOP:
        if direction_mat[i][j] == DIAG:
            aseqq.append(seqq[i - 1])
            aseqt.append(seqt[j - 1])
            i -= 1
            j -= 1
        elif direction_mat[i][j] == UP:
            aseqq.append(seqq[i - 1])
            aseqt.append('-')
            i -= 1
        elif direction_mat[i][j] == LEFT:
            aseqq.append('-')
            aseqt.append(seqt[j - 1])
            j -= 1

    aseqq = ''.join(reversed(aseqq))
    aseqt = ''.join(reversed(aseqt))

    return aseqq, aseqt, max_score


def join(double_hit, seqq, seqt):
    poss1, posq1, poss2, posq2 = double_hit
    length = poss2-poss1+WORD_LEN
    score = sum([aa_sub_score(seqq[posq1+i], seqt[poss1+i])
                for i in range(0, length)])
    return seqq[posq1:length], seqt[poss1:length], score, posq1, poss1


def extend(joined_alignment, seqq, seqt):
    jaseqq, jaseqt, score, posq1, poss1 = joined_alignment
    # to the left
    rlseqq = seqq[posq1::-1]
    rlseqt = seqt[poss1::-1]
    laseqq, laseqt, lascore = alignment(rlseqq, rlseqt)
    # to the right
    rseqq = seqq[posq1+len(jaseqq):]
    rseqt = seqt[poss1+len(jaseqq):]
    raseqq, raseqt, rascore = alignment(rseqq, rseqt)
    # final alignement
    aseqq = laseqq[::-1] + jaseqq + raseqq
    aseqt = laseqt[::-1] + jaseqt + raseqt
    score = score + lascore + rascore
    posq = posq1 - len(laseqq.replace('_', ''))
    poss = poss1 - len(laseqt.replace('_', ''))
    return aseqq, aseqt, posq, poss, score


def are_similar(aa1, aa2):
    """Returns True if a substitution is conservative, False otherwise."""
    for group in similar_groups:
        if aa1 in group and aa2 in group:
            return True
    return False


def print_blast_alignment(query, subject, query_start, subject_start):
    """
    Prints the alignment of two protein sequences.

    Parameters:
    query -- Query protein sequence with alignment gaps ("_").
    subject -- Subject protein sequence with alignment gaps.
    query_start -- Starting position of the query in the alignment.
    subject_start -- Starting position of the subject (db sequence) in the alignment.
    """
    line_length = 60  # Maximum number of characters per line for display
    query_pos = query_start
    subject_pos = subject_start

    for i in range(0, len(query), line_length):
        # Get current segment of the alignment
        query_slice = query[i:i + line_length]
        subject_slice = subject[i:i + line_length]

        # Generate match line: "|" for exact matches, ":" for conservative, " " for mismatches or gaps
        match_line = ''.join(
            ['|' if query_slice[j] == subject_slice[j]
             else ':' if query_slice[j] != '_' and subject_slice[j] != '_' and are_similar(query_slice[j], subject_slice[j])
             else ' '
             for j in range(len(query_slice))])

        # Print query line with one space at the beginning for alignment
        print(f"Query   {query_pos:<4}  {query_slice}\
              {query_pos + query_slice.replace('_', '').count('') - 1}")
        query_pos += len(query_slice.replace('_', ''))

        # Print match line with one less space to align with sequences
        print(f"              {match_line}")

        # Print subject line, showing gaps and using the actual position for the non-gap characters
        print(f"Subject {subject_pos:<4}  {subject_slice}\
              {subject_pos + subject_slice.replace('_', '').count('') - 1}")
        subject_pos += len(subject_slice.replace('_', ''))

        print()  # Blank line between alignment blocks


def run_blast(fastaf, dbf):
    records = list(SeqIO.parse(fastaf, "fasta"))
    if len(records) > 2:
        print("The fasta file provided as query contains more than one sequence. The first sequence of the file will be used for the blast search.")
    seq = str(records[0].seq)
    words = seq_to_words(seq)
    db_hits = db_search(dbf, words)
    del db_hits['index']
    alignments = {}
    nb_seq_aligned = 0
    for seqid, values in db_hits.items():
        max_score = 0
        for dbh in values[1:]:
            seqt = values[0]
            jalignment = join(dbh, seq, seqt)
            alignment = extend(jalignment, seq, seqt)
            score = alignment[4]
            norm_score = normalized_score(score)
            if norm_score > max_score:
                max_score = norm_score
                evalue = e_value(norm_score, len(alignment[0]))
                best_alignment = alignment + (evalue,)
        alignments[seqid] = best_alignment
        nb_seq_aligned += 1
        if nb_seq_aligned % 20 == 0:
            print(nb_seq_aligned)
    # show results
    for key, values in alignments.items():
        aseqq, aseqt, posq, poss, score, evalue = values
        print(f"Sequence : {key}, score : {score}, evalue : {evalue}")
        print_blast_alignment(aseqq, aseqt, posq, poss)
    return alignments


if __name__ == "__main__":
    # run_blast("data/P06007.fasta", "data/photosystem_viruses_archea.fasta")
    run_blast("data/smallseq.fasta", "data/smalldb.fasta")
