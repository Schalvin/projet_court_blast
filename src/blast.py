import sys
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
GAP_PENALTY = -11
GAP_EXT_PENALTY = -1

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
    """Get the score of substitution depending on the provided MATRIX

    Args:
        aa1 (str): first amino acid
        aa2 (str): second amino acid

    Returns:
        float: the substitution score
    """
    i = MATRIX.alphabet.index(aa1)
    j = MATRIX.alphabet.index(aa2)
    score = MATRIX[i][j]
    return score

# def close_words_gen(word):
#    word


def seq_to_words(seq):
    """Decomposition of a sequence into words of previously specified length

    Args:
        fastaf (str): path to fasta file containing the sequence to be searched

    Returns:
        dict: keys are the words obtained from the sequence, and value is the position in 
        the sequence of the first letter of the word
    """

    words = {}
    for i in range(0, len(seq)-WORD_LEN+1):
        word = seq[i:i+WORD_LEN]
        if word not in words:
            words[word] = []
        score = sum([aa_sub_score(aa, aa) for aa in word])
        if score > 11:
            words[word].append(i)
    return words


def double_hit_finder(words, seq):
    """Get all non over-lapping, double hits with maximum distance 
    between first position of hit of MAX_HITS_GAP_LEN.

    Args:
        words (dict): keys are words of length constant WORD_LEN, 
        values are positions in the query sequence. 
        generated using the seq_to_words() function
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
    for i in range(0, len(seq)):
        wordseq = seq[i:i+WORD_LEN]
        if wordseq in words:
            diags = [j-i for j in words[wordseq]]
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


def double_hit_fuser(double_hits):
    """Fuses together clusters of double hits on same diagonal

    Args:
        double_hits (tuple): output from double_hit_finder()

    Returns:
        tuple: (
            - int : position of first hit on subject sequence
            - int : position of first hit on query sequence
            - int : position of last hit on subject sequence
            - int : position of last hit on query sequence)
    """
    i = 2
    while i < len(double_hits):
        poss1_1, posq1_1, poss2_1, posq2_1 = double_hits[i-1]
        poss1_2, posq1_2, poss2_2, posq2_2 = double_hits[i]
        if (posq2_1 == posq1_2) & (poss2_1 == poss1_2):
            double_hits[i-1] = (poss1_1, posq1_1, poss2_2, posq2_2)
            del double_hits[i]
        else:
            i += 1
    return double_hits


def db_search(dbf, words):
    """Search of words in a sequence database in fasta format

    Args:
        dbf (str): path to fasta file containing the db
        words (dict): keys are words of length constant WORD_LEN, 
        values are positions in the query sequence.
        generated using the seq_to_words() function

    Returns:
        dict: keys are the sequence IDs, 
        values are lists containing the target sequence, 
        followed by double hits found in tuples with
            poss1 - position on sequence in first hit
            posq1 - position on query in first hit
            poss2 - position on sequence in second hit
            posq2 - position on query in second hit
    """
    results = {'index': ('poss1', 'posq1', 'poss2', 'posq2')}
    for dbrecord in SeqIO.parse(dbf, "fasta"):
        dbh = double_hit_finder(words, str(dbrecord.seq))
        result = double_hit_fuser(dbh)
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
    """Aligns provided sequence based on Smith and Waterman algorithm
    The best local alignment is kept

    Args:
        seqq (str): query sequence
        seqt (str): subject sequence

    Returns:
        tuple: contains the query sequence alignment sequence (gaps as '_'), 
                the subject sequence alignment sequence (gaps as '_')
                the score of the alignment
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

            score_mat[i][j] = max(sub, dele, ins)

            if score_mat[i][j] == sub:
                direction_mat[i][j] = DIAG
            elif score_mat[i][j] == dele:
                direction_mat[i][j] = UP
            elif score_mat[i][j] == ins:
                direction_mat[i][j] = LEFT

            # highest score
            if score_mat[i][j] >= max_score:
                max_score = score_mat[i][j]
                max_pos = (i, j)
    aseqq = []
    aseqt = []
    score = score_mat[i][j]
    if max_score > 0:
        # get the optimal alignment

        i, j = max_pos

        while (i != 0) & (j != 0):
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

    return aseqq, aseqt, score


def join(double_hit, seqq, seqt):
    """Joins the hits of a double hit and gives their alignment score

    Args:
        double_hit (tuple): output from double_hit_finder()
        seqq (str): query sequence
        seqt (str): target sequence

    Returns:
        tuple: contains 
                - str : aligned sequence of query sequence
                - str : aligned sequence of subject sequence
                - int : alignment score
                - int : position of alignment on query sequence
                - int : position of alignment on subject sequence
    """
    poss1, posq1, poss2, posq2 = double_hit
    length = poss2-poss1+WORD_LEN
    score = sum([aa_sub_score(seqq[posq1+i], seqt[poss1+i])
                for i in range(0, length)])
    return seqq[posq1:posq1+length], seqt[poss1:poss1+length], score, posq1, poss1


def extend(joined_alignment, seqq, seqt):
    """Extends an existing alignment with the Smith and Waterman algorithm

    Args:
        joined_alignment (tuple): output from join() function
        seqq (str): query sequence
        seqt (str): target sequence

    Returns:
        tuple: contains 
                - str : aligned sequence of query sequence
                - str : aligned sequence of subject sequence
                - int : position of alignment on query sequence
                - int : position of alignment on subject sequence
                - int : alignment score
    """
    jaseqq, jaseqt, score, posq1, poss1 = joined_alignment
    # to the left
    rlseqq = seqq[max(posq1-1, 0)::-1]
    rlseqt = seqt[max(poss1-1, 0)::-1]
    if (len(rlseqq) < 2) or (len(rlseqt) < 2):
        laseqq, laseqt, lascore = '', '', 0
    else:
        laseqq, laseqt, lascore = alignment(rlseqq, rlseqt)
    # to the right
    rseqq = seqq[posq1+len(jaseqq):]
    rseqt = seqt[poss1+len(jaseqq):]
    if (len(rseqq) < 2) or (len(rseqt) < 2):
        raseqq, raseqt, rascore = '', '', 0
    else:
        raseqq, raseqt, rascore = alignment(rseqq, rseqt)
    # final alignement
    aseqq = laseqq[::-1] + jaseqq + raseqq
    aseqt = laseqt[::-1] + jaseqt + raseqt
    score = score + lascore + rascore
    posq = posq1 - len(laseqq.replace('_', ''))
    poss = poss1 - len(laseqt.replace('_', ''))
    return aseqq, aseqt, posq, poss, score


def are_similar(aa1, aa2):
    """Returns True if a substitution is conservative, False otherwise.

    Args:
        aa1 (str): first amino acid
        aa2 (str): second amino acid

    Returns:
        boolean : True or False depending on similarity
    """
    for group in similar_groups:
        if aa1 in group and aa2 in group:
            return True
    return False


def order_by_score(alignments):
    """Sort alignments depending on score (decreasing order) 

    Args:
        alignments (dict): keys are sequence IDs, values are a tuple containing :
                - str : aligned sequence of query sequence
                - str : aligned sequence of subject sequence
                - int : position of alignment on query sequence
                - int : position of alignment on subject sequence
                - int : alignment score
                - float : normalized alignment score
                - float : evalue
    Returns:
        dict: same dict as given but ordered in decreasing score values
    """
    trans_list = []
    for key, values in alignments.items():
        trans_list.append(tuple((key, values[4])))
    trans_list.sort(key=lambda a: a[1], reverse=True)
    return trans_list


def print_blast_alignment(query, subject, query_start, subject_start):
    """Prints the alignment of two protein sequences.

    Args:
        query (str) : query protein sequence with alignment gaps ("_").
        subject (str) : subject protein sequence with alignment gaps.
        query_start (int) : starting position of the query in the alignment.
        subject_start (int) : starting position of the subject (db sequence) in the alignment.
    """
    line_length = 60  # Maximum number of characters per line for display
    query_pos = query_start+1
    subject_pos = subject_start+1

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
    """Main function of the programm that runs the whole BLAST process

    Args:
        fastaf (str): path to the query sequence fasta file
        dbf (str): path to the database in fasta format

    Returns:
        dict: alignments with sequence IDs as keys an values contain a tuple with :
                - str : aligned sequence of query sequence
                - str : aligned sequence of subject sequence
                - int : position of alignment on query sequence
                - int : position of alignment on subject sequence
                - int : alignment score
                - float : normalized alignment score
                - float : evalue
    """

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
                best_alignment = alignment + (norm_score, evalue)
        alignments[seqid] = best_alignment
        nb_seq_aligned += 1
        if nb_seq_aligned % 20 == 0:
            print(f"Please be patient! Already \
                {nb_seq_aligned} sequences aligned!")
    # show results
    ordered_ids = order_by_score(alignments)
    print()
    for ids in ordered_ids[0:99]:
        aseqq, aseqt, posq, poss, score, norm_score, evalue = alignments[ids[0]]
        print(f"Sequence : {ids[0]}, score : \
            {norm_score:.2f} bits ({score}), evalue : {evalue:.2f}")
        print_blast_alignment(aseqq, aseqt, posq, poss)
    print()
    for ids in ordered_ids:
        print(f"{ids[0]},{alignments[ids[0]][5]:.2f},{alignments[ids[0]][6]}")

    return alignments


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Too little or too many arguments provided,\n\
            Please provide the fasta file of the sequence to search, \
                followed by the  fasta file of th db")
    else:
        fastaf = sys.argv[1]
        dbf = sys.argv[2]
        run_blast(fastaf, dbf)
