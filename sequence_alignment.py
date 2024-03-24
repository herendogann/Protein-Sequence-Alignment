import blosum as bl
import sys

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap_open=-10, gap_extend=-5):
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1

    score_matrix = [[0] * len_seq2 for _ in range(len_seq1)]

    for i in range(len_seq1):
        score_matrix[i][0] = i * (gap_extend if i > 0 else gap_open)

    for j in range(len_seq2):
        score_matrix[0][j] = j * (gap_extend if j > 0 else gap_open)


    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            match_mismatch_score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            diagonal_score = score_matrix[i - 1][j - 1] + match_mismatch_score
            horizontal_score = score_matrix[i][j - 1] + (gap_extend if j > 1 else gap_open)
            vertical_score = score_matrix[i - 1][j] + (gap_extend if i > 1 else gap_open)

            score_matrix[i][j] = max(diagonal_score, horizontal_score, vertical_score)


    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = len_seq1 - 1, len_seq2 - 1

    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + (gap_extend if i > 1 else gap_open):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1

        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, score_matrix[len_seq1 - 1][len_seq2 - 1]

def calculate_alignment_score(aligned_seq1, aligned_seq2, scoring_matrix, gap_open, gap_extend):
    score = 0
    in_gap = False
    gap_length = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != "-" and b != "-":
            if in_gap:
                score += gap_open + gap_extend * (gap_length - 1)
                in_gap = False
                gap_length = 0
            score += scoring_matrix[a][b]
        elif a == "-" and b == "-":
            continue
        elif a == "-":
            if in_gap:
                gap_length += 1
            else:
                in_gap = True
                gap_length = 1
        elif b == "-":
            if in_gap:
                gap_length += 1
            else:
                in_gap = True
                gap_length = 1

    if in_gap:
        score += gap_open + gap_extend * (gap_length - 1)

    return score

def calculate_identity(aligned_seq1, aligned_seq2):
    matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2))
    total = len(aligned_seq1)
    identity = (matches / total) * 100
    return matches, total, identity

def print_alignment(aligned_seq1, aligned_seq2):
    print(aligned_seq1)
    print(''.join('|' if a == b else ' ' for a, b in zip(aligned_seq1, aligned_seq2)))
    print(aligned_seq2)

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return lines[0].strip(), lines[1].strip()

def initialize_matrix(rows, cols):
    return [[0] * cols for _ in range(rows)]

def create_scoring_matrix(matrix):
    return {key1: {key2: matrix[key1][key2] for key2 in matrix[key1]} for key1 in matrix}

if len(sys.argv) != 5:
        print("Usage: python alignment.py <sequence_file> <blosumXX> <gap_open> <gap_extend>")
        sys.exit(1)

sequence_file = sys.argv[1]
blosum_number = int(sys.argv[2][6:])
gap_open = int(sys.argv[3])
gap_extend = int(sys.argv[4])

seq1, seq2 = read_fasta(sequence_file)


matrix = bl.BLOSUM(blosum_number)
scoring_matrix = create_scoring_matrix(matrix)
alignment_result = needleman_wunsch(seq1, seq2)
aligned_seq1 = alignment_result[0]
aligned_seq2 = alignment_result[1]

print_alignment(alignment_result[0], alignment_result[1])
matches, total, identity = calculate_identity(aligned_seq1, aligned_seq2)
alignment_score = calculate_alignment_score(aligned_seq1, aligned_seq2, scoring_matrix, gap_open, gap_extend)
print(f"\nAlignment score: {alignment_score}")
print(f"Identity value: {matches}/{total} ({identity:.1f}%)")
print()
