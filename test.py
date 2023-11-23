import sys
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

def similarity_score(a, b):
    return p if a == b else -q

def fill_scoring_matrix(seq1, seq2, W1):
    def s(i, j):
        return similarity_score(seq1[i], seq2[j])

    n = len(seq1)
    m = len(seq2)
    H = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = H[i - 1, j - 1] + s(i - 1, j - 1)
            delete = H[i - 1, j] - W1
            insert = H[i, j - 1] - W1
            H[i, j] = max(match, delete, insert, 0)

    return H

def find_max_and_traceback(scoring_matrix):
    max_value = np.max(scoring_matrix)
    max_indices = np.argwhere(scoring_matrix == max_value)[0]
    i, j = max_indices
    traceback_path = [(i, j)]

    match = scoring_matrix[i - 1, j - 1] if i > 0 and j > 0 else 0
    delete = scoring_matrix[i - 1, j] if i > 0 else 0
    insert = scoring_matrix[i, j - 1] if j > 0 else 0

    while match != 0 and delete != 0 and insert != 0:
        if match >= delete and match >= insert:
            i, j = i - 1, j - 1
        elif delete >= insert:
            i -= 1
        else:
            j -= 1

        match = scoring_matrix[i - 1, j - 1] if i > 0 and j > 0 else 0
        delete = scoring_matrix[i - 1, j] if i > 0 else 0
        insert = scoring_matrix[i, j - 1] if j > 0 else 0

        traceback_path.append((i, j))

    traceback_path.reverse()

    return max_value, traceback_path

# Example sequences
seq2 = "GTTGCTGTAGCATTCTGGTAGG"
seq1 = "AGGCATCTACCCTCCTCACGGA"

# Gap penalty
W1 = 2
p = 4
q = 4
scoring_matrix = fill_scoring_matrix(seq1, seq2, W1)

# Display the scoring matrix as a heatmap
plt.imshow(scoring_matrix, cmap='viridis', origin='lower')
plt.colorbar(label='Score')
plt.title('Scoring Matrix Heatmap')
plt.xlabel('Sequence 2')
plt.ylabel('Sequence 1')
plt.xticks(np.arange(len(seq2) + 1), [''] + list(seq2))
plt.yticks(np.arange(len(seq1) + 1), [''] + list(seq1))

# Highlight the traceback path
max_value, traceback_path = find_max_and_traceback(scoring_matrix)
traceback_path = np.array(traceback_path)
plt.plot(traceback_path[:, 1], traceback_path[:, 0], marker='o', color='red', label='Traceback Path')

plt.legend()
plt.show()
