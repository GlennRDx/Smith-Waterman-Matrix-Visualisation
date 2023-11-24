import numpy as np
import matplotlib.pyplot as plt
"""
Author: Glenn Ross-Dolan
Date Created: 24-11-2023
Description: This file contains a number of functions which produces a 
Smith-Waterman scoring matrix from two sequences using numpy and visualises 
it through matplotlib.
"""


#-----------------------------------------------------------------------------

def smith_waterman(seq1, seq2, gap_pen, match_score, mismatch_score):

    """
    Perform Smith-Waterman local sequence alignment.

    Args:
        seq1: Reference sequence
        seq2: Comparison sequence
        gap_pen:   Gap penalty

    Returns:
        np.ndarray: Scoring matrix
    """

    row_len = len(seq2)
    col_len = len(seq1)
    scoring_matrix = np.zeros((row_len + 1, col_len + 1), dtype=int)

    for row in range(1, row_len + 1):
        for col in range(1, col_len + 1):
            if seq2[row - 1] == seq1[col - 1]:
                match = scoring_matrix[row - 1, col - 1] + match_score
            else:
                match = scoring_matrix[row - 1, col - 1] + mismatch_score
            delete = scoring_matrix[row - 1, col] - gap_pen
            insert = scoring_matrix[row, col - 1] - gap_pen
            scoring_matrix[row, col] = max(match, delete, insert, 0)

    return scoring_matrix

#-----------------------------------------------------------------------------

def traceback(scoring_matrix):

    """
    A function which calculates the indices of each point during
    the traceback step.

    Args: 
        np.ndarray: Scoring matrix
    
    Returns:
        int: max_value
        list: traceback_path

    """

    max_value = np.max(scoring_matrix)
    max_indices = np.argwhere(scoring_matrix == max_value)[0]
    row, col = max_indices
    traceback_path = [(row, col)]

    
    match = scoring_matrix[row - 1, col - 1] if row > 0 and col > 0 else 0
    delete = scoring_matrix[row - 1, col] if row > 0 else 0
    insert = scoring_matrix[row, col - 1] if col > 0 else 0

    while match != 0 and delete != 0 and insert != 0:
        if match >= delete and match >= insert:
            row, col = row - 1, col - 1
        elif delete >= insert:
            row -= 1
        else:
            col -= 1

        match = scoring_matrix[row - 1, col - 1] if row > 0 and col > 0 else 0
        delete = scoring_matrix[row - 1, col] if row > 0 else 0
        insert = scoring_matrix[row, col - 1] if col > 0 else 0

        traceback_path.append((row, col))

    traceback_path.reverse()

    return max_value, traceback_path

#-----------------------------------------------------------------------------

def visualise_matrix(scoring_matrix):

    """
    This function displays the scoring matrix as a heat map
    Args:
        scoring_matrix
    Returns:
        Nothing but it produces a matplotlib plot
    """

    plt.imshow(scoring_matrix, cmap='magma', origin='upper')
    plt.colorbar(label='Score')
    plt.title('Scoring Matrix Heatmap')
    plt.xlabel('Query Sequence (seq2)')
    plt.ylabel('Reference Sequence (seq1)')
    plt.xticks(np.arange(0, len(seq1),step = int(0.249999*len(seq1))))
    plt.yticks(np.arange(0, len(seq2), step = int(0.249999*len(seq2))))

    # Highlight the traceback path
    traceback_path = traceback(scoring_matrix)[1]
    traceback_path = np.array(traceback_path)
    plt.plot(traceback_path[:, 1], traceback_path[:, 0], marker='.', \
             color='white', label='Traceback Path')


    plt.legend()
    plt.show()

#-----------------------------------------------------------------------------

# Reference Sequence
seq1 = "TTTAGCATGCGCAT"

# Query Sequence
seq2 = "ATAGACGACAT"

# Call the functions
matrix_1 = smith_waterman(seq1, seq2, gap_pen = 2, \
                          match_score = 3, mismatch_score = -3)

visualise_matrix(matrix_1)
