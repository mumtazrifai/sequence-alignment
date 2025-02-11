import numpy as np
import pandas as pd

amino_acid_order = list("ARNDCQEGHILKMFPSTWYV")

score_mat = pd.read_csv("src\\score_matrices\\BLOSUM50.txt", sep=" ", dtype=float, header=None)
score_mat.index = amino_acid_order
score_mat.columns = amino_acid_order

seqx = "PAWHEAE"
seqy = "HEAGAWGHEE"

def global_allgnment(seqx, seqy, score_mat, d):
    alignment_mat = np.empty(shape=(len(seqx) + 1, len(seqy) + 1))
    pointer_mat = np.empty(shape=(len(seqx) + 1, len(seqy) + 1, 2), dtype=int)

    # Create function to recursively fill the alignment score matrix

    def calc_score(i, j):
        if i == 0 and j == 0:
            alignment_mat[0, 0] = 0
            return 0
        
        if i == 0:
            best_score = calc_score(0, j - 1) - d
            alignment_mat[0, j] = best_score
            pointer_mat[0, j, :] = (0, -1)
            return best_score
        
        if j == 0:
            best_score = calc_score(i - 1, 0) - d
            alignment_mat[i, 0] = best_score
            pointer_mat[i, 0, :] = (-1, 0)
            return best_score

        both_residues_aligned = calc_score(i-1, j-1) + score_mat.loc[seqx[i - 1], seqy[j - 1]]
        x_aligned_with_gap = calc_score(i-1, j) - d
        y_aligned_with_gap = calc_score(i, j-1) - d

        best_score = max(both_residues_aligned, x_aligned_with_gap, y_aligned_with_gap)

        alignment_mat[i, j] = best_score
        
        if best_score == both_residues_aligned:
            pointer_mat[i, j, :] = (-1, -1)
        elif best_score == x_aligned_with_gap:
            pointer_mat[i, j, :] = (-1, 0)
        else:
            pointer_mat[i, j, :] = (0, -1)

        return best_score

    calc_score(len(seqx), len(seqy))

    # Traceback

    aligned_seqx = ""
    aligned_seqy = ""

    i = pointer_mat.shape[0] - 1
    j = pointer_mat.shape[1] - 1

    while not (i == 0 and j == 0):
        index_change = tuple(pointer_mat[i, j, :])

        if index_change == (-1, -1):
            aligned_seqx = seqx[i - 1] + aligned_seqx
            aligned_seqy = seqy[j - 1] + aligned_seqy
        elif index_change == (-1, 0):
            aligned_seqx = seqx[i - 1] + aligned_seqx
            aligned_seqy = "-" + aligned_seqy
        else:
            aligned_seqx = "-" + aligned_seqx
            aligned_seqy = seqy[j - 1] + aligned_seqy

        i += index_change[0]
        j += index_change[1]

    return (aligned_seqx, aligned_seqy)


print(global_allgnment(seqx, seqy, score_mat, 8))