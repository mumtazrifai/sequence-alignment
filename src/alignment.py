import numpy as np
import pandas as pd
import sys
import argparse
import pathlib

def load_substitution_matrix(path):
    amino_acid_order = list("ARNDCQEGHILKMFPSTWYV")

    sub_mat = pd.read_csv(path, sep=" ", dtype=float, header=None)
    sub_mat.index = amino_acid_order
    sub_mat.columns = amino_acid_order

    return sub_mat

def parse_fasta(path):
    with open(path) as file:
        lines = file.readlines()
        seq_name_indicies = []

        for i in range(len(lines)):
            if ">" in lines[i]:
                seq_name_indicies.append(i)
        
        seqx_name = lines[seq_name_indicies[0]].replace(">", "").strip()
        seqx = "".join(lines[seq_name_indicies[0] + 1 : seq_name_indicies[1]]).replace("\n", "").replace("\r", "")

        seqy_name = lines[seq_name_indicies[1]].replace(">", "").strip()
        seqy = "".join(lines[seq_name_indicies[1] + 1:]).replace("\n", "").replace("\r", "")

        return (seqx, seqy, seqx_name, seqy_name)

def global_allgnment(seqx, seqy, sub_mat, d):
    alignment_mat = np.full([len(seqx) + 1, len(seqy) + 1], np.nan)
    pointer_mat = np.full([len(seqx) + 1, len(seqy) + 1, 2], np.nan, dtype=int)

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

        if not np.isnan(alignment_mat[i, j]):
            return alignment_mat[i, j] 

        both_residues_aligned = calc_score(i-1, j-1) + sub_mat.loc[seqx[i - 1], seqy[j - 1]]
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

def get_output_text(seqx, seqy, namex, namey):
    seq_start_length = max(len(namex), len(namey)) + 1

    line1 = namex + (" " * (seq_start_length - len(namex))) + "| " + seqx + "\n"
    line2 = namey + (" " * (seq_start_length - len(namey))) + "| " + seqy

    return (line1, line2)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("fasta_path", metavar="input_fasta", action="store", type=pathlib.Path, help="The path of the file containing the two sequences to be aligned in FASTA format.")
    parser.add_argument("sub_mat_path", metavar="substitution_matrix", action="store", type=pathlib.Path, help="The path of the subsitution matrix file with space-separated values. The order of the rows and columns is: \"ARNDCQEGHILKMFPSTWYV\".")
    parser.add_argument("-o", "--output", metavar="output", action="store", type=pathlib.Path, help="The path of the alignment file to be outputted. If not set, will simply print the results")

    args = parser.parse_args()

    sub_mat = load_substitution_matrix(args.sub_mat_path)
    seqx, seqy, seqx_name, seqy_name = parse_fasta(args.fasta_path)

    sys.setrecursionlimit(len(seqx) + len(seqy) + 3)

    aligned_seqx, aligned_seqy = global_allgnment(seqx, seqy, sub_mat, 8)

    output_line1, output_line2 = get_output_text(aligned_seqx, aligned_seqy, seqx_name, seqy_name)

    if args.output:
        with open(args.output, "w") as file:
            file.writelines([output_line1, output_line2])
    else:
        print(output_line1, output_line2, sep="")