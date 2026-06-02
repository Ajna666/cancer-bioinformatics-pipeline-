import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram


# ----------------------------
# FASTA READER
# ----------------------------
def read_fasta(file_path):
    seq = []

    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())

    return "".join(seq).upper()


# ----------------------------
# SIMILARITY
# ----------------------------
def similarity(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / max(len(seq1), len(seq2))


# ----------------------------
# BUILD MATRIX
# ----------------------------
def build_matrix(folder):
    names = []
    sequences = []

    for file in os.listdir(folder):
        if file.endswith(".fasta"):
            names.append(file.replace(".fasta", ""))
            sequences.append(read_fasta(os.path.join(folder, file)))

    n = len(sequences)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            matrix[i][j] = 1 - similarity(sequences[i], sequences[j])

    return names, matrix


# ----------------------------
# PLOT TREE
# ----------------------------
def plot_tree(names, matrix):
    linked = linkage(matrix, method="average")

    plt.figure(figsize=(8, 5))

    dendrogram(linked, labels=names)

    plt.title("Phylogenetic Tree of Cancer Genes")
    plt.tight_layout()

    plt.savefig("plots/phylogenetic_tree.png", dpi=300)
    plt.close()


# ----------------------------
# RUN
# ----------------------------
if __name__ == "__main__":
    names, matrix = build_matrix("data")
    plot_tree(names, matrix)

    print("Phylogenetic tree completed ✔")
