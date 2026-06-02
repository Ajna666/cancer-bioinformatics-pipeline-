import os
import matplotlib.pyplot as plt


# ----------------------------
# FASTA READER
# ----------------------------
def read_fasta(file_path):
    sequence = []

    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())

    return "".join(sequence).upper()


# ----------------------------
# GC CONTENT
# ----------------------------
def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


# ----------------------------
# CODON ANALYSIS
# ----------------------------
def get_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq)-2, 3)]


def codon_stats(codons):
    return {
        "ATG": codons.count("ATG"),
        "TAA": codons.count("TAA"),
        "TAG": codons.count("TAG"),
        "TGA": codons.count("TGA"),
    }


# ----------------------------
# ANALYZE GENE
# ----------------------------
def analyze_gene(file_path):
    seq = read_fasta(file_path)
    codons = get_codons(seq)

    stats = codon_stats(codons)

    return {
        "length": len(seq),
        "gc": gc_content(seq),
        "start_codon": stats["ATG"],
        "stop_codons": stats["TAA"] + stats["TAG"] + stats["TGA"]
    }


# ----------------------------
# RUN PIPELINE
# ----------------------------
def run_pipeline(folder):
    results = {}

    for file in os.listdir(folder):
        if file.endswith(".fasta"):
            gene = file.replace(".fasta", "")
            path = os.path.join(folder, file)

            results[gene] = analyze_gene(path)

            print(f"{gene} analyzed ✔")

    return results


# ----------------------------
# VISUALIZATION
# ----------------------------
def plot_gc(results):
    genes = list(results.keys())
    gc_values = [results[g]["gc"] for g in genes]

    plt.figure()
    plt.bar(genes, gc_values)
    plt.title("GC Content Across Cancer Genes")
    plt.xlabel("Genes")
    plt.ylabel("GC %")

    plt.tight_layout()
    plt.savefig("plots/gc_content.png", dpi=300)
    plt.close()


# ----------------------------
# MAIN
# ----------------------------
if __name__ == "__main__":
    results = run_pipeline("data")

    print("\n=== SUMMARY ===")
    for gene, r in results.items():
        print(gene, "| GC:", round(r["gc"], 2), "| Length:", r["length"])

    plot_gc(results)

    print("\nGC analysis completed ✔")
