
Cancer Gene Bioinformatics Pipeline

A Python-based computational biology project for analyzing cancer-related genes using real DNA sequences.

Features

- GC content analysis
- Codon-level statistics
- Multi-gene comparison
- Phylogenetic tree construction
- Data visualization

Data Source

Sequences are retrieved from the :contentReference[oaicite:0]{index=0} in FASTA format.

:contentReference[oaicite:1]{index=1}  

Outputs

- GC content bar chart
- Phylogenetic tree (sequence similarity-based clustering)

Phylogenetic Analysis

![Tree](plots/phylogenetic_tree.png)

Genes are clustered based on sequence similarity to infer structural relationships.

How to Run

```bash
pip install -r requirements.txt
python src/multi_gene_system.py
python src/phylo_tree.py
