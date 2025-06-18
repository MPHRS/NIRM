#!/usr/bin/env python3
"""
compute_rf_distances.py

Given:
  - A tab-delimited combos list: each line “combo_name\tgene1 gene2 … geneK”
  - A directory where each combo has its IQ-TREE output (treefile named: tree.treefile)
  - A reference treefile

Outputs:
  distances.csv with columns:
    combo_name, genes…, rf_distance, max_rf, norm_rf
"""

import csv
import os
import sys
from ete3 import Tree

def load_reference(ref_tree_path):
    """Load the reference tree from the given file path."""
    ref = Tree(ref_tree_path, format=1)
    return ref

def compute_distance(tree_path, ref_tree):
    """Compute the Robinson-Foulds distance between a tree and the reference."""
    t = Tree(tree_path, format=1)
    # ETE3 returns a tuple with at least 6 elements. We only need RF distance and max RF.
    rf_results = ref_tree.robinson_foulds(t, unrooted_trees=True)
    rf = rf_results[0]
    max_parts = rf_results[1]
    return rf, max_parts

def main(combos_file, combos_dir, ref_tree_path, out_csv="distances.csv"):
    # Load reference tree
    ref = load_reference(ref_tree_path)

    # Open combos list and output CSV
    with open(combos_file) as fin, open(out_csv, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout)
        # Write header
        writer.writerow(["combo_name", "genes", "rf_distance", "max_rf", "norm_rf"]);

        # Iterate through each combo
        for combo_name, genes_str in reader:
            genes = genes_str.split()
            # Construct expected treefile path
            treefile = os.path.join(combos_dir, combo_name, "tree.treefile")
            if not os.path.isfile(treefile):
                print(f"[WARN] missing tree for {combo_name} at {treefile}", file=sys.stderr)
                continue

            # Compute RF distance
            rf, max_rf = compute_distance(treefile, ref)
            norm = rf / max_rf if max_rf > 0 else 0.0
            # Write results
            writer.writerow([combo_name, ";".join(genes), rf, max_rf, f"{norm:.4f}"])

    print(f"Written distances to {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: compute_rf_distances.py <combos_list.txt> <combos_root_dir> <reference.treefile>")
        sys.exit(1)

    combos_list = sys.argv[1]
    combos_root = sys.argv[2]
    ref_tree = sys.argv[3]
    main(combos_list, combos_root, ref_tree)
