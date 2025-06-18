#!/usr/bin/env python3

import csv
import os
import sys
import argparse
import dendropy
from dendropy.calculate import treecompare

# Функция загрузки дерева
def load_tree(path, taxon_namespace):
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace
    )

# Вычисление RF
def compute_rf(tree, ref):
    return treecompare.symmetric_difference(ref, tree)

# Вычисление KF
def compute_kf(tree, ref):
    tree.encode_bipartitions()
    ref.encode_bipartitions()
    ref_splits = {edge.bipartition: edge.length for edge in ref.postorder_edge_iter() if edge.bipartition and not edge.bipartition.is_trivial()}
    tree_splits = {edge.bipartition: edge.length for edge in tree.postorder_edge_iter() if edge.bipartition and not edge.bipartition.is_trivial()}
    common = set(ref_splits.keys()) & set(tree_splits.keys())
    return sum((ref_splits[b] - tree_splits[b])**2 for b in common)

# Основной скрипт
def main(combos_file, combos_dir, ref_tree, out_csv):
    ns = dendropy.TaxonNamespace()
    ref = load_tree(ref_tree, ns)
    ref.encode_bipartitions()

    with open(combos_file, encoding="utf-8") as fin, \
         open(out_csv, "w", newline="", encoding="utf-8") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout)
        writer.writerow(["combo_name", "genes", "rf_distance", "norm_rf", "kf_distance"])

        for combo_name, genes_str in reader:
            genes = genes_str.split()
            treefile = os.path.join(combos_dir, combo_name, "tree.treefile")
            if not os.path.isfile(treefile):
                sys.stderr.write(f"[WARN] отсутствует {treefile}\n")
                continue

            t = load_tree(treefile, ns)
            rf = compute_rf(t, ref)
            n_leaves = len(t.leaf_nodes())
            max_rf = 2 * (n_leaves - 3) if n_leaves >= 3 else 0
            norm_rf = rf / max_rf if max_rf > 0 else 0.0
            kf = compute_kf(t, ref)

            writer.writerow([combo_name, ";".join(genes), rf, f"{norm_rf:.4f}", kf])

    print(f"Готово — результаты записаны в {out_csv}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Вычисляет RF и KF расстояния до референсного дерева")
    p.add_argument("combos_list", help="файл со списком комбо (tab-delimited)")
    p.add_argument("combos_dir", help="корневая папка с поддиректориями комбо")
    p.add_argument("ref_tree", help="файл reference.treefile")
    p.add_argument("-o", "--output", default="distances.csv", help="путь и имя выходного CSV")
    args = p.parse_args()

    main(
        combos_file=args.combos_list,
        combos_dir=args.combos_dir,
        ref_tree=args.ref_tree,
        out_csv=args.output
    )