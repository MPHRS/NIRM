#!/usr/bin/env python3
"""
NGMAST_qual.py: Оценивает дерево по bootstrap (grep-style) и RF-distance к референсному.
Requirements: ete3
Usage:
  python3 NGMAST_qual.py <treefile> <reference_treefile>
"""

import sys
import re
from ete3 import Tree


def strip_prefix(tree, prefix="strain_"):
    """Убирает префикс из имён листьев."""
    for leaf in tree.iter_leaves():
        if leaf.name.startswith(prefix):
            leaf.name = leaf.name[len(prefix):]


def load_tree_dynamic(treefile):
    """
    Загружает дерево с ETE3, автоматически подбирая формат и возвращает (Tree, has_support)
    has_support=True, если в Newick найдены bootstrap-значения ')<digits>'.
    """
    with open(treefile) as f:
        newick = f.read()
    has_supp = True if re.search(r"\)\d+", newick) else False
    fmt = 3 if has_supp else 1
    return Tree(newick, format=fmt), has_supp


def compute_bootstrap_from_newick(treefile):
    """Парсит Newick-файл и находит все bootstrap-значения по grep-способу."""
    with open(treefile) as f:
        newick = f.read()
    matches = re.findall(r"\)([0-9]+)", newick)
    values = [int(x) for x in matches]
    return (sum(values) / len(values)) if values else 0.0


def main(treefile, reference_treefile):
    # 1) Загружаем оба дерева с динамическим выбором формата
    tree, has_supp_tree = load_tree_dynamic(treefile)
    ref_tree, has_supp_ref = load_tree_dynamic(reference_treefile)

    # 2) Убираем префикс 'strain_'
    strip_prefix(tree)
    strip_prefix(ref_tree)

    # 3) Bootstrap: grep-style parsing (если есть)
    mean_bs = compute_bootstrap_from_newick(treefile) if has_supp_tree else None

    # 4) RF-distance
    rf, max_rf, common_leaves, *_ = tree.robinson_foulds(ref_tree, unrooted_trees=True)
    norm_rf = rf / max_rf if max_rf > 0 else 0.0

    # 5) Вывод результатов
    print(f"Input tree              : {treefile}")
    print(f"Reference tree          : {reference_treefile}\n")
    print(f"Bootstrap present in input tree? {has_supp_tree}")
    print(f"Bootstrap present in ref tree?   {has_supp_ref}")
    if has_supp_tree:
        print(f"Average bootstrap support (grep-style): {mean_bs:.2f}")
    else:
        print("No bootstrap values found in input tree.")
    print(f"Common leaves           : {len(common_leaves)}")
    print(f"Raw RF distance         : {rf}")
    print(f"Maximum possible splits : {max_rf}")
    print(f"Normalized RF distance  : {norm_rf:.3f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
