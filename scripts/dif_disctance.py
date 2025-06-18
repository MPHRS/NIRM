#!/usr/bin/env python3
"""
compute_distances_py.py

Вычисляет для каждого дерева:
  - RF (Robinson-Foulds)
  - norm_RF
  - KF (Kuhner–Felsenstein, sum of squared branch-length differences for common splits)
  - BHV (геодезическое дерево-пространство), если установлен tree_geodist

Зависимости:
  pip install dendropy
  # (по желанию для BHV) pip install tree-geodist
"""

import csv
import os
import sys
import argparse

import dendropy
from dendropy.calculate import treecompare

# Попытка импортировать geodesic_distance для BHV
try:
    from tree_geodist import geodesic_distance
    _BHV_AVAILABLE = True
except ImportError:
    _BHV_AVAILABLE = False
    sys.stderr.write("[WARN] tree_geodist не найден: BHV-расстояние вычисляться не будет\n")


def load_tree(path, taxon_namespace):
    """Читает дерево из Newick-файла с общим таксономическим пространством."""
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace
    )


def compute_rf(tree, ref):
    """Возвращает RF расстояние (число несовпадающих сплитов)."""
    return treecompare.symmetric_difference(ref, tree)


def compute_kf(tree, ref):
    """
    Вычисляет KF расстояние как сумму квадратов разностей длин ветвей
    для общих сплитов.
    """
    # Убедимся, что бипартиции закодированы
    tree.encode_bipartitions()
    ref.encode_bipartitions()
    # Собираем мапы: bipartition -> branch length
    ref_splits = {
        edge.bipartition: edge.length
        for edge in ref.postorder_edge_iter()
        if edge.bipartition and not edge.bipartition.is_trivial()
    }
    tree_splits = {
        edge.bipartition: edge.length
        for edge in tree.postorder_edge_iter()
        if edge.bipartition and not edge.bipartition.is_trivial()
    }
    # Общие сплиты
    common = set(ref_splits.keys()) & set(tree_splits.keys())
    # Сумма квадратов разницы длин
    kf = sum((ref_splits[b] - tree_splits[b])**2 for b in common)
    return kf


def main(combos_file, combos_dir, ref_path, out_csv):
    # Создаём единое таксономическое пространство
    ns = dendropy.TaxonNamespace()

    # Загрузка и подготовка опорного дерева
    ref = load_tree(ref_path, ns)
    ref.encode_bipartitions()

    with open(combos_file, encoding="utf-8") as fin, \
         open(out_csv, "w", newline="", encoding="utf-8") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout)
        # Заголовок
        header = [
            "combo_name", "genes",
            "rf_distance", "norm_rf",
            "kf_distance"
        ]
        if _BHV_AVAILABLE:
            header.append("bhv_distance")
        writer.writerow(header)

        for combo_name, genes_str in reader:
            genes = genes_str.split()
            treefile = os.path.join(combos_dir, combo_name, "tree.treefile")
            if not os.path.isfile(treefile):
                sys.stderr.write(f"[WARN] отсутствует {treefile}\n")
                continue

            # загрузка дерева с тем же TaxonNamespace
            t = load_tree(treefile, ns)

            # 1) RF
            rf = compute_rf(t, ref)

            # нормированный RF
            n_leaves = len(t.leaf_nodes())
            max_rf = 2 * (n_leaves - 3) if n_leaves >= 3 else 0
            norm_rf = rf / max_rf if max_rf > 0 else 0.0

            # 2) KF
            kf = compute_kf(t, ref)

            row = [
                combo_name,
                ";".join(genes),
                rf,
                f"{norm_rf:.4f}",
                kf
            ]

            # 3) BHV, если доступно
            if _BHV_AVAILABLE:
                try:
                    bhv = geodesic_distance(ref, t)
                except Exception as e:
                    sys.stderr.write(f"[ERROR] BHV не посчитан для {combo_name}: {e}\n")
                    bhv = ""
                row.append(bhv)

            writer.writerow(row)

    print(f"Готово — результаты записаны в {out_csv}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Вычисляет RF, KF и (опционально) BHV расстояния до референсного дерева"
    )
    p.add_argument("combos_list", help="файл со списком комбо (tab-delimited)")
    p.add_argument("combos_dir",  help="корневая папка с поддиректориями комбо")
    p.add_argument("ref_tree",    help="файл reference.treefile")
    p.add_argument(
        "-o", "--output",
        default="distances.csv",
        help="путь и имя выходного CSV (default: distances.csv)"
    )
    args = p.parse_args()

    main(
        combos_file=args.combos_list,
        combos_dir=args.combos_dir,
        ref_path=args.ref_tree,
        out_csv=args.output
    )
