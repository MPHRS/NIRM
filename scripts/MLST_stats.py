#!/usr/bin/env python3
import numpy as np
import dendropy
from dendropy.calculate import treecompare
from pathlib import Path

# Пути к файлам (от корня проекта)
REF_TREE     = Path("final_trees/tree_all_genes.treefile")
MLST_TREE    = Path("mlst_extracted/mlst_iqtree.treefile")
RESULTS_FILE = Path("mlst_extracted/results_stats.txt")

# Создаем общее пространство таксонов
ns = dendropy.TaxonNamespace()

# Загружаем деревья с единой TaxonNamespace
ref  = dendropy.Tree.get(path=str(REF_TREE),
                         schema="newick",
                         taxon_namespace=ns,
                         preserve_underscores=True)
mlst = dendropy.Tree.get(path=str(MLST_TREE),
                         schema="newick",
                         taxon_namespace=ns,
                         preserve_underscores=True)

# 1) RF-расстояние
rf = treecompare.symmetric_difference(ref, mlst)

# 2) Средняя bootstrap-поддержка
bs_values = []
for nd in mlst.internal_nodes():
    if nd.label is not None:
        try:
            bs_values.append(float(nd.label))
        except ValueError:
            pass
mean_bs = np.mean(bs_values) if bs_values else 0.0

# Пишем результаты
RESULTS_FILE.parent.mkdir(exist_ok=True, parents=True)
with open(RESULTS_FILE, "w") as f:
    f.write(f"RF distance: {rf}\n")
    f.write(f"Mean bootstrap support: {mean_bs:.1f}%\n")

print("Done. See", RESULTS_FILE)
