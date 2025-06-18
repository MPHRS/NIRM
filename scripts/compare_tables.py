#!/usr/bin/env python3
import pandas as pd
from ete3 import Tree

def load_bootstrap(path):
    return pd.read_csv(path, sep='\t')

def compute_norm_rf(ref, tree_path):
    t = Tree(tree_path, format=1)
    rf, max_rf, *_ = ref.robinson_foulds(t, unrooted_trees=True)
    return rf / max_rf if max_rf > 0 else 0.0

def load_top_variable(sorted_entropy_file):
    ogs = []
    with open(sorted_entropy_file) as f:
        for line in f:
            path, _ = line.split()
            og = path.rsplit('/',1)[1].replace('.aln','')
            ogs.append(og)
    return ogs

def load_top_ml(top_ml_file):
    with open(top_ml_file) as f:
        return [line.strip() for line in f if line.strip()]

def main():
    # 1) Загрузка данных
    bs_var = load_bootstrap('bootstrap_summary.txt')
    bs_ml  = load_bootstrap('bootstrap_summary_ml.txt')
    ref    = Tree('final_trees/tree_all_genes.treefile', format=1)
    top_var = load_top_variable('tmp_files/sorted_entropy.txt')
    top_ml  = load_top_ml('tmp_files/top_ml_genes.txt')

    # 2) Добавляем расстояния
    bs_var['Distance_Var'] = bs_var['Genes'].apply(
        lambda n: compute_norm_rf(ref, f'final_trees/tree_{n}genes.treefile')
    )
    bs_ml['Distance_ML'] = bs_ml['Genes'].apply(
        lambda n: compute_norm_rf(ref, f'final_trees_ml/tree_{n}.treefile')
    )

    # 3) Подготовка словарей N→список OG
    var_dict = {n: top_var[:n] for n in bs_var['Genes']}
    ml_dict  = {n: top_ml[:n]  for n in bs_ml['Genes']}

    # 4) Объединение таблиц
    bs_var = bs_var.rename(columns={'Average_Bootstrap':'Bootstrap_Var'})
    bs_ml  = bs_ml.rename(columns={'Average_Bootstrap':'Bootstrap_ML'})
    df = pd.merge(
        bs_var[['Genes','Bootstrap_Var','Distance_Var']],
        bs_ml[['Genes','Bootstrap_ML','Distance_ML']],
        on='Genes', how='outer'
    ).sort_values('Genes').reset_index(drop=True)

    # 5) Добавляем списки ортогрупп
    df['OGs_Var'] = df['Genes'].map(lambda n: ";".join(var_dict.get(n, [])))
    df['OGs_ML']  = df['Genes'].map(lambda n: ";".join(ml_dict.get(n, [])))

    # 6) Сохраняем и выводим
    df.to_csv('comparison_summary_with_OGs.csv', index=False)
    print(df.to_string(index=False))
    print("\nСохранено в comparison_summary_with_OGs.csv")

if __name__ == '__main__':
    main()
