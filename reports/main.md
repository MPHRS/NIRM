
## 1. Сбор и загрузка данных

* **Команды:**

  * `ncbi-genome-download bacteria -g "Neisseria gonorrhoeae" --assembly-levels complete,chromosome --formats genbank,cds-fasta -n … > tmp_files/list_assemblies`
  * `grep "WHO" tmp_files/list_assemblies > tmp_files/list_who_assemblies`
  * Формирование `tmp_files/who_accessions.txt` и повторный `ncbi-genome-download` с `--assembly-accessions` → `data/raw/`
  * Распаковка: `find data/raw -name "*.gz" -exec gunzip {} \;`

---

## 2. Парсинг метаданных и маппинг штаммов

* **Скрипт:** `scripts/parse_metadata.py <assembly_summary.txt> <data/genomes/refseq/bacteria/>`

  * Читает `assembly_summary.txt`, фильтрует нужные GCF\_… сборки и сохраняет `mapping.csv` с колонками `assembly_accession,biosample,strain`.

---

## 3. Аннотация FASTA-заголовков

* **Скрипт:** `scripts/update_headers.py <cds_fasta_file.fna> mapping.csv`

  * По `mapping.csv` добавляет к каждому `SeqRecord.id` поля `|biosample|strain`, результат в `*_annotated.fna`.

---

## 4. Инференс ортогрупп

* **Команда:**

  ```bash
  mkdir orthofinder_input_annotated
  find … -name "*_annotated.fna" | while read f; do
    mv "$f" orthofinder_input_annotated/
  done
  orthofinder -f orthofinder_input_annotated/ --dna
  ```

  * Результат: каталог `OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences` с .fa каждого OG.

---

## 5. Выравнивание и расчёт вариабельности

* **Мультишаговый цикл:**

  ```bash
  mkdir -p alignments variability
  find … -name "*.fa" | parallel '
    og=$(basename {} .fa)
    mafft --auto {} > alignments/${og}.aln
    python3 scripts/entropy.py alignments/${og}.aln > variability/${og}.entropy
  '
  cat variability/*.entropy | sort -k2 -nr > tmp_files/sorted_entropy.txt
  ```
* **Скрипт:** `scripts/entropy.py <alignment.fasta>` — вычисляет среднюю энтропию (Шеннона).

---

## 6. Конкатенация топ‑N самых вариабельных генов

* **Шелл‑скрипт:** `scripts/conc_var.sh`

  * Считывает `tmp_files/sorted_entropy.txt`, создаёт папки `concat_align/variant_{1..9}`, копирует в них первые N `.aln`.
* **Скрипт:** `scripts/concat_alignments.py <concat_dir> mapping.csv`

  * Для каждого `variant_*` конкатенирует выравнивания в `concat.fasta`.

---

## 7. Построение филогенетических деревьев

* **Шелл‑скрипты:**

  * `trees.sh`: для каждого N генов запускает `iqtree -s concat.fasta -m MFP -bb 1000 … -pre trees_temp/…` и копирует `.treefile` в `final_trees/`.
  * `calc_bootstrap.sh`: собирает `bootstrap_summary.txt` с усреднёнными значениями бутстрэпа.

* **Референс:** также строится дерево по всем генам (`concat_align/all/concat.fasta` → `final_trees/tree_all_genes.treefile`).

---

## 8. ML‑аналитика: случайные комбинации и фичи

1. **Генерация случайных комбинаций OG:**

   * `scripts/random_combos_pipeline.sh` → папка `random_combos/{combo_name}/concat.fasta` + IQ‑TREE + лог.
2. **Расчёт расстояний до референсного дерева:**

   * `scripts/compute_rf_distances.py combos_list.txt random_combos final_trees/tree_all_genes.treefile` → `distances.csv` (RF, норм. RF).
3. **Генерация признаков (feature engineering):**

   * `scripts/features.py --align-dir alignments --out-csv features.csv` → таблица OG×признаки (энтропия, SNP‑плотность, #аллелей, gap‑доля и др.).
4. **Запуск ML‑пайплайна:**

   * `scripts/run_ml_pipeline.sh -t tmp_files/top_ml_genes.txt -i alignments -c concat_ml -r final_trees_ml -m trees_ml_temp -b bootstrap_summary_ml.txt`
   * Аналогично для K-F distance: `-t tmp_files/top_ml_kf.txt … bootstrap_summary_ml_kf.txt`.

---

## 9. NG‑MAST (porB + tbpB)

1. **Экстракция генов:**

   * **Скрипт:** `scripts/NGMAST.py <root_dir> <output_dir>` → для каждого штамма `porB`+`tbpB` → `porB_tbpB_per_strain/*.fasta`.
2. **Конкатенация, выравнивание и дерево:**

   * **Скрипт:** `scripts/align_and_tree.py <input_dir> porB_tbpB.fasta porB_tbpB.aln porB_tbpB_tree`

     * делает `mafft`, `iqtree`.

