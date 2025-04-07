# Скачивание данных


```{Bash}
ncbi-genome-download bacteria -g "Neisseria gonorrhoeae" --assembly-levels complete,chromosome --formats genbank,fasta,protein-fasta --output-folder ./genomes -v
```

## Распакуем

```{Bash}
gunzip -r genomes/
```

# Orthofinder
## подготовка данных 

если у нас не скачиваются с ncbi .faa protein-fasta файлы, то можно сделать при помощи prodigal

---
надо подумать какие именно входные файлы у нас и если что оптимизирвать предобработку

---


	а пока запустим orthofinder


и соберём их в одну папку, так как orthofinder требует тчобы они лежали в одной директории 

```
mkdir -p orthofinder_input

find genomes/ -type f -name "*.faa" -exec cp {} orthofinder_input/ \;
```


```
 orthofinder -f orthofinder_input/
```

сейчас я запустил для 25 сборок

заняло около 21 минуты

---

почитать как именно работает orthofinder

---


# Выравнивание ортогрупп

```{Bash}
#!/bin/bash
ALIGN_DIR="alignments"
VAR_DIR="variability"
SCOS="./orthofinder_input_2/OrthoFinder/Results_Apr07/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
SEQ_DIR="./orthofinder_input_2/OrthoFinder/Results_Apr07/Single_Copy_Orthologue_Sequences"

mkdir -p ${ALIGN_DIR} ${VAR_DIR}

parallel -j $(nproc) --bar --colsep '\t' '
    og={1};
    input="'${SEQ_DIR}'/${og}.fa";
    aln="'${ALIGN_DIR}'/${og}.aln";
    entropy="'${VAR_DIR}'/${og}.entropy";

    if [ -f "${input}" ]; then
        mafft --auto --thread 2 "${input}" > "${aln}" 2> /dev/null
        python3 entropy.py "${aln}" > "${entropy}"
    else
        echo "Файл ${input} не найден! Пропускаем..."
    fi
' :::: "${SCOS}"  # Исправлено здесь

echo "Готово! Результаты:"
echo "Alignments: ${ALIGN_DIR}"
echo "Variability: ${VAR_DIR}"
```

тут присутсвует вызов скрипта entropy.py

```{Python}
# File: entropy.py
from Bio import AlignIO
import numpy as np
import sys

def calculate_entropy(alignment):
    entropy = []
    for col in range(alignment.get_alignment_length()):
        counts = {}
        for rec in alignment[:, col]:
            aa = rec.upper()
            if aa not in ["-", "X", "J", "Z", "B"]:
                counts[aa] = counts.get(aa, 0) + 1
        total = sum(counts.values())
        if total == 0:
            entropy.append(0.0)
            continue
        s = -sum((c/total) * np.log2(c/total) for c in counts.values())
        entropy.append(s)
    return np.mean(entropy)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <alignment.fasta>")
        sys.exit(1)

    aln = AlignIO.read(sys.argv[1], "fasta")
    mean_entropy = calculate_entropy(aln)
    print(f"{sys.argv[1]}\t{mean_entropy:.4f}")
```


---
вопрос по энтропии шеннона -  у нас составная характеристика, а вот именно последовательность не учитывается, подумать как можно сортировать иначе с учётом последовтельностей или локальных вариаций или такого достаточно будет

---



после этого соберем отсортированно все ортогруппы
```
cat variability/*.entropy | sort -k2 -nr > sorted_entropy.txt
```


	(GR_data) mhprs@Seraph:~/PISH/grant/pipeline_final$ head sorted_entropy.txt
		alignments/OG0000127.aln        0.3013
		alignments/OG0001315.aln        0.2599
		alignments/OG0001314.aln        0.2058
		alignments/OG0000134.aln        0.2056
		alignments/OG0001094.aln        0.2041
		alignments/OG0001280.aln        0.2012
		alignments/OG0001416.aln        0.2008
		alignments/OG0001485.aln        0.1771
		alignments/OG0000990.aln        0.1731
		alignments/OG0000140.aln        0.1461

и их отдельно соберем для построения дерева 

```{bash}

N=5
mkdir -p top_variable_alignments



head -n 5 sorted_entropy.txt | while read -r line; do
    full_path=$(echo "$line" | awk '{print $1}')  # Берём полный путь (alignments/XXX.aln)
    og=$(basename "$full_path" .aln)              # Извлекаем только имя ортогруппы (XXX)
    entropy=$(echo "$line" | awk '{print $2}')    # Отдельно сохраняем энтропию
    
    # Копируем, используя полный путь из файла
    cp "$full_path" top_variable_alignments/
    
    # Проверяем, успешно ли скопировалось
    if [ -f "top_variable_alignments/$(basename "$full_path")" ]; then
        echo "Успешно: $og.aln (Энтропия: $entropy)"
    else
        echo "Ошибка: $full_path не найден!" >&2
    fi
done

```

	(GR_data) mhprs@Seraph:~/PISH/grant/pipeline_final/top_variable_alignments$ ls
		OG0000127.aln  
		OG0000134.aln 
		OG0001094.aln  
		OG0001314.aln  
		OG0001315.aln

теперь нам надо их собрать воедино для построения дерева


```{Bash}
AMAS.py concat -i top_variable_alignments/*.aln -f fasta -d aa -t concat_align/concatenated_alignment.fa -p concat_align/partitions.txt
```
	/concat_align$ ls
		concatenated_alignment.fa  partitions.txt
	
	- `concatenated_alignment.fa` — итоговое выравнивание в формате FASTA.
	- `partitions.txt` — файл с разметкой (какие гены откуда взяты).

и строим дерево 

```{Bash}
iqtree -s concatenated_alignment.fa -m TEST -bb 1000
```


получим 
	concat_align$ ls
	concatenated_alignment.fa          concatenated_alignment.fa.iqtree    concatenated_alignment.fa.splits.nex
	concatenated_alignment.fa.bionj    concatenated_alignment.fa.log       concatenated_alignment.fa.treefile
	concatenated_alignment.fa.ckp.gz   concatenated_alignment.fa.mldist    partitions.txt
	concatenated_alignment.fa.contree  concatenated_alignment.fa.model.gz




	cat concatenated_alignment.fa.treefile
	(WP_002214231.1:1.1961568151,(((WP_002217735.1:0.0028596312,(((WP_003692822.1:0.0709516245,(WP_003694978.1:0.0170596594,WP_047923833.1:0.0000010000)82:0.0367383018)92:0.0620938966,((((WP_010359920.1:0.0033742464,WP_105211411.1:0.0000010000)100:0.0231694635,WP_010951046.1:0.0348280149)91:0.0165037870,WP_047920986.1:0.0279616656)99:0.2179785135,(WP_047916984.1:0.0000010000,WP_048339587.1:0.0033651589)97:0.0276601858)48:0.0045916397)67:0.0298812248,WP_047949523.1:0.0037348774)26:0.0087497270)14:0.0000161576,((((((WP_002237787.1:0.1505153417,WP_047917697.1:0.0000010000)82:0.0401271430,(WP_047919925.1:0.0000010000,WP_124723901.1:0.0375818369)95:0.1098034759)39:0.0000010000,WP_017147127.1:0.0000010000)58:0.0068974802,((WP_003690584.1:0.0955813238,WP_003704890.1:0.0000010000)44:0.0000904707,WP_003690899.1:0.0056893878)9:0.0000080932)11:0.0061335142,WP_010359820.1:0.0000010000)4:0.0000010000,WP_003704449.1:0.0397705756)14:0.0091071979)47:0.0448822053,WP_047921846.1:0.0030523104)35:0.0440671907,(((((WP_002232153.1:0.0101621505,(((WP_003690898.1:0.0286620365,(((((WP_003692820.1:0.0596389020,WP_010951045.1:0.0768372327)88:0.1021553589,WP_082277597.1:0.0519982605)66:0.0120416357,(WP_003694976.1:0.0108385476,WP_047923834.1:0.0122833862)91:0.0386340819)25:0.0068779179,WP_047920984.1:0.0458129065)43:0.0208703766,(WP_012503479.1:0.0000010000,WP_124693319.1:0.0055448298)86:0.0310285603)33:0.0121740327)33:0.0157644632,(WP_010359918.1:0.0000010000,WP_082298766.1:0.0110857425)85:0.0295226131)7:0.0000482256,WP_082277595.1:0.0171666658)9:0.0000598183)8:0.0002020002,WP_047949522.1:0.0425939229)5:0.0000010000,WP_226888059.1:0.3852078773)12:0.0000025542,WP_012503480.1:0.0000010000)36:0.0104786721,WP_071198137.1:0.0000010000)37:0.0000022621);



---
посмотреть mafft как вариабельность там можно посмотреть 

---


- в гитхаб пайплайны
- по одинаковым данным для orthofinder и protein ortho 
- деревья 