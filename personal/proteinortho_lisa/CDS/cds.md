## STEP 1 - downloading (mannually, using hands)
filters:    
✔ Neisseria gonorrhoeae (taxon 485)    
✔ Annotated genomes only   
✔ GenBank assemblies (GCA_)   
✔ Chromosome and complete levels   
✔ Search for "WHO" in metadata   
✔ Include Genomic coding sequences (FASTA)   

Got data like this `/home/ktsel/nirma/cds/*/*.fna` Applying script "copy_and_rename.sh" to make them all in one folder with different naming `GCA_xxxxxx_cds_from_genomic.fna` like.

## STEP 2 - run proteinortho

```
proteinortho6.pl -project=cds_ortho \
    -cpus=8 \
    -p=blastn+ \
    -e=1e-10 \
    filtered/*.fna
``` 


*Parameter-info:* 
```
Parameter-vector : (version=6.3.4,step=0,verbose=1,debug=0,synteny=0,duplication=2,cs=3,alpha=0.5,connectivity=0.1,cpus=8,evalue=1e-10,purity=-1,coverage=50,identity=25,blastmode=blastn+,sim=0.95,report=3,keep=0,force=0,selfblast=0,twilight=0,core=0,coreMinSpecies=0,coreMaxProts=10,pseudo=1,omni=0,identical=0,range=-1,singles=0,clean=0,blastOptions=,makeBlastOptions=,nograph=0,xml=0,desc=0,tmp_path=./proteinortho_cache_cds_ortho/,blastversion=unknown,binpath=,makedb=,blast=,jobs_todo=0,project=cds_ortho,inproject=cds_ortho,po_path=/home/ktsel/miniconda3/envs/proteinortho/bin/,run_id=,threads_per_process=1,um=0)--------
```
Extract core proteome  
`proteinortho_grab_proteins.pl -core cds_ortho.proteinortho.tsv filtered/*.fna`

check that each core file contains one protein for each of the groups that span across all species:    
`for f in *core; do echo $f; grep ">" $f | wc -l; done`
>(1864)   

Meaning that we ve got 1864 core orthogroups finally.

## Step 3. Working with tsv and core files 

Parse the tsv file obtained after proteinortho.
Orthologous groups are selected where the number of species (# Species) and the number of genes (Genes) coincide with the number of samples used, and the parameter "Alg.-Conn." is not taken into account. For each such orthogroup, extract proteins from core files (using indices). Make a mafft alignment. Evaluate the variability of the orthogroup using the Shannon index and another simple variabily calculation in comparison.

### Shannon script

We will use python script for that - [align+shannon.py](<align+shannon.py>), there is also one extra parametr to check - write out logs for each alignment.

Firstly create this kind of structure:
```
.
├── align+shannon.py
├── cds/
│   ├── GCA_XXXX_cds_from_genomic.fna.core
│   ├── GCA_YYYY_cds_from_genomic.fna.core
│   └──  etc.
├── cds_ortho.proteinortho.tsv
```
(workig in mafft env)

`python3 align+shannon.py` (time 7:37 - 10:00)
> number of orthogroup for analysis  1769   
> mv *.fna cds/mafft_out/    
> mv mafft_logs/ cds/mafft_out/

As an output we got bunch (1769) of aligned fna files ([example](<132_aligned.fna>)), tsv file with shannon indexes and same bunch of mafft logs ([example](<132_mafft.log>)).  They provide information about the alignment process itself, but unfortunately they do not directly report the percentage of sequence identity or variability. We will now work with tsv file to get statistically significant data.

**Full Summary Stats (Mean, Median, Min, Max, Std Dev)**
```
awk '
NR > 1 {
    sum += $2; 
    vals[NR] = $2; 
    if ($2 > max || NR == 2) max = $2; 
    if ($2 < min || NR == 2) min = $2; 
} 
END { 
    mean = sum / (NR - 1); 
    asort(vals); 
}' shannon_indices.tsvstdev;  - 1)); [i] - mean)^2 } int(NR/2)+1]) / 2 : vals[int(NR/2)+1]; 
Mean: 0.0145891
Median: 0.003
Min: 0.0000
Max: 0.5414
Std Dev: 0.0458642  
```
Make another simple python [script](<shannon_filt.py>) - Reads data from the file shannon_indices.tsv. Calculates and displays general statistics. Applies a filter (Mean + 3 × Std Dev) to the Shannon index values. Saves the filtered data to the file [shannon_indices_filt.tsv](<shannon_indices_filt.tsv>):

```
Mean: 0.0145891
Median: 0.0030000
Min: 0.0000000
Max: 0.5414000
Std Dev: 0.0458772
Highly variable orthogroups selected: 36
Filtered groups saved to 'shannon_indices_filt.tsv' (threshold: 0.1522208)
```
Afterwords we recieve filtered groups with Shannon index as filtering parametr, you can brifely look at sorted data below:   
`head -n1 shannon_indices_filt.tsv && tail -n +2 shannon_indices_filt.tsv | sort -k2,2 -gr | column -t `

```
Orthogroup      ShannonIndex
1861  0.5414
1858  0.4783
909   0.4628
1859  0.4455
764   0.4243
1833  0.3959
1603  0.3811
1863  0.3707
1857  0.3638
126   0.3486
1854  0.3433
1419  0.3413
132   0.3387
185   0.3286
1034  0.328
1362  0.3166
1714  0.3023
498   0.288
1596  0.2594
1409  0.2579
225   0.2505
949   0.2326
511   0.2277
1855  0.2139
405   0.2082
1856  0.2075
812   0.2047
1715  0.2018
1351  0.1931
1161  0.1909
1333  0.1908
118   0.1894
1738  0.1779
1727  0.1575
776   0.1568
1350  0.1534
```

### Identity/variability script 

MAFFT itself does not directly output percent identity or variability in the standard output. However, there is a simple approach that allows you to obtain percent similarity or variability based on the aligned sequences after running MAFFT (as an alternative to Shannon index). The interpitation is simple as it is - A high percentage of variability (>50%) indicates high variability of the group.
Where is the difference tho? The Shannon index takes into account the frequency and diversity of all characters in the alignment columns and measures the "entropy". The percent identity (or variability) simply estimates the proportion of identical positions among all aligned sequences.

First steps in the script ([align+variab.py](<align+variab.py>)) are the same, the only part that differs:
```
    # Calculate identity and variability percentages

Initialize identical_positions counter to 0
For each position (column) in the alignment:
    If all characters in this position are identical:
        Increment identical_positions counter by 1

Calculate total_positions as the total number of columns in the alignment

Compute identity_percent as:
    (identical_positions / total_positions) * 100

Compute variability_percent as:
    100 - identity_percent

Store identity_percent and variability_percent for the current orthogroup

```

> `python3 align+variab.py` (time 10:40 - 2:00)

The output will be safed at cds/mafft_out_v

**Full Summary Stats (Mean, Median, Min, Max, Std Dev)**

```
awk '
BEGIN {
    FS="\t"; 
    print "Statistics for Variability(%):"
}
NR > 1 {  # Skip header
    sum += $3; 
    vals[NR] = $3; 
    if (NR == 2 || $3 < min) min = $3; 
    if (NR == 2 || $3 > max) max = $3; 
} 
END { 
    n = NR - 1; 
    mean = sum / n; 
    
    # Calculate Median
    asort(vals); 
    if (n % 2 == 1) median = vals[int(n/2)+2]; 
    else median = (vals[n/2+1] + vals[n/2]) / 2; 
}' identity_variability_cds.tsv 
Statistics for Variability(%):
Mean: 2.63%
Median: 0.60%
Min: 0.00%
Max: 75.57%
Std Dev: 7.05%
```

Make another simple python [script](<variab_filt.py>) - Reads data from the file identity_variability_cds.tsv. Calculates and displays general statistics. Applies a filter (Mean + 3 × Std Dev) to the values. Saves the filtered data to the file [variable_orthogroups_filt.tsv](<variable_orthogroups_filt.tsv>):

```
Mean Variability: 2.6289486%
Median Variability: 0.6000000%
Min Variability: 0.0000000%
Max Variability: 75.5700000%
Std Dev of Variability: 7.0556013%
Highly variable orthogroups selected: 42
Filtered groups saved to 'variable_orthogroups_filt.tsv' (threshold: 23.7957525%)

```
Afterwords we recieved filtered groups with manual calculation as filtering parametr, you can brifely look at sorted data below:    
`head -n1 variable_orthogroups_filt.tsv && tail -n +2 variable_orthogroups_filt.tsv | sort -k3,3 -gr | column -t `

```
Orthogroup      Identity(%)     Variability(%)
1858  24.43  75.57
1596  34.49  65.51
1857  36.54  63.46
1034  38.03  61.97
909   39.16  60.84
1861  41.29  58.71
1863  42.86  57.14
126   47.86  52.14
1856  48.98  51.02
1419  49.02  50.98
1859  49.53  50.47
185   50.74  49.26
1805  51.2   48.8
1854  53.33  46.67
132   54.76  45.24
498   55.15  44.85
764   55.49  44.51
1833  57.76  42.24
776   58.4   41.6
1855  59.09  40.91
1362  60.0   40.0
1598  60.3   39.7
1603  61.09  38.91
511   63.01  36.99
1714  63.81  36.19
1599  63.94  36.06
873   64.37  35.63
1409  66.36  33.64
209   67.97  32.03
1194  67.98  32.02
1586  69.48  30.52
1715  69.5   30.5
849   69.85  30.15
812   71.97  28.03
969   72.25  27.75
1547  73.1   26.9
1061  73.32  26.68
622   73.78  26.22
225   74.0   26.0
1351  74.22  25.78
453   75.11  24.89
949   75.8   24.2

```

## Step 4. Aligning and tree-ing

### 4.1 Running [script](<tree_top5.py>) for top5 highly variable orthogroups (less then 10 mins)

The following will be created in the ./trees_5_shannon/ folder:

Files of the concatenated_{N}.fna type — concatenated sequences for N orthogroups.

Files of the aligned_{N}.fna type — aligned sequences.

IQ-TREE tree files (.treefile, .log, .contree, etc.).

output for snannon tree:
```
(iqtree) ktsel@LPT-pc1:~/nirma/cds$ iqtree -s trees_5_shannon/aligned_1.fna -m GTR+G -bb 1000 -nt AUTO
IQ-TREE multicore version 2.4.0 for Linux x86 64-bit built Feb 12 2025
Developed by Bui Quang Minh, Nguyen Lam Tung, Olga Chernomor, Heiko Schmidt,
Dominik Schrempf, Michael Woodhams, Ly Trong Nhan, Thomas Wong

Host:    LPT-pc1 (AVX2, FMA3, 7 GB RAM)
Command: iqtree -s trees_5_shannon/aligned_1.fna -m GTR+G -bb 1000 -nt AUTO
Seed:    37255 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Fri Apr 25 09:14:40 2025
Kernel:  AVX+FMA - auto-detect threads (20 CPU cores detected)

Reading alignment file trees_5_shannon/aligned_1.fna ... Fasta format detected
Reading fasta file: done in 0.00466819 secs
Alignment most likely contains DNA/RNA sequences
Alignment has 29 sequences with 7333 columns, 158 distinct patterns
154 parsimony-informative, 27 singleton sites, 7152 constant sites
                                     Gap/Ambiguity  Composition  p-value
Analyzing sequences: done in 0.00669291 secs using 1356% CPU
   1  lcl|CP145015.1_cds_V5G07_02935_567    57.66%    failed      1.39%
   2  lcl|CP145017.1_cds_V5G14_03135_603    57.66%    failed      1.39%
   3  lcl|CP145019.1_cds_V5G05_02930_566    57.66%    failed      1.39%
   4  lcl|CP145021.1_cds_V5G12_02920_565    57.66%    failed      1.39%
   5  lcl|CP145024.1_cds_V5G17_06195_1213   20.96%    passed     73.33%
   6  lcl|CP145026.1_cds_XCB86095.1_565      0.42%    passed     18.94%
   7  lcl|CP145028.1_cds_V5G19_05795_1139   20.63%    passed     83.92%
   8  lcl|CP145032.1_cds_V5F95_03520_682    57.66%    failed      1.14%
   9  lcl|CP145035.1_cds_XCB90816.1_1147     0.38%    passed     17.97%
  10  lcl|CP145037.1_cds_XCB94074.1_502      0.42%    passed     26.87%
  11  lcl|CP145041.1_cds_V5G13_02925_566    57.66%    failed      1.39%
  12  lcl|CP145045.1_cds_V5G18_06615_1289   20.63%    passed     82.28%
  13  lcl|CP145048.1_cds_V5G06_05825_1138   20.96%    passed     73.33%
  14  lcl|CP145050.1_cds_V5F92_03280_631    57.66%    failed      1.39%
  15  lcl|CP145052.1_cds_V5F98_06040_1188   20.63%    passed     69.70%
  16  lcl|CP145053.1_cds_XCC06720.1_1077     0.42%    passed     26.87%
  17  lcl|CP145057.1_cds_V5F97_06170_1208   20.96%    passed     77.35%
  18  lcl|CP145059.1_cds_V5F96_06210_1216   20.96%    passed     73.33%
  19  lcl|CP145062.1_cds_XCC12972.1_1148     0.26%    passed     24.97%
  20  lcl|CP145064.1_cds_V5F91_05705_1120    0.33%    passed     22.32%
  21  lcl|CP145067.1_cds_V5G04_02945_570    57.66%    failed      1.04%
  22  lcl|CP145071.1_cds_V5G11_05850_1145   20.96%    passed     74.59%
  23  lcl|CP145074.1_cds_V5F94_02570_501    57.66%    failed      0.46%
  24  lcl|CP145077.1_cds_V5G03_02920_565    57.66%    failed      0.94%
  25  lcl|CP145080.1_cds_V5G00_05785_1137   20.63%    passed     79.75%
  26  lcl|CP145083.1_cds_XCC26611.1_501      0.42%    passed     26.87%
  27  lcl|CP145086.1_cds_V5G10_05480_1076   20.63%    passed     83.92%
  28  lcl|CP145088.1_cds_XCC29306.1_1139     0.34%    passed     12.08%
  29  lcl|CP145090.1_cds_V5G15_02575_502    57.66%    failed      1.91%
WARNING: 11 sequences contain more than 50% gaps/ambiguity
****  TOTAL                                 29.14%  11 sequences failed composition chi2 test (p-value<5%; df=3)
NOTE: lcl|CP145017.1_cds_V5G14_03135_603 is identical to lcl|CP145015.1_cds_V5G07_02935_567 but kept for subsequent analysis
NOTE: lcl|CP145048.1_cds_V5G06_05825_1138 is identical to lcl|CP145024.1_cds_V5G17_06195_1213 but kept for subsequent analysis
NOTE: lcl|CP145086.1_cds_V5G10_05480_1076 is identical to lcl|CP145028.1_cds_V5G19_05795_1139 but kept for subsequent analysis
NOTE: lcl|CP145053.1_cds_XCC06720.1_1077 is identical to lcl|CP145037.1_cds_XCB94074.1_502 but kept for subsequent analysis
NOTE: 6 identical sequences (see below) will be ignored for subsequent analysis
NOTE: lcl|CP145019.1_cds_V5G05_02930_566 (identical to lcl|CP145015.1_cds_V5G07_02935_567) is ignored but added at the end
NOTE: lcl|CP145021.1_cds_V5G12_02920_565 (identical to lcl|CP145015.1_cds_V5G07_02935_567) is ignored but added at the end
NOTE: lcl|CP145041.1_cds_V5G13_02925_566 (identical to lcl|CP145015.1_cds_V5G07_02935_567) is ignored but added at the end
NOTE: lcl|CP145050.1_cds_V5F92_03280_631 (identical to lcl|CP145015.1_cds_V5G07_02935_567) is ignored but added at the end
NOTE: lcl|CP145059.1_cds_V5F96_06210_1216 (identical to lcl|CP145024.1_cds_V5G17_06195_1213) is ignored but added at the end
NOTE: lcl|CP145083.1_cds_XCC26611.1_501 (identical to lcl|CP145037.1_cds_XCB94074.1_502) is ignored but added at the end
Alignment was printed to trees_5_shannon/aligned_1.fna.uniqueseq.phy

For your convenience alignment with unique sequences printed to trees_5_shannon/aligned_1.fna.uniqueseq.phy

Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.007 seconds
Generating 1000 samples for ultrafast bootstrap (seed: 37255)...

NOTE: 1 MB RAM (0 GB) is required!
Measuring multi-threading efficiency up to 20 CPU cores
Increase to 10 rounds for branch lengths
4047 trees examined
Threads: 1 / Time: 20.004 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -13985
Threads: 2 / Time: 17.207 sec / Speedup: 1.163 / Efficiency: 58% / LogL: -13985
Threads: 3 / Time: 18.019 sec / Speedup: 1.110 / Efficiency: 37% / LogL: -13985
BEST NUMBER OF THREADS: 2

Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -11876.702
2. Current log-likelihood: -11725.664
3. Current log-likelihood: -11725.227
Optimal log-likelihood: -11725.221
Rate parameters:  A-C: 1.26084  A-G: 3.28282  A-T: 0.82953  C-G: 1.48923  C-T: 6.05228  G-T: 1.00000
Base frequencies:  A: 0.191  C: 0.232  G: 0.383  T: 0.194
Gamma shape alpha: 0.020
Parameters optimization took 3 rounds (0.030 sec)
Wrote distance file to... 
Computing ML distances based on estimated model parameters...
Calculating distance matrix: done in 0.00215569 secs using 261.5% CPU
Computing ML distances took 0.002288 sec (of wall-clock time) 0.005775 sec (of CPU time)
Setting up auxiliary I and S matrices: done in 0.000337722 secs using 127% CPU
Computing RapidNJ tree took 0.006759 sec (of wall-clock time) 0.005338 sec (of CPU time)
Log-likelihood of RapidNJ tree: -11817.102
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 0.081 second
Computing log-likelihood of 98 initial trees ... 0.098 seconds
Current best score: -11724.614

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -11722.799
BETTER TREE FOUND at iteration 3: -11722.798
BETTER TREE FOUND at iteration 10: -11722.797
Iteration 10 / LogL: -11722.797 / Time: 0h:0m:55s
Iteration 20 / LogL: -11722.799 / Time: 0h:0m:55s
Finish initializing candidate tree set (15)
Current best tree score: -11722.797 / CPU time: 0.499
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
UPDATE BEST LOG-LIKELIHOOD: -11722.795
Iteration 30 / LogL: -11722.812 / Time: 0h:0m:56s (0h:2m:34s left)
Iteration 40 / LogL: -11722.807 / Time: 0h:0m:56s (0h:1m:40s left)
Iteration 50 / LogL: -11722.810 / Time: 0h:0m:56s (0h:1m:9s left)
Log-likelihood cutoff on original alignment: -11771.439
Iteration 60 / LogL: -11722.807 / Time: 0h:0m:56s (0h:0m:48s left)
Iteration 70 / LogL: -11722.820 / Time: 0h:0m:56s (0h:0m:33s left)
Iteration 80 / LogL: -11722.830 / Time: 0h:0m:57s (0h:0m:21s left)
Iteration 90 / LogL: -11722.837 / Time: 0h:0m:57s (0h:0m:12s left)
UPDATE BEST LOG-LIKELIHOOD: -11722.795
Iteration 100 / LogL: -11722.830 / Time: 0h:0m:57s (0h:0m:5s left)
Log-likelihood cutoff on original alignment: -11772.264
NOTE: Bootstrap correlation coefficient of split occurrence frequencies: 0.994
Iteration 110 / LogL: -11725.338 / Time: 0h:0m:57s (0h:0m:47s left)
TREE SEARCH COMPLETED AFTER 111 ITERATIONS / Time: 0h:0m:57s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -11722.795
Optimal log-likelihood: -11722.792
Rate parameters:  A-C: 1.30829  A-G: 3.36371  A-T: 0.76788  C-G: 1.65099  C-T: 6.24190  G-T: 1.00000
Base frequencies:  A: 0.191  C: 0.232  G: 0.383  T: 0.194
Gamma shape alpha: 0.020
Parameters optimization took 1 rounds (0.006 sec)
BEST SCORE FOUND : -11722.792
Creating bootstrap support values...
Split supports printed to NEXUS file trees_5_shannon/aligned_1.fna.splits.nex
Total tree length: 0.045

Total number of iterations: 111
CPU time used for tree search: 4.889 sec (0h:0m:4s)
Wall-clock time used for tree search: 2.587 sec (0h:0m:2s)
Total CPU time used: 114.002 sec (0h:1m:54s)
Total wall-clock time used: 58.015 sec (0h:0m:58s)

Computing bootstrap consensus tree...
Reading input file trees_5_shannon/aligned_1.fna.splits.nex...
23 taxa and 171 splits.
Consensus tree written to trees_5_shannon/aligned_1.fna.contree
Reading input trees file trees_5_shannon/aligned_1.fna.contree
Log-likelihood of consensus tree: -11722.883

Analysis results written to: 
  IQ-TREE report:                trees_5_shannon/aligned_1.fna.iqtree
  Maximum-likelihood tree:       trees_5_shannon/aligned_1.fna.treefile
  Likelihood distances:          trees_5_shannon/aligned_1.fna.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          trees_5_shannon/aligned_1.fna.splits.nex
  Consensus tree:                trees_5_shannon/aligned_1.fna.contree
  Screen log file:               trees_5_shannon/aligned_1.fna.log

Date and Time: Fri Apr 25 09:15:38 2025

```
resultant tree for 1st highly variable orthogroup (according shannon index) - [aligned_1.fna.treefile](<aligned_1.fna.treefile>) 

### 4.2 Gene ID → GCA genome accession → Strain name
This file contain gene ID in trees, we would like to see, what are the strains that hided under these IDs. 
[Script](<get_strain.py>) parses the Newick tree from IQ-TREE output; matches gene IDs (lcl|...) to sample IDs (GCA_* ) via cds_ortho.tsv.; maps sample IDs (GCA_* ) to strain names via data_summary.tsv.; creates a new Newick-format tree with strains instead of genes: 


Figure represents two trees (on the left handside - 1st highly variable orthogroup (aligned_1.fna.treefile) and on the right handside - first 5 highly variable ones transfered to strains)

![alt text](<Shannon_1_5_strain.png>)

tree-strains-shannon files are in folder [modified_trees](modified_trees/)

Repeating step 4.1 (run 7 mins) and 4.2 (run several seconds) for variability via simle calculation.

tree-strains-variab files are in folder [modified_trees_variab](modified_trees_variab/)


#################################
```
Orthogroup      ShannonIndex    CDS
1861  0.5414    lcl|CP145015.1_cds_V5G07_02935_567 [locus_tag=V5G07_02935] [protein=PLxRFG domain-containing protein] [pseudo=true] [location=550831..553935] [gbkey=CDS]
1858  0.4783    lcl|CP145015.1_cds_XCB76284.1_977 [locus_tag=V5G07_05020] [protein=baseplate J/gp47 family protein] [protein_id=XCB76284.1] [location=941730..942785] [gbkey=CDS]
909   0.4628    lcl|CP145015.1_cds_XCB75446.1_2238 [locus_tag=V5G07_11550] [protein=hypothetical protein] [protein_id=XCB75446.1] [location=2188079..2188225] [gbkey=CDS]
1859  0.4455    lcl|CP145015.1_cds_XCB76563.1_2175 [locus_tag=V5G07_11235] [protein=glycosyltransferase family 4 protein] [protein_id=XCB76563.1] [location=complement(2125622..2126320)] [gbkey=CDS]
764   0.4243    lcl|CP145015.1_cds_XCB75290.1_2049 [locus_tag=V5G07_10560] [protein=hypothetical protein] [protein_id=XCB75290.1] [location=2006082..2006435] [gbkey=CDS]
```

```
Orthogroup     Variability(%)   CDS
1858  75.57 lcl|CP145015.1_cds_XCB76284.1_977 [locus_tag=V5G07_05020] [protein=baseplate J/gp47 family protein] [protein_id=XCB76284.1] [location=941730..942785] [gbkey=CDS]     
1596  65.51 lcl|CP145015.1_cds_XCB76297.1_1002 [locus_tag=V5G07_05145] [protein=hypothetical protein] [protein_id=XCB76297.1] [location=complement(961954..962322)] [gbkey=CDS]     
1857  63.46 lcl|CP145015.1_cds_XCB75055.1_1758 [locus_tag=V5G07_09020] [protein=hypothetical protein] [protein_id=XCB75055.1] [location=complement(1714244..1714906)] [gbkey=CDS]
1034  61.97 lcl|CP145015.1_cds_XCB75585.1_122 [locus_tag=V5G07_00620] [protein=hypothetical protein] [protein_id=XCB75585.1] [location=134186..134458] [gbkey=CDS]
909  60.84 lcl|CP145015.1_cds_XCB75446.1_2238 [locus_tag=V5G07_11550] [protein=hypothetical protein] [protein_id=XCB75446.1] [location=2188079..2188225] [gbkey=CDS] 
```