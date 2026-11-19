```python
print("Day 1 - Dengue Genomic Analysis Project")
print("Starting bioinformatics pipeline")
```

    Day 1 - Dengue Genomic Analysis Project
    Starting bioinformatics pipeline



```python
import subprocess
subprocess.run(["pip", "install", "biopython"])

```




    CompletedProcess(args=['pip', 'install', 'biopython'], returncode=0)




```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt
print("All libraries loaded successfully!")

```

    All libraries loaded successfully!


# Dengue Virus Genomic Analysis

## Background

Dengue fever is a mosquito-borne viral disease caused by the dengue 
virus (DENV), transmitted primarily by Aedes aegypti mosquitoes. 
It is a major public health concern globally, with an estimated 
400 million infections occurring annually.

The dengue virus exists as four genetically distinct serotypes: 
DENV-1, DENV-2, DENV-3, and DENV-4. While all four serotypes 
cause dengue fever, they are sufficiently different at the genomic 
level that infection with one serotype does not confer immunity 
against the others. This phenomenon — known as the lack of 
cross-protective immunity — means that individuals in endemic 
regions like Tamil Nadu can be infected up to four times in 
their lifetime.

A second infection by a different serotype carries a significantly 
higher risk of severe dengue (dengue haemorrhagic fever), which 
can be fatal. This makes understanding the genomic differences 
between serotypes a critical research question.

## Research Question
How genomically distinct are the four dengue virus serotypes, 
and can these differences be visualised through phylogenetic 
analysis and GC content comparison?

## Step 1 - Loading and Exploring the Sequences



```python
from Bio import SeqIO
import os

data_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data"

print("=== DENV-1 Sequences ===")
for record in SeqIO.parse(os.path.join(data_path, "DENV1.fasta"), "fasta"):
    print(f"ID: {record.id} | Length: {len(record.seq)} bp")
```

    === DENV-1 Sequences ===
    ID: MW315195.1 | Length: 10734 bp
    ID: PX297362.1 | Length: 10734 bp
    ID: MW362474.1 | Length: 10735 bp
    ID: KR052012.1 | Length: 10714 bp
    ID: AF226685.2 | Length: 10735 bp



```python
import os
print(os.path.expanduser("~"))
```

    C:\Users\Gaayathri PS



```python
serotypes = ["DENV1", "DENV2", "DENV3", "DENV4"]

for serotype in serotypes:
    print(f"\n=== {serotype} Sequences ===")
    for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
        print(f"ID: {record.id} | Length: {len(record.seq)} bp")
```

    
    === DENV1 Sequences ===
    ID: MW315195.1 | Length: 10734 bp
    ID: PX297362.1 | Length: 10734 bp
    ID: MW362474.1 | Length: 10735 bp
    ID: KR052012.1 | Length: 10714 bp
    ID: AF226685.2 | Length: 10735 bp
    
    === DENV2 Sequences ===
    ID: PV752130.1 | Length: 10723 bp
    ID: PQ657766.1 | Length: 10723 bp
    ID: MH822944.2 | Length: 10723 bp
    ID: MZ008467.1 | Length: 10724 bp
    ID: KP406804.1 | Length: 10712 bp
    
    === DENV3 Sequences ===
    ID: PX297385.1 | Length: 10709 bp
    ID: PV918782.1 | Length: 10707 bp
    ID: MH051731.1 | Length: 10766 bp
    ID: MH891766.1 | Length: 10764 bp
    ID: JX669494.1 | Length: 10710 bp
    
    === DENV4 Sequences ===
    ID: MG601754.1 | Length: 10772 bp
    ID: KF041260.1 | Length: 10652 bp
    ID: JQ922560.1 | Length: 10739 bp
    ID: JN983813.1 | Length: 10649 bp
    ID: AY618993.1 | Length: 10653 bp



```python
import pandas as pd

serotypes = ["DENV1", "DENV2", "DENV3", "DENV4"]

data = []

for serotype in serotypes:
    for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
        data.append({
            "Accession": record.id,
            "Serotype": serotype,
            "Length": len(record.seq)
        })

df = pd.DataFrame(data)
print(df)
```

         Accession Serotype  Length
    0   MW315195.1    DENV1   10734
    1   PX297362.1    DENV1   10734
    2   MW362474.1    DENV1   10735
    3   KR052012.1    DENV1   10714
    4   AF226685.2    DENV1   10735
    5   PV752130.1    DENV2   10723
    6   PQ657766.1    DENV2   10723
    7   MH822944.2    DENV2   10723
    8   MZ008467.1    DENV2   10724
    9   KP406804.1    DENV2   10712
    10  PX297385.1    DENV3   10709
    11  PV918782.1    DENV3   10707
    12  MH051731.1    DENV3   10766
    13  MH891766.1    DENV3   10764
    14  JX669494.1    DENV3   10710
    15  MG601754.1    DENV4   10772
    16  KF041260.1    DENV4   10652
    17  JQ922560.1    DENV4   10739
    18  JN983813.1    DENV4   10649
    19  AY618993.1    DENV4   10653



```python
from Bio.SeqUtils import gc_fraction

data = []

for serotype in serotypes:
    for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
        gc = round(gc_fraction(record.seq) * 100, 2)
        data.append({
            "Accession": record.id,
            "Serotype": serotype,
            "Length": len(record.seq),
            "GC Content (%)": gc
        })

df = pd.DataFrame(data)
print(df)
```

         Accession Serotype  Length  GC Content (%)
    0   MW315195.1    DENV1   10734           46.28
    1   PX297362.1    DENV1   10734           46.62
    2   MW362474.1    DENV1   10735           46.62
    3   KR052012.1    DENV1   10714           46.76
    4   AF226685.2    DENV1   10735           46.44
    5   PV752130.1    DENV2   10723           45.73
    6   PQ657766.1    DENV2   10723           45.77
    7   MH822944.2    DENV2   10723           45.83
    8   MZ008467.1    DENV2   10724           45.63
    9   KP406804.1    DENV2   10712           46.07
    10  PX297385.1    DENV3   10709           46.61
    11  PV918782.1    DENV3   10707           46.74
    12  MH051731.1    DENV3   10766           46.38
    13  MH891766.1    DENV3   10764           46.66
    14  JX669494.1    DENV3   10710           46.50
    15  MG601754.1    DENV4   10772           47.34
    16  KF041260.1    DENV4   10652           47.44
    17  JQ922560.1    DENV4   10739           47.56
    18  JN983813.1    DENV4   10649           47.19
    19  AY618993.1    DENV4   10653           47.17



```python
import matplotlib.pyplot as plt

# Calculate average GC per serotype
avg_gc = df.groupby("Serotype")["GC Content (%)"].mean().round(2)

# Plot
plt.figure(figsize=(8, 5))
bars = plt.bar(avg_gc.index, avg_gc.values, 
               color=["#4C9BE8", "#E8734C", "#4CE87B", "#E8C84C"],
               edgecolor="black", width=0.5)

# Add value labels on top of each bar
for bar, val in zip(bars, avg_gc.values):
    plt.text(bar.get_x() + bar.get_width()/2, 
             bar.get_height() + 0.05, 
             f"{val}%", ha="center", fontsize=11)

plt.title("Average GC Content by Dengue Virus Serotype", fontsize=13)
plt.xlabel("Serotype", fontsize=11)
plt.ylabel("GC Content (%)", fontsize=11)
plt.ylim(44, 48.5)
plt.tight_layout()
plt.savefig(r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\images\GC_content.png", dpi=300)
plt.show()

print("Chart saved!")
```


    
![png](Dengue_analysis_final_files/Dengue_analysis_final_10_0.png)
    


    Chart saved!



```python
import os

data_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data"
output_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data\all_dengue.fasta"

serotypes = ["DENV1", "DENV2", "DENV3", "DENV4"]

with open(output_path, "w") as outfile:
    for serotype in serotypes:
        for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
            outfile.write(f">{serotype}_{record.id}\n{str(record.seq)}\n")

print("Combined file created!")
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[1], line 10
          8 with open(output_path, "w") as outfile:
          9     for serotype in serotypes:
    ---> 10         for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
         11             outfile.write(f">{serotype}_{record.id}\n{str(record.seq)}\n")
         13 print("Combined file created!")


    NameError: name 'SeqIO' is not defined



```python
from Bio import SeqIO
import os

data_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data"
```


```python
import os

data_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data"
output_path = r"C:\Users\Gaayathri PS\OneDrive\Documents\DENGUE_PROJECT\Data\all_dengue.fasta"

serotypes = ["DENV1", "DENV2", "DENV3", "DENV4"]

with open(output_path, "w") as outfile:
    for serotype in serotypes:
        for record in SeqIO.parse(os.path.join(data_path, f"{serotype}.fasta"), "fasta"):
            outfile.write(f">{serotype}_{record.id}\n{str(record.seq)}\n")

print("Combined file created!")
```

    Combined file created!


## Methods

### Data collection
Complete genome sequences for all four dengue virus serotypes were retrieved from NCBI GenBank. Five genomes per serotype (20 in total, each about 10,650-10,770 bp) were selected to cover different countries and years, including isolates from South Asia, Southeast Asia, East Asia, the Americas and the Pacific.

### Analysis
- Sequences were read and GC content calculated in Python using BioPython, with results organised in pandas and plotted with matplotlib.
- The 20 genomes were combined into one FASTA file and aligned with MUSCLE in MEGA 12.
- A Neighbour-Joining tree (Maximum Composite Likelihood distances, pairwise deletion, 1000 bootstrap replicates) was built in MEGA 12.

## Results

### GC Content Analysis
GC content was calculated for all 20 dengue genome sequences (5 per serotype). Average values per serotype were:

| Serotype | Average GC Content (%) |
|----------|------------------------|
| DENV-1   | 46.54                  |
| DENV-2   | 45.81                  |
| DENV-3   | 46.58                  |
| DENV-4   | 47.34                  |

DENV-4 had the highest average GC content (47.34%) and DENV-2 the lowest (45.81%), a difference of about 1.5 percentage points. Values were similar within each serotype, and every DENV-4 genome was higher than every DENV-2 genome. With only 5 genomes per serotype and no statistical test, this is an observed trend rather than a tested difference.

![GC Content Chart](images/GC_content.png)

### Phylogenetic Analysis
A Neighbour-Joining tree was built from 20 aligned complete genomes in MEGA 12 (1000 bootstrap replicates). The sequences separated into four clusters matching the four serotypes, each supported by 100% bootstrap values, regardless of country or year of isolation. DENV-1 and DENV-3 grouped most closely, DENV-2 joined them next, and DENV-4 was the most divergent in this analysis. Lower bootstrap values (69-99%) occurred only within serotypes, indicating less certain relationships among individual strains.

The tree is unrooted and based on 5 genomes per serotype, so conclusions are limited to serotype-level separation.

![Phylogenetic Tree](images/phylogenetic_tree.png)

## Discussion

The results show that the four dengue serotypes are distinct at the genome level. In this dataset they differ slightly in average GC content, and in the phylogenetic tree every genome grouped with its own serotype, with 100% bootstrap support for each of the four clusters.

Genomes from different countries and years (for example India, Brazil, Thailand and China within DENV-1) clustered by serotype rather than by place. This suggests that serotype is the main source of genetic difference in this dataset, larger than geographic variation within a serotype.

The clear genetic separation between serotypes is consistent with the well-documented clinical observation that infection with one serotype does not give lasting protection against the others. Because the serotypes differ in the proteins the immune system recognises, such as the Envelope (E) protein, antibodies made against one serotype may not neutralise another. This analysis looked at whole genomes and did not test the E protein or immune response directly, so it supports this explanation but does not prove it.

These findings are relevant to Tamil Nadu, India, where dengue is endemic and more than one serotype can circulate. They also fit the difficulty of building a single vaccine that protects against all four serotypes.

### Limitations
- Only 5 genomes per serotype were analysed, chosen for geographic and temporal variety rather than at random.
- The GC content comparison was not statistically tested.
- Neighbour-Joining is a distance-based method and the tree is unrooted, so no serotype can be called "ancestral". A maximum likelihood tree would be a stronger follow-up.
- Sequences from Tamil Nadu specifically were not included.

## Conclusion

This project compared 20 complete dengue virus genomes (5 per serotype, DENV-1 to DENV-4) from NCBI GenBank using BioPython for GC content analysis and MEGA 12 for alignment and phylogenetic reconstruction. Two findings emerged:

1. Average GC content was highest in DENV-4 (47.34%) and lowest in DENV-2 (45.81%), with similar values within each serotype.

2. A Neighbour-Joining tree (1000 bootstraps) separated the genomes into four clusters matching the four serotypes, each with 100% bootstrap support. DENV-1 and DENV-3 grouped most closely, and DENV-4 was the most divergent.

Together these results show clear genomic separation between serotypes, consistent with the lack of cross-protective immunity that makes dengue a continuing public health problem in endemic regions such as South and Southeast Asia. The small sample size and unrooted tree limit how far these conclusions can be taken.

## Tools and Data Sources
- NCBI GenBank: sequence data
- BioPython: reading sequences and GC content analysis
- MEGA 12: multiple sequence alignment (MUSCLE) and Neighbour-Joining tree construction
- Python 3 (pandas, matplotlib): data organisation and visualisation
- Jupyter Notebook: analysis environment
- GitHub: project repository
