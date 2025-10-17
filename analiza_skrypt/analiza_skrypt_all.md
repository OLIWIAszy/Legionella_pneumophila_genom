# Legionella_pneumophila_genom
Analiza genomu bakterii *Legionella pneumophila*, po sekwencjonowaniu wysokoprzepustowym.


## Sekwencjonowanie NGS na platformie Illumina 

### Dane wejściowe 

```shell

# Dane wejściowe
conda activate bioinfo 

BACTERIA_FOLDER="/Legionella_pneumophila/short_reads" 
SAMPLE_NAME='L7_S252'


# Definiowanie ścieżek do folderów
RAW_READS_DIR="${BACTERIA_FOLDER}/raw_reads"
FASTP_DIR="${BACTERIA_FOLDER}/fastp"
UNICYCLER_DIR="${BACTERIA_FOLDER}/unicycler/unicycler_${SAMPLE_NAME}"
REFERENCE_DIR="${BACTERIA_FOLDER}/ref_genome_legionella"
REF_FILE="${REFERENCE_DIR}/GCF_001941585.1_ASM194158v1_genomic.fna"
RAGTAG_DIR="${BACTERIA_FOLDER}/ragtag/ragtag_${SAMPLE_NAME}"
PROKKA_DIR="${BACTERIA_FOLDER}/prokka/prokka_${SAMPLE_NAME}"
QUAST_DIR="${BACTERIA_FOLDER}/quast/quast_${SAMPLE_NAME}"
CHECKM_DIR="${BACTERIA_FOLDER}/checkm/checkm_${SAMPLE_NAME}"





```

### Fastp analiza jakości odczytów

```shell

mkdir -p "${FASTP_DIR}"

fastp -i "${RAW_READS_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz" \
      -I "${RAW_READS_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz" \
      -o "${FASTP_DIR}/fastp_${SAMPLE_NAME}_R1_001.fastq.gz" \
      -O "${FASTP_DIR}/fastp_${SAMPLE_NAME}_R2_001.fastq.gz" \
      --json "${FASTP_DIR}/${SAMPLE_NAME}_report.json" \
      --html "${FASTP_DIR}/${SAMPLE_NAME}_report.html"

```
### Unicycler składanie genomu

```shell

unicycler -1 "${FASTP_DIR}/fastp_${SAMPLE_NAME}_R1_001.fastq.gz" \
          -2 "${FASTP_DIR}/fastp_${SAMPLE_NAME}_R2_001.fastq.gz" \
          -o "${UNICYCLER_DIR}" \
          --min_fasta_length 1000

```

### Ragtag ustawianie kontigów

```shell

mkdir -p "${RAGTAG_DIR}"

ragtag.py scaffold \
          "${REF_FILE}"\
          "${UNICYCLER_DIR}/assembly.fasta" \
          -o ragtag

```

### Prokka adnotacja genomu

```shell

mkdir -p "${PROKKA_DIR}"

prokka --outdir "${PROKKA_DIR}" --force \
       --prefix "${SAMPLE_NAME}" --addgenes \
       --kingdom "Bacteria" \
       --genus "Legionella" \
       --species "pneumophila" \
       --cpus 4 \
       "${RAGTAG_DIR}/ragtag.scaffold.fasta"

```

### Quast analiza jakości złożenia genomu

```shell

mkdir -p "${QUAST_DIR}"

quast.py -o "${QUAST_DIR}" \
        "${RAGTAG_DIR}/ragtag.scaffold.fasta" \
        -r "${REF_FILE}" \
        -1 "${RAW_READS_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz" \
        -2 "${RAW_READS_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz" \
        -g "${PROKKA_DIR}/${SAMPLE_NAME}.gff"

```

### CheckM analiza jakości adnotacji genomu

```shell

mkdir -p "${CHECKM_DIR}"

checkm taxonomy_wf genus Legionella \
                   --genes \
                   -t 10 \
                   -x faa \
                   --file ${BACTERIA_FOLDER}/checkm_stats_${SAMPLE_NAME}.txt \
                   "${PROKKA_DIR}" \
                   "${CHECKM_DIR}"
                  
```

### Eggnog-mapper funkcjonalna adnotacja sekwencji genomowych
#### http://eggnog-mapper.embl.de/ i użycie pliku z Prokka .faa

### Łączenie plików eggmapper i prokka files w jedną adnotację genomu
```python
import pandas as pd

simple_name = "R32_S257"

# Ścieżki dostępu do plików
prokka_path = f"/Legionella_pneumophila/short_reads/prokka/prokka_{simple_name}/{simple_name}.tsv"
eggmaper_path = f"/Legionella_pneumophila/short_reads/eggmaper/{simple_name}/eggmaper_{simple_name}.tsv"
merged_path = f"/Legionella_pneumophila/short_reads/prokkaxeggmaper/merge_{simple_name}.tsv"
filtered_path = f"/Legionella_pneumophila/short_reads/prokkaxeggmaper/allright_{simple_name}.tsv"

# Wczytywanie plików
prokkatsv = pd.read_csv(prokka_path, sep='\t')
eggmapertsv = pd.read_csv(eggmaper_path, sep='\t')

# Łączenie
Output_df2 = pd.merge(prokkatsv, eggmapertsv, on='locus_tag', how='outer')
Output_df2.to_csv(merged_path, sep="\t", header=True, index=False)

# Filtrowanie
outputgenes = pd.read_csv(merged_path, sep='\t')
allright = outputgenes[outputgenes['ftype'] != 'gene']
allright.to_csv(filtered_path, sep="\t", header=True, index=False)


```


## Sekwencjonowanie trzeciej generacji oparte na nanoporach

### Dane wejściowe
```shell

# Dane wejściowe
conda activate bioinfo 

BACTERIA_FOLDER="/Legionella_pneumophila/long_reads" 
SAMPLE_NAME='L40'
NAME='legionella'

# Definiowanie ścieżek do folderów
RAW_READS_DIR="${BACTERIA_FOLDER}/raw_reads"

FLYE_DIR="${BACTERIA_FOLDER}/flye/flye_${SAMPLE_NAME}"
CIRCLATOR_DIR="${BACTERIA_FOLDER}/circlator/circlator_${SAMPLE_NAME}"
BAKTA_DIR="${BACTERIA_FOLDER}/bakta/bakta_${SAMPLE_NAME}.40"
QUAST_DIR="${BACTERIA_FOLDER}/quast/quast_${SAMPLE_NAME}"
CHECKM2_DIR="${BACTERIA_FOLDER}/checkm2/checkm2_${SAMPLE_NAME}"
DNAA_DIR="${BACTERIA_FOLDER}/dnaa_seq"
BAKTA_DB="/db"
CHECKM_INPUT_DIR="${BACTERIA_FOLDER}/checkm2/imput_checkm2/checkm_input_${SAMPLE_NAME}"
REF_FILE="${BACTERIA_FOLDER}/ref_genome_legionella/GCF_001941585.1_ASM194158v1_genomic.fna"

```

### Flye składanie sekwencji nukleotydowych długich pojedyńczych cząsteczek

```shell

mkdir -p "${FLYE_DIR}"

flye --nano-raw "${RAW_READS_DIR}/${SAMPLE_NAME}.fastq.gz" \
     --out-dir "${FLYE_DIR}" \
     --genome-size 4m \
     --asm-coverage 70 \
     --threads 4

```

### Circlator cyrkulacja genomu

```shell

mkdir -p "${BACTERIA_FOLDER}/circlator"

circlator all \
    --merge_min_id 85 \
    --threads 4 \
    --verbose \
    --data_type nanopore-raw \
    --merge_breaklen 1000 \
    "${FLYE_DIR}/assembly.fasta" \
    "${RAW_READS_DIR}/${SAMPLE_NAME}.fastq.gz" \
    "${CIRCLATOR_DIR}"

```

### CheckM2 analiza jakości genomu

```shell

conda activate checkm2

mkdir -p "${CHECKM2_DIR}"
mkdir -p "${CHECKM_INPUT_DIR}"

cp "${CIRCLATOR_DIR}/06.fixstart.fasta" "${CHECKM_INPUT_DIR}/${SAMPLE_NAME}.fasta"

checkm2 predict \
  --input "${CHECKM_INPUT_DIR}/${SAMPLE_NAME}.fasta" \
  --output-directory "${CHECKM2_DIR}" \
  --threads 4 \
  --force


```

### Bakta szybka i standaryzowana adnotacja genomów bakteryjnych i plazmidów

```shell

conda activate bakta-env

mkdir -p "${BAKTA_DIR}"

bakta --db "${BAKTA_DB}" \
      --verbose \
      --output "${BAKTA_DIR}" \
      --prefix legionella \
      --genus Legionella \
      --species pneumophila \
      --complete \
      --locus-tag legio \
      --threads 4 \
      --force \
      "${CIRCLATOR_DIR}/06.fixstart.fasta" 

```

### Quast analiza jakości złożenia genomu

```shell

mkdir -p "${QUAST_DIR}"

quast.py "${CHECKM_INPUT_DIR}/${SAMPLE_NAME}.fasta" \
  -o "${QUAST_DIR}" \
  -r "${REF_FILE}" \
  -g "${BAKTA_DIR}/${NAME}.gff3" \
  -t 4

```

### Eggnog-mapper funkcjonalna adnotacja sekwencji genomowych
#### http://eggnog-mapper.embl.de/ i użycie pliku z Bakta .faa

### Łączenie plików eggmapper i prokka files w jedną adnotację genomu
```python

import pandas as pd

simple_name = "R32_S257"
NAME = "legionella"

# Ścieżki dostępu do plików
prokka_path = f"/Legionella_pneumophila/long_reads/bakta/bakta_{simple_name}/{NAME}.tsv"
eggmaper_path = f"/Legionella_pneumophila/long_reads/eggmaper/{simple_name}/eggmaper_{simple_name}.tsv"
merged_path = f"/Legionella_pneumophila/long_reads/baktaxeggmaper/merge_{simple_name}.tsv"

# Wczytywanie plików
prokkatsv = pd.read_csv(prokka_path, sep='\t')
eggmapertsv = pd.read_csv(eggmaper_path, sep='\t')

# Łączenie
Output_df2 = pd.merge(prokkatsv, eggmapertsv, on='Locus_Tag', how='outer')
Output_df2.to_csv(merged_path, sep="\t", header=True, index=False)

```

## Wizualizacja danych

### Tworzenie boxplota w celu porównania kompletności oraz kontaminacji genomów powstałych w wyniku sekwencjonowania drugiej i trzeciej generacji

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Wprowadzenie danych
data = {
    "short_reads": ["L5", "L6", "L7", "L8", "L9", "L10", "L11", "L39", "L40", "R32"],
    "short_complitnes": [99.93, 99.93, 99.93, 99.93, 99.93, 99.93, 99.93, 99.59, 99.59, 99.93],
    "short_contamination": [0.11]*10,
    "long_reads": ["L5", "L6", "L7", "L8", "L9", "L10", "L11", "L39", "L40", "R32"],
    "long_complitnes": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
    "long_contamination": [0.69, 0.6, 0.6, 0.6, 0.6, 0.63, 0.57, 0.15, 0.15, 0.18]
}

df = pd.DataFrame(data)

# Przekształcenie do formatu długiego
df_long = pd.DataFrame({
    "Sample": df["short_reads"].tolist() + df["long_reads"].tolist(),
    "Read_Type": ["Short"] * 10 + ["Long"] * 10,
    "Contamination": df["short_contamination"].tolist() + df["long_contamination"].tolist()
})


# Tworzenie boxplota
plt.figure(figsize=(8, 6))
sns.boxplot(data=df_long, x="Read_Type", y="Contamination", palette="Set2")
sns.stripplot(data=df_long, x="Read_Type", y="Contamination", color="black", alpha=0.5, jitter=True)

plt.title("Porównanie kontaminacji między short_reads a long_reads")
plt.ylabel("Contamination (%)")
plt.xlabel("Typ odczytu")
plt.grid(True)
plt.tight_layout()
plt.show()

```

### Tworzenie wykresu COG
#### Rozdzielanie kolumny zawierającej dane do wykresu COG

```python

import pandas as pd
from collections import Counter

# Wczytaj plik .tsv
df = pd.read_csv("/Legionella_pneumophila/long_reads/eggmaper/cog_l5.tsv", sep="\t")

# Usuń puste wartości 
df = df.dropna(subset=["COG_category"])

# Rozdziel litery i zlicz wszystkie kategorie
all_categories = []
for entry in df["COG_category"]:
    all_categories.extend(list(entry.strip())) 

# Zlicz ile wystąpień ma każda litera
counter = Counter(all_categories)

# Zamień na DataFrame dla lepszego podglądu lub zapisu
category_counts = pd.DataFrame(counter.items(), columns=["Category", "Count"]).sort_values(by="Count", ascending=False)

# Wyświetl wynik
print(category_counts)

# (opcjonalnie) Zapisz do pliku
category_counts.to_csv("/Legionella_pneumophila/long_reads/eggmaper/L5_cog_Coun2.csv", index=False, sep="\t")
```

#### Tworzenie wykresu dla danych COG

```python

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# wczytywanie danych wejściowych
df = pd.read_csv("/Legionella_pneumophila/long_reads/eggmaper/L5_cog_Coun2.csv", sep="\t")

# Mapowanie kategorii COG na grupy funkcjonalne
grupa_map = {
    # Procesy i sygnały komórkowe
    'D': 'Procesy i sygnały komórkowe', 'M': 'Procesy i sygnały komórkowe',
    'N': 'Procesy i sygnały komórkowe', 'O': 'Procesy i sygnały komórkowe',
    'T': 'Procesy i sygnały komórkowe', 'U': 'Procesy i sygnały komórkowe',
    'V': 'Procesy i sygnały komórkowe', 'W': 'Procesy i sygnały komórkowe',
    'Z': 'Procesy i sygnały komórkowe',

    # Przechowywanie i procesowanie informacji
    'A': 'Przechowywanie i procesowanie informacji', 'J': 'Przechowywanie i procesowanie informacji',
    'K': 'Przechowywanie i procesowanie informacji', 'L': 'Przechowywanie i procesowanie informacji',
    'B': 'Przechowywanie i procesowanie informacji',

    # Metabolizm
    'C': 'Metabolizm', 'E': 'Metabolizm', 'F': 'Metabolizm',
    'G': 'Metabolizm', 'H': 'Metabolizm', 'I': 'Metabolizm',
    'P': 'Metabolizm', 'Q': 'Metabolizm',

    # Słabo scharakteryzowane
    'S': 'Słabo scharakteryzowane', 
}

# Dodanie kolumny z grupą
df['Grupa funkcjonalna'] = df['Category'].map(grupa_map)

# Sortowanie danych
df_sorted = df.sort_values(by=['Grupa funkcjonalna', 'Count'], ascending=[False, True])

# Ustawienia stylu wykresu
plt.figure(figsize=(10, 8))

# Wykres słupkowy poziomy
sns.barplot(
    data=df_sorted,
    y='Category',
    x='Count',
    hue='Grupa funkcjonalna',
    dodge=False,
    palette='Set2'
)

# Etykiety 
plt.xlabel('Liczba genów')
plt.ylabel('Kategoria COG')
plt.title('Liczba genów według kategorii COG i grupy funkcjonalnej')
plt.legend(title='Grupa funkcyjna', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Zapis
plt.savefig("/Legionella_pneumophila/long_reads/eggmaper/L5_cog_wykres.png", dpi=300)
plt.show()
```

### Tworzenie heatmapy dla genów systemów sekrecji typu II, IVA, IVB

#### Łączenie plików w celu stworzenia heatmapy

```python

import pandas as pd
import os

# Ścieżka do katalogu z plikami
folder_path = "/Legionella_pneumophila/long_reads/VFAnalyzer/sekretion_system"

# Lista nazw plików
file_names = [
    "wstep_L5.tsv", "wstep_L6.tsv", "wstep_L7.tsv", "wstep_L8.tsv", "wstep_L9.tsv",
    "wstep_L10.tsv", "wstep_L11.tsv", "wstep_L39.tsv", "wstep_L40.tsv", "wstep_R32.tsv"
]

dfs = []

for file_name in file_names:
    file_path = os.path.join(folder_path, file_name)
    
    # Wczytaj plik
    df = pd.read_csv(file_path, sep="\t")
    
    # Ustaw pierwszą kolumnę jako indeks
    df.set_index(df.columns[0], inplace=True)
    
    # Zachowaj drugą kolumnę pod jej oryginalną nazwą
    dfs.append(df.iloc[:, [0]])  # tylko druga kolumna (jako DataFrame, nie Series)

# Połącz wszystkie dane po indeksie
combined_df = pd.concat(dfs, axis=1)

# Dodaj wspólną kolumnę z pierwszego pliku
reference_df = pd.read_csv(os.path.join(folder_path, file_names[0]), sep="\t")
reference_df.set_index(reference_df.columns[0], inplace=True)
combined_df[reference_df.columns[1]] = reference_df.iloc[:, 1]

# Zapisz
combined_df.to_csv("/Legionella_pneumophila/long_reads/VFAnalyzer/sekretion_system/polaczone.tsv", sep="\t")

```

#### Tworzenie heatmapy 

```python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# wczytaj dane
df = pd.read_csv("/Legionella_pneumophila/short_reads/VFAnalizer/FILE/sekretion_system/polaczone.tsv", sep="\t", index_col=0)

# macierz obecności: 1 = obecność, 0 = brak ("-" lub puste)
presence_matrix = (df.notna() & (df != "-")).astype(int)

# 0 = niebieski, 1 = różowy
colors = ["#4A90E2", "#F78FB3"]
cmap = sns.color_palette(colors)

# rysowanie heatmapy
plt.figure(figsize=(13, 20))
ax = sns.heatmap(presence_matrix, cmap=cmap, cbar=False, linewidths=0.5, linecolor='grey', square=True, )
ax.set_yticks(range(len(presence_matrix.index)))
ax.set_yticklabels(presence_matrix.index, rotation=0)

plt.title("Heatmapa obecności genów")
plt.xlabel("Genomy")
plt.ylabel("Geny")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()

# wygląd legendy i żeby była po prawej
legend_labels = [
    mpatches.Patch(color=colors[0], label="Nie"),
    mpatches.Patch(color=colors[1], label="Tak")
]
ax.legend(
    handles=legend_labels,
    title="Obecność genu",
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),
    fontsize=10,
    title_fontsize=11,
    frameon=True
)

# bez ucinania legendy 
plt.tight_layout(rect=[0, 0, 0.85, 1])

# zapisywanie pliku w odpowiednim folderze
plt.savefig("/Legionella_pneumophila/short_reads/VFAnalizer/FILE/sekretion_system/heatmapa_short_.png", dpi=300, bbox_inches="tight")

plt.show()

```
