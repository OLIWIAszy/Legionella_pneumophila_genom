# Legionella_pneumophila_genom
Analiza genomu bakterii Legionella pneumophila, po sekwencjonowaniu wysokoprzepustowym NGS na platformie Illumina.

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

### Łączenie plików eggmaster i prokka files w jedną adnotację genomu
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
