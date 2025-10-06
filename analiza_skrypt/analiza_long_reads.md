# Legionella_pneumophila_genom
Analiza genomu bakterii Legionella pneumophila, po sekwencjonowaniu wysokoprzepustowym trzeciej generacji opartym na nanoporach.


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

### Łączenie plików eggmaster i prokka files w jedną adnotację genomu
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
