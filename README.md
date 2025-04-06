# Legionella_pneumophila_genom
Analiza genomu bakterii Legionella pneumophila, po sekwencjonowaniu wysokoprzepustowym.


## Sekwencjonowanie NGS na platformie Illumina 

### Dane wejściowe 

```shell

# Input data
conda activate bioinfo 



BACTERIA_FOLDER="/media/oli/Nowy/bioinformatyka/analiza_DNA/Legionella_pneumophila/short_reads" 
SAMPLE_NAME='L7_S252'


# Define the directories
RAW_READS_DIR="$BACTERIA_FOLDER/raw_reads"
FASTP_DIR="$BACTERIA_FOLDER/fastp"
UNICYCLER_DIR="$BACTERIA_FOLDER/unicycler/unicycler_${SAMPLE_NAME}"
REFERENCE_DIR="$BACTERIA_FOLDER/refgenom"
PROKKA_DIR="$BACTERIA_FOLDER/prokka"
QUAST_DIR="$BACTERIA_FOLDER/quast"
CHECKM_DIR="$BACTERIA_FOLDER/checkm"


REF_FILE="$REFERENCE_DIR/GCF_002252505.1_ASM225250v1_genomic.fna




```

### Fastp analiza jakości odczytów

```shell

mkdir /media/oli/Nowy/bioinformatyka/analiza_DNA/Legionella_pneumophila/short_reads/fastp

fastp -i "$RAW_READS_DIR/${SAMPLE_NAME}_R1_001.fastq.gz" \
      -I "$RAW_READS_DIR/${SAMPLE_NAME}_R2_001.fastq.gz" \
      -o "$FASTP_DIR/fastp_${SAMPLE_NAME}_R1_001.fastq.gz" \
      -O "$FASTP_DIR/fastp_${SAMPLE_NAME}_R2_001.fastq.gz" \
      --json "$FASTP_DIR/${SAMPLE_NAME}_report.json" \
      --html "$FASTP_DIR/${SAMPLE_NAME}_report.html"

```
### Unicycler składanie genomu

```shell

unicycler -1 "$FASTP_DIR/fastp_${SAMPLE_NAME}_R1_001.fastq.gz" \
          -2 "$FASTP_DIR/fastp_${SAMPLE_NAME}_R2_001.fastq.gz" \
          -o "$UNICYCLER_DIR" \
          --min_fasta_length 1000

```

### Ragtag ustawianie kontigów

```shell

ragtag.py scaffold \
          "$REF_FILE"\
          "$UNICYCLER_DIR/assembly.fasta" \
          -o ragtag

```

### Prokka adnotacja genomu

```shell

prokka --outdir "$PROKKA_DIR" --force \
       --prefix "$SAMPLE_NAME" --addgenes \
       --kingdom "Bacteria" \
       --genus "Brucella" \
       --species "grignonensis" \
       --cpus 4 \
       "$UNICYCLER_DIR/assembly.fasta"

```

### Quast analiza jakości złożenia genomu

```shell

quast.py -o "$QUAST_DIR" \
        "$UNICYCLER_DIR/assembly.fasta" \
        -r "$REF_FILE" \
        -1 "$INPUT_FILE1" \
        -2 "$INPUT_FILE2" \
        -g "$PROKKA_DIR/${SAMPLE_NAME}.gff"

```

### CheckM analiza jakości adnotacji genomu

```shell

checkm taxonomy_wf genus Brucella \
                   --genes \
                   -t 10 \
                   -x faa \
                   --file ${BACTERIA_FOLDER}/checkm_stats.txt \
                   "$PROKKA_DIR" \
                   "$CHECKM_DIR"
                   
```

### Eggnog-mapper funkcjonalna adnotacja sekwencji genomowych
#### http://eggnog-mapper.embl.de/ i użycie pliku z Prokka .faa

### Łączenie plików eggmaster i prokka files w jedną adnotację genomu
```python

import pandas as pd

prokkatsv = pd.read_csv("/home/oli/Documents/analiza_DNA/Brucella_grignonensis/prokka/OA139_S10.tsv", sep='\t')
eggmapertsv = pd.read_csv("/home/oli/Documents/analiza_DNA/Brucella_grignonensis/MM_4sp0uox_.emapper.annotations(1).tsv", sep='\t')

Output_df2 = pd.merge(prokkatsv, eggmapertsv, on='locus_tag', how='outer')

Output_df2.to_csv("/home/oli/Documents/analiza_DNA/Brucella_grignonensis/Outputgenes2.tsv",
                 sep="\t", header=True,
                 index=False)

outputgenes = pd.read_csv("/home/oli/Documents/analiza_DNA/Brucella_grignonensis/Outputgenes2.tsv", sep='\t')



allright = outputgenes[outputgenes['ftype'] != 'gene']

allright.to_csv("/home/oli/Documents/analiza_DNA/Brucella_grignonensis/prokkaandeggmaper.tsv",
                 sep="\t", header=True,
                 index=False)

```




## Sekwencjonowanie trzeciej generacji oparte na nanoporach

### Dane wejściowe
```shell

# Input data
conda activate bioinfo 



BACTERIA_FOLDER="/home/oli/Documents/analiza_DNA/Legionella_pneumophila/long_reads" 
SAMPLE_NAME='L39'

mkdir /home/oli/Documents/analiza_DNA/Legionella_pneumophila/long_reads/flye/flye_${SAMPLE_NAME}

# Define the directories
RAW_READS_DIR="$BACTERIA_FOLDER/raw_reads"

FLYE_DIR="$BACTERIA_FOLDER/flye/flye_${SAMPLE_NAME}"
CIRCLATOR_DIR="$BACTERIA_FOLDER/circlator"
BAKTA_DIR="$BACTERIA_FOLDER/bakta"
QUAST_DIR="$BACTERIA_FOLDER/quast"
CHECKM_DIR="$BACTERIA_FOLDER/checkm"

# Define the files
INPUT_FILE="$RAW_READS_DIR/${SAMPLE_NAME}



```

### Flye składanie sekwencji nukleotydowych długich pojedyńczych cząsteczek

```shell

flye --nano-raw "$RAW_READS_DIR/${SAMPLE_NAME}.fastq.gz" \
     --out-dir "$FLYE_DIR" \
     --threads 4

```

### Circlator cyrkulacja genomu

```shell

circlator all --merge_min_id 85 \
          --merge_breaklen 1000 \
          --uniprot_search dnaa \
          assembly.fasta \
          reads "$CIRCLATOR_DIR" \

```

### Bakta szybka i standaryzowana adnotacja genomów bakteryjnych i plazmidów

```shell

bakta --db <db-path> \
      --verbose \
      --output results/ \
      --prefix ecoli123 \
      --locus-tag eco634 \
      --prodigal-tf eco.tf \
      --replicons replicon.tsv 
      --threads 8 genome.fasta



```

### Quast analiza jakości złożenia genomu

```shell

quast.py -o "$QUAST_DIR" \
        "$UNICYCLER_DIR/assembly.fasta" \
        -r "$REF_FILE" \
        -1 "$INPUT_FILE1" \
        -g "$BAKTA_DIR/${SAMPLE_NAME}.gff"

```

### CheckM analiza jakości adnotacji genomu

```shell

checkm taxonomy_wf genus Brucella \
                   --genes \
                   -t 10 \
                   -x faa \
                   --file ${BACTERIA_FOLDER}/checkm_stats.txt \
                   "$BAKTA_DIR" \
                   "$CHECKM_DIR"
```

### Eggnog-mapper funkcjonalna adnotacja sekwencji genomowych
#### http://eggnog-mapper.embl.de/ i użycie pliku z Bakta .faa

