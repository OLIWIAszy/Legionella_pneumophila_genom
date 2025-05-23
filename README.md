# Analiza genomu *Legionella pneumophila*

## Opis analizy
Celem analizy jest porównanie wyników sekwencjonowania wysokoprzepustowego, drugiej (Illumina) i trzeiej (Oxford Nanopore) generacji dla genomów *L. pneumophila*

Analiza obejmuje:

- **Kontrolę jakości i obróbkę odczytów**
- **Składanie genomu**
- **Ocena jakości genomów**
- **Adnotacja funkcjonalna**
- **Identyfikacja genów wirulencji i systemów sekrecji**
- **Wizualizacje (heatmapy, wykresy completeness/contamination)**

## Narzędzia
- `fastp`, `Filtlong`, `Nanoplot`
- `Flye`, `Unicycler`
- `CheckM`, `CheckM2`, `QUAST`
- `Prokka`
- `eggNOG-mapper`
- `VFAnalyzer`
- `Python`, `pandas`, `matplotlib`, `seaborn`


## Wymagania
- Python 3.8+, zainstalowane biblioteki: `pandas`, `matplotlib`, `seaborn`, `Bio`

## Wizualizacje
Pliki `heatmap_*.png` i `completeness_plot.png` znajdują się w katalogu `schematy/`.

