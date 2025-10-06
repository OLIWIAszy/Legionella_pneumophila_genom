# Analiza genomu *Legionella pneumophila*

## Opis analizy
Celem analizy jest porównanie wyników sekwencjonowania wysokoprzepustowego, [drugiej (Illumina)](./analiza_skrypt/analiza_short_reads.md) i [trzeciej (Oxford Nanopore)](./analiza_skrypt/analiza_long_reads.md) generacji dla genomów *L. pneumophila*. Analizy wykonane były przy użyciu skryptu

Analiza obejmuje:

- **Kontrolę jakości i obróbkę odczytów**
- **Składanie genomu**
- **Ocena jakości genomów**
- **Adnotacja funkcjonalna**
- **Identyfikacja genów wirulencji i systemów sekrecji**
- **[Wizualizacje](./analiza_skrypt/Wizualizacje_danych.md) (heatmapy [krótkich](./Wizualizacje/heatmapa_genow_short_reads.png) i [długich](./Wizualizacje/heatmapa_genow_long_reads.png) odczytów, wykresy [completeness](./Wizualizacje/porownanie_kompletnosci.png)[contamination](./Wizualizacje/porownanie_kontaminacji.png), [wykres grup funkcjonalnych COG](./Wizualizacje/L5_cog_wykres.png))**

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


