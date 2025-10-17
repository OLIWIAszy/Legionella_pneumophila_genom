# Legionella_pneumophila_genom 
Wizualizacja danych genomu bakterii *Legionella pneumophila*, po sekwencjonowaniu wysokoprzepustowym.


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
