import argparse
import pandas as pd
import matplotlib.pyplot as plt
import re

def plot_admixture(df, output_file):
    # Estrai i nomi degli individui dalla prima colonna (ignorando l'header)
    individuals = df.iloc[:, 0].tolist()

    # Estrai i dati di ancestral population
    ancestries = df.columns[1:]

    # Estrai le proporzioni per ciascuna ancestral population
    proportions = df.iloc[:, 1:].map(lambda x: float(re.match(r'^([0-9.]+)', x).group(1))).values

    # Inizializza un dizionario per le corrispondenze ancestry-color
    ancestry_colors = {}

    for ancestry in ancestries:
        # Estrai il colore dalla colonna e rimuovi i caratteri non numerici
        color = df[ancestry].str.extract(r'#([0-9a-fA-F]+)').iloc[0, 0]
        ancestry_colors[ancestry] = f'#{color}'

    # Crea un grafico a barre per ogni individuo
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = [0] * len(individuals)

    for i, ancestry in enumerate(ancestries):
        proportion = proportions[:, i]
        ax.bar(individuals, proportion, bottom=bottom, color=ancestry_colors[ancestry], label=ancestry)
        bottom = [sum(x) for x in zip(bottom, proportion)]

    ax.set_ylabel('Global Ancestry Proportion')
    ax.set_title('Global Ancestry Admixture')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=len(ancestries))  # Sposta la leggenda in alto a destra
    plt.xticks(rotation=90, fontsize=10)
    
    # Aggiungi una leggenda personalizzata sopra il grafico
    handles, labels = ax.get_legend_handles_labels()
    legend_labels = [label for label in labels]
    plt.legend(handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=len(ancestries))
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf')
    plt.close()

def main():
    # Definisci gli argomenti da linea di comando
    parser = argparse.ArgumentParser(description="Crea un grafico di Admixture dai dati di ancestral population")
    parser.add_argument("--input", help="File di input .bed", required=True)
    parser.add_argument("--output", help="Nome del file di output del grafico", required=True)
    args = parser.parse_args()

    # Carica il file di input .bed utilizzando pandas
    df = pd.read_csv(args.input, sep='\t')

    # Crea il grafico di Admixture
    plot_admixture(df, args.output)

if __name__ == "__main__":
    main()
