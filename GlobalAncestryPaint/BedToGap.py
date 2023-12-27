import argparse
import pandas as pd

# Funzione per aggiungere i colori agli ancestral populations
def add_colors(df, ancestry_colors):
    for ancestry, color in ancestry_colors.items():
        df[ancestry] = df[ancestry].apply(lambda x: f"{x:.5f} ({color})")
    return df

def main():
    # Definisci gli argomenti da linea di comando
    parser = argparse.ArgumentParser(description="Aggiungi colori alle ancestral populations in un file .bed")
    parser.add_argument("--input", help="File di input .bed", required=True)
    parser.add_argument("--output", help="Nome del file di output .bed", required=True)
    parser.add_argument("--ancestry1", nargs=2, help="Nome e colore per ancestry1", required=True)
    parser.add_argument("--ancestry2", nargs=2, help="Nome e colore per ancestry2", required=True)
    parser.add_argument("--ancestry3", nargs=2, help="Nome e colore per ancestry3", required=True)
    args = parser.parse_args()

    # Carica il file di input .bed utilizzando pandas
    df = pd.read_csv(args.input, sep='\t')

    # Crea un dizionario con le corrispondenze ancestry-color
    ancestry_colors = {
        args.ancestry1[0]: args.ancestry1[1],
        args.ancestry2[0]: args.ancestry2[1],
        args.ancestry3[0]: args.ancestry3[1]
    }

    # Aggiungi i colori alle ancestral populations
    df = add_colors(df, ancestry_colors)

    # Salva il DataFrame risultante nel file di output .bed
    df.to_csv(args.output, sep='\t', index=False)

    print(f"File {args.output} creato con successo.")

if __name__ == "__main__":
    main()
