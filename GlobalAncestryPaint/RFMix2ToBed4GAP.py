import os
import argparse

# Funzione per ottenere la lista di terze righe dai file .rfmix.Q
def get_third_lines(file_list):
    lines_to_write = []
    for filename in file_list:
        with open(filename, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 3:
                lines_to_write.append(lines[2])
    return lines_to_write

def main():
    # Definisci gli argomenti da linea di comando
    parser = argparse.ArgumentParser(description="Combina i file .rfmix.Q in un file .bed")
    parser.add_argument("--input", nargs='+', help="Lista dei file .rfmix.Q di input", required=True)
    parser.add_argument("--output", help="Nome del file di output .bed", required=True)
    args = parser.parse_args()

    # Ottieni le terze righe dai file di input
    lines_to_write = get_third_lines(args.input)

    # Scrivi le prime due righe e le terze righe nel file di output .bed
    if len(lines_to_write) > 0:
        header1 = ""
        with open(args.input[0], 'r') as file:
            header1 = file.readline()
            header2 = file.readline()
        
        with open(args.output, 'w') as output:
            output.write(header2)
            for line in lines_to_write:
                output.write(line)
        
        print(f"File {args.output} creato con successo.")
    else:
        print("Nessuna terza riga trovata nei file di input.")

if __name__ == "__main__":
    main()
