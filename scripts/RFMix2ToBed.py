import pandas as pd
import argparse
import os

# Constants for columns to keep or drop
COLUMNS_TO_KEEP = [0, 1, 2, 6, 7]
COLUMNS_TO_DROP = [3, 4, 5]

# Mapping of ancestries
ANCESTRY_MAPPING = {0: "ancestry0", 1: "ancestry1", 2: "ancestry2", 3: "ancestry3", 4: "ancestry4", 5: "ancestry5", 6: "ancestry6", 7: "ancestry7"}

def process_file(input_file, output_file_prefix):
    # Check the input file extension
    if not input_file.endswith(".msp.tsv"):
        print("Error: The input file must have the extension .msp.tsv.")
        return

    # Check the output file name extension
    if output_file_prefix.endswith(".bed"):
        print("Error: The output file name should not have the extension .bed.")
        return

    # Read the file, skip the first two rows
    msp = pd.read_csv(input_file, header=None, sep='\t', skiprows=2)

    # Select only the first 8 columns
    msp2 = msp.iloc[:, :8]

    # Drop unwanted columns
    msp2.drop(columns=COLUMNS_TO_DROP, inplace=True)

    # Rename columns
    msp2.columns = ["chm", "spos", "epos", "ind1", "ind2"]

    # Create 2 dataframes
    ch1 = msp2.drop(columns=["ind2"]).copy()
    ch2 = msp2.drop(columns=["ind1"]).copy()

    # Replace values in columns "ind1" and "ind2" with corresponding ancestry values
    ch1['ind1'] = ch1['ind1'].replace(ANCESTRY_MAPPING)
    ch2['ind2'] = ch2['ind2'].replace(ANCESTRY_MAPPING)

    # Save the two files for each individual
    ch1.to_csv(f"{output_file_prefix}_hap1.bed", index=False, header=False, sep='\t', quoting=3)
    ch2.to_csv(f"{output_file_prefix}_hap2.bed", index=False, header=False, sep='\t', quoting=3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to process MSP file and generate haplotype BED files.")
    parser.add_argument("--input", help="Input file. Extension .msp.tsv must be specified.", required=True)
    parser.add_argument("--output", help="Output file name. Extension .bed must NOT be specified.", required=True)
    
    args = parser.parse_args()

    process_file(args.input, args.output)
