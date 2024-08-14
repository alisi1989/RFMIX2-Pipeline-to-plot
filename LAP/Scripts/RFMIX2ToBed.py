import argparse
import glob
import os
import pandas as pd

# Mapping of ancestries extended to ancestry10
ANCESTRY_MAPPING = {
    0: "ancestry0", 1: "ancestry1", 2: "ancestry2", 3: "ancestry3", 4: "ancestry4",
    5: "ancestry5", 6: "ancestry6", 7: "ancestry7", 8: "ancestry8", 9: "ancestry9", 10: "ancestry10"
}

def combine_chromosomes(prefix, chromosomes):
    """
    Combine chromosome files into a single file for each individual.
    """
    individuals = set()
    
    # Create a list of files to combine and identify individuals
    files_to_combine = {}
    for chrom in chromosomes:
        # Create two patterns to search for files with and without underscore
        pattern1 = f"{prefix}*chr{chrom}.msp.tsv"
        pattern2 = f"{prefix}_*chr{chrom}.msp.tsv"
        matching_files1 = glob.glob(pattern1)
        matching_files2 = glob.glob(pattern2)
        matching_files = matching_files1 + matching_files2

        for file in matching_files:
            individual = file.split('chr')[0]
            individuals.add(individual)
            if individual not in files_to_combine:
                files_to_combine[individual] = []
            files_to_combine[individual].append(file)

    return files_to_combine, sorted(individuals)

def process_file(input_file, output_file_prefix):
    # Check the input file extension
    if not input_file.endswith(".msp.tsv"):
        print("Error: The input file must have the extension .msp.tsv.")
        return

    # Initialize a list to store header lines
    headers = []

    # Read the file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Process each line
    data_lines = []
    for line in lines:
        if line.startswith("#"):
            headers.append(line.strip())
        else:
            data_lines.append(line.strip().split('\t'))

    # Convert the processed data lines to a DataFrame
    msp = pd.DataFrame(data_lines, columns=["chm", "spos", "epos", "sgpos", "egpos", "n_snps", "ind1", "ind2"])

    # Apply ancestry mapping to both columns "ind1" and "ind2"
    msp['ind1'] = msp['ind1'].astype(int).map(ANCESTRY_MAPPING).fillna(msp['ind1'])
    msp['ind2'] = msp['ind2'].astype(int).map(ANCESTRY_MAPPING).fillna(msp['ind2'])

    # Create 2 dataframes for haplotypes
    ch1 = msp.drop(columns=["ind2", "sgpos", "egpos", "n_snps"]).copy()
    ch2 = msp.drop(columns=["ind1", "sgpos", "egpos", "n_snps"]).copy()

    # Rename columns for output
    ch1.columns = ["chm", "spos", "epos", "ancestry"]
    ch2.columns = ["chm", "spos", "epos", "ancestry"]

    # Save the two files for each individual, including the headers
    with open(f"{output_file_prefix}_hap1.bed", 'w') as file:
        file.write("\n".join(headers) + '\n')
        ch1.to_csv(file, index=False, header=False, sep='\t', quoting=3)

    with open(f"{output_file_prefix}_hap2.bed", 'w') as file:
        file.write("\n".join(headers) + '\n')
        ch2.to_csv(file, index=False, header=False, sep='\t', quoting=3)

    print(f"Haplotype files created: {output_file_prefix}_hap1.bed and {output_file_prefix}_hap2.bed")

def main():
    parser = argparse.ArgumentParser(description="Combine RFMix v.2 files and generate haplotype BED files.")
    parser.add_argument('--prefix', required=True, help="Prefix of the RFMix v.2 output file names.")
    parser.add_argument('--chr', nargs='+', default=[str(i) for i in range(1, 23)], help="Chromosomes of interest. Default is all 22.")
    parser.add_argument('--output', required=True, help="Output directory for the generated files.")
    args = parser.parse_args()

    files_to_combine, individuals = combine_chromosomes(args.prefix, args.chr)

    for individual in individuals:
        combined_output_file = os.path.join(args.output, f"{os.path.basename(individual).strip('_')}.msp.tsv")
        
        # Combine files for the current individual
        with open(combined_output_file, 'w') as outfile:
            header_written = False
            for file in files_to_combine[individual]:
                with open(file, 'r') as infile:
                    lines = infile.readlines()
                    if not header_written:
                        outfile.write(lines[0])  # Write the first header line (#Subpopulation)
                        outfile.write(lines[1])  # Write the second header line (#chr)
                        header_written = True
                    for line in lines[2:]:  # Start writing from the third line onwards
                        outfile.write(line)
        
        # Process the combined file to generate haplotype BED files
        process_file(combined_output_file, combined_output_file.replace('.msp.tsv', ''))

        # Remove the combined file
        os.remove(combined_output_file)
        print(f"Combined file removed: {combined_output_file}")

if __name__ == "__main__":
    main()
