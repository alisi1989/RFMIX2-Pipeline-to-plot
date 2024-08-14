import argparse
import glob
import numpy as np
import os

def combine_files(prefix, chromosomes):
    """
    Combine files for each individual and chromosome specified.
    """
    individuals = set()
    
    # Create a list of files to combine and identify individuals
    files_to_combine = {}
    for chrom in chromosomes:
        # Create two patterns to search for files with and without underscore
        pattern1 = f"{prefix}*chr{chrom}.rfmix.Q"
        pattern2 = f"{prefix}_*chr{chrom}.rfmix.Q"
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

def calculate_mean(input_files, individual_name):
    """
    Calculate the mean values from the combined files and return the result as a string.
    """
    values = []
    header_lines = ""

    for i, input_file in enumerate(input_files):
        try:
            with open(input_file, 'r') as file:
                if i == 0:  # Take the first two lines from the first file
                    header_lines = next(file).strip() + "\n" + next(file).strip() + "\n"
                else:
                    next(file)  # Skip the first line
                    next(file)  # Skip the second line
                line = next(file).strip()  # Read the third line with the values
                values.append([float(val) for val in line.split("\t")[1:]])
        except FileNotFoundError:
            print(f"File {input_file} not found. Ensure that all files for the specified chromosomes are present.")
            return None, None
    
    # Calculate the mean per column
    mean_values = np.mean(values, axis=0)
    
    # Prepare the data for output
    data_line = individual_name + "\t" + "\t".join([f"{mean:.5f}" for mean in mean_values])
    
    return header_lines, data_line

def combine_rfmix_to_bed(prefix, chromosomes, output_path, sort_ancestry=None):
    """
    Combine .rfmix.Q files into a single .bed file, keeping the header from the second line.
    """
    combined_files, individuals = combine_files(prefix, chromosomes)
    all_lines_to_write = []
    combined_output_file_prefix = os.path.join(output_path, os.path.basename(prefix).strip('_'))

    ancestry_index = None
    for individual in individuals:
        header_lines, data_line = calculate_mean(combined_files[individual], os.path.basename(individual.strip('_')))
        if header_lines and data_line:
            if not all_lines_to_write:  # Only add header once
                all_lines_to_write.append(header_lines.strip())
                # Determine the index of the specified ancestry
                if sort_ancestry:
                    ancestry_columns = header_lines.strip().split("\t")[1:]  # Exclude the 'sample' column
                    if sort_ancestry in ancestry_columns:
                        ancestry_index = ancestry_columns.index(sort_ancestry) + 1  # Adjust for zero-based index
                    else:
                        print(f"Warning: Specified ancestry '{sort_ancestry}' not found in header. Sorting will not be applied.")
            all_lines_to_write.append(data_line)

    # Sort the data lines by the specified ancestry if provided
    if sort_ancestry and ancestry_index is not None:
        all_lines_to_write[1:] = sorted(all_lines_to_write[1:], key=lambda x: -float(x.split("\t")[ancestry_index]))

    # Write the combined .bed file
    bed_output_file = f"{combined_output_file_prefix}.bed"
    with open(bed_output_file, 'w') as output:
        for line in all_lines_to_write:
            output.write(line + '\n')
    
    print(f"File {bed_output_file} successfully created.")

def main():
    parser = argparse.ArgumentParser(description="Combine RFMix files and process them to generate .bed files.")
    parser.add_argument('--prefix', required=True, help="Prefix of the RFMix output file names.")
    parser.add_argument('--chr', nargs='+', default=[str(i) for i in range(1, 23)], help="Chromosomes of interest. Default is all 22.")
    parser.add_argument('--output', required=True, help="Output directory for the generated files.")
    parser.add_argument('--sort-ancestry', help="Column name to use for sorting in the .bed file.")
    args = parser.parse_args()

    combine_rfmix_to_bed(args.prefix, args.chr, args.output, args.sort_ancestry)

if __name__ == "__main__":
    main()
