#!/usr/bin/env python3
# Authors: Alessandro Lisi & Michael C. Campbell
# Make sure to install click beforehand (pip install click)

# Use "python BedToLAP.py --help" to see options

import os
import click

# Function to remove the 'chr' prefix, if present
def remove_chr_prefix(chromosome):
    if chromosome.startswith("chr"):
        return chromosome[3:]
    return chromosome

@click.command()
@click.option("--bed1", default=None, type=click.File("r"), help="First BED file for RFMix2 painting")
@click.option("--bed2", default=None, type=click.File("r"), help="Second BED file for RFMix2 painting")
@click.option("--ancestry0", default="#a32e2e", help="Color for ancestry0")
@click.option("--ancestry1", default="#0a0ae0", help="Color for ancestry1")
@click.option("--ancestry2", default="#bfa004", help="Color for ancestry2")
@click.option("--ancestry3", default="#d18311", help="Color for ancestry3")
@click.option("--ancestry4", default="#22ba9d", help="Color for ancestry4")
@click.option("--ancestry5", default="#839dfc", help="Color for ancestry5")
@click.option("--ancestry6", default="#9a5dc1", help="Color for ancestry6")
@click.option("--ancestry7", default="#26962b", help="Color for ancestry7")
@click.option("--ancestry8", default="#707070", help="Color for ancestry8")
@click.option("--ancestry9", default="#00cfff", help="Color for ancestry9")
@click.option("--ancestry10", default="#790ee0", help="Color for ancestry10")
@click.option("--unknown", default="#808080", help="Color for unknown regions, default gray #808080")
@click.option("--from-bp", "start", type=int, help="Start position (in bp) of your gene to highlight")
@click.option("--to-bp", "end", type=int, help="End position (in bp) of your gene to highlight")
@click.option("--chr", "chromosome", type=int, help="Chromosome number for the --start and --end positions")
@click.option("--feature-type", type=click.Choice(['line', 'triangle']), default='line', help="Feature type: 'line' for a dashed line or 'triangle' for a triangle marker to represent your target gene")
@click.option("--out", "output", default=None, type=click.File("w"), help="Output file name (extension .bed must be specified)")
def main(bed1, bed2, ancestry0, ancestry1, ancestry2, ancestry3, ancestry4, ancestry5, ancestry6, ancestry7, ancestry8, ancestry9, ancestry10, unknown, start, end, chromosome, feature_type, output):
    # Define the default colors
    default_colors = {
        "ancestry0": "#a32e2e",
        "ancestry1": "#0a0ae0",
        "ancestry2": "#bfa004",
        "ancestry3": "#d18311",
        "ancestry4": "#22ba9d",
        "ancestry5": "#839dfc",
        "ancestry6": "#9a5dc1",
        "ancestry7": "#26962b",
        "ancestry8": "#707070",
        "ancestry9": "#00cfff",
        "ancestry10": "#790ee0"
    }

    # Update the colors dictionary with the provided or default values
    colors = {
        "UNK": unknown,
        "ancestry0": ancestry0 or default_colors["ancestry0"],
        "ancestry1": ancestry1 or default_colors["ancestry1"],
        "ancestry2": ancestry2 or default_colors["ancestry2"],
        "ancestry3": ancestry3 or default_colors["ancestry3"],
        "ancestry4": ancestry4 or default_colors["ancestry4"],
        "ancestry5": ancestry5 or default_colors["ancestry5"],
        "ancestry6": ancestry6 or default_colors["ancestry6"],
        "ancestry7": ancestry7 or default_colors["ancestry7"],
        "ancestry8": ancestry8 or default_colors["ancestry8"],
        "ancestry9": ancestry9 or default_colors["ancestry9"],
        "ancestry10": ancestry10 or default_colors["ancestry10"]
    }

    output_lines = []
    
    # Read the header from bed1 and bed2
    header = bed1.readline().strip()
    output.write(header + "\n")

    # Process bed1 and bed2
    for line in bed1:
        if line.startswith('#'):  # Skip header/comment lines
            continue
        line = line.strip().split("\t")
        line[0] = remove_chr_prefix(line[0])  # Remove 'chr' prefix, if present
        ancestry = line[3]
        bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\tgeom_rect\t{colors.get(ancestry, unknown)}\t1"
        output_lines.append(bedLine)

    for line in bed2:
        if line.startswith('#'):  # Skip header/comment lines
            continue
        line = line.strip().split("\t")
        line[0] = remove_chr_prefix(line[0])  # Remove 'chr' prefix, if present
        ancestry = line[3]
        bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\tgeom_rect\t{colors.get(ancestry, unknown)}\t2"
        output_lines.append(bedLine)

    # Add lines with specified start and end positions and chromosome only if all are provided
    if start is not None and end is not None and chromosome is not None:
        feature_line_1 = f"{chromosome}\t{start}\t{end}\tgeom_{feature_type}\t#000000\t1"
        feature_line_2 = f"{chromosome}\t{start}\t{end}\tgeom_{feature_type}\t#000000\t2"
    
        # Add the lines representing the gene or the specific feature to the output_lines
        output_lines.extend([feature_line_1, feature_line_2])

    # Sort and write the output lines
    output_lines.sort(key=lambda x: int(x.split('\t')[1]))
    output.write("\n".join(output_lines) + "\n")

if __name__ == "__main__":
    main()
