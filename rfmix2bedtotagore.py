#!/usr/bin/env python3
__author__ = "Aroon Chande edit by Alessandro Lisi"
__copyright__ = "Copyright 2019, Applied Bioinformatics Lab edit by Michael Campbell Lab"
__license__ = "MIT" 

import os

import click

commentBlock = """
# Each column is explained below:
# 1. chr - The chromosome on which a feature has to be drawn
# 2. start - Start position (in bp) for feature
# 3. stop - Stop position (in bp) for feature
# 4. feature - The shape of the feature to be drawn
# 	 0 - will draw a rectangle
# 	 1 - will draw a circle
# 	 2 - will draw a triangle pointing to the genomic location
# 	 3 - will draw a line at that genomic location
# 5. size - The horizontal size of the feature. Should range between 0 and 1.
# 6. color - Specify the color of the genomic feature with a hex value (#FF0000 for red, etc.)
# 7. chrCopy - Specify the chromosome copy on which the feature should be drawn (1 or 2).  To draw the same feature on both chromosomes, you must specify the feature twice
"""

@click.command()
@click.option(
    "--chr1",
    "-1",
    "chr1",
    default=None,
    type=click.File("r"),
    help="Chromosome 1 RFMix painting",
)
@click.option(
    "--chr2",
    "-2",
    "chr2",
    default=None,
    type=click.File("r"),
    help="Chromosome 2 RFMix painting",
)
@click.option("--anc0", "ANC0", default="#0000ff", help="Color for ANC1 blocks")
@click.option("--anc1", "ANC1", default="#0000ff", help="Color for ANC1 blocks")
@click.option("--anc2", "ANC2", default="#F4A500", help="Color for ANC2 blocks")
@click.option("--anc3", "ANC3", default="#F4A500", help="Color for ANC3 blocks")
@click.option("--anc4", "ANC4", default="#F4A500", help="Color for ANC4 blocks")

@click.option("--unk", "Unknown", default="#808080", help="Color for Unknown regions")
@click.option(
    "--out",
    "-o",
    "output",
    default=None,
    type=click.File("w"),
    help="Output da Vinci bed",
)
def main(chr1, chr2, ANC0, ANC1, ANC2, ANC3, ANC4, Unknown, output):
    colors = {"UNK": Unknown, "ANC0": ANC0, "ANC1": ANC1, "ANC2": ANC2, "ANC3": ANC3, "ANC4": ANC4}
    output.write("#chr\tstart\tstop\tfeature\tsize\tcolor\tchrCopy")
    output.write(commentBlock)
    for line in chr1:
        line = line.strip().split("\t")
        bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t1"
        output.write(bedLine + "\n")
    for line in chr2:
        line = line.strip().split("\t")
        bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t2"
        output.write(bedLine + "\n")

if __name__ == "__main__":
    main()
