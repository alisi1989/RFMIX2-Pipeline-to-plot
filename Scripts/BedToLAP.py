#!/usr/bin/env python3
#Autors Alessandro Lisi & Michael C. Campbell
#be sure to priorly install click (pip install click) 

#use python BedToLAP.py --help to see the option

import os
import click

@click.command()
@click.option(
    "--bed1",
    default=None,
    type=click.File("r"),
    help="Bed 1 RFMix2 painting",
)
@click.option(
    "--bed2",
    default=None,
    type=click.File("r"),
    help="Bed TO RFMix2 painting",
)
@click.option("--ancestry0", default="#0000ff", help="Color for ancestry0, default blue #0000ff")
@click.option("--ancestry1", default="#850b39", help="Color for ancestry1, default magenta #850b39")
@click.option("--ancestry2", default="#F4A500", help="Color for ancestry2, default gold #F4A500")
@click.option("--ancestry3", default="#292518", help="Color for ancestry3, default dark gold #292518")
@click.option("--ancestry4", default="#076940", help="Color for ancestry4, default green #076940")
@click.option("--unknown", default="#808080", help="Color for Unknown regions, default gray #808080")
@click.option("--from-bp", "start", type=int, help="Start position (in bp) of your Gene to Highlight")
@click.option("--to-bp", "end", type=int, help="End position (in bp) of your Gene to Highlight")
@click.option("--chr", "chromosome", type=int, help="Chromosome number to use for --start and --end positions")
@click.option("--feature-type", type=click.Choice(['line', 'triangle']), default='line', help="Feature type: 'line' for dashed line to represent your target Gene, 'triangle' for triangle marker to represent your target Gene")
@click.option(
    "--out",
    "output",
    default=None,
    type=click.File("w"),
    help=".bed must be specified",
)
def main(bed1, bed2, ancestry0, ancestry1, ancestry2, ancestry3, ancestry4, unknown, start, end, chromosome, feature_type, output):
    colors = {"UNK": unknown, "ancestry0": ancestry0, "ancestry1": ancestry1, "ancestry2": ancestry2, "ancestry3": ancestry3, "ancestry4": ancestry4}
    output_lines = []
    
        # Leggi l'intestazione da bed1 e bed2
    header = bed1.readline().strip()
    output.write(header + "\n")
    

    # Process bed1 and bed2
    for line in bed1:
        if line.startswith('#'):  # Salta le righe di intestazione/commento
            continue
        line = line.strip().split("\t")

        if feature_type == 'line':
            bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t1"
        else:
            bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t1"
        output_lines.append(bedLine)

    for line in bed2:
        if line.startswith('#'):  # Salta le righe di intestazione/commento
            continue
        line = line.strip().split("\t")
        if feature_type == 'line':
            bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t2"
        else:
            bedLine = f"{line[0]}\t{line[1]}\t{line[2]}\t0\t1\t{colors[line[3]]}\t2"
        output_lines.append(bedLine)

    # Add lines with specified start and end positions and chromosome
    if start is not None and end is not None and chromosome is not None:
        if feature_type == 'line':
            extra_lines = [
                f"{chromosome}\t{start}\t{end}\t1\t1\t#000000\t1",
                f"{chromosome}\t{start}\t{end}\t1\t1\t#000000\t2"
            ]
        else:
            extra_lines = [
                f"{chromosome}\t{start}\t{end}\t2\t1\t#000000\t1",
                f"{chromosome}\t{start}\t{end}\t2\t1\t#000000\t2"
            ]
        output_lines.extend(extra_lines)

    # Ordina e scrivi le linee di output
    output_lines.sort(key=lambda x: int(x.split('\t')[-1]))
    output.write("\n".join(output_lines))

if __name__ == "__main__":
    main()