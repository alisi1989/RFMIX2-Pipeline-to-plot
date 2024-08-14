#!/usr/bin/env python3
#Autors Alessandro Lisi & Michael C. Campbell
#be sure to priorly install librsvg (brew install rsvg) or for Linux machine "https://manpages.ubuntu.com/manpages/trusty/man1/rsvg-convert.1.html"

#use python LAP.py --help to see the option

import argparse
import xml.etree.ElementTree as ET
import os
import subprocess

def insert_colored_regions(svg_file, bed_file, output_file, individual_name, ancestries, chromosome_length):

    chromosome_coordinates = {
        'chromosome1_hap1': {'x': 224.900, 'y': 103.20, 'width': 21, 'height': 666.300},
        'chromosome1_hap2': {'x': 245.900, 'y': 103.20, 'width': 21, 'height': 666.300},
        'chromosome2_hap1': {'x': 287.900, 'y': 121.300, 'width': 21, 'height': 648.200},
        'chromosome2_hap2': {'x': 308.900, 'y': 121.300, 'width': 21, 'height': 648.200},
        'chromosome3_hap1': {'x': 350.800, 'y': 238.600, 'width': 21, 'height': 530.900},
        'chromosome3_hap2': {'x': 371.700, 'y': 238.600, 'width': 21, 'height': 530.900},
        'chromosome4_hap1': {'x': 413.500, 'y': 260.200, 'width': 21, 'height': 509.300},
        'chromosome4_hap2': {'x': 434.200, 'y': 260.200, 'width': 21, 'height': 509.300},
        'chromosome5_hap1': {'x': 476.300, 'y': 283.300, 'width': 21, 'height': 486.200},
        'chromosome5_hap2': {'x': 497.200, 'y': 283.300, 'width': 21, 'height': 486.200},
        'chromosome6_hap1': {'x': 539.100, 'y': 312.000, 'width': 21, 'height': 457.500},
        'chromosome6_hap2': {'x': 560.000, 'y': 312.000, 'width': 21, 'height': 457.500},
        'chromosome7_hap1': {'x': 601.900, 'y': 342.600, 'width': 21, 'height': 426.900},
        'chromosome7_hap2': {'x': 622.800, 'y': 342.600, 'width': 21, 'height': 426.900},
        'chromosome8_hap1': {'x': 664.600, 'y': 380.600, 'width': 21, 'height': 388.900},
        'chromosome8_hap2': {'x': 685.500, 'y': 380.600, 'width': 21, 'height': 388.900},
        'chromosome9_hap1': {'x': 727.200, 'y': 398.600, 'width': 21, 'height': 370.900},
        'chromosome9_hap2': {'x': 748.100, 'y': 398.600, 'width': 21, 'height': 370.900},
        'chromosome10_hap1': {'x': 790.200, 'y': 410.800, 'width': 21, 'height': 358.600},
        'chromosome10_hap2': {'x': 811.100, 'y': 410.800, 'width': 21, 'height': 358.600},
        'chromosome11_hap1': {'x': 853.000, 'y': 407.300, 'width': 21, 'height': 362.000},
        'chromosome11_hap2': {'x': 873.900, 'y': 407.300, 'width': 21, 'height': 362.00},
        'chromosome12_hap1': {'x': 915.700, 'y': 412.300, 'width': 21, 'height': 357.200},
        'chromosome12_hap2': {'x': 936.700, 'y': 412.300, 'width': 21, 'height': 357.200},
        'chromosome13_hap1': {'x': 978.300, 'y': 462.600, 'width': 21, 'height': 306.700},
        'chromosome13_hap2': {'x': 999.200, 'y': 462.600, 'width': 21, 'height': 306.700},
        'chromosome14_hap1': {'x': 1041.300, 'y': 482.400, 'width': 21, 'height': 287.100},
        'chromosome14_hap2': {'x': 1062.200, 'y': 482.400, 'width': 21, 'height': 287.100},
        'chromosome15_hap1': {'x': 1104.300, 'y': 495.700, 'width': 21, 'height': 273.600},
        'chromosome15_hap2': {'x': 1125.300, 'y': 495.700, 'width': 21, 'height': 273.600},
        'chromosome16_hap1': {'x': 1166.300, 'y': 527.100, 'width': 21, 'height': 242.400},
        'chromosome16_hap2': {'x': 1187.300, 'y': 527.100, 'width': 21, 'height': 242.400},
        'chromosome17_hap1': {'x': 1229.600, 'y': 545.700, 'width': 21, 'height': 223.500},
        'chromosome17_hap2': {'x': 1250.500, 'y': 545.700, 'width': 21, 'height': 223.500},
        'chromosome18_hap1': {'x': 1292.200, 'y': 553.500, 'width': 21, 'height': 215.800},
        'chromosome18_hap2': {'x': 1313.100, 'y': 553.500, 'width': 21, 'height': 215.800},
        'chromosome19_hap1': {'x': 1355.100, 'y': 611.800, 'width': 21, 'height': 157.700},
        'chromosome19_hap2': {'x': 1376.100, 'y': 611.800, 'width': 21, 'height': 157.700},
        'chromosome20_hap1': {'x': 1417.900, 'y': 596.200, 'width': 21, 'height': 173.300},
        'chromosome20_hap2': {'x': 1438.700, 'y': 596.200, 'width': 21, 'height': 173.300},
        'chromosome21_hap1': {'x': 1480.700, 'y': 643.300, 'width': 21, 'height': 125.900},
        'chromosome21_hap2': {'x': 1501.600, 'y': 643.300, 'width': 21, 'height': 125.900},
        'chromosome22_hap1': {'x': 1543.500, 'y': 632.700, 'width': 21, 'height': 136.800},
        'chromosome22_hap2': {'x': 1564.400, 'y': 632.700, 'width': 21, 'height': 136.800},
        
    }

    # Chromosomes Lenghts
    chromosome_lengths = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]

    # Parsing the original SVG file
    original_svg_tree = ET.parse(svg_file)
    original_svg_root = original_svg_tree.getroot()

    # Add the individual's name as a title in the SVG file
    title_element = ET.Element('text', {'x': "800", 'y': "30", 'fill': "black", 'font-size': "32"})
    title_element.text = individual_name
    original_svg_root.insert(0, title_element)  # Insert the title at the beginning of the document

# Create a legend for specified ancestries
    if ancestries:
        y_offset = 40  # Start drawing the legend from this y coordinate
        for name, color in ancestries.items():
            rect_element = ET.Element('rect', {'x': "1450", 'y': str(y_offset), 'width': "25", 'height': "15", 'fill': color})
            original_svg_root.insert(0, rect_element)

            # Text for the ancestry
            text_element = ET.Element('text', {'x': "1480", 'y': str(y_offset + 12), 'fill': "black", 'font-size': "16"})
            text_element.text = name
            original_svg_root.insert(0, text_element)

            y_offset += 30 # Update the y coordinate for the next entry in the legend

            
    rectangle_elements = []
    line_elements = []
    
    # Extract regions from the BED file and add them to the SVG
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip()
            # Ignora le righe di commento che iniziano con "#"
            if line.startswith("#"):
                continue
                
            parts = line.split('\t')
            if len(parts) < 6:
                continue

            chromosome, start, end, line_type, color, haplotype = parts[:6]
            chromosome = int(chromosome)
            start, end, haplotype = int(start), int(end), int(haplotype)

            chromosome_id = f'chromosome{chromosome}_hap{haplotype}'
            if chromosome_id in chromosome_coordinates:
                chromosome_data = chromosome_coordinates[chromosome_id]
                chromosome_length = chromosome_lengths[chromosome - 1]

                x = chromosome_data['x']
                width = chromosome_data['width']
                start_y = chromosome_data['y'] + (start / chromosome_length) * chromosome_data['height']
                end_y = chromosome_data['y'] + (end / chromosome_length) * chromosome_data['height']

                if line_type == 'geom_line':
                    offset=15
                    line_element_start = ET.Element('line', {
                        'x1': str(x - offset),
                        'y1': str(start_y),
                        'x2': str(x + width + offset),
                        'y2': str(start_y),
                        'stroke': color,
                        'stroke-width': '2',
                        'stroke-dasharray': '10'
                    })
                    line_element_end = ET.Element('line', {
                        'x1': str(x - offset),
                        'y1': str(end_y),
                        'x2': str(x + width + offset),
                        'y2': str(end_y),
                        'stroke': color,
                        'stroke-width': '2',
                        'stroke-dasharray': '10'
                    })
                    original_svg_root.insert(0, line_element_start)
                    original_svg_root.insert(0, line_element_end)
                else:
                    rect_element = ET.Element('rect', {
                        'x': str(x),
                        'y': str(start_y),
                        'width': str(width),
                        'height': str(end_y - start_y),
                        'fill': color
                    })
                    rectangle_elements.append(rect_element)
    for elem in rectangle_elements + line_elements:
        original_svg_root.insert(0, elem)

    # Save the modified result in the SVG file
    original_svg_tree.write(output_file)

    # Call rsvg-convert to convert the SVG file to a different format (e.g., PDF)
    output_format = "pdf"
    input_svg = output_file
    output_file_pdf = os.path.splitext(output_file)[0] + f'.{output_format}'
    subprocess.run(["rsvg-convert", "-f", output_format, "-o", output_file_pdf, input_svg])
    print(f"File converted to {output_format}: {output_file_pdf}")

def parse_svg_filename(filename):
    if not filename.endswith('.svg'):
        filename += '.svg'
    return filename

def main():
    parser = argparse.ArgumentParser(description='Insert colored regions from a BED file into an SVG file.')
    parser.add_argument('-B', choices=['hg38', 'hg37'], default='hg38', help='Input build "hg38" or "hg37"')
    parser.add_argument('-I', type=str, required=True, help='Input BED file')
    parser.add_argument('-O', type=str, required=True, help='Output SVG file')
    
    # Colori di default per le ancestries
    default_ancestry_colors = [
        ("ancestry0", "#b35806"),
        ("ancestry1", "#b2182b"),
        ("ancestry2", "#fdb863"),
        ("ancestry3", "#fddbc7"),
        ("ancestry4", "#f4a582"),
        ("ancestry5", "#d6604d"),
        ("ancestry6", "#4393c3"),
        ("ancestry7", "#2166ac"),
        ("ancestry8", "#e08214"),
        ("ancestry9", "#92c5de"),
        ("ancestry10", "#ffffff"),
    ]

    # Aggiungiamo gli argomenti --ancestry come opzionali
    for i, (ancestry_name, default_color) in enumerate(default_ancestry_colors):
        parser.add_argument(f'--{ancestry_name}', nargs='*', help=f'Name and color for {ancestry_name}, default color {default_color}')

    args = parser.parse_args()
    
        # Determina le lunghezze dei cromosomi e il file SVG in base alla build specificata
    if args.B == 'hg37':
        chromosome_lengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
        svg_file = 'hg37.svg'
    else:  # default to hg38 if not specified or specified as hg38
        chromosome_lengths = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
        svg_file = 'hg38.svg'


    # Mappa degli ancestries con colori di default se non specificati
# Mappa degli ancestries basata sugli input dell'utente
    ancestries = {}
    for i, (ancestry_name, default_color) in enumerate(default_ancestry_colors):
        ancestry_arg = getattr(args, ancestry_name, None)
        if ancestry_arg:
            # Assumi che il primo elemento sia il nome e il secondo (se presente) sia il colore
            name = ancestry_arg[0] if ancestry_arg else ancestry_name  # Il nome dell'ancestry è obbligatorio
            color = ancestry_arg[1] if len(ancestry_arg) > 1 else default_color  # Usa il colore di default se non specificato
            ancestries[name] = color  # Usa il nome fornito come chiave

    
    individual_name = os.path.basename(args.I).split('.')[0]

    insert_colored_regions(svg_file, args.I, args.O, individual_name, ancestries, chromosome_lengths)

if __name__ == "__main__":
    main()
