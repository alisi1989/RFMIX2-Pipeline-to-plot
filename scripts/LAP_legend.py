#Autors Alessandro Lisi & Micahel C. Campbell
#be sure to priorly install Matplotlib (pip install Matplotlib)

#use python LAP.py --help to see the option

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Create a custom legend with colored rectangles. Color code must be inside quotations mark")
    parser.add_argument("--anc0", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc0 Africa ''#0000ff'' or others country and color code")
    parser.add_argument("--anc1", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc1 Europe  ''#850b39'' magenta or others country and color code")
    parser.add_argument("--anc2", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc2 MiddleEast  ''#F4A500'' or others country and color code")
    parser.add_argument("--anc3", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc3 Asia  ''#292518'' or others country and color code")
    parser.add_argument("--anc4", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc4 America  ''#076940'' or others country and color code")
    parser.add_argument("--anc5", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc4 Country1  ''#00000'' or others country and color code")
    parser.add_argument("--anc6", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc4 Country2  ''#00000'' or others country and color code")
    parser.add_argument("--anc7", nargs=2, metavar=("NAME", "COLOR"), help="Example --anc4 Country3  ''#00000'' or others country and color code")
    parser.add_argument("--output", metavar=("FILENAME"), help="Specify the output file name")
    parser.add_argument("--format", metavar=("FORMAT"), help="Specify the output file format (e.g., svg, png, pdf). Default is svg.")

    args = parser.parse_args()
    ancestry_info = {}

    for i in range(8):
        arg_name = f"anc{i}"
        if hasattr(args, arg_name) and getattr(args, arg_name):
            ancestry_info[getattr(args, arg_name)[0]] = getattr(args, arg_name)[1]

    fig, ax = plt.subplots(figsize=(2, len(ancestry_info) * 0.7))
    ax.set_axis_off()

    legend_handles = []
    legend_height = 4
    font_size = 42

    for name, color in ancestry_info.items():
        rect = patches.Rectangle((0, i * 0.7), 1, legend_height, linewidth=0, edgecolor='none', facecolor=color, label=name)
        legend_handles.append(rect)

    legend = ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1, 1), prop={'size': font_size})

    for spine in ax.spines.values():
        spine.set_visible(False)

    legend.set_bbox_to_anchor((1.05, 1))

    for text in legend.get_texts():
        text.set_fontsize(font_size)

    output_format = args.format or "svg"
    output_filename = args.output or "Ancestry_Legends"
    
    # Estrai l'estensione dal nome del file di output
    _, ext = os.path.splitext(output_filename)
    
    # Aggiungi l'estensione solo se non è stata specificata dall'utente
    if not ext:
        output_filename += "." + output_format
    
    plt.savefig(output_filename, format=output_format, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

if __name__ == "__main__":
    main()
