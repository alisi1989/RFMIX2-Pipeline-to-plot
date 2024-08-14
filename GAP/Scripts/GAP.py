import argparse
import pandas as pd
import matplotlib.pyplot as plt
import re

def plot_admixture(df, output_file):
   # Extract individual names from first column (ignoring header)
    individuals = df.iloc[:, 0].tolist()

  # Extract ancestral population data
    ancestries = df.columns[1:]

  # Extract the proportions for each ancestral population
    proportions = df.iloc[:, 1:].map(lambda x: float(re.match(r'^([0-9.]+)', x).group(1))).values

  # Initialize a dictionary for ancestry-color matches
    ancestry_colors = {}

    for ancestry in ancestries:
       # Extract color from column and remove non-numeric characters
        color = df[ancestry].str.extract(r'#([0-9a-fA-F]+)').iloc[0, 0]
        ancestry_colors[ancestry] = f'#{color}'

   # Create a bar graph for each individual
    fig, ax = plt.subplots(figsize=(15, 6))
    bottom = [0] * len(individuals)

    for i, ancestry in enumerate(ancestries):
        proportion = proportions[:, i]
        ax.bar(individuals, proportion, bottom=bottom, color=ancestry_colors[ancestry], label=ancestry)
        bottom = [sum(x) for x in zip(bottom, proportion)]

    ax.set_ylabel('Global Ancestry Proportion')
    ax.set_title('Global Ancestry Admixture')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=len(ancestries))  # Move the legend to the top right
    plt.xticks(rotation=90, fontsize=10)
    
    # Add a custom legend above the chart
    handles, labels = ax.get_legend_handles_labels()
    legend_labels = [label for label in labels]
    plt.legend(handles, legend_labels, loc='upper left', bbox_to_anchor=(1, 1), ncol=1, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf')
    plt.close()

def main():
# Define command line arguments
    parser = argparse.ArgumentParser(description="Create an Admixture chart from ancestral population data")
    parser.add_argument("--input", help="Input File .bed", required=True)
    parser.add_argument("--output", help="Name of the graph output file", required=True)
    args = parser.parse_args()

   # Load the .bed input file using pandas
    df = pd.read_csv(args.input, sep='\t')

 # Create Admixture plot
    plot_admixture(df, args.output)

if __name__ == "__main__":
    main()
