import argparse
import pandas as pd

# Function to add colors to ancestral populations
def add_colors(df, ancestry_colors):
    for ancestry, color in ancestry_colors.items():
        if ancestry in df.columns:
            df[ancestry] = df[ancestry].apply(lambda x: f"{x:.5f} ({color})")
    return df

def main():
   # Define command line arguments
    parser = argparse.ArgumentParser(description="Add colors to ancestral populations in a .bed file. Each ancestry can have a default color or a user-specified color.")
    
   # List of predefined colors
    default_colors = ["#a32e2e", "#0a0ae0", "#f4a500", "#0d6601", "#22ba9d", "#839dfc", "#9a5dc1", "#26962b", "#707070", "#00cfff", "#790ee0"]
    
    parser.add_argument("--input", help="Input File .bed", required=True)
    parser.add_argument("--output", help="Name of the output file .bed", required=True)
    
   # Generate --ancestry arguments dynamically
    for i in range(11):
        parser.add_argument(f"--ancestry{i}", nargs='*', help=f"Name for ancestry{i} and optionally its color. If no color is specified, a default color will be assigned.")

    args = parser.parse_args()

   # Load the .bed input file using pandas, skipping the first line
    df = pd.read_csv(args.input, sep='\t', skiprows=1)
    
    # Get the list of ancestry columns
    ancestry_columns = df.columns[1:]

    # Create a dictionary with ancestry-color matches
    ancestry_colors = {}
    for i, ancestry_name in enumerate(ancestry_columns):
        ancestry_arg = getattr(args, f"ancestry{i}", None)
        if ancestry_arg:
            ancestry_name = ancestry_arg[0]
            if len(ancestry_arg) > 1:
               # If the user specified a color
                ancestry_color = ancestry_arg[1]
            else:
               # Otherwise, assign a default color
                ancestry_color = default_colors[i % len(default_colors)]  # Usa i colori predefiniti in loop
        else:
         # Use default name and default color
            ancestry_color = default_colors[i % len(default_colors)]  # Use default colors in loop
        ancestry_colors[ancestry_name] = ancestry_color

    # Add colors to ancestral populations
    df = add_colors(df, ancestry_colors)

 # Save the resulting DataFrame to the .bed output file
 # Add the original comment lines except the first line
    with open(args.output, 'w') as output:
        with open(args.input, 'r') as input_file:
            lines = input_file.readlines()
         #output.write(lines[1]) # Write the second header line
        df.to_csv(output, sep='\t', index=False)

    print(f"File {args.output} successfully created.")

if __name__ == "__main__":
    main()
