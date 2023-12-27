#!/bin/bash

# Output file name
output_file="Mozabite2_allchr.rfmix.Q"

# Initialize counters for the average
sum_col2=0
sum_col3=0
sum_col4=0
count=0

# Loop through chromosome files
for i in {1..22}; do
    file="Mozabite2_chr${i}.rfmix.Q"

    # Extract the line with Mozabite2 and calculate the average
    if grep -q "Mozabite1" "$file"; then
        value_col2=$(awk '$1 == "Mozabite1" {print $2}' "$file")
        value_col3=$(awk '$1 == "Mozabite1" {print $3}' "$file")
        value_col4=$(awk '$1 == "Mozabite1" {print $4}' "$file")

        sum_col2=$(awk "BEGIN {print $sum_col2 + $value_col2}")
        sum_col3=$(awk "BEGIN {print $sum_col3 + $value_col3}")
        sum_col4=$(awk "BEGIN {print $sum_col4 + $value_col4}")
        count=$((count + 1))
    fi
done

# Calculate the average
avg_col2=$(awk "BEGIN {print $sum_col2 / $count}")
avg_col3=$(awk "BEGIN {print $sum_col3 / $count}")
avg_col4=$(awk "BEGIN {print $sum_col4 / $count}")

# Create the output file with the average
echo -e "#rfmix diploid global ancestry .Q format output\n#sample\tAfrica\tEurope\tMiddle_East" > "$output_file"
echo -e "Mozabite2\t$avg_col2\t$avg_col3\t$avg_col4" >> "$output_file"
