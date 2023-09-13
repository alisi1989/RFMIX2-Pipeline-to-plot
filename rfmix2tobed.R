# Read the file
msp <- read.table("~/RFMIX2_pipeline_master/Output/Mozabite1_ind1_allchr.msp.tsv", header = FALSE, sep = "\t")

# Select only the first 8 columns
msp2 <- msp[, 1:8]

# Remove columns 4, 5, and 6
msp2 <- msp2[, -4]
msp2 <- msp2[, -4]
msp2 <- msp2[, -4]

# Rename the columns
names(msp2) <- c("chm", "spos", "epos", "ind1", "ind1")

# Create 2 dataframes
ch1 <- msp2[, -5]
ch2 <- msp2[, -4]

# Replace values in the "ind1" column with corresponding ancestry values
ch1$ind1 <- ifelse(ch1$ind1 == 0, "ANC0",
                   ifelse(ch1$ind1 == 1, "ANC1",
                   ifelse(ch1$ind1 == 2, "ANC2",
                   ifelse(ch1$ind1 == 3, "ANC3",
                   ifelse(ch1$ind1 == 4, "ANC4", ch1$ind1)))))
                   
# Replace values in the "ind1" column with corresponding ancestry values
ch2$ind1 <- ifelse(ch2$ind1 == 0, "ANC0",
                   ifelse(ch2$ind1 == 1, "ANC1",
                   ifelse(ch2$ind1 == 2, "ANC2",
                   ifelse(ch2$ind1 == 3, "ANC3",
                   ifelse(ch2$ind1 == 4, "ANC4", ch2$ind1)))))

# Save the two files for each individual
write.table(ch2, "~RFMIX2_pipeline_master/Output/Mozabite1_ind1_hap2.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(ch1, "~/RFMIX2_pipeline_master/Output/Mozabite1_ind1_hap1.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
