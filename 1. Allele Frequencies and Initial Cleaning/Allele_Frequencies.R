# This script filters our initial allele data frame

# Load packages
library(tidyverse)
library(stringr)

# Set directory
setwd("~/Desktop")

#############################################
# FILTER OUT DUPLICATE SNPS
############################################

# load data
SNP_data <- read.csv("SNP_numbers.csv", header = TRUE)

# Set the names as rownames
row.names(SNP_data) <- SNP_data$X
SNP <- SNP_data[-1]

# Remove duplicate rows with less data present
SNP <- SNP[-4]
SNP <- SNP[-7]
SNP <- SNP[-7]
SNP <- SNP[-8]
SNP <- SNP[-11]
SNP <- SNP[-12]
SNP <- SNP[-13]
SNP <- SNP[-13]
SNP <- SNP[-15]
SNP <- SNP[-16]
SNP <- SNP[-17]
SNP <- SNP[-17]
SNP <- SNP[-18]

df <- SNP

#############################################
# CONVERT MINOR, HETEROZYGOUS, AND MAJOR TO 0,1,2
############################################

# Create data frame showing major and minor alleles for each SNP
ID <- colnames(df)
Maj <- c("C","T","C","A","A","G","T","C","T","G","T","T","T","G","T","T","T")
Min <- c("T","C","T","C","G","A","C","T","C","A","A","G","A","A","C","C","C")

Maj_Min_alleles <- data.frame(ID, Maj, Min)

# Find indices of matching SNPs in Maj_Min_alleles
match_indices <- match(ID, Maj_Min_alleles$ID)

# Extract major alleles from matching SNPs
Major <- Maj_Min_alleles$Maj[match_indices]

# Create data frame with SNP names and major alleles
maj_df <- data.frame(ID, Major)

# Create a matrix to store major allele frequencies
allele_freq <- matrix(ncol = ncol(df), nrow = nrow(df))

# Iterate through each column of the SNP dataframe
for (i in 1:ncol(df)) {
  snp <- maj_df[i, 1]  
  major_allele <- as.character(maj_df[i, 2])  # Major allele for the SNP
  
  # Loop through each row of the allele_freq data frame
  for (j in 1:length(allele_freq[,i])) {
    # Check if there is missing genotype data for the SNP in the individual
    if (df[j, i] == "--") {
      allele_freq[j, i] <- "--" # Preserve missing genotype data
    } else {
      # Calculate the count of major alleles in the individuals genotype
      allele_freq[j, i] <- str_count(df[j, i], major_allele)
    }
  }
}

# Convert allele_freq to a data frame and set row and column names
allele_freq <- as.data.frame(allele_freq)
rownames(allele_freq) <- rownames(df)
colnames(allele_freq) <- c("rs1801278", "rs35874116", "rs5400", "rs5417", "rs5418", "rs12970134", "rs1421085", "rs2943641","rs7102710","rs17683430", "rs9939609", "rs17817449","rs1558902","rs17300539","rs17782313","rs4994","rs1801260")

# Convert NAs to 1 (mean)
allele_freq[is.na(allele_freq)] <- 1 

# Calculate frequencies for each column
column_frequencies <- sapply(allele_freq, function(x) table(factor(x, levels = 0:2)))

# Print the frequencies
print(column_frequencies)


#############################################
# SAVES
############################################

write.csv(column_frequencies, file = "SNP_Frequencies.csv", row.names = TRUE)
write.csv(allele_freq, file = "SNP_Alleles.csv", row.names = TRUE)

