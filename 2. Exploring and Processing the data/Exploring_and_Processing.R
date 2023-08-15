# This script explores the data, allelic proportions

# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(stringr)
library(gridExtra)
library(grid)
library(ggpubr)
library(readxl)
library(ggcorrplot)
library(magick)

# Set directory
setwd("~/Desktop")

# Read in frequency and allele file
SNP_Frequencies <- read.csv("SNP_Frequencies.csv", header=TRUE)
SNP_Alleles <- read.csv("SNP_Alleles.csv", header=TRUE)

# Turn the O's to 0's so the IDs match
SNP_Alleles$X <- gsub("O", "0", SNP_Alleles$X)
rownames(SNP_Alleles) <- SNP_Alleles$X
SNP_Alleles$X <- NULL

##################################################
# Children
##################################################
filtered_children_major_allele_frequencies<-SNP_Alleles %>% 
  filter(grepl("A|B|C", rownames(.)))
# IDs beginning with B or C are the younger siblings.

# Removing any columns that have only "0"s in the entire column
filtered_children_major_allele_frequencies <- filtered_children_major_allele_frequencies %>%
  select(where(~ any(. != 0)))
filtered_children_major_allele_frequencies$X.1 <- NULL

lapply(names(filtered_children_major_allele_frequencies), function(col) {
  # Filter out the "--" values from each column
  filtered_data <- filtered_children_major_allele_frequencies[filtered_children_major_allele_frequencies[[col]] != "--", ]
  
  # Convert the fill variable to a factor
  filtered_data[[col]] <- factor(filtered_data[[col]])
  
  plot<-ggplot(filtered_data, aes(.data[[col]], ..count..)) + 
    geom_bar(aes(fill = .data[[col]]), position = "dodge") +
    scale_fill_manual(values=c("darkorange2", "ivory3", "slateblue")) +
    ggtitle(paste("Major Allele Frequency of", col)) + 
    scale_x_discrete(labels = str_wrap(c('Minor Alleles', 'Heterozygous', 'Major Alleles'), width = 10)) +
    ylab("Percent") +
    theme(legend.position = "none", axis.ticks.x = element_blank())
  
  
}) -> list_plots
# Can use list_plots[[1]] to view each plot one at a time

# Calculate frequencies for each column
child_column_frequencies <- sapply(filtered_children_major_allele_frequencies, function(x) table(factor(x, levels = 0:2))) %>%
  `rownames<-`(c("Minor Alleles", "Heterozygous", "Major Alleles")) %>%
  t() %>%
  as.data.frame()

# Get percentages from the frequencies
total_counts <- rowSums(child_column_frequencies)
child_column_percentages <- child_column_frequencies / total_counts * 100
child_column_percentages <- round(child_column_percentages, 1)

# Creating df that has percentages, rather than frequencies
child_allelic_distribution_percentages <- data.frame()  # Create an empty data frame

for (i in 1:nrow(child_column_frequencies)) {
  row_sum <- sum(child_column_frequencies[i, ])  # Calculate the sum of the row
  normalized_row <- child_column_frequencies[i, ] / row_sum  # Divide each cell by the row sum
  child_allelic_distribution_percentages <- rbind(child_allelic_distribution_percentages, normalized_row)  # Append the normalized row to the new data frame
}

child_column_percentages %>%
  mutate(across(everything(), ~ round(., 1))) %>%
  kbl() %>% 
  kable_classic() %>% 
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(1:4, width = "1em") #%>% #Uncomment to save table 
  #as_image(., width = 7, file="~/Desktop/Allele_Percentage_Distribution_Table.png", zoom=15)


child_allelic_distribution_percentages_long<-data.frame(rsID = rownames(child_allelic_distribution_percentages),
                                                        child_allelic_distribution_percentages,
                                                        row.names=NULL) %>%
  tidyr::gather(., Allele, Percent, -rsID)

# Convert the 'Allele' column to a factor with the desired order
child_allelic_distribution_percentages_long$Allele <- factor(child_allelic_distribution_percentages_long$Allele, 
                                                             levels = c("Minor.Alleles", "Heterozygous", "Major.Alleles"))

# Plot the grouped barplot
children_grouped_plot<-ggplot(child_allelic_distribution_percentages_long, aes(x = rsID, y = Percent, fill = Allele)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.4, width=0.8) +
  labs(x = "rsID", y = "Percent", fill = "Allele") +
  scale_fill_manual(values = c("darkorange2", "ivory3", "slateblue"),  labels = c("Minor Alleles", "Heterozygous", "Major Alleles")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "right",
        
        plot.title = element_text(size=15, face="bold", hjust= 0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.99)) +
  ggtitle("Allelic Proportion Distribution")


children_grouped_plot
# Specify the output file path and dimensions, then save plot
output_file <- "~/Desktop/Allele_Proportion_Distribution.png"

ggsave(filename = output_file, plot=children_grouped_plot, width = 8, height = 6, dpi = 900)

##################################################
# PARENTS
##################################################
filtered_parent_major_allele_frequencies<- SNP_Alleles %>% 
  filter(grepl("P|S", rownames(.)))

# Removing any columns that have only "0"s in the entire column
filtered_parent_major_allele_frequencies <- filtered_parent_major_allele_frequencies %>%
  select(where(~ any(. != 0)))
filtered_children_major_allele_frequencies$X.1 <- NULL

lapply(names(filtered_parent_major_allele_frequencies), function(col) {
  # Filter out the "--" values from each column
  filtered_parent_data <- filtered_parent_major_allele_frequencies[filtered_parent_major_allele_frequencies[[col]] != "--", ]
  
  # Convert the fill variable to a factor
  filtered_parent_data[[col]] <- factor(filtered_parent_data[[col]])
  
  plot<-ggplot(filtered_parent_data, aes(.data[[col]], ..count..)) + 
    geom_bar(aes(fill = .data[[col]]), position = "dodge") +
    scale_fill_manual(values=c("darkorange2", "ivory3", "slateblue")) +
    ggtitle(paste("Major Allele Frequency of", col)) + 
    scale_x_discrete(labels = str_wrap(c("Minor Alleles", "Heterozygous", "Major Alleles"), width = 10)) +
    ylab("Percent") +
    theme(legend.position = "none", axis.ticks.x = element_blank())
  
}) -> list_plots
# Can use list_plots[[1]] to view each plot one at a time


# # Calculate frequencies for each column
parent_column_frequencies <- sapply(filtered_parent_major_allele_frequencies, function(x) table(factor(x, levels = 0:2))) %>%
  `rownames<-`(c("Minor Alleles", "Heterozygous", "Major Alleles")) %>%
  t() %>%
  as.data.frame() 

p_total_counts <- rowSums(parent_column_frequencies)
parent_column_percentages <- parent_column_frequencies / p_total_counts * 100
parent_column_percentages <- round(parent_column_percentages, 1)
# Creating a df that has percentages, rather than frequencies
parent_allelic_distribution_percentages <- data.frame()  # Create an empty data frame

for (i in 1:nrow(parent_column_frequencies)) {
  row_sum <- sum(parent_column_frequencies[i, ])  # Calculate the sum of the row
  normalized_row <- parent_column_frequencies[i, ] / row_sum  # Divide each cell by the row sum
  parent_allelic_distribution_percentages <- rbind(parent_allelic_distribution_percentages, normalized_row)  # Append the normalized row to the new data frame
}

parent_column_percentages %>%
  mutate(across(everything(), ~ round(., 3))) %>%
  kbl() %>% 
  kable_classic() %>% 
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(1:4, width = "1em") #%>%
  #as_image(., width = 7, file="~/Desktop/P_Allele_Percentage_Distribution_Table.png", zoom=15)

parent_allelic_distribution_percentages_long<-data.frame(rsID = rownames(parent_allelic_distribution_percentages),
                                                         parent_allelic_distribution_percentages,
                                                         row.names=NULL) %>%
  tidyr::gather(., Allele, Percent, -rsID) 

# Convert the 'category' column to a factor with the desired order
parent_allelic_distribution_percentages_long$Allele <- factor(parent_allelic_distribution_percentages_long$Allele, 
                                                              levels = c("Minor.Alleles", "Heterozygous", "Major.Alleles"))

# Plot the grouped barplot
parent_grouped_plot<-ggplot(parent_allelic_distribution_percentages_long, aes(x = rsID, y = Percent, fill = Allele)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.4, width=0.8) +
  labs(x = "rsID", y = "Percent", fill = "Allele") +
  scale_fill_manual(values = c("darkorange2", "ivory3", "slateblue"),  labels = c("Minor Alleles", "Heterozygous", "Major Alleles")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "right",
        
        plot.title = element_text(size=15, face="bold", hjust= 0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.99)) +
  ggtitle("Allelic Proportion Distribution")

parent_grouped_plot

# Specify the output file path and dimensions, then save plot
output_file <- "~/Desktop/P_Allele_Proportion_Distribution.png"

ggsave(filename = output_file, plot=parent_grouped_plot, width = 8, height = 6, dpi = 900)

##################################################
# Combining tables
##################################################
parent_column_percentages
child_column_percentages
combined_proportion_df<-cbind(child_column_percentages,
                              parent_column_percentages)

combined_proportion_df %>%
  `colnames<-`(c("Minor Alleles", "Heterozygous", "Major Alleles","Minor Alleles", "Heterozygous", "Major Alleles")) %>%
  kable(format = "html",escape = FALSE, caption = "Allele Proportion Within Group")%>% 
  kable_classic() %>%
  kable_styling(full_width = TRUE) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1:7, width = "1em") %>% 
  column_spec(1, bold = TRUE) %>% 
  column_spec (4, border_right = T) %>% 
  add_header_above(., c(" ", "Children" = 3, "Parent" = 3)) %>% 
  as_image(., width = 7, file="~/Desktop/All_allelic_proportion_table.png", zoom=15) 

##################################################
# Combine the grouped plots
##################################################

children_grouped_plot_1 <- children_grouped_plot + theme(plot.title = element_blank()) +
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0, vjust = 1, 
                                  margin = margin(0, 0, 10, 0))) +
  theme(plot.background = element_blank()) # removes line between plots+


parent_grouped_plot_2 <- parent_grouped_plot + theme(plot.title = element_blank()) +
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0, vjust = 1, 
                                  margin = margin(0, 0, 10, 0))) +
  theme(plot.background = element_blank()) # removes line between plots

g <- ggarrange(children_grouped_plot_1, parent_grouped_plot_2, 
               labels = c("Children", "Parents"),
               ncol = 2, nrow = 1)+
  theme(plot.background = element_rect(fill = "white", 
                                       color = NA))  # Set color = NA to remove the border line
# Create a custom title plot
title_plot <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.95, label = "Allelic Proportion Distribution"),
            hjust = 0.5, vjust = 0, size = 10) +
  theme_void()

# Combine the title plot and the grouped plot vertically
final_plot <- grid.arrange(title_plot, g, ncol = 1, heights = c(0.1, 1))


ggsave(filename = "~/Desktop/All_allelic_proportion.png", final_plot,dpi = 900, width = 11, height = 8)
