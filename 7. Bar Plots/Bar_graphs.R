# Create Barplots to convey data in a more universally understable way
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)

# Set directory
setwd("~/Desktop")

children_final_models<-read.csv("children_intersection_refit_models.csv",check.names = FALSE,row.names = 1)
children_combined_unscaled<-read.csv("children_anthropometric_data_and_allele_freq_combined.csv", check.names = FALSE, row.names = 1)
parent_final_models<-read.csv("parent_intersection_refit_models.csv",check.names = FALSE,row.names = 1)
parent_combined_unscaled<-read.csv("parent_anthropometric_data_and_allele_freq_combined.csv", check.names = FALSE, row.names = 1)

# Remove . that had replaced spaces
children_final_models$Measure <- gsub("\\.", " ", children_final_models$Measure)
children_final_models$Measure <- gsub("BMI Z Score", "BMI Z-Score", children_final_models$Measure)
parent_final_models$Measure <- gsub("\\.", " ", parent_final_models$Measure)

GEE_final_models<-read.csv("gee_coefficients_long.csv",check.names = FALSE,row.names = 1)
GEE_final_models <- GEE_final_models %>% rename(Measure = colname)
GEE_final_models <- GEE_final_models %>% rename(SNP = rowname)

# Combine parent and children data
GEE_Parent <- parent_combined_unscaled[,-2] #Remove BMI column
GEE_Child <- children_combined_unscaled[, -2] #Remove BMI-Z column
GEE_combined_unscaled<-rbind(GEE_Parent, GEE_Child)

##################################################
# Add alleles and their frequency to each table
##################################################

# Creating a function where, when inputting the rsID, it returns the gene name
find_allele_from_SNP <- function(rsID) {
  matching_rows <- SugarGenes_SNPs$SNP == rsID
  
  # Check if any matching rows are found
  if (any(matching_rows)) {
    # Extract the corresponding Alleles
    Alleles <- SugarGenes_SNPs$Major.MinorAlleles[matching_rows]
  } else {
    Alleles <- "No matching alleles found"
  }
  
  return(Alleles)
}

SugarGenes_SNPs<-read.csv("Sugar_locations.csv")
SugarGenes_SNPs$SNP<-tolower(SugarGenes_SNPs$SNP)


for (i in 1:nrow(children_final_models)){
  children_final_models$Alleles[i]<-find_allele_from_SNP(children_final_models$SNP[i])
}
for (i in 1:nrow(parent_final_models)){
  parent_final_models$Alleles[i]<-find_allele_from_SNP(parent_final_models$SNP[i])
}
for (i in 1:nrow(GEE_final_models)){
  GEE_final_models$Alleles[i]<-find_allele_from_SNP(GEE_final_models$SNP[i])
}


# Split the alleles into a major and minor column
children_final_models <- separate(children_final_models, col = Alleles, into = c("major", "minor"), sep = "/")
parent_final_models <- separate(parent_final_models, col = Alleles, into = c("major", "minor"), sep = "/")
GEE_final_models <- separate(GEE_final_models, col = Alleles, into = c("major", "minor"), sep = "/")

##################################################
# Add Children frequency to the table (manually)
##################################################

# Columns 1 to 17 that you want to get value counts for
columns_to_count <- 11:27

# Initialize an empty list to store the value count tables for each column
value_count_tables <- list()

# Loop through each column and apply the table() function to get value counts
for (col_index in columns_to_count) {
  col_name <- colnames(children_combined_unscaled)[col_index]
  value_count_tables[[col_name]] <- table(children_combined_unscaled[[col_index]])
}
# Make individual tables for each column..... probably dont need to do this if you don't have alleles with 0
# Only extract tables that have been selected to make your life a bit easier
children_final_models$SNP
names(value_count_tables)

child_rs5400 <- value_count_tables[[3]] # 0 = 5 so you need to collapse it with 2
child_rs1801278 <- value_count_tables[[1]] # No 0
child_rs17683430 <- value_count_tables[[10]] # No 0

# collapse rs5400 manually
child_rs5400
child_rs5400 <- c(0, 50, 159)

child_rs5400 <- as.data.frame(child_rs5400)
rownames(child_rs5400) <- c("min","het","maj")
child_rs5400 <- t(child_rs5400)
child_rs5400 <- child_rs5400 / sum(child_rs5400) * 100 # This gives us percentage instead of count
# Collapse this manually as well
child_rs1801278
child_rs1801278 <- c(0, 36, 173)

child_rs1801278 <- as.data.frame(child_rs1801278)
rownames(child_rs1801278) <- c("min","het","maj")
child_rs1801278 <- t(child_rs1801278)
child_rs1801278 <- child_rs1801278 / sum(child_rs1801278) * 100

# Last one
child_rs17683430
child_rs17683430 <- c(0, 17, 192)

child_rs17683430 <- as.data.frame(child_rs17683430)
rownames(child_rs17683430) <- c("min","het","maj")
child_rs17683430 <- t(child_rs17683430)
child_rs17683430 <- child_rs17683430 / sum(child_rs17683430) * 100

# Combine all the allele frequencies
child_frequencies <- rbind(child_rs5400, child_rs1801278, child_rs17683430)
rownames(child_frequencies) <- c("rs5400", "rs1801278", "rs17683430")
child_frequencies <- as.data.frame(child_frequencies)


# Add frequencies to their matching rsID in the final_models table
child_frequencies
children_final_models

# Function to add frequencies to table
update_frequencies <- function(rsID) {
  matching_rows <- rownames(child_frequencies) == rsID
  
  # Check if any matching rows are found
  if (any(matching_rows)) {
    # Extract the corresponding frequencies
    minor_freq <- round(child_frequencies$min[matching_rows], 1)
    het_freq <- round(child_frequencies$het[matching_rows], 1)
    major_freq <- round(child_frequencies$maj[matching_rows], 1)
    
    # Return a list containing the frequencies
    return(list(minor_freq = minor_freq, het_freq = het_freq, major_freq = major_freq))
  } else {
    print("No matching frequencies found")
    return(NULL)
  }
}

for (i in 1:nrow(children_final_models)) {
  snp_id <- children_final_models$SNP[i]
  freqs <- update_frequencies(snp_id)
  
  if (!is.null(freqs)) {
    # Update the frequencies in children_final_models
    children_final_models$minor_freq[i] <- freqs$minor_freq
    children_final_models$het_freq[i] <- freqs$het_freq
    children_final_models$major_freq[i] <- freqs$major_freq
  } else {
    # Handle the case where no matching frequencies are found
    children_final_models$minor_freq[i] <- NA
    children_final_models$het_freq[i] <- NA
    children_final_models$major_freq[i] <- NA
  }
}


##################################################
# Add Parent frequency to the table (manually)
##################################################

# Columns 1 to 17 that you want to get value counts for
columns_to_count <- 11:27

# Initialize an empty list to store the value count tables for each column
p_value_count_table <- list()

# Loop through each column and apply the table() function to get value counts
for (col_index in columns_to_count) {
  col_name <- colnames(parent_combined_unscaled)[col_index]
  p_value_count_table[[col_name]] <- table(parent_combined_unscaled[[col_index]])
}
# Make individual tables for each column..... probably dont need to do this if you don't have alleles with 0
# Only extract tables that have been selected to make your life a bit easier
parent_final_models$SNP
names(p_value_count_table)

parent_rs12970134 <- p_value_count_table[[6]] # This one is fine
parent_rs7102710 <- p_value_count_table[[9]] # No 0
parent_rs35874116 <- p_value_count_table[[2]] # This one is fine
parent_rs1801260 <- p_value_count_table[[17]] # This one is fine

# Add 0
parent_rs7102710
parent_rs7102710 <- c(0, 30, 231)

parent_rs7102710 <- as.data.frame(parent_rs7102710)
rownames(parent_rs7102710) <- c("min","het","maj")
parent_rs7102710 <- t(parent_rs7102710)
parent_rs7102710 <- parent_rs7102710 / sum(parent_rs7102710) * 100 # This gives us percentage instead of count

# percentages for the two good ones
parent_rs12970134 <- parent_rs12970134 / sum(parent_rs12970134) * 100
parent_rs35874116 <- parent_rs35874116 / sum(parent_rs35874116) * 100
parent_rs1801260 <- parent_rs1801260 / sum(parent_rs1801260) * 100

# Ok now combine all of them

parent_frequencies <- rbind(parent_rs7102710, parent_rs12970134, parent_rs35874116, parent_rs1801260)
rownames(parent_frequencies) <- c("rs7102710", "rs12970134", "rs35874116", "rs1801260")
parent_frequencies <- as.data.frame(parent_frequencies)


# Now, add these to their matching rsID in the final_models table
parent_frequencies
parent_final_models

# Repeat same process done on children for parents
p_update_frequencies <- function(rsID) {
  matching_rows <- rownames(parent_frequencies) == rsID
  
  # Check if any matching rows are found
  if (any(matching_rows)) {
    # Extract the corresponding frequencies
    minor_freq <- round(parent_frequencies$min[matching_rows], 1)
    het_freq <- round(parent_frequencies$het[matching_rows], 1)
    major_freq <- round(parent_frequencies$maj[matching_rows], 1)
    
    # Return a list containing the frequencies
    return(list(minor_freq = minor_freq, het_freq = het_freq, major_freq = major_freq))
  } else {
    print("No matching frequencies found")
    return(NULL)
  }
}

for (i in 1:nrow(parent_final_models)) {
  snp_id <- parent_final_models$SNP[i]
  freqs <- p_update_frequencies(snp_id)
  
  if (!is.null(freqs)) {
    # Update the frequencies in children_final_models
    parent_final_models$minor_freq[i] <- freqs$minor_freq
    parent_final_models$het_freq[i] <- freqs$het_freq
    parent_final_models$major_freq[i] <- freqs$major_freq
  } else {
    # Handle the case where no matching frequencies are found
    parent_final_models$minor_freq[i] <- NA
    parent_final_models$het_freq[i] <- NA
    parent_final_models$major_freq[i] <- NA
  }
}

parent_final_models


##################################################
# Add GEE frequency to the table (manually)
##################################################

# Columns 1 to 17 that you want to get value counts for
columns_to_count <- 10:26

# Initialize an empty list to store the value count tables for each column
gee_value_count_table <- list()

# Loop through each column and apply the table() function to get value counts
for (col_index in columns_to_count) {
  col_name <- colnames(GEE_combined_unscaled)[col_index]
  gee_value_count_table[[col_name]] <- table(GEE_combined_unscaled[[col_index]])
}
# Make individual tables for each column..... probably dont need to do this if you don't have alleles with 0
# Only extract tables that have been selected
GEE_final_models$SNP
names(gee_value_count_table)

gee_rs5417 <- p_value_count_table[[4]] # This one is fine
gee_rs5418 <- p_value_count_table[[5]] # This one is fine
gee_rs9939609 <- p_value_count_table[[11]] # This one is fine
gee_rs17817449 <- p_value_count_table[[12]] # This one is fine
gee_rs1558902 <- p_value_count_table[[13]] # This one is fine
gee_rs1421085 <- p_value_count_table[[7]] # This one is fine

gee_rs5417
gee_rs5417 <- c(68, 214, 188)
gee_rs5417 <- as.data.frame(gee_rs5417)
rownames(gee_rs5417) <- c("min","het","maj")
gee_rs5417 <- t(gee_rs5417)

# percentages for the the good ones
gee_rs5417 <- gee_rs5417 / sum(gee_rs5417) * 100
gee_rs5418 <- gee_rs5418 / sum(gee_rs5418) * 100
gee_rs9939609 <- gee_rs9939609 / sum(gee_rs9939609) * 100
gee_rs17817449 <- gee_rs17817449 / sum(gee_rs17817449) * 100
gee_rs1558902 <- gee_rs1558902 / sum(gee_rs1558902) * 100
gee_rs1421085 <- gee_rs1421085 / sum(gee_rs1421085) * 100

# Ok now combine all of them

gee_frequencies <- rbind(gee_rs5417, gee_rs5418, gee_rs9939609, gee_rs17817449, gee_rs1558902, gee_rs1421085)
rownames(gee_frequencies) <- c("rs5417", "rs5418", "rs9939609", "rs17817449", "rs1558902", "rs1421085")
gee_frequencies <- as.data.frame(gee_frequencies)


# Now, add these to their matching rsID in the final_models table
gee_frequencies
GEE_final_models

# Same process done on children and parents

gee_update_frequencies <- function(rsID) {
  matching_rows <- rownames(gee_frequencies) == rsID
  
  # Check if any matching rows are found
  if (any(matching_rows)) {
    # Extract the corresponding frequencies
    minor_freq <- round(gee_frequencies$min[matching_rows], 1)
    het_freq <- round(gee_frequencies$het[matching_rows], 1)
    major_freq <- round(gee_frequencies$maj[matching_rows], 1)
    
    # Return a list containing the frequencies
    return(list(minor_freq = minor_freq, het_freq = het_freq, major_freq = major_freq))
  } else {
    print("No matching frequencies found")
    return(NULL)
  }
}

for (i in 1:nrow(GEE_final_models)) {
  snp_id <- GEE_final_models$SNP[i]
  freqs <- gee_update_frequencies(snp_id)
  
  if (!is.null(freqs)) {
    # Update the frequencies in children_final_models
    GEE_final_models$minor_freq[i] <- freqs$minor_freq
    GEE_final_models$het_freq[i] <- freqs$het_freq
    GEE_final_models$major_freq[i] <- freqs$major_freq
  } else {
    # Handle the case where no matching frequencies are found
    GEE_final_models$minor_freq[i] <- NA
    GEE_final_models$het_freq[i] <- NA
    GEE_final_models$major_freq[i] <- NA
  }
}

GEE_final_models


##################################################
# Bargraphs of significant things
##################################################

# This Function works for plots with no alleles missing
plot_bargraphs_for_final_models<-function(final_models_df, combined_unscaled_df, folder){
  model_bargraphs<-list()
  model_names<-unique(final_models_df$Measure)
  
  # My input dfs don't have the units of measurement. This new df will allow us to get the label for the y-axis
  y_axis_labels_df <- cbind(
    colnames(parent_combined_unscaled)[!grepl("rs", colnames(parent_combined_unscaled))],
    c("Sex ","BMI","Body Mass (kg)", "Height (cm)", "Waist Circumference (cm)", "Total Sugar (mg)", "Added Sugar (mg)", "Caloric Intake (KCAL)", "Age (Years)",  "Waist Circumference-to-Height Ratio")) %>%
    as.data.frame() %>%
    `colnames<-`(c("measure", "label"))
  
  for (i in 1:length(model_names)){
    max_y_value<-NULL # Used so that all plots for a given measure have the same ylim
    measure<-model_names[i]
    y_label<-y_axis_labels_df$label[y_axis_labels_df$measure == paste0(measure)] 
    individual_model_bargraph_list<-list()
    SNP_IDs<-final_models_df$SNP[final_models_df$Measure == model_names[i]]
    Gene<-final_models_df$Gene[final_models_df$Measure == model_names[i]]
    
    for (SNP in SNP_IDs) {
      df <- combined_unscaled_df[, c(model_names[i], SNP)] %>% as.data.frame()
      mean_value_by_allele <- aggregate(df[, 1] ~ df[, 2], 
                                        data = df, 
                                        FUN = mean, 
                                        na.rm = TRUE) %>%
        `colnames<-`(c("Allele", "Mean_Measure")) %>%
        mutate(Allele = factor(Allele))
      
      # Generate x_labels based on SNP IDs
      row_index <- match(SNP, final_models_df$SNP) # Find the corresponding row in final_models_df
      major_allele <- final_models_df$major[row_index]
      minor_allele <- final_models_df$minor[row_index]
      minor_freq <- final_models_df$minor_freq[row_index]
      het_freq <- final_models_df$het_freq[row_index]
      major_freq <- final_models_df$major_freq[row_index]
      x_labels <- c(paste0(minor_allele, minor_allele, "\n", minor_freq, "%"), paste0(major_allele, minor_allele, "\n", het_freq, "%"), paste0(minor_allele, major_allele, "\n", major_freq, "%"))
      
      plot <- ggplot(mean_value_by_allele, aes(x = Allele, y = Mean_Measure, fill = Allele)) +
        geom_col() +
        labs(title = paste0(SNP, ": ", Gene),
             x = "SNP Allele",
             y = paste0("Mean ", y_label)) +
        scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
        scale_x_discrete(labels = x_labels) +  # Use the correct x_labels
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "none") +
        geom_text(aes(label = round(Mean_Measure, 2)), nudge_y = max(mean_value_by_allele$Mean_Measure)*0.05)
      
      individual_model_bargraph_list[[SNP]]<-plot
      max_y_value <- max(max_y_value, max(mean_value_by_allele$Mean_Measure)) # Used so that all plots for a given measure have the same ylim; update the maximum y-axis value if the current plot's maximum is higher
      
    }
    
    # Use the maximum y-axis value as the common ylim for all plots in the list
    for (j in 1:length(individual_model_bargraph_list)) {
      individual_model_bargraph_list[[j]] <- individual_model_bargraph_list[[j]] +
        coord_cartesian(ylim = c(0, max_y_value * 1.1))   
    }
    
    # Df to determine how many cols the final figure should have
    fig_col_numbers_df<-rbind(c("1","1"),
                              c("2","2"),
                              c("3","3"),
                              c("4","2"),
                              c("5","3"),
                              c("6","3"),
                              c("7","3"),
                              c("8","4"),
                              c("9","3")
    ) %>%
      as.data.frame() %>%
      `colnames<-` (c("NumFigs","NumCol"))
    
    ncol<-fig_col_numbers_df$NumCol[fig_col_numbers_df$NumFigs == length(individual_model_bargraph_list)] %>% as.numeric() # Using the fig_col_numbers_df to determine how many columns the model-specific figure should have
    
    individual_model_bargraph_figure <- grid.arrange(
      grobs = individual_model_bargraph_list,
      ncol = ncol,
      top = paste0("SNPs in final ",measure, " Model"))
    
    ggsave(filename = paste0(folder, measure,".png"),
           individual_model_bargraph_figure,
           dpi = 1200,
           width = 10, height = 7,
           units = c("in"))
    
    print(individual_model_bargraph_figure)
    model_bargraphs[[measure]]<-individual_model_bargraph_figure
    
  }
  
  
}

plot_bargraphs_for_final_models(children_final_models,children_combined_unscaled, "~/Desktop/")
plot_bargraphs_for_final_models(parent_final_models,parent_combined_unscaled, "~/Desktop/")
plot_bargraphs_for_final_models(GEE_final_models,GEE_combined_unscaled, "~/Desktop/")

##################################################
# Function if minor allele is missing...
##################################################

# This will still run through everything, but it will at least properly model the plots where the minor allele is 0, so make sure to pay attention you are using the correct saved plot
plot_bargraphs_for_no_minor<-function(final_models_df, combined_unscaled_df, folder){
  model_bargraphs<-list()
  model_names<-unique(final_models_df$Measure)
  
  # My input dfs don't have the units of measurement. This new df will allow us to get the label for the y-axis
  y_axis_labels_df <- cbind(
    colnames(parent_combined_unscaled)[!grepl("rs", colnames(parent_combined_unscaled))],
    c("Sex ","BMI","Body Mass (kg)", "Height (cm)", "Waist Circumference (cm)", "Total Sugar (mg)", "Added Sugar (mg)", "Caloric Intake (KCAL)", "Age (Years)",  "Waist Circumference-to-Height Ratio")) %>%
    as.data.frame() %>%
    `colnames<-`(c("measure", "label"))
  
  for (i in 1:length(model_names)){
    max_y_value<-NULL # Used so that all plots for a given measure have the same ylim
    measure<-model_names[i]
    y_label<-y_axis_labels_df$label[y_axis_labels_df$measure == paste0(measure)] 
    individual_model_bargraph_list<-list()
    SNP_IDs<-final_models_df$SNP[final_models_df$Measure == model_names[i]]
    Gene<-final_models_df$Gene[final_models_df$Measure == model_names[i]]
    
    for (SNP in SNP_IDs) {
      df <- combined_unscaled_df[, c(model_names[i], SNP)] %>% as.data.frame()
      mean_value_by_allele <- aggregate(df[, 1] ~ df[, 2], 
                                        data = df, 
                                        FUN = mean, 
                                        na.rm = TRUE) %>%
        `colnames<-`(c("Allele", "Mean_Measure")) %>%
        mutate(Allele = factor(Allele))
      
      # Generate x_labels based on SNP IDs
      row_index <- match(SNP, final_models_df$SNP) # Find the corresponding row in final_models_df
      major_allele <- final_models_df$major[row_index]
      minor_allele <- final_models_df$minor[row_index]
      het_freq <- final_models_df$het_freq[row_index]
      major_freq <- final_models_df$major_freq[row_index]
      x_labels <- c(paste0(major_allele, minor_allele, "\n", het_freq, "%"), paste0(minor_allele, major_allele, "\n", major_freq, "%"), paste0("Use other function \n for this plot"))
      
      plot <- ggplot(mean_value_by_allele, aes(x = Allele, y = Mean_Measure, fill = Allele)) +
        geom_col() +
        labs(title = paste0(SNP, ": ", Gene),
             x = "SNP Allele",
             y = paste0("Mean ", y_label)) +
        scale_fill_manual(values = c("#56B4E9", "#009E73", "white")) +
        scale_x_discrete(labels = x_labels) +  # Use the correct x_labels
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "none") +
        geom_text(aes(label = round(Mean_Measure, 2)), nudge_y = max(mean_value_by_allele$Mean_Measure)*0.05)
      # coord_cartesian(ylim = c(0, max(mean_value_by_allele$Mean_Measure) * 1.1))
      
      individual_model_bargraph_list[[SNP]]<-plot
      max_y_value <- max(max_y_value, max(mean_value_by_allele$Mean_Measure)) # Used so that all plots for a given measure have the same ylim; update the maximum y-axis value if the current plot's maximum is higher
      
    }
    
    # Use the maximum y-axis value as the common ylim for all plots in the list
    for (j in 1:length(individual_model_bargraph_list)) {
      individual_model_bargraph_list[[j]] <- individual_model_bargraph_list[[j]] +
        coord_cartesian(ylim = c(0, max_y_value * 1.1))   
    }
    
    # Df to determine how many cols the final figure should have
    fig_col_numbers_df<-rbind(c("1","1"),
                              c("2","2"),
                              c("3","3"),
                              c("4","2"),
                              c("5","2"),
                              c("6","3"),
                              c("7","3"),
                              c("8","4"),
                              c("9","3")
    ) %>%
      as.data.frame() %>%
      `colnames<-` (c("NumFigs","NumCol"))
    
    ncol<-fig_col_numbers_df$NumCol[fig_col_numbers_df$NumFigs == length(individual_model_bargraph_list)] %>% as.numeric() # Using the fig_col_numbers_df to determine how many columns the model-specific figure should have
    
    individual_model_bargraph_figure <- grid.arrange(
      grobs = individual_model_bargraph_list,
      ncol = ncol,
      top = paste0("SNPs in final ",measure, " Model"))
    
    ggsave(filename = paste0(folder, measure,".png"),
           individual_model_bargraph_figure,
           dpi = 1200,
           width = 10, height = 7,
           units = c("in"))
    
    print(individual_model_bargraph_figure)
    model_bargraphs[[measure]]<-individual_model_bargraph_figure
    
  }
  
  
}

plot_bargraphs_for_no_minor(children_final_models,children_combined_unscaled, "~/Desktop/")
plot_bargraphs_for_no_minor(parent_final_models,parent_combined_unscaled, "~/Desktop/")

##################################################
# Function if minor allele was collapsed
##################################################

# This will still run through everything, but it will at least properly model the plots the minor allele needs to be collapsed
plot_bargraphs_for_collapsed<-function(final_models_df, combined_unscaled_df, folder){
  model_bargraphs<-list()
  model_names<-unique(final_models_df$Measure)
  
  # My input dfs don't have the units of measurement. This new df will allow us to get the label for the y-axis
  y_axis_labels_df <- cbind(
    colnames(parent_combined_unscaled)[!grepl("rs", colnames(parent_combined_unscaled))],
    c("Sex ","BMI","Body Mass (kg)", "Height (cm)", "Waist Circumference (cm)", "Total Sugar (mg)", "Added Sugar (mg)", "Caloric Intake (KCAL)", "Age (Years)",  "Waist Circumference-to-Height Ratio")) %>%
    as.data.frame() %>%
    `colnames<-`(c("measure", "label"))
  
  for (i in 1:length(model_names)){
    max_y_value<-NULL # Used so that all plots for a given measure have the same ylim
    measure<-model_names[i]
    y_label<-y_axis_labels_df$label[y_axis_labels_df$measure == paste0(measure)] 
    individual_model_bargraph_list<-list()
    SNP_IDs<-final_models_df$SNP[final_models_df$Measure == model_names[i]]
    Gene<-final_models_df$Gene[final_models_df$Measure == model_names[i]]
    
    for (SNP in SNP_IDs) {
      df <- combined_unscaled_df[, c(model_names[i], SNP)] %>% as.data.frame()
      
      # Recode the alleles to combine 0's and 1's in one bar together
      df[, 2] <- ifelse(df[, 2] %in% c(0, 1), "0-1", "2")
      
      mean_value_by_allele <- aggregate(df[, 1] ~ df[, 2], 
                                        data = df, 
                                        FUN = mean, 
                                        na.rm = TRUE) %>%
        setNames(c("Allele", "Mean_Measure")) %>%
        mutate(Allele = factor(Allele, levels = c("0-1", "2"))) # Use factor levels to ensure correct ordering      
      # Generate x_labels based on SNP IDs
      row_index <- match(SNP, final_models_df$SNP) # Find the corresponding row in final_models_df
      major_allele <- final_models_df$major[row_index]
      minor_allele <- final_models_df$minor[row_index]
      het_freq <- final_models_df$het_freq[row_index]
      major_freq <- final_models_df$major_freq[row_index]
      x_labels <- c(paste0(minor_allele, minor_allele, "/",major_allele, minor_allele, "\n", het_freq, "%"), paste0(major_allele, major_allele, "\n", major_freq, "%"), paste0(" "))
      
      plot <- ggplot(mean_value_by_allele, aes(x = Allele, y = Mean_Measure, fill = Allele)) +
        geom_col() +
        labs(title = paste0(SNP, ": ", Gene),
             x = "SNP Allele",
             y = paste0("Mean ", y_label)) +
        scale_fill_manual(values = c("#56B4E9", "#009E73", "white")) +
        scale_x_discrete(labels = x_labels) +  # Use the correct x_labels
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "none") +
        geom_text(aes(label = round(Mean_Measure, 2)), nudge_y = max(mean_value_by_allele$Mean_Measure)*0.05)
      # coord_cartesian(ylim = c(0, max(mean_value_by_allele$Mean_Measure) * 1.1))
      
      individual_model_bargraph_list[[SNP]]<-plot
      max_y_value <- max(max_y_value, max(mean_value_by_allele$Mean_Measure)) # Used so that all plots for a given measure have the same ylim; update the maximum y-axis value if the current plot's maximum is higher
      
    }
    
    # Use the maximum y-axis value as the common ylim for all plots in the list
    for (j in 1:length(individual_model_bargraph_list)) {
      individual_model_bargraph_list[[j]] <- individual_model_bargraph_list[[j]] +
        coord_cartesian(ylim = c(0, max_y_value * 1.1))   
    }
    
    # Df to determine how many cols the final figure should have
    fig_col_numbers_df<-rbind(c("1","1"),
                              c("2","2"),
                              c("3","3"),
                              c("4","2"),
                              c("5","2"),
                              c("6","3"),
                              c("7","3"),
                              c("8","4"),
                              c("9","3")
    ) %>%
      as.data.frame() %>%
      `colnames<-` (c("NumFigs","NumCol"))
    
    ncol<-fig_col_numbers_df$NumCol[fig_col_numbers_df$NumFigs == length(individual_model_bargraph_list)] %>% as.numeric() # Using the fig_col_numbers_df to determine how many columns the model-specific figure should have
    
    individual_model_bargraph_figure <- grid.arrange(
      grobs = individual_model_bargraph_list,
      ncol = ncol,
      top = paste0("SNPs in final ",measure, " Model"))
    
    ggsave(filename = paste0(folder, measure,".png"),
           individual_model_bargraph_figure,
           dpi = 1200,
           width = 10, height = 7,
           units = c("in"))
    
    print(individual_model_bargraph_figure)
    model_bargraphs[[measure]]<-individual_model_bargraph_figure
    
  }
  
  
}

plot_bargraphs_for_collapsed(children_final_models,children_combined_unscaled, "~/Desktop/")
