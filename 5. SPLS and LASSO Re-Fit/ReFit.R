# This script refits the SPLS and lasso results to a simple linear regression
setwd("~/Desktop")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(stringr)
library(gridExtra)
library(grid)

##################################################
# Selecting SNPs based on frequency
##################################################

# Load in results
frequency_of_child_spls<-read.csv("SPLS_Child_Frequency.csv",check.names = FALSE,row.names = 1)
frequency_of_child_lasso<-read.csv("Lasso_Child_Frequency.csv",check.names = FALSE,row.names = 1)
frequency_of_parent_spls<-read.csv("SPLS_Parent_Frequency.csv",check.names = FALSE,row.names = 1)
frequency_of_parent_lasso<-read.csv("Lasso_Parent_Frequency.csv",check.names = FALSE,row.names = 1)

# Creating a function where, when inputting the rsID, it returns the gene name
find_gene_from_SNP <- function(rsID) {
  matching_rows <- SugarGenes_SNPs$SNP == rsID
  
  # Check if any matching rows are found
  if (any(matching_rows)) {
    # Extract the corresponding gene(s)
    gene <- SugarGenes_SNPs$Gene[matching_rows]
  } else {
    gene <- "No matching SNP found"
  }
  
  return(gene)
}

SugarGenes_SNPs<-read.csv("Sugar_locations.csv")
SugarGenes_SNPs$SNP<-tolower(SugarGenes_SNPs$SNP)

# Function to create SNP and gene list from model results
find_SNP_and_gene_list_from_model_frequencies <- function(frequency_df, frequency_threshold) {
  # Creating a list with the SNP rsIDs that were above 10% (out of 10 runs) for SPLS
  # Initialize an empty list to store the results
  list_of_SNPs <- list()
  list_of_genes <- list()
  
  # Iterate through each column
  for (col_index in 1:ncol(frequency_df)) {
    # Get column name
    col_name <- colnames(frequency_df)[col_index]
    
    # Collect row names with values above the frequency_threshold for the column
    snps_above_frequency_threshold <- rownames(frequency_df)[frequency_df[, col_index] > frequency_threshold]
    
    # Assign the collected row names as a single string to the nested list
    list_of_SNPs[[col_name]] <- paste(snps_above_frequency_threshold, collapse = ",")
    
    genes <- c()
    for (i in snps_above_frequency_threshold) {
      gene <- find_gene_from_SNP(i)
      genes <- c(genes, gene)
    }
    
    # Assign the collected genes as a single string to the list
    list_of_genes[[col_name]] <- paste(genes, collapse = ",")
    rm(col_name, snps_above_frequency_threshold, genes)
  }
  
  # Print the complete list of row names
  output_lists <- list(list_of_SNPs, list_of_genes)
  return(output_lists)
}

# The 1 and 10 represent the threshold we chose of 10%
Child_SNPs_for_SPLS<-find_SNP_and_gene_list_from_model_frequencies(frequency_of_child_spls,1)[[1]]
Child_SNPs_for_LASSO<-find_SNP_and_gene_list_from_model_frequencies(frequency_of_child_lasso,10)[[1]]
Parent_SNPs_for_SPLS<-find_SNP_and_gene_list_from_model_frequencies(frequency_of_parent_spls,1)[[1]]
Parent_SNPs_for_LASSO<-find_SNP_and_gene_list_from_model_frequencies(frequency_of_parent_lasso,10)[[1]]

# Find the union or intersection of SNPs found (liberal or conservative approach)
identify_SNPs_for_model<-function(SNP_list1,SNP_list2, union_or_intersection) {
  SNPs_df <- data.frame(Measure = character(), Group_of_SNPs = character(),
                        stringsAsFactors = FALSE)
  check_if_models_identical<-identical(names(SNP_list1),names(SNP_list2)) 
  if (check_if_models_identical=="FALSE") {
    print("Models are not the same in input or are not in the same order")
    break
  }
  
  for (i in 1:length(SNP_list1)){
    model_name<-names(SNP_list1)[i]
    
    # Getting SNPs from both algorithms for each model/measure
    list1_SNPs<-SNP_list1[[i]] %>%
      strsplit(., ",")%>%
      unlist
    list2_SNPs<-SNP_list2[[i]] %>%
      strsplit(., ",") %>%
      unlist
    
    # For intersection
    if (union_or_intersection=="intersection"){
      intersection_of_SNPs <- intersect(list1_SNPs, list2_SNPs) 
      
      if (length(intersection_of_SNPs) > 0) {
        intersection_of_SNPs <- intersection_of_SNPs %>%
          paste(collapse = ", ")
        temp_df <- data.frame(Measure =model_name,
                              Group_of_SNPs = intersection_of_SNPs,
                              stringsAsFactors = FALSE)
        
        SNPs_df <- rbind(SNPs_df, temp_df)}
    }
    # For union
    else if (union_or_intersection=="union"){
      union_of_SNPs <- union(list1_SNPs, list2_SNPs) 
      
      if (length(union_of_SNPs) > 0) {
        union_of_SNPs <- union_of_SNPs %>%
          paste(collapse = ", ")
        
        temp_df <- data.frame(Measure = model_name,
                              Group_of_SNPs = union_of_SNPs,
                              stringsAsFactors = FALSE)
        
        SNPs_df <- rbind(SNPs_df, temp_df)
      }
    }
  }
  
  return(SNPs_df)
}

# Run function to get results
Child_union<-identify_SNPs_for_model(Child_SNPs_for_LASSO,Child_SNPs_for_SPLS,"union")
#Child_union
Child_intersection<-identify_SNPs_for_model(Child_SNPs_for_LASSO,Child_SNPs_for_SPLS,"intersection")
#Child_intersection

Parent_union<-identify_SNPs_for_model(Parent_SNPs_for_LASSO,Parent_SNPs_for_SPLS,"union")
#Parent_union
Parent_intersection<-identify_SNPs_for_model(Parent_SNPs_for_LASSO,Parent_SNPs_for_SPLS,"intersection")
#Parent_intersection

################################################
# Re-fitting linear models with selected SNPs
################################################
Parent_combined_df <- read.csv("p_anthropometric_data_and_allele_freq_combined_scaled.csv", row.names = 1,check.names = FALSE)
Child_combined_df <- read.csv("anthropometric_data_and_allele_freq_combined_scaled.csv", row.names = 1,check.names = FALSE)
colnames(Child_combined_df) <- c("Sex","BMI Z-Score","Body Mass","Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Age", "Waist Circumference to Height Ratio", "rs1801278", "rs35874116", "rs5400", "rs5417", "rs5418", "rs12970134", "rs1421085", "rs2943641", "rs7102710", "rs17683430", "rs9939609", "rs17817449", "rs1558902", "rs17300539", "rs17782313", "rs4994", "rs1801260")
colnames(Parent_combined_df) <- c("Sex","BMI","Body Mass","Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Age", "Waist Circumference to Height Ratio", "rs1801278", "rs35874116", "rs5400", "rs5417", "rs5418", "rs12970134", "rs1421085", "rs2943641", "rs7102710", "rs17683430", "rs9939609", "rs17817449", "rs1558902", "rs17300539", "rs17782313", "rs4994", "rs1801260")

# Function to refit results to simple linear regression
refit_models <- function(groups_of_SNPs_df, data_type) {
  list_of_refit_models <- list()
  for (i in 1:nrow(groups_of_SNPs_df)) {
    measure <- groups_of_SNPs_df$Measure[i]
    SNP_group <- groups_of_SNPs_df[i, 2] %>%
      strsplit(., ", ")
    
    if (data_type == "child") {
      combined_data <- Child_combined_df %>%
        select(c(all_of(measure), unlist(SNP_group)))
    } else if (data_type == "parent") {
      combined_data <- Parent_combined_df %>%
        select(c(all_of(measure), unlist(SNP_group)))
    } else {
      stop("Invalid data_type. Please specify 'child' or 'parent'.")
    }
    
    formula <- as.formula(paste0("`", measure, "` ~ ."))  # Enclose variable name with ``
    model <- lm(formula, data = na.omit(combined_data))
    list_of_refit_models[[i]] <- model  # Use double square brackets [[]] for list assignment
  }
  names(list_of_refit_models) <- groups_of_SNPs_df$Measure
  
  # Create an empty vector to store the MAE values
  mae_values <- numeric(length(list_of_refit_models))
  
  # Calculate MAE for each model
  for (i in 1:length(list_of_refit_models)) {
    measure <- names(list_of_refit_models)[i]
    model <- list_of_refit_models[[i]]
    predictions <- predict(model)
    actual <- model$model[[measure]]
    
    mae_val <- mean(abs(predictions - actual))  # Calculate MAE
    print(mae_val)
    mae_values[i] <- mae_val
  }
  
  df_with_MAE <- groups_of_SNPs_df %>%
    cbind(., MAE = mae_values) %>%
    setNames(c("Measure", "Group of SNPs", "MAE"))
  
  return(df_with_MAE)
}

# Function to collapse rows for our tables
collapse_rows_df <- function(df, variable){
  
  group_var <- enquo(variable)
  
  df %>%
    group_by(!! group_var) %>%
    mutate(groupRow = 1:n()) %>%
    ungroup() %>%
    mutate(!!quo_name(group_var) := ifelse(groupRow == 1, as.character(!! group_var), "")) %>%
    select(-c(groupRow))
}

# Function to add our coefficients to the table
create_coefficient_dfs <- function(SNP_groups_with_MAE_df, data_type) {
  # The input df is the df created in the last section
  
  SNP_groups_with_lowest_MAE <- SNP_groups_with_MAE_df %>%
    group_by(Measure) %>%
    filter(MAE == min(MAE)) %>% 
    as.data.frame()
  
  coefficient_summary_dfs_list <- list()
  for (i in 1:nrow(SNP_groups_with_lowest_MAE)) {
    
    new_df <- SNP_groups_with_lowest_MAE[i,]
    measure <- new_df[1, 1]
    SNPs <- new_df[1, 2] %>%
      strsplit(., ", ")
    
    num_snps <- length(SNPs[[1]])
    MAE <- new_df[1, 3] %>%
      round(., 5) %>%
      c(., rep("", (num_snps - 1)))
    
    if (data_type == "child") {
      combined_data <- Child_combined_df %>%
        select(c(all_of(measure), unlist(SNPs)))
    } else if (data_type == "parent") {
      combined_data <- Parent_combined_df %>%
        select(c(all_of(measure), unlist(SNPs)))
    } else {
      stop("Invalid data_type. Please specify 'child' or 'parent'.")
    }
    
    formula <- as.formula(paste0("`", measure, "` ~ ."))  # Enclose variable name with ``
    model <- lm(formula, data = na.omit(combined_data))
    summary(model)
    
    r_sq <- summary(model)$r.squared %>%
      round(., 5) %>%
      c(., rep("", (num_snps - 1)))
    
    adj_r_sq <- summary(model)$adj.r.squared %>%
      round(., 4) %>%
      c(., rep("", (num_snps - 1)))
    
    coefficients_df <- summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(Coefficient = paste(round(Estimate, 3), round(`Std. Error`, 3), sep = " ± ")) %>%
      select(rowname, Coefficient) %>%
      filter(rowname != "(Intercept)")
    
    model_summary_df <- cbind(MAE, r_sq, adj_r_sq, coefficients_df) %>%
      `colnames<-`(c("MAE", "R. Sq.", "Adj. R. Sq.", "SNP", "Coefficient ± Std. Error"))
    
    coefficient_summary_dfs_list[[i]] <- model_summary_df
    
  }
  
  names(coefficient_summary_dfs_list) <- SNP_groups_with_lowest_MAE$Measure
  
  # Create an empty list to store the modified data frames
  combined_dfs <- list()
  
  # Iterate over the list of data frames
  for (i in seq_along(coefficient_summary_dfs_list)) {
    # Get the current data frame
    df <- coefficient_summary_dfs_list[[i]]
    
    # Add a new column with the data frame name
    df$Measure <- names(coefficient_summary_dfs_list)[i]
    
    # Append the modified data frame to the combined list
    combined_dfs[[i]] <- df
  }
  
  # Combine the modified data frames into a single data frame
  combined_df <- do.call(rbind, combined_dfs) %>%
    select(Measure, everything())
  
  return(combined_df)
}

# Models with SNP Intersection
children_intersection_refit_models <- refit_models(Child_intersection, data_type = "child") %>% 
  create_coefficient_dfs(data_type = "child") %>% 
  mutate(Gene = map_chr(SNP, find_gene_from_SNP)) %>%
  select(Measure, MAE, 'R. Sq.', 'Adj. R. Sq.', SNP, Gene, everything()) %>%
  mutate_all(~ifelse(is.na(.), "", as.character(.))) %>%
  mutate(coef = as.numeric(gsub(" ±.*", "", .$`Coefficient ± Std. Error`)))

children_intersection_refit_models
write.csv(children_intersection_refit_models,"~/Desktop/children_intersection_refit_models.csv")

parent_intersection_refit_models <- refit_models(Parent_intersection, data_type = "parent") %>% 
  create_coefficient_dfs(data_type = "parent") %>% 
  mutate(Gene = map_chr(SNP, find_gene_from_SNP)) %>%
  select(Measure, MAE, 'R. Sq.', 'Adj. R. Sq.', SNP, Gene, everything()) %>%
  mutate_all(~ifelse(is.na(.), "", as.character(.))) %>%
  mutate(coef = as.numeric(gsub(" ±.*", "", .$`Coefficient ± Std. Error`)))

parent_intersection_refit_models
write.csv(parent_intersection_refit_models,"~/Desktop/parent_intersection_refit_models.csv")

# Create kable of our intersection results
# Create a gradient of colors for the coef column
color_scale <- scales::col_numeric(palette = c("indianred1","white", "seagreen1"), domain = range(children_intersection_refit_models$coef))
p_color_scale <- scales::col_numeric(palette = c("indianred1","white", "seagreen1"), domain = range(parent_intersection_refit_models$coef))

children_SNP_intersection_kable<-children_intersection_refit_models%>%
  collapse_rows_df(Measure)%>%
  select(-coef) %>%
  `colnames<-`(c("Linear Model Response Variable","MAE","R. Sq.","Adj. R. Sq.", "SNP", "Gene", "Coefficient ± Std. Error")) %>%
  kable(caption = "Children Data Linear Regression Models", align= "l") %>% 
  kable_classic() %>%
  column_spec(column = 7, background = color_scale(children_intersection_refit_models$coef)) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE, hline_after = TRUE) %>%
  column_spec(1,bold=TRUE, border_right = TRUE) %>%
  row_spec(c(0), extra_css = "border-bottom: 4px solid")%>%
  row_spec(c(2,3,4,5,6), extra_css = "border-bottom: 1px solid") #%>% # Change this line for wherever you need lines spliting rows
  #as_image(., width = 7, file="~/Desktop/children_SNP_intersection_kable.png", zoom=15) # Uncomment this to save as png

children_SNP_intersection_kable

parent_SNP_intersection_kable<-parent_intersection_refit_models%>%
  collapse_rows_df(Measure)%>%
  select(-coef) %>%
  `colnames<-`(c("Linear Model Response Variable","MAE","R. Sq.","Adj. R. Sq.", "SNP", "Gene", "Coefficient ± Std. Error")) %>%
  kable(caption = "Parent Data Linear Regression Models", align= "l") %>% 
  kable_classic() %>%
  column_spec(column = 7, background = p_color_scale(parent_intersection_refit_models$coef)) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE, hline_after = TRUE) %>%
  column_spec(1,bold=TRUE, border_right = TRUE) %>%
  row_spec(c(0), extra_css = "border-bottom: 4px solid")%>%
  row_spec(c(2,3,4,6,7), extra_css = "border-bottom: 1px solid") #%>% # Change this line for wherever you need lines splitting rows
  #as_image(., width = 7, file="~/Desktop/parent_SNP_intersection_kable.png", zoom=15) # Uncomment this to save as png

parent_SNP_intersection_kable

#####################################
# Union.... Good to note everything flagged, but will use conservative approach for report
#####################################

# Function to get coefficients for union approach
create_coefficient_dfs_Union <- function(SNP_groups_with_MAE_df, data_type) {
  # The input df is the df created in the last section
  
  SNP_groups_with_lowest_MAE <- SNP_groups_with_MAE_df %>%
    group_by(Measure) %>%
    filter(MAE == min(MAE)) %>% 
    as.data.frame()
  
  coefficient_summary_dfs_list <- list()
  for (i in 1:nrow(SNP_groups_with_lowest_MAE)) {
    
    new_df <- SNP_groups_with_lowest_MAE[i,]
    measure <- new_df[1, 1]
    SNPs <- new_df[1, 2] %>%
      strsplit(., ", ")
    
    num_snps <- length(SNPs[[1]])
    MAE <- new_df[1, 3] %>%
      round(., 5) %>%
      c(., rep("", (num_snps - 1)))
    
    if (data_type == "child") {
      combined_data <- Child_combined_df %>%
        select(c(all_of(measure), unlist(SNPs)))
    } else if (data_type == "parent") {
      combined_data <- Parent_combined_df %>%
        select(c(all_of(measure), unlist(SNPs)))
    } else {
      stop("Invalid data_type. Please specify 'child' or 'parent'.")
    }
    
    formula <- as.formula(paste0("`", measure, "` ~ ."))  # Enclose variable name with ``
    model <- lm(formula, data = na.omit(combined_data))
    summary(model)
    
    r_sq <- summary(model)$r.squared %>%
      round(., 5) %>%
      c(., rep("", (num_snps - 1)))
    
    adj_r_sq <- summary(model)$adj.r.squared %>%
      round(., 4) %>%
      c(., rep("", (num_snps - 1)))
    
    # Create vectors of equal length
    max_len <- max(num_snps, length(MAE), length(r_sq), length(adj_r_sq))
    MAE <- rep(MAE, length.out = max_len)
    r_sq <- rep(r_sq, length.out = max_len)
    adj_r_sq <- rep(adj_r_sq, length.out = max_len)
    
    coefficients_df <- summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(Coefficient = paste(round(Estimate, 3), round(`Std. Error`, 3), sep = " ± ")) %>%
      select(rowname, Coefficient) %>%
      filter(rowname != "(Intercept)")
    
    # Pad the coefficients_df with NAs to have the same number of rows as max_len
    coefficients_df <- rbind(coefficients_df, data.frame(rowname = character(max_len - nrow(coefficients_df)), Coefficient = character(max_len - nrow(coefficients_df))))
    
    model_summary_df <- cbind(MAE, r_sq, adj_r_sq, coefficients_df) %>%
      `colnames<-`(c("MAE", "R. Sq.", "Adj. R. Sq.", "SNP", "Coefficient ± Std. Error"))
    
    coefficient_summary_dfs_list[[i]] <- model_summary_df
    
  }
  
  names(coefficient_summary_dfs_list) <- SNP_groups_with_lowest_MAE$Measure
  
  # Create an empty list to store the modified data frames
  combined_dfs <- list()
  
  # Iterate over the list of data frames
  for (i in seq_along(coefficient_summary_dfs_list)) {
    # Get the current data frame
    df <- coefficient_summary_dfs_list[[i]]
    
    # Add a new column with the data frame name
    df$Measure <- names(coefficient_summary_dfs_list)[i]
    
    # Append the modified data frame to the combined list
    combined_dfs[[i]] <- df
  }
  
  # Combine the modified data frames into a single data frame
  combined_df <- do.call(rbind, combined_dfs) %>%
    select(Measure, everything())
  
  return(combined_df)
}


# children SNP Union models
children_union_refit_models <- refit_models(Child_union, data_type = "child") %>% 
  create_coefficient_dfs_Union(data_type = "child") %>% 
  mutate(Gene = map_chr(SNP, find_gene_from_SNP)) %>%
  select(Measure, MAE, 'R. Sq.', 'Adj. R. Sq.', SNP, Gene, everything()) %>%
  mutate_all(~ifelse(is.na(.), "", as.character(.))) %>%
  mutate(coef = as.numeric(gsub(" ±.*", "", .$`Coefficient ± Std. Error`)))

children_union_refit_models
write.csv(children_union_refit_models,"~/Desktop/children_union_refit_models.csv")

# parent SNP Union models
parent_union_refit_models <- refit_models(Parent_union, data_type = "parent") %>% 
  create_coefficient_dfs_Union(data_type = "parent") %>% 
  mutate(Gene = map_chr(SNP, find_gene_from_SNP)) %>%
  select(Measure, MAE, 'R. Sq.', 'Adj. R. Sq.', SNP, Gene, everything()) %>%
  mutate_all(~ifelse(is.na(.), "", as.character(.))) %>%
  mutate(coef = as.numeric(gsub(" ±.*", "", .$`Coefficient ± Std. Error`)))

parent_union_refit_models
write.csv(parent_union_refit_models,"~/Desktop/parent_union_refit_models.csv")

# Remove extra rows with NA values
children_union_refit_models <- na.omit(children_union_refit_models)
parent_union_refit_models <- na.omit(parent_union_refit_models)

# Create kable of our union results
# Create a gradient of colors for the coef column
color_scale <- scales::col_numeric(palette = c("indianred1","white", "seagreen1"), domain = range(children_union_refit_models$coef))
p_color_scale <- scales::col_numeric(palette = c("indianred1","white", "seagreen1"), domain = range(parent_union_refit_models$coef))


children_SNP_union_kable<-children_union_refit_models%>%
  collapse_rows_df(Measure)%>%
  select(-coef) %>%
  `colnames<-`(c("Linear Model Response Variable","MAE","R. Sq.","Adj. R. Sq.", "SNP", "Gene", "Coefficient ± Std. Error")) %>%
  kable(caption = "Children Data Linear Regression Models", align= "l") %>% 
  kable_classic() %>%
  column_spec(column = 7, background = color_scale(children_union_refit_models$coef)) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE, hline_after = TRUE) %>%
  column_spec(1,bold=TRUE, border_right = TRUE) %>%
  row_spec(c(0), extra_css = "border-bottom: 4px solid")%>%
  row_spec(c(2,3,13,14,26,27,28), extra_css = "border-bottom: 1px solid") #%>% # Change this line for wherever you need to split rows
  #as_image(., width = 7, file="~/Desktop/children_SNP_union_kable.png", zoom=15) # Uncomment this to save as png

children_SNP_union_kable

parent_SNP_union_kable<-parent_union_refit_models%>%
  collapse_rows_df(Measure)%>%
  select(-coef) %>%
  `colnames<-`(c("Linear Model Response Variable","MAE","R. Sq.","Adj. R. Sq.", "SNP", "Gene", "Coefficient ± Std. Error")) %>%
  kable(caption = "Children Data Linear Regression Models", align= "l") %>% 
  kable_classic() %>%
  column_spec(column = 7, background = p_color_scale(parent_union_refit_models$coef)) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE, hline_after = TRUE) %>%
  column_spec(1,bold=TRUE, border_right = TRUE) %>%
  row_spec(c(0), extra_css = "border-bottom: 4px solid")%>%
  row_spec(c(13,15,23,29,31,33,49), extra_css = "border-bottom: 1px solid") #%>% # Change this line for wherever you need to split rows
  #as_image(., width = 7, file="~/Desktop/parent_SNP_union_kable.png", zoom=15) # Uncomment this to save as png

parent_SNP_union_kable
