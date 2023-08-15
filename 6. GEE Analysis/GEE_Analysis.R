# This script performs the GEE analysis
library(tidyverse)
library(dplyr)
library(stringr)
library(geepack)
library(IRanges)
library(lubridate)
library(magrittr)
library(kableExtra)
library(gridExtra)
library(grid)

# Set directory
setwd("~/Desktop")

# Read in the anthropometric and allele files (GEE can work with non-normalized data so we don't need scaled data) 
Parent_combined_df <- read.csv("parent_anthropometric_data_and_allele_freq_combined.csv", header=TRUE)
Child_combined_df <- read.csv("children_anthropometric_data_and_allele_freq_combined.csv", header=TRUE)

# Set Rownames
rownames(Parent_combined_df) <- Parent_combined_df$X
Parent_combined_df$X <- NULL
rownames(Child_combined_df) <- Child_combined_df$X
Child_combined_df$X <- NULL

# Function to scale the children and parent anthropometric data separately
scale_data<-function(anthropometric_data_and_allele_freq_combined_df){
  columns_to_scale <- 2:10
  scaled_df<-anthropometric_data_and_allele_freq_combined_df
  # Apply scaling to the selected columns
  scaled_df[, columns_to_scale] <- scale(scaled_df[, columns_to_scale])
  
  return(scaled_df)
}

# Scale children and adults separately
Child_combined_df <- scale_data(Child_combined_df)
Parent_combined_df <- scale_data(Parent_combined_df)

# Create FamID
Parent_combined_df$Fam_ID <- sub("[A-Z]", "F", rownames(Parent_combined_df))
Child_combined_df$Fam_ID <- sub("[A-Z]", "F", rownames(Child_combined_df))

# Remove BMI and BMI-z columns as there's no validated way to get BMI Z-score in adults
Parent_combined_df$BMI <- NULL
Child_combined_df$BMI.Z.Score <- NULL

# Combine both parents and children data frames
comb_df <- rbind(Parent_combined_df, Child_combined_df)

##################################
# GEE Function
##################################
rsIDs <- colnames(comb_df)[grepl("^rs", colnames(comb_df))] %>%
  paste(collapse = " + ")

# Function to perform GEE analysis
get_gee_results <- function(anthropometric_measurement) {
  df_for_variable <- comb_df %>%
    select(Fam_ID, anthropometric_measurement, matches("^rs"))
  
  formula <- reformulate(rsIDs, response = anthropometric_measurement)
 
  gee_model <- geeglm(formula,
                      data = df_for_variable,
                      family = gaussian(),
                      id = Fam_ID,
                      corstr = "exchangeable")
  
  summary_results <- summary(gee_model)
  pvalues_df <- summary_results$coefficients[, "Pr(>|W|)"] %>%
    as.data.frame() %>%
    setNames(anthropometric_measurement)
  
  coefficients_df <- coef(gee_model) %>%
    as.data.frame() %>%
    setNames(anthropometric_measurement)
  
  return(list(pvalues_df = pvalues_df, coefficients_df = coefficients_df))
}

##################################
# Run GEE
##################################

# Order by family ID
comb_df <- comb_df %>% arrange(Fam_ID)
variables_list <- c("Body.Mass", "Height", "Waist.Circumference", "Total.Sugar", "Added.Sugar", "Caloric.Intake", "Waist.Circumference.to.Height.Ratio")

# Remove F0 from family ID to just keep number
comb_df$Fam_ID <- substr(comb_df$Fam_ID, start = nchar(comb_df$Fam_ID) - 2, stop = nchar(comb_df$Fam_ID))

# Initialize the number of tests for Bonferroni correction (number of measurements in variables_list)
num_tests <- length(variables_list)

# Create an empty list to store results
results_list <- list()

# Loop through each measurement in variables_list
for (measurement in variables_list) {
  # Call the function to get the p-values and coefficients for the current measurement
  results <- get_gee_results(measurement)
  
  # Apply Bonferroni correction to p-values
  results$pvalues_df <- results$pvalues_df * num_tests
  results$pvalues_df[results$pvalues_df > 1] <- 1
  
  # Replace p-values greater than 0.05 with NA and corresponding coefficients with NA
  results$pvalues_df[results$pvalues_df > 0.05] <- NA
  results$coefficients_df[is.na(results$pvalues_df), ] <- NA
  
  # Combine p-values and coefficients together in a data frame
  result_df <- cbind(results$pvalues_df, results$coefficients_df)
  
  # Set the column names of the result data frame
  colnames(result_df) <- c("p-value", "coefficient")
  
  # Add the result data frame to the list with the measurement name as the key
  results_list[[measurement]] <- result_df
}

# Coefficients df
# Initialize an empty list to store the coefficient vectors
coefficients_list <- list()

# Iterate over the anthropometric measurements
for (measurement in variables_list) {
  
  # Extract the coefficient column from the result data frame for the current measurement
  coefficients <- results_list[[measurement]]$coefficient
  
  # Set 0 for rows with no coefficient
  coefficients[is.na(coefficients)] <- 0
  
  # Add the coefficient vector to the list with the measurement name as the key
  coefficients_list[[measurement]] <- coefficients
}

# Create the coefficients dataframe
coefficients_df <- do.call(cbind, coefficients_list)

# Set the column names of the coefficients dataframe
colnames(coefficients_df) <- c("Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio")
rownames(coefficients_df) <- rownames(result_df)
coefficients_df <- coefficients_df[-1, ] #Remove intercept row

#################################
# Create Coefficients Heatmap (Not used in report as table shows the same results)
#################################
# Code created with Nishita Sharif to keep plots consistent throughout both our research

coefficients_df <- as.data.frame(coefficients_df)
coefficients_df_long <- coefficients_df %>% rownames_to_column() %>% gather(colname, value, -rowname)
coefficients_df_long$order<-factor(coefficients_df_long$colname, levels=c("Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio"))

gee_coef_heatmap <- ggplot(coefficients_df_long, 
                                   aes(x = order, # Anthropocentric measure
                                       y = rowname, # rsIDs
                                       fill = value, # Coefficient
                                       label = round(value, 2))) + 
  geom_tile(colour = "black") + # Add grid lines between tiles
  geom_text(colour = "black", # Add labels for the values
            size = 2, 
            hjust = 0.5, 
            vjust = 0.5,
            data = subset(coefficients_df_long, value != 0) # shows values not equal to 0
  ) +  
  scale_fill_gradient2(low = "indianred1",mid = "white", high = "seagreen1", limits = c(-6, 6)) +
  labs(y = "SNP rsID", x = "Anthropometric Measurement") +
  ggtitle("GEE Coefficients") +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(angle = 25, hjust = 1, size = 5),
    axis.title.y = element_text(size = 8),
    plot.title = element_text( size = 12, hjust = 0.5)
  ) +
  scale_x_discrete(expand=c(0,0),labels = function(x) stringr::str_wrap(x, width = 15))+
  scale_y_discrete(expand=c(0,0))

gee_coef_heatmap

output_file <- paste0("~/Desktop/gee_coef_heatmap.png")
ggsave(filename = output_file, 
       gee_coef_heatmap,
       dpi = 600, width = 6, height = 4)

#################################
# Significant SNPs table
#################################
SugarGenes_SNPs<-read.csv("Sugar_locations.csv")
SugarGenes_SNPs$SNP<-tolower(SugarGenes_SNPs$SNP)

# Find gene from SNP function
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

coefficients_df_long <- filter(coefficients_df_long, value != 0)

for (i in 1:nrow(coefficients_df_long)){
  coefficients_df_long$Gene[i]<-find_gene_from_SNP(coefficients_df_long$rowname[i])
}

# Create a gradient of colors for the Score column
color_scale <- scales::col_numeric(palette = c("indianred1","white", "seagreen1"), domain = range(coefficients_df_long$value))

coefficients_df_long <- coefficients_df_long %>%
  select(colname, rowname, Gene, everything()) %>%
  mutate_at(vars(-value), ~ ifelse(is.na(.), "", as.character(.)))
coefficients_df_long$value <- as.numeric(coefficients_df_long$value)

# Function to collapse rows
collapse_rows_df <- function(df, variable){
  
  group_var <- enquo(variable)
  
  df %>%
    group_by(!! group_var) %>%
    mutate(groupRow = 1:n()) %>%
    ungroup() %>%
    mutate(!!quo_name(group_var) := ifelse(groupRow == 1, as.character(!! group_var), "")) %>%
    select(-c(groupRow))
}

# Create GEE kable
Gee_coefficients_kable <- coefficients_df_long %>%
  select(-order) %>%
  collapse_rows_df(colname) %>%
  `colnames<-`(c("Measure", "SNP", "Gene", "Coefficient")) %>%
  mutate(Coefficient = sprintf("%.3f", Coefficient))  # Format coefficient to three decimal places

Gee_coefficients_kable <- Gee_coefficients_kable %>%
  kable(caption = "Significantly Associated SNPs in Family Unit Analysis", align = "l") %>% 
  kable_classic() %>%
  column_spec(column = 4, background = color_scale(coefficients_df_long$value)) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE) %>%
  row_spec(c(5, 11), extra_css = "border-bottom: 1px solid") %>%
  row_spec(c(0), extra_css = "border-bottom: 4px solid") %>%
  as_image(., width = 7, file="~/Desktop/GEE_SNP_kable.png", zoom=15) # Uncomment this to save as png

Gee_coefficients_kable

#Saves
write.csv(coefficients_df_long,"~/Desktop/gee_coefficients_long.csv")

