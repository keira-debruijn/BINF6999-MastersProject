# This script performs the Lasso analysis
library(tidyverse)
library(dplyr)
library(stringr)
library(kableExtra)
library(gridExtra)
library(grid)

# Set directory
setwd("~/Desktop")

# Read in combined scaled anthropometric and allele files
Parent_combined_df <- read.csv("p_anthropometric_data_and_allele_freq_combined_scaled.csv", header=TRUE)
Child_combined_df <- read.csv("anthropometric_data_and_allele_freq_combined_scaled.csv", header=TRUE)

# Set Rownames
rownames(Parent_combined_df) <- Parent_combined_df$X
Parent_combined_df$X <- NULL
rownames(Child_combined_df) <- Child_combined_df$X
Child_combined_df$X <- NULL

Child_combined_df <- na.omit(Child_combined_df)
Parent_combined_df <- na.omit(Parent_combined_df)

# Keep only oldest children
Child_combined_df <- Child_combined_df %>% filter(grepl("A", rownames(.)))

##################################
# Lasso Function
##################################

library(glmnet)

lasso_analysis <- function(x, y) {
  # Create empty list to store models
  lasso_results <- list()
  
  # Empty list to store all coef data frames 
  list_of_lasso_coefficient_dfs <- list()
  
  for (i in 1:100) {
    # Fit cross-validated LASSO models for each response variable
    cv_models <- list()
    response_vars <- colnames(y)
    
    for (response_var in response_vars) {
      cv_models[[response_var]] <- cv.glmnet(x, y[, response_var], alpha = 1)
    }
    
    # Run LASSO for each response variable using the optimal lambda from cross-validation
    for (response_var in response_vars) {
      lasso_results[[response_var]] <- glmnet(x, y[, response_var], alpha = 1, lambda = cv_models[[response_var]]$lambda.min)
    }
    
    # Empty list for row names
    Identified_SNPs <- list()
    
    for (model_name in names(lasso_results)) {
      coefs <- as.matrix(coef(lasso_results[[model_name]]))
      rowname_list <- rownames(coefs)
      rowname_list <- rowname_list[coefs != 0 & !grepl("Intercept", rowname_list)]
      Identified_SNPs[[model_name]] <- rowname_list
    }
    
    # Create an empty data frame
    lasso_coefficients <- matrix(ncol = length(response_vars), nrow = length(colnames(x)))
    
    # Add column names based on the response variable names
    colnames(lasso_coefficients) <- response_vars
    
    # Add row names
    rownames(lasso_coefficients) <- colnames(x)
    
    for (model_name in names(lasso_results)) {
      coefs <- as.matrix(coef(lasso_results[[model_name]])[-1, , drop = FALSE])
      lasso_coefficients[, model_name] <- coefs
    }
    
    list_of_lasso_coefficient_dfs[[i]] <- as.data.frame(lasso_coefficients)
  }
  
  return(list_of_lasso_coefficient_dfs)
}

#################################
# Run Lasso
#################################
set.seed(100)

# Split combined data frame into x and y variables
Child_y <- as.data.frame(Child_combined_df[,2:10])
Child_y[, 1:9] <- lapply(Child_y[, 1:9], as.numeric)
Child_y$Age <- NULL
Child_x <- as.data.frame(Child_combined_df[, 11:27])
Child_x[, 1:17] <- lapply(Child_x[, 1:17], as.numeric)
Child_x <- as.matrix(Child_x)
Parent_y <- as.data.frame(Parent_combined_df[,2:10])
Parent_y[, 1:9] <- lapply(Parent_y[, 1:9], as.numeric)
Parent_y$Age <- NULL
Parent_x <- as.data.frame(Parent_combined_df[, 11:27])
Parent_x[, 1:17] <- lapply(Parent_x[, 1:17], as.numeric)
Parent_x <- as.matrix(Parent_x)

#Run Lasso
children_lasso_results <- lasso_analysis(Child_x, Child_y)
parents_lasso_results <- lasso_analysis(Parent_x, Parent_y)

# Count how many times each SNP is significant
# Initialize the count dataframe
frequency_of_child_lasso <- data.frame(matrix(0, nrow = nrow(children_lasso_results[[1]]), ncol = ncol(children_lasso_results[[1]])))
rownames(frequency_of_child_lasso) <- rownames(children_lasso_results[[1]])
colnames(frequency_of_child_lasso) <- c("BMI Z-Score", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio")

frequency_of_parent_lasso <- data.frame(matrix(0, nrow = nrow(parents_lasso_results[[1]]), ncol = ncol(parents_lasso_results[[1]])))
rownames(frequency_of_parent_lasso) <- rownames(parents_lasso_results[[1]])
colnames(frequency_of_parent_lasso) <- c("BMI", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio")

# Update the count dataframe
for (df in children_lasso_results) {
  frequency_of_child_lasso <- frequency_of_child_lasso + (df != 0)
}

for (df in parents_lasso_results) {
  frequency_of_parent_lasso <- frequency_of_parent_lasso + (df != 0)
}
#frequency_of_child_lasso
#frequency_of_parent_lasso

# Save results
write.csv(children_lasso_results, file = "Children_Lasso_Results.csv", row.names = TRUE)
write.csv(parents_lasso_results, file = "Parent_Lasso_Results.csv", row.names = TRUE)
write.csv(frequency_of_child_lasso, file = "Lasso_Child_Frequency.csv", row.names = TRUE)
write.csv(frequency_of_parent_lasso, file = "Lasso_Parent_Frequency.csv", row.names = TRUE)

# Creating a table figure with kable
frequency_of_child_lasso %>%
  kable(caption = "Frequency of significant factors with Lasso") %>% 
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE, width = "2em") %>%
  column_spec(2:9, width = "6em") %>%
  as_image(., width = 7, file="~/Desktop/Child_Lasso_frequency_table.png", zoom=7)

# Creating a table figure with kable
frequency_of_parent_lasso %>%
  kable(caption = "Frequency of significant factors with Lasso") %>% 
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE, width = "2em") %>%
  column_spec(2:9, width = "6em") %>%
  as_image(., width = 7, file="~/Desktop/Parent_Lasso_frequency_table.png", zoom=7)

#################################
# Create Heatmap
#################################
# Code created with Nishita Sharif to keep plots consistent throughout both our research

# Creating a heatmap with the freuqency that a SNP was found to be significant when doing SPLS analysis
frequency_of_child_lasso_long <- frequency_of_child_lasso %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

frequency_of_parent_lasso_long <- frequency_of_parent_lasso %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

frequency_of_child_lasso_long$order<-factor(frequency_of_child_lasso_long$colname, levels=c("BMI Z-Score", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio"))
frequency_of_parent_lasso_long$order<-factor(frequency_of_parent_lasso_long$colname, levels=c("BMI", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio"))

child_lasso_freq_heatmap <- ggplot(frequency_of_child_lasso_long, 
                                  aes(x = order, # Anthropocentric measure
                                      y = rowname, # rsIDs
                                      fill = value, # Coefficient
                                      label = round(value, 2))) + 
  geom_tile(colour = "black") + # Add grid lines between tiles
  geom_text(colour = "black", # Add labels for the values
            size = 2, 
            hjust = 0.5, 
            vjust = 0.5,
            data = subset(frequency_of_child_lasso_long, value >= 1)
  ) +  
  scale_fill_gradient2(low = "white",high = "darkorange", limits = c(0, 100)) +
  labs(y = "SNP rsID", x = "Anthropometric Measurement") +
  ggtitle("Lasso Child") +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(angle = 25, hjust = 1, size = 5),
    axis.title.y = element_text(size = 8),
    plot.title = element_text( size = 12, hjust = 0.5)
  ) +
  scale_x_discrete(expand=c(0,0),labels = function(x) stringr::str_wrap(x, width = 15))+
  scale_y_discrete(expand=c(0,0))

parent_lasso_freq_heatmap <- ggplot(frequency_of_parent_lasso_long, 
                                   aes(x = order, # Anthropocentric measure
                                       y = rowname, # rsIDs
                                       fill = value, # Coefficient
                                       label = round(value, 2))) + 
  geom_tile(colour = "black") + # Add grid lines between tiles
  geom_text(colour = "black", # Add labels for the values
            size = 2, 
            hjust = 0.5, 
            vjust = 0.5,
            data = subset(frequency_of_parent_lasso_long, value >= 1) 
  ) +  
  scale_fill_gradient2(low = "white",high = "darkorange", limits = c(0, 100)) +
  labs(y = "SNP rsID", x = "Anthropometric Measurement") +
  ggtitle("Lasso Parent") +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(angle = 25, hjust = 1, size = 5),
    axis.title.y = element_text(size = 8),
    plot.title = element_text( size = 12, hjust = 0.5)
  ) +
  scale_x_discrete(expand=c(0,0),labels = function(x) stringr::str_wrap(x, width = 15))+
  scale_y_discrete(expand=c(0,0))

#child_lasso_freq_heatmap
#parent_lasso_freq_heatmap

output_file <- paste0("~/Desktop/child_lasso_freq_heatmap.png")
output_file_parent <- paste0("~/Desktop/parent_lasso_freq_heatmap.png")

ggsave(filename = output_file, 
       child_lasso_freq_heatmap,
       dpi = 600, width = 6, height = 4)
ggsave(filename = output_file_parent, 
       parent_lasso_freq_heatmap,
       dpi = 600, width = 6, height = 4)
