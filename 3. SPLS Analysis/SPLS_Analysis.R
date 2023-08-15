# This script performs the SPLS analysis
library(spls)
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

# Keep only oldest children
Child_combined_df <- Child_combined_df %>% filter(grepl("A", rownames(.)))

##################################
# SPLS Function
##################################
library(spls)

spls_analysis <- function(x, y, anthropometric_measurements) {
  # Initializing empty list; each item in this list will be a dataframe with the SPLS beta coefficients
  list_of_spls_coefficient_dfs <- list()
  
  for (p in 1:10) {
    # Create an empty dataframe
    spls_coefficients_df <- matrix(ncol = length(anthropometric_measurements), 
                                   nrow = length(colnames(x)[grep("^rs", colnames(x))])) %>% 
      as.data.frame()
    
    # Add column names
    colnames(spls_coefficients_df) <- anthropometric_measurements
    
    # Add row names
    rownames(spls_coefficients_df) <- colnames(x)[grep("^rs", colnames(x))]
    
    # Create an empty list to store the results
    results_list <- list()
    
    for (i in seq_along(anthropometric_measurements)) {
      # Find optimal parameters for k and eta
      cvs_model <- cv.spls(x, y[, i], K = c(1:15), eta = seq(0.01, 0.99, 0.01))
      k_value <- cvs_model$K.opt
      eta <- cvs_model$eta.opt
      
      spls_model <- spls(x, y[, i], K = k_value, eta = eta) # Fitting the SPLS model with the optimal K and eta values
      spls_output <- capture.output(print(spls_model)) # Assigning the output to an object; will use regex to get the selected variables in the next line
      elements_vector <- unlist(strsplit(gsub("\t", " ", capture.output(print(spls_model)))[length(spls_output)], "\\s+")) # The vector created here has the selected variables (i.e., SNPs)
      
      spls_coefficients_df[i] <- spls_model$betahat
      # Assign the elements_vector as a list item with the corresponding name
      results_list[[i]] <- elements_vector
    }
    
    # Assign column names to the result list
    names(results_list) <- anthropometric_measurements
    
    # Add the results to the list of coefficient dataframes
    list_of_spls_coefficient_dfs[[p]] <- spls_coefficients_df
  }
  
  return(list_of_spls_coefficient_dfs)
}

#################################
# Run SPLS
#################################
set.seed(100)

# Remove rows with missing values
Child_combined_df <- na.omit(Child_combined_df)
Parent_combined_df <- na.omit(Parent_combined_df)

child_anthropometric_measurements <- c("BMI.Z.Score", "Body.Mass", "Height", "Waist.Circumference", "Total.Sugar", "Added.Sugar", "Caloric.Intake", "Waist.Circumference.to.Height.Ratio")
parent_anthropometric_measurements <- c("BMI", "Body.Mass", "Height", "Waist.Circumference", "Total.Sugar", "Added.Sugar", "Caloric.Intake", "Waist.Circumference.to.Height.Ratio")

# Split combined data frame into x and y variables
Child_y <- as.data.frame(Child_combined_df[,2:10])
Child_y[, 1:9] <- lapply(Child_y[, 1:9], as.numeric)
Child_x <- as.data.frame(Child_combined_df[, 11:27])
Child_x[, 1:17] <- lapply(Child_x[, 1:17], as.numeric)
Child_x <- as.data.frame(Child_x)
Parent_y <- as.data.frame(Parent_combined_df[,2:10])
Parent_y[, 1:9] <- lapply(Parent_y[, 1:9], as.numeric)
Parent_x <- as.data.frame(Parent_combined_df[, 11:27])
Parent_x[, 1:17] <- lapply(Parent_x[, 1:17], as.numeric)
Parent_x <- as.data.frame(Parent_x)

#Run SPLS
children_spls_results <- spls_analysis(Child_x, Child_y, child_anthropometric_measurements)
parents_spls_results <- spls_analysis(Parent_x, Parent_y, parent_anthropometric_measurements)

# Count how many times each SNP is significant
# Initialize the count dataframe
frequency_of_child_spls <- data.frame(matrix(0, nrow = nrow(children_spls_results[[1]]), ncol = ncol(children_spls_results[[1]])))
rownames(frequency_of_child_spls) <- rownames(children_spls_results[[1]])
colnames(frequency_of_child_spls) <- c("BMI Z-Score", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio")

frequency_of_parent_spls <- data.frame(matrix(0, nrow = nrow(parents_spls_results[[1]]), ncol = ncol(parents_spls_results[[1]])))
rownames(frequency_of_parent_spls) <- rownames(parents_spls_results[[1]])
colnames(frequency_of_parent_spls) <- c("BMI", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio")

# Update the count dataframe
for (df in children_spls_results) {
  frequency_of_child_spls <- frequency_of_child_spls + (df != 0)
}

for (df in parents_spls_results) {
  frequency_of_parent_spls <- frequency_of_parent_spls + (df != 0)
}
#frequency_of_child_spls
#frequency_of_parent_spls

# Save results
write.csv(children_spls_results, file = "Children_SPLS_Results.csv", row.names = TRUE)
write.csv(parents_spls_results, file = "Parent_SPLS_Results.csv", row.names = TRUE)
write.csv(frequency_of_child_spls, file = "SPLS_Child_Frequency.csv", row.names = TRUE)
write.csv(frequency_of_parent_spls, file = "SPLS_Parent_Frequency.csv", row.names = TRUE)

# Creating a table figure with kable
frequency_of_child_spls %>%
  kable(caption = "Frequency of significant factors with SPLS") %>% 
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE, width = "2em") %>%
  column_spec(2:9, width = "6em") %>%
  as_image(., width = 7, file="~/Desktop/Child_SPLS_frequency_table.png", zoom=7)

# Creating a table figure with kable
frequency_of_parent_spls %>%
  kable(caption = "Frequency of significant factors with SPLS") %>% 
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(0,bold=TRUE) %>%
  column_spec(1, bold = TRUE, width = "2em") %>%
  column_spec(2:9, width = "6em") %>%
  as_image(., width = 7, file="~/Desktop/Parent_SPLS_frequency_table.png", zoom=7)


#################################
# Create Heatmap
#################################
# Code created with Nishita Sharif to keep plots consistent throughout both our research

# Creating a heatmap with the frequency that a SNP was found to be significant when doing SPLS analysis
frequency_of_child_spls_long <- frequency_of_child_spls %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

frequency_of_parent_spls_long <- frequency_of_parent_spls %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

frequency_of_child_spls_long$order<-factor(frequency_of_child_spls_long$colname, levels=c("BMI Z-Score", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio"))
frequency_of_parent_spls_long$order<-factor(frequency_of_parent_spls_long$colname, levels=c("BMI", "Body Mass", "Height", "Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake", "Waist Circumference to Height Ratio"))


child_spls_freq_heatmap <- ggplot(frequency_of_child_spls_long, 
                          aes(x = order, # Anthropocentric measure
                              y = rowname, # rsIDs
                              fill = value, # Coefficient
                              label = round(value, 2))) + 
  geom_tile(colour = "black") + # Add grid lines between tiles
  geom_text(colour = "black", # Add labels for the values
            size = 2, 
            hjust = 0.5, 
            vjust = 0.5,
            data = subset(frequency_of_child_spls_long, value >= 1) # Can remove this line if you want the 0s to appear in the tiles
  ) +  
  scale_fill_gradient2(low = "white",high = "darkorange", limits = c(0, 10)) +
  labs(y = "SNP rsID", x = "Anthropometric Measurement") +
  ggtitle("SPLS Child") +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(angle = 25, hjust = 1, size = 5),
    axis.title.y = element_text(size = 8),
    plot.title = element_text( size = 12, hjust = 0.5)
  ) +
  scale_x_discrete(expand=c(0,0),labels = function(x) stringr::str_wrap(x, width = 15))+
  scale_y_discrete(expand=c(0,0))

parent_spls_freq_heatmap <- ggplot(frequency_of_parent_spls_long, 
                                  aes(x = order, # Anthropocentric measure
                                      y = rowname, # rsIDs
                                      fill = value, # Coefficient
                                      label = round(value, 2))) + 
  geom_tile(colour = "black") + # Add grid lines between tiles
  geom_text(colour = "black", # Add labels for the values
            size = 2, 
            hjust = 0.5, 
            vjust = 0.5,
            data = subset(frequency_of_parent_spls_long, value >= 1) # Can remove this line if you want the 0s to appear in the tiles
  ) +  
  scale_fill_gradient2(low = "white",high = "darkorange", limits = c(0, 10)) +
  labs(y = "SNP rsID", x = "Anthropometric Measurement") +
  ggtitle("SPLS Parent") +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(angle = 25, hjust = 1, size = 5),
    axis.title.y = element_text(size = 8),
    plot.title = element_text( size = 12, hjust = 0.5)
  ) +
  scale_x_discrete(expand=c(0,0),labels = function(x) stringr::str_wrap(x, width = 15))+
  scale_y_discrete(expand=c(0,0))

child_spls_freq_heatmap
parent_spls_freq_heatmap

output_file <- paste0("~/Desktop/child_spls_freq_heatmap.png")
output_file_parent <- paste0("~/Desktop/parent_spls_freq_heatmap.png")

ggsave(filename = output_file, 
       child_spls_freq_heatmap,
       dpi = 600, width = 6, height = 4)
ggsave(filename = output_file_parent, 
       parent_spls_freq_heatmap,
       dpi = 600, width = 6, height = 4)


