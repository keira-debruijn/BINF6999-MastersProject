setwd("~/Desktop")

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
library(knitr)
library(forecast)
library(bestNormalize)
library(purrr)

##################################################
# Functions
##################################################

# Function to filter and clean data
filter_and_clean_anthropometric_data<-function(group){
  SNP_Alleles <- read.csv("SNP_Alleles.csv", header=TRUE)
  major_allele_freq_df <- data.frame(lapply(SNP_Alleles, function(x) ifelse(x == "--", NA, x)))
  
  major_allele_freq_df$X <- gsub("O", "0", major_allele_freq_df$X)
  rownames(major_allele_freq_df)<-major_allele_freq_df$X
  major_allele_freq_df$X <- NULL
  
  
  if (group=="children"){
    major_allele_frequencies<-major_allele_freq_df %>% 
      filter(grepl("A|B|C", rownames(.)))
    anthropometric_data<-read_excel("Children_Study_Data.xlsx")
    
    
    anthropometric_data$age <- round(interval(
      anthropometric_data$dob,
      anthropometric_data$date_ha
    ) / years(1), 5) %>%
      as.numeric()
    anthropometric_data<-anthropometric_data%>%
      select(-dob, -date_ha)
    
    selected_columns <- c("pid","age", "sex", "bmi_z", "bm_kg", "ht_cm", "wc_cm", "SUGR", "ADD_SUGARS", "KCAL")
  }
  else{
    major_allele_frequencies<-major_allele_freq_df %>% 
      filter(grepl("P|S", rownames(.)))
    anthropometric_data<-read.csv("Parent_Study_Data.csv", header=TRUE) %>%
      { .[!grepl("pregnant", .$bm_kg, ignore.case = TRUE), ] }
    selected_columns <- c("pid", "age_years_ha", "sex", "bmi", "bm_kg", "ht_cm", "wc_cm", "SUGR", "ADD_SUGARS", "KCAL")
  }
  # Removing any columns that have only "0"s in the entire column
  filtered_major_allele_frequencies <- major_allele_frequencies %>%
    select(where(~ any(. != 0)))
  
  rownames_filtered<-rownames(filtered_major_allele_frequencies)
  
  # Apply as.character() to each cell (necessary for creating boxplots later)
  filtered_major_allele_frequencies <- as.data.frame(apply(filtered_major_allele_frequencies, 2, as.character), stringsAsFactors = FALSE)
  rownames(filtered_major_allele_frequencies)<-rownames_filtered
  
  
  p_IDs<-rownames(major_allele_frequencies)
  
  # Selecting for pids that are also found in the pids
  updated_anthropometric_data <- anthropometric_data[anthropometric_data$pid %in% p_IDs, ] 
  
  
  # Create a new dataframe with selected columns
  updated_anthropometric_data <- updated_anthropometric_data[, selected_columns] %>% 
    as.data.frame() %>%
    mutate_all(~ ifelse(. == "NaN" | . == "", NA, .)) %>%
    mutate(across(-c(pid, sex), ~ as.numeric(as.character(.)))) %>%
    mutate(wc_ht_ratio = wc_cm/ht_cm)
  
  # Changing rownames to the IDs (and removing pid column)
  rownames(updated_anthropometric_data)<-updated_anthropometric_data$pid
  updated_anthropometric_data<-updated_anthropometric_data[,-1]
  
  updated_anthropometric_data$wc_ht_ratio<-as.numeric(updated_anthropometric_data$wc_cm)/as.numeric(updated_anthropometric_data$ht_cm)
  
  
  if (group=="children"){
    updated_anthropometric_data <- updated_anthropometric_data %>%
      select(sex, bmi_z, bm_kg, ht_cm, wc_cm, SUGR, ADD_SUGARS, KCAL, age, wc_ht_ratio) %>%
      `colnames<-`(c("Sex","BMI Z-Score", "Body Mass", "Height","Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake","Age", "Waist Circumference-to-Height Ratio"))
    
  }
  
  else{
    
    updated_anthropometric_data <-  updated_anthropometric_data %>%
      select(sex, bmi, bm_kg, ht_cm, wc_cm, SUGR, ADD_SUGARS, KCAL, age_years_ha, wc_ht_ratio) %>%
      `colnames<-`(c("Sex","BMI", "Body Mass", "Height","Waist Circumference", "Total Sugar", "Added Sugar", "Caloric Intake","Age", "Waist Circumference-to-Height Ratio"))
  }
  
  
  return(updated_anthropometric_data)
  
  
}

# Function to replace outliers with NA
replace_outliers_with_na <- function(x) {
  # Calculate mean and standard deviation of the column
  mean_val <- mean(x, na.rm = TRUE)
  std_dev <- sd(x, na.rm = TRUE)
  
  # Set a threshold for outlier detection (e.g., 3 standard deviations)
  threshold <- 3 * std_dev
  
  # Replace outliers with NA
  x[abs(x - mean_val) > threshold] <- NA
  return(x)
}

# Function to normalize added and total sugar per 1000 calories consumed
normalize_total_sugar_added_sugar <- function(anthropometric_data_and_allele_freq_combined_df) {
  raw_data <- anthropometric_data_and_allele_freq_combined_df
  
  # Normalize "Total Sugar" and "Added Sugar" by the amount of calories consumed (per 1000 cal)
  if ("Caloric Intake" %in% colnames(raw_data)) {
    total_sugar_normalized <- (raw_data[["Total Sugar"]] / raw_data[["Caloric Intake"]]) * 1000
    added_sugar_normalized <- (raw_data[["Added Sugar"]] / raw_data[["Caloric Intake"]]) * 1000
    raw_data[["Total Sugar"]] <- total_sugar_normalized
    raw_data[["Added Sugar"]] <- added_sugar_normalized
  } else {
    # If "Caloric Intake" column is not present, return the original data
    print("Warning: 'Caloric Intake' column not found. 'Total Sugar' and 'Added Sugar' will not be normalized.")
  }
  
  return(raw_data)
}

# Function to normalize all data
normalize_data<-function(anthropometric_data_and_allele_freq_combined_df){
  raw_data<-anthropometric_data_and_allele_freq_combined_df[2:10]
  
  outliers_removed <- raw_data %>%
    map_df(~replace_outliers_with_na(.))
  
  
  list_of_all_transformed<-list()
  list_chosen_transformations<-list()
  for (i in 1:ncol(outliers_removed)){
    # if (i==1) {next}
    # if (i==3) {next}
    
    print(colnames(outliers_removed)[i])
    print(class(outliers_removed)[i])
    # if (i==3){break}
    
    data<-outliers_removed[[i]]
    
    normalized<-bestNormalize(data, r=100)
    
    # df_blah<-cbind(df_blah,normalized)
    list_of_all_transformed[[i]]<-normalized$x.t
    list_chosen_transformations[[i]]<-normalized$chosen_transform
    # print(normalized)
    print(i)
    
    
  }
  
  transformed_df<-data.frame(list_of_all_transformed) %>%
    `colnames<-`(colnames(anthropometric_data_and_allele_freq_combined_df)[2:10])
  # plot_distributions(transformed_df)
  transformed_df_normality<-create_normality_df(transformed_df)
  
  transformed_df_normality$col_name[transformed_df_normality$normal==""] %>% # Things that aren't normally distributed
    print()
  print(transformed_df_normality)
  print(list_chosen_transformations)
  anthropometric_data_and_allele_freq_combined_df[2:10]<-transformed_df
  return(anthropometric_data_and_allele_freq_combined_df)
  
  
}

# Function to create Shipiro-Wilks normality data frame
create_normality_df<-function(anthropometric_data_df){
  
  # Initialize an empty dataframe to store the results
  normality_results <- data.frame(col_name = character(),
                                  p_value = numeric(),
                                  normal = character(),
                                  stringsAsFactors = FALSE)
  
  # Iterate through each column in the dataframe
  for (col_name in colnames(anthropometric_data_df)) {
    # Check if the column is numeric
    if (is.numeric(anthropometric_data_df[[col_name]])) {
      # Perform Shapiro-Wilk test on the numeric column to determine if the data is normally distributed
      shapiro_result <- shapiro.test(anthropometric_data_df[[col_name]])
      
      formatted_p_value <- format(shapiro_result$p.value, scientific = TRUE, digits = 3)
      
      # Create a new row for the result dataframe
      new_row <- data.frame(col_name = col_name,
                            p_value = formatted_p_value,
                            normal = ifelse(shapiro_result$p.value <= 0.05, "", "Normal"),
                            stringsAsFactors = FALSE)
      
      # Append the new row to the result dataframe
      normality_results <- rbind(normality_results, new_row)
    }
    # rm(col_name)
  }
  # rm(new_row,formatted_p_value,shapiro_result)
  return(normality_results)
}

# Function to plot all the distributions
plot_distributions<-function(anthropometric_data_df){
  if (grepl("parent", deparse(substitute(anthropometric_data_df)), ignore.case = TRUE)) {
    variable_x_axes<-c("Sex","BMI", "Body Mass (kg)", "Height (cm)","Waist Circumference (cm)", "Total Sugar (mg)", "Added Sugar (mg)", "Calories (KCal)","Age", "Waist Circumference-to-Height Ratio")
    fill_colour<-"#F0E442"
  }
  else{
    variable_x_axes<-c("Sex","BMI Z-Score", "Body Mass (kg)", "Height (cm)","Waist Circumference (cm)", "Total Sugar (mg)", "Added Sugar (mg)", "Calories (KCal)","Age", "Waist Circumference-to-Height Ratio")
    fill_colour<-"#0072B2"
  }
  
  distribution_plots <- list()
  for (i in 1:length(colnames(anthropometric_data_df))) {
    variable=colnames(anthropometric_data_df)[i]
    
    if (variable=="Sex"){
      next
    }
    plot <- ggplot(anthropometric_data_df, aes(x = !!as.name(variable))) +
      geom_histogram(color = "black", fill = fill_colour, bins=25, size = 0.1) +
      ylab("Frequency") +
      xlab(paste0(variable_x_axes[i]))+
      ggtitle(paste0(variable))+
      theme(
        plot.title = element_text(size=5, hjust = 0.5, vjust=10,margin = margin(b = -8)),
        axis.text.x = element_text(size = 3),  # Set the size of x-axis labels to 8
        axis.text.y = element_text(size = 3),
        axis.title.x = element_text(size = 4),  # Set the size of x-axis title to 8
        axis.title.y = element_text(size = 4),
        # plot.margin = margin(t = 0.1),
        axis.ticks = element_line(size = 0.1),
        panel.grid = element_line(size = 0.1)
      ) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0))
    distribution_plots <- c(distribution_plots, list(plot))
    
  }
  
  distribution_plots_figure<-grid.arrange(distribution_plots[[1]],distribution_plots[[2]],
                                          distribution_plots[[3]],distribution_plots[[4]],
                                          distribution_plots[[5]],distribution_plots[[6]],
                                          distribution_plots[[7]],distribution_plots[[8]],
                                          distribution_plots[[9]],
                                          ncol=3)
  
  
  print(distribution_plots_figure)
  return(distribution_plots_figure)
}

# Function to plot all the correlations
create_correlations_plot<-function(anthropometric_data_df){
  
  correlations<-cor(anthropometric_data_df[,c(2:10)],use = "complete.obs")
  variable_x_axes<-c("Body Mass", "Height", "Waist Circumference", "Sugar Total", "Added Sugar", "Caloric Intake", "Age", "Waist Circumference-to-Height Ratio")
  if (grepl("children", deparse(substitute(anthropometric_data_df)), ignore.case = TRUE)) {
    
    variable_y_axes<-c("BMI Z-Score", "Body Mass", "Height", "Waist Circumference", "Sugar Total", "Added Sugar", "Caloric Intake", "Age")
    
  }
  else{
    
    variable_y_axes<-c("BMI", "Body Mass", "Height", "Waist Circumference", "Sugar Total", "Added Sugar", "Caloric Intake", "Age")
    
  }
  
  ## Correlation plot
  # "#0072B2", "#D55E00"
  corr_plot<-ggcorrplot(correlations, lab = TRUE, type = "lower", colors = c("#0072B2", "white", "#D55E00"),) +
    # labs(title = "Correlations Between Anthropometric Measures") +
    theme(panel.grid.major=element_blank(),
          axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 8, angle = 20, hjust = 1),
          plot.background = element_rect(fill = "white", color =NA)) +
    scale_x_discrete(labels=str_wrap(variable_x_axes, width=16)) +
    scale_y_discrete(labels=str_wrap(variable_y_axes,width=16)) +
    geom_hline(yintercept = seq(0.5, 9, by = 1), color = "gray90", linetype = "solid") +
    geom_vline(xintercept = seq(0.5, 9, by = 1), color = "gray90", linetype = "solid")
  
  # corr_plot <- corr_plot +
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "green", midpoint = 0, na.value = "white")
  
  
  return(corr_plot)
  
}

# Function to combine the anthropomteric and allelic data frames
combine_anthropometric_and_allelic_data<-function(anthropometric_data_df){
  major_allele_freq_df<-read.csv("SNP_Alleles.csv", header=TRUE)
  major_allele_freq_df <- data.frame(lapply(major_allele_freq_df, function(x) ifelse(x == "--", NA, x)))
  
  major_allele_freq_df$X <- gsub("O", "0", major_allele_freq_df$X)
  rownames(major_allele_freq_df)<-major_allele_freq_df$X
  major_allele_freq_df$X <- NULL
  
  anthropometric_data_and_allele_freq_combined <- merge(anthropometric_data_df,major_allele_freq_df, by = "row.names", all = FALSE)
  
  rownames(anthropometric_data_and_allele_freq_combined)<-anthropometric_data_and_allele_freq_combined$Row.names
  anthropometric_data_and_allele_freq_combined<-anthropometric_data_and_allele_freq_combined %>%
    select(-Row.names)
  
  anthropometric_data_and_allele_freq_combined <- anthropometric_data_and_allele_freq_combined %>%
    mutate_at(vars(starts_with("rs")), ~ ifelse(is.na(.), 1, .))
  
  return(anthropometric_data_and_allele_freq_combined)
  
}

# Creates the demographics table
demographics_table<-function(anthropometrics_data_df){
  anthropometrics_data_df<- anthropometrics_data_df %>%
    select(-contains("rs"))
  num_total<-dim(anthropometrics_data_df)[1]
  
  male_df <- subset(anthropometrics_data_df, Sex == "M") %>%
    select_if(is.numeric)
  female_df <- subset(anthropometrics_data_df, Sex == "F")    %>%
    select_if(is.numeric)
  
  num_male<-dim(male_df)[1]
  num_female<-dim(female_df)[1]
  
  ## Male ##
  # Calculate the mean and standard deviation for each column
  mean_values <- colMeans(male_df, na.rm = TRUE) %>% 
    round(., 1)
  sd_values <- apply(male_df, 2, sd, na.rm = TRUE) %>% 
    round(., 1)
  
  # Create a new dataframe with mean ± sd format (rounded to 1 decimal point)
  male_summary_df <- data.frame(Characteristic = colnames(male_df),
                                M = paste(mean_values, "±", sd_values)) %>%
    rbind(c("N",num_male),.)
  
  
  ## Female ##
  # Calculate the mean and standard deviation for each column
  mean_values <- colMeans(female_df, na.rm = TRUE) %>% 
    round(., 1)
  sd_values <- apply(female_df, 2, sd, na.rm = TRUE) %>% 
    round(., 1)
  
  
  # Create a new dataframe with mean ± sd format (rounded to 1 decimal point)
  female_summary_df <- data.frame(Characteristic = colnames(female_df),
                                  F = paste(mean_values, "±", sd_values)) %>%
    rbind(c("N",num_female),.)
  
  
  ## Total ##
  
  anthropometrics_data_df<-anthropometrics_data_df %>%
    select_if(is.numeric)
  # Calculate the mean and standard deviation for each column
  mean_values <- colMeans(anthropometrics_data_df, na.rm = TRUE) %>% 
    round(., 1)
  sd_values <- apply(anthropometrics_data_df, 2, sd, na.rm = TRUE) %>% 
    round(., 1)
  
  
  # Create a new dataframe with mean ± sd format (rounded to 1 decimal point)
  total_summary_df <- data.frame(bleh = colnames(anthropometrics_data_df),
                                 Total = paste(mean_values, "±", sd_values)) %>%
    rbind(c("N",num_total),.)
  
  combined_df<- cbind(male_summary_df,female_summary_df,total_summary_df) %>%   select(-3, -5)
  return(combined_df)
}

# Function to get QQ Plots
plot_qq_plots<-function(df) {
  all_qqplots<-list()
  # Columns 2 to 11
  selected_columns <- 2:10 
  
  # Loop through each column and create the QQ plot
  for (column_index in selected_columns) {
    column_name <- names(df)[column_index]
    data <- df[[column_index]]
    
    # Remove NA values from the data
    data <- data[!is.na(data)]
    
    p <- ggplot(data.frame(Quantiles = qqnorm(data, plot.it = FALSE)$x, Values = qqnorm(data, plot.it = FALSE)$y),
                aes(x = Quantiles, y = Values)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(title = paste("QQ Plot of", column_name), x = "Theoretical Quantiles", y = "Sample Quantiles")
    all_qqplots[[column_index - 1]] <- p
    
  }
  
  
  qqplot_figure <- grid.arrange(
    grobs = all_qqplots,
    ncol = 3)
  return(qqplot_figure)
  
}


##################################################


##################################################
# Children
##################################################

#Run and save all the functions previously created
children_anthropometric_data <- filter_and_clean_anthropometric_data("children")
children_anthropometric_data_sug <- normalize_total_sugar_added_sugar(children_anthropometric_data)
children_normality_df<-create_normality_df(children_anthropometric_data_sug)
children_distribution_figure<-plot_distributions(children_anthropometric_data_sug)

children_combined_anthropometric_and_allelic_data<-combine_anthropometric_and_allelic_data(children_anthropometric_data_sug)
children_correlation_plot<-create_correlations_plot(children_anthropometric_data_sug)
children_demopgraphics_table<-demographics_table(children_combined_anthropometric_and_allelic_data)

set.seed(5)
normalized_children_anthropometric_data<-normalize_data(children_anthropometric_data_sug)
children_normality_df<-create_normality_df(normalized_children_anthropometric_data) 
children_distribution_figure<-plot_distributions(normalized_children_anthropometric_data)
children_qqplots<-plot_qq_plots(normalized_children_anthropometric_data)
children_combined_anthropometric_and_allelic_data_scaled<-combine_anthropometric_and_allelic_data(normalized_children_anthropometric_data)

# Tried to manually normalize height, still did not normalize, luckily height being normalized is not critical in this analysis

# Saves
write.csv(children_combined_anthropometric_and_allelic_data,"~/Desktop/children_anthropometric_data_and_allele_freq_combined.csv")
write.csv(children_combined_anthropometric_and_allelic_data_scaled,"~/Desktop/anthropometric_data_and_allele_freq_combined_scaled.csv")

##################################################
# Parents
##################################################
#Run and save all the functions previously created
parent_anthropometric_data<-filter_and_clean_anthropometric_data("parent")
parent_anthropometric_data_sug <- normalize_total_sugar_added_sugar(parent_anthropometric_data)
parent_normality_df<-create_normality_df(parent_anthropometric_data_sug)
parent_distribution_figure<-plot_distributions(parent_anthropometric_data_sug)

parent_correlation_plot<-create_correlations_plot(parent_anthropometric_data_sug)
parent_combined_anthropometric_and_allelic_data<-combine_anthropometric_and_allelic_data(parent_anthropometric_data_sug)
parent_demopgraphics_table<-demographics_table(parent_combined_anthropometric_and_allelic_data)

set.seed(14)
normalized_parent_anthropometric_data<-normalize_data(parent_anthropometric_data_sug)
parent_normality_df<-create_normality_df(normalized_parent_anthropometric_data) 
parent_distribution_figure<-plot_distributions(normalized_parent_anthropometric_data)
parent_qqplots<-plot_qq_plots(normalized_parent_anthropometric_data)
parent_combined_anthropometric_and_allelic_data_scaled<-combine_anthropometric_and_allelic_data(normalized_parent_anthropometric_data)

# transform bmi manually
bmi_transformed <- bestNormalize(parent_anthropometric_data_sug$BMI, r = 100) %>%
  .$x.t
hist(bmi_transformed)
shapiro.test(bmi_transformed)
normalized_parent_anthropometric_data$BMI <- bmi_transformed

# Saves
write.csv(parent_combined_anthropometric_and_allelic_data,"~/Desktop/parent_anthropometric_data_and_allele_freq_combined.csv")
write.csv(parent_combined_anthropometric_and_allelic_data_scaled,"~/Desktop/p_anthropometric_data_and_allele_freq_combined_scaled.csv")

##################################################
#      Creating Demographics Tables
##################################################
children_demopgraphics_table%>%
  kable(format = "html",escape = FALSE, caption = "Children Demographics")%>%
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  column_spec (1, border_right = T, width = "15em") %>%
  row_spec(0, bold = TRUE) %>%
  as_image(., width = 7, file="~/Desktop/children_demographics_table.png", zoom=15)

parent_demopgraphics_table%>%
  kable(format = "html",escape = FALSE, caption = "Parent Demographics")%>%
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  column_spec (1, border_right = T, width = "15em") %>%
  row_spec(0, bold = TRUE) %>%
  as_image(., width = 7, file="~/Desktop/parent_demographics_table.png", zoom=15)

##################################################
# Combining Distribution figures
##################################################

# Create the grid arrangement with extra white space between plots
combined_anthropometric_distirbution <- grid.arrange(
  arrangeGrob(children_distribution_figure, top = "Children"),
  arrangeGrob(parent_distribution_figure, top = "Parents"),
  ncol = 2,
  top = textGrob("Anthropometric Distribution", gp = gpar(fontsize = 20, font = 2)))

ggsave(filename = "~/Desktop/All_Measures_Distribution.png",
       combined_anthropometric_distirbution,
       dpi = 900, 
       width = 16, height = 13,
       units = c("in"))

##################################################
# Combining Normality Tables
##################################################

combined_normality_df<-cbind(children_normality_df,parent_normality_df) %>%
  setNames(NULL) %>%
  rbind(c("Measure", "P-Value", "Normality","Measure", "P-Value", "Normality"),.)

combined_normality_df %>% 
  kable(format = "html",escape = FALSE, caption = "Shapiro-Wilk's Test Results")%>%
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  column_spec (3, border_right = T) %>%
  row_spec(1, bold = TRUE) %>%
  row_spec(1, extra_css = "border-bottom: 1px solid")  %>%
  add_header_above(., c("Children" = 3, "Parents" = 3)) %>%
  as_image(., width = 7, file="~/Desktop/Normality_Tables.png", zoom=15)

##################################################
# Combining Q-Q Plots
##################################################

combined_correlation_plots <- grid.arrange(
  arrangeGrob(children_qqplots, top = "Children"),
  arrangeGrob(parent_qqplots, top = "Parents"),
  ncol = 2,
  top = textGrob("Correlations Between Anthropometric Measures", gp = gpar(fontsize = 20, font = 2))
)

ggsave(filename = "~/Desktop/QQ_Plots.png",
       combined_correlation_plots,
       dpi = 900, 
       width = 17, height = 6,
       units = c("in"))

##################################################
# Combining Correlation Plots
##################################################

children_correlation_plot<-create_correlations_plot(children_anthropometric_data_sug)
parent_correlation_plot<-create_correlations_plot(parent_anthropometric_data_sug)

combined_correlation_plots <- grid.arrange(
  arrangeGrob(children_correlation_plot, top = "Children"),
  arrangeGrob(parent_correlation_plot, top = "Parents"),
  ncol = 2,
  top = textGrob("Correlations Between Anthropometric Measures", gp = gpar(fontsize = 20, font = 2))
)

ggsave(filename = "~/Desktop/Correlation_Plots.png",
       combined_correlation_plots,
       dpi = 900, 
       width = 12, height = 6,
       units = c("in"))

##################################################
# Family, Income, and education demographics
##################################################

pids<-c(rownames(children_anthropometric_data_sug),rownames(parent_anthropometric_data_sug))
length(pids) #470

all_fids<-substr(pids, start = nchar(pids) - 2, stop = nchar(pids))

num_families<-length(unique(all_fids))

num_children<-pids %>%
  as.data.frame() %>%
  filter(str_detect(., "A|B|C")) %>%
  nrow() 

num_parents<-pids %>%
  as.data.frame() %>%
  filter(str_detect(., "P|S")) %>%
  nrow()

female_count <- sum(children_combined_anthropometric_and_allelic_data$Sex == 'F')
male_count <- sum(children_combined_anthropometric_and_allelic_data$Sex == 'M')
p_female_count <- sum(parent_combined_anthropometric_and_allelic_data$Sex == "F", na.rm = TRUE)
p_male_count <- sum(parent_combined_anthropometric_and_allelic_data$Sex == "M", na.rm = TRUE)

num_families # 168
num_children # 209
num_parents # 261
female_count #113
male_count #96
p_female_count #154
p_male_count #106

children_fid <- substr(rownames(children_anthropometric_data_sug), start = nchar(rownames(children_anthropometric_data_sug)) - 2, stop = nchar(rownames(children_anthropometric_data_sug))) 

parent_fid <- substr(rownames(parent_anthropometric_data_sug), start = nchar(rownames(parent_anthropometric_data_sug)) - 2, stop = nchar(rownames(parent_anthropometric_data_sug)))

not_in_children <- parent_fid[!parent_fid %in% children_fid]
not_in_parent <- children_fid[!children_fid %in% parent_fid]
not_in_children # None
not_in_parent # 471, 703, 703

# Finding number of dyads
children_selected <- table(children_fid) %>%
  {names(.)[. == 1]} # Select items that occur once
children_selected # Display the selected items
parent_selected <- table(parent_fid) %>%
  {names(.)[. == 1]}# Select items that occur once
parent_selected  # Display the selected items

num_dyads <- length(intersect(children_selected, parent_selected))
num_dyads # 57

# Finding number of triads
children_selected <- table(children_fid) %>%
  {names(.)[. == 2]} # Select items that occur once
children_selected # Display the selected items
parent_selected <- table(parent_fid) %>%
  {names(.)[. == 1]}# Select items that occur once
parent_selected # Display the selected items

num_triads <- length(intersect(children_selected, parent_selected))
num_triads #14
#two-parents or three-siblings???

income_household_summary <- read.csv("Parent_Study_Data.csv", header = TRUE) %>%
  filter(fid %in% parent_fid) %>%
  select(fid, income_household) %>%
  `colnames<-`(c("fid","group")) %>%
  distinct(fid, .keep_all = TRUE) %>%
  arrange(fid) %>%
  group_by(group) %>%
  summarize(count = n()) %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  as.data.frame() %>%
  .[match(c("$10,000 to $19,999", "$20,000 to $29,999", "$30,000 to $39,999", "$40,000 to $49,999",
            "$50,000 to $59,999", "$60,000 to $69,999", "$70,000 to $79,999", "$80,000 to $89,999",
            "$90,000 to $99,999", "$100,000 to $149,999", "$150,000 or more", "I am not comfortable answering this question", "I don’t know"),
          .$group), ] %>%
  `rownames<-`(NULL)

income_household_summary<-data.frame(
  group=c("< $50 000", 
          "$50 000-79 999", 
          "$80 000-100 000", 
          ">$100 000", 
          "Don't know or am not comfortable answering"),
  count=c(sum(income_household_summary$count[1:4]),
          sum(income_household_summary$count[5:7]),
          sum(income_household_summary$count[8:9]),
          sum(income_household_summary$count[10:11]),
          sum(income_household_summary$count[12:13])),
  percentage=c(sum(income_household_summary$percentage[1:4]),
               sum(income_household_summary$percentage[5:7]),
               sum(income_household_summary$percentage[8:9]),
               sum(income_household_summary$percentage[10:11]),
               sum(income_household_summary$percentage[12:13]))
)

income_household_summary

df <- c("Don't know / < $50 000", 32, 19.27)
df <- t(df)
colnames(df) <- c("group", "count", "percentage")
income_household_summary <- rbind(df, income_household_summary)
income_household_summary <- income_household_summary[-c(2,6),]


# Read the data and get the complete parent_education_summary
parent_education_summary <-
  read.csv("Parent_Study_Data.csv", header = TRUE) %>%
  # filter(fid %in% parent_fid) %>%
  filter(pid %in% rownames(parent_anthropometric_data)) %>%
  select(fid, school_completed) %>%
  `colnames<-`(c("fid","group")) %>%
  arrange(fid) %>%
  group_by(group) %>%
  summarize(count = n()) %>%
  # mutate(proportion = count / sum(count)) %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  as.data.frame() %>%
  .[match(c("Some high school", "High School graduate", "Some college or technical school", "Some university", "College graduate", "University graduate", "Postgraduate training or degree", ""), .$group), ] %>%
  mutate(group = ifelse(group == "", NA, group)) %>%
  `rownames<-`(NULL)

parent_education_summary

df2 <- c("High School or less / Some college or technical school", 30, 11.5)
df2 <- t(df2)
colnames(df2) <- c("group", "count", "percentage")
parent_education_summary <- rbind(df2, parent_education_summary)
parent_education_summary <- parent_education_summary[-c(2,3,4,9),]

df3 <- c("Some university / University graduate", 94, 36.01)
df3 <- t(df3)
colnames(df3) <- c("group", "count", "percentage")
parent_education_summary <- rbind(df3, parent_education_summary)
parent_education_summary <- parent_education_summary[-c(3,5),]

# Manually specify the desired order of rows using base row numbers
desired_order <- c(2,3,1,4)
# Add row numbers to the dataframe
parent_education_summary <- parent_education_summary %>%
  mutate(row_num = row_number())
# Sort the rows based on the desired order using base row numbers
parent_education_summary <- parent_education_summary %>%
  arrange(match(row_num, desired_order))
# Remove the row_num column
parent_education_summary <- parent_education_summary %>%
  select(-row_num)


education_and_income_df<-rbind(c("Group"," ","N"),
                               c("Parents"," ",num_parents),
                               c("Male"," ",p_male_count),
                               c("Female"," ",p_female_count),
                               c("Children"," ",num_children),
                               c("Male"," ",male_count),
                               c("Female"," ",female_count),
                               c("Families"," ",num_families),
                               c("Dyads"," ",num_dyads),
                               c("Level of Education", "Count", "Percentage"),
                               parent_education_summary,
                               c("Income", "Count", "Percentage"),
                               income_household_summary) %>%
  select(group, count, percentage)


print(education_and_income_df)

# Just education and income data
education_and_income_df %>%
  select(group, count, percentage) %>%
  mutate(group = cell_spec(group, italic = ifelse(row_number() %in% c(2:9,11:14,16:19), T, F))) %>% 
  kable(booktabs = TRUE, escape = FALSE, row.names = FALSE, col.names = NULL) %>%
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(column=1,width="22em") %>%
  row_spec(row=c(1,10,15), bold=TRUE)  %>%
  group_rows("Group", 1, 9) %>% # change this to the number of rows of parents
  group_rows("Parent Education", 10, 14) %>% # change this to the number of rows in the education df
  group_rows("Household Income per Family", 15, 19) %>%# change this to the number of rows in the income df
  add_indent(positions = c(2:9,11:14,16:19)) %>%
  add_indent(positions = c(3:4,6:7), level_of_indent = 2) #%>%
  #as_image(., width = 4, file="~/Desktop/Table1_demographics.png", zoom=15)

# Get demographics table with all anthropometric data
children_demopgraphics_table<-demographics_table(children_combined_anthropometric_and_allelic_data)
parent_demopgraphics_table<-demographics_table(parent_combined_anthropometric_and_allelic_data)

parent_demopgraphics_table<-rbind(c("Characteristic", "M", "F", "Total"),parent_demopgraphics_table)
children_demopgraphics_table<-rbind(c("Characteristic", "M", "F", "Total"),children_demopgraphics_table)

total_demographics<-rbind(parent_demopgraphics_table,children_demopgraphics_table)

total_demographics %>%
  mutate(Characteristic = cell_spec(Characteristic, italic = ifelse(row_number() %in% c(2:12,14:24), T, F))) %>%
  kable(booktabs = TRUE, escape = FALSE, row.names = FALSE, col.names = NULL) %>%
  kable_classic() %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(column=1,width="22em") %>%
  row_spec(row=c(1,13), bold=TRUE)%>% 
  group_rows("Parent Demographics", 1, 11) %>% # change this to the rows of parent demographics
  group_rows("Children Demographics", 12, 22) %>%
  add_indent(positions = c(2:11,13:22)) %>%
  as_image(., width = 4, file="~/Desktop/Table2_demographics.png", zoom=15)

