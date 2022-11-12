### functions to compute bat counts
# written by Allan Roberts on Mar 29, 2017
# modified by Joanna Burgar on Nov 11, 2022 
# now works as functions sourced by "Predicting_ABSpecies.R"


# Define functions for computing the category with the highest output probability, and the second and third highest probability. 
# In the case of ties, an NA value is returned:  

get_most_probable <- function(species, rank){
	species <- as.character(species)
	mp <- species[rank == 1];
	if (length(mp) != 1) {mp <- NA};
	mp <- as.character(mp)
	return(mp)
}

get_second_probable <- function(species, rank){
	species <- as.character(species)
	mp <- species[rank == 2];
	if (length(mp) != 1) {mp <- NA};
	mp <- as.character(mp)
	return(mp)
}

get_third_probable <- function(species, rank){
	species <- as.character(species)
	mp <- species[rank == 3];
	if (length(mp) != 1) {mp <- NA};
	mp <- as.character(mp)
	return(mp)
}

### Function to create Alberta eBat output from random forest model input
compute_bat_counts_fn <- function(bat_data=bat_data){

  # Add an index for the row number:
  bat_data$row <- 1:nrow(bat_data)
  
  # Convert the `Date` column to character strings:
  bat_data$Date <- as.character(bat_data$Date)
  
  # Add a column for the date in POSIX format:
  bat_data$Date_POSIX <- as.POSIXct(strptime(bat_data$Date, format = "%Y%m%d"))

  # Add a new column that pastes together the site name and the file name; this is intended to help provide unique identifiers: 
  bat_data$Site_Filename <- paste(bat_data$Site, bat_data$Filename)
  
  # Add a column indicating the number of rows for each call sequence. Rows are counted for each file name, within each site:
  bat_data <- mutate(group_by(bat_data, Site_Filename), n_calls = length(Site_Filename))
  
  ###Counting Species by Highest Average Probability  
  # Convert the data frames to a long format:  
  bat_data_long <- melt(bat_data, variable.name = "species", value.name = "probability", id.vars = c("Site", "Date", "Filename", "Site_Filename", "row", "Date_POSIX", "n_calls"))
  
  bat_data_long$species <- as.character(bat_data_long$species)
  
  # Get the average probability for each species (averaged across calls) within each call sequence (i.e. `Filename` within `Site`):
  bat_data_mean_probabilities <- summarize(group_by(bat_data_long, Site, Date, Date_POSIX, Filename, species, n_calls), probability = mean(probability))
  
  # Rank the mean probabilities for each category in each call sequence, from highest to lowest, with the highest probability having a rank of 1:
  bat_data_mean_probabilities <- mutate(group_by(bat_data_mean_probabilities, Site, Date, Date_POSIX, Filename),  rank = rank(1 - probability))
  
  # Get the species with the most weight within each category. When species are not combined into aggregated categories, the value in the "category" column will generally be the same as the value in the "sp1" column, or "unknown" if the "prob1" value, corresponding to the category given by "sp1" is not above the required threshold given by `threshold_bat`. If the number of rows in a call sequence is less than 3, the category column is set to "unknown".  
  bat_data_most_weight_df <- summarize(group_by(bat_data_mean_probabilities, Site, Date, Date_POSIX, Filename, n_calls), prob1 = max(probability), sp1 = get_most_probable(species, rank), prob2 = rev(sort(probability))[2], sp2 = get_second_probable(species, rank), prob3 = rev(sort(probability))[3], sp3 = get_third_probable(species, rank))
  
  bat_data_output <- bat_data_most_weight_df
  
  bat_data_output <- rename(bat_data_output, "Site_Combined" = Site)
  
  bat_data_output$Site <- word(str_replace(bat_data_output$Site_Combined, "/_", " "), 1)
  
  bat_data_output$Detector <- word(str_replace(bat_data_output$Site_Combined, "/_", " "), 2)
  bat_data_output$Detector <- case_when(is.na(bat_data_output$Detector) ~ bat_data_output$Site,
                                        TRUE ~ as.character(bat_data_output$Detector))
  
  ###Setting Thresholds for the Classification of Bats and Noise  
  # Now thresholds are set for the positive identification of noise, and of bat species. These thresholds are used to determine if the probability for the category is high enough for a classification to be made.  
  
  threshold_noise <- 0.8
  threshold_bat   <- 0.5
  
  # Start with all call sequences categorized as "unknown":  
  bat_data_most_weight_df$category <- "unknown"
  
  # In the first phase of categorization, all call sequences that have a `sp1` value of `noise`, and a `prob1` value over the threshold for the positive identification of noise (here, 0.80) are categorized as noise:  
  # Index cases to be categorized as noise:
  index_noise <- with(bat_data_most_weight_df, (sp1 == "noise") & (prob1 > threshold_noise))
  
  bat_data_most_weight_df %>% filter(n_calls>2)
  # Set the indexed categories to "noise":
  bat_data_most_weight_df$category[index_noise] <- "noise"
  
  # If sp1 == "noise" and `prob1` is below the threshold value for noise, filter out the noise value, and standardize the remaining probabilities:
  # Index the noise values to be filtered out:
  index_remove_noise <- with(bat_data_most_weight_df, (category == "unknown") & (sp1 == "noise"))
  
  # This `NA` value is due to an `NA` in the `sp1` column which is the result of a tie in the random forest probabilities:  
  bat_data_most_weight_df[is.na(index_remove_noise),]
  
  # Set any `NA` value to `FALSE` (i.e. not a noise value to be filtered out):
  index_remove_noise[is.na(index_remove_noise)] <- FALSE 
  
  # For the noise values to be filtered out, get the corresponding probabilities:
  prob_noise <- bat_data_most_weight_df[index_remove_noise, "prob1"]
  
  # Record the state of the data frame before the filtering out of noise:
  bat_data_most_weight_df_before_noise_removal <- bat_data_most_weight_df
  
  # For the call sequences that are not categorized as noise, the following is done: shift the values from the `sp2` column into the `sp1`; shift the values from the `sp3` column into the `sp2` column; shift the values from the `prob2` column into the `prob1` column; shift the values from the `prob3` column into the `prob2` column; set the values in the `sp3` and `prob3` columns to `NA`:
  bat_data_most_weight_df[index_remove_noise, "sp1"] <- bat_data_most_weight_df[index_remove_noise, "sp2"]
  
  bat_data_most_weight_df[index_remove_noise,"sp2"] <- bat_data_most_weight_df[index_remove_noise, "sp3"]
  
  bat_data_most_weight_df[index_remove_noise, "prob1"] <- bat_data_most_weight_df[index_remove_noise, "prob2"]
  
  bat_data_most_weight_df[index_remove_noise,"prob2"] <- bat_data_most_weight_df[index_remove_noise, "prob3"]
  
  # Now that the `prob3` value has been shifted to the `prob2` column, replace the `prob3` column with `NA` values, for the cases where noise has been filtered out:
  bat_data_most_weight_df[index_remove_noise, "prob3"] <- NA
  bat_data_most_weight_df[index_remove_noise, "sp3"] <- NA
  
  # Now standardize `prob1` and `prob2` to account for the removal of the noise category. For example, if a "noise" value in the "sp1" column with a "prob1" value of 0.5 has been filtered out, the remaining probabilities will be multiplied by 2 so that the remaining probabilities sum to 1.  
  bat_data_most_weight_df[index_remove_noise, "prob1"] <- bat_data_most_weight_df[index_remove_noise, "prob1"] / (1 - prob_noise)
  
  bat_data_most_weight_df[index_remove_noise, "prob2"] <- bat_data_most_weight_df[index_remove_noise, "prob2"] / (1 - prob_noise)
  
  # Categorize as a single species if `prob1` is above the value for `threshold_bat`, and prob2 is not too high. For a categorization to take place the call sequence must have a length of at least 3, the `prob1` value must be above the `threshold` bat column, and the `prob2` value divided by the `prob1` value must be less than or equal to 0.80:
  index_bat_species <- with(bat_data_most_weight_df, (category == "unknown") & (n_calls >= 3) & (prob1 > threshold_bat) & (prob2 / prob1 <= 0.80))
  
  bat_data_most_weight_df$category[index_bat_species] <- bat_data_most_weight_df$sp1[index_bat_species]
 
  # If the category is "unknown" consider the possibility of aggregating sp1 and sp2 into an aggregated classification:
  
  LABO.MYLU_index <- (bat_data_most_weight_df$sp1 %in% c("LABO", "MYLU")) & (bat_data_most_weight_df$sp2 %in% c("LABO", "MYLU")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)
  
  EPFU.LANO_index <- (bat_data_most_weight_df$sp1 %in% c("EPFU", "LANO")) & (bat_data_most_weight_df$sp2 %in% c("EPFU", "LANO")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)
  
  MYEV.MYSE_index <- (bat_data_most_weight_df$sp1 %in% c("MYEV", "MYSE")) & (bat_data_most_weight_df$sp2 %in% c("MYEV", "MYSE")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)
  
  # If the category is "unknown" and the conditions for aggregrations of two species have not been satisfied, consider an aggregation of three species:  
  
  My_40k_index <- (bat_data_most_weight_df$sp1 %in% c("MYCI", "MYLU", "MYVO")) & (bat_data_most_weight_df$sp2 %in% c("MYCI", "MYLU", "MYVO")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)
  
  # Define a new data frame that includes the aggregated categories:  
  
  bat_data_most_weight_aggregated_df <- bat_data_most_weight_df
  
  bat_data_most_weight_aggregated_df$category[LABO.MYLU_index] <- "LABO.MYLU"
  bat_data_most_weight_aggregated_df$category[EPFU.LANO_index] <- "EPFU.LANO"
  bat_data_most_weight_aggregated_df$category[My_40k_index] <- "My_40k"
  bat_data_most_weight_aggregated_df$category[MYEV.MYSE_index] <- "MYEV.MYSE"
  
  # In the final phase of categorization, remove all calls with less than 3 calls
  bat_data_most_weight_df <- bat_data_most_weight_df %>% filter(n_calls >2)
  
  # Add a column for the ratio of the second-to-largest average probability to the maximum average probability, and make a data plot:  
  bat_data_most_weight_df$prob2_over_prob1 <- (bat_data_most_weight_df$prob2 / bat_data_most_weight_df$prob1)
  
  bat_data_most_weight_df_wide <- dcast(bat_data_most_weight_df, Site + Date + Date_POSIX ~ category, fun.aggregate = length, value.var = "Filename")
  
  bat_data_counts <- melt(bat_data_most_weight_df_wide, id.vars = c("Site","Date", "Date_POSIX"), variable.name = "category", value.name = "count")
  
  bat_data_counts_by_site <- summarize(group_by(bat_data_counts, Site, category), mean_count = mean(count), max = max(count), min = min(count), n_days = length(Date), SD = sd(count))
  
  bat_data_counts_by_site <- as.data.frame(bat_data_counts_by_site)
  
  bat_data_counts_by_site$SE <- bat_data_counts_by_site$SD / sqrt(bat_data_counts_by_site$n_days)
  
  # Get a count within each Site and Date for Aggregated Categories:
  
  bat_data_most_weight_aggregated_df_wide <- dcast(bat_data_most_weight_aggregated_df, Site + Date + Date_POSIX ~ category, fun.aggregate = length, value.var = "Filename")
  
  # Add columns of zeros for categories that are absent, so that absent categories can be shown:  
  
  if (!("EPFU" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$EPFU <- 0
  }
  
  if (!("LANO" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$LANO <- 0
  }
  
  if (!("EPFU.LANO" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$EPFU.LANO <- 0
  }
  
  if (!("LACI" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$LACI <- 0
  }
  
  if (!("LABO" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$LABO <- 0
  }
  
  if (!("LABO.MYLU" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$LABO.MYLU <- 0
  }
  
  if (!("MYLU" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYLU <- 0
  }
  
  if (!("My_40k" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$My_40k <- 0
  }
  
  if (!("MYEV" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYEV <- 0
  }
  
  if (!("MYEV" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYEV <- 0
  }
  
  if (!("MYEV.MYSE" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYEV.MYSE <- 0
  }
  
  if (!("MYSE" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYSE <- 0
  }
  
  if (!("MYVO" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYVO <- 0
  }
  
  if (!("MYCA" %in% names(bat_data_most_weight_aggregated_df_wide))){
    bat_data_most_weight_aggregated_df_wide$MYCA <- 0
  }
  
  names(bat_data_most_weight_aggregated_df_wide)
  
  bat_data_counts_aggregated <- melt(bat_data_most_weight_aggregated_df_wide, id.vars = c("Site","Date", "Date_POSIX"), variable.name = "category", value.name = "count")
  
  bat_data_counts_aggregated_output <- bat_data_counts_aggregated
  
  # Make a new version of the data frame that puts the column names into a standard form, that can be used as the basis for a final output:
  bat_data_counts_aggregated_output <- bat_data_counts_aggregated
  
  bat_data_counts_aggregated_output <- rename(bat_data_counts_aggregated_output, "Site_Combined" = Site, "Classification" = category, "Count" = count)
  
  bat_data_counts_aggregated_output$Site <- word(str_replace(bat_data_counts_aggregated_output$Site_Combined, "\\_", " "), 1)
  
  bat_data_counts_aggregated_output$Detector <- word(str_replace(bat_data_counts_aggregated_output$Site_Combined, "\\_", " "), 2)
  bat_data_counts_aggregated_output$Detector <- case_when(is.na(bat_data_counts_aggregated_output$Detector) ~ bat_data_counts_aggregated_output$Site,
                                                          TRUE ~ as.character(bat_data_counts_aggregated_output$Detector))
  
  bat_data_counts_aggregated_output <- bat_data_counts_aggregated_output[ ,c("Site", "Date", "Detector", "Classification", "Count")]
 
  # The two data frames with the ending "output" provide summaries of the classification probabilities, and the species counts. These were shown earlier in this document, but are being re-displayed here for the sake of clarity.  
  # Select columns to include in the output file:
  bat_data_output <- bat_data_output[ ,c("Site", "Date", "Filename", "n_calls", "prob1", "sp1", "prob2", "sp2", "prob3", "sp3")]
  
  # Note that the single case of `NA` values in `sp1` and `sp2` has resulted in a categorization of "unknown" as would be expected:  
  bat_data_output[is.na(bat_data_output$sp1), ]
  bat_data_most_weight_df[is.na(bat_data_most_weight_df$sp1),c("Site", "Date","sp1", "sp2", "category")]
  
  # Display the output file out bat species counts:
  bat_data_class_counts_agg <- bat_data_counts_aggregated_output[bat_data_counts_aggregated_output$Count>0,]

  return(list(bat_counts = bat_data_class_counts_agg, bat_summary = bat_data_output))
}
