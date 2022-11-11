#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#' ---
#' title: "Compute Bat Species Counts"
#' date: "Mar 29, 2017"
#' ---

#' R Script by: Allan Roberts  
#' Modified by QReserve Inc. for Alberta eBat Pipeline Integration
#'  


#' ###Introduction  
#'  

#' The `R` script in this document computes bat species counts, based on the `M1_Edson_2016.csv` data file. The goal is to have an R script that can be conveniently adapted to other datasets. The method used is to average the probabilities for each category (species or noise) within each call sequence; each call sequence is then categorized according to the overall greatest average probability.  
#'  

#' For a call sequence to be positively categorized, a call sequence must have a minimum of three calls (i.e. rows); for call sequences with fewer than 3 rows, the species is categorized as "unknown". Minimum probability thresholds are required for the categorization of both bat calls, and noise. Here, the average probability allocated to a bat species, by the random forest algorithm, must be greater than 0.50 for a call sequence to be categorized as that bat species; for a call sequence to be categorized as "noise" the average probability for the noise category must be greater than 0.80. These thresholds are set by the values for `threshold_bat` and threshold_noise`.  
#'  

#' The following aggregated categories are considered: "EPFU-LANO" (combined EPFU and LANO calls), "LABO-MYLU" (LABO and MYLU calls), "MYEV-MYSE" (MYEV and MYSE calls), and "My_40k" (MYCI, MYLU and MYVO calls).  
#'  

#' This script requires that the `Date` column be in a standard format, so that it can be converted into a POSIX time format; however, if each date is represented by a unique character string, it should be fairly easy to adapt this script so that it can compute mean site counts, by simply removing any use of POSIX time.  
#'  


#' Load some R packages:  
#'  

library(reshape2)	# for formatting data frames
library(dplyr)		# for applying functions to subsets of data frames
library(ggplot2)	# for data visualization
library(stringr)	# for working with character strings
library(tidyr)		# for data formatting functions
library(knitr)		# for the "kable" function for formatting tables

#' $\pagebreak$  
#'  





#' ###Define Some Functions  
#'  


#' Define functions for computing the category with the highest output probability, and the second and third highest probability. In the case of ties, an NA value is returned:  
#'  

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



#' ###Loading a dataset
#'  

#setwd("C://Users/JBurgar/Documents/R/Analysis/BatID_Automation/data")


#' Load random forest outputs:  

bat_data <- read.csv(paste(args[1], "M1_analook.csv", sep="/"))
names(bat_data)
head(bat_data)

#' Add an index for the row number:
bat_data$row <- 1:nrow(bat_data)

#' Check that all of the random forest categorization probabilities sum to 1.

# Not run:
# n <- nrow(bat_data)
# S <- numeric(n)
# for (i in 1:n){
#	S[i] <- sum(bat_data[i, 1:11])
# }

# print(min(S), digits = 20)

#' $\pagebreak$  
#'  





#' ###The Date Column  
#'  


#' Convert the `Date` column to character strings:
bat_data$Date <- as.character(bat_data$Date)

#' Check that length of the character strings is the same for all values in the `Date` column:  

table(nchar(bat_data$Date))

#' Add a column for the date in POSIX format:
bat_data$Date_POSIX <- as.POSIXct(strptime(bat_data$Date, format = "%Y%m%d"))

#' $\pagebreak$  
#'  






#' Check the number of rows for each site:
with(bat_data, kable(data.frame(table(Site, useNA = "always"))))

#' Check the number of rows for each date:
with(bat_data, kable(data.frame(table(Date_POSIX, useNA = "always"))))

#' $\pagebreak$  
#'  







#' ###The Number of Rows for Each Filename  
#'  

#' Make a dotchart of the number of rows for each file:  

file_lengths <- summarize(group_by(bat_data, Site, Date_POSIX, Filename), file_length = length(Filename))

range(file_lengths$file_length)

#' Make a histogram of the file lengths:

hist(file_lengths$file_length, breaks = seq(-0.5, 275.5, 1), col = "darkgrey", border = "lightgrey", main = "Number of rows for call sequences", font.main = 1, las = 1, xlab = "Number of rows in sequence", ylab = "Number of sequences")

file_lengths$row <- 1:nrow(file_lengths)

p <- ggplot(file_lengths)
p <- p + theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white"))
p <- p + geom_point(aes(file_length, row), size = 1)
p <- p + facet_wrap(~Site, ncol = 1, scales = "free_y")
p <- p + xlab("Length of Call Sequence in Rows") + ylab("")
p <- p + theme(axis.text.y = element_blank())
p <- p + theme(axis.ticks.y = element_blank())
print(p)

#' $\pagebreak$  
#'  





#' ###Add a Column for the Number of Rows per Call Sequence 
#'  

#' Here a column is added to the data frame of bat data, indicating the number of calls (i.e. rows) for each call sequence; the information in this column is later used to set the `species` column to "unknown" for call sequences with fewer than 3 rows.  
#'  

# Check the number of distinct file names: 
n_distinct(bat_data$Filename)

# Add a new column that pastes together the site name and the file name; this is intended to help provide unique identifiers: 

bat_data$Site_Filename <- paste(bat_data$Site, bat_data$Filename)


#' Add a column indicating the number of rows for each call sequence. Rows are counted for each file name, within each site:
bat_data <- mutate(group_by(bat_data, Site_Filename), n_calls = length(Site_Filename))



#' $\pagebreak$  
#'  







#' ###Counting Species by Highest Average Probability  
#'    

#' Convert the data frames to a long format:  
 
bat_data_long <- melt(bat_data, variable.name = "species", value.name = "probability", id.vars = c("Site", "Date", "Filename", "Site_Filename", "row", "Date_POSIX", "n_calls"))


#' Abbreviate the `species` column by removing the substring `rf4.`:

bat_data_long$species <- as.character(bat_data_long$species)
bat_data_long$species <- str_replace(bat_data_long$species, "\\.", " ")
bat_data_long$species <- word(bat_data_long$species, 2)


#' Get the average probability for each species (averaged across calls) within each call sequence (i.e. `Filename` within `Site`):
bat_data_mean_probabilities <- summarize(group_by(bat_data_long, Site, Date, Date_POSIX, Filename, species, n_calls), probability = mean(probability))




#' Rank the mean probabilities for each category in each call sequence, from highest to lowest, with the highest probability having a rank of 1:

bat_data_mean_probabilities <- mutate(group_by(bat_data_mean_probabilities, Site, Date, Date_POSIX, Filename),  rank = rank(1 - probability))


#' Get the species with the most weight within each category. When species are not combined into aggregated categories, the value in the "category" column will generally be the same as the value in the "sp1" column, or "unknown" if the "prob1" value, corresponding to the category given by "sp1" is not above the required threshold given by `threshold_bat`. If the number of rows in a call sequence is less than 3, the category column is set to "unknown".  
#'  


bat_data_most_weight_df <- summarize(group_by(bat_data_mean_probabilities, Site, Date, Date_POSIX, Filename, n_calls), prob1 = max(probability), sp1 = get_most_probable(species, rank), prob2 = rev(sort(probability))[2], sp2 = get_second_probable(species, rank), prob3 = rev(sort(probability))[3], sp3 = get_third_probable(species, rank))

#' View the first few rows of the data frame:
head(bat_data_most_weight_df[ ,c("Site", "Date", "n_calls", "prob1", "sp1", "prob2", "sp2", "prob3", "sp3")])

#' Get some summaries of the data frame:
table(bat_data_most_weight_df$sp1, useNA = "always")

#' ###Ties in Classification Values can Result in NA Values  
#'  

#' This section shows that a tie in the values for `prob1` and `prob2` can result in an `NA` value in the `sp1` and `sp2` columns. An `NA` value would also result if there were a tie between the `prob2` and `prob3` columns. The number of such cases is expected to be very small, and to result in a classification of "unknown".  
#'   


#' There is one row where a tied value between `prob1` and `prob2` has resulted in `NA` values for `sp1` and `sp2`:
bat_data_most_weight_df[is.na(bat_data_most_weight_df$sp1), c("Site", "Date", "n_calls", "prob1", "sp1", "prob2", "sp2", "prob3", "sp3") ]  
#'  

#' The `NA` values in the `sp2` column are associated with `prob1` and `prob2` values that are tied at 0:
bat_data_most_weight_df[is.na(bat_data_most_weight_df$sp2), c("Site", "Date", "n_calls", "prob1", "sp1", "prob2", "sp2", "prob3", "sp3") ]


#' $\pagebreak$  
#'  







#' ###Create an Output File with the Classification Probabilities   
#'  

bat_data_output <- bat_data_most_weight_df

bat_data_output <- rename(bat_data_output, "Site_Combined" = Site)

bat_data_output$Site <- word(str_replace(bat_data_output$Site_Combined, "/_", " "), 1)

bat_data_output$Detector <- word(str_replace(bat_data_output$Site_Combined, "/_", " "), 2)

#' The resulting data frame, `bat_data_output` is displayed in the section at the end of this document, after a selection of the columns to be included has been made.  
#'  



#' $\pagebreak$  
#'  







#' ###Setting Thresholds for the Classification of Bats and Noise  
#'  

#' Now thresholds are set for the positive identification of noise, and of bat species. These thresholds are used to determine if the probability for the category is high enough for a classification to be made.  
#'  
  

threshold_noise <- 0.8
threshold_bat   <- 0.5


#' Start with all call sequences categorized as "unknown":  
bat_data_most_weight_df$category <- "unknown"

#' In the first phase of categorization, all call sequences that have a `sp1` value of `noise`, and a `prob1` value over the threshold for the positive identification of noise (here, 0.80) are categorized as noise:  

# Index cases to be categorized as noise:
index_noise <- with(bat_data_most_weight_df, (sp1 == "noise") & (prob1 > threshold_noise))

# Set the indexed categories to "noise":
bat_data_most_weight_df$category[index_noise] <- "noise"

#' After the first phase of categorization, all call sequences are categorized as `noise` or as `unknown`:
table(bat_data_most_weight_df$category, useNA = "always")


#' $\pagebreak$  
#'  


#' If sp1 == "noise" and `prob1` is below the threshold value for noise, filter out the noise value, and standardize the remaining probabilities:

# Index the noise values to be filtered out:
index_remove_noise <- with(bat_data_most_weight_df, (category == "unknown") & (sp1 == "noise"))

# Note that there is one `NA` value:
sum(is.na(index_remove_noise))

# This `NA` value is due to an `NA` in the `sp1` column which is the result of a tie in the random forest probabilities:  
bat_data_most_weight_df[is.na(index_remove_noise),]

# Set any `NA` value to `FALSE` (i.e. not a noise value to be filtered out):
index_remove_noise[is.na(index_remove_noise)] <- FALSE 

# For the noise values to be filtered out, get the corresponding probabilities:
prob_noise <- bat_data_most_weight_df[index_remove_noise, "prob1"]

# Record the state of the data frame before the filtering out of noise:
bat_data_most_weight_df_before_noise_removal <- bat_data_most_weight_df

#' For the call sequences that are not categorized as noise, the following is done: shift the values from the `sp2` column into the `sp1`; shift the values from the `sp3` column into the `sp2` column; shift the values from the `prob2` column into the `prob1` column; shift the values from the `prob3` column into the `prob2` column; set the values in the `sp3` and `prob3` columns to `NA`:

bat_data_most_weight_df[index_remove_noise, "sp1"] <- bat_data_most_weight_df[index_remove_noise, "sp2"]

bat_data_most_weight_df[index_remove_noise,"sp2"] <- bat_data_most_weight_df[index_remove_noise, "sp3"]

bat_data_most_weight_df[index_remove_noise, "prob1"] <- bat_data_most_weight_df[index_remove_noise, "prob2"]

bat_data_most_weight_df[index_remove_noise,"prob2"] <- bat_data_most_weight_df[index_remove_noise, "prob3"]

# Now that the `prob3` value has been shifted to the `prob2` column, replace the `prob3` column with `NA` values, for the cases where noise has been filtered out:
bat_data_most_weight_df[index_remove_noise, "prob3"] <- NA
bat_data_most_weight_df[index_remove_noise, "sp3"] <- NA


#' Now standardize `prob1` and `prob2` to account for the removal of the noise category. For example, if a "noise" value in the "sp1" column with a "prob1" value of 0.5 has been filtered out, the remaining probabilities will be multiplied by 2 so that the remaining probabilities sum to 1.  
#'  


bat_data_most_weight_df[index_remove_noise, "prob1"] <- bat_data_most_weight_df[index_remove_noise, "prob1"] / (1 - prob_noise)

bat_data_most_weight_df[index_remove_noise, "prob2"] <- bat_data_most_weight_df[index_remove_noise, "prob2"] / (1 - prob_noise)

# Compare "sp1" columns before and after filtering out noise:
table(bat_data_most_weight_df_before_noise_removal$sp1)
table(bat_data_most_weight_df$sp1)

#' Next following tables check the effect of removing the `noise` values that were below the threshold for noise detection.

#' This plot shows that `noise` values removed from the `sp1` column were usually replaced with `LACI` (267 cases):
table(bat_data_most_weight_df_before_noise_removal$sp1, bat_data_most_weight_df$sp1, useNA = "always")

#' In the following table, the values in the row labelled `LACI` show which species classifications were moved into the `sp2` column when the `LACI` values were moved into the `sp1` column:  

table(bat_data_most_weight_df_before_noise_removal$sp2, bat_data_most_weight_df$sp2, useNA = "always")

#' Removing the sub-threshold noise values from `sp1` results in new `NA` values in the `sp3` column, but otherwise leaves the `sp3` column unchanged:
table(bat_data_most_weight_df_before_noise_removal$sp3, bat_data_most_weight_df$sp3, useNA = "always")



#' For this dataset, the most obvious effect of filtering out noise observations below the threshold for noise has been to classify more call sequences as `LACI`.  
#'  

#' Categorize as a single speces if `prob1` is above the value for `threshold_bat`, and prob2 is not too high. For a categorization to take place the call sequence must have a length of at least 3, the `prob1` value must be above the `threshold` bat column, and the `prob2` value divided by the `prob1` value must be less than or equal to 0.80:

index_bat_species <- with(bat_data_most_weight_df, (category == "unknown") & (n_calls >= 3) & (prob1 > threshold_bat) & (prob2 / prob1 <= 0.80))

bat_data_most_weight_df$category[index_bat_species] <- bat_data_most_weight_df$sp1[index_bat_species]


#' $\pagebreak$  
#'  

#' Check the number of files that have been categorized as "unknown":

table(bat_data_most_weight_df$category == "unknown", useNA = "always")


#' Check the number of call sequence files with fewer than 3 rows:

table(bat_data_most_weight_df$n_calls < 3, useNA = "always")





#' Check the resulting data frame:
print(as.data.frame(bat_data_most_weight_df)[1:10 ,c("Site", "Date_POSIX","category", "sp1", "prob1", "sp2", "prob2", "sp3", "prob3")], digits = 2)

#' $\pagebreak$  
#'  





#' ###Aggregate Species Categories

#' This section uses logical conditions to create a new data frame called `most_weight_aggregated_df`, which aggregates certain combinations in the `sp1` and `sp2` columns of `most_weight_aggregated`.  
#'  

#' If the category is "unknown" consider the possibility of aggregating sp1 and sp2 into an aggregated classification:

LABO.MYLU_index <- (bat_data_most_weight_df$sp1 %in% c("LABO", "MYLU")) & (bat_data_most_weight_df$sp2 %in% c("LABO", "MYLU")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)

EPFU.LANO_index <- (bat_data_most_weight_df$sp1 %in% c("EPFU", "LANO")) & (bat_data_most_weight_df$sp2 %in% c("EPFU", "LANO")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)

MYEV.MYSE_index <- (bat_data_most_weight_df$sp1 %in% c("MYEV", "MYSE")) & (bat_data_most_weight_df$sp2 %in% c("MYEV", "MYSE")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)

#' If the category is "unknown" and the conditions for aggregrations of two species have not been satisfied, consider an aggregation of three species:  

My_40k_index <- (bat_data_most_weight_df$sp1 %in% c("MYCI", "MYLU", "MYVO")) & (bat_data_most_weight_df$sp2 %in% c("MYCI", "MYLU", "MYVO")) & (bat_data_most_weight_df$category == "unknown") & ((bat_data_most_weight_df$prob1 + bat_data_most_weight_df$prob2) > threshold_bat)

#' Check the number of call sequences belonging to each aggregated category:
sum(LABO.MYLU_index)
sum(EPFU.LANO_index)
sum(MYEV.MYSE_index)
sum(My_40k_index)

#' $\pagebreak$  
#'  







#' Define a new data frame that includes the aggregated categories:  

bat_data_most_weight_aggregated_df <- bat_data_most_weight_df

bat_data_most_weight_aggregated_df$category[LABO.MYLU_index] <- "LABO.MYLU"
bat_data_most_weight_aggregated_df$category[EPFU.LANO_index] <- "EPFU.LANO"
bat_data_most_weight_aggregated_df$category[My_40k_index] <- "My_40k"
bat_data_most_weight_aggregated_df$category[MYEV.MYSE_index] <- "MYEV.MYSE"

#' Compare the aggregated with the non-aggregated categorizations:  

kable(table(bat_data_most_weight_aggregated_df$category, bat_data_most_weight_df$category))

#' $\pagebreak$  
#'  






#' Do a visual check of the `My_40k` categorizations:

print(subset(as.data.frame(bat_data_most_weight_aggregated_df), category == "My_40k")[,c("n_calls", "sp1", "prob1", "sp2", "prob2", "category")], digits = 2)

#' Do a visual check of the `unknown` categorizations:

#print(subset(as.data.frame(bat_data_most_weight_aggregated_df), category == "unknown")[,c("n_calls", "sp1", "prob1", "sp2", "prob2", "category")], digits = 2)


#' $\pagebreak$  
#'  


#' In the final phase of categorization, make sure that the value for `category` to "unknown" whenever the number of rows of data for a call sequence is less than 3:

bat_data_most_weight_df$category[bat_data_most_weight_df$n_calls < 3] <- "unknown"




#' ###Uncertainty of Categorizations 
#'  

#' Make a boxplot of the confidence (average probability) that has been used for the categorization of each call sequence. This table shows the probability value that was used for each classification. For the bat species, these probabilities are greater than 0.50 because that was the value that was used for `threshold_bat`; for "noise" the probabilities are greater than 0.80 because that was the value that was used for `threshold_noise`. Cases that did not meet these thresholds were categorized as "unknown". Some "unknown" classifications can also be due to call sequences with fewer than 2 rows, or to cases where the "prob1" and "prob2" values are tied, or close to tied (e.g. if prob2/prob1 > 0.8).

par(mar = c(7, 7, 3, 3))

with(bat_data_most_weight_df, plot(factor(category), prob1, las = 2, cex = 0, cex.axis = 0.7))

with(bat_data_most_weight_df, points(jitter(as.numeric(factor(category))), prob1, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.5)))

mtext(side = 2, line = 3, "Maximum Average Probability\nwithin Call Sequence", cex = 1)

mtext(side = 1, line = 4, "Categorization of Call Sequence", cex = 1)

abline(h = 0.5, lty = "dashed")	# add a reference line

#' $\pagebreak$  
#'  





#' Add a column for the ratio of the second-to-largest average probability to the maximum average probability, and make a data plot:  
#'  



bat_data_most_weight_df$prob2_over_prob1 <- (bat_data_most_weight_df$prob2 / bat_data_most_weight_df$prob1)


par(mar = c(9, 7, 3, 3))
with(bat_data_most_weight_df, plot(factor(category), prob2_over_prob1, las = 2, cex = 0, cex.axis = 0.6))
with(bat_data_most_weight_df, points(jitter(as.numeric(factor(category))), prob2_over_prob1, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.5)))
mtext(side = 2, line = 3, "prob2_over_prob1", cex = 1)
mtext(side = 1, line = 5, "Categorization of Call Sequence", cex = 1)

#' The previous plot shows the value of the second-to-highest probability, divided by the highest probability, for each categorization that has been made.
#' $\pagebreak$  
#'  




#' ###Computing Counts by Day and Site  
#'  

#' Get species counts within each Site and Date. For each day within each site, the count for a species is equal to the number of call sequence files that have been categorized as representing that species:


bat_data_most_weight_df_wide <- dcast(bat_data_most_weight_df, Site +             Date + Date_POSIX ~ category, fun.aggregate = length, value.var = "Filename")


#' Display the data frame:  
head(bat_data_most_weight_df_wide[ ,-2])

bat_data_counts <- melt(bat_data_most_weight_df_wide, id.vars = c("Site",             "Date", "Date_POSIX"), variable.name = "category", value.name = "count")

#' Check that the total of counts for all categories matches the number of filenames:

sum(bat_data_counts$count)
nrow(bat_data_most_weight_df)
n_distinct(bat_data$Filename)

p <- ggplot(bat_data_counts)
p <- p + theme_bw() + theme(strip.background = element_blank())
p <- p + geom_point(aes(category, count))
p <- p + facet_wrap(~Site, ncol = 1, scales = "free_y")
p <- p + theme(axis.text.x = element_text(angl = 90))
print(p)

#' $\pagebreak$  
#'  






bat_data_counts_by_site <- summarize(group_by(bat_data_counts, Site, category), mean_count = mean(count), max = max(count), min = min(count), n_days = length(Date), SD = sd(count))

bat_data_counts_by_site <- as.data.frame(bat_data_counts_by_site)

bat_data_counts_by_site$SE <- bat_data_counts_by_site$SD / sqrt(bat_data_counts_by_site$n_days)

#' As a quality check, confirm that the total count for all categories has remained unchanged:
sum(bat_data_counts_by_site$mean_count * bat_data_counts_by_site$n_days)

#' Check that the counts for each categories match between the `bat_data_counts_by_site` data frame and the `bat_data_most_weight_df` data frame:

as.data.frame(summarize(group_by(bat_data_counts_by_site, category), total = sum(n_days * mean_count)))

as.data.frame(table(bat_data_most_weight_df$category))

#' $\pagebreak$  
#'  





#' Get a count within each Site and Date for Aggregated Categories:

bat_data_most_weight_aggregated_df_wide <- dcast(bat_data_most_weight_aggregated_df, Site + Date + Date_POSIX ~ category, fun.aggregate = length, value.var = "Filename")

head(bat_data_most_weight_aggregated_df_wide)[ ,1:min(7, ncol(bat_data_most_weight_aggregated_df_wide))] # possible replacement

#tail(bat_data_most_weight_aggregated_df_wide)[ ,1:7]


#' Add columns of zeros for categories that are absent, so that absent categories can be shown:  
  
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

bat_data_counts_aggregated <- melt(bat_data_most_weight_aggregated_df_wide, id.vars = c("Site",             "Date", "Date_POSIX"), variable.name = "category", value.name = "count")

p <- ggplot(bat_data_counts_aggregated)
p <- p + theme_bw() + theme(strip.background = element_blank())
p <- p + geom_point(aes(category, count))
p <- p + facet_wrap(~Site, ncol = 1, scales = "free_y")
p <- p + theme(axis.text.x = element_text(angl = 90))
print(p)

#' $\pagebreak$  
#'  

bat_data_counts_aggregated_output <- bat_data_counts_aggregated

#' Make a new version of the data frame that puts the column names into a standard form, that can be used as the basis for a final output:

bat_data_counts_aggregated_output <- bat_data_counts_aggregated

bat_data_counts_aggregated_output <- rename(bat_data_counts_aggregated_output, "Site_Combined" = Site, "Classification" = category, "Count" = count)

bat_data_counts_aggregated_output$Site <- word(str_replace(bat_data_counts_aggregated_output$Site_Combined, "\\_", " "), 1)

bat_data_counts_aggregated_output$Detector <- word(str_replace(bat_data_counts_aggregated_output$Site_Combined, "\\_", " "), 2)

bat_data_counts_aggregated_output <- bat_data_counts_aggregated_output[ ,c("Site", "Date", "Detector", "Classification", "Count")]

head(bat_data_counts_aggregated_output)






bat_data_counts_aggregated_by_site <- summarize(group_by(bat_data_counts_aggregated, Site, category), mean_count = mean(count), max = max(count), min = min(count), n_days = length(Date), SD = sd(count))

bat_data_counts_aggregated_by_site <- as.data.frame(bat_data_counts_aggregated_by_site)

bat_data_counts_aggregated_by_site$SE <- bat_data_counts_aggregated_by_site$SD / sqrt(bat_data_counts_aggregated_by_site$n_days)


#' As a quality check, confirm that the total counts for each category have remain unchanged:
summarize(group_by(bat_data_counts_aggregated_by_site, category), total = sum(n_days * mean_count))

as.data.frame(table(bat_data_most_weight_aggregated_df$category))

#' $\pagebreak$  
#'  









#' Display the data as a table:

kable(bat_data_counts_by_site, digits = 2)

#' $\pagebreak$  
#'  


kable(bat_data_counts_aggregated_by_site, digits = 2)

#' $\pagebreak$  
#'  







#' ###Plot Counts by Site  
#'  

bat_data_counts_by_site$category <- factor(bat_data_counts_by_site$category, levels = c("EPFU", "LANO", "EPFU.LANO", "LACI", "LABO", "LABO.MYLU", "MYLU", "My_40k", "MYEV", "MYEV.MYSE", "MYSE", "MYVO", "MYCA", "unknown", "noise"))

bat_data_counts_by_site$colour <- "black"
bat_data_counts_by_site$colour[bat_data_counts_by_site$max == 0] <- "transparent"

bat_data_counts_by_site$colour <- "black"
bat_data_counts_by_site$colour[bat_data_counts_by_site$max == 0] <- "transparent"

bat_data_counts_by_site$shape <- 19
bat_data_counts_by_site$shape[bat_data_counts_by_site$max == 0] <- 1

#' Plot the mean counts and standard errors by site:  
#p <- ggplot(subset(bat_data_counts_by_site, !(category %in% c("noise", "unknown"))))
#p <- p + theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white"))
#p <- p + theme(panel.grid = element_blank())
#p <- p + geom_point(aes(category, mean_count, shape = shape))
#p <- p + facet_wrap(~Site, ncol = 1, scales = "free_y")
#p <- p + geom_errorbar(aes(category, ymin = mean_count - SE, ymax = mean_count + SE, colour = colour), width = 0.25)
#p <- p + scale_colour_identity()
#p <- p + scale_shape_identity()
#p <- p + theme(axis.title.x = element_blank())
#p <- p + theme(axis.text.x = element_text(angle = 90))
#p <- p + ylab("Mean ±1 SE Nightly Bat Activity")
#print(p)

#' Summarize information that could be useful for a caption:
#summarize(group_by(bat_data_most_weight_df, Site), n_days = n_distinct(Date), n_files = n_distinct(paste(Site, Filename)), n_bat = sum(!(category %in% c("noise", "unknown"))), n_unknown = sum(category == "unknown", na.rm = TRUE), n_noise = sum(category == "noise", na.rm = TRUE), n_short = sum(n_calls < 3))

#' An overall summary that is not broken down by `Unit`:
#summarize(group_by(bat_data_most_weight_df), n_days = n_distinct(Date), n_files = n_distinct(paste(Site, Filename)), n_bat = sum(!(category %in% c("noise", "unknown"))), n_unknown = sum(category == "unknown", na.rm = TRUE), n_noise = sum(category == "noise", na.rm = TRUE), n_short = sum(n_calls < 3))


#' $\pagebreak$  
#'  



#' Plot the mean counts and standard errors for the aggregated categories:  


bat_data_counts_aggregated_by_site$category <- factor(bat_data_counts_aggregated_by_site$category, levels = c("EPFU", "LANO", "EPFU.LANO", "LACI", "LABO", "LABO.MYLU", "MYLU", "My_40k", "MYEV","MYEV.MYSE", "MYSE", "MYVO", "MYCA", "unknown", "noise"))

bat_data_counts_aggregated_by_site$category <- as.character(bat_data_counts_aggregated_by_site$category)

bat_data_counts_aggregated_by_site$category <- str_replace(bat_data_counts_aggregated_by_site$category,"\\.","-")

bat_data_counts_aggregated_by_site$category <- str_replace(bat_data_counts_aggregated_by_site$category,"My_40k","Myotis 40k")

bat_data_counts_aggregated_by_site$category <- factor(bat_data_counts_aggregated_by_site$category, levels = c("EPFU", "LANO", "EPFU-LANO", "LACI", "LABO", "LABO-MYLU", "MYLU", "Myotis 40k", "MYEV", "MYEV-MYSE", "MYSE", "MYVO", "MYCA", "unknown", "noise"))

bat_data_counts_aggregated_by_site$colour <- "black"
bat_data_counts_aggregated_by_site$colour[bat_data_counts_aggregated_by_site$max == 0] <- "transparent"

bat_data_counts_aggregated_by_site$shape <- 19
bat_data_counts_aggregated_by_site$shape[bat_data_counts_aggregated_by_site$max == 0] <- 1


#' Define a version of the data frame to use for the data plot, which doesn't display the values for "noise" and "unknown", and which removes the underscore symbol from the site name: 

bat_data_counts_aggregated_by_site_display <- bat_data_counts_aggregated_by_site

bat_data_counts_aggregated_by_site_display$Site <- str_replace(bat_data_counts_aggregated_by_site_display$Site, "_", " ")

bat_data_counts_aggregated_by_site_display <- subset(bat_data_counts_aggregated_by_site_display, !(category %in% c("noise", "unknown", NA)))
str(bat_data_counts_aggregated_by_site_display)


#print(p %+% bat_data_counts_aggregated_by_site_display)
#ggsave(paste(args[1], "classification.png", sep="/"), width=8, height = 8)
#dev.off()

#' Summarize information that could be useful for a caption:
summarize(group_by(bat_data_most_weight_aggregated_df, Site), n_days = n_distinct(Date), n_files = n_distinct(paste(Site, Filename)), n_bat = sum(!(category %in% c("noise", "unknown"))), n_unknown = sum(category == "unknown", na.rm = TRUE), n_noise = sum(category == "noise", na.rm = TRUE), n_short = sum(n_calls < 3), n_NA = sum(is.na(category)))

#' An overall summary that is not broken down by `Unit`:
summarize(group_by(bat_data_most_weight_aggregated_df), n_days = n_distinct(Date), n_files = n_distinct(paste(Site, Filename)), n_bat = sum(!(category %in% c("noise", "unknown", na.rm = TRUE))), n_unknown = sum(category == "unknown", na.rm = TRUE), n_noise = sum(category == "noise", na.rm = TRUE), n_short = sum(n_calls < 3), n_NA = sum(is.na(category)))



#' $\pagebreak$  
#'  

#' ###View the Output files  
#'  

#' The two data frames with the ending "output" provide summaries of the classification probabilities, and the species counts. These were shown earlier in this document, but are being re-displayed here for the sake of clarity.  

# Select columns to include in the output file:
bat_data_output <- bat_data_output[ ,c("Site", "Date", "Filename", "n_calls", "prob1", "sp1", "prob2", "sp2", "prob3", "sp3")]

# Display the output file of classification probabilities:
glimpse(bat_data_output, width = 55)

write.csv(bat_data_output,paste(args[1], "anabat_class_prob.csv", sep="/"), row.names=F)

# Note that the single case of `NA` values in `sp1` and `sp2` has resulted in a categorization of "unknown" as would be expected:  

bat_data_output[is.na(bat_data_output$sp1), ]
bat_data_most_weight_df[is.na(bat_data_most_weight_df$sp1),c("Site", "Date","sp1", "sp2", "category")]

# Display the output file out bat species counts:
head(bat_data_counts_aggregated_output)
bat_data_class_counts_agg <- bat_data_counts_aggregated_output[bat_data_counts_aggregated_output$Count>0,]

write.csv(bat_data_class_counts_agg, paste(args[1], "anabat_class_counts_agg.csv", sep="/"), row.names=F)


# End of R script.
