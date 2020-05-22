# This script demonstrates how to read
# an exported file from Nuclear Morphology Analysis,
# filter and select columns of interest, and create
# a basic plot of some data

library(tidyverse)

# Assign a file name in the working directory to 
# the variable 'data.file'
data.file = "plot_multiple_ggplots_demo.csv"

# Read the contents of the file into the variable 'data'
# header = T tells R that the first line contains column names
data = read.csv(data.file, header=T, stringsAsFactors = F)

# The 'select' function in the dplyr package allows specific
# columns to be selected
# The '%>%' operator is a pipe: the contents of 'data' is piped into 
# the select function
data.subset = data %>% dplyr::select(Dataset, Area_square_microns, Perimeter_microns)


# The 'filter' function in the dplyr package allows rows to be 
# excluded based on their values - e.g. only rows with an area
# greater than 22 square microns
# The '%>%' operator is a pipe: the contents of 'data' is piped into 
# the filter function
data.filtered = data %>% dplyr::filter(Area_square_microns>20)

# Piping allows functions to be chained. We can filter and select in 
# one statement:
data.filter.subset = data %>% 
  dplyr::select(Dataset, Area_square_microns, Perimeter_microns) %>%
  dplyr::filter(Area_square_microns>20)

# ggplot allows simple plotting of charts
# E.g. a violin chart as seen in Nuclear Morphology Analysis
# In ggplot syntax, the '+' at the end of each line is a pipe,
# like the '%>%' we used in dplyr
# In this case, we plot the entire dataset ('data'). Change this
# to one of the other filtered or subset datasets to see how the 
# plot changes
ggplot(data, aes(x=Dataset, y=Area_square_microns))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  ylab("Area (square microns)")

