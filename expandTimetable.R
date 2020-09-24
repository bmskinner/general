# Reformat the School timetable into something useful

library(tidyverse)
library(lubridate)

# Copy the School timetable from the website and paste into a text file. 
# It should be tab separated by default
in.file = "2020-21_Life_Sciences_AU_full.txt"

# Set the first day of term
Week1Day1 = lubridate::ymd("2020-10-05")

data = read.csv(file=in.file, sep="\t", header = T, stringsAsFactors = F)

# Split the columns of type 'Week 2-9' into separate rows for each week
trans = data %>% tidyr::separate_rows(Weeks, sep = ", " ) %>%
  dplyr::mutate(StartWeek = as.numeric(stringr::str_extract(Weeks, "^(\\d+)")),
                EndWeek = as.numeric(stringr::str_extract(Weeks, "(\\d+$)"))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Week = paste(seq(StartWeek, EndWeek), collapse=", ")) %>%
  tidyr::separate_rows(Week, sep = ", " ) %>%
  dplyr::select(-StartWeek, -EndWeek, -Weeks) %>%
  dplyr::mutate(Day = factor(Day, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")),
                Week = as.numeric(Week))

# Find the date of each event by counting from week 1
trans$Date = Week1Day1 + lubridate::dweeks(trans$Week-1) + lubridate::ddays(as.integer(trans$Day)-1)

# Order by date
trans = trans %>% dplyr::arrange(Date, Start)

# Export the new table
write.table(trans, file = paste0(in.file, ".converted.txt"), sep="\t", col.names = T, row.names = F)



