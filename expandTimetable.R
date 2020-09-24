# Reformat the School timetable into something useful

library(tidyverse)

in.file = "2020-21_Life_Sciences_AU_full.txt"

data = read.csv(file=in.file, sep="\t", header = T, stringsAsFactors = F)

Week1Day1 = lubridate::ymd("2020-10-05") # Start of term

trans = data %>% tidyr::separate_rows(Weeks, sep = ", " ) %>%
  dplyr::mutate(StartWeek = as.numeric(stringr::str_extract(Weeks, "^(\\d+)")),
                EndWeek = as.numeric(stringr::str_extract(Weeks, "(\\d+$)"))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Week = paste(seq(StartWeek, EndWeek), collapse=", ")) %>%
  tidyr::separate_rows(Week, sep = ", " ) %>%
  dplyr::select(-StartWeek, -EndWeek, -Weeks) %>%
  dplyr::mutate(Day = factor(Day, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")),
                Week = as.numeric(Week))

trans$Date = Week1Day1 + lubridate::dweeks(trans$Week-1) + lubridate::ddays(as.integer(trans$Day)-1)

trans = trans %>% dplyr::arrange(Date, Start)

write.table(trans, file = paste0(in.file, ".converted.txt"), sep="\t", col.names = T, row.names = F)



