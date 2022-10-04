# Reformat the School timetable into something useful

library(tidyverse)
library(xlsx)

# Create an Excel file with column filtering
create.xlsx = function(data, file.name){
  oldOpt = options()
  options(xlsx.date.format="yyyy-mm-dd") # change date format
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb, sheetName = "Timetable")
  xlsx::addDataFrame(data, sh, row.names = F)
  xlsx::addAutoFilter(sh, "A1:L1")
  xlsx::createFreezePane(sh, 2, 1, 2, 1) # freeze top row
  xlsx::autoSizeColumn(sh, 6) # hardcoded col 6 as Date
  xlsx::saveWorkbook(wb, file=file.name)
  options(oldOpt)
}

in.file = commandArgs(trailingOnly=TRUE)[1] #"2022-07-25_draft_timetable.tsv"

data = read.csv(file=in.file, sep="\t", header = T, stringsAsFactors = F)

Week1Day1 = lubridate::ymd("2022-10-03") # Start of term

# Transform input to long format
trans = data %>% tidyr::separate_rows(Weeks, sep = ", " ) %>%
  dplyr::mutate(StartWeek = as.numeric(stringr::str_extract(Weeks, "^(\\d+)")),
                EndWeek = as.numeric(stringr::str_extract(Weeks, "(\\d+$)"))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Week = paste(seq(StartWeek, EndWeek), collapse=", ")) %>%
  tidyr::separate_rows(Week, sep = ", " ) %>%
  dplyr::select(-StartWeek, -EndWeek, -Weeks) %>%
  dplyr::mutate(Day = factor(Day, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")),
                Week = as.numeric(Week)) %>%
  dplyr::mutate(Date = Week1Day1 + lubridate::dweeks(Week - 1) + lubridate::ddays(as.integer(Day)-1))  %>%
  dplyr::arrange(Date, Start) %>%
  dplyr::select(Term, Module.s., Group.s., Type, Module.Name, Date, Week, Day, Start, Finish, Lecturer.s., Location.s.) %>%
  as.data.frame # needed for write.xlsx

# trans$Date = Week1Day1 + lubridate::dweeks(trans$Week - 1) + lubridate::ddays(as.integer(trans$Day)-1)
# trans = trans %>% dplyr::arrange(Date, Start) %>%
#   dplyr::select(Term, Module.s, Group.s, Type, Module.Name, Date, Week, Day, Start, Finish, Lecturer.s, Location.s)
# write.table(trans, file = paste0(in.file, ".converted.txt"), sep="\t", col.names = T, row.names = F)

# Get the modules I teach on
my.modules = trans %>% dplyr::filter(grepl("Skinner", Lecturer.s.)) %>%
  dplyr::select(Module = Module.s.) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Module)

mine = trans %>% dplyr::filter(Module.s. %in% my.modules$Module) %>%
  as.data.frame

# Export the data to Excel files
create.xlsx(trans, file.name = paste0(in.file, ".converted.xlsx"))
create.xlsx(mine, file.name = paste0(in.file, ".converted.mine.xlsx"))
