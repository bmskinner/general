# Read Moodle quiz responses and save into separate docx files for each student.

library(tidyverse)
library(officer)

read.data = function(file) read.csv(file, header = T, sep = ",", stringsAsFactors = F)

# This should be the file of quiz responses downloaded from Moodle
file = "BS102-4-__-CO-__-__-Practical 3 Assessment-responses.csv"

# Read in the responses
data = read.data(file)

# Find the number of questions in the assessment
n.questions = data %>% dplyr::select(dplyr::starts_with("Response")) %>% ncol

# Ensure the output directory exists
out.dir = "report"
if(!dir.exists(out.dir))
  dir.create(out.dir)

# Write a single student response to a word file
# row - a single row from the input data
write.doc = function(row){
  
  doc = officer::read_docx() %>%
    officer::body_add_par(paste(row['First.name'], row['ï..Surname']))
  
  write.question = function(i){
    doc %>%
      officer::body_add_par(paste("Question", i, ":")) %>%
      officer::body_add_par(row[ paste0("Response.", i)]) %>%
      officer::body_add_par("")
  }
  
  lapply(1:n.questions, write.question)
  
  path = paste0(out.dir, "/", row['ID.number'],".docx")
  print(doc, target=path)
}

# Write each row in turn
apply(data, 1, write.doc)





