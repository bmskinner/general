# Fetch the Entrez IDs for a list of entrez accessions using rentrez
library(dplyr)
library(rentrez)

gene.ids = c(76074, 99633, 99167)

# Split into parts if too many genes for one web call
genes.split = split(gene.ids, ceiling(seq_along(gene.ids)/100), drop=TRUE)

GetSummary = function(id.list){  
  
  summ    = entrez_summary(db="gene", id=id.list)
  
  # Get the basic values into a data frame
  values  = extract_from_esummary(summ, c("uid", "name", "chromosome"), simplify = T)
  tvalues = as.data.frame(t(values))
  tvalues = as.data.frame(lapply( tvalues, unlist))
  
  # Turn the genomic info into a data frame
  # This holds the current chromosome and position if available
  gen.list   = lapply(summ,  "[[", "genomicinfo")
  gen.df     = do.call(rbind.data.frame, gen.list)
  gen.df$uid = row.names(gen.df)
  rownames(gen.df) = c()
  
  # Merge in the location info
  tvalues = merge(x=tvalues, y=gen.df, by = c("uid") , all=T)
  
  # Genomic coordiates are in relation to the current assembly.
  # The location history can fetch previous mappings.
  # Included here to demonstrate how to extract this data
  
  # Turn the location history into a data frame
  # loc.list   = lapply(summ,  "[[", "locationhist")
  # loc.df     = do.call(rbind.data.frame, loc.list)
  # loc.df$uid = row.names(loc.df)
  # loc.df$uid = gsub("\\.\\d+", "", loc.df$uid)
  # rownames(loc.df) = c()
  
  # Filter for desired build, or comment out to fetch all builds
  # loc.df = loc.df %>% filter(annotationrelease == "105")
  
  # Merge into the values data frame
  # tvalues = merge(x=tvalues, y=loc.df, by = c("uid"), all=T)
  
  # return the annotated data
  tvalues
  
}

summary.list  = lapply(genes.split, GetSummary)
summary.table = as.data.frame(do.call(rbind, summary.list))
