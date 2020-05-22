# Custom gene annotation of Affy mouse arrays
# using Ensembl and MGI sources


# source("https://bioconductor.org/biocLite.R")
# biocLite("AnnotationDbi")
# biocLite("mogene21sttranscriptcluster.db")
# biocLite("oligo")
# biocLite("limma")
# biocLite("genefilter")
# biocLite("Biobase")
# biocLite("biomaRt")
library(dplyr)
library(AnnotationDbi)
library(mogene21sttranscriptcluster.db) # Annotations for mouse array
library(pd.mogene.2.1.st) # Array probes
library(tidyr)
library(oligo)
library(biomaRt)
library(reshape2)

analysis.folder = "//path/to/data/using/windows/network/share/"
setwd(analysis.folder)
cat("Analysis working directory: ", analysis.folder, "\n")

####################
#
# Define gene annotation
# sources
#
####################

# Get the MGI annotations
annotation.url  = "http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"
annotation.file = paste0(analysis.folder, "/MGI_EntrezGene.rpt")
if(!file.exists(annotation.file)){
  download.file(annotation.url, annotation.file) # ensure we are up to date
}

# Read the annotations file. It does not come with a header.
mgi.annotations = read.csv(annotation.file, header=F, sep="\t")
colnames(mgi.annotations) = c("MGI_Marker_Accession_ID","Marker_Symbol","Status","Marker_Name","cM_Position","Chromosome","Type","Secondary_Accession_IDs","Entrez_Gene_ID","Synonyms","Feature_Types","Genome_Coordinate_Start","Genome_Coordinate_End","Strand","BioTypes")

# Fetch the transcript annotations for the array
k <- keys(mogene21sttranscriptcluster.db,keytype="PROBEID")
dbi.annotations = AnnotationDbi::select(mogene21sttranscriptcluster.db, keys=k, columns=c("SYMBOL","GENENAME", "ENTREZID"), keytype="PROBEID") %>% dplyr::distinct()

# Fetch ensembl annotations from BioMart
ens.mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl") # Where to look
ens.filters    = c("affy_mogene_2_1_st_v1") # Stuff to look up
ens.attributes = c("affy_mogene_2_1_st_v1", "mgi_id", "ensembl_gene_id", "external_gene_name", "description", "gene_biotype") # Columns to get back
# Use: listAttributes(ens.mart) to see all available columns to fetch

####################
#
# Add annotations
#
####################

addExtraGeneAnnotations = function(data){
  # Add ensembl data
  probes = unique(data$probesetid)
  result = getBM(mart=ens.mart, attributes=ens.attributes, filters=ens.filters, values = probes )
  data   = merge(data, result, by.x = "probesetid", by.y = "affy_mogene_2_1_st_v1", all.x = T, all.y=F)

  # Add annotations direct from the mogene21sttranscriptcluster.db
  data = merge(data, dbi.annotations, by.x = "probesetid", by.y = "PROBEID", all.x = T, all.y=F)

  # Add the MGI information
  data = merge(data, mgi.annotations, by.x = "mgi_id", by.y = "MGI_Marker_Accession_ID", all.x = T, all.y=F)

  data
}

# Example: to perform annotation:
annotated.data.frame = addExtraGeneAnnotations(data.frame.with.probesetid.column)

# Example: to filter for protein coding genes:
annotated.data.frame = annotated.data.frame %>%
dplyr::filter( gene_biotype == "protein_coding")

# Example: to select specific columns for export
export.cols = list("probesetid", "log2_average", "log2_diff", "ensembl_gene_id", "TranscriptID", "Symbol", "Desc", "Chromosome", "Genome_Coordinate_Start","Genome_Coordinate_End","Strand" )
desired.columns = annotated.data.frame %>%
dplyr::select_( .dots = export.cols )


# Example: to select probesets that hit only one ensembl gene
single.hit.genes =  annotated.data.frame %>%
dplyr::distinct() %>%
dplyr::group_by(probesetid) %>%
dplyr::mutate( Hits = n() ) %>%
dplyr::filter( Hits==1) %>%
dplyr::select(-Hits)
