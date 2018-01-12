# This script demonstrates how to add a custom
# chromosome ideogram to a Gviz plot, when the ideogram
# is not available from the UCSC genome.
library(Gviz)
library(GenomicRanges)
library(tidyverse)

genome="susScr11" 
chr   = "chr1"
chr.length = 274330532
start = 263148525
end   = 263717204

range = GRanges(seqnames = chr, ranges=IRanges(start=start, end=end))

# Example ideogram table format; these are not real bands
# Use FLpters with the chromosome length to create band locations
ideogramData = tribble(
  ~chrom, ~chromStart,  ~chromEnd, ~name, ~gieStain,
  "chr1", 1,         50000000, "test1", "gpos25",
  "chr1", 50000000,  99000000, "test2", "gpos100",
  "chr1", 100000000, 110000000, "cent", "acen",
  "chr1", 111000000, 200000000, "test3", "gneg",
  "chr1", 200000000, 250000000, "test4", "gpos25",
  "chr1", 250000000, 300000000, "test5", "gpos75"
)

# Set the options to use the faster European server
options(ucscChromosomeNames=TRUE, Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")

# Custom ideogram track using the data defined above
iTrack = Gviz::IdeogramTrack(genome=genome, chromosome=chr, bands=ideogramData)

# Make a genome axis over the region
gAxis = Gviz::GenomeAxisTrack(range, lwd=2, fontsize=10)

# Examples of annotations
# A gene region track: Genscan genes from UCSC
gTrack = Gviz::UcscTrack(genome=genome, chromosome=chr, track="Genscan Genes", from = start, to = end, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", fill = "#8282d2", name = "Genscan Genes")

# An annotation track: CPGs from UCSC
cpgTrack = Gviz::UcscTrack(genome=genome, chromosome=chr, track="cpgIslandExt", from = start, to = end, trackType = "AnnotationTrack",  start = "chromStart", end = "chromEnd", id = "name", shape = "box", fill = "#006400", name = "CpG Islands")

# Plot the tracks in the desired order
plotTracks(list(iTrack, gAxis, gTrack, cpgTrack ))


# Alternative: fetch gene data from BioMart - but cannot be used with UCSC tracks, including ideogram tracks

# host = "www.ensembl.org" 
# mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="sscrofa_gene_ensembl", host=host) #connect
# fm   = Gviz:::.getBMFeatureMap()
# fm["symbol"] = "external_gene_id"
# dTrack = Gviz::BiomartGeneRegionTrack(genome= "Sscrofa11.1", 
#                                       chromosome = chr, start = start, end = end,
#                                       name = "ENSEMBL", biomart = mart, showId = TRUE)
