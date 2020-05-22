# This script demonstrates how to add a custom
# chromosome ideogram to a Gviz plot, when the ideogram
# is not available from the UCSC genome.
library(Gviz)
library(GenomicRanges)
library(tidyverse)
library(biomaRt)

genome="susScr11" 
chr   = "chr1"
chr.length = 274330532
start = 263148525
end   = 263717204

range = GRanges(seqnames = chr, ranges=IRanges(start=start, end=end))

# Example ideogram table format; these are not real bands
# Use FLpters with the chromosome length to create band locations
# ideogramData = tribble(
#   ~chrom, ~chromStart,  ~chromEnd, ~name, ~gieStain,
#   "chr1", 1,         50000000, "test1", "gpos25",
#   "chr1", 50000000,  99000000, "test2", "gpos100",
#   "chr1", 100000000, 110000000, "cent", "acen",
#   "chr1", 111000000, 200000000, "test3", "gneg",
#   "chr1", 200000000, 250000000, "test4", "gpos25",
#   "chr1", 250000000, 300000000, "test5", "gpos75"
# )

# Example ideogram table using recorded FLpter measurements
# for SSC1. Coordinates are calculated from chromosome length.
ideogramData = tribble(
  ~chrom,	~name,	~flpter_start,	~flpter_end,	~gieStain,
  "chr1",	"2.5",	0,	0.033101045,	"gneg",
  "chr1",	"2.4",	0.033101045,	0.059233449,	"gpos100",
  "chr1",	"2.3",	0.059233449,	0.104529617,	"gneg",
  "chr1",	"2.2",	0.104529617,	0.173344948,	"gpos100",
  "chr1",	"2.1",	0.173344948,	0.216898955,	"gneg",
  "chr1",	"1.4",	0.216898955,	0.25,	"gpos100",
  "chr1",	"1.3",	0.25,	0.288327526,	"gneg",
  "chr1",	"1.2",	0.288327526,	0.318815331,	"gpos100",
  "chr1",	"1.1",	0.318815331,	0.341463415,	"gneg",
  "chr1",	"cent",	0.341463415,	0.348432056,	"acen",
  "chr1",	"1.1",	0.348432056,	0.363240418,	"gneg",
  "chr1",	"1.2",	0.363240418,	0.385888502,	"gpos100",
  "chr1",	"1.3",	0.385888502,	0.404181185,	"gneg",
  "chr1",	"1.4",	0.404181185,	0.433797909,	"gpos100",
  "chr1",	"1.5",	0.433797909,	0.452090592,	"gneg",
  "chr1",	"1.6",	0.452090592,	0.482578397,	"gpos100",
  "chr1",	"1.7",	0.482578397,	0.519163763,	"gneg",
  "chr1",	"1.8",	0.519163763,	0.568815331,	"gpos100",
  "chr1",	"2.1",	0.568815331,	0.613240418,	"gneg",
  "chr1",	"2.2",	0.613240418,	0.657665505,	"gpos100",
  "chr1",	"2.3",	0.657665505,	0.709930314,	"gneg",
  "chr1",	"2.4",	0.709930314,	0.730836237,	"gpos100",
  "chr1",	"2.5",	0.730836237,	0.770034843,	"gneg",
  "chr1",	"2.6",	0.770034843,	0.799651568,	"gpos100",
  "chr1",	"2.7",	0.799651568,	0.825783972,	"gneg",
  "chr1",	"2.8",	0.825783972,	0.867595819,	"gpos100",
  "chr1",	"2.9",	0.867595819,	0.905052265,	"gneg",
  "chr1",	"2.10",	0.905052265,	0.919860627,	"gpos100",
  "chr1",	"2.11",	0.919860627,	0.93989547,	"gneg",
  "chr1",	"2.12",	0.93989547,	0.960801394,	"gpos100",
  "chr1",	"2.13",	0.960801394,	1,	"gneg"
) %>% mutate(chromStart = chr.length * flpter_start, chromEnd = chr.length * flpter_end)

# Set the options to use the faster European server
options(ucscChromosomeNames=TRUE, Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")

# Custom ideogram track using the data defined above
iTrack = Gviz::IdeogramTrack(genome=genome, chromosome=chr, bands=ideogramData)

# Ideogram track using band data from UCSC; not present for pig
# iTrack = Gviz::IdeogramTrack(genome=genome, chromosome=chr)

# Make a genome axis over the region
gAxis = Gviz::GenomeAxisTrack(range, lwd=2, fontsize=10)

# Examples of annotations
# A gene region track: Genscan genes from UCSC
gTrack = Gviz::UcscTrack(genome=genome, chromosome=chr, track="Genscan Genes", from = start, to = end, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", fill = "#8282d2", name = "Genscan Genes")

# An annotation track: CPGs from UCSC
cpgTrack = Gviz::UcscTrack(genome=genome, chromosome=chr, track="cpgIslandExt", from = start, to = end, trackType = "AnnotationTrack",  start = "chromStart", end = "chromEnd", id = "name", shape = "box", fill = "#006400", name = "CpG Islands")

# Plot the tracks in the desired order
plotTracks(list(iTrack, gAxis, gTrack, cpgTrack ))

# Alternative: fetch gene data from BioMart
host   = "www.ensembl.org"
mart   = useMart("ENSEMBL_MART_ENSEMBL", dataset="sscrofa_gene_ensembl", host=host) #connect
fm     = Gviz:::.getBMFeatureMap()
dTrack = Gviz::BiomartGeneRegionTrack(genome= "Sscrofa11.1",
                                      chromosome = chr, start = start, end = end,
                                      name = "Ensembl", biomart = mart, showId = TRUE,                                               featureMap = fm)

plotTracks(list(iTrack, gAxis, dTrack, cpgTrack ), stacking="squish")
