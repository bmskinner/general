# general
This repo contains generic scripts that are not tied to any particular project. They have been useful for demonstrations or used in general purpose across many projects.

# rentrez_demo.r

This shows how to use the rentrez package to access the NCBI Entrez api and fetch summary information from the gene database.


# plot_multiple_ggplots_demo.r

This shows how to grid and align ggplots generated via a function. It uses the sample dataset in plot_multiple_ggplots_demo.csv

# addGeneAnnotations.r

This uses BioMart and MGI annotations to improve standard Affymetrix microarray gene annotation, allowing probesets that are incorrectly annotated to be detected for downstream scripts. An example is probesetid 17215370, which is (as of 2017-08) annotated as the protein-coding gene Atg16l1, but actually maps to a scaRNA within an intron (Gm25395). BioMart contains Ensembl's own mappings of the sequence features associated with the oligos, so the scaRNA id maps to the probesetid. Mismatches between the expected and Ensembl biotype can be used as an extra filter. 

# Gviz_demo.R

This shows how to generate a custom ideogram track for use when plotting UCSC gene data
