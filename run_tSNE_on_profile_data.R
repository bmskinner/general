# Demonstration of how to read a sample tsv of morphology data, extract the relevant
# columns for the angle profile, run a tSNE, cluster, and display the result
library(tidyverse)
library(Rtsne)
library(cluster)
library(dendextend)

# Read in the input data file
data = read.csv("NMA_full_export_demo.csv", sep="\t", header=T, stringsAsFactors = F)

# Take just the angle profile columns
profiles = data %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,1))))
  
# Set the rownames of the data using the cell id and the strain for convenience
rownames(profiles) = data$CellID

# Set the random number generator with a seed for reproducilble results
set.seed(42)

# Run a tSNE on the profile data - the Rtsne function requires a matrix
rtsne_out = Rtsne(as.matrix(profiles), perplexity=20, max_iter=1000)
tsne.values = as.data.frame(rtsne_out$Y)

# Display the results as a scatter plot
ggplot(tsne.values, aes(x=V1, y=V2))+
  geom_point()

#Cluster the tSNE values using hierarchical clustering
# The agnes function is agglomerative nesting
hc = cluster::agnes(tsne.values, method = "ward")

# Make a dendrogram (tree) from the cluster data
dend = as.dendrogram(hc)

# Cut the dendrogram to get 4 separate groups
clusters = dendextend::cutree(dend, 4)

# Assign the group names to the original data
data$agnes = clusters

# Assign the group names to the tSNE results
tsne.values$agnes = clusters

# Plot the tSNE results coloured by cluster
ggplot(tsne.values, aes(x=V1, y=V2, col=agnes))+
  geom_point()
