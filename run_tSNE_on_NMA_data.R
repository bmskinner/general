# Demonstration of how to read a sample tsv of morphology data, extract the relevant
# columns for the angle profile, run a tSNE, cluster, and display the result

# Load the required packages
library(tidyverse)
library(Rtsne)
library(cluster)
library(dendextend)

# Define the input data file
data.file = "plot_multiple_ggplots_demo.csv"

# Read in the input data file
data = read.csv(data.file, sep="\t", header=T, stringsAsFactors = F)

# Take just the columns whose name starts with "Angle_profile_"
profiles = data %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,1))))
  
  
# Example of making a chart
# Plot the areas as a volin chart plus with a bar-chart
# Rotate the x-axis labels to 45 degrees to avoid overlaps
ggplot(data, aes(x=Dataset, y=Area_square_pixels))+
  geom_violin()+
  geom_boxplot(width=0.2, outlier.shape = NA)+ # Outliers are shown as dots by default - can make chart messier
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Theme functions allow all aspects of the plot to be changed
  

# Set the rownames of the data using the cell id for convenience
rownames(profiles) = data$CellID

# Set the random number generator with a seed for reproducilble results
set.seed(42)

# Run a tSNE on the profile data - the Rtsne function requires a matrix
rtsne_out = Rtsne(as.matrix(profiles), perplexity=100, max_iter=1000)
tsne.values = as.data.frame(rtsne_out$Y)

# Display the results as a scatter plot
ggplot(tsne.values, aes(x=V1, y=V2))+
  geom_point()

# Cluster the tSNE values using hierarchical clustering
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

