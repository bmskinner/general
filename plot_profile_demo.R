# Demonstration of how to read a sample tsv of morphology data, extract the relevant
# columns for the angle profile, and plot all profiles based on their dataset
library(tidyverse)

# Read in the input data file
data = read.csv("NMA_full_export_demo.csv", sep="\t", header=T, stringsAsFactors = F)

# Extract the profile data and format it suitable for plotting
profile.data = data %>% ungroup() %>% dplyr::select(Dataset, matches("Angle_profile_\\d+$")) %>%
  tidyr::gather(Profile_position, value, -Dataset) %>%
  dplyr::group_by(Profile_position, Dataset) %>%
  dplyr::mutate(Position = as.numeric(gsub("Angle_profile_", "", Profile_position)),
                Median = median(value),
                Q25 = quantile(value, 0.25),
                Q75 = quantile(value, 0.75)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-value, -Profile_position) %>%
  dplyr::distinct()

# Make a chart showing the profiles for each sample
ggplot(profile.data, aes(x=Position, y=Median, fill=Dataset))+
  geom_hline(yintercept=180)+
  geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.5)+
  geom_line(aes(col=Dataset))+
  xlab("Position")+
  ylab("Angle")+
  ylim(50,250)+
  facet_wrap(~Dataset)+
  theme_classic() +
  theme(legend.position = "top")
