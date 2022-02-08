# Pilocytic astrocytoma transcriptomic analysis
# Last updated: 08/02/2022

########## Initialize ##########
library(Rtsne)

########## Analysis ##########
# 1) Select 10,000 most variable probes by standard deviation
# 2) The input for the t-SNE calculation is 1-Pearson correlation, weighted by variance