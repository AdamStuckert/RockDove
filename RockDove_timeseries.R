R

source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
require("dplyr")
require("foreach")
require("data.table")
require("splines")
require("ggthemes")
require("scales")
require("gridExtra")

# load packages required for analyses


library(sleuth)
library(dplyr)
#library(foreach)
#library(doParallel)
library(data.table)
library(splines)
library(ggthemes)
library(scales)
library(gridExtra)


sessionInfo()

memory.limit(32000000)

# get the directory/path for each sample in this study
base_dir <- getwd()
# name of each sample
sample_id <- dir(file.path(base_dir, "kallisto_mappings")) 
# append them to get the location of each samples' quantification
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_mappings", id))
# read in a tab-delimited file with information about samples and treatments
samples <- read.csv("Parental_data_for_analyses.csv", header = TRUE)
# order them to match up with sample paths
# samples <- samples[order(samples$sample),] ##not sure this is needed?
# append paths to the sample dataframe
samples <- dplyr::mutate(samples, path = kal_dirs)

samples

# First, prepare the model by building the design; these need to be in numeric format (they should be already but just in case)
samples$dayas.numeric(samples$day)

spline_design <- model.matrix(formula( ~ ns(samples$day, df = 3) + samples$sex + samples$tissue))
spline_design

##############add annotation data here
ann <- fread() #### add in annotation data once I have it.

# import everything into a sleuth object using the Sleuth package
so <- sleuth_prep(samples, num_cores = 18) #, target_mapping = ann


