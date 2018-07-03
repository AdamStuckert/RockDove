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


# get the directory/path for each sample in this study
base_dir <- getwd()
# name of each sample
sample_id <- dir(file.path(base_dir, "kallisto_mappings")) 
# append them to get the location of each samples' quantification
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_mappings", id))
# read in a tab-delimited file with information about samples and treatments
samples <- read.csv("AdamTimeSeriesDocs/Parental_data_for_analyses.csv", header = TRUE)

# order them to match up with sample paths
samples <- samples[order(samples$sample),] 
# append paths to the sample dataframe
samples <- dplyr::mutate(samples, path = kal_dirs)

head(samples)

# First, prepare the model by building the design; these need to be in numeric format (they should be already but just in case)
samples$day <- as.numeric(samples$day)

full_spline_design <- model.matrix(formula( ~ ns(samples$day, df = 3) + samples$sex + samples$tissue))
notime_spline_design <- model.matrix(formula( ~  samples$sex + samples$tissue))
nosex_spline_design <- model.matrix(formula( ~ ns(samples$day, df = 3) + samples$tissue))
notimesex_spline_design <- model.matrix(formula( ~ samples$tissue))

##############add annotation data here
ann <- fread("AdamTimeSeriesDocs/name.mapping.withcodes.parental.csv", header = FALSE) #### add in annotation data once I have it.
colnames(ann) <- c("target_id", "peptide_id", "gene_number", "gene_symbol")

# import everything into a sleuth object using the Sleuth package
so <- sleuth_prep(samples, target_mapping = ann, num_cores = 24)
so <- sleuth_fit(so, formula = full_spline_design, fit_name = "full") 
so <- sleuth_fit(so, formula = notime_spline_design, fit_name = "notime")
so <- sleuth_fit(so, formula = nosex_spline_design, fit_name = "nosex")
so <- sleuth_fit(so, formula = notimesex_spline_design, fit_name = "notimesex")

# print the models
models(so)

# run a likelihood ratio test between the full and the model without time. This basically tests whether the inclusion of the time points explains the data better than the reduced model with just the lane. Transcripts better explained by the inclusion of time should be considered differentially expressed over time. Qvalues are corrected for multiple comparisons.
so_lrt <- sleuth_lrt(so, "notime", "full")

# save the model results to a data frame
lrt_results <- sleuth_results(so_lrt, 'notime:full', test_type = 'lrt')

# how many transcripts are differentially expressed if we use a cut off of a = 0.05?
table(lrt_results[,"qval"] < 0.05)

# We are interested in time series figures for a number of key parental care genes. 

genes <- read.csv("AdamTimeSeriesDocs/parental_care_genes_of_interest.csv", header = TRUE)

## search through the results for anything that annotated in the full or xenopus annotation to a color gene
lrt_par_genes <- lrt_results %>% filter(gene_number %in% genes$gene_number)

dir.create("colorgenespline-figures")

# escape X11!
options(bitmapType='cairo')

# Now extract only the data of interest.
for (i in 1:nrow(lrt_par_genes)){
  transcript <- lrt_par_genes[i,"target_id"]
  gene <- lrt_par_genes[i,"gene_symbol"]
  qval <- lrt_par_genes[i,"qval"]
  tmp <- so$obs_norm %>% dplyr::filter(target_id == transcript)  ### These are normalized values I think!
  tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
  tmp
  
# subset out the three tissues for the graph...

tmp.hyp <- tmp %>% dplyr::filter(tissue == "hypo")
tmp.pit <- tmp %>% dplyr::filter(tissue == "pit")
tmp.gonad <- tmp %>% dplyr::filter(tissue == "gonad")
  
a <- ggplot(tmp.hyp, aes(x=day, y=est_counts)) + geom_point(aes(size = 3, color = sex)) + scale_colour_manual(values = c("pink", "blue"), guide = guide_legend(title = "Sex", override.aes = list(size=4))) +  geom_smooth(method = loess, aes(group = sex, color = sex)) + ggtitle(paste0(gene)) + guides(size=FALSE) + theme_bw()
   ggsave(paste0("colorgenespline-figures/", gene, "_hypothalamus.png"), width = 6.81, height = 3.99)

b <- ggplot(tmp.pit, aes(x=day, y=est_counts)) + geom_point(aes(size = 3, color = sex)) + scale_colour_manual(values = c("pink", "blue"), guide = guide_legend(title = "Sex", override.aes = list(size=4))) +  geom_smooth(method = loess, aes(group = sex, color = sex)) + ggtitle(paste0(gene)) + guides(size=FALSE) + theme_bw()
   ggsave(paste0("colorgenespline-figures/", gene, "_pituitary.png"), width = 6.81, height = 3.99)

c <- ggplot(tmp.gonad, aes(x=day, y=est_counts)) + geom_point(aes(size = 3, color = sex)) + scale_colour_manual(values = c("pink", "blue"), guide = guide_legend(title = "Sex", override.aes = list(size=4))) +  geom_smooth(method = loess, aes(group = sex, color = sex)) + ggtitle(paste0(gene)) + guides(size=FALSE) + theme_bw()
   ggsave(paste0("colorgenespline-figures/", gene, "_gonad.png"), width = 6.81, height = 3.99)
}

# I need to make 


# saved the slueth run:
sleuth_save(so, "Preliminary_sleuth_run")