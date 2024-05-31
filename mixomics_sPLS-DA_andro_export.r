########################################################################################
### MixOmics sPLS-DA - Androgen-microbes-ENS GUS Gene Intensity
########################################################################################

##### Usage
# Rscript mixomics_sPLS-DA_andro_export.r

##### prep environment
### if auto install fails:
# install mixomics at https://github.com/mixOmicsTeam/mixOmics
# manually load packages instead if needed
# suppressMessages(library(mixOmics))
# suppressMessages(library(dplyr))
# suppressMessages(library(stringr))

### auto install

rm(list=ls()) # Clean workspace 
options(repos = "https://cloud.r-project.org") # Set the mirror

# Define required packages
required_packages <- c("mixOmics", "dplyr","stringr")

# Identify packages that are and are not installed
installed <- installed.packages()[,"Package"]
not_installed <- setdiff(required_packages, installed)

# Install missing packages
if(length(not_installed) > 0) {
  if('mixOmics' %in% not_installed){ # install mixOmics
    # install BiocManager if not installed
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::valid()
    BiocManager::install('mixOmics')
  installed <- installed.packages()[,"Package"]
  not_installed <- setdiff(required_packages, installed)
  install.packages(not_installed, dependencies = TRUE)
  }
}

# Load packages invisibly
invisible(suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE)))

# Return an error if required packages are not loaded, invisible above hides any such errors
loaded_packages <- character(length(required_packages))
for (i in seq_along(required_packages)) {
  loaded_packages[i] <- tryCatch({
    library(required_packages[i], character.only = TRUE)
    required_packages[i]
  }, error = function(e) {
    message(paste("Error loading package:", required_packages[i], "-", conditionMessage(e)))
    NA
  })
}

########################################################################################
##### main
########################################################################################

################################################################################
# mixomics example
# data(srbct) # extract the small round bull cell tumour data
# X <- srbct$gene # use the gene expression data as the X matrix
# Y <- srbct$class # use the class data as the Y matrix
# print(X)
# print(dim(X))
# print(Y)
# print(length(Y))
# PLS-DA
# result.plsda.srbct <- plsda(X, Y) # run the method
# plotIndiv(result.plsda.srbct) # plot the samples
# plotVar(result.plsda.srbct) # plot the variables 
################################################################################

# set paths
in_csv <- "andro_hms_rep_summary_ALL_tax_GUS_full70.csv"
outdir <- "./andro_MGX_full_md_ALL_GUS_genes/"

md_df <- read.csv(in_csv,check.names=FALSE) # read data

row_size <- .47 # set rowsize
marginVector <- c(7, 10) # set margin vector
legend_size <- 1 

# Find the index of the split column
split_column <- "14_2_MGG44437_Loop_1"
split_column_index <- which(names(md_df) == split_column)

md_col <- "treatment" # get md col of interest

# create outdir
if (!dir.exists(outdir)){
  dir.create(outdir)
}

md_df <- md_df %>% arrange(str_to_lower(!!sym(md_col))) # arrange df based on alphabetical order of entries in md col (important for labeling)

X <- md_df[, split_column_index:ncol(md_df)] # use split_column to get X matrix
Y <- factor(md_df[[md_col]]) # use the class data as the Y matrix
pdf(paste0(outdir,md_col,".pdf")) # specify ofnm and open pdf

# run method - PLS-DA
result.splsda <- splsda(X, Y)

### CIM
# set the styling of the legend
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = md_col, # legend title
        cex = legend_size) # legend size

cim(result.splsda,
    row.cex=row_size,
    row.names = FALSE, 
    legend = legend,
    margins=marginVector, # first controls x axis dist from edge, 2nd is y (or opp if trans == TRUE)
    row.sideColors = color.mixo(Y),
    transpose = TRUE)
