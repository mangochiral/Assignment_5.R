# Create R package
library(devtools)
devtools::install_local("Assignment_5")

# Package documentation
#' @title myDEGPkg: Identify Differentially Expressed Genes
#' @description This package provides functions for identifying differentially expressed genes (DEGs) from gene expression data.
#' @author Your Name
#' @import GEOdata, limma, cluster, ggplot2, ReactomePA, fgsea
#' @examples
#' \dontrun{
#' # Load data
#' library(GEOdata)
#' gset <- getGEO("GSE33126")
#'
#' # Identify DEGs
#' library(Assignment_5)
#'
#'}
#' @export

# Function for identifying Log Normalize
norm_log <- function(data) {
  # Log-normalize data
  exprs_norm <- log10(gene_exprs)+1
 
 # Perform differential expression analysis
 library(limma)
design <- model.matrix(~0+metadata_modified$group)
design
colnames(design) <- c(colnames(metadata_modified))
fit <- lmFit(gene_exprs, design)
contrast_matrix <- makeContrasts(normal - tumor, levels = design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)
DEGs <- topTable(fit_ebayes, coef = 1, number = Inf)

  # Return DEGs
  return(DEGs)
}

# Vignette for the package
