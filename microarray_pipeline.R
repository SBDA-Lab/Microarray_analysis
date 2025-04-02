# Microarray Data Analysis Pipeline
# This script processes and analyzes microarray data from .cel files

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
required_packages <- c("affy", "limma", "affyPLM", "AnnotationDbi", 
                      "hgu133plus2.db", "genefilter", "RColorBrewer",
                      "pheatmap", "shiny")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE))
    BiocManager::install(package)
}

# Load required libraries
library(affy)
library(limma)
library(affyPLM)
library(AnnotationDbi)
library(hgu133plus2.db)
library(genefilter)
library(RColorBrewer)
library(pheatmap)

# Function to create output directories
create_output_dirs <- function(base_dir) {
  dirs <- c("QC", "Normalized_Data", "DE_Analysis", "Plots")
  for (dir in dirs) {
    dir.path <- file.path(base_dir, dir)
    if (!dir.exists(dir.path)) {
      dir.create(dir.path, recursive = TRUE)
    }
  }
  return(list(
    qc = file.path(base_dir, "QC"),
    norm = file.path(base_dir, "Normalized_Data"),
    de = file.path(base_dir, "DE_Analysis"),
    plots = file.path(base_dir, "Plots")
  ))
}

# Function to read CEL files
read_cel_files <- function(cel_dir) {
  # Read all CEL files in the directory
  cel_files <- list.files(cel_dir, pattern = "\\.cel$|\\.CEL$", full.names = TRUE)
  if (length(cel_files) == 0) {
    stop("No CEL files found in the specified directory")
  }
  
  # Read the CEL files
  raw_data <- ReadAffy(filenames = cel_files)
  return(raw_data)
}

# Function for quality control
perform_qc <- function(raw_data, output_dirs) {
  # RNA degradation plot
  pdf(file.path(output_dirs$qc, "RNA_degradation.pdf"))
  deg <- AffyRNAdeg(raw_data)
  plotAffyRNAdeg(deg)
  dev.off()
  
  # Box plots of raw data
  pdf(file.path(output_dirs$qc, "Raw_Intensity_Boxplot.pdf"))
  boxplot(raw_data, main="Raw Intensity Values", las=2)
  dev.off()
  
  # MA plots
  pdf(file.path(output_dirs$qc, "MA_plots.pdf"))
  MAplot(raw_data, pairs=TRUE)
  dev.off()
  
  # Density plots
  pdf(file.path(output_dirs$qc, "Density_plots.pdf"))
  hist(raw_data, main="Density Plot of Raw Intensity Values")
  dev.off()
}

# Function for data normalization
normalize_data <- function(raw_data, output_dirs) {
  # RMA normalization
  eset_rma <- rma(raw_data)
  
  # Save normalized expression values
  exprs_data <- exprs(eset_rma)
  write.csv(exprs_data, 
            file = file.path(output_dirs$norm, "normalized_expression_values.csv"))
  
  # Create post-normalization QC plots
  pdf(file.path(output_dirs$qc, "Post_Normalization_Boxplot.pdf"))
  boxplot(exprs_data, main="Normalized Expression Values", las=2)
  dev.off()
  
  return(eset_rma)
}

# Function to create heatmap of top 25 genes
create_heatmap <- function(eset, results, groups, output_dirs) {
  # Get top 25 genes by adjusted p-value
  top_genes <- head(order(results$adj.P.Val), 25)
  
  # Extract expression data for top genes
  expr_data <- exprs(eset)[top_genes, ]
  
  # Get gene symbols for labeling
  gene_symbols <- mapIds(hgu133plus2.db,
                        keys=rownames(expr_data),
                        column="SYMBOL",
                        keytype="PROBEID",
                        multiVals="first")
  
  # Create annotation for samples
  sample_anno <- data.frame(
    Group = groups,
    row.names = colnames(expr_data)
  )
  
  # Create color schemes
  ann_colors <- list(
    Group = c(control = "#1B9E77", treatment = "#D95F02")
  )
  
  # Create and save heatmap
  pdf(file.path(output_dirs$plots, "Top25_Genes_Heatmap.pdf"))
  pheatmap(expr_data,
           labels_row = gene_symbols,
           annotation_col = sample_anno,
           annotation_colors = ann_colors,
           scale = "row",
           main = "Top 25 Differentially Expressed Genes",
           fontsize_row = 8,
           fontsize_col = 8)
  dev.off()
  
  # Return the data for potential interactive plotting
  return(list(
    expr_data = expr_data,
    gene_symbols = gene_symbols,
    sample_anno = sample_anno
  ))
}

# Function for differential expression analysis
perform_de_analysis <- function(eset, groups, output_dirs) {
  # Create design matrix
  design <- model.matrix(~0 + factor(groups))
  colnames(design) <- levels(factor(groups))
  
  # Fit linear model
  fit <- lmFit(eset, design)
  
  # Create contrasts for all pairwise comparisons
  group_levels <- levels(factor(groups))
  contrasts <- makeContrasts(
    contrasts = combn(group_levels, 2, paste, collapse="-"),
    levels = design
  )
  
  # Compute contrasts
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  # Get results for each comparison
  for(i in 1:ncol(contrasts)) {
    comparison <- colnames(contrasts)[i]
    results <- topTable(fit2, coef=i, number=Inf)
    
    # Add gene symbols
    results$Symbol <- mapIds(hgu133plus2.db,
                           keys=rownames(results),
                           column="SYMBOL",
                           keytype="PROBEID",
                           multiVals="first")
    
    # Save results
    write.csv(results,
              file = file.path(output_dirs$de,
                             paste0("DE_results_", comparison, ".csv")))
    
    # Create volcano plot
    pdf(file.path(output_dirs$plots,
                  paste0("Volcano_plot_", comparison, ".pdf")))
    plot(results$logFC, -log10(results$adj.P.Val),
         main=paste("Volcano Plot -", comparison),
         xlab="log2 Fold Change",
         ylab="-log10 adjusted p-value",
         pch=20)
    dev.off()
    
    # Create heatmap for top 25 genes
    heatmap_data <- create_heatmap(eset, results, groups, output_dirs)
  }
}

# Main pipeline function
run_microarray_pipeline <- function(cel_dir, output_dir, groups) {
  # Create output directories
  output_dirs <- create_output_dirs(output_dir)
  
  # Read CEL files
  message("Reading CEL files...")
  raw_data <- read_cel_files(cel_dir)
  
  # Perform quality control
  message("Performing quality control...")
  perform_qc(raw_data, output_dirs)
  
  # Normalize data
  message("Normalizing data...")
  eset_normalized <- normalize_data(raw_data, output_dirs)
  
  # Perform differential expression analysis
  message("Performing differential expression analysis...")
  perform_de_analysis(eset_normalized, groups, output_dirs)
  
  message("Pipeline completed successfully!")
}

# Example usage:
# cel_dir <- "path/to/cel/files"
# output_dir <- "path/to/output"
# groups <- factor(c("control", "control", "treatment", "treatment"))  # Must match the order of CEL files
# run_microarray_pipeline(cel_dir, output_dir, groups) 