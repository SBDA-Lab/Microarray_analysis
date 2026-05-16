# Microarray Analysis Pipeline - Shiny Interface

## Overview

The **Microarray Analysis Pipeline** is a Shiny-based web application that provides an interactive interface for analyzing microarray data. It allows users to:

* **Upload** CEL files and assign sample groups.
* **Perform** quality control (QC) analysis with various diagnostic plots.
* **Run** differential expression analysis and visualize results.
* **Download** results for further downstream analysis.

This pipeline is designed to simplify the microarray data analysis process for researchers and bioinformaticians alike.

---

## Features

* **Data Upload:** Supports uploading multiple CEL files and assigning experimental groups.
* **Quality Control:** Generates RNA degradation plots, boxplots, MA plots, and density plots.
* **Differential Expression Analysis:** Displays results in an interactive table format alongside dynamic volcano plots and heatmaps.
* **Downloadable Results:** Users can download processed data and analysis results in standard CSV format.
* **User-Friendly Interface:** Built using `shiny` and `shinydashboard`, making it fully accessible via any modern web browser.

---

## Installation & Dependencies

The application requires **R** and several Bioconductor and CRAN packages.

### Required R Packages

The following core UI and visualization packages are utilized:

* `shiny`
* `shinydashboard`
* `DT`
* `plotly`

> [!IMPORTANT]
> Ensure you have `BiocManager` installed before running the application to handle any genomic dependencies.

To install the necessary dependencies manually, execute the following script in your R console:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("shiny", "shinydashboard", "DT", "plotly")

for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package)
    }
}
```

---

## How to Run the Application

1. Install dependencies using the script provided above if they are not already present.
2. Ensure the core logic script `microarray_pipeline.R` is in the exact same directory as `app.R`.
3. Open R or RStudio and set the working directory to the folder containing the app files.
4. Run the following commands in your R console:

```R
library(shiny)
runApp("app.R")
```

5. The application will automatically launch in your default web browser.

---

## Usage Guide

### 1. Upload Data

* Navigate to the **Upload Data** tab in the sidebar.
* Upload your raw `.CEL` files.
* Assign corresponding experimental groups (e.g., *Control* vs. *Treatment*).
* Specify your preferred output directory path.
* Click **"Run Analysis"** to initialize the pipeline.

### 2. Quality Control

* Navigate to the **Quality Control** tab.
* Inspect the automatically generated diagnostic plots:
  * RNA degradation plots
  * Intensity box plots
  * MA plots
  * Density plots

### 3. View and Download Results

* Head over to the **Analysis Results** tab.
* Select your target contrast/comparison framework to view:
  * An interactive **Differential Expression Table** (filterable via `DT`).
  * A dynamic **Volcano Plot** to identify statistically significant fold changes.
  * A **Heatmap** clustering the top 25 differentially expressed genes.
* Click the download button to export your results matrix as a `.csv` file.

---

## File Structure

```text
Microarray-Analysis-Pipeline/
│
├── app.R                       # Main Shiny user interface & server execution logic
├── microarray_pipeline.R       # Core bioinformatics analysis functions
│
├── data/                       # (Optional) Directory for storing example raw datasets
└── microarray_results/          # Generated automatically upon pipeline completion
```

---

## Troubleshooting

* **Application doesn't launch?**
  * Verify that all packages in the dependency block installed successfully without compilation errors.
  * Run `sessionInfo()` in your console to check your current R version and package attachments.

* **No results appear after clicking "Run Analysis"?**
  * Check your R console layout for underlying errors.
  * Verify that all `.CEL` files are valid, uncorrupted, and properly mapped to a specific experimental group.
  * Double-check that `microarray_pipeline.R` resides in the same root directory as `app.R`.

* **Plots not displaying or interactive features broken?**
  * Ensure your web browser has JavaScript enabled.
  * Re-verify your `plotly` installation by calling `library(plotly)` directly in your console to check for missing system dependencies.
