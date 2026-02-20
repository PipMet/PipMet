
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PipMet

<!-- badges: start -->

The **PipMet** package was developed to perform end-to-end processing of
metabolomic-based GC-MS data, with automated generation of high-quality
figures throughout the workflow. All user inputs are obtained through
pop-up windows.  


<!-- badges: end -->

## Installation

### Dependencies installation

Before installing `PipMet`, make sure that all required dependencies are
installed.  
The following code will automatically check for and install all
necessary **Bioconductor** and **CRAN** packages:

```r
# 1. Install BiocManager and align with Bioconductor 3.22 (Required for R 4.5)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.22", ask = FALSE, update = TRUE)

# 2. Install Devtools
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# 3. Install Bioconductor specific dependencies
bioc_pkgs <- c("xcms", "MSnbase", "CluMSID", "metaMS", "BiocParallel", 
               "Biobase", "ProtGenerics", "CAMERA", "NormalyzerDE")
BiocManager::install(bioc_pkgs, ask = FALSE, update = TRUE)

# 4. Install CRAN dependencies
cran_pkgs <- c("svDialogs", "pheatmap", "ddpcr", "webchem", "fritools", "pracma")
new_pkgs <- cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

```
### PipMet installation

You can install the released version of PipMet from
[GitHub](https://github.com) with:

``` r {eval=FALSE}
remotes::install_github("PipMet/PipMet", force = TRUE)
```

## Example

The package is constituted of two main functions with pre-set parameters
and algorithms for GC-MS data processing. The `workData()` function
reads, treats and process GC-MS sample data, with metabolite
identification and quantification.

The package was thought to be as friendly-user as possible. Therefore, when information is needed,
pop-ups will appear to collect input.

To run an example use

``` r {eval=FALSE}
library(PipMet)
result <- workData(example = TRUE)
```

To process your dataset use

``` r {eval=FALSE}
library(PipMet)
result <- workData()
```

For more information, see the manual: https://github.com/PipMet/PipMet/blob/main/Manual_PipMet.pdf.
