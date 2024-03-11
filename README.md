# MED 263 Final Project: Microbiome Analysis in R

> The goal of this project is to learn basic techniques for analyzing microbiome data and investigate if human gut, palm, and tongue microbiomes differ in terms of composition and microbial diversity. 

## Introduction:
Over 90% of the human body is made up of non-human cells, mostly consisting of bacteria (~99%)<sup>1,2</sup>. These microorganisms, collectively termed the human microbiota<sup>3</sup>, have been shown to impact various host processes, including malnutrition and obesity<sup>4</sup>, inflammatory bowel disease<sup>5</sup>, and drug metabolism<sup>6</sup>. Therefore, understanding microbes and their composition in microbial communities at a genome level has important applications in the mechanisms of human disease and research involving bacteria. Using 16S ribosomal RNA (rRNA) sequencing, profiling microbiome data can identify species present in the samples and determine their relative abundance, enabling the measurement of microbial diversity of samples and the comparison of microbial diversity across samples<sup>7,8,9</sup>. We plan to host a tutorial covering bioinformatic methods in microbiome analysis. Using bioinformatic tools we can classify samples based on their microbial diversity and look at changes of aggregated measures under different treatments to the sampled microbe community.

# Methods:
This tutorial uses a subset of publicly available hypervariable region 4 (V4) 16S rRNA dataset from Caporaso et al. (2011)<sup>10</sup>. The subset data includes human microbial samples taken from gut, left palm, right palm, and tongue. Three tsv files (feature-table.tsv, sample-metadata.tsv, taxonomy.tsv) will be made available through GitHub. R packages required in this tutorial include phyloseq<sup>11</sup>, tidyverse<sup>12</sup>, vegan<sup>13</sup>, DESeq2<sup>14</sup>, and ComplexHeatmap<sup>15</sup>.
The class will first load data, then perform basic data inspection and cleaning to get familiar with high dimensional, compositional, and sparse microbiome data. They will then build a “phyloseq” object to perform rarefaction to minimize the bias from varying sampling depths. Furthermore, they will use “phyloseq” built-in functions to calculate alpha and beta diversity of samples, make principal coordinates analysis (PCoA) plots, and conduct statistical tests to compare these metrics between samples taken from different environments. To visualize the microbial composition of each sample, they will learn how to create relative abundance plots. Finally, the class will learn how to identify differentially abundant species using DESeq2. Through the use of ComplexHeatmap, they can visualize microbes that are differentially abundant in each environment.

# Task List:

* Data cleaning & formatting
* Rarefaction
* Calculation of alpha & beta diversity
* PCoA plots
* PERMANOVA to compare samples from different environments
* Relative abundance profiles on genus and species level
* Differential abundance analysis & heatmaps

# Data Availability
Access to the Data files to follow along the tutorial can be accessed here: [Google Drive Folder Link](https://drive.google.com/drive/folders/1L4oR2CkqmCIcGuHUSuyaQnW9GlauV3yo?usp=drive_link).

# Prerequisites and Installation
## If using RStudio:
The packages can be accessed by running this code chunk:

```{r}
p1 <- c("tidyverse", "vegan", "BiocManager")
p2 <- c("phyloseq", "DESeq2", "ComplexHeatmap")
load_package <- function(p) {
if (!requireNamespace(p, quietly = TRUE)) {
ifelse(p %in% p1,
install.packages(p, repos = "http://cran.us.r-project.aorg/"),
BiocManager::install(p))
}
library(p, character.only = TRUE, quietly = TRUE)
}
invisible(lapply(c(p1,p2), load_package))
```

Est. run time for package installation: .5-1 hours

## If using R in Jupyter Notebook:
Use the following code chunk:
```{r}
install.packages("R.utils")
p1 <- c("tidyverse", "vegan", "BiocManager")
p2 <- c("phyloseq", "DESeq2", "ComplexHeatmap")
load_package <- function(p) {
if (!requireNamespace(p, quietly = TRUE)) {
1
ifelse(p %in% p1,
install.packages(p, repos = "http://cran.us.r-project.org/"),
BiocManager::install(p))
}
library(p, character.only = TRUE, quietly = TRUE)
}
invisible(lapply(c(p1,p2), load_package))
```

You then should be able to access these libraries without issue:
```{r}
library(IRdisplay)
library(tidyverse)
library(vegan)
library(BiocManager)
library(DESeq2)
library(ComplexHeatmap)
library(phyloseq)
```

# Troubleshooting
## Using Rstudio:

1. If time is an issue, remove ANCOMBC as a package to be installed from the vector p2
2. If asked to update packages, select all (“a”)

## Using Jupyter Notebook:

1. If running into issues with the phyloseq package, run the following command to load dependencies of
the package:

```{r}
install.packages("devtools")
library(devtools)
install_github("igraph/rigraph")
BiocManager::install("phyloseq")
```

# Contributors:
This tutorial by Yan Hui was adapted by Sherlyn Weng, Kate Johnson, Shengjia Tu, and Avery Williams for the MED 263 Final Project.

# Acknowledgments:
We thank the MED 263 course instructors and guest lecturers for their support and guidance throughout the winter quarter of 2024.

# References:

1. Grice, E. A., & Segre, J. A. (2012). The human microbiome: our second genome. Annual review of genomics and human genetics, 13, 151-170.

2. Sender, R., Fuchs, S., & Milo, R. (2016). Revised estimates for the number of human and bacteria cells in the body. PLoS biology, 14(8), e1002533.

3. Ursell, L. K., Metcalf, J. L., Parfrey, L. W., & Knight, R. (2012). Defining the human microbiome. Nutrition reviews, 70(suppl_1), S38-S44.

4. Fan, Y., & Pedersen, O. (2021). Gut microbiota in human metabolic health and disease. Nature Reviews Microbiology, 19(1), 55-71.

5. Halfvarson, J. et al. (2017). Dynamics of the human gut microbiome in inflammatory bowel disease. Nat Microbiol 2: 17004.

6. Enright, E. F., Joyce, S. A., Gahan, C. G., & Griffin, B. T. (2017). Impact of gut microbiota-mediated bile acid metabolism on the solubilization capacity of bile salt micelles and drug solubility. Molecular Pharmaceutics, 14(4), 1251-1263.

7. Caporaso, J. G. et al. Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. Proceedings of the national academy of sciences, 108(supplement_1), 4516-4522.

8. Kuczynski, J., Lauber, C. L., Walters, W. A., Parfrey, L. W., Clemente, J. C., Gevers, D., & Knight, R. (2012). Experimental and analytical tools for studying the human microbiome. Nature Reviews Genetics, 13(1), 47-58.

9. Morgan, X. C., & Huttenhower, C. (2012). Chapter 12: Human microbiome analysis. PLoS computational biology, 8(12), e1002808.

10. Caporaso, J. G., Lauber, C. L., Costello, E. K., Berg-Lyons, D., Gonzalez, A., Stombaugh, J., ... & Knight, R. (2011). Moving pictures of the human microbiome. Genome biology, 12(5), 1-8.

11. McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217

12. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686.

13. Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2022). _vegan: Community Ecology Package_. R package version 2.6-4, <https://CRAN.R-project.org/package=vegan>.

14. Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

15. Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. https://doi.org/10.1093/bioinformatics/btw313.
