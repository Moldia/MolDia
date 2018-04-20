---
title: "<center>__MolDia__</center><center>Single cell In-Situ RCA data analysis</center>"
author: <center>*Mats Nilsson*</center><center>Department of Biochemistry and Biophysics</center><center>Stockholm
  University, Sweden</center>
date: '<center> Date: *`r Sys.Date()`*</center><center>__*Version: 1.0.0*__</center>'
output:
  BiocStyle::html_document:
    number_sections: yes
    toc: yes
    toc_depth: 6
  BiocStyle::pdf_document:
    number_sections: yes
    toc: yes
    toc_depth: 6
bibliography: bibliography.bib
vignette: |
  %\VignetteIndexEntry{<center>__MolDia__</center><center>Single cell In-Situ RCA data analysis</center>}   %\VignetteEngine{knitr::rmarkdown}   %\VignetteEncoding{UTF-8}   %\usepackage[utf8x]{inputenc}
---

# MolDia
R package to analyze single cell in-situ sequencing (ISS) data

# Install package 
```{r}
install.packages("devtools")
library(devtools)
install_github("mashranga/MolDia")
```
# Scope of the package
The scope of the package can be devided into following broad category
## Read and manipulate RCA ISS data
1. __readRCA__ : Function to read RCA ISS data. THis function has different parameter, how to read and what to read. See detail in function help.
2. __ISS_rotate__ : Rotate a tissue by specific angle or flip  by x or y axis.

## ISS data summary and vizualization
1. __genesSummary__: Summary of Specific gene of interest in ISS data. This include % of gene of interest positive cell with their distribution in other positive cells and visa-versa.
2. __ISS_pieplot__:
3. __ISS_ratiocor__:

## Multiple ISS data comparision
1. __ISS_compare__: Relation of multiple ISS data interms of total reads per gene. The result is presented in fitted regression line with R^2.



## Need to updrade 

1. In function RCA_map plot empty cells foe specific genes.
2. legend problem In function RCA_map plot
3. Add violon plot in function RCA_map
4. Add tsne subplot for selected gene 