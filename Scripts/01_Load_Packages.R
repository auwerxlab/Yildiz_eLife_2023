# Utilities
source("./Scripts/00_Helper_Functions.R")
# RSeQC retrieved from http://rseqc.sourceforge.net/#download-rseqc
RSeQC <- "/home/elalam/Project_RNAseq_CD_HFD_epcam/Tools/RSeQC-4.0.0"
# BedOps retrieved from https://github.com/bedops/bedops/releases/tag/v2.4.40
BedOps <- "/home/elalam/Project_RNAseq_CD_HFD_epcam/Tools/bedops_linux_x86_64-v2.4.40"

# Genomic Databases
library("org.Mm.eg.db")
library("AnnotationDbi")

# Differential Expression
library("edgeR")

# Enrichment Analysis
library("enrichplot")
library("clusterProfiler")

# Cell Type Deconvolution
library("ggpubr")
library("stringr")
library("GEOquery")
library("data.table")

# Implementation
library("parallel")

# Plotting & Downstream analysis
library("ggplot2")
library("ggrepel")
library("cowplot")
library("VennDiagram")
library("RColorBrewer")
library("FactoMineR")
library("ComplexHeatmap")

# Data Export
library("openxlsx")
library("readxl")
library("Rsamtools")
library("dplyr")
