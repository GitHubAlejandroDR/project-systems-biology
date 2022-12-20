#ESTABLECEMOS DIRECTORIO DONDE INSTALAMOS LAS LIBRERIAS
args=commandArgs(trailingOnly = TRUE)
lib_dir <- args[1]
.libPaths(lib_dir)

#LIBRERIAS
if (!require("BiocManager",quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(version="3.16")

if (!require("clusterProfiler",quietly=TRUE))
  BiocManager::install("clusterProfiler")

if (!require("org.Hs.eg.db",quietly=TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!require("DOSE",quietly=TRUE))
  BiocManager::install("DOSE")

if (!require("tidyverse",quietly=TRUE))
  BiocManager::install("tidyverse")

if (!require("igraph",quietly=TRUE))
  BiocManager::install("igraph")

if (!require("linkcomm",quietly=TRUE))
  BiocManager::install("linkcomm")

if (!require("dplyr",quietly=TRUE))
  BiocManager::install("dplyr")

if (!require("STRINGdb",quietly=TRUE))
  BiocManager::install("STRINGdb")

if (!require("DT",quietly=TRUE))
  BiocManager::install("DT")

if (!require("AnnotationDbi",quietly=TRUE))
  BiocManager::install("AnnotationDbi",force=T)








