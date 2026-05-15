library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# Query al repositorio TCGA-BRCA
query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

# Ver cuántas muestras hay por tipo
getResults(query)
  dplyr::count(sample_type)

  # Ver todos los barcodes disponibles por tipo
  resultados <- getResults(query)
  
  # Separar tumor y normal
  tumor   <- resultados |> dplyr::filter(sample_type == "Primary Tumor")
  normal  <- resultados |> dplyr::filter(sample_type == "Solid Tissue Normal")
  
  cat("Tumores disponibles:", nrow(tumor), "\n")
  cat("Normales disponibles:", nrow(normal), "\n")
  
  # Seleccionar subconjunto balanceado (40 tumor + 20 normal)
  set.seed(42)
  barcodes_seleccionados <- c(
    sample(tumor$cases, 40),
    sample(normal$cases, 20)
  )
  
  cat("Total muestras seleccionadas:", length(barcodes_seleccionados), "\n")
  