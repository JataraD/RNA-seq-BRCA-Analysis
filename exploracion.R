# ============================================================
# exploracion.R
# Exploración y control de calidad — TCGA-BRCA RNA-seq
# Autor: Juan David Atará Delgado
# Fecha: Mayo 2026
# ============================================================

library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)

# Cargar objeto guardado (no necesitas volver a descargar nunca)
se <- readRDS("RNA-seq-BRCA-Analysis/data/processed/se_brca_60samples.rds")

# Exploración básica
cat("=== ESTRUCTURA DEL DATASET ===\n")
cat("Genes totales:  ", nrow(se), "\n")
cat("Muestras totales:", ncol(se), "\n")

# Ver tipos de muestra
table(se$sample_type)

#2

# Extraer conteos crudos (unstranded es el estándar para DESeq2)
counts_matrix <- assay(se, "unstranded")

cat("Dimensiones matriz de conteos:", dim(counts_matrix), "\n")

# Ver primeras filas y columnas
counts_matrix[1:5, 1:3]

# Estadísticas básicas
cat("\nConteos cero (genes no expresados):", 
    sum(rowSums(counts_matrix) == 0), "\n")
cat("Conteos totales por muestra (primeras 5):\n")
print(colSums(counts_matrix)[1:5])

#3

# Filtro estándar: eliminar genes con < 10 conteos en total
keep <- rowSums(counts_matrix) >= 10
counts_filtrado <- counts_matrix[keep, ]

cat("Genes antes del filtro:", nrow(counts_matrix), "\n")
cat("Genes después del filtro:", nrow(counts_filtrado), "\n")
cat("Genes eliminados:", nrow(counts_matrix) - nrow(counts_filtrado), "\n")

#4

# Metadata de las muestras
metadata <- data.frame(
  sample     = colnames(counts_filtrado),
  condition  = ifelse(se$sample_type == "Primary Tumor", "tumor", "normal"),
  row.names  = colnames(counts_filtrado)
)

# Verificar balance
table(metadata$condition)

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtrado,
  colData   = metadata,
  design    = ~ condition
)

# Guardar
saveRDS(dds, "RNA-seq-BRCA-Analysis/data/processed/dds_brca.rds")
cat("✅ Objeto DESeq2 creado y guardado\n")

#5

# Normalización rápida para visualización
vsd <- vst(dds, blind = TRUE)

# PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Gráfico
png("RNA-seq-BRCA-Analysis/results/figures/01_PCA_calidad.png", 
    width = 800, height = 600)

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("tumor" = "#E74C3C", "normal" = "#2ECC71")) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  ggtitle("PCA — TCGA-BRCA: Tumor vs. Tejido Normal",
          subtitle = "RNA-seq | 60 muestras | STAR Counts") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

dev.off()
cat("✅ PCA guardado en results/figures/\n")
