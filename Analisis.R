# ============================================================
# 03_DESeq2_analysis.R
# Análisis de Expresión Diferencial — TCGA-BRCA
# Tumor vs. Tejido Normal | DESeq2
# Autor: Juan David Atará Delgado
# Fecha: Mayo 2026
# ============================================================

library(DESeq2)
library(tidyverse)

# Cargar objeto guardado
dds <- readRDS("RNA-seq-BRCA-Analysis/data/processed/dds_brca.rds")

# Establecer referencia (normal = control)
dds$condition <- relevel(dds$condition, ref = "normal")

# Correr DESeq2 — este paso tarda 2–5 minutos
cat("Corriendo DESeq2...\n")
dds <- DESeq(dds)
cat("✅ DESeq2 completado\n")

# Guardar
saveRDS(dds, "RNA-seq-BRCA-Analysis/data/processed/dds_brca_fitted.rds")

#2

# Resultados: tumor vs normal
res <- results(dds, 
               contrast  = c("condition", "tumor", "normal"),
               alpha     = 0.05)

# Resumen
summary(res)

# Convertir a dataframe limpio
res_df <- as.data.frame(res) |>
  rownames_to_column("gene_id") |>
  arrange(padj) |>
  filter(!is.na(padj))

cat("\nTotal genes analizados:", nrow(res_df), "\n")
cat("Genes significativos (padj < 0.05):", 
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated en tumor:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated en tumor:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("Downregulated en tumor:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE), "\n")

#3

# Guardar tabla completa
write.csv(res_df, 
          "RNA-seq-BRCA-Analysis/results/tables/DEG_tumor_vs_normal.csv",
          row.names = FALSE)

# Top 20 genes más significativos
top20 <- res_df |>
  filter(padj < 0.05, abs(log2FoldChange) > 1) |>
  slice_head(n = 20)

print(top20[, c("gene_id", "log2FoldChange", "padj")])

cat(" Resultados guardados en results/tables/\n")


#4

library(org.Hs.eg.db)
library(AnnotationDbi)

# Limpiar IDs (quitar versión: .10, .12, etc.)
res_df$gene_id_clean <- gsub("\\..*", "", res_df$gene_id)

# Mapear Ensembl ID → símbolo de gen
res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys      = res_df$gene_id_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Ver top 20 con nombres — usando subset() en lugar de select()
top20_anotado <- res_df |>
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1) |>
  subset(select = c(gene_symbol, gene_id, log2FoldChange, padj)) |>
  head(20)

print(top20_anotado)

# Guardar tabla anotada
write.csv(res_df,
          "RNA-seq-BRCA-Analysis/results/tables/DEG_anotado.csv",
          row.names = FALSE)

cat(" Tabla anotada guardada\n")


#5

library(EnhancedVolcano)

# Resolver duplicados: quedarse con el gen más significativo de cada símbolo
res_volcano <- res_df |>
  dplyr::filter(!is.na(gene_symbol)) |>
  dplyr::arrange(padj) |>
  dplyr::distinct(gene_symbol, .keep_all = TRUE) |>  # elimina duplicados
  tibble::column_to_rownames("gene_symbol")

# Genes a etiquetar
genes_label <- rownames(res_volcano)[
  res_volcano$padj < 0.05 & 
    abs(res_volcano$log2FoldChange) > 4
][1:15]

# Guardar
png("RNA-seq-BRCA-Analysis/results/figures/02_Volcano_DEG.png",
    width = 1000, height = 800)

EnhancedVolcano(res_volcano,
                lab            = rownames(res_volcano),
                x              = "log2FoldChange",
                y              = "padj",
                title          = "Expresión Diferencial — TCGA-BRCA",
                subtitle       = "Tumor vs. Tejido Normal | DESeq2 | padj < 0.05 | |LFC| > 1",
                pCutoff        = 0.05,
                FCcutoff       = 1,
                pointSize      = 2.5,
                labSize        = 4,
                selectLab      = genes_label,
                drawConnectors = TRUE,
                col            = c("grey70", "#3498DB", "#2ECC71", "#E74C3C"),
                legendLabels   = c("NS", "LFC", "p-adj", "p-adj & LFC"),
                legendPosition = "right",
                caption        = paste0("Total genes: ", nrow(res_volcano),
                                        " | Significativos: ",
                                        sum(res_volcano$padj < 0.05, na.rm = TRUE))
)

dev.off()
cat(" Volcano plot guardado\n")

#6

library(pheatmap)
library(DESeq2)
library(SummarizedExperiment)

# Cargar objeto fitted
dds_fitted <- readRDS("RNA-seq-BRCA-Analysis/data/processed/dds_brca_fitted.rds")

# VST normalización
vsd <- DESeq2::vst(dds_fitted, blind = FALSE)

# Top 50 genes más significativos con símbolo conocido
top50_ids <- res_df |>
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1, !is.na(gene_symbol)) |>
  dplyr::arrange(padj) |>
  dplyr::slice_head(n = 50) |>
  dplyr::pull(gene_id)

# Limpiar IDs para hacer match
top50_ids_clean <- gsub("\\..*", "", top50_ids)

# Matriz de expresión normalizada
mat <- assay(vsd)
rownames(mat) <- gsub("\\..*", "", rownames(mat))
mat <- mat[top50_ids_clean, ]

# Reemplazar IDs por símbolos
simbolos <- res_df$gene_symbol[match(top50_ids_clean, res_df$gene_id_clean)]
rownames(mat) <- ifelse(is.na(simbolos), top50_ids_clean, simbolos)

# Anotación columnas
anno_col <- data.frame(
  Condicion = ifelse(vsd$condition == "tumor", "Tumor", "Normal"),
  row.names = colnames(mat)
)

anno_colors <- list(
  Condicion = c(Tumor = "#E74C3C", Normal = "#2ECC71")
)

# Guardar
png("RNA-seq-BRCA-Analysis/results/figures/03_Heatmap_top50.png",
    width = 1000, height = 1200)

pheatmap(mat,
         annotation_col   = anno_col,
         annotation_colors = anno_colors,
         scale            = "row",
         show_colnames    = FALSE,
         fontsize_row     = 8,
         color            = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
         main             = "Top 50 DEGs — TCGA-BRCA\nTumor vs. Normal | VST normalizado",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean"
)

dev.off()
cat(" Heatmap guardado\n")

