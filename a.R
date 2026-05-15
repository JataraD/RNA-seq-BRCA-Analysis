# Crear carpetas en la ubicación correcta
dir.create("~/Juan/RNA-seq-BRCA-Analysis/results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("~/Juan/RNA-seq-BRCA-Analysis/results/tables",  recursive = TRUE, showWarnings = FALSE)

# Copiar desde Documents hacia Juan
file.copy(
  from = "C:/Users/lucho/OneDrive/Documents/RNA-seq-BRCA-Analysis/results/figures/01_PCA_calidad.png",
  to   = "~/Juan/RNA-seq-BRCA-Analysis/results/figures/01_PCA_calidad.png",
  overwrite = TRUE
)
file.copy(
  from = "C:/Users/lucho/OneDrive/Documents/RNA-seq-BRCA-Analysis/results/figures/02_Volcano_DEG.png",
  to   = "~/Juan/RNA-seq-BRCA-Analysis/results/figures/02_Volcano_DEG.png",
  overwrite = TRUE
)
file.copy(
  from = "C:/Users/lucho/OneDrive/Documents/RNA-seq-BRCA-Analysis/results/figures/03_Heatmap_top50.png",
  to   = "~/Juan/RNA-seq-BRCA-Analysis/results/figures/03_Heatmap_top50.png",
  overwrite = TRUE
)
file.copy(
  from = "C:/Users/lucho/OneDrive/Documents/RNA-seq-BRCA-Analysis/results/tables/DEG_anotado.csv",
  to   = "~/Juan/RNA-seq-BRCA-Analysis/results/tables/DEG_anotado.csv",
  overwrite = TRUE
)

# Verificar
list.files("~/Juan/RNA-seq-BRCA-Analysis/results/figures/")

knitr::include_graphics("results/figures/01_PCA_calidad.png")
knitr::include_graphics("results/figures/02_Volcano_DEG.png")
knitr::include_graphics("results/figures/03_Heatmap_top50.png")