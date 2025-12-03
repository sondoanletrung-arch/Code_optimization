# Automatic cell type annotation by SingleR
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

annotate_cells <- function(sr_obj, ref_data = NULL) {
  #install HumanPrimaryCellAtlasData for annotation
  if (is.null(ref_data)) {
    ref_data <- HumanPrimaryCellAtlasData()
  }
  
  #Seurat object -> SingleCellExperiment object
  sce <- as.SingleCellExperiment(sr_obj)
  
  #SingleR
  pred.main <- SingleR(test = sce, ref = ref_data, labels = ref_data$label.main)
  
  #Add SingleR labels to Seurat
  sr_obj$SingleR.labels <- pred.main$labels
  
  return(sr_obj)
}
