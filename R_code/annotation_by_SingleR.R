# annotation_by_SingleR.R
# Automatic cell type annotation using SingleR and celldex references


library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)


annotate_cells <- function(sr_obj, ref_data = NULL) {
# Use HumanPrimaryCellAtlasData if no reference provided
if (is.null(ref_data)) {
ref_data <- HumanPrimaryCellAtlasData()
}


# Convert Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(sr_obj)


# Run SingleR (assign main labels)
pred.main <- SingleR(test = sce, ref = ref_data, labels = ref_data$label.main)


# Attach SingleR labels to Seurat metadata
sr_obj$SingleR.labels <- pred.main$labels


return(sr_obj)
}
