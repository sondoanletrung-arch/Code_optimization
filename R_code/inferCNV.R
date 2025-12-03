# inferCNV.R
# Prepare inputs for inferCNV and helper utilities


library(Seurat)
library(dplyr)


prepare_infercnv_inputs <- function(sr_obj, work_dir) {
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)


# Create annotation dataframe: expected inferCNV format is: cell_name<tab>cell_type (no header)
annot_df <- data.frame(
cell_id = colnames(sr_obj),
cell_type = sr_obj$SingleR.labels
)


annot_path <- file.path(work_dir, "infercnv_annotations.txt")


write.table(annot_df,
file = annot_path,
sep = "\t",
quote = FALSE,
row.names = FALSE,
col.names = FALSE)


# Extract raw counts matrix. Use saved layer if present, otherwise counts.
if ("raw_count" %in% Layers(sr_obj, assay = "RNA")) {
raw_matrix <- GetAssayData(sr_obj, assay = "RNA", layer = "raw_count")
} else {
raw_matrix <- GetAssayData(sr_obj, assay = "RNA", slot = "counts")
}


return(list(
annotation_file = annot_path,
raw_counts_matrix = raw_matrix
))
}
