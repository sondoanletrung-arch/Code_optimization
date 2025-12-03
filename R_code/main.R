# main.R
# Pipeline driver: load/merge, preprocess, annotate, and run inferCNV


# Source helper modules
source("functions.R")
source("annotation_by_SingleR.R")
source("inferCNV.R")


# Required libraries
library(infercnv)
library(ggplot2)


# Input / output paths (edit to your local environment)
data_dir <- "C:/Users/Administrator/Downloads/GSE282701_RAW"
work_dir <- "C:/gse"


if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)


# MODULE 1: Load and merge 10X datasets
sr <- load_and_merge_data(data_directory = data_dir, project_name = "GSE282701")


# Project-specific: determine tissue type if orig.ident encodes this information
if ("orig.ident" %in% colnames(sr@meta.data)) {
sr$TissueType <- ifelse(grepl("T$", sr$orig.ident), "Tumor", "Paracancerous")
print(table(sr$orig.ident, sr$TissueType))
}


# MODULE 2: Preprocessing (QC, normalization, PCA, UMAP)
sr <- process_seurat_object(sr)


# MODULE 3: Cell type annotation with SingleR
sr <- annotate_cells(sr)


# MODULE 4: Prepare inferCNV inputs
gene_order_file <- file.path(work_dir, "gene_ordering_file.txt")
out_dir_cnv <- file.path(work_dir, "infercnv_output_HMM")


infercnv_data <- prepare_infercnv_inputs(sr, work_dir)
raw_matrix <- infercnv_data$raw_counts_matrix
annot_file <- infercnv_data$annotation_file


# Create infercnv object
# NOTE: ensure gene_order_file exists and contains the required gene ordering information
infercnv_obj <- CreateInfercnvObject(
raw_counts_matrix = raw_matrix,
annotations_file = annot_file,
delim = "\t",
gene_order_file = gene_order_file,
ref_group_names = c("T_cells")
)


# Run inferCNV with reasonable defaults; adjust cutoff and threads as needed
infercnv_obj <- infercnv::run(
infercnv_obj,
cutoff = 0.1,
out_dir = out_dir_cnv,
cluster_by_groups = TRUE,
denoise = TRUE,
HMM = TRUE,
}, silent = TRUE)
