# functions.R
# Read 10X files, create Seurat object, and perform preprocessing (QC, normalization, dimension reduction)


library(Seurat)
library(dplyr)


# MODULE 1: Load and merge 10X directories into a single Seurat object
load_and_merge_data <- function(data_directory, project_name = "SeuratProject") {
sample_folders <- list.dirs(data_directory, full.names = FALSE, recursive = FALSE)
seurat_list <- list()


# Read and create Seurat object for each 10X folder
for (folder in sample_folders) {
data_path <- file.path(data_directory, folder)
if (!dir.exists(data_path)) next
data_counts <- Read10X(data.dir = data_path)
obj <- CreateSeuratObject(counts = data_counts, project = folder, min.cells = 3, min.features = 200)
seurat_list[[folder]] <- obj
}


if (length(seurat_list) == 0) stop("No valid 10X folders found in the provided data_directory")


if (length(seurat_list) > 1) {
sr_merged <- merge(x = seurat_list[[1]],
y = seurat_list[2:length(seurat_list)],
add.cell.ids = names(seurat_list),
project = project_name)
} else {
sr_merged <- seurat_list[[1]]
}


# If using Seurat v5 and layers/assays are present, you can join layers
if (exists("JoinLayers")) {
try({
sr_merged <- JoinLayers(sr_merged)
}, silent = TRUE)
}


return(sr_merged)
}


# MODULE 2: QC, normalization, PCA, UMAP, clustering
process_seurat_object <- function(sr_obj,
min_nFeature = 200,
max_nFeature = 6000,
max_mt = 15,
n_variable_features = 2000,
resolution = 0.5) {


# Save raw counts into a layer to preserve original counts
if (!"raw_count" %in% Layers(sr_obj, assay = "RNA")) {
sr_obj[["RNA"]]$raw_count <- sr_obj[["RNA"]]$counts
}


# Compute mitochondrial percentage (human genes usually "MT-" prefix)
sr_obj[["percent.mt"]] <- PercentageFeatureSet(sr_obj, pattern = "^MT-")


# Quality control filtering
sr_obj <- subset(sr_obj, subset = nFeature_RNA > min_nFeature &
}
