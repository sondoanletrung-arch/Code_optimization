# functions.R
# Read 10X file and create Seurat object
# Preprocess data: QC, normalize, dimension Reduction
library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

#MODULE 1
#Read folder
load_and_merge_data <- function(data_directory, project_name = "SeuratProject") {
  
  sample_folders <- list.dirs(data_directory, full.names = FALSE, recursive = FALSE)
  seurat_list <- list()

  #Read and create Seurat object by each 10X file in the folder
  for (folder in sample_folders) {
    data_counts <- Read10X(data.dir = file.path(data_directory, folder))
    obj <- CreateSeuratObject(counts = data_counts, project = folder, min.cells = 3, min.features = 200)
    seurat_list[[folder]] <- obj
  }

  if (length(seurat_list) > 1) {
    sr_merged <- merge(x = seurat_list[[1]], 
                       y = seurat_list[2:length(seurat_list)],
                       add.cell.ids = names(seurat_list), 
                       project = project_name)
  } else {
    sr_merged <- seurat_list[[1]]
  }

  sr_merged <- JoinLayers(sr_merged) #Join layer in Seuratv5
  return(sr_merged)
}

#MODULE 2
#QC, normalize, PCA, UMAP
process_seurat_object <- function(sr_obj, 
                                  min_nFeature = 200, 
                                  max_nFeature = 6000, 
                                  max_mt = 15, 
                                  n_variable_features = 2000,
                                  resolution = 0.5) {

  
  # Save raw count
  sr_obj[["RNA"]]$raw_count <- sr_obj[["RNA"]]$counts
  
  #MT percentage
  sr_obj[["percent.mt"]] <- PercentageFeatureSet(sr_obj, pattern = "^MT-")
  
  #Quality control
  sr_obj <- subset(sr_obj, subset = nFeature_RNA > min_nFeature & 
                                    nFeature_RNA < max_nFeature & 
                                    percent.mt < max_mt)

  #Normalization
  sr_obj <- NormalizeData(sr_obj)
  sr_obj <- FindVariableFeatures(sr_obj, selection.method = "vst", nfeatures = n_variable_features)
  
  #Scale data, PCA
  all.genes <- rownames(sr_obj)
  sr_obj <- ScaleData(sr_obj, features = all.genes)
  sr_obj <- RunPCA(sr_obj, features = VariableFeatures(object = sr_obj), verbose = FALSE)

  #Cluster, UMAP
  sr_obj <- FindNeighbors(sr_obj, dims = 1:20)
  sr_obj <- FindClusters(sr_obj, resolution = resolution)
  sr_obj <- RunUMAP(sr_obj, dims = 1:20)
  
  return(sr_obj)
}
