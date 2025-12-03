# functions.R
# Read 10X file and create Seurat object

library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

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
