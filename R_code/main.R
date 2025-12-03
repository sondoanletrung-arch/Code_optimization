#Running Seurat and inferCNV
#Load module
source("functions.R")
source("annotation_by_SingleR.R")
source("infercnv_utils.R")

#Load necessary libraries for execution, other libraries also load
library(infercnv)
library(ggplot2)

#input - output
data_dir <- "C:/Users/Administrator/Downloads/GSE282701_RAW-20251128T111238Z-1-001"
work_dir <- "C:/gse" 

if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
# 3. EXECUTE SINGLE-CELL ANALYSIS PIPELINE (MOD 1, 2, 3)
# ------------------------------------------------------------------------------

#MODULE 1: Load and Merge Data
sr <- load_and_merge_data(data_directory = data_dir, project_name = "GSE282701")

sr$TissueType <- ifelse(grepl("T$", sr$orig.ident), "Tumor", "Paracancerous") # Project-specific logic: Assign TissueType (Tumor/Paracancerous)
print(table(sr$orig.ident, sr$TissueType))

#MODULE 2: Preprocessing (QC, Normalize, PCA, UMAP)
# Using default parameters defined in functions.R (e.g., resolution=0.5, max_mt=15)
sr <- process_seurat_object(sr)

#MODULE 3: Cell Type Annotation
sr <- annotate_cells(sr)


#MODULE 4: INFERCNV ANALYSIS
gene_order_file <- file.path(work_dir, "gene_ordering_file.txt") 
out_dir_cnv <- file.path(work_dir, "infercnv_output_HMM") # InferCNV output folder

# 4.1. Prepare inputs (Annotation file and Raw Matrix extraction))
infercnv_data <- prepare_infercnv_inputs(sr, work_dir)

raw_matrix <- infercnv_data$raw_counts_matrix
annot_file <- infercnv_data$annotation_file # Path to the generated annotation file

# 4.2. Create InferCNV Object
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = raw_matrix,
    annotations_file = annot_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = c("T_cells") # Reference group for CNV calculation
)

# 4.3. Run InferCNV (with HMM)
message(">>> Running InferCNV analysis (HMM = TRUE). This may consume significant RAM.")
infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = out_dir_cnv,
    cluster_by_groups = T,
    denoise = T,
    HMM = T, # Run with Hidden Markov Model (As requested)
    num_threads = 4
)


#Save

# Save the final annotated Seurat object
saveRDS(sr, file = file.path(work_dir, "GSE282701_Annotated.rds"))

# Plot UMAP based on cell type annotation
DimPlot(sr, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
    ggtitle("Cell Type Annotation Results")

# Plot UMAP based on tissue type
DimPlot(sr, reduction = "umap", group.by = "TissueType", label = TRUE, repel = TRUE) + 
    ggtitle("Tissue Type (Tumor vs. Paracancerous)")
