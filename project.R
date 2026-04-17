options(repos = c(CRAN = "https://cloud.r-project.org"))

# install if needed
# install.packages(c(
#   "BiocManager", "dplyr", "Seurat", "Matrix", "uwot", "harmony",
#   "ggplot2", "mclust", "aricode", "cluster", "FNN", "igraph", "NMF"
# ))
# BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(dplyr)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(uwot)
library(harmony)
library(ggplot2)
library(mclust)
library(aricode)
library(cluster)
library(FNN)
library(igraph)
library(NMF)


set.seed(100)

### counts matrix (already normalized) ###
download.file(
  "https://datasets.cellxgene.cziscience.com/d0919452-3528-46e9-bc53-130f81ab0126.h5ad",
  "spectrum_tcells.h5ad",
  method = "libcurl"
)

sce <- readH5AD("spectrum_tcells.h5ad")

mat <- assay(sce, "X")
if (!inherits(mat, "dgCMatrix")) {
  mat <- as(mat, "dgCMatrix")
}

meta <- as.data.frame(colData(sce))

# inspecting 
dim(mat)
colnames(meta)
unique(meta$cell_type)

#### seurat object ####

seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta)
LayerData(seurat_obj, assay = "RNA", layer = "data") <- mat

#### pca, with subset of 30000 cells ####

cells_use <- sample(colnames(seurat_obj), 30000)
seurat_obj_small <- subset(seurat_obj, cells = cells_use)

seurat_obj_small <- FindVariableFeatures(
  seurat_obj_small,
  selection.method = "vst",
  nfeatures = 2000
)

seurat_obj_small <- ScaleData(
  seurat_obj_small,
  features = VariableFeatures(seurat_obj_small),
  verbose = FALSE
)

seurat_obj_small <- RunPCA(
  seurat_obj_small,
  features = VariableFeatures(seurat_obj_small),
  npcs = 50,
  verbose = FALSE
)

Reductions(seurat_obj_small)

### louvain with no harmony ###

louvain_noharm <- seurat_obj_small

louvain_noharm <- FindNeighbors(
  louvain_noharm,
  reduction = "pca",
  dims = 1:50,
  verbose = FALSE
)

louvain_noharm <- FindClusters(
  louvain_noharm,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 1,
  verbose = FALSE
)

louvain_noharm <- RunUMAP(
  louvain_noharm,
  reduction = "pca",
  dims = 1:50,
  verbose = FALSE
)

louvain_noharm$RNA_snn_res.0.1 <- factor(
  louvain_noharm$RNA_snn_res.0.1,
  levels = sort(as.numeric(levels(louvain_noharm$RNA_snn_res.0.1)))
)

louvain_noharm$RNA_snn_res.0.2 <- factor(
  louvain_noharm$RNA_snn_res.0.2,
  levels = sort(as.numeric(levels(louvain_noharm$RNA_snn_res.0.2)))
)

louvain_noharm$RNA_snn_res.0.3 <- factor(
  louvain_noharm$RNA_snn_res.0.3,
  levels = sort(as.numeric(levels(louvain_noharm$RNA_snn_res.0.3)))
)

p_louvain_noharm_01 <- DimPlot(
  louvain_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain without Harmony, resolution 0.1")

p_louvain_noharm_02 <- DimPlot(
  louvain_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain without Harmony, resolution 0.2")

p_louvain_noharm_03 <- DimPlot(
  louvain_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain without Harmony, resolution 0.3")

print(p_louvain_noharm_01)
print(p_louvain_noharm_02)
print(p_louvain_noharm_03)

ggsave("louvain_no_harmony_res_0.1.png", plot = p_louvain_noharm_01, width = 6, height = 5, dpi = 300)
ggsave("louvain_no_harmony_res_0.2.png", plot = p_louvain_noharm_02, width = 6, height = 5, dpi = 300)
ggsave("louvain_no_harmony_res_0.3.png", plot = p_louvain_noharm_03, width = 6, height = 5, dpi = 300)


### louvain with harmony ###

louvain_harm <- seurat_obj_small

louvain_harm <- RunHarmony(
  object = louvain_harm,
  group.by.vars = "patient_id",
  reduction.use = "pca",
  dims.use = 1:50,
  verbose = FALSE
)

louvain_harm <- FindNeighbors(
  louvain_harm,
  reduction = "harmony",
  dims = 1:50,
  verbose = FALSE
)

louvain_harm <- FindClusters(
  louvain_harm,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 1,
  verbose = FALSE
)

louvain_harm <- RunUMAP(
  louvain_harm,
  reduction = "harmony",
  dims = 1:50,
  verbose = FALSE
)

louvain_harm$RNA_snn_res.0.1 <- factor(
  louvain_harm$RNA_snn_res.0.1,
  levels = sort(as.numeric(levels(louvain_harm$RNA_snn_res.0.1)))
)

louvain_harm$RNA_snn_res.0.2 <- factor(
  louvain_harm$RNA_snn_res.0.2,
  levels = sort(as.numeric(levels(louvain_harm$RNA_snn_res.0.2)))
)

louvain_harm$RNA_snn_res.0.3 <- factor(
  louvain_harm$RNA_snn_res.0.3,
  levels = sort(as.numeric(levels(louvain_harm$RNA_snn_res.0.3)))
)

p_louvain_harm_01 <- DimPlot(
  louvain_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain with Harmony, resolution 0.1")

p_louvain_harm_02 <- DimPlot(
  louvain_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain with Harmony, resolution 0.2")

p_louvain_harm_03 <- DimPlot(
  louvain_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("Louvain with Harmony, resolution 0.3")

print(p_louvain_harm_01)
print(p_louvain_harm_02)
print(p_louvain_harm_03)

ggsave("louvain_harmony_res_0.1.png", plot = p_louvain_harm_01, width = 6, height = 5, dpi = 300)
ggsave("louvain_harmony_res_0.2.png", plot = p_louvain_harm_02, width = 6, height = 5, dpi = 300)
ggsave("louvain_harmony_res_0.3.png", plot = p_louvain_harm_03, width = 6, height = 5, dpi = 300)


### leiden no harmony ###

leiden_noharm <- seurat_obj_small

leiden_noharm <- FindNeighbors(
  leiden_noharm,
  reduction = "pca",
  dims = 1:50,
  verbose = FALSE
)

leiden_noharm <- FindClusters(
  leiden_noharm,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 4,
  verbose = FALSE
)

leiden_noharm <- RunUMAP(
  leiden_noharm,
  reduction = "pca",
  dims = 1:50,
  verbose = FALSE
)

leiden_noharm$RNA_snn_res.0.1 <- factor(
  leiden_noharm$RNA_snn_res.0.1,
  levels = sort(as.numeric(levels(leiden_noharm$RNA_snn_res.0.1)))
)

leiden_noharm$RNA_snn_res.0.2 <- factor(
  leiden_noharm$RNA_snn_res.0.2,
  levels = sort(as.numeric(levels(leiden_noharm$RNA_snn_res.0.2)))
)

leiden_noharm$RNA_snn_res.0.3 <- factor(
  leiden_noharm$RNA_snn_res.0.3,
  levels = sort(as.numeric(levels(leiden_noharm$RNA_snn_res.0.3)))
)

p_leiden_noharm_01 <- DimPlot(
  leiden_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden without Harmony, resolution 0.1")

p_leiden_noharm_02 <- DimPlot(
  leiden_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden without Harmony, resolution 0.2")

p_leiden_noharm_03 <- DimPlot(
  leiden_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden without Harmony, resolution 0.3")

print(p_leiden_noharm_01)
print(p_leiden_noharm_02)
print(p_leiden_noharm_03)

ggsave("leiden_no_harmony_res_0.1.png", plot = p_leiden_noharm_01, width = 7, height = 6, dpi = 300)
ggsave("leiden_no_harmony_res_0.2.png", plot = p_leiden_noharm_02, width = 7, height = 6, dpi = 300)
ggsave("leiden_no_harmony_res_0.3.png", plot = p_leiden_noharm_03, width = 7, height = 6, dpi = 300)

### leiden with harmony ###

leiden_harm <- seurat_obj_small

leiden_harm <- RunHarmony(
  object = leiden_harm,
  group.by.vars = "patient_id",
  reduction.use = "pca",
  dims.use = 1:50,
  verbose = FALSE
)

leiden_harm <- FindNeighbors(
  leiden_harm,
  reduction = "harmony",
  dims = 1:50,
  verbose = FALSE
)

leiden_harm <- FindClusters(
  leiden_harm,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 4,
  verbose = FALSE
)

leiden_harm <- RunUMAP(
  leiden_harm,
  reduction = "harmony",
  dims = 1:50,
  verbose = FALSE
)

leiden_harm$RNA_snn_res.0.1 <- factor(
  leiden_harm$RNA_snn_res.0.1,
  levels = sort(as.numeric(levels(leiden_harm$RNA_snn_res.0.1)))
)

leiden_harm$RNA_snn_res.0.2 <- factor(
  leiden_harm$RNA_snn_res.0.2,
  levels = sort(as.numeric(levels(leiden_harm$RNA_snn_res.0.2)))
)

leiden_harm$RNA_snn_res.0.3 <- factor(
  leiden_harm$RNA_snn_res.0.3,
  levels = sort(as.numeric(levels(leiden_harm$RNA_snn_res.0.3)))
)

p_leiden_harm_01 <- DimPlot(
  leiden_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden with Harmony, resolution 0.1")

p_leiden_harm_02 <- DimPlot(
  leiden_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden with Harmony, resolution 0.2")

p_leiden_harm_03 <- DimPlot(
  leiden_harm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("Leiden with Harmony, resolution 0.3")


print(p_leiden_harm_01)
print(p_leiden_harm_02)
print(p_leiden_harm_03)

ggsave("leiden_harmony_res_0.1.png", plot = p_leiden_harm_01, width = 7, height = 6, dpi = 300)
ggsave("leiden_harmony_res_0.2.png", plot = p_leiden_harm_02, width = 7, height = 6, dpi = 300)
ggsave("leiden_harmony_res_0.3.png", plot = p_leiden_harm_03, width = 7, height = 6, dpi = 300)

######## metrics###############

#### silhouette scores ####

sil_cells <- sample(colnames(seurat_obj_small), 5000)
# pca embeddings for no harmony
emb_noharm <- Embeddings(seurat_obj_small, "pca")[sil_cells, 1:20]
# pca embeddings for harmony
emb_harm <- Embeddings(louvain_harm, "harmony")[sil_cells, 1:20]

dist_noharm <- dist(emb_noharm)
dist_harm <- dist(emb_harm)

### louvain no harmony ###
sil_louvain_noharm_01 <- silhouette(
  as.numeric(louvain_noharm@meta.data[sil_cells, "RNA_snn_res.0.1"]),
  dist_noharm
)
sil_louvain_noharm_02 <- silhouette(
  as.numeric(louvain_noharm@meta.data[sil_cells, "RNA_snn_res.0.2"]),
  dist_noharm
)
sil_louvain_noharm_03 <- silhouette(
  as.numeric(louvain_noharm@meta.data[sil_cells, "RNA_snn_res.0.3"]),
  dist_noharm
)

#print(sil_louvain_noharm_01)

### louvain harmony ###
sil_louvain_harm_01 <- silhouette(
  as.numeric(louvain_harm@meta.data[sil_cells, "RNA_snn_res.0.1"]),
  dist_harm
)
sil_louvain_harm_02 <- silhouette(
  as.numeric(louvain_harm@meta.data[sil_cells, "RNA_snn_res.0.2"]),
  dist_harm
)
sil_louvain_harm_03 <- silhouette(
  as.numeric(louvain_harm@meta.data[sil_cells, "RNA_snn_res.0.3"]),
  dist_harm
)

### leiden no harmony ###
sil_leiden_noharm_01 <- silhouette(
  as.numeric(leiden_noharm@meta.data[sil_cells, "RNA_snn_res.0.1"]),
  dist_noharm
)
sil_leiden_noharm_02 <- silhouette(
  as.numeric(leiden_noharm@meta.data[sil_cells, "RNA_snn_res.0.2"]),
  dist_noharm
)
sil_leiden_noharm_03 <- silhouette(
  as.numeric(leiden_noharm@meta.data[sil_cells, "RNA_snn_res.0.3"]),
  dist_noharm
)

### leiden harmony ###
sil_leiden_harm_01 <- silhouette(
  as.numeric(leiden_harm@meta.data[sil_cells, "RNA_snn_res.0.1"]),
  dist_harm
)
sil_leiden_harm_02 <- silhouette(
  as.numeric(leiden_harm@meta.data[sil_cells, "RNA_snn_res.0.2"]),
  dist_harm
)
sil_leiden_harm_03 <- silhouette(
  as.numeric(leiden_harm@meta.data[sil_cells, "RNA_snn_res.0.3"]),
  dist_harm
)

### silhouette results table ###

silhouette_results <- data.frame(
  method = c(
    "Louvain_noHarmony_0.1", "Louvain_noHarmony_0.2", "Louvain_noHarmony_0.3",
    "Louvain_Harmony_0.1",   "Louvain_Harmony_0.2",   "Louvain_Harmony_0.3",
    "Leiden_noHarmony_0.1",  "Leiden_noHarmony_0.2",  "Leiden_noHarmony_0.3",
    "Leiden_Harmony_0.1",    "Leiden_Harmony_0.2",    "Leiden_Harmony_0.3"
  ),
  silhouette = c(
    mean(sil_louvain_noharm_01[, "sil_width"]),
    mean(sil_louvain_noharm_02[, "sil_width"]),
    mean(sil_louvain_noharm_03[, "sil_width"]),
    mean(sil_louvain_harm_01[, "sil_width"]),
    mean(sil_louvain_harm_02[, "sil_width"]),
    mean(sil_louvain_harm_03[, "sil_width"]),
    mean(sil_leiden_noharm_01[, "sil_width"]),
    mean(sil_leiden_noharm_02[, "sil_width"]),
    mean(sil_leiden_noharm_03[, "sil_width"]),
    mean(sil_leiden_harm_01[, "sil_width"]),
    mean(sil_leiden_harm_02[, "sil_width"]),
    mean(sil_leiden_harm_03[, "sil_width"])
  )
)

print(silhouette_results)
write.csv(silhouette_results, "silhouette_results.csv", row.names = FALSE)

#### ARI and NMI ####

ari_nmi_results <- data.frame(
  comparison = c(
    "Leiden_Harmony_vs_Leiden_noHarmony_0.1",
    "Leiden_Harmony_vs_Leiden_noHarmony_0.2",
    "Leiden_Harmony_vs_Leiden_noHarmony_0.3",
    
    "Leiden_Harmony_vs_Louvain_Harmony_0.1",
    "Leiden_Harmony_vs_Louvain_Harmony_0.2",
    "Leiden_Harmony_vs_Louvain_Harmony_0.3",
    
    "Leiden_noHarmony_vs_Louvain_noHarmony_0.1",
    "Leiden_noHarmony_vs_Louvain_noHarmony_0.2",
    "Leiden_noHarmony_vs_Louvain_noHarmony_0.3",
    
    "Louvain_Harmony_vs_Louvain_noHarmony_0.1",
    "Louvain_Harmony_vs_Louvain_noHarmony_0.2",
    "Louvain_Harmony_vs_Louvain_noHarmony_0.3"
  ),
  ARI = c(
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.1, leiden_noharm$RNA_snn_res.0.1),
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.2, leiden_noharm$RNA_snn_res.0.2),
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.3, leiden_noharm$RNA_snn_res.0.3),
    
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.1, louvain_harm$RNA_snn_res.0.1),
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.2, louvain_harm$RNA_snn_res.0.2),
    adjustedRandIndex(leiden_harm$RNA_snn_res.0.3, louvain_harm$RNA_snn_res.0.3),
    
    adjustedRandIndex(leiden_noharm$RNA_snn_res.0.1, louvain_noharm$RNA_snn_res.0.1),
    adjustedRandIndex(leiden_noharm$RNA_snn_res.0.2, louvain_noharm$RNA_snn_res.0.2),
    adjustedRandIndex(leiden_noharm$RNA_snn_res.0.3, louvain_noharm$RNA_snn_res.0.3),
    
    adjustedRandIndex(louvain_harm$RNA_snn_res.0.1, louvain_noharm$RNA_snn_res.0.1),
    adjustedRandIndex(louvain_harm$RNA_snn_res.0.2, louvain_noharm$RNA_snn_res.0.2),
    adjustedRandIndex(louvain_harm$RNA_snn_res.0.3, louvain_noharm$RNA_snn_res.0.3)
  ),
  NMI = c(
    NMI(leiden_harm$RNA_snn_res.0.1, leiden_noharm$RNA_snn_res.0.1),
    NMI(leiden_harm$RNA_snn_res.0.2, leiden_noharm$RNA_snn_res.0.2),
    NMI(leiden_harm$RNA_snn_res.0.3, leiden_noharm$RNA_snn_res.0.3),
    
    NMI(leiden_harm$RNA_snn_res.0.1, louvain_harm$RNA_snn_res.0.1),
    NMI(leiden_harm$RNA_snn_res.0.2, louvain_harm$RNA_snn_res.0.2),
    NMI(leiden_harm$RNA_snn_res.0.3, louvain_harm$RNA_snn_res.0.3),
    
    NMI(leiden_noharm$RNA_snn_res.0.1, louvain_noharm$RNA_snn_res.0.1),
    NMI(leiden_noharm$RNA_snn_res.0.2, louvain_noharm$RNA_snn_res.0.2),
    NMI(leiden_noharm$RNA_snn_res.0.3, louvain_noharm$RNA_snn_res.0.3),
    
    NMI(louvain_harm$RNA_snn_res.0.1, louvain_noharm$RNA_snn_res.0.1),
    NMI(louvain_harm$RNA_snn_res.0.2, louvain_noharm$RNA_snn_res.0.2),
    NMI(louvain_harm$RNA_snn_res.0.3, louvain_noharm$RNA_snn_res.0.3)
  )
)

print(ari_nmi_results)
write.csv(ari_nmi_results, "ari_nmi_results.csv", row.names = FALSE)



##### NMF #####

# checking for nonnegative inputs
range(LayerData(seurat_obj_small, assay = "RNA", layer = "data"))

nmf_mat <- LayerData(seurat_obj_small, assay = "RNA", layer = "data")
var_genes <- VariableFeatures(seurat_obj_small)
nmf_mat <- nmf_mat[var_genes, ]
nmf_cells <- sample(colnames(nmf_mat), 5000)
nmf_genes <- sample(rownames(nmf_mat), 1000)

nmf_mat <- nmf_mat[nmf_genes, nmf_cells]
nmf_mat <- as.matrix(nmf_mat)

# remove rows that are all zero
nmf_mat <- nmf_mat[complete.cases(nmf_mat), ]
nmf_mat <- nmf_mat[rowSums(nmf_mat) > 0, ]
nmf_mat <- nmf_mat[apply(nmf_mat, 1, var) > 0, ]


dim(nmf_mat)
range(nmf_mat)
sum(is.na(nmf_mat))

length(unique(louvain_noharm$RNA_snn_res.0.3))
length(unique(leiden_noharm$RNA_snn_res.0.3))

nmf_rank <- 9

nmf_fit <- nmf(
  nmf_mat,
  rank = nmf_rank,
  method = "brunet",
  nrun = 2,
  seed = 100
)

nmf_embed <- t(coef(nmf_fit))   # cells x factors
rownames(nmf_embed) <- nmf_cells
colnames(nmf_embed) <- paste0("NMF_", 1:ncol(nmf_embed))

nmf_obj <- subset(seurat_obj_small, cells = nmf_cells)

nmf_reduction <- CreateDimReducObject(
  embeddings = nmf_embed,
  key = "NMF_",
  assay = "RNA"
)

nmf_obj[["nmf"]] <- nmf_reduction


### clustering in NMF space with louvain ###
nmf_obj <- FindNeighbors(
  nmf_obj,
  reduction = "nmf",
  dims = 1:nmf_rank,
  verbose = FALSE
)

nmf_obj <- FindClusters(
  nmf_obj,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 1,   # Louvain
  verbose = FALSE
)

nmf_obj <- RunUMAP(
  nmf_obj,
  reduction = "nmf",
  dims = 1:nmf_rank,
  verbose = FALSE
)

p_nmf_louvain_01 <- DimPlot(
  nmf_obj,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Louvain, resolution 0.1")

p_nmf_louvain_02 <- DimPlot(
  nmf_obj,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Louvain, resolution 0.2")

p_nmf_louvain_03 <- DimPlot(
  nmf_obj,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Louvain, resolution 0.3")

print(p_nmf_louvain_01)
print(p_nmf_louvain_02)
print(p_nmf_louvain_03)

ggsave("nmf_louvain_res_0.1.png", plot = p_nmf_louvain_01, width = 6, height = 5, dpi = 300)
ggsave("nmf_louvain_res_0.2.png", plot = p_nmf_louvain_02, width = 6, height = 5, dpi = 300)
ggsave("nmf_louvain_res_0.3.png", plot = p_nmf_louvain_03, width = 6, height = 5, dpi = 300)


### clustering in NMF space with leiden ###

nmf_obj_leiden <- nmf_obj

nmf_obj_leiden <- FindNeighbors(
  nmf_obj_leiden,
  reduction = "nmf",
  dims = 1:nmf_rank,
  verbose = FALSE
)

nmf_obj_leiden <- FindClusters(
  nmf_obj_leiden,
  resolution = c(0.1, 0.2, 0.3),
  algorithm = 4,
  verbose = FALSE
)

p_nmf_leiden_01 <- DimPlot(
  nmf_obj_leiden,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Leiden, resolution 0.1")

p_nmf_leiden_02 <- DimPlot(
  nmf_obj_leiden,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Leiden, resolution 0.2")

p_nmf_leiden_03 <- DimPlot(
  nmf_obj_leiden,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF Leiden, resolution 0.3")

nmf_obj_leiden <- RunUMAP(
  nmf_obj_leiden,
  reduction = "nmf",
  dims = 1:nmf_rank,
  verbose = FALSE
)

print(p_nmf_leiden_01)
print(p_nmf_leiden_02)
print(p_nmf_leiden_03)

ggsave("nmf_leiden_res_0.1.png", plot = p_nmf_leiden_01, width = 6, height = 5, dpi = 300)
ggsave("nmf_leiden_res_0.2.png", plot = p_nmf_leiden_02, width = 6, height = 5, dpi = 300)
ggsave("nmf_leiden_res_0.3.png", plot = p_nmf_leiden_03, width = 6, height = 5, dpi = 300)

### nmf metrics ###
nmf_embed_all <- Embeddings(nmf_obj, "nmf")

# subset cells for silhouette to reduce computation
sil_cells_nmf <- sample(rownames(nmf_embed_all), min(3000, nrow(nmf_embed_all)))

emb_nmf_sub <- nmf_embed_all[sil_cells_nmf, 1:nmf_rank]
dist_nmf <- dist(emb_nmf_sub)

### Louvain on NMF ###
sil_nmf_louvain_01 <- silhouette(
  as.numeric(nmf_obj@meta.data[sil_cells_nmf, "RNA_snn_res.0.1"]),
  dist_nmf
)

sil_nmf_louvain_02 <- silhouette(
  as.numeric(nmf_obj@meta.data[sil_cells_nmf, "RNA_snn_res.0.2"]),
  dist_nmf
)

sil_nmf_louvain_03 <- silhouette(
  as.numeric(nmf_obj@meta.data[sil_cells_nmf, "RNA_snn_res.0.3"]),
  dist_nmf
)

### Leiden on NMF ###
sil_nmf_leiden_01 <- silhouette(
  as.numeric(nmf_obj_leiden@meta.data[sil_cells_nmf, "RNA_snn_res.0.1"]),
  dist_nmf
)

sil_nmf_leiden_02 <- silhouette(
  as.numeric(nmf_obj_leiden@meta.data[sil_cells_nmf, "RNA_snn_res.0.2"]),
  dist_nmf
)

sil_nmf_leiden_03 <- silhouette(
  as.numeric(nmf_obj_leiden@meta.data[sil_cells_nmf, "RNA_snn_res.0.3"]),
  dist_nmf
)

nmf_silhouette_results <- data.frame(
  method = c(
    "NMF_Louvain_0.1",
    "NMF_Louvain_0.2",
    "NMF_Louvain_0.3",
    "NMF_Leiden_0.1",
    "NMF_Leiden_0.2",
    "NMF_Leiden_0.3"
  ),
  silhouette = c(
    mean(sil_nmf_louvain_01[, "sil_width"]),
    mean(sil_nmf_louvain_02[, "sil_width"]),
    mean(sil_nmf_louvain_03[, "sil_width"]),
    mean(sil_nmf_leiden_01[, "sil_width"]),
    mean(sil_nmf_leiden_02[, "sil_width"]),
    mean(sil_nmf_leiden_03[, "sil_width"])
  )
)

print(nmf_silhouette_results)
write.csv(nmf_silhouette_results, "nmf_silhouette_results.csv", row.names = FALSE)


### NMF ARI and NMI ###

nmf_ari_nmi_results <- data.frame(
  comparison = c(
    "NMF_Leiden_vs_NMF_Louvain_0.1",
    "NMF_Leiden_vs_NMF_Louvain_0.2",
    "NMF_Leiden_vs_NMF_Louvain_0.3"
  ),
  ARI = c(
    adjustedRandIndex(nmf_obj_leiden$RNA_snn_res.0.1, nmf_obj$RNA_snn_res.0.1),
    adjustedRandIndex(nmf_obj_leiden$RNA_snn_res.0.2, nmf_obj$RNA_snn_res.0.2),
    adjustedRandIndex(nmf_obj_leiden$RNA_snn_res.0.3, nmf_obj$RNA_snn_res.0.3)
  ),
  NMI = c(
    NMI(nmf_obj_leiden$RNA_snn_res.0.1, nmf_obj$RNA_snn_res.0.1),
    NMI(nmf_obj_leiden$RNA_snn_res.0.2, nmf_obj$RNA_snn_res.0.2),
    NMI(nmf_obj_leiden$RNA_snn_res.0.3, nmf_obj$RNA_snn_res.0.3)
  )
)

print(nmf_ari_nmi_results)
write.csv(nmf_ari_nmi_results, "nmf_ari_nmi_results.csv", row.names = FALSE)

#### nmf vs pca results ####

common_cells <- colnames(nmf_obj)


nmf_vs_pca_results <- data.frame(
  comparison = c(
    "NMF_Louvain_vs_PCA_Louvain_noHarmony_0.3",
    "NMF_Louvain_vs_PCA_Leiden_noHarmony_0.3",
    "NMF_Leiden_vs_PCA_Louvain_noHarmony_0.3",
    "NMF_Leiden_vs_PCA_Leiden_noHarmony_0.3"
  ),
  ARI = c(
    adjustedRandIndex(
      nmf_obj@meta.data[common_cells, "RNA_snn_res.0.3"],
      louvain_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    adjustedRandIndex(
      nmf_obj@meta.data[common_cells, "RNA_snn_res.0.3"],
      leiden_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    adjustedRandIndex(
      nmf_obj_leiden@meta.data[common_cells, "RNA_snn_res.0.3"],
      louvain_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    adjustedRandIndex(
      nmf_obj_leiden@meta.data[common_cells, "RNA_snn_res.0.3"],
      leiden_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    )
  ),
  NMI = c(
    NMI(
      nmf_obj@meta.data[common_cells, "RNA_snn_res.0.3"],
      louvain_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    NMI(
      nmf_obj@meta.data[common_cells, "RNA_snn_res.0.3"],
      leiden_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    NMI(
      nmf_obj_leiden@meta.data[common_cells, "RNA_snn_res.0.3"],
      louvain_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    ),
    NMI(
      nmf_obj_leiden@meta.data[common_cells, "RNA_snn_res.0.3"],
      leiden_noharm@meta.data[common_cells, "RNA_snn_res.0.3"]
    )
  )
)

print(nmf_vs_pca_results)
write.csv(nmf_vs_pca_results, "nmf_vs_pca_results.csv", row.names = FALSE)

### plotting pca and nmf side by side ###

### louvain ###
p_pca <- DimPlot(
  louvain_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("PCA + Louvain (0.3)")

p_nmf <- DimPlot(
  nmf_obj,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF + Louvain (0.3)")

p_nmf <- p_nmf +
  annotate("text", x = Inf, y = Inf,
           label = "Silhouette = 0.289",
           hjust = 1.1, vjust = 1.5,
           size = 4)

p_pca <- p_pca +
  annotate("text", x = Inf, y = Inf,
           label = "Silhouette = 0.092",
           hjust = 1.1, vjust = 1.5,
           size = 4)

p_compare <- p_pca | p_nmf
print(p_compare)
ggsave("pca_vs_nmf_louvain.png", p_compare, width = 10, height = 5, dpi = 300)

### leiden ###

p_pca_leiden <- DimPlot(
  leiden_noharm,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("PCA + Leiden (0.3)")


p_nmf_leiden <- DimPlot(
  nmf_obj_leiden,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  raster = TRUE
) + ggtitle("NMF + Leiden (0.3)")

# silhouette scores
p_nmf_leiden <- p_nmf_leiden +
  annotate("text", x = Inf, y = Inf,
           label = "Silhouette = 0.294",
           hjust = 1.1, vjust = 1.5,
           size = 4)

p_pca_leiden <- p_pca_leiden +
  annotate("text", x = Inf, y = Inf,
           label = "Silhouette = 0.080",
           hjust = 1.1, vjust = 1.5,
           size = 4)

p_compare_leiden <- p_pca_leiden | p_nmf_leiden

print(p_compare_leiden)

ggsave(
  "pca_vs_nmf_leiden.png",
  p_compare_leiden,
  width = 10,
  height = 5,
  dpi = 300
)