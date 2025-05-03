# ------------------------------------------
# Assignment 3: Single-Cell RNA-Seq Analysis
# ------------------------------------------

# Install required packages (only if not already installed)
packages <- c("Seurat", "patchwork", "dplyr", "ggplot2", "hdf5r", "caret", "e1071", "randomForest")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(caret)
library(e1071)
library(randomForest)

# --------------------
# Load & Preprocess Data
# --------------------

# Path to HDF5 file
h5_file <- "C:/Users/tooba/Downloads/pbmc_1k_v2_filtered_feature_bc_matrix.h5"

# Load data
sc_data <- Read10X_h5(h5_file)
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "PBMC1k", min.cells = 3, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC plots
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize & find variable features
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# ------------------------
# Dimensionality Reduction & Clustering
# ------------------------
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# PCA plots with legend and axes
VizDimLoadings(seurat_obj, dims = 1, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca", label = TRUE) + theme_minimal()

ElbowPlot(seurat_obj)

set.seed(42)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# UMAP with legend and axes
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  theme_minimal() + 
  ggtitle("UMAP Clustering") + 
  theme(legend.position = "right")

# ------------------------
# Differential Expression
# ------------------------
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# Visualize marker genes
VlnPlot(seurat_obj, features = c("CD3D", "MS4A1"), pt.size = 0.1) + 
  theme_minimal()

DoHeatmap(seurat_obj, features = unique(top_markers$gene[1:30])) + 
  theme_minimal() +
  labs(title = "Top Marker Gene Expression")

# ------------------------
# Machine Learning Classification
# ------------------------
pca_features <- Embeddings(seurat_obj, "pca")[, 1:10]
cell_labels <- as.factor(Idents(seurat_obj))
ml_data <- data.frame(pca_features)
ml_data$Cluster <- cell_labels

set.seed(123)
train_idx <- createDataPartition(ml_data$Cluster, p = 0.7, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data <- ml_data[-train_idx, ]

# Train kNN classifier
knn_model <- train(Cluster ~ ., data = train_data, method = "knn", tuneLength = 5)
predictions <- predict(knn_model, newdata = test_data)
cm <- confusionMatrix(predictions, test_data$Cluster)
print(cm)
print(cm$overall['Accuracy'])
