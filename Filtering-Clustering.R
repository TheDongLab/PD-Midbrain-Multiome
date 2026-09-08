# =============================================================================
#  QC filtering, normalization, Harmony integration, WNN clustering
# -----------------------------------------------------------------------------
# Purpose : Apply cell/sample QC filters, add donor metadata, run RNA (SCT) and
#           ATAC (LSI/TF-IDF) processing with Harmony batch integration,
#           weighted-nearest-neighbor (WNN) joint clustering, doublet removal
#           (scDblFinder), TSS-enrichment refinement, manual cell-type
#           annotation, and finally add the CD8 T-cell subcluster.
# Inputs  : merged Seurat object from step 01 (here: ./01-Clusterings4/merged_seurat.rda)
#           donor metadata spreadsheet (METADATA_XLSX)
# Outputs : <DIR_FINALWNN>/merged_seurat4.rda,
#           merged_seurat4+Tcells.rda and a series of QC / UMAP / marker PDFs.
# Run     : Rscript Filtering-Clustering.R
# Notes   : Package-install and non-English comments were removed/translated;
#           analysis steps and output files are preserved.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(ggplot2)
  library(hdf5r)
  library(openxlsx)
  library(harmony)
  library(cluster)
  library(scDblFinder)
  library(BiocParallel)
  library(mbkmeans)
  library(RColorBrewer)
})

# Donor metadata spreadsheet (edit path for your environment).
METADATA_XLSX <- "/metadata/2023-07-26-midbrain-scATAC-metadata.xlsx"

# Input merged object from step 01.
load("./merged_seurat.rda")

## ---- QC visualization -------------------------------------------------------

p1 <- DensityScatter(merged_seurat, x = "nCount_ATAC", y = "TSS.enrichment",
                     log_x = TRUE, quantiles = c(2, 5, 10, 90, 95, 98))
pdf(file.path(DIR_FINALWNN, "QC_DensityScatter.pdf"), width = 10, height = 10); print(p1); dev.off()

Idents(merged_seurat) <- "All"
p2 <- VlnPlot(merged_seurat,
              features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
                           "nucleosome_signal", "percent.mt"),
              ncol = 5, pt.size = 0, combine = TRUE)
pdf(file.path(DIR_FINALWNN, "QC_VlnPlot.pdf"), width = 15, height = 6); print(p2); dev.off()

p3 <- DensityScatter(merged_seurat, x = "nCount_ATAC", y = "nCount_RNA",
                     log_x = TRUE, quantiles = c(2, 5, 10, 90, 95, 98))
pdf(file.path(DIR_FINALWNN, "QC_DensityScatter_ATACvsRNA.pdf"), width = 10, height = 10); print(p3); dev.off()

p4 <- DensityScatter(merged_seurat, x = "nCount_RNA", y = "percent.mt",
                     log_x = TRUE, quantiles = c(2, 5, 10, 90, 95, 98))
pdf(file.path(DIR_FINALWNN, "QC_DensityScatter_percent.mtvsRNA.pdf"), width = 10, height = 10); print(p4); dev.off()

## ---- Cell- and sample-level QC filtering ------------------------------------

merged_seurat2 <- subset(x = merged_seurat, subset =
                           nCount_RNA > 1e2 & nCount_RNA < 1e5 &
                           nCount_ATAC > 5e1 & nCount_ATAC < 1e5 &
                           percent.mt < 15 &
                           nucleosome_signal < 5 &
                           TSS.enrichment > 1)


## ---- Add donor metadata -----------------------------------------------------
metadata <- read.xlsx(METADATA_XLSX)
merged_seurat2 <- AddMetaData(
  object = merged_seurat2,
  metadata = metadata[match(merged_seurat2$sample_id, metadata$Sample), ]
)
# Keep samples with a valid RIN value.
sample2remain <- merged_seurat2@meta.data$sample_id[which(merged_seurat2@meta.data$RIN != "N/A")]
merged_seurat2 <- subset(merged_seurat2, subset = sample_id %in% sample2remain)

# Derive age / PMI / RIN brackets.
age_bracket <- as.vector(cut(merged_seurat2@meta.data$Age, c(50, 60, 70, 80, 90, 100, 110)))
age_bracket <- gsub("\\(|]", "", age_bracket)
age_bracket <- gsub(",", "-", age_bracket)
age_bracket <- gsub("90-100", "90-110", gsub("100-110", "90-110", age_bracket))
merged_seurat2@meta.data$age_bracket <- age_bracket
merged_seurat2@meta.data$PMI_bracket <- droplevels(cut(as.numeric(merged_seurat2@meta.data$PMI), seq(0, 7, by = 1)))
merged_seurat2@meta.data$RIN_bracket <- droplevels(cut(as.numeric(merged_seurat2@meta.data$RIN), seq(0, 10, by = 1)))

RIN_bracket <- merged_seurat2@meta.data$RIN_bracket
RIN_bracket <- gsub(",", "-", RIN_bracket); RIN_bracket <- gsub("\\(", "", RIN_bracket); RIN_bracket <- gsub("]", "", RIN_bracket)
merged_seurat2@meta.data$RIN_bracket <- as.factor(RIN_bracket)

PMI_bracket <- merged_seurat2@meta.data$PMI_bracket
PMI_bracket <- gsub(",", "-", PMI_bracket); PMI_bracket <- gsub("\\(", "", PMI_bracket); PMI_bracket <- gsub("]", "", PMI_bracket)
merged_seurat2@meta.data$PMI_bracket <- as.factor(PMI_bracket)

## ---- Rebuild RNA assay by merging layers; drop low-signal features -----------

rna_layers <- names(merged_seurat2[["RNA"]]@layers)
matrix_list <- list()
for (layer in rna_layers) {
  matrix_list[[layer]] <- GetAssayData(object = merged_seurat2[["RNA"]], slot = layer)
}
merged_matrix <- do.call(cbind, matrix_list)
merged_matrix <- Matrix::Matrix(merged_matrix, sparse = TRUE)
merged_seurat2[["RNA"]] <- CreateAssayObject(counts = merged_matrix)

# Keep features detected in >= 10 cells (RNA and ATAC).
tmp <- Matrix::rowSums(merged_seurat2[["RNA"]]@counts > 0)
merged_seurat2[["RNA"]] <- subset(merged_seurat2[["RNA"]], features = names(which(tmp >= 10)))
tmp <- Matrix::rowSums(merged_seurat2[["ATAC"]]@counts > 0)
merged_seurat2[["ATAC"]] <- subset(merged_seurat2[["ATAC"]], features = names(which(tmp >= 10)))

save(merged_seurat2, file = file.path(DIR_FINALWNN, "merged_seurat2.rda"))

## ---- Normalization, dimensional reduction, Harmony, WNN ---------------------

load(file.path(DIR_FINALWNN, "merged_seurat2.rda"))
n.pc <- 30; n.lsi <- 10
vars_to_correct1 <- c("BatchLib", "Sample", "Sex", "age_bracket", "PMI_bracket", "RIN_bracket")
vars_to_correct2 <- c("BatchLib", "Sample", "Sex", "age_bracket", "PMI_bracket")

options(future.globals.maxSize = 20 * 1024^3)

# RNA: SCT -> PCA -> Harmony -> UMAP / clustering.
DefaultAssay(merged_seurat2) <- "RNA"
merged_seurat2 <- SCTransform(merged_seurat2)
merged_seurat2 <- RunPCA(merged_seurat2)
merged_seurat2 <- RunUMAP(merged_seurat2, dims = 1:n.pc, reduction = "pca",
                          reduction.name = "umap.rna", assay = "SCT", reduction.key = "rnaUMAP_")
merged_seurat2 <- RunHarmony(merged_seurat2, group.by.vars = vars_to_correct1,
                             assay.use = "SCT", reduction = "pca", reduction.save = "pca.harmony")
merged_seurat2 <- RunUMAP(merged_seurat2, dims = 1:n.pc, reduction = "pca.harmony", reduction.name = "pca.harmony.umap")
merged_seurat2 <- FindNeighbors(merged_seurat2, reduction = "pca.harmony", dims = 1:n.pc, assay = "SCT")
merged_seurat2 <- FindClusters(merged_seurat2, graph.name = "SCT_snn", algorithm = 3, resolution = 0.2)

# ATAC: TF-IDF -> SVD (LSI) -> Harmony -> UMAP / clustering.
DefaultAssay(merged_seurat2) <- "ATAC"
merged_seurat2 <- RunTFIDF(merged_seurat2)
merged_seurat2 <- FindTopFeatures(merged_seurat2, min.cutoff = "q0")
merged_seurat2 <- RunSVD(merged_seurat2)
merged_seurat2 <- RunUMAP(merged_seurat2, reduction = "lsi", dims = 2:n.lsi,
                          reduction.name = "umap.atac", reduction.key = "atacUMAP_")
merged_seurat2 <- RunHarmony(object = merged_seurat2, group.by.vars = vars_to_correct2,
                             assay.use = "ATAC", reduction = "lsi", reduction.save = "lsi.harmony", project.dim = FALSE)
merged_seurat2 <- RunUMAP(merged_seurat2, dims = 2:n.lsi, reduction = "lsi.harmony", reduction.name = "lsi.harmony.umap")
merged_seurat2 <- FindNeighbors(merged_seurat2, reduction = "lsi.harmony", dims = 2:n.lsi, assay = "ATAC")
merged_seurat2 <- FindClusters(merged_seurat2, graph.name = "ATAC_snn", algorithm = 3, resolution = 0.2)

# WNN joint graph.
merged_seurat2 <- FindMultiModalNeighbors(merged_seurat2, prune.SNN = 1/20,
                                          reduction.list = list("pca.harmony", "lsi.harmony"),
                                          dims.list = list(1:n.pc, 2:n.lsi),
                                          modality.weight.name = c("RNA.weight", "ATAC.weight"))
merged_seurat2 <- FindClusters(merged_seurat2, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.2)
merged_seurat2 <- RunUMAP(merged_seurat2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
## ---- Choose clustering resolution by silhouette score -----------------------

load(file.path(DIR_FINALWNN, "merged_seurat2.rda"))

calculate_silhouette_downsample <- function(seurat_obj, num_subsamples = 5, subsample_size = 10000) {
  silhouette_scores <- numeric(num_subsamples)
  for (i in 1:num_subsamples) {
    set.seed(42 + i)
    subset_cells <- sample(Cells(seurat_obj), size = subsample_size)
    subset_seurat <- subset(seurat_obj, cells = subset_cells)
    clusters <- Idents(subset_seurat)
    print(table(clusters))
    dist_matrix <- dist(Embeddings(subset_seurat, "wnn.umap"))
    sil <- silhouette(as.numeric(clusters), dist_matrix)
    silhouette_scores[i] <- mean(sil[, 3])
    print(silhouette_scores[i])
  }
  mean(silhouette_scores)
}

resolutions <- seq(0.1, 1.0, by = 0.1)
avg_silhouette_scores <- numeric(length(resolutions))
pdf(file.path(DIR_FINALWNN, "UMAP_plots_all_Resolutions2.pdf"))
for (i in seq_along(resolutions)[5:10]) {
  resolution <- resolutions[i]
  merged_seurat2 <- FindClusters(merged_seurat2, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
  merged_seurat2 <- RunUMAP(merged_seurat2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  group <- paste0("wsnn_res.", resolution)
  avg_silhouette_scores[i] <- calculate_silhouette_downsample(merged_seurat2, num_subsamples = 2, subsample_size = 10000)
  umap_plot <- DimPlot(merged_seurat2, reduction = "wnn.umap", group.by = group, label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle(paste("Resolution:", resolution))
  print(umap_plot)
}
dev.off()

Idents(merged_seurat2) <- merged_seurat2$wsnn_res.0.2
pdf(file.path(DIR_FINALWNN, "resolutions_vs_silhouetteScores.pdf"))
plot(resolutions, avg_silhouette_scores, type = "b", xlab = "Resolution",
     ylab = "Average Silhouette Score", main = "Silhouette Scores for Different Resolutions")
dev.off()
best_resolution <- resolutions[which.max(avg_silhouette_scores)]
cat("Best resolution based on silhouette score:", best_resolution, "\n")

# Final resolution = 0.2.
merged_seurat2 <- FindClusters(merged_seurat2, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)
merged_seurat2 <- RunUMAP(merged_seurat2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

## ---- Doublet detection (scDblFinder) on RNA and ATAC ------------------------

load(file.path(DIR_FINALWNN, "merged_seurat2.rda"))
DefaultAssay(merged_seurat2) <- "RNA"
seurat_object <- merged_seurat2
sce <- as.SingleCellExperiment(seurat_object)
sce <- scDblFinder(sce, samples = "sample_id", BPPARAM = MulticoreParam(3))
table(sce$scDblFinder.class)
merged_seurat2$scDblFinder_class_RNA <- sce$scDblFinder.class

cluster_counts <- table(seurat_object$wsnn_res.0.2)
doublet_counts <- table(seurat_object$wsnn_res.0.2, merged_seurat2$scDblFinder_class_RNA)
doublet_proportions <- doublet_counts[, "doublet"] / cluster_counts
doublet_prop_df <- data.frame(Cluster = names(doublet_proportions),
                              Doublet_Proportion = as.numeric(doublet_proportions))
write.table(doublet_prop_df, file = file.path(DIR_FINALWNN, "doublet_prop_df.RNA.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t", col.names = TRUE)

DefaultAssay(merged_seurat2) <- "ATAC"
seurat_object2 <- merged_seurat2
sce2 <- as.SingleCellExperiment(seurat_object2)
sce2 <- scDblFinder(sce2, artificialDoublets = 1, aggregateFeatures = TRUE, nfeatures = 25,
                    processing = "normFeatures", samples = "sample_id")
table(sce2$scDblFinder.class)
merged_seurat2$scDblFinder_class_ATAC <- sce2$scDblFinder.class

cluster_counts <- table(seurat_object2$wsnn_res.0.2)
doublet_counts <- table(seurat_object2$wsnn_res.0.2, merged_seurat2$scDblFinder_class_ATAC)
doublet_proportions <- doublet_counts[, "doublet"] / cluster_counts
doublet_prop_df_ATAC <- data.frame(Cluster = names(doublet_proportions),
                                   Doublet_Proportion = as.numeric(doublet_proportions))
write.table(doublet_prop_df_ATAC, file = file.path(DIR_FINALWNN, "doublet_prop_df.ATAC.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t", col.names = TRUE)

# Mark a cell a doublet if flagged in either modality.
merged_seurat2$doublets <- "singlet"
merged_seurat2$doublets[which(merged_seurat2$scDblFinder_class_ATAC == "doublet" |
                                merged_seurat2$scDblFinder_class_RNA == "doublet")] <- "doublet"

## ---- Remove doublets + cluster 14, then re-cluster --------------------------

merged_seurat3 <- subset(merged_seurat2, idents = 14, invert = TRUE)
merged_seurat3 <- subset(merged_seurat3, subset = scDblFinder_class_RNA == "singlet" & scDblFinder_class_ATAC == "singlet")

n.pc <- 30; n.lsi <- 10
DefaultAssay(merged_seurat3) <- "RNA"
merged_seurat3 <- FindNeighbors(merged_seurat3, reduction = "pca.harmony", dims = 1:n.pc, assay = "SCT", graph.name = "SCT_snn_new")
merged_seurat3 <- FindClusters(merged_seurat3, graph.name = "SCT_snn_new", algorithm = 3, resolution = 0.2)
merged_seurat3 <- RunUMAP(merged_seurat3, dims = 1:n.pc, reduction = "pca.harmony", reduction.name = "pca.harmony.umap_new")

DefaultAssay(merged_seurat3) <- "ATAC"
merged_seurat3 <- FindNeighbors(merged_seurat3, reduction = "lsi.harmony", dims = 2:n.lsi, assay = "ATAC", graph.name = "ATAC_snn_new")
merged_seurat3 <- FindClusters(merged_seurat3, graph.name = "ATAC_snn_new", algorithm = 3, resolution = 0.2)
merged_seurat3 <- RunUMAP(merged_seurat3, dims = 2:n.lsi, reduction = "lsi.harmony", reduction.name = "lsi.harmony.umap_new")

merged_seurat3 <- FindMultiModalNeighbors(merged_seurat3, prune.SNN = 1/20,
                                          reduction.list = list("pca.harmony", "lsi.harmony"),
                                          dims.list = list(1:n.pc, 2:n.lsi),
                                          modality.weight.name = c("RNA.weight", "ATAC.weight"),
                                          knn.graph.name = "wknn_new", snn.graph.name = "wsnn_new")
merged_seurat3 <- FindClusters(merged_seurat3, graph.name = "wsnn_new", algorithm = 3, verbose = FALSE, resolution = 0.2)
merged_seurat3 <- RunUMAP(merged_seurat3, nn.name = "weighted.nn", reduction.name = "wnn.umap_new", reduction.key = "new.wnnUMAP_")

## ---- Cluster markers and manual annotation ----------------------------------
Idents(merged_seurat3) <- merged_seurat3$wsnn_new_res.0.2
markers <- FindAllMarkers(merged_seurat3, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

wb <- createWorkbook()
for (cluster in unique(top_markers$cluster)) {
  cluster_markers <- top_markers %>% filter(cluster == !!cluster)
  addWorksheet(wb, sheetName = paste0("Cluster_", cluster))
  writeData(wb, sheet = paste0("Cluster_", cluster), x = cluster_markers)
}
saveWorkbook(wb, file = file.path(DIR_FINALWNN, "top_markers_by_cluster_res0.2_new.xlsx"), overwrite = TRUE)

## ---- Final re-clustering (res = 0.2) ----------------------------

DefaultAssay(merged_seurat4) <- "RNA"
merged_seurat4 <- FindClusters(merged_seurat4, graph.name = "SCT_snn_new", algorithm = 3, resolution = 0.2)
merged_seurat4 <- RunUMAP(merged_seurat4, dims = 1:n.pc, reduction = "pca.harmony", reduction.name = "pca.harmony.umap_new")

DefaultAssay(merged_seurat4) <- "ATAC"
merged_seurat4 <- FindClusters(merged_seurat4, graph.name = "ATAC_snn_new", algorithm = 3, resolution = 0.2)
merged_seurat4 <- RunUMAP(merged_seurat4, dims = 2:n.lsi, reduction = "lsi.harmony", reduction.name = "lsi.harmony.umap_new")

merged_seurat4 <- FindMultiModalNeighbors(merged_seurat4, prune.SNN = 1/20,
                                          reduction.list = list("pca.harmony", "lsi.harmony"),
                                          dims.list = list(1:n.pc, 2:n.lsi),
                                          modality.weight.name = c("RNA.weight", "ATAC.weight"),
                                          knn.graph.name = "wknn_new", snn.graph.name = "wsnn_new")
merged_seurat4 <- FindClusters(merged_seurat4, graph.name = "wsnn_new", algorithm = 3, verbose = FALSE, resolution = 0.2)
merged_seurat4 <- RunUMAP(merged_seurat4, nn.name = "weighted.nn", reduction.name = "wnn.umap_new", reduction.key = "newwnnUMAP_")

p1 <- DimPlot(merged_seurat4, reduction = "pca.harmony.umap_new", group.by = "SCT_snn_new_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(merged_seurat4, reduction = "lsi.harmony.umap_new", group.by = "ATAC_snn_new_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(merged_seurat4, reduction = "wnn.umap_new", group.by = "wsnn_new_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.Final_after_res0.222.pdf"), width = 18, height = 5)
print(p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))); dev.off()

p4 <- FeaturePlot(merged_seurat4, features = "RNA.weight", reduction = "wnn.umap_new", cols = c("blue", "red")) + ggtitle("RNA.weight") & theme(plot.title = element_text(size = 10))
p5 <- FeaturePlot(merged_seurat4, features = "ATAC.weight", reduction = "wnn.umap_new", cols = c("blue", "red")) + ggtitle("ATAC.weight") & theme(plot.title = element_text(size = 10))
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.RNAWeights_Final_after_res0.2.pdf"), width = 6, height = 5); print(p4); dev.off()
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.ATACWeights_Final_after_res0.2.pdf"), width = 6, height = 5); print(p5); dev.off()

DefaultAssay(merged_seurat4) <- "RNA"
features <- paste0("sct_", markerGenes)
Idents(merged_seurat4) <- merged_seurat4$wsnn_new_res.0.2
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.Markers_WNN_Final_after_res0.2.pdf"), width = 25, height = 80)
print(FeaturePlot(merged_seurat4, features = features, ncol = 3, reduction = "wnn.umap_new")); dev.off()

markers <- FindAllMarkers(merged_seurat4, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
h <- DoHeatmap(merged_seurat4, features = top_markers$gene, assay = "SCT")
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.top_markers_Final_after_res0.2.pdf"), width = 30, height = 15); print(h); dev.off()

wb <- createWorkbook()
for (cluster in unique(top_markers$cluster)) {
  cluster_markers <- top_markers %>% filter(cluster == !!cluster)
  addWorksheet(wb, sheetName = paste0("Cluster_", cluster))
  writeData(wb, sheet = paste0("Cluster_", cluster), x = cluster_markers)
}
saveWorkbook(wb, file = file.path(DIR_FINALWNN, "top_markers_by_cluster_Final_after_res0.2.xlsx"), overwrite = TRUE)

# WNN annotation (final-after).
Idents(merged_seurat4) <- merged_seurat4$wsnn_new_res.0.2
merged_seurat4 <- RenameIdents(merged_seurat4, "4" = "Astrocytes")
merged_seurat4 <- RenameIdents(merged_seurat4, "0" = "Oligodendrocytes", "1" = "Oligodendrocytes", "2" = "Oligodendrocytes", "3" = "Oligodendrocytes", "5" = "Oligodendrocytes", "8" = "Oligodendrocytes", "11" = "Oligodendrocytes", "14" = "Oligodendrocytes", "12" = "Oligodendrocytes", "13" = "Oligodendrocytes")
merged_seurat4 <- RenameIdents(merged_seurat4, "7" = "OPCs")
merged_seurat4 <- RenameIdents(merged_seurat4, "10" = "Endothelial cell")
merged_seurat4 <- RenameIdents(merged_seurat4, "6" = "Microglia")
merged_seurat4 <- RenameIdents(merged_seurat4, "9" = "Neurons")
merged_seurat4$celltype.WNN.after <- Idents(merged_seurat4)
=
## ---- Sub-clustering to resolve neuron / immune / vascular subtypes ----------

Idents(merged_seurat4) <- "wsnn_new_res.0.2"
merged_seurat4 <- FindSubCluster(merged_seurat4, cluster = 9, graph.name = "wsnn_new", algorithm = 3, res = 0.05)
Idents(merged_seurat4) <- "sub.cluster"
merged_seurat4 <- FindSubCluster(merged_seurat4, cluster = 10, graph.name = "wsnn_new", algorithm = 3, res = 0.05)
Idents(merged_seurat4) <- "sub.cluster"
merged_seurat4 <- FindSubCluster(merged_seurat4, cluster = 6, graph.name = "wsnn_new", algorithm = 3, res = 0.05)
Idents(merged_seurat4) <- "sub.cluster"

p3 <- DimPlot(merged_seurat4, reduction = "wnn.umap_new", group.by = "sub.cluster", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pdf(file.path(DIR_FINALWNN, "UMAP.Clustering.Final_after_res0.2-subclusters.pdf"), width = 7, height = 6); print(p3 & theme(plot.title = element_text(hjust = 0.5))); dev.off()

Idents(merged_seurat4) <- merged_seurat4$sub.cluster
merged_seurat4 <- RenameIdents(merged_seurat4, "4" = "Astrocytes")
merged_seurat4 <- RenameIdents(merged_seurat4, "0" = "Oligodendrocytes", "1" = "Oligodendrocytes", "2" = "Oligodendrocytes", "3" = "Oligodendrocytes", "5" = "Oligodendrocytes", "8" = "Oligodendrocytes", "11" = "Oligodendrocytes", "14" = "Oligodendrocytes", "12" = "Oligodendrocytes", "13" = "Oligodendrocytes")
merged_seurat4 <- RenameIdents(merged_seurat4, "7" = "OPCs")
merged_seurat4 <- RenameIdents(merged_seurat4, "10_0" = "Endothelial cell", "10_1" = "Pericytes")
merged_seurat4 <- RenameIdents(merged_seurat4, "6_0" = "Microglia", "6_1" = "macrophages","6_2" = "T cells")
merged_seurat4 <- RenameIdents(merged_seurat4, "9_1" = "GABA neurons", "9_0" = "Glu neurons", "9_2" = "DA neurons")
merged_seurat4$celltype.WNN.after.sub <- Idents(merged_seurat4)
save(merged_seurat4, file = file.path(DIR_FINALWNN, "merged_seurat4.rda"))

## ---- Cell composition by diagnosis ------------------------------------------
metadata <- merged_seurat4@meta.data
cell_composition <- metadata %>%
  group_by(Diagnosis, celltype.WNN.after.combine.update2) %>%
  summarise(count = n()) %>% ungroup() %>%
  group_by(Diagnosis) %>% mutate(proportion = count / sum(count)) %>% ungroup()

pdf(file.path(DIR_FINALWNN, "celltype_vs_diagnosis+Tcells.pdf"), width = 8, height = 6)
print(ggplot(cell_composition, aes(x = Diagnosis, y = proportion, fill = celltype.WNN.after.combine.update2)) +
        geom_bar(stat = "identity", position = "fill") +
        labs(title = "Cell Composition by Diagnosis", x = "Diagnosis", y = "Proportion", fill = "Cell Type") +
        theme_minimal() + scale_fill_brewer(palette = "Set3") +
        scale_y_continuous(labels = scales::percent) + coord_flip())
dev.off()

## ---- Covariate UMAPs --------------------------------------------------------

for (cov in c("age_bracket", "Sex", "Diagnosis", "BatchLib")) {
  show_legend <- !(cov == "BatchLib")
  p <- DimPlot(merged_seurat4, reduction = "wnn.umap_new", group.by = cov, label = (cov != "BatchLib"),
               label.size = 0, repel = TRUE, raster = FALSE, pt.size = 0.05) + ggtitle(cov)
  base <- file.path(DIR_FINALWNN, paste0("UMAP.Clustering.Manual.WNN_Final_", cov))
  pdf(paste0(base, ".pdf"), width = 6, height = 5)
  print(p & theme(plot.title = element_text(hjust = 0.5),
                  legend.position = if (show_legend) "right" else "none")); dev.off()
  ggsave(filename = paste0(base, ".png"), plot = p & theme(plot.title = element_text(hjust = 0.5)),
         width = 6, height = 5, dpi = 300)
}

## ---- ATAC weight distribution across cell types -----------------------------

DefaultAssay(merged_seurat4) <- "RNA"
Idents(merged_seurat4) <- "celltype.WNN.after.combine.update"
meta <- merged_seurat4@meta.data
p <- ggplot(meta, aes(x = celltype.WNN.after.combine.update, y = ATAC.weight,
                      fill = celltype.WNN.after.combine.update)) +
  geom_violin(scale = "width") +
  stat_summary(fun = median, geom = "point", size = 1, color = "black") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.minor = element_blank(), legend.position = "none") +
  labs(title = "Distribution of ATAC Weight Across Cell Types", x = "Cell Type", y = "ATAC Weight")
ggsave(file.path(DIR_DIFF, "dars/nebula_sc_renew/ATAC_weight_by_celltype_violin.pdf"),
       plot = p, width = 9, height = 6, dpi = 300)
