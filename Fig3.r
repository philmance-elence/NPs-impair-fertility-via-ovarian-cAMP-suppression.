
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)


control.data <- Read10X(data.dir = "./data/Control/filtered_feature_bc_matrix/")
psnps.data <- Read10X(data.dir = "./data/PS_NPs/filtered_feature_bc_matrix/")


control <- CreateSeuratObject(counts = control.data, project = "Control", min.cells = 3, min.features = 200)
control$group <- "Control"
psnps <- CreateSeuratObject(counts = psnps.data, project = "PS-NPs", min.cells = 3, min.features = 200)
psnps$group <- "PS-NPs"
sc_obj <- merge(control, y = psnps, add.cell.ids = c("Ctrl", "PS"))

sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-") 
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

sc_obj <- ScaleData(sc_obj, features = rownames(sc_obj))
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))


sc_obj <- FindNeighbors(sc_obj, dims = 1:30)
sc_obj <- FindClusters(sc_obj, resolution = 0.4)
sc_obj <- RunUMAP(sc_obj, dims = 1:30)


markers_B <- c("Cd34", "Pecam1", "Cdh5", "Epcam", "Krt19", "Krt18", "Amh", "Fshr", 
               "Hsd17b1", "Ptprc", "Cd74", "Ctss", "Ptgfr", "Hsd17b7", "Ssu2", 
               "Col1a1", "Dcn", "Mgp", "Aldh1a1", "Ephx1", "Srd5a1")
level_order_B <- c("Theca", "Stroma", "Luteal", "Immune", "Granulosa", "Epithelium", "Endothelium")
sc_obj$cell_type <- factor(sc_obj$cell_type, levels = level_order_B)

DotPlot(sc_obj, features = markers_B, group.by = "cell_type") + 
  coord_flip() +
  scale_color_gradientn(colors = c("#e0e0e0", "#ff7f50", "#ff0000")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank())

DimPlot(sc_obj, reduction = "umap", group.by = "cell_type", split.by = "group", label = FALSE) + 
  theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_blank())

ratio_df <- sc_obj@meta.data %>%
  group_by(group, cell_type) %>%
  tally() %>%
  mutate(percent = n / sum(n) * 100)

ggplot(ratio_df, aes(x = group, y = percent, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5, position = "stack") +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Cells (%)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold"))


ggplot(ratio_df, aes(x = cell_type, y = percent, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.9, color = "white") +
  coord_polar(start = 0) + 
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = element_blank())


granulosa_obj <- subset(sc_obj, idents = "Granulosa")
granulosa_obj <- FindVariableFeatures(granulosa_obj, nfeatures = 1500)
granulosa_obj <- ScaleData(granulosa_obj)
granulosa_obj <- RunPCA(granulosa_obj)
granulosa_obj <- FindNeighbors(granulosa_obj, dims = 1:15)
granulosa_obj <- FindClusters(granulosa_obj, resolution = 0.3)
granulosa_obj <- RunUMAP(granulosa_obj, dims = 1:15)

markers_F <- c("Inhbb", "Fst", "Gja1", "Itih5", "Ghr", "Pik3ip1", "Igfbp5", 
               "Gatm", "Col18a1", "Top2a", "Cenpa", "Racgap1")

DotPlot(granulosa_obj, features = markers_F, group.by = "sub_cluster") + 
  scale_color_gradientn(colors = c("#f0f0f0", "#ff7f50")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DimPlot(granulosa_obj, reduction = "umap", group.by = "sub_cluster", split.by = "group") + 
  theme_bw() + 
  theme(panel.grid = element_blank())

ratio_sub_df <- granulosa_obj@meta.data %>%
  group_by(group, sub_cluster) %>%
  tally() %>%
  mutate(percent = n / sum(n) * 100)

ggplot(ratio_sub_df, aes(x = group, y = percent, fill = sub_cluster)) +
  geom_bar(stat = "identity", width = 0.5, position = "stack") +
  theme_classic()

library(monocle)

# 1. 表达矩阵与元数据准备
matrix_counts <- as(as.matrix(granulosa_obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = granulosa_obj@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(matrix_counts), 
                                                 row.names = row.names(matrix_counts)))

cds <- newCellDataSet(matrix_counts, phenoData = pd, featureData = fd, 
                       expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
ordering_genes <- VariableFeatures(granulosa_obj)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5) + 
  facet_wrap(~group) +
  scale_color_gradient(low = "#000033", high = "#66ccff") +
  theme_bw()
library(clusterProfiler)
library(org.Mm.eg.db)

Idents(granulosa_obj) <- "sub_cluster"
preantral_cells <- subset(granulosa_obj, idents = "Preantral")
Idents(preantral_cells) <- "group"
preantral_deg <- FindMarkers(preantral_cells, ident.1 = "PS-NPs", ident.2 = "Control", logfc.threshold = 0.1)

pre_genes <- rownames(preantral_deg[preantral_deg$p_val < 0.05, ])
pre_ids <- bitr(pre_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
kegg_j <- enrichKEGG(gene = pre_ids$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

antral_cells <- subset(granulosa_obj, idents = "Antral")
Idents(antral_cells) <- "group"
antral_deg <- FindMarkers(antral_cells, ident.1 = "PS-NPs", ident.2 = "Control", logfc.threshold = 0.1)

ant_genes <- rownames(antral_deg[antral_deg$p_val < 0.05, ])
ant_ids <- bitr(ant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
kegg_k <- enrichKEGG(gene = ant_ids$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

dotplot(kegg_k, x = "GeneRatio", showCategory = 10, color = "p.adjust") +
  scale_size_continuous(range = c(3, 8)) +
  labs(title = "Antral granulosa (PS-NPs vs Control)", x = "KEGG Enrichment") +
  theme_bw()
