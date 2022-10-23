# reading required_packages
library(Seurat)
library(patchwork)
library(dplyr) 
library(Seurat)
library(AUCell)
library(GSEABase)
library(stringr)
library(writexl)

#reading all the data and subsetting only the mammary organoid specific cells
Mammary2.data <- Read10X(data.dir = "./../filtered_feature_bc_matrix/")
x <- colnames(Mammary2.data)
y1 <- str_detect(x, "-2$|-3$")
df1 <- Mammary2.data[, y1]
y2 <- str_detect(x,"-4$|-5$")
df2 <- Mammary2.data[, y2]
Mammary2 <- CreateSeuratObject(counts = df1, project = "Mammary2", min.cells = 5, min.features = 500, max.features = 7000)
Skin2 <- CreateSeuratObject(counts = df2, project = "Skin2", min.cells = 5, min.features = 500, max.feature = 7000)
Mammary2[["biotype"]] <- "Mammary2"
Skin2[["biotype"]] <- "Skin2"
add.cell.ids <- c("Mammary2", "Skin2")
gcdata <- merge(x = Mammary2, y = Skin2, add.cell.ids = add.cell.ids, merge.data = FALSE)
ifnb.list <- SplitObject(gcdata, split.by = "biotype")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "biotype")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "biotype")

#write.table(immune.combined@assays$integrated@scale.data,"top_2000_scaled_data_both_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

#FeaturePlot(immune.combined, features = c("Adipoq","Mitf","Acta2", "Krt18","Krt5", "Krt14", "Krt17"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("Prrx1"), min.cutoff = "q9", pt.size = 1)
VlnPlot(object = immune.combined, features = 'Ccnd1', split.by = 'biotype', cols=c('#D93128','#097BED'))

rm(df1,df2,Mammary2.data,gcdata,ifnb.list,immune.anchors)

#write.table(immune.combined@meta.data,"./meta_data_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

#write.table(immune.combined@reductions$umap@cell.embeddings,"./umap_data_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

comp.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

comp.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(immune.combined, features = top10$gene, group.colors = c("#db6d00","#004949","#ffb6db","#920000","#009292","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#ff6db6","#924900","#000000","#24ff24","#ffdf00")) + NoLegend()


comp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
df <- comp.markers
#write_xlsx(df, "combined_markers_clust_0.5_res.xlsx")

new.cluster.ids <- c("Mesenchyme-1", "Mammary Mesenchyme-1", "Mesenchyme-2", "Neuronal Derivatives", "Mammary Mesenchyme-2", "Chondrocytes", "Dermal Mesenchyme", "Neuronal cells","Myocytes","Dividing cells","Keratinocytes","Adipocytes","Neural crest cells","Myeloid cells","Melanocytes")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined<- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 1)+ NoLegend()

markers.to.plot <- c("Dct","Pmel","Mlana","Mmp12","Lyz2","Ccl9","Ngfr","Plp1","Sox10","Lpl","Plin2","Adipoq","Krt14","Krt8","Krt6a","Krt7","Krt18",
                     "Top2a","Mki67","Stmn1","Acta1","Mylpf","Myl1","Ptn","Ckb","Nrxn1","Col8a1","Ogn","Abi3bp","Col2a1","Col11a1","Acan","Rarres2","Serping1","Ifitm1",
                     "Crabp1","Nptx2","Vcan","S100a4","Tgfbi","Tm4sf1","Gpc6","Bmp2","Ero1l","Hspa5","Ndrg1")
DotPlot(immune.combined, features = markers.to.plot, cols = c("red", "blue"), dot.scale = 8, split.by = "biotype") + RotatedAxis()


