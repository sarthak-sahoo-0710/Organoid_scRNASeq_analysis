# reading required_packages
library(Seurat)
library(patchwork)
library(dplyr) 
library(Seurat)
library(AUCell)
library(GSEABase)
library(stringr)
library(writexl)
library(ggplot2)

#reading all the data and subsetting only the mammary organoid specific cells
Mammary2.data <- Read10X(data.dir = "./../filtered_feature_bc_matrix/")
x <- colnames(Mammary2.data)
y1 <- str_detect(x, "-2$")
df1 <- Mammary2.data[, y1]
y2 <- str_detect(x,"-3$")
df2 <- Mammary2.data[, y2]
Mammary2 <- CreateSeuratObject(counts = df1, project = "Mammary1", min.cells = 5, min.features = 500, max.features = 7000)
Skin2 <- CreateSeuratObject(counts = df2, project = "Mammary2", min.cells = 5, min.features = 500, max.feature = 7000)
Mammary2[["biotype"]] <- "Mammary1"
Skin2[["biotype"]] <- "Mammary2"
add.cell.ids <- c("Mammary1", "Mammary2")
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
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "biotype")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "biotype")

#write.table(immune.combined@assays$integrated@scale.data,"top_2000_scaled_data_both_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

FeaturePlot(immune.combined, features = c("Lbh"), min.cutoff = "q9", pt.size = 0.6)
VlnPlot(object = immune.combined, features = c("Krt5","Krt8","Krt14","Krt18","Krt6a","Cdh1","Epcam","Foxa1","Cd24a","Sox9","Nrp2","Areg"),pt.size=0.8,cols=c("#db6d00","#004949","#ffb6db","#920000","#009292","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#ff6db6","#924900","#24ff24","#454545","#ffdf00"))+NoLegend()
VlnPlot(object = immune.combined, features = "Lbh",pt.size=0.2,cols=c("#db6d00","#004949","#ffb6db","#920000","#009292","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#ff6db6","#924900","#24ff24","#454545","#ffdf00"))+NoLegend()

rm(df1,df2,Mammary2.data,gcdata,ifnb.list,immune.anchors)

#write.table(immune.combined@meta.data,"./meta_data_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

#write.table(immune.combined@reductions$umap@cell.embeddings,"./umap_data_combined.txt",sep="\t",row.names = TRUE, quote = FALSE)

comp.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

comp.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(immune.combined, features = top10$gene,group.colors =c("#db6d00","#004949","#ffb6db","#920000","#009292","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#ff6db6","#924900","#24ff24","#454545","#ffdf00")) + NoLegend()

comp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
df <- comp.markers
write_xlsx(df, "combined_markers_clust_0.5_res.xlsx")

new.cluster.ids <- c("Mesenchyme-1","Mammary Mesenchyme-1","Mammary Mesenchyme-2","Mesenchyme-2","Mammary Mesenchyme-3","Neuronal cells-1","Mitochondria-specific","Myocytes-1","Adipocytes","Myocytes-2","Dividing cells","Chondrocytes-1","Chondrocytes-2","Mammary epithelium","Neuronal cells-2")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined<- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "umap", label = TRUE)+ NoLegend()

markers.to.plot <- c("Mmp12","Lyz2","Ccl9","Ngfr","Plp1","Sox10","Lpl","Plin2","Adipoq","Krt14","Krt8","Krt6a","Krt7","Krt18",
                     "Top2a","Mki67","Stmn1","Acta1","Mylpf","Myl1","Ptn","Ckb","Nrxn1","Col8a1","Ogn","Abi3bp","Col2a1","Col11a1","Acan","Rarres2","Serping1","Ifitm1",
                     "Crabp1","Nptx2","Vcan","S100a4","Tgfbi","Tm4sf1","Gpc6","Bmp2","Ero1l","Hspa5","Ndrg1")
markers.to.plot <- c("Krt5","Krt6a","Krt8","Krt14","Krt18","Krt19","Cldn4","Cldn7","Epcam","Cdh1","Cd24a","Col2a1","Col11a1","Acan","Top2a","Mki67","Ccnb1","Acta1","Mylpf","Myl1","Lpl","Rbp4","Adipoq","Ptn","Ckb","Nrxn1","Sox9","Sox4","Bmp2","Ccnd1","Tgfbi","Col1a1","Vcan","S100a4")
markers.to.plot <- c("Lad1","Cxcl15","Pkp3","Dsg2","Itgb4","Jup")
markers.to.plot <- c("Pax3","Pax6", "Otx2", "Gbx2", "Bmp4", "Hand1")

#celltypes.to.plot <- c('Neuronal cells-1', 'Neuronal cells-2', 'Mammary epithelium', 'Chondrocytes-1','Dividing cells', 'Myocytes-2','Adipocytes',"Mesenchyme-1","Mammary Mesenchyme-1","Mammary Mesenchyme-2","Mesenchyme-2","Mammary Mesenchyme-3")
celltypes.to.plot <- c('Mammary epithelium','Dividing cells', 'Chondrocytes-1', 'Myocytes-2','Adipocytes', 'Neuronal cells-1',"Mesenchyme-1","Mammary Mesenchyme-1","Mammary Mesenchyme-2")
DotPlot(immune.combined, idents = celltypes.to.plot,features = markers.to.plot,cluster.idents=FALSE, dot.scale = 8) + RotatedAxis() + scale_colour_gradient2(low = "#D3D3D3",mid="#87CEEB", high = "#0000FF")

plots <- VlnPlot(immune.combined, features = c("Adipoq", "Cldn7"),idents = celltypes.to.plot, pt.size = 0.5, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


library(ggplot2)
data <- read.csv("GOdata_mouse_mam_mes_inte.txt", sep ="\t", header = TRUE, stringsAsFactors = FALSE)
df2=data[order(data$Fold_Enrichment),]
df2$Fold_Enrichment=factor(df2$Fold_Enrichment,levels=df2$Fold_Enrichment)
S1<- ggplot(df2, aes(x= Fold_Enrichment, y=Pathways, size=Genes, color=p_value)) + geom_point(alpha = 0.8) + 
  theme_classic() + scale_colour_gradient2(low = "#FF0000",mid="#0000FF", high = "#FF0000") + theme(legend.position="bottom")
S1
