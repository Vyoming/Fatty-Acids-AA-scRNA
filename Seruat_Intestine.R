# Load in Packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
# Load in the Data
Vyom_color_Umap <- CustomPalette(low = "darkblue", high = "yellow", mid = "#3cb371", k = 100)

contr4.data <- Read10X(data.dir ="~/data/intestine/SB04/Control/filtered_feature_bc_matrix/")
contr4 <- CreateSeuratObject(counts = contr4.data, project = "Control4")
contr4

AA4.data <- Read10X(data.dir = "~/data/intestine/SB04/AA/filtered_feature_bc_matrix/")
AA4 <- CreateSeuratObject(counts = AA4.data, project = "AA4")
AA4

Data4 <- merge(contr4, y = c(AA4), add.cell.ids = c("Control", "AA"), project = "SB04")
Data4
Data4[["percent.mt"]] <- PercentageFeatureSet(Data4, pattern = "mt-")
VlnPlot(Data4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Data4 <- subset(Data4, subset = nCount_RNA > 10000 & nCount_RNA < 100000 & nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 15)
Data4
VlnPlot(Data.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'type', pt.size = 0)
AA5.data <- Read10X(data.dir = "~/data/intestine/SB05/AA/")
AA5 <- CreateSeuratObject(counts = AA5.data, project = "AA5")
AA5

Control5.data <- Read10X(data.dir = "~/data/intestine/SB05/Control/")
C5 <- CreateSeuratObject(counts = Control5.data, project = "C5")
C5

P5.data <- Read10X(data.dir = "~/data/intestine/SB05/Pge2/")
P5 <- CreateSeuratObject(counts = P5.data, project = "P5")
P5

Data5 <- merge(C5, y = c(AA5, P5), add.cell.ids = c("Control", "AA", "PGE2"), project = "SB05")
Data5
Data5[["percent.mt"]] <- PercentageFeatureSet(Data5, pattern = "mt-")

VlnPlot(Data5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Data5 <- subset(Data5, subset = nCount_RNA > 3000 & nCount_RNA < 35000 & nFeature_RNA > 1500 & nFeature_RNA < 7000 & percent.mt < 15)
Data5

Data <- merge(Data4, y = c(Data5), add.cell.ids = c("SB04", "SB05"), project = "intest")
Data
# Begin Preprocessing and Filtering
Data[["percent.mt"]] <- PercentageFeatureSet(Data, pattern = "mt-")
VlnPlot(Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(Data@meta.data, 5)
plot1 <- FeatureScatter(Data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Data

# Batch correct data
Data.list <- SplitObject(Data, split.by = "ident")
Data.list <- Data.list[c("AA4", "AA5", "C5", "Control4", "P5")]
for (i in 1:length(Data.list)) {
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
options (future.globals.maxSize = 4000 * 1024^5)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
rm(AA4, AA5, C5, contr4, P5, Data4, Data5,AA4.data, AA5.data, C5.data, contr4.data, P5.data)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, Data)
Data.integrated <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)
rm(Data.anchors)
# Visulization and Clustering
Data.integrated <- RunPCA(Data.integrated, verbose = FALSE)

Data.integrated <- RunUMAP(Data.integrated, dims = 1:10)

levels(factor(Data.integrated@meta.data$orig.ident))
Idents(Data.integrated) <- Data.integrated$orig.ident
Data.integrated[["Type"]] <- Idents(Data.integrated)
new.cluster.ids <- c("AA", "AA", "Control", "Control", "Pge2")
names(new.cluster.ids) <- levels(Data.integrated)
Data.integrated <- RenameIdents(Data.integrated, new.cluster.ids)

Data.integrated[["Type"]] <- Idents(Data.integrated)
head(Idents(Data.integrated))
plots <- DimPlot(Data.integrated, group.by = c("Type"))
plots 

Data.integrated <- FindNeighbors(Data.integrated, dims = 1:10)
Data.integrated <- FindClusters(Data.integrated, resolution = .3)#.3
DimPlot(Data.integrated, reduction = "umap", label = TRUE)
#dotplots
b <- DotPlot(Data.integrated, features = c("Lgr5","Ascl2", "Hopx", "Clu", "S100a6","Hopx", "Defa24", "Tubb5", "Fabp1","Ube2c", "Muc2","Slc39a2", "Adh1"))
b


#Export <- CreateSeuratObject(counts = Data.integrated)
#Export



markers <- FindAllMarkers(Data.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(markers,'Sigs_Per_Clust.csv')

FeaturePlot(Data.integrated, features = c("S100a6", "Gata3"), blend = TRUE)
DE_AA_Control <- FindMarkers(Data.integrated, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", logfc.threshold = log(2), min.pct = 0.5)

VlnPlot(Data.integrated, features = c("S100a6"), ncol = 3, group.by = "Type", pt.size = 0, assay = "RNA", y.max = 5)
RidgePlot(Data.integrated, features = c("Wnt_targets1", "Yap_targets1", "Creb_targets1"), ncol = 3,  group.by= "Cell_Type")
VlnPlot(Data.integrated, features = c("Wnt_targets1", "Yap_targets1", "Creb_targets1"), ncol = 3,  group.by= "Cell_Type", pt.size = 0,y.max = NULL)& geom_boxplot(width=0.2, fill="#A4A4A4",color="#585858", outlier.shape = NA)
Data.integrated$Cell_Type

DotPlot(Data.integrated, features = yap.targets, group.by = "Type", assay = "RNA")
DotPlot(Data.integrated, features = Wnt.targets,  group.by = "Type")
DotPlot(Data.integrated, features = Creb.targets,  group.by = "Type")
names(x = Data.integrated[[]])
CustomPalette(low = "blue", high = "red", mid = NULL, k = 50)
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Wnt.targets), name = 'Wnt_targets')
FeaturePlot(object = Data.integrated, features = 'Wnt_Express1', cols = Vyom_color_Umap)
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(yap.targets), name = 'Yap_targets')
FeaturePlot(object = Data.integrated, features = 'yap.targets1', cols = Vyom_color_Umap)
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Wnt.targets), name = 'Creb_targets')
FeaturePlot(object = Data.integrated, features = 'Creb.targets1', cols = Vyom_color_Umap)

pdf(file = "~/CREB_Targets_UMAP.pdf", width = 7.5, height = 5)
FeaturePlot(object = Data.integrated, features = 'Creb.targets1', cols = Vyom_color_Umap)
dev.off()
FeaturePlot(object = Data.integrated, features = 'Type', cols = Vyom_color_Umap)

#annotate cluster decision
stem.markers <- c("Lgr5","Ascl2","Slc12a2","Axin2","Olfm4","Gkn3")
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(stem.markers), name = 'stem')
FeaturePlot(object = Data.integrated, features = 'stem1', cols = Vyom_color_Umap)


revival.markers <- c("Clu","Areg","Anxa1","Anxa2","Krt19","Ccnd2","Cxadr","Basp1","Ly6d","Lamc2","Ly6a","Krt7","S100a6","Reg3g","Krt8","Cldn4","Cxcl16","Malat1","Cdkn1a","Cldn3")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(revival.markers), name = 'revival.markers',assay =  'RNA')
FeaturePlot(object = Data.integrated, features = 'revival.markers1', cols = Vyom_color_Umap)


reserve.markers <- c("Bmi1","Lrig1","Hopx","Mex3a","Tert")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(reserve.markers), name = 'reserve.markers')
FeaturePlot(object = Data.integrated, features = 'reserve.markers1', cols = Vyom_color_Umap)


TA.markers <- c("Tubb5","Hmgb2","Stmn1","H2afz1","Tuba1b","Hmgb11","Hmgn22","Ptma1","Kiaa0101","Tk1","Cenpw","Tyms","Ranbp11","Idh21","Ran1","Dtymk","Nucks11","Dut1","Cks1b")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(TA.markers), name = 'TA.markers')
FeaturePlot(object = Data.integrated, features = 'TA.markers1', cols = Vyom_color_Umap)

paneth.markers <- c("Lyz1","Defa17","Defa22","Defa24","Ang4")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(paneth.markers), name = 'paneth.markers')
FeaturePlot(object = Data.integrated, features = 'paneth.markers1', cols = Vyom_color_Umap)

Enterocyte.markers <- c("Alpi","Apoa1","Apoa4","Fabp1","Adh6a")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
FeaturePlot(object = Data.integrated, features = 'Enterocyte.markers1', cols = Vyom_color_Umap)

EP.markers <- c("Ccnb1","Cdc20","Cenpa","Cdkn3","Ccnb2","Cdc25c","Kif22","Ube2c","Sapcd2","Rbp7","Aurka","Ccna2","Cdkn2d","Kif23","Nek2","Slc16a1")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(EP.markers), name = 'EP.markers')
FeaturePlot(object = Data.integrated, features = 'EP.markers1', cols = Vyom_color_Umap)

Enteroendocrine.markers <- c("Chga","Chgb","Tac1","Tph1","Neurog3","Gch1","Slc39a2")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
FeaturePlot(object = Data.integrated, features = 'Enteroendocrine.markers1', cols = Vyom_color_Umap)

Goblet.markers <- c("Muc2","Clca3","Tff3","Agr2","Spink4","Fcgbp","Zg16","Ccl9","Atoh1")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Goblet.markers), name = 'Goblet.markers')
FeaturePlot(object = Data.integrated, features = 'Goblet.markers1', cols = Vyom_color_Umap)

Tuft.markers <- c("Trpm5","Gfi1b","Il25","Alox5ap","Lrmp","Rgs13","Ltc4s","Adh1")

Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Tuft.markers), name = 'Tuft.markers')
FeaturePlot(object = Data.integrated, features = 'Tuft.markers1', cols = Vyom_color_Umap)

Sigs_per_clust <- DotPlot(Data.integrated, group.by = "Cell_Type", features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))
Sigs_per_clust
VlnPlot(Data.integrated, pt.size = 0,group.by = "Cell_Type", features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))
VlnPlot(Data.integrated, pt.size = 0,group.by = "Cell_Type", features = c("Tubb5","Hmgb2","Stmn1","H2afz1","Tuba1b","Hmgb11","Hmgn22","Ptma1","Kiaa0101","Tk1","Cenpw","Tyms","Ranbp11","Idh21","Ran1","Dtymk","Nucks11","Dut1","Cks1b"), assay= "RNA")
AverageExpression(Data.integrated,assay ='RNA', slot = "data", features = check)

sanity_check<- VlnPlot(Data.integrated, pt.size = 0,ncol = 4,group.by = "Cell_type", features = check, assay= "RNA")

stem_score <- Data.integrated$revival.markers1

Cell_Type_Sig_Score <- data.frame(Stem=Data.integrated$stem1, Revival=Data.integrated$revival.markers1, Reserve=Data.integrated$reserve.markers1, Transit_Amplifying =Data.integrated$TA.markers1, Paneth=Data.integrated$paneth.markers1, Enterocyte=Data.integrated$Enterocyte.markers1, EP=Data.integrated$EP.markers1,Enteroendocrine= Data.integrated$Enteroendocrine.markers1,Goblet=Data.integrated$Goblet.markers1,Tuft=Data.integrated$Tuft.markers1)

write.csv(Cell_Type_Sig_Score, "Cell_Type_Score.csv", row.names = TRUE)

# module score distribution
modulescores <- Cell_Type_Sig_Score %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
   facet_wrap(~celltype) +
   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

 p <- ggplot(modulescores)
 p <- ggplot(onescore)
 p + geom_histogram(aes(x=score)) 
 
# Cowplots for Single Cell Analysis
 
# Umap
cluster_umap <-DimPlot(Data.integrated, reduction = "umap")
cluster_umap
#pathway sigs
rev.mark <- c('S100a6',	'Clu',	'Ly6a',	'S100a11',	'Ly6d',	'Sprr2a3',	'Mif',	'Krt19',	'Krt7',	'Npm1',	'Krt18',	'Rps12-ps3',	'Fxyd3',	'Xcl1',	'Anxa2',	'Ccnd1',	'Ncl',	'Ccnd2',	'Cd44',	'Ptma',	'Anxa3',	'Lmna',	'Ranbp1',	'Rplp0',	'Rpl29',	'Nme1',	'Ran',	'Gm10073',	'Rpl35',	'Dynll1',	'Sprr1a',	'Rps2',	'Rps12',	'Anxa1',	'Rps27l',	'Erh',	'Gm10076',	'Rplp1',	'Pglyrp1',	'Rps28',	'Nhp2',	'Snrpg',	'Rps26',	'Gm8730',	'Sat1',	'Eif2s2',	'Wdr89',	'Hsp90ab1',	'Rpl41',	'Rpl27')
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(rev.mark), name = 'rev_mark')
FeaturePlot(object = Data.integrated, features = 'rev_mark1', cols = Vyom_color_Umap)

names(x = Data.integrated[[]])
CustomPalette(low = "blue", high = "red", mid = NULL, k = 50)
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Wnt.targets), name = 'Wnt_Express')
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(yap.targets), name = 'yap.targets')
Data.integrated <- AddModuleScore(object = Data.integrated, features = list(Wnt.targets), name = 'Creb.targets')

wnt_umap <- FeaturePlot(object = Data.integrated, features = 'Wnt_Express1', cols = Vyom_color_Umap)
yap_umap <- FeaturePlot(object = Data.integrated, features = 'yap.targets1', cols = Vyom_color_Umap)
creb_umap <- FeaturePlot(object = Data.integrated, features = 'Creb.targets1', cols = Vyom_color_Umap)

FeaturePlot(object = Data.integrated, features = 'Pou2f3', cols = Vyom_color_Umap)

VlnPlot(Data.integrated, features = c(Wnt.targets),pt.size = .01, ncol = 1)
VlnPlot(Data.integrated, features = c("Wnt_Express1"),pt.size = .01, group.by = "orig.ident", ncol = 1)

#supplementary figures
#QC
Sigs_per_clust <- DotPlot(Data.integrated, features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))
Sigs_per_clust

#dot plots
Sigs_per_clust <- DotPlot(Data.integrated, features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))
Sigs_per_clust
yap_dot <-DotPlot(Data.integrated, features = yap.targets)
wnt_dot <- DotPlot(Data.integrated, features = Wnt.targets, assay = 'RNA')
Creb_dot <- DotPlot(Data.integrated, features = Creb.targets, assay = 'RNA')


#differential expression analysis for all cells
Idents(Data.integrated) = Data.integrated$type
AAvC1 <-FindMarkers(Data.integrated, ident.1 = "AA", ident.2 = "Control", test.use = "MAST")
  write.csv(AAvC,'AA_vs_Contr_DE.csv')
EnhancedVolcano(PGEvC,
                lab = rownames(PGEvC),
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-5, 5))
AAvC <-FindMarkers(Data.integrated, ident.1 = "Control", ident.2 = "AA", test.use = "MAST")
write.csv(AAvC,'AA_vs_Contr_DE.csv')

PGEvC <-FindMarkers(Data.integrated, ident.1 = "Control", ident.2 = "Pge2", test.use = "MAST")
write.csv(PGEvC,'PGE2_vs_Contr_DE.csv')
check_sanity = FALSE

Rev1vsRev2 <-FindMarkers(Data.integrated, ident.1 = "Revival Stem 1", ident.2 = "Revival Stem 2", test.use = "MAST")
write.csv(Rev1vsRev2,'Rev1_vs_Rev2_DE.csv')

Data.integrated1 <-Data.integrated
Data.integrated1 <- NormalizeData(object = Data.integrated1,normalization.method = "LogNormalize", assay = "RNA")
AAvC_YAP <-FindMarkers(Data.integrated1, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", features = yap.targets, assay = "RNA")
write.csv(AAvC_YAP,'AA_vs_Contr_DESeq2_YAP.csv')

AAvC_YAP <-FindMarkers(Data.integrated1, ident.1 = "Pge2", ident.2 = "Control", test.use = "DESeq2", features = yap.targets, assay = "RNA")
write.csv(AAvC_YAP,'Pge2_vs_Contr_DESeq2_YAP.csv')

AAvC_wnt <-FindMarkers(Data.integrated1, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", features = Wnt.targets, assay = "RNA")
write.csv(AAvC_wnt,'AA_vs_Contr_DESeq2_WNT.csv')

AAvC_wnt <-FindMarkers(Data.integrated1, ident.1 = "Pge2", ident.2 = "Control", test.use = "DESeq2", features = Wnt.targets, assay = "RNA")
write.csv(AAvC_wnt,'Pge2_vs_Contr_DESeq2_WNT.csv')

AAvC_CREB <-FindMarkers(Data.integrated1, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", features = Creb.targets, assay = "RNA")
write.csv(AAvC_CREB,'AA_vs_Contr_DESeq2_WNT.csv')

AAvC_CREB <-FindMarkers(Data.integrated1, ident.1 = "Pge2", ident.2 = "Control", test.use = "MASDESeq2T", features = Creb.targets, assay = "RNA")
write.csv(AAvC_CREB,'Pge2_vs_Contr_DESeq2_WNT.csv')
Idents(Data.integrated)
DimPlot(Data.integrated)
#Further subcluster 5 to get paneth tuft and goblet

Cluster_7 <- subset(Data.integrated, idents = 7)

Cluster_7
Cluster_7 <- RunPCA(Cluster_7, verbose = FALSE)

Cluster_7 <- RunUMAP(Cluster_7, dims = 1:10)

Cluster_7 <- FindNeighbors(Cluster_7, dims = 1:10)
Cluster_7 <- FindClusters(Cluster_7, resolution = .3)
DimPlot(Cluster_7, reduction = "umap", label = TRUE)


Cluster_7 <- AddModuleScore(object = Cluster_7, features = list(reserve.markers), name = 'reserve.markers', assay = 'RNA')
FeaturePlot(object = Cluster_7, features = 'reserve.markers1', cols = Vyom_color_Umap)
Cluster_7 <- AddModuleScore(object = Cluster_7, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
FeaturePlot(object = Cluster_7, features = 'Enteroendocrine.markers1', cols = Vyom_color_Umap)

new.cluster.ids <- c("Enteroendocrine","Stem 4","Stem 4","Enteroendocrine","Enteroendocrine")
Cluster_7[["Cell_Type"]] <- Idents(Cluster_7)
names(new.cluster.ids) <- levels(Cluster_7)
Cluster_7 <- RenameIdents(Cluster_7, new.cluster.ids)
Cluster_7[["Cell_Type"]] <- Idents(Cluster_7)
DimPlot(Cluster_7, reduction = "umap", label = TRUE)

Cluster_5 <- subset(Data.integrated, idents = 5)
Cluster_5
Cluster_5 <- RunPCA(Cluster_5, verbose = FALSE)

Cluster_5 <- RunUMAP(Cluster_5, dims = 1:10)

Cluster_5 <- FindNeighbors(Cluster_5, dims = 1:10)
Cluster_5 <- FindClusters(Cluster_5, resolution = .3)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Goblet.markers), name = 'Goblet.markers')
FeaturePlot(object = Cluster_5, features = 'Goblet.markers1', cols = Vyom_color_Umap)
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Tuft.markers), name = 'Tuft.markers')
FeaturePlot(object = Cluster_5, features = 'Tuft.markers1', cols = Vyom_color_Umap)
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(paneth.markers), name = 'paneth.markers')
FeaturePlot(object = Cluster_5, features = 'paneth.markers1', cols = Vyom_color_Umap)
new.cluster.ids <- c("Goblet","Goblet","Stem Like Secretory Progenitor","Paneth","Paneth","Tuft")
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
names(new.cluster.ids) <- levels(Cluster_5)
Cluster_5 <- RenameIdents(Cluster_5, new.cluster.ids)
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


#pre annotate real clusters
new.cluster.ids <- c("Stem 1","Transit Amplifying","Enterocyte Progenitor","Stem 1","Enterocyte", "","Stem 2", "", "Stem 3", "Stem Like Progenitor")
Data.integrated[["Cell_type"]] <- Idents(Data.integrated)
names(new.cluster.ids) <- levels(Data.integrated)
Data.integrated <- RenameIdents(Data.integrated, new.cluster.ids)
Data.integrated[["Cell_type"]] <- Idents(Data.integrated)
DimPlot(Data.integrated, reduction = "umap", label = TRUE)


# Generate a new column called Cell_type in the metadata
Data.integrated$Cell_Type <- as.character(Idents(Data.integrated))

# Change the information of cells containing sub-cluster information
Data.integrated$Cell_Type[WhichCells(Cluster_5)] <- paste("",Idents(Cluster_5))
Data.integrated$Cell_Type[WhichCells(Cluster_7)] <- paste("",Idents(Cluster_7))
DimPlot(Data.integrated, group.by = "Cell_Type", label = TRUE, pt.size=1, label.size = 3) + NoLegend()
Idents(Data.integrated) <- Data.integrated$Cell_Type
levels(factor(Data.integrated@meta.data$Type))
Data.integrated@meta.data$type <- factor(Data.integrated@meta.data$Type, levels = c("Control", "AA", "Pge2")) 
prop.table(x = table(Idents(Data.integrated), Data.integrated$Type), margin = 2)

Idents(Data.integrated) = Data.integrated$Cell_Type

Data.integrated$integrated_snn_res.0.3
DimPlot(Data.integrated, group.by = "integrated_snn_res.0.3", label = TRUE, pt.size=1, label.size = 3) + NoLegend()
VlnPlot(Data.integrated, group.by = 'integrated_snn_res.0.3', split.by ='Type', features = c('Xpo1'), assay = 'RNA')
FeaturePlot(Data.integrated,features = '')
Idents(Data.integrated) <- Data.integrated$Cell_Type
stem_cells <- subset(Data.integrated, idents = 'Stem')

stem_cells
stem_cells <- RunPCA(stem_cells, verbose = FALSE)

stem_cells <- RunUMAP(stem_cells, dims = 1:10)

stem_cells <- FindNeighbors(stem_cells, dims = 1:10)
stem_cells <- FindClusters(stem_cells, resolution = .15)
DimPlot(stem_cells, reduction = "umap", label = TRUE)
stem_cells <- NormalizeData(object = stem_cells,normalization.method = "LogNormalize", assay = "RNA")

stem.markers <- c("Lgr5","Ascl2","Axin2","Olfm4","Gkn3",'Birc5', 'Hmgb2', 'Ccnb2','Hoxb7', 'Lyzl4', "S100a6",'Mki67','Xpo1', 'Lyz1','Stmn1')
stem_cells <- NormalizeData(object = stem_cells,normalization.method = "LogNormalize", assay = "RNA")
VlnPlot(stem_cells, group.by = 'integrated_snn_res.0.15', split.by = "type", ncol = 5, slot = "data", features = stem.markers, pt.size = 0, assay = "RNA", log = TRUE)
stem_cells$integrated_snn_res.0.15
stem.markers <- c("Lgr5","Ascl2","Olfm4","","Olfm4","Gkn3")
stem_cells <- AddModuleScore(object = stem_cells, features = list(stem.markers), name = 'stem', assay = 'RNA')
stem_cells <- AddModuleScore(object = stem_cells, features = 'Xpo1', name = 'gene', assay = 'RNA')
stem_cells$type
FeaturePlot(object = stem_cells, features = 'gene1', cols = Vyom_color_Umap)

markers_stem <- FindAllMarkers(stem_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_stem %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# cell based DE analysis
Data.integrated1 <-Data.integrated
Data.integrated1 <- NormalizeData(object = Data.integrated1,normalization.method = "LogNormalize", assay = "RNA")
gene.level.list <- c('Lgr5', 'Ascl2', 'Olfm4', 'Clu', 'Ly6a', 'S100a6', 'Bmi1', 'Hopx', 'Lrig1', 'Mex3a')
Idents(Data.integrated1) <- Data.integrated1$Type
Sub_cluster <- subset(Data.integrated1, idents = c('Control','AA'))
Gene_level_DE <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene_level_list, assay = "RNA", logfc.threshold = 0)
AverageExpressio
Gene_level_DE
Idents(Sub_cluster) <- Sub_cluster$Cell_Type
VlnPlot(Sub_cluster, group.by = "Type", features = gene.level.list, pt.size = 0, assay = "RNA")
DotPlot(Sub_cluster, group.by = "Type", features = gene_level_list, dot.scale = 6, assay = "RNA")
DotPlot(Sub_cluster, group.by = "Cell_Type", features = gene_level_list, dot.scale = 6, assay = "RNA")
ascl_promoters <- ASCL2_targets_promoter$...1

Sub_cluster <- subset(Data.integrated1, idents = c('Control', 'Pge2'))
Idents(Sub_cluster) <- Sub_cluster$Cell_Type
VlnPlot(Sub_cluster, split.by = "Type", slot = "data", features = 'ascl_promoters1', pt.size = 0, assay = "RNA")
Idents(Data.integrated) <- Data.integrated$Type
Sub_cluster <- subset(Data.integrated, idents = c('Control','AA'))
Idents(Sub_cluster) <- Sub_cluster$Cell_Type
VlnPlot(Sub_cluster, split.by = "Type", slot = "data", features = 'Cd55', pt.size = 0, assay = "RNA")

DotPlot(Sub_cluster, split.by = "Type", group.by = "Cell_Type", cols= c('darkblue','Red'),features = "Cd55", dot.scale = 6, assay = "RNA")

AverageExpression(Gene_de,assay ='RNA', slot = "data", features = gene.level.list)
VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 3, split.by = "type", slot = "data", features = gene.level.list, pt.size = 0, assay = "RNA", log = TRUE, cols = c('#00BFC4','#F8766D'))
Idents(Data.integrated1) <- Data.integrated1$Cell_type
Data.integrated1 <- AddModuleScore(object = Data.integrated1, features = list(ascl_promoters), name = 'ascl_promoters')
VlnPlot(stem_cells, features = 'Cd55', group.by = 'integrated_snn_res.0.15', assay = "RNA")
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AverageExpression(arasco.obj, assay ='RNA', slot = "data", features = List)
VlnPlot(Sub_cluster, group.by = "Cell_Type", ncol = 3, split.by = "Type", slot = "data", features = gene.level.list, pt.size = 0, assay = "RNA", cols = c('#00BFC4','#F8766D'))

List <- c('Ptges', 'Ptgs1', 'Ptgs2', 'Ptgis', 'Ptger4', 'Ptgds', 'Cbr1', 'Ptges2', 'Ptges3', 'Hpgds', 'Pla2g4a')
my_levels <- c("Stem 1", "Stem 2", "Stem 3", " Stem 4", "Stem Like Progenitor", " Stem Like Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", " Enteroendocrine", " Goblet", " Tuft", " Paneth" )
Data.integrated$Cell_type <- factor(x = Data.integrated$Cell_Type, levels = my_levels)

VlnPlot(Sub_cluster, slot = "data",group.by = "Cell_type", split.by = "Type", features = gene.level.list, pt.size = 0,assay = "RNA")

DimPlot(Sub_cluster)
p = VlnPlot(arasco.obj, slot = "data",group.by = "Cell_Type", split.by = "Type", features = List, pt.size = 0,assay = "RNA")
p
p1 <- list() 
for (i in seq_along(p)){
  #Change x and y tick label font size.
  p1[[i]] = p[[i]] + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) 
}
plot_grid(plotlist = p1, ncol = 3) 
Vyom_color_Umap <- CustomPalette(low = "purple", high = "yellow", mid = "#3cb371", k = 100)
Vyom_color_Umap1 <- CustomPalette(low = "Darkblue", high = "red", mid = "purple", k = 100)

Gene_level_DE <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene_level_list, assay = "RNA", logfc.threshold = 0)

Sub_cluster <- subset(Data.integrated1, idents = c('Control'))
Dotc = DotPlot(Sub_cluster, group.by = "Cell_Type", cols= c('darkblue','Red'),features = 'Cd55', dot.scale = 6, assay = "RNA" ) & NoLegend()
Sub_cluster <- subset(Data.integrated1, idents = c('AA'))
Dota = DotPlot(Sub_cluster, group.by = "Cell_Type", cols= c('darkblue','Red'),features = "Cd55", dot.scale = 6, assay = "RNA") + theme(axis.title.y = element_blank())
plot_grid(Dotc,
          Dota,
          labels = c('Control','AA'),
          label_x = 0.5,
          ncol = 2)
# DE per cell type

Idents(Data.integrated) <- Data.integrated$Cell_Type
Sub_cluster <- subset(Data.integrated,  idents = "Stem")
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'STEM_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'STEM_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c("Revival Stem 1", "Revival Stem 2"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,' combrSTEM_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'combrSTEM_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c("Enterocyte Progenitor"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'EP_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'EP_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c("Transit Amplifying"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'TA_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'TA_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c("Enterocyte"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Ent_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Ent_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,idents =c(" Tuft"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Tuft_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Tuft_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c(" Paneth"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Paneth_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Paneth_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c(" Goblet"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Goblet_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Goblet_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c(" Enteroendocrine"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'entend_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'entend_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c("Stem Like Progenitor Cell"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Proj_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Proj_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c(" Stem & Tuft"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'ST_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'ST_PGE2_vs_Contr_DE_gene_level.csv')

Sub_cluster <- subset(Data.integrated,  idents =c(" Reserve Stem"))
Idents(Sub_cluster) = Sub_cluster$Type
Sub_cluster <- NormalizeData(object = Sub_cluster,normalization.method = "LogNormalize", assay = "RNA")
AAvCs <-FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(AAvCs,'Reserve_AA_vs_Contr_DE_gene_level.csv')
PGEvCs <-FindMarkers(Sub_cluster, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", features = gene.level.list, logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")
write.csv(PGEvCs,'Reserve_PGE2_vs_Contr_DE_gene_level.csv')

UMAP_Coords = Embeddings(Data.integrated@reductions$umap)

write.csv(UMAP_Coords,'WHALE_Coords.csv')

pbmc.loom <- as.loom(arasco.obj, filename = "/Users/vyom/data/Intestine/Arasco.loom", verbose = FALSE)
pbmc.loom
pbmc.loom$close_all()



yap.targets <-c('Cat',	'Ctgf',	'Cyr61',	'Gpatch4',	'Lmnb2',	'Ptgs2',	'Wsb2',	'Amotl2',	'Ankrd1',	'Igfbp3',	'F3',	'Fjx1',	'Nuak2',	'Lats2',	'Crim1',	'Gadd45a',	'Tgfb2',	'Ptpn14',	'Nt5e',	'Foxf2',	'Axl',	'Dock5',	'Asap1',	'Rbms3',	'Myof',	'Arhgef17',	'Ccdc80',	'Tnfrsf12a',	'Areg',	'Clu',	'Anxa3',	'Tinagl1',	'Gprc5a',	'Msln',	'Il33',	'Edn1',	'Plaur',	'S100a6',	'Birc5',	'Sox9',	'Ccnb2',	'Cdkn3',	'Tubb6',	'Anxa5',	'Elf3',	'Bcl2',	'Ccl4',	'Cldn4',	'Lyz2',	'Hkdc1',	'Cygb',	'Abi3',	'Slc4a9',	'C1qb',	'Lcp2',	'Lst1',	'Tacstd2',	'Csf1r',	'Pip4k2a',	'Pkib',	'Was',	'Cldn6',	'Evl',	'Fgl2',	'Prtn3',	'Ddit4l',	'Aif1',	'Mdk',	'Cyp1b1',	'Fkbp7',	'Arhgdib',	'Ncf4',	'Irf8',	'Omd',	'Wipf1',	'Arhgap30',	'Col6a1',	'S100a8',	'Ntrk2',	'Vcam1',	'Tyrobp')
Wnt.targets <-c(	'Adam10',	'Ascl2',	'Axin2',	'Bambi',	'Bcl2l2',	'Bcl2l1',	'Birc5',	'Bmi1',	'Bmp4',	'Ccnd1',	'Cd44',	'Cdkn2a',	'Cdx1',	'Cldn1',	'Dkk1',	'Dkk4',	'Dnmt1',	'Edn1',	'Efnb1',	'Enc1',	'Ephb2',	'Ephb3',	'Fgf18',	'Fgfbp1',	'Fgfbp3',		'Fscn1',	'Gast',		'Nedd9',	'Hes1',	'Id2',	'Tcf4',	'Jag1',	'Jun',	'L1cam',	'Lamc2',	'Lef1',	'Lgr5',	'Met',	'Mmp14',	'Mmp7',	'Myb',	'Myc',	'Mycbp',	'Nos2',	'Notch2',	'Nrcam',	'Plaur',	'Ppard',	'S100a4',	'Sgk1',	'Smc3',	'Sox9',	'Sp5',	'Srsf3',	'Suz12',	'Tiam1',	'Tnc',	'Yap1',	'Ccnd3',	'Fzd1',	'Ctnnb1',		'Rnf43',	'Lrp6',	'Lrp5',	'Olfm4',	'Prex1',	'Rac1',	'Cldn15',	'Edn3',	'Efnb2',	'Ephb4',	'Fgfr4',	'Nedd8',	'Tcf7',	'Jag2',	'Smc2',	'Zeb1',	'Zeb2',	'Cdh2',	'Cdh1',	'Actb',	'Gapdh',	'Znrf3',	'Lgr4',	'Mycn',	'Ccnd2',	'Hnf1a',	'Fosl1',	'Ephb1',	'Ephb6',	'Msl1',	'Dkk2',	'Dkk3',	'Dkkl1',	'Tert',	'Fgf9',	'Lbh',	'Fgf20',	'Sox17',	'Runx2',	'Grem1',	'Sall4',	'Sox2',	'Dll1',	'Foxn1',	'Nanog',	'Snai1',	'Fn1',	'Fzd7',	'Fst',	'Wnt3a',	'Isl1',	'Mmp2')
Creb.targets <-c('Hoxa1',	'Hoxa2',	'Hoxa3',	'Hoxa5',	'Hoxa11',	'Hoxb4',	'Hoxb5',	'Hoxc4',	'Hoxc5',	'Hoxc8',	'Hoxc9',	'Hoxc10',	'Hoxd8',	'Pitx2',	'Zic1',	'Dlx6',	'Hhex',	'Lhx1',	'Nkx2-2',	'Nr4a2',	'Otx2',	'Pax3',	'Pax5',	'Pax6',	'Pou4f2',	'Ppard',	'Shox2',	'Stat3',	'Nr2e1',	'Phox2b',	'Lhx2',	'Hey2',	'Heyl',	'Lhx5',	'Foxp2',	'Hira',	'Klf4',	'Gsc',	'Tbx5',	'Hand1',	'Hand2',	'E2f1',	'E2f3',	'E2f5',	'Erf',	'Fosb',	'Junb',	'Nfyc',	'Rb1',	'Thra',	'Ctcf',	'Bcl6',	'Zfp36l1',	'Zfp36l2',	'Btg1',	'Egr4',	'Irf2',	'Pou1f1',	'Btg2',	'Fosl1',	'Ccnl1',	'Brf1',	'Snapc5',	'Hmgb2',	'Sox1',	'Hmg20a',	'Sox21',	'Cebpe',	'Gata3',	'Nr3c1',	'Nfx1',	'Xbp1',	'Foxn1')
ASCL2.targets <- c('1190003M12Rik', 'Adfp', 'Adora1', 'Adra2a', 'AI173486', 'AI428936', 'Aqp4', 'Ascl2', 'BC057022', 'Cdc42ep1', 'Cdk6', 'Ces3', 'Cfi', 'Crlf1', 'E030025L21Rik', 'Ell3', 'Ephb3', 'Kif12', 'Lgr5', 'Mia1', 'NAP027049-1', 'Nr2e3', 'Olfm4', 'Psrc1', 'Ptpro', 'Slc12a2', 'Slc14a1', 'Slco3a1', 'Soat1', 'Sox9', 'Tnfrsf19', 'Zbtb12')
ASCL2.enhancers<- c('Agpat2',	'Ptges',	'Arhgap40',	'Pard6b',	'Pear1',	'Lmna',	'Bcar3',	'Coro2a',	'Klf4',	'Tek',	'Bmp8b',	'Tinagl1',	'Fam46b',	'Ldlrap1',	'Id3',	'Epha2',	'Efhd2',	'Espn',	'Mmp23',	'Hgfac',	'Por',	'Adap1',	'Tead4',	'Bcl3',	'Dedd2',	'4732471J01Rik',	'Sertad1',	'Sult2b1',	'Trim30d',	'Ano1',	'Ccnd1',	'Plat',	'Fcho1',	'Palm3',	'Apoa1',	'Paqr5',	'Rbp1',	'Als2cl',	'Vill',	'Plcd1',	'Csrnp1',	'Hk1',	'Celf5',	'Suox',	'Myl7',	'Wwc1',	'Pitpnm3',	'Itga3',	'Tns4',	'Krt19',	'Plcd3',	'Myo15b',	'Itgb4',	'Socs3',	'Slc16a3',	'5730507C01Rik',	'Id2',	'Fos',	'Slc9a3',	'Lhfpl2',	'Arhgap22',	'Ajuba',	'Pdlim2',	'Fam83h',	'Wnt7b',	'Dhh',	'Fignl2',	'Krt4',	'Kalrn',	'Pim1',	'Rnf39',	'Stap2',	'Ston1',	'Celf4',	'Smad7',	'Sh3pxd2a',	'Trex2')
ASCL2.promoters <- c('Tcf24',	'Dst',	'Tmem237',	'Epha4',	'Cd55',	'Camsap2',	'Ptpn14',	'Casc4',	'Slc4a11',	'4930402H24Rik',	'Atrn',	'Slc23a2',	'Wfdc2',	'Sulf2',	'Lama5',	'Anxa5',	'Abhd18',	'Wwtr1',	'Pear1',	'S100a6',	'Ppm1j',	'Npnt',	'Ttll7',	'Ak4',	'Bmp8b',	'Ncmap',	'Nbl1',	'Tmem51os1',	'Slc5a1',	'Dok7',	'Sh3tc1',	'Sgcb',	'Ereg',	'Prkg2',	'Ptpn13',	'Vps37d',	'Nyap1',	'Gpr146',	'Lncpint',	'Cyp26b1',	'B4galnt3',	'Tspan9',	'Tead4',	'Cox6b2',	'Kcnk6',	'Wtip',	'Ppfibp2',	'Mical2',	'Far1',	'Tspan4',	'Syt8',	'Atp11a',	'Micu3',	'Stox2',	'Wwc2',	'Pbx4',	'Junb',	'Kifc3',	'Maml2',	'Adamts15',	'Gramd1b',	'Grik4',	'Gnb5',	'Plod2',	'Ppp2r3a',	'Als2cl',	'Tmie',	'Slc6a20a',	'Samd5',	'Hebp2',	'Vsir',	'Hk1',	'Adamtsl5',	'Tmem263',	'Nab2',	'Sowaha',	'Chd3',	'Spns2',	'Mmp28',	'Srcin1',	'Itga2b',	'Sphk1',	'Slc16a3',	'Tc2n',	'Tnfaip2',	'Crip2',	'Foxq1',	'Pdlim7',	'Slc25a48',	'Fbp1',	'Lysmd3',	'Plau',	'Ogdhl',	'Arhgap22',	'Cdhr1',	'Fam213a',	'Ero1l',	'Lgals3',	'Ltb4r1',	'Neil2',	'Pdlim2',	'Sox21',	'Abcc4',	'Ndrg1',	'Ccdc134',	'Rarg',	'Zfp385a',	'Ppl',	'Emp2',	'Abcc5',	'Nxpe3',	'Dcbld2',	'Btg3',	'Bace2',	'Neurl1b',	'Cdkn1a',	'H2-Q1',	'Ppp1r18',	'Abhd3',	'Epb41l4a',	'Synpo',	'Tubb6',	'Clcf1',	'Fosl1',	'Ahnak',	'Gcnt1',	'Tcf7l2',	'Ablim1',	'Acsl4',	'Ap1s2')
Sub_cluster <- AddModuleScore(object = Sub_cluster, features = list(Wnt.targets), name = 'Wnt_targets', assay = 'RNA')
FeaturePlot(object = Data.integrated, features = 'Wnt_Express1', cols = Vyom_color_Umap)
Sub_cluster <- AddModuleScore(object = Sub_cluster, features = list(yap.targets), name = 'Yap_targets', assay = 'RNA')
FeaturePlot(object = Data.integrated, features = 'yap.targets1', cols = Vyom_color_Umap)
Sub_cluster <- AddModuleScore(object = Sub_cluster, features = list(Creb.targets), name = 'Creb_targets', assay = 'RNA')
FeaturePlot(object = Data.integrated, features = 'Creb.targets1', cols = Vyom_color_Umap)
Sub_cluster <- AddModuleScore(object = Sub_cluster, features = list(ASCL2.promoters), name = 'ASCL2.targets', assay = 'RNA')

# add decoy cell identity
Idents(Data.integrated) <- Data.integrated$Type

# compute average gene expression
Sub_cluster@meta.data$yap.targets <- AverageExpression(Sub_cluster, features = yap.targets, assays = 'RNA')

# add it to meta.data
Sub_cluster@meta.data$yap.targets <- gene.set.exp

#Gene Lebel Analysis
# average expression of gene set and create feature for violin plot
mean.exp <- colMeans(x = Gene_de@assays$RNA@data[yap.targets, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Yap_Targets <- mean.exp
}

mean.exp <- colMeans(x = Gene_de@assays$RNA@data[Wnt.targets, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$ßCatenin_Targets <- mean.exp
}

mean.exp <- colMeans(x = Gene_de@assays$RNA@data[Creb.targets, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$CREB_Targets <- mean.exp
}

mean.exp <- colMeans(x = Gene_de@assays$RNA@data[ASCL2.enhancers, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$ASCL2_Enhancers <- mean.exp
}

mean.exp <- colMeans(x = Gene_de@assays$RNA@data[ASCL2.promoters, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$ASCL2_Promoters <- mean.exp
}
VlnPlot(Gene_de, group.by = "Cell_Type", ncol = 2, split.by = "Type", slot = "data", features = c('Yap_Targets','ßCatenin_Targets','CREB_Targets','ASCL2_Targets'), pt.size = 0, assay = "RNA", cols = c('#00BFC4','#F8766D'))
VlnPlot(Gene_de, group.by = "Cell_type", split.by = "type", ncol = 3, slot = "data", features = c('Yap_Targets', 'ßCatenin_Targets','CREB_Targets','ASCL2_Enhancers', 'ASCL2_Promoters'), pt.size = 0, assay = "RNA", log = TRUE,cols = c('#00BFC4','#F8766D'))

VlnPlot(Gene_de, group.by = "Cell_type", ncol = 2, split.by = "Type", slot = "data", features = c('ASCL2_Enhnacers', 'ASCL2_Promoters' ), pt.size = 0, assay = "RNA", log = TRUE)

AAvC_path <-FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control",  features = c('Yap_Targets','ßCatenin_Targets','CREB_Targets','ASCL2_Enhnacers', 'ASCL2_Promoters'), logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")

All_prominant_genes <- c(yap.targets, Wnt.targets, Creb.targets,ASCL2.enhancers, ASCL2.promoters, gene.level.list)
All_prominant_genes <- unique(All_prominant_genes) 

AAvC_gene_level <-FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control",  features = All_prominant_genes , logfc.threshold = -Inf, min.pct = -Inf, assay = "RNA")

Gene_de <- NormalizeData(object = Gene_de,normalization.method = "LogNormalize", assay = "RNA")

Sub_cluster <- subset(Data.integrated,  idents =c(" Stem & Tuft"))
Data.integrated1 <- Data.integrated
Gene_de <- subset(Data.integrated,  idents =c('AA', 'Control'))
DimPlot(Data.integrated, group.by = "Cell_type", pt.size = 2)
Idents(Data.integrated) <- Data.integrated$Type
Data.integrated$Type <- Data.integrated$type

Data.integrated <- NormalizeData(object = Data.integrated,normalization.method = "LogNormalize", assay = "RNA")
Idents(Sub_cluster) <- Sub_cluster$Type
DE_AA_Control <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "DESeq2", logfc.threshold = .1, min.pct = .1, assay = 'RNA')


