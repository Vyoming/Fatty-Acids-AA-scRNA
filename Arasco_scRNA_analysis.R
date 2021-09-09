library(Rmagic)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)

AA_10x <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB10/Beyaz_10_10xgex_araasco/filtered_feature_bc_matrix.h5")

Control_10X <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB10/Beyaz_10_10xgex_control/filtered_feature_bc_matrix.h5")

AA_11x <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB11/Beyaz_11_10xgex_arasco/filtered_feature_bc_matrix.h5")

Control_11X <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB11/Beyaz_11_10xgex_control/filtered_feature_bc_matrix.h5")

AA_10x <- CreateSeuratObject(AA_10x, project = "AA_10x")
Control_10X <- CreateSeuratObject(Control_10X, project = "Control_10X")
AA_11x <- CreateSeuratObject(AA_11x, project = "AA_11x")
Control_11X <- CreateSeuratObject(Control_11X, project = "Control_11X")

d <- merge(Control_10X, y = c(AA_10x,Control_11X, AA_11x ), add.cell.ids = c("Control_10", "Arasco_10","Control_11", "Arasco_11"), project = "Arasco")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 2000 & nCount_RNA < 40000 & nFeature_RNA > 1000 & nFeature_RNA < 9000 & percent.mt < 15)
d
d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control_10X", "AA_10x","Control_11X", "AA_11x")]
for (i in 1:length(Data.list)) {
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
options (future.globals.maxSize = 4000 * 1024^5)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, d)
arasco_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)
rm(Data.anchors)
# Visulization and Clustering
arasco_obj <- RunPCA(arasco_obj, verbose = FALSE)

arasco_obj <- RunUMAP(arasco_obj, dims = 1:10)
DimPlot(arasco_obj, group.by = 'Cell_type')
levels(factor(arasco_obj@meta.data$orig.ident))
arasco_obj[["Type"]] <- Idents(arasco_obj)
new.cluster.ids <- c("Arasco", "Arasco", "Control", "Control")
names(new.cluster.ids) <- levels(arasco_obj)
arasco_obj <- RenameIdents(arasco_obj, new.cluster.ids)

arasco_obj[["Type"]] <- Idents(arasco_obj)
head(Idents(arasco_obj))
DimPlot(arasco_obj, group.by = c("Type"))


FeaturePlot(arasco_obj, features = c("Lgr5", "Ascl2", "Clu", "S100a6","Hopx", "Defa24", "Tubb5", "Fabp1","Ube2c", "Muc2","Slc39a2", "Adh1"), cols = Vyom_color_Umap)

x<-FeaturePlot(arasco_obj, features = c("Lgr5", "Clu","S100a6"))
x

arasco_obj <- FindNeighbors(arasco_obj, dims = 1:10)
arasco_obj <- FindClusters(arasco_obj, resolution = .6)#.3
DimPlot(arasco_obj, reduction = "umap", label = TRUE)

#Annotation decisions

stem.markers <- c("Lgr5" ,"Ascl2" ,"Slc12a2" ,"Axin2" ,"Olfm4" ,"Gkn3")
TA.markers <- c("Tubb5" ,"Hmgb2" ,"Stmn1" ,"H2afz1" ,"Tuba1b" ,"Hmgb11" ,"Hmgn22" ,"Ptma1" ,"Kiaa0101" ,"Tk1" ,"Cenpw" ,"Tyms" ,"Ranbp11" ,"Idh21" ,"Ran1" ,"Dtymk" ,"Nucks11" ,"Dut1" ,"Cks1b")
paneth.markers <- c("Lyz1" ,"Defa17" ,"Defa22" ,"Defa24" ,"Ang4")
Enterocyte.markers <- c("Alpi" ,"Apoa1" ,"Apoa4" ,"Fabp1" ,"Adh6a")
EP.markers <- c("Ccnb1" ,"Cdc20" ,"Cenpa" ,"Cdkn3" ,"Ccnb2" ,"Cdc25c" ,"Kif22" ,"Ube2c" ,"Sapcd2" ,"Rbp7" ,"Aurka" ,"Ccna2" ,"Cdkn2d" ,"Kif23" ,"Nek2" ,"Slc16a1")
Enteroendocrine.markers <- c("Chga" ,"Chgb" ,"Tac1" ,"Tph1" ,"Neurog3" ,"Gch1" ,"Slc39a2")
Goblet.markers <- c("Muc2" ,"Clca3" ,"Tff3" ,"Agr2" ,"Spink4" ,"Fcgbp" ,"Zg16" ,"Ccl9" ,"Atoh1")
Tuft.markers <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1")
library(readxl)
Roulis_gene_lists_for_metagenes <- read_excel("Desktop/Roulis_gene_lists_for_metagenes.xlsx")
revival.markers <- c('S100a6','Ly6a','Clu','Anxa3','Areg')

arasco.obj <- AddModuleScore(object = arasco.obj, features = list(pathways$LGR5_stem_cell_signature), name = 'Stem')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(revival.markers), name = 'revival.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(TA.markers), name = 'TA.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(paneth.markers), name = 'paneth.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(EP.markers), name = 'EP.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Goblet.markers), name = 'Goblet.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Tuft.markers), name = 'Tuft.markers')

FeaturePlot(object = arasco.obj, features = 'stem1', pt.size = .001) + scale_color_viridis(option = 'B') + 
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Stem.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Goblet.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Goblet.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'paneth.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_paneth.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Enterocyte.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Enterocyte.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Enteroendocrine.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Enteroendocrine.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'revival.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_RevStem.pdf', width=2, height=2, units="in")

DotPlot(arasco.obj, features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))

Cell_Type_Sig_Score <- data.frame(Stem=arasco_obj$stem1, Revival=arasco_obj$revival.markers1, Reserve=arasco_obj$reserve.markers1, Transit_Amplifying =arasco_obj$TA.markers1, Paneth=arasco_obj$paneth.markers1, Enterocyte=arasco_obj$Enterocyte.markers1, EP=arasco_obj$EP.markers1,Enteroendocrine= arasco_obj$Enteroendocrine.markers1,Goblet=arasco_obj$Goblet.markers1,Tuft=arasco_obj$Tuft.markers1)
FeaturePlot(object = arasco_obj, features = c('Goblet.markers1','paneth.markers1','S100a6'), cols = Vyom_color_Umap)

# module score distribution
modulescores <- Cell_Type_Sig_Score %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

Idents(arasco_obj) <- arasco_obj$Type

z_Stem_cluster<-WhichCells(object = arasco_obj, expression = stem1 > 1)
z_RevSC_cluster<-WhichCells(object = arasco_obj,idents = c('Arasco') ,expression = revival.markers1 > 1, )
z_reserve_cluster<-WhichCells(object = arasco_obj,idents = c('Arasco') , expression = reserve.markers1 > 1)
z_TA_cluster<-WhichCells(object = arasco_obj, expression = TA.markers1 > .5)
z_paneth_cluster<-WhichCells(object = arasco_obj, expression = paneth.markers1 > 5)
z_Enterocyte_cluster<-WhichCells(object = arasco_obj, expression = Enterocyte.markers1 > 1)
z_EP_cluster<-WhichCells(object = arasco_obj, expression = EP.markers1 > .5)
z_Enteroendocrine_cluster<-WhichCells(object = arasco_obj, expression = Enteroendocrine.markers1 > 11)
z_Goblet_cluster<-WhichCells(object = arasco_obj, expression = Goblet.markers1 > 1)
z_Tuft_cluster<-WhichCells(object = arasco_obj, expression = Tuft.markers1 > 1)

table(arasco_obj$Cell_type)

arasco_obj <- FindClusters(arasco_obj, resolution = .6)#.3
DimPlot(arasco_obj, reduction = "umap", label = TRUE)

DimPlot(arasco_obj, label=T,group.by = 'Cell_type', cells.highlight= list( z_Stem_cluster, z_RevSC_cluster,  z_TA_cluster, z_paneth_cluster, z_Enterocyte_cluster, z_EP_cluster, z_Enteroendocrine_cluster, z_Goblet_cluster, z_Tuft_cluster),  cols.highlight = c("darkblue", "darkred", 'green','lightgreen', 'darkgreen', 'yellow', 'pink', 'purple', 'orange','lightblue'),cols= "grey")
DimPlot(arasco_obj, label=T,group.by = 'Cell_type', cells.highlight= list(z_paneth_cluster),  cols.highlight = c("darkblue", "darkred"),cols= "grey")
VlnPlot(arasco_obj, features = c('paneth.markers1','Goblet.markers1'), pt.size = 0)

z_RevSC_clusterx<-WhichCells(object = arasco_obj ,expression = revival.markers1 > 1)
x_rev <- WhichCells(object = arasco_obj,idents = c('x'), z_RevSC_clusterx )
write.csv(x_rev)
#______________________________Clustering________________________________________#
new.cluster.ids <- c('Enterocyte (Proximal)', 'Stem 1', 'Goblet', 'Enterocyte Progenitor', 'Enterocyte (Proximal)', 'Enterocyte (Proximal)', 'Stem 1', 'Transit Amplifying','Enterocyte Progenitor','Transit Amplifying', 'Stem 1', 'Enterocyte (Proximal)', 'Enterocyte (Distal)', 'Enteroendocrine', 'Enterocyte (Distal)', 'Enterocyte Progenitor', 'Tuft', 'CD8+ intraepithelial cells', 'Macrophage', 'Goblet', 'Enterocyte (Proximal)', 'x', 'Paneth')
arasco_obj[["Cell_type"]] <- Idents(arasco_obj)
names(new.cluster.ids) <- levels(arasco_obj)
arasco_obj <- RenameIdents(arasco_obj, new.cluster.ids)
arasco_obj[["Cell_type"]] <- Idents(arasco_obj)
DimPlot(arasco_obj, reduction = "umap", group.by= 'Cell_type')

Idents(arasco_obj) <- arasco_obj$Cell_type
arasco_obj$Cell_Type <- as.character(Idents(arasco_obj))
stem2_cells1 <- WhichCells(object = arasco_obj,idents = c('x'), z_RevSC_cluster )
stem2_cells <- c(stem2_cells1, z_Enteroendocrine_cluster)
arasco_obj$Cell_Type[stem2_cells] <- paste('Stem 2')
DimPlot(arasco_obj, reduction = "umap", group.by= 'Cell_Type')
arasco_obj$Cell_type


levels(factor(arasco_obj@meta.data$Type))
arasco_obj@meta.data$type <- factor(arasco_obj@meta.data$Type, levels = c("Control", "Arasco")) 
Idents(arasco_obj) <- arasco_obj$Cell_Type
prop.table(x = table(Idents(arasco_obj), arasco_obj$type), margin = 2)
tail(arasco_obj$type)
arasco.obj <- subset(arasco_obj, idents = c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft'))
DimPlot(arasco.obj)
my_levels <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
arasco.obj$Cell_type <- factor(x = arasco.obj$Cell_Type, levels = my_levels)
Idents(arasco.obj) <- arasco.obj$Cell_type
Proportion_Arasco<- prop.table(x = table(arasco.obj$Cell_type, arasco.obj$type), margin = 2)
DimPlot(arasco.obj, group.by = "Cell_type", label = FALSE, pt.size=1, label.size = 3)
arasco.obj@meta.data$type <- factor(arasco.obj@meta.data$Type, levels = c("Control", "Arasco")) 

arasco.obj <- NormalizeData(object = arasco.obj,normalization.method = "LogNormalize", assay = "RNA")
Arasco_DE <- FindMarkers(arasco.obj, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = .1, min.pct = .05, assay = 'RNA')

# final Plots
DefaultAssay(arasco.obj) <- "RNA"
arasco.obj <- magic(arasco.obj)

All_Genes <- arasco.obj@assays$RNA@data@Dimnames[[1]]
Fetal_Genes <- intersect(All_Genes, Fetal$Gene_Name)
Radiation_Genes <- intersect(All_Genes, Radiation$Gene_Name)
Homeostatic_Genes <- intersect(All_Genes, Homeostatic$Gene_Name)
Regeneration_Genes <- intersect(All_Genes, Regeneration$Gene_Name)
Granuloma_Genes <- intersect(All_Genes, Granuloma$Gene_Name)
CREB_Targets_Genes <- intersect(All_Genes, CREB_Targets$Gene_Name)
BCatenin_Enhanced_Genes <- intersect(All_Genes, BCatenin_Enhanced$Gene_Name)
YAP_Targets_Genes <- intersect(All_Genes, YAP_Targets$Gene_Name)
ECM_Induced_Genes <- intersect(All_Genes, ECM_Induced$Gene_Name)
ASCL_Targets_Genes <- intersect(All_Genes, ASCL_Targets$Gene_Name)

Scores <- list(Fetal = Fetal_Genes, Radiation = Radiation_Genes, Homeostatic = Homeostatic_Genes, Regeneration = Regeneration_Genes, Granuloma = Granuloma_Genes, ECM_Induced = ECM_Induced_Genes, CREB_Targets = CREB_Targets_Genes, BCatenin_Enhanced = BCatenin_Enhanced_Genes, YAP_Targets = YAP_Targets_Genes, ASCL_Targets = ASCL_Targets_Genes)

writeGmtPathways(Scores, 'Beyaz_AA_final.gmt')

pathways = getPathways()
getPathways = function(genesOrth = NULL)
{
  pathwaysM = c(gmtPathways("~/analysis/scRNA_Intestine/Beyaz_AA_final.gmt"))
  return(pathwaysM)
}
BiocManager::install("fgsea")
mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[pathways[['Fetal']], ], na.rm = TRUE), dist = 'norm')


mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[Radiation_Genes, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$Radiation_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[Fetal_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$Fetal_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[Granuloma_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$Granuloma_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[Regeneration_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$Regeneration_Score <- mean.exp
}
mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[Homeostatic_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$Homeostatic_Stem_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[CREB_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$CREB_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[BCatenin_Enhanced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$BCatenin_Enhanced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[YAP_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$YAP_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[ECM_Induced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$ECM_Induced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = arasco.obj@assays$MAGIC_RNA@data[ASCL_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = arasco.obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  arasco.obj@meta.data$ASCL_Targets_Score <- mean.exp
}

Radiation_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Radiation_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')
Fetal_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Fetal_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Granuloma_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Granuloma_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Regen_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Regeneration_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Homeostatic_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Homeostatic_Stem_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
YAP_graph <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('YAP_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
wnt_graph <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('BCatenin_Enhanced_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Creb_graph <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('CREB_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
ASCL_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('ASCL_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
ECM_grpah <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('ECM_Induced_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')

Score_final <- plot_grid(Radiation_grpah, Fetal_grpah, Granuloma_grpah, Regen_grpah, ECM_grpah, ncol = 5, labels = c('A','B','C','D','E'), label_size = 12)
Score_final

ggsave(file = 'Radiation_grpah.svg', plot=Radiation_grpah, width=10, height=10)
ggsave(file = 'fetal.svg', plot=Fetal_grpah, width=10, height=10)
ggsave(file = 'granuloma.svg', plot=Granuloma_grpah, width=10, height=10)
ggsave(file = 'regen.svg', plot=Regen_grpah, width=10, height=10)
ggsave(file = 'Homeostatic_Stem_Score.svg', plot=Homeostatic_grpah, width=10, height=10)


gene.level.list <- c('Radiation_Score', 'Fetal_Score', 'Granuloma_Score', 'Regeneration_Score', 'Homeostatic_Score')
All <- VlnPlot(arasco.obj, group.by = 'Cell_type' ,  split.by = "type", features = gene.level.list, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'),log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
ggsave(file = paste0('Arasco_Score_vln.pdf'), plot=All,  width=4, height=3, units="in")


Lgr5 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Lgr5'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Ascl2 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ascl2'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
S100a6 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('S100a6'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('')
Ly6a <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ly6a'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Cd55 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Cd55'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ptger4'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')

ggsave(file = 'Lgr5.svg', plot=Lgr5, width=10, height=10)
ggsave(file = 'Ascl2.svg', plot=Ascl2, width=10, height=10)
ggsave(file = 'S100a6.svg', plot=S100a6, width=10, height=10)
ggsave(file = 'Ly6a.svg', plot=Ly6a, width=10, height=10)
ggsave(file = 'Cd55.svg', plot=Cd55, width=10, height=10)

PGE2_genes <-  c('Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')
Gene_de$type
for(i in PGE2_genes){
  Gene_plot_pge <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c(i), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('') 
  ggsave(file = paste0('ARAsco_', i, '.svg'), plot=Gene_plot_pge, width=10, height=10)
}

Gene_final <- plot_grid(Lgr5, Ascl2, S100a6, Ly6a, Cd55, ncol = 5, labels = c('A','B','C','D','E'), label_size = 12)
Gene_final

gene.level.list <- c('Lgr5', 'Ascl2', 'Olfm4', 'Clu', 'Ly6a', 'S100a6', 'Bmi1', 'Hopx', 'Lrig1', 'Mex3a')
VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 5, split.by = "type", features = gene.level.list, pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE)

check <- c('Lgr5','S100a6','Hopx',"Ccnb1","Fabp1", "Agr2", "Defa17", "Trpm5", 'Chgb', 'Tk1', 'Stmn1','Atoh1'  )
Lgr5 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Lgr5'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr51 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('S100a6'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr52 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Hopx'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr53 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Ccnb1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr54 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Fabp1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr55 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Agr2'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr56 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Defa17'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr57 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Trpm5'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr58 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Chgb'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr59 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Tk1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr510 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Stmn1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr511 <- VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Atoh1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
alt_final <- plot_grid(Lgr5,Lgr51, Lgr52, Lgr53, Lgr54, Lgr55, Lgr56, Lgr57, Lgr58, Lgr59, Lgr510, Lgr511, ncol = 4, labels = c('A','B','C','D','E','F','G','H','I','J', 'K','L'), label_size = 12)
 

DimPlot(arasco.obj, pt.size = 2, group.by = 'Cell_type')

#proportion analysis
Prop_table<- prop.table(x = table(arasco.obj$Cell_type, arasco.obj$type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
#Prop_Table1$Freq <- 1/log10(Prop_Table$Freq)*-1
Prop_Table1$Freq[2] <- c(0)
my_levels <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','Arasco'))

Prop_Table1 %>% anova_test(dv = Freq, wid = rownames(Prop_Table1), between = Var1)
glm.Prop_Table1 <- glm(Freq ~ Var1 + Var2, family = poisson(), data = Prop_Table1)
anova.test <- anova(glm.Prop_Table1, test = "Chisq")
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="ARAsco_prop.svg", plot=plot, width=10, height=8)

prop_sig<- table(arasco.obj$Cell_type, arasco.obj$orig.ident)
Prop_Sig_Table <- as.data.frame.matrix(prop_sig)
sum(Prop_Sig_Table$AA_10x)
num <- length(rownames(Prop_Sig_Table))

for(i in 1:num){
  dat = data.frame(sampleID = c(colnames(Prop_Sig_Table)), current_cluster=c(Prop_Sig_Table[i,]$AA_10x, Prop_Sig_Table[i,]$AA_11x, Prop_Sig_Table[i,]$Control_10X, Prop_Sig_Table[i,]$Control_11X), all_other_clusters=c(sum(Prop_Sig_Table$AA_10x) - Prop_Sig_Table[i,]$AA_10x, sum(Prop_Sig_Table$AA_11x) - Prop_Sig_Table[i,]$AA_11x, sum(Prop_Sig_Table$Control_10X) - Prop_Sig_Table[i,]$Control_10X, sum(Prop_Sig_Table$Control_11X) - Prop_Sig_Table[i,]$Control_11X), Treatment = c("Arasco", "Arasco", "Control", "Control"))
  glm_output = glm(data = dat, cbind(current_cluster, all_other_clusters) ~ Treatment, family = "binomial")
  print(coef(summary(glm_output))[,'Pr(>|z|)'])
  }


#markers
Idents(arasco.obj) <- arasco.obj$Cell_type
Acluster.markers1 <- FindAllMarkers(arasco.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'SCT')
top10 <- Acluster.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
check_heatmap <- DoHeatmap(arasco.obj, features = top10$gene, group.by = 'Cell_type', assay = 'SCT',size = 2 , angle = 60) + NoLegend() + scale_fill_viridis(option="inferno")#top10$Gene
ggsave(file = 'ARAsco_top10_cluster_Heatmap.svg', plot=check_heatmap, width=10, height=15)
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
check_heatmap <- DoHeatmap(arasco.obj, features = DotPlot_Sig, group.by = 'Cell_type', assay = 'SCT',size = 2 , angle = 60) + NoLegend() + scale_fill_viridis(option="inferno")#top10$Gene
ggsave(file = 'ARAsco_haber_cluster_Heatmap.svg', plot=check_heatmap, width=10, height=10)

write.csv(Acluster.markers1,'ARAsco_Sigs_Per_Clust.csv')
my_levels <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine','Tuft', 'Goblet', 'Paneth')
arasco.obj$Cell_type <- factor(x = arasco.obj$Cell_type, levels = my_levels)
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Pou2f3","Avil","Tuba1a","Adh1","Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4") 
DotPlot(arasco.obj, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text= element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file= 'Arasco_Haber_check_dotplot.pdf', width=4.5, height=2.6, units="in")


stem.markers <- c("Lgr5" ,"Ascl2" ,"Slc12a2" ,"Axin2" ,"Olfm4" ,"Gkn3")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(stem.markers), name = 'stem')
FeaturePlot(object = arasco.obj, features = 'stem1') + ggtitle("Stem 1 Signature") & scale_color_viridis_c()

revival.markers <- c(Roulis_gene_lists_for_metagenes$Gene)
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(revival.markers), name = 'revival.markers' ,assay =  'RNA')
FeaturePlot(object = arasco.obj, features = 'revival.markers1') + ggtitle("Stem 2 Signature") & scale_color_viridis_c()

TA.markers <- c("Tubb5" ,"Hmgb2" ,"Stmn1" ,"H2afz1" ,"Tuba1b" ,"Hmgb11" ,"Hmgn22" ,"Ptma1" ,"Kiaa0101" ,"Tk1" ,"Cenpw" ,"Tyms" ,"Ranbp11" ,"Idh21" ,"Ran1" ,"Dtymk" ,"Nucks11" ,"Dut1" ,"Cks1b")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(TA.markers), name = 'TA.markers')
FeaturePlot(object = arasco.obj, features = 'TA.markers1') + ggtitle("Transit Amplifying Signature") & scale_color_viridis_c()

paneth.markers <- c("Lyz1" ,"Defa17" ,"Defa22" ,"Defa24" ,"Ang4")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(paneth.markers), name = 'paneth.markers')
FeaturePlot(object = arasco.obj, features = 'paneth.markers1') + ggtitle("Paneth Signature") & scale_color_viridis_c()

Enterocyte.markers <- c("Alpi" ,"Apoa1" ,"Apoa4" ,"Fabp1" ,"Adh6a")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
FeaturePlot(object = arasco.obj, features = 'Enterocyte.markers1') + ggtitle("Enterocyte Signature") & scale_color_viridis_c()

EP.markers <- c("Ccnb1" ,"Cdc20" ,"Cenpa" ,"Cdkn3" ,"Ccnb2" ,"Cdc25c" ,"Kif22" ,"Ube2c" ,"Sapcd2" ,"Rbp7" ,"Aurka" ,"Ccna2" ,"Cdkn2d" ,"Kif23" ,"Nek2" ,"Slc16a1")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(EP.markers), name = 'EP.markers')
FeaturePlot(object = arasco.obj, features = 'EP.markers1')  + ggtitle("Enterocyte Progenitor Signature") & scale_color_viridis_c()

Enteroendocrine.markers <- c("Chga" ,"Chgb" ,"Tac1" ,"Tph1" ,"Neurog3" ,"Gch1" ,"Slc39a2")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
FeaturePlot(object = arasco.obj, features = 'Enteroendocrine.markers1') + ggtitle("Enteroendocrine Signature") & scale_color_viridis_c()

Goblet.markers <- c("Muc2" ,"Clca3" ,"Tff3" ,"Agr2" ,"Spink4" ,"Fcgbp" ,"Zg16" ,"Ccl9" ,"Atoh1")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Goblet.markers), name = 'Goblet.markers')
FeaturePlot(object = arasco.obj, features = 'Goblet.markers1') + ggtitle("Goblet Signature") & scale_color_viridis_c()

Tuft.markers <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1")
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Tuft.markers), name = 'Tuft.markers')
FeaturePlot(object = arasco.obj, features = 'Tuft.markers1') + ggtitle("Tuft Signature") & scale_color_viridis_c()


# Pseudotime Analysis
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(arasco.obj) <- arasco.obj$Cell_type
pseudo_stem <- WhichCells(object = arasco.obj, idents = 'Stem 1')

arasco.obj1 <- arasco.obj
cds <- as.cell_data_set(arasco.obj1)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(arasco.obj1[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)

cds <- order_cells(cds, reduction_method = "UMAP")
saveRDS(cds, file = "ARAsco_pseudotime.rds")
plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE)
ggsave(file= 'Arasco_Pseudotime_umap.svg',  width=10, height=10)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file= 'Arasco_Pseudotime_umap.svg',  width=10, height=10)



cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.000000001))

plot_cells(cds, genes = c("S100a6", "Ly6a", "Olfm4", "Ascl2"),
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE)
write.csv(cds_pr_test_res, 'Pseudotime_prominent.csv')

Psuedo_genes <- c("S100a6", "Ly6a", "Olfm4", "Ascl2")
Gene_cds <- cds[rowData(cds)$gene_short_name %in% Psuedo_genes,]

heatmap <- plot_pseudotime_heatmap(Gene_cds, cluster_rows = FALSE, show_rownames = TRUE, return_heatmap = T)
heatmap
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("S100a6", "Ly6a", "Olfm4", "Ascl2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="Cell_type",
                         min_expr=0.1)
cdscore = cds@colData$Regeneration_Score
new_df <- data.frame(pseudotime=pseudotime(cds),score= cds@colData$Regeneration_Score)#be sure the oxfos score is in the same order as the cells in colData

ggplot(new_df, aes(x = pseudotime,y = score)) + geom_point() + geom_smooth(method = 'lm')


plot_genes_violin(Gene_cds, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
cds$Gm37381

#visualize Seurat object using the pseudotime variable
arasco.obj <- AddMetaData(
  object = arasco.obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

Idents(arasco.obj) <- arasco.obj$Cell_type
colScale = viridis::magma(40)
Pseudotime_umap <- FeaturePlot(arasco.obj, c("Pseudotime"), pt.size = 2, label = FALSE, cols = colScale) 
Pseudotime_umap + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime") + scale_fill_viridis(option="inferno")

ggsave(file= 'Arasco_Pseudotime_umap.svg', plot= Pseudotime_umap, width=10, height=10)

Idents(arasco.obj) <- arasco.obj$type
ARA_Control <- subset(arasco.obj,  idents =c('Control'))
ARA_ARA <- subset(arasco.obj,  idents =c('Arasco'))

arasco.obj@meta.data$quantile <- as.factor(ntile(arasco.obj$Pseudotime, 10))
Idents(arasco.obj) <- arasco.obj$quantile
{
Sub_cluster <- subset(arasco.obj, idents = c('10'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_10 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_10 <- ARA_psuedo_gene_10[order(rownames(ARA_psuedo_gene_10)),]
ARA_psuedo_gene_10 <- subset(ARA_psuedo_gene_10, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('1'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_1 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_1 <- ARA_psuedo_gene_1[order(rownames(ARA_psuedo_gene_1)),]
ARA_psuedo_gene_1 <- subset(ARA_psuedo_gene_1, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('2'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_2 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_2 <- ARA_psuedo_gene_2[order(rownames(ARA_psuedo_gene_2)),]
ARA_psuedo_gene_2 <- subset(ARA_psuedo_gene_2, select = c('p_val_adj'))
Sub_cluster <- subset(arasco.obj, idents = c('3'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_3 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_3 <- ARA_psuedo_gene_3[order(rownames(ARA_psuedo_gene_3)),]
ARA_psuedo_gene_3 <- subset(ARA_psuedo_gene_3, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('4'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_4 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_4 <- ARA_psuedo_gene_4[order(rownames(ARA_psuedo_gene_4)),]
ARA_psuedo_gene_4 <- subset(ARA_psuedo_gene_4, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('5'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_5 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_5 <- ARA_psuedo_gene_5[order(rownames(ARA_psuedo_gene_5)),]
ARA_psuedo_gene_5 <- subset(ARA_psuedo_gene_5, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('6'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_6 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_6 <- ARA_psuedo_gene_6[order(rownames(ARA_psuedo_gene_6)),]
ARA_psuedo_gene_6 <- subset(ARA_psuedo_gene_6, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('7'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_7 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_7 <- ARA_psuedo_gene_7[order(rownames(ARA_psuedo_gene_7)),]
ARA_psuedo_gene_7 <- subset(ARA_psuedo_gene_7, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('8'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_8 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_8 <- ARA_psuedo_gene_8[order(rownames(ARA_psuedo_gene_8)),]
ARA_psuedo_gene_8 <- subset(ARA_psuedo_gene_8, select = c('p_val_adj'))

Sub_cluster <- subset(arasco.obj, idents = c('9'))
Idents(Sub_cluster) <- Sub_cluster$type
ARA_psuedo_gene_9 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
ARA_psuedo_gene_9 <- ARA_psuedo_gene_9[order(rownames(ARA_psuedo_gene_9)),]
ARA_psuedo_gene_9 <- subset(ARA_psuedo_gene_9, select = c('p_val_adj'))

S100a6 <- c(ARA_psuedo_gene_1['S100a6',], ARA_psuedo_gene_2['S100a6',], ARA_psuedo_gene_3['S100a6',], ARA_psuedo_gene_4['S100a6',], ARA_psuedo_gene_5['S100a6',], ARA_psuedo_gene_6['S100a6',], ARA_psuedo_gene_7['S100a6',], ARA_psuedo_gene_8['S100a6',], ARA_psuedo_gene_9['S100a6',], ARA_psuedo_gene_10['S100a6',])
Ly6a <- c(ARA_psuedo_gene_1['Ly6a',], ARA_psuedo_gene_2['Ly6a',], ARA_psuedo_gene_3['Ly6a',], ARA_psuedo_gene_4['Ly6a',], ARA_psuedo_gene_5['Ly6a',], ARA_psuedo_gene_6['Ly6a',], ARA_psuedo_gene_7['Ly6a',], ARA_psuedo_gene_8['Ly6a',], ARA_psuedo_gene_9['Ly6a',], ARA_psuedo_gene_10['Ly6a',])
Lgr5 <- c(ARA_psuedo_gene_1['Lgr5',], ARA_psuedo_gene_2['Lgr5',], ARA_psuedo_gene_3['Lgr5',], ARA_psuedo_gene_4['Lgr5',], ARA_psuedo_gene_5['Lgr5',], ARA_psuedo_gene_6['Lgr5',], ARA_psuedo_gene_7['Lgr5',], ARA_psuedo_gene_8['Lgr5',], ARA_psuedo_gene_9['Lgr5',], ARA_psuedo_gene_10['Lgr5',])
Cd55 <- c(ARA_psuedo_gene_1['Cd55',], ARA_psuedo_gene_2['Cd55',], ARA_psuedo_gene_3['Cd55',], ARA_psuedo_gene_4['Cd55',], ARA_psuedo_gene_5['Cd55',], ARA_psuedo_gene_6['Cd55',], ARA_psuedo_gene_7['Cd55',], ARA_psuedo_gene_8['Cd55',], ARA_psuedo_gene_9['Cd55',], ARA_psuedo_gene_10['Cd55',])
Ascl2 <- c(ARA_psuedo_gene_1['Ascl2',], ARA_psuedo_gene_2['Ascl2',], ARA_psuedo_gene_3['Ascl2',], ARA_psuedo_gene_4['Ascl2',], ARA_psuedo_gene_5['Ascl2',], ARA_psuedo_gene_6['Ascl2',], ARA_psuedo_gene_7['Ascl2',], ARA_psuedo_gene_8['Ascl2',], ARA_psuedo_gene_9['Ascl2',], ARA_psuedo_gene_10['Ascl2',])
pseudo_PVal <- data.frame(S100a6, Ly6a, Lgr5, Cd55, Ascl2)
write.csv(pseudo_PVal, 'Arasco_Psuedo_gene_signficance.csv')
}
#make plots
pseudo_PVal <- read.csv('./data/AA_ARA_Intestine_project/Arasco_Psuedo_gene_signficance.csv')
colnames(pseudo_PVal)<- c('num' ,"S100a6", "Ly6a",  "Lgr5",  "Cd55", "Ascl2")

ARA_rownames <- rownames(ARA_ARA@meta.data)
Control_rownames <- rownames(ARA_Control@meta.data)

Cell_type_frame <- data.frame(arasco.obj$type,arasco.obj$orig.ident )
all <-  cbind(type = Cell_type_frame, Pseudotime = arasco.obj$Pseudotime)
colnames(all) <- c("Type",'Replicate',"Pseudotime")

all <- all %>% mutate(quantile = ntile(all$Pseudotime, 10))
quart <- data.frame(var1 = tapply(all$Pseudotime, all$quantile, mean), var2 = tapply(all$Pseudotime, all$quantile, max))

features_pseudo_scores <- c("Fetal_Score", "Granuloma_Score", "Radiation_Score", "Homeostatic_Score", "Regeneration_Score")
for(i in features_pseudo_scores) {
  score_iter = i
  quart$maxi = max(arasco.obj@meta.data[,i])
  ARAsco <- subset(all, Type=='Arasco')
  Control <- subset(all, Type=='Control')
  
  ARAsco <-  cbind(ARAsco, ARA_ARA@meta.data[,i])
  Control <-  cbind(Control, ARA_Control@meta.data[,i])
  
  colnames(ARAsco) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(ARAsco, quantile== i)
    Cinqunat  <- subset(Control, quantile== i)

    p.vals.prop = c(p.vals.prop, wilcox.test(AAinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  Pseudotimedot <- ggplot() +
    # Points
    geom_jitter(data=Control, aes(x=Pseudotime, y=score), size=.25, alpha=0.1, width=1,stroke = .3, colour = "#1b9e77" ) + 
    geom_jitter(data=ARAsco, aes(x=Pseudotime, y=score), size=.25, alpha=0.1, width=1,stroke = .3, colour = "#d95f02") +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=1) + 
    geom_smooth(data=ARAsco, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=1) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"), position=position_dodge2(0.75), vjust=.5, inherit.aes = FALSE)
  nameoooo <-  'Arasco_'
  name1 <- score_iter
  Name2 <- '.svg'
  File_name <-  paste(nameoooo, name1, Name2, sep = "")
  ggsave(file= File_name, plot=Pseudotimedot, width=5, height=2.5)
}

features_pseudo_genes <- c( "S100a6", "Ascl2", "Lgr5", "Ly6a", "Cd55")
for(i in features_pseudo_genes) {
  score_iter = i
  quart$maxi = max(arasco.obj@assays$MAGIC_RNA@data[i, ])
  ARAsco <- subset(all, Type=='Arasco')
  Control <- subset(all, Type=='Control')
  
  ARAsco <-  cbind(ARAsco, ARA_ARA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, ARA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(ARAsco) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(ARAsco, quantile== i)
    Cinqunat  <- subset(Control, quantile== i)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(AAinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.05] = "*"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.01] = "**"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  Pseudotimedot <- ggplot() +
    # Points
    ggrastr::rasterise(geom_jitter(data=Control, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#1b9e77" ), dpi = 300) + 
    ggrastr::rasterise(geom_jitter(data=ARAsco, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=ARAsco, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) + theme_cem +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + theme(text = element_text(size=5),axis.text = element_text(size = 6) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 0))
  nameoooo <-  'Arasco_'
  name1 <- score_iter
  Name2 <- '.pdf'
  File_name <-  paste(nameoooo, name1, Name2, sep = "")
  ggsave(file= File_name, plot=Pseudotimedot, width = 2.2, height = 1, units="in")
}

{
  #main figure pseudotime plot
  #Lgr5
  i = 'Lgr5'
  score_iter = i
  quart$maxi = max(arasco.obj@assays$MAGIC_RNA@data[i, ])
  ARAsco <- subset(all, Type=='Arasco')
  Control <- subset(all, Type=='Control')
  
  ARAsco <-  cbind(ARAsco, ARA_ARA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, ARA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(ARAsco) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(ARAsco, quantile== i)
    Cinqunat  <- subset(Control, quantile== i)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(AAinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.05] = "*"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.01] = "**"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  Pseudotimedot_Lgr5 <- ggplot() +
    # Points
    ggrastr::rasterise(geom_jitter(data=Control, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#1b9e77" ), dpi = 300) + 
    ggrastr::rasterise(geom_jitter(data=ARAsco, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=ARAsco, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  #Ascl2
  i = 'Ascl2'
  score_iter = i
  quart$maxi = max(arasco.obj@assays$MAGIC_RNA@data[i, ])
  ARAsco <- subset(all, Type=='Arasco')
  Control <- subset(all, Type=='Control')
  
  ARAsco <-  cbind(ARAsco, ARA_ARA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, ARA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(ARAsco) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(ARAsco, quantile== i)
    Cinqunat  <- subset(Control, quantile== i)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(AAinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.05] = "*"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.01] = "**"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  Pseudotimedot_Ascl2 <- ggplot() +
    # Points
    ggrastr::rasterise(geom_jitter(data=Control, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#1b9e77" ), dpi = 300) + 
    ggrastr::rasterise(geom_jitter(data=ARAsco, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=ARAsco, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  #Ly6a
  i = 'Ascl2'
  score_iter = i
  quart$maxi = max(arasco.obj@assays$MAGIC_RNA@data[i, ])
  ARAsco <- subset(all, Type=='Arasco')
  Control <- subset(all, Type=='Control')
  
  ARAsco <-  cbind(ARAsco, ARA_ARA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, ARA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(ARAsco) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(ARAsco, quantile== i)
    Cinqunat  <- subset(Control, quantile== i)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(AAinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.05] = "*"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.01] = "**"
  quart$p_asterik[pseudo_PVal[[score_iter]] < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  Pseudotimedot_Ly6a <- ggplot() +
    # Points
    ggrastr::rasterise(geom_jitter(data=Control, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#1b9e77" ), dpi = 300) + 
    ggrastr::rasterise(geom_jitter(data=ARAsco, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=ARAsco, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  plot_grid(Density_diff_stem_1, Pseudotimedot_Lgr5, Pseudotimedot_Ascl2, labels = c('Enterocyte Progenitor','Lgr5', 'Ascl2'), label_size = 6, ncol = 1, align = "v", hjust = -.005) #'All Cells', 'Stem 1', 'Stem 2'
  ggsave(paste0("Arasco_pseudo_gene.pdf"), width = 2.2, height = 3.1, units="in")
}

#DE analysis
Idents(arasco.obj) <- arasco.obj$type
arasco.obj1 <- arasco.obj
arasco.obj <- arasco.obj1
arasco.obj1 <- NormalizeData(object = arasco.obj1, normalization.method = "LogNormalize", assay = "RNA")
all.genes <- rownames(arasco.obj1)
arasco.obj1 <- ScaleData(arasco.obj1, features = all.genes, assay = 'RNA')
arasco.obj1@assays$RNA@scale.data
DE_valX <- FindMarkers(arasco.obj, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'RNA')

DE_valX<- DE_valX[!grepl("mt-", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Rpl", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Rps", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Atp", rownames(DE_valX)),]
DE_valX <- DE_valX[order(-DE_valX$avg_logFC),]
DE_valX['S100a6',]
EnhancedVolcano(DE_valX, lab = rownames(DE_valX), x = 'avg_logFC', y = 'p_val_adj', title = 'Arasco vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-1.2, 1.2), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_valX, 'ARAsco_DE_All.csv')
cluster.ids <- c("Stem 1", "Stem 2", "Transit Amplifying")
Idents(arasco.obj) <- arasco.obj$Cell_type
Sub_cluster <- subset(arasco.obj,  idents = cluster.ids)
Idents(Sub_cluster) = Sub_cluster$type
Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
DE_val1 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control",  slot='scale.data', test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'RNA')
DE_val1<- DE_val1[!grepl("mt-", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Mt", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Rpl", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Rps", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Atp", rownames(DE_val1)),]
DE_val1 <- DE_val1[order(-DE_val1$avg_diff),]
DE_val1['S100a6',]
EnhancedVolcano(DE_val1, lab = rownames(DE_val1), x = 'avg_diff', y = 'p_val_adj', title = 'Arasco vs Control', pCutoff = 10e-25, FCcutoff = 0.5, xlim = c(-5, 5), ylim = c(0,310), subtitle = 'Pluripotent Cells' )
write.csv(DE_val1, 'ARAsco_DE_pluripotent.csv')

#DE check

# STem1v2
Idents(arasco.obj) <- arasco.obj$Cell_type
Stem2v1 <- FindMarkers(arasco.obj, ident.1 = "Stem 2", ident.2 = "Stem 1", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

Stem2v1<- Stem2v1[!grepl("mt-", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Rpl", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Rps", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Atp", rownames(Stem2v1)),]
Stem2v1 <- Stem2v1[order(-Stem2v1$avg_logFC),]
volc <- EnhancedVolcano(Stem2v1, lab = rownames(Stem2v1), x = 'avg_logFC', y = 'p_val_adj', title = 'Stem 2 vs Stem 1', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-4, 4), ylim = c(0,330),  shade = c(''), shadeAlpha = 0 )
ggsave(file = 'DE_Stem2v1_ARAsco.svg', plot=volc, width=10, height=10)
write.csv(Stem2v1, 'DE_Stem2v1_ARAsco.csv')

#ara v c
AA_induced_Genes <- read_csv("AA_induced_Genes.csv")$gene_name
arasco.obj <- AddModuleScore(object = arasco.obj, features = 'Ptger4', name = 'ptger4')
VlnPlot(arasco.obj, features = 'ptger41',group.by = 'type')
ggsave('plot.svg')
arasco.obj1 <-  arasco.obj
arasco.obj$AA_induced <- arasco.obj$Cell_type
new.cluster.ids <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
Idents(arasco.obj) <- arasco.obj$AA_induced
current.cluster.ids <- c('Not Treated','AA Induced', 'Not Treated', 'Not Treated','Not Treated','Not Treated', 'Not Treated', 'Not Treated', 'Not Treated','Not Treated')
arasco.obj$AA_induced <- plyr::mapvalues(x = arasco.obj$AA_induced, from = new.cluster.ids, to = current.cluster.ids)
DimPlot(arasco.obj$AA_induced)


AA_induced_cells <-WhichCells(object = arasco.obj, expression = ptger41 > .01)
Idents(arasco.obj) <- arasco.obj$AA_induced
arasco.obj$AA_induced[AA_induced_cells] <- paste('AA Induced')
table(arasco.obj$AA_induced, arasco.obj$type)

Idents(arasco.obj) = arasco.obj$AA_induced
AA_induced_1 <- subset(arasco.obj,  idents =c('AA Induced'))
AA_induced_not <- subset(arasco.obj,  idents =c('Not Treated'))
Idents(AA_induced_1) <- AA_induced_1$Type
Idents(AA_induced_not) <- AA_induced_not$Type
AA_induced_1@assays$SCT@scale.data
DE_AA_induced_vs_All <- FindMarkers(AA_induced_1, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("mt-", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Rpl", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Rps", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Atp", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All <- DE_AA_induced_vs_All[order(-DE_AA_induced_vs_All$avg_logFC),]

volc <- EnhancedVolcano(DE_AA_induced_vs_All, lab = rownames(DE_AA_induced_vs_All), x = 'avg_logFC', y = 'p_val_adj', title = 'Arasco vs Control', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-2, 2), ylim = c(0,100), subtitle = 'Ptger4 Expressing Cells',  shade = c('Ly6d','Ly6a','Lyz2', 'Fabp1', 'Tyrobp', 'S100a6', 'Cd74', 'Rac2'), shadeAlpha = 0 )
ggsave(file = 'DE_PTGER4_expressing_ARAsco.svg', plot=volc, width=10, height=10)
write.csv(DE_AA_induced_vs_All, 'DE_PTGER4_expressing_ARAsco.csv')


DE_AA_not_induced_vs_All <- FindMarkers(AA_induced_not, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST",  logfc.threshold = .05, min.pct = 0, assay = 'RNA')

DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("mt-", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Rpl", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Rps", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Atp", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All <- DE_AA_not_induced_vs_All[order(-DE_AA_not_induced_vs_All$avg_logFC),]

volc <- EnhancedVolcano(DE_AA_not_induced_vs_All, lab = rownames(DE_AA_not_induced_vs_All), x = 'avg_logFC', y = 'p_val_adj', title = 'Arasco vs Control', pCutoff = 10e-25, FCcutoff = .15, xlim = c(-2, 2), ylim = c(0,275), subtitle = 'Ptger4 Null Cells' )
ggsave(file = 'DE_PTGER4_NULL_ARAsco.svg', plot=volc, width=10, height=10)
write.csv(DE_AA_not_induced_vs_All, 'DE_PTGER4_NULL_ARAsco.csv')
#DE heatmap
Idents(arasco.obj) = arasco.obj$type
DE_ARA_Control <- FindMarkers(arasco.obj, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = .1, assay = 'RNA')

DE_ARA_Control1 <- DE_ARA_Control
DE_ARA_Control1<- DE_ARA_Control1[!grepl("mt-", rownames(DE_ARA_Control1)),]
DE_ARA_Control1<- DE_ARA_Control1[!grepl("Rpl", rownames(DE_ARA_Control1)),]
DE_ARA_Control1<- DE_ARA_Control1[!grepl("Rps", rownames(DE_ARA_Control1)),]
DE_ARA_Control1<- DE_ARA_Control1[!grepl("Atp", rownames(DE_ARA_Control1)),]
DE_ARA_Control1 <- DE_ARA_Control1[order(-DE_ARA_Control1$avg_logFC),]
write.csv(DE_ARA_Control1 ,'ARAsco_DE.csv')
ARAVC_top_genes <- rownames(DE_ARA_Control1)[1:100]


Idents(arasco.obj) = arasco.obj$Cell_type
Sub_cluster <- subset(arasco.obj,  idents = 'Stem 1')
Idents(Sub_cluster) = Sub_cluster$type
Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
DE_val1 <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", features = ARAVC_top_genes,  test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
DE_val1 <- DE_val1[order(rownames(DE_val1)),]
Cluster_ARAvC <- subset(DE_val1, select = c('avg_logFC'))
Cluster_ARAvC$gene_name <- rownames(Cluster_ARAvC)
clusters <- c('Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')

for(i in clusters){
  Sub_cluster <- subset(arasco.obj,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
  DE_val <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", features = ARAVC_top_genes, test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
  DE_val <- DE_val[order(rownames(DE_val)),]
  DE_val <- subset(DE_val, select = c('avg_logFC'))
  Cluster_ARAvC <- cbindX(Cluster_ARAvC, DE_val)
  #DE_val$gene_name <- rownames(DE_val)
  #Cluster_ARAvC <- merge(Cluster_ARAvC, DE_val, by="gene_name", all = T)
}
Cluster_ARAvC
Cluster_ARAvC_genes <- Cluster_ARAvC$gene_name
write.csv(Cluster_ARAvC, "Cluster_ARAvC.csv")
Cluster_ARAvC <- read_csv("Cluster_ARAvC.csv")
FeaturePlot(arasco.obj, features = 'Ptger4')
ggsave('ptger4.svg')
prominent_genes <-  c(Radiation_Genes,Fetal_Genes,Homeostatic_Genes,Regeneration_Genes,CREB_Targets_Genes,Granuloma_Genes,BCatenin_Enhanced_Genes,YAP_Targets_Genes,ECM_Induced_Genes,ASCL_Targets_Genes)

Prom_in_DE <- intersect(Cluster_ARAvC_genes, prominent_genes)
YAP_in_DE <- intersect(Cluster_ARAvC_genes, YAP_Targets_Genes)
WNT_in_DE <- intersect(Cluster_ARAvC_genes, BCatenin_Enhanced_Genes)

library(viridis)
heatmap1 <-Cluster_ARAvC
heatmap1 <- subset(heatmap1, select = -c(Gene_Name))
rownames(heatmap1) <- Cluster_ARAvC$Gene_Name
dt2 <- heatmap1 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)

genes_DE <- Cluster_ARAvC$Gene_Name
a <-'remove this'
for(i in genes_DE){
  if(length(intersect(i,WNT_in_DE)) == 1){
    a <- c(a, "red")}
  else{
    a <- c(a, "black")}}
a <- a[2:101]

ggplot(dt2, aes(x = rowname, y = colname, fill = value)) + geom_tile(color="white", size=0.25) + scale_fill_viridis(discrete=FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = a), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Gene") + ylab("Cell Type") + labs(fill = "Log FC") + scale_y_discrete(limits = c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft'))

#diff express all clust
names <- c('Stem 1', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
Idents(arasco.obj) <- arasco.obj$Cell_type
for(i in names){
  Sub_cluster <- subset(arasco.obj,  idents = i)
  Idents(Sub_cluster) <- Sub_cluster$type
  Gene_pval <- FindMarkers(Sub_cluster, ident.1 = "Arasco", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
  write.csv(Gene_pval, file = paste0(i, '_Arasco_gene_pval.csv'))
}
Enterocyte_Proximal_Arasco_gene_pval
Enterocyte_Progenitor_Arasco_gene_pval
Enteroendocrine_Arasco_gene_pval
Goblet_Arasco_gene_pval
Paneth_Arasco_gene_pval
Stem_1_Arasco_gene_pval
Transit_Amplifying_Arasco_gene_pval
Tuft_Arasco_gene_pval

rownames(Enterocyte_Distal_Arasco_gene_pval) <- Enterocyte_Distal_Arasco_gene_pval$X1
rownames(Enterocyte_Proximal_Arasco_gene_pval) <- Enterocyte_Proximal_Arasco_gene_pval$X1
rownames(Enterocyte_Progenitor_Arasco_gene_pval) <- Enterocyte_Progenitor_Arasco_gene_pval$X1
rownames(Enteroendocrine_Arasco_gene_pval) <- Enteroendocrine_Arasco_gene_pval$X1
rownames(Goblet_Arasco_gene_pval) <- Goblet_Arasco_gene_pval$X1
rownames(Paneth_Arasco_gene_pval) <- Paneth_Arasco_gene_pval$X1
rownames(Stem_1_Arasco_gene_pval) <- Stem_1_Arasco_gene_pval$X1
rownames(Transit_Amplifying_Arasco_gene_pval) <- Transit_Amplifying_Arasco_gene_pval$X1
rownames(Tuft_Arasco_gene_pval) <- Tuft_Arasco_gene_pval$X1

Enterocyte_Distal <- ''
Enterocyte_Proximal <- ''
Enterocyte_Progenitor <- ''
Enteroendocrine <- ''
Goblet <- ''
Paneth <- ''
Stem_1 <- ''
Transit_Amplifying <- ''
Tuft <- ''

gene_list <-  c('S100a6','Lgr5','Ascl2','Ly6a','Cd55','Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')

for(i in gene_list){
  Enterocyte_Distal <- c(Enterocyte_Distal,Enterocyte_Distal_Arasco_gene_pval[i,]$p_val_adj)
  Enterocyte_Proximal <- c(Enterocyte_Proximal,Enterocyte_Proximal_Arasco_gene_pval[i,]$p_val_adj)
  Enterocyte_Progenitor <- c(Enterocyte_Progenitor,Enterocyte_Progenitor_Arasco_gene_pval[i,]$p_val_adj)
  Enteroendocrine <- c(Enteroendocrine,Enteroendocrine_Arasco_gene_pval[i,]$p_val_adj)
  Goblet <- c(Goblet,Goblet_Arasco_gene_pval[i,]$p_val_adj)
  Paneth <- c(Paneth,Paneth_Arasco_gene_pval[i,]$p_val_adj)
  Stem_1 <- c(Stem_1,Stem_1_Arasco_gene_pval[i,]$p_val_adj)
  Transit_Amplifying <- c(Transit_Amplifying,Transit_Amplifying_Arasco_gene_pval[i,]$p_val_adj)
  Tuft <- c(Tuft,Tuft_Arasco_gene_pval[i,]$p_val_adj)
}

Enterocyte_Distal <- Enterocyte_Distal[2:29]
Enterocyte_Proximal <- Enterocyte_Proximal[2:29]
Enterocyte_Progenitor <- Enterocyte_Progenitor[2:29]
Enteroendocrine <- Enteroendocrine[2:29]
Goblet <- Goblet[2:29]
Paneth <- Paneth[2:29]
Stem_1 <- Stem_1[2:29]
Transit_Amplifying <- Transit_Amplifying[2:29]
Tuft <- Tuft[2:29]

Stem_2 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
PGE2_heatmap <- data.frame(Stem_1, Stem_2, Transit_Amplifying, Enterocyte_Progenitor, Enterocyte_Proximal, Enterocyte_Distal, Enteroendocrine, Goblet, Paneth, Tuft, stringsAsFactors=FALSE)

rownames(PGE2_heatmap) <- gene_list
PGE2_heatmap[is.na(PGE2_heatmap)] <- 0
PGE2_heatmap <- data.frame(sapply(PGE2_heatmap, as.numeric))

PGE2_heatmap1 <- data.frame(PGE2_heatmap)
rownames(PGE2_heatmap1) <- gene_list
colnames(PGE2_heatmap1) <- c('Stem 1', 'Stem 2' ,'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
write.csv(PGE2_heatmap1, 'ARAsco_MAST_PVALS.csv')

dt2 <- PGE2_heatmap1 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
geneColorsDiff = brewer.pal(n=5, "RdYlBu")
ggplot(dt2, aes(x = rowname, y = colname, fill = value)) + scale_fill_viridis_c(option = "magma") + geom_tile(color="white", size=0.25) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Gene") + ylab("Cell Type") + labs(fill = "Log FC") + scale_y_discrete(limits = c('Stem 1', 'Stem 2' ,'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft'))



#housekeeping
saveRDS(arasco.obj, file = "arasco_Sobj.rds")

arasco.obj <- readRDS('./data/Seurat_Objects/arasco_Sobj.rds', refhook = NULL)

Idents(arasco.obj) <- arasco.obj$Cell_type
current.cluster.ids <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
new.cluster.ids <- c('Homeostatic Stem','Revival Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
arasco.obj@active.ident <- plyr::mapvalues(x = arasco.obj@active.ident, from = new.cluster.ids, to = current.cluster.ids)
DimPlot(arasco.obj)
arasco.obj$Cell_type <- Idents(arasco.obj)
DimPlot(arasco.obj, group.by = "Cell_type")

gene.level.list <- c('Lgr5', 'Ascl2', 'Olfm4', 'Clu', 'Ly6a', 'S100a6', 'Bmi1', 'Hopx', 'Lrig1', 'Mex3a')
VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 5, split.by = "type", features = gene.level.list, pt.size = 0, assay = "RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE)
VlnPlot(arasco.obj1, group.by = 'Cell_type' , ncol = 5, split.by = "type", features = gene.level.list, pt.size = 0, assay = "RNA",slot = 'scale.data', log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE)

identities <- levels(arasco.obj$Cell_type)
my_color_palette <- hue_pal()(length(identities))
