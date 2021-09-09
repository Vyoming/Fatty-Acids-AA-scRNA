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
# Load in the Data
Vyom_color_Umap <- CustomPalette(low = "purple", high = "yellow", mid = "black", k = 100)

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
VlnPlot(organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'type', pt.size = 0)
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
options (future.globals.maxSize = 4000 * 1024^10)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
rm(AA4, AA5, C5, contr4, P5, Data4, Data5 ,AA4.data, AA5.data, C5.data, contr4.data, P5.data)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, Data)
organoid <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)
rm(Data.anchors)
# Visulization and Clustering
organoid <- RunPCA(organoid, verbose = FALSE)

organoid <- RunUMAP(organoid, dims = 1:10)

levels(factor(organoid@meta.data$orig.ident))
Idents(organoid) <- organoid$orig.ident
organoid[["Type"]] <- Idents(organoid)
new.cluster.ids <- c("AA", "AA", "Control", "Control", "Pge2")
names(new.cluster.ids) <- levels(organoid)
organoid <- RenameIdents(organoid, new.cluster.ids)

organoid[["Type"]] <- Idents(organoid)
head(Idents(organoid))
plots <- DimPlot(organoid, group.by = c("Type"))
plots 

organoid <- FindNeighbors(organoid, dims = 1:10)
organoid <- FindClusters(organoid, resolution = .3)#.3
DimPlot(organoid, reduction = "umap", label = TRUE)

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
revival.markers <- c(Roulis_gene_lists_for_metagenes$Gene)

organoid <- AddModuleScore(object = organoid, features = list(pathways$LGR5_stem_cell_signature), name = 'Stem')
organoid <- AddModuleScore(object = organoid, features = list(revival.markers), name = 'revival.markers', assay =  'RNA')
organoid <- AddModuleScore(object = organoid, features = list(TA.markers), name = 'TA.markers')
organoid <- AddModuleScore(object = organoid, features = list(paneth.markers), name = 'paneth.markers')
organoid <- AddModuleScore(object = organoid, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
organoid <- AddModuleScore(object = organoid, features = list(EP.markers), name = 'EP.markers')
organoid <- AddModuleScore(object = organoid, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
organoid <- AddModuleScore(object = organoid, features = list(Goblet.markers), name = 'Goblet.markers')
organoid <- AddModuleScore(object = organoid, features = list(Tuft.markers), name = 'Tuft.markers')
plasma <- viridis(10, direction = 1, option = "C")
FeaturePlot(object = organoid, features = 'stem1', pt.size = .001)  + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Stem.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Goblet.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Goblet.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'paneth.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_paneth.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Enterocyte.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Enterocyte.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Enteroendocrine.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Enteroendocrine.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'revival.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_RevStem.pdf', width=2, height=2, units="in")

DotPlot(organoid, group.by = "Cell_Type", features = c('stem1' ,'revival.markers1' ,'reserve.markers1' ,'TA.markers1' ,'paneth.markers1' ,'Enterocyte.markers1' ,'EP.markers1' ,'Enteroendocrine.markers1' ,'Goblet.markers1' ,'Tuft.markers1'))

Cell_Type_Sig_Score <- data.frame(Stem=organoid$stem1, Revival=organoid$revival.markers1, Reserve=organoid$reserve.markers1, Transit_Amplifying =organoid$TA.markers1, Paneth=organoid$paneth.markers1, Enterocyte=organoid$Enterocyte.markers1, EP=organoid$EP.markers1 ,Enteroendocrine= organoid$Enteroendocrine.markers1 ,Goblet=organoid$Goblet.markers1 ,Tuft=organoid$Tuft.markers1)

# module score distribution
modulescores <- Cell_Type_Sig_Score %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p
Idents(organoid) <- organoid$Type

write.csv(z_RevSC_cluster)
z_Stem_cluster<-WhichCells(object = organoid, expression = stem1 > .5)
z_RevSC_cluster<-WhichCells(object = organoid ,idents = c('AA' ,'Pge2'), expression = revival.markers1 > 11) #idents = c('AA' ,'Pge2')  ,
z_reserve_cluster<-WhichCells(object = organoid ,idents = c('AA' ,'Pge2') , expression = reserve.markers1 > 1)
z_TA_cluster<-WhichCells(object = organoid, expression = TA.markers1 > 3.1)
z_paneth_cluster<-WhichCells(object = organoid, expression = paneth.markers1 > 13.5)
z_Enterocyte_cluster<-WhichCells(object = organoid, expression = Enterocyte.markers1 > 1)
z_EP_cluster<-WhichCells(object = organoid, expression = EP.markers1 > 3.7)
z_Enteroendocrine_cluster<-WhichCells(object = organoid, expression = Enteroendocrine.markers1 > 2)
z_Goblet_cluster<-WhichCells(object = organoid, expression = Goblet.markers1 > 12.9)
z_Tuft_cluster<-WhichCells(object = organoid, expression = Tuft.markers1 > 14)

table(organoid$Cell_type)

organoid <- FindClusters(organoid, resolution = .3)#.3
DimPlot(organoid, reduction = "umap", label = TRUE)

DimPlot(organoid, label=T, group.by="Cell_Type", cells.highlight= list(z_Tuft_cluster, z_Stem_cluster, z_RevSC_cluster, z_reserve_cluster, z_TA_cluster, z_paneth_cluster, z_Enterocyte_cluster, z_EP_cluster, z_Enteroendocrine_cluster, z_Goblet_cluster  ),  cols.highlight = c("darkblue", "darkred", 'green' ,'lightgreen', 'darkgreen', 'yellow', 'pink', 'purple', 'orange' ,'lightblue') ,cols= "grey")
z_RevSC_cluster1<-WhichCells(object = organoid ,idents = c('Pge2'), expression = revival.markers1 > 5)
z_RevSC_cluster2<-WhichCells(object = organoid ,idents = c('AA'), expression = revival.markers1 > 5)
z_RevSC_cluster3<-WhichCells(object = organoid ,idents = c('Control'), expression = revival.markers1 > 5)

DimPlot(organoid, label=T, cells.highlight= list(z_RevSC_cluster),  cols.highlight = c("darkblue") ,cols= "grey")

#Annotate Clusters
Cluster_5 <- subset(organoid, idents = 5)
Cluster_5
Cluster_5 <- RunPCA(Cluster_5, verbose = FALSE)

Cluster_5 <- RunUMAP(Cluster_5, dims = 1:10)

Cluster_5 <- FindNeighbors(Cluster_5, dims = 1:10)
Cluster_5 <- FindClusters(Cluster_5, resolution = .3)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Goblet.markers), name = 'Goblet.markers')
FeaturePlot(object = Cluster_5, features = 'Goblet.markers1')
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Tuft.markers), name = 'Tuft.markers')
FeaturePlot(object = Cluster_5, features = 'Tuft.markers1')
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(paneth.markers), name = 'paneth.markers')
FeaturePlot(object = Cluster_5, features = 'paneth.markers1')
new.cluster.ids <- c("Goblet" ,"Goblet" ,"Secretory Progenitor" ,"Paneth" ,"Paneth" ,"Tuft")
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
names(new.cluster.ids) <- levels(Cluster_5)
Cluster_5 <- RenameIdents(Cluster_5, new.cluster.ids)
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


#pre annotate clusters
new.cluster.ids <- c("Stem 1" ,"Transit Amplifying" ,"Enterocyte Progenitor" ,"Stem 1" ,"Enterocyte", "" ,"Stem 1", "Enteroendocrine", "Transit Amplifying", "Stem Like Progenitor")
organoid[["Cell_type"]] <- Idents(organoid)
names(new.cluster.ids) <- levels(organoid)
organoid <- RenameIdents(organoid, new.cluster.ids)
organoid[["Cell_type"]] <- Idents(organoid)

# Generate a new column called Cell_type in the metadata
organoid$Cell_Type <- as.character(Idents(organoid))

# Change the information of cells containing sub-cluster information z_reserve_cluster
organoid$Cell_Type[WhichCells(Cluster_5)] <- paste(Idents(Cluster_5))
Idents(organoid) <- organoid$Cell_type
Idents(organoid) <- organoid$integrated_snn_res.0.3
#6,8
Stem_2_x <- WhichCells(object = organoid ,idents = c('6'), z_RevSC_cluster)

Stem3_x <- WhichCells(object = organoid ,idents = c('8'), z_RevSC_cluster)

Stem2_cells <- c(Stem_2_x, z_EP_cluster)
Stem3_cells <- c(Stem3_x, z_Tuft_cluster)

Idents(organoid) <- organoid$Cell_Type
organoid$Cell_Type[Stem2_cells] <- paste('Stem 2')
organoid$Cell_Type[Stem3_cells] <- paste('Stem 3')

table(organoid$Cell_Type, organoid$Type)
#organoid$Cell_Type[WhichCells(object = organoid ,idents = c('AA' ,'Pge2') , expression = reserve.markers1 > 1.5)] <- paste('Stem 3')

DimPlot(organoid, group.by = "Cell_Type", label = FALSE, pt.size=1, label.size = 3)
Idents(organoid) <- organoid$Cell_Type
organoid$Cell_type <- organoid$Cell_Type
my_levels <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
organoid$Cell_type <- factor(x = organoid$Cell_Type, levels = my_levels)
DimPlot(organoid, group.by = "Cell_type", label = FALSE, pt.size=1, label.size = 2)
Idents(organoid) <- organoid$Type
my_levels <- c('Control' ,'AA' ,'Pge2')
organoid$type <- factor(x = organoid$Type, levels = my_levels)

#Cluster Markers
Idents(organoid) <- organoid$Cell_type

markers_stem <- FindAllMarkers(organoid, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
markers_stem %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(markers_stem,'organoid_Sigs_Per_Clust.csv')

Idents(organoid) <- organoid$Cell_type
organoid <- NormalizeData(object = organoid ,normalization.method = "LogNormalize", assay = "RNA")

Stem1V2 <- FindMarkers(organoid, ident.1 = "Stem 2", ident.2 = "Stem 3", test.use = "MAST", logfc.threshold = .05, min.pct = 0.15, assay = 'RNA')
Stem1V2<- Stem1V2[!grepl("mt-", rownames(Stem1V2)),]
Stem1V2<- Stem1V2[!grepl("Rpl", rownames(Stem1V2)),]
Stem1V2<- Stem1V2[!grepl("Rps", rownames(Stem1V2)),]
Stem1V2<- Stem1V2[!grepl("Atp", rownames(Stem1V2)),]
Stem1V2 <- Stem1V2[order(-Stem1V2$avg_logFC),]
write.csv(Stem1V2 ,'Stem2_vs_Stem3.csv')

EnhancedVolcano(Stem1V2, lab = rownames(Stem1V2), x = 'avg_logFC', y = 'p_val_adj', title = 'Stem 2 vs Stem 3', pCutoff = 10e-8, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,60), subtitle = '' )

write.csv(Stem1V2, 'stem2v3.csv')
#Proportions per Cell Type
Prop_table<- prop.table(x = table(Gene_de$Cell_type, Gene_de$type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
#Prop_Table1$Freq <- 1/log10(Prop_Table$Freq)*-1
Prop_Table1$Freq[2:3] <- c(0,0)
my_levels <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','AA'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Gene_de_prop.svg", plot=plot, width=10, height=8)

#prop sigs
prop_sig <- table(Gene_de$Cell_type, Gene_de$orig.ident)
Prop_Sig_Table <- as.data.frame.matrix(prop_sig)
num <- length(rownames(Prop_Sig_Table))

for(i in 1:num){
  dat = data.frame(sampleID = c(colnames(Prop_Sig_Table)), current_cluster=c(Prop_Sig_Table[i,]$AA4, Prop_Sig_Table[i,]$AA5, Prop_Sig_Table[i,]$C5, Prop_Sig_Table[i,]$Control4), all_other_clusters=c(sum(Prop_Sig_Table$AA4) - Prop_Sig_Table[i,]$AA4, sum(Prop_Sig_Table$AA5) - Prop_Sig_Table[i,]$AA5, sum(Prop_Sig_Table$C5) - Prop_Sig_Table[i,]$C5, sum(Prop_Sig_Table$Control4) - Prop_Sig_Table[i,]$Control4), Treatment = c("AA", "AA", "Control", "Control"))
  glm_output = glm(data = dat, cbind(current_cluster, all_other_clusters) ~ Treatment, family = "binomial")
  print(coef(summary(glm_output))[,'Pr(>|z|)'])
}

#Gene level analysis
Idents(Gene_de) <- Gene_de$type
Gene_de <- subset(organoid,  idents =c('Control', 'AA'))
Gene_de <- NormalizeData(object = Gene_de ,normalization.method = "LogNormalize", assay = "RNA")
# do some abra kadabra magic

my_levels <- c('Control' ,'AA')
Gene_de$type <- factor(x = Gene_de$Type, levels = my_levels)
DefaultAssay(Gene_de) <- "RNA"
Gene_de <- magic(Gene_de)


#zscore testing:
library(readxl)
CREB_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Creb_Targets")
YAP_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Yap_Targets")
BCatenin_Enhanced <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "bcatenin_enhanced")
Fetal <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Fetal_Spheroid")
Granuloma <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Granuloma_Induced")
Radiation <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Radiation_induced")
ECM_Induced <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "ECM_Induced")
ASCL_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "ASCL_Targets")
Regeneration <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Regeneration_Induced")
Homeostatic <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Homeostatic_Signature")

All_Genes <- Gene_de@assays$RNA@data@Dimnames[[1]]
Radiation_Genes <- intersect(All_Genes, Radiation$Gene_Name)
Fetal_Genes <- intersect(All_Genes, Fetal$Gene_Name)
Homeostatic_Genes <- intersect(All_Genes, Homeostatic$Gene_Name)
Regeneration_Genes <- intersect(All_Genes, Regeneration$Gene_Name)
Granuloma_Genes <- intersect(All_Genes, Granuloma$Gene_Name)
CREB_Targets_Genes <- intersect(All_Genes, CREB_Targets$Gene_Name)
BCatenin_Enhanced_Genes <- intersect(All_Genes, BCatenin_Enhanced$Gene_Name)
YAP_Targets_Genes <- intersect(All_Genes, YAP_Targets$Gene_Name)
ECM_Induced_Genes <- intersect(All_Genes, ECM_Induced$Gene_Name)
ASCL_Targets_Genes <- intersect(All_Genes, ASCL_Targets$Gene_Name)

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[Radiation_Genes, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Radiation_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[Fetal_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Fetal_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[Granuloma_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Granuloma_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[Regeneration_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Regeneration_Score <- mean.exp
}
mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[Homeostatic_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$Homeostatic_Stem_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[CREB_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$CREB_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[BCatenin_Enhanced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$BCatenin_Enhanced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[YAP_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$YAP_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[ECM_Induced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$ECM_Induced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = Gene_de@assays$MAGIC_RNA@data[ASCL_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = Gene_de@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Gene_de@meta.data$ASCL_Targets_Score <- mean.exp
}

Radiation_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Radiation_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')
Fetal_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Fetal_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Granuloma_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Granuloma_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Regen_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Regeneration_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Homeostatic_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Homeostatic_Stem_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
YAP_graph <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('YAP_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
wnt_graph <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('BCatenin_Enhanced_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Creb_graph <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('CREB_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
ASCL_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('ASCL_Targets_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
ECM_grpah <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('ECM_Induced_Score'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')

Score_final <- plot_grid(Radiation_grpah, Fetal_grpah, Granuloma_grpah, Regen_grpah, ECM_grpah, ncol = 5, labels = c('A','B','C','D','E'), label_size = 12)
Score_final

ggsave(file = 'Radiation_grpah.svg', plot=Radiation_grpah, width=10, height=10)
ggsave(file = 'fetal.svg', plot=Fetal_grpah, width=10, height=10)
ggsave(file = 'granuloma.svg', plot=Granuloma_grpah, width=10, height=10)
ggsave(file = 'regen.svg', plot=Regen_grpah, width=10, height=10)
ggsave(file = 'Homeostatic_stem.svg', plot=Homeostatic_grpah, width=10, height=10)


gene.level.list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'Cd55')
gene.level.list <- c('Radiation_Score', 'Fetal_Score', 'Granuloma_Score', 'Regeneration_Score', 'Homeostatic_Score')
gene.level.list <- c('Areg', 'Ereg', 'Sapcd2')
gene.level.list <- c('Ptger1', 'Ptger2', 'Ptger3', 'Ptges', 'Ptgs1', 'Ptgs2')
All <- VlnPlot(arasco.obj, group.by = 'Cell_type' ,  split.by = "type", features = gene.level.list, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + 
  scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + 
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
ggsave(file = paste0('arasco_PGE_Gene.pdf'), plot=All, width=5, height=3, units="in")


Gene_Level_PVal_adj <- data.frame(Gene_Level_PVal1)
Gene_Level_PVal_adj$Lgr5_text = ""
Gene_Level_PVal_adj$Lgr5_text[Gene_Level_PVal_adj$Lgr5 < 0.05] = "*"
Gene_Level_PVal_adj$Lgr5_text[Gene_Level_PVal_adj$Lgr5 < 0.01] = "**"
Gene_Level_PVal_adj$Lgr5_text[Gene_Level_PVal_adj$Lgr5 < 0.001] = "***"
ggplot(curRes, aes(x=pathway, y=NES, fill=Tissue)) + ....... + geom_text(aes(y=NES+0.5*sign(NES), label=ptext, color=Tissue, fontface="bold"), position=position_dodge2(0.9), vjust=0.75) + ....


Lgr5 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Lgr5'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('') 
Ascl2 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ascl2'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
S100a6 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('S100a6'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('')
Ly6a <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ly6a'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Cd55 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Cd55'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')


PGE2_genes <-  c('Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')

for(i in PGE2_genes){
  Gene_plot_pge <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c(i), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('') 
  ggsave(file = paste0('Organoid_', i, '.svg'), plot=Gene_plot_pge, width=10, height=10)
}


ggsave(file = 'Lgr5.svg', plot=Lgr5, width=10, height=10)
ggsave(file = 'Ascl2.svg', plot=Ascl2, width=10, height=10)
ggsave(file = 'S100a6.svg', plot=S100a6, width=10, height=10)
ggsave(file = 'Ly6a.svg', plot=Ly6a, width=10, height=10)
ggsave(file = 'Cd55.svg', plot=Cd55, width=10, height=10)

VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Areg'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('') + 
  theme(legend.position = 'none', text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
ggsave(file = 'Organoid_areg.pdf',  width=2, height=3, units="in")
VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ereg'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('')+ 
  theme(legend.position = 'none', text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
ggsave(file = 'Organoid_ereg.pdf', width=2, height=3, units="in")

VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Areg'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('') + 
  theme(legend.position = 'none', text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
ggsave(file = 'Arasco_areg.pdf', width=2, height=3, units="in")
VlnPlot(arasco.obj, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ereg'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE, combine = TRUE) + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('') + 
  theme(legend.position = 'none', text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
ggsave(file = 'Arasco_ereg.pdf', width=2, height=3, units="in")

Gene_final <- plot_grid(Lgr5, Ascl2, S100a6, Ly6a, Cd55, ncol = 5, labels = c('A','B','C','D','E'), label_size = 12)
Gene_final
Lgr5 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Lgr5'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr51 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('S100a6'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr52 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Hopx'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr53 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Ccnb1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr54 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Fabp1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr55 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Agr2'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr56 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Defa17'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr57 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Trpm5'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr58 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Chgb'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr59 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Tk1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr510 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Stmn1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Lgr511 <- VlnPlot(Gene_de, group.by = 'Cell_type' , ncol = 1,  slot = "data", features = c('Atoh1'), pt.size = 0, assay = "RNA", log = TRUE) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
alt_final <- plot_grid(Lgr5,Lgr51, Lgr52, Lgr53, Lgr54, Lgr55, Lgr56, Lgr57, Lgr58, Lgr59, Lgr510, Lgr511, ncol = 4, labels = c('A','B','C','D','E','F','G','H','I','J', 'K','L'), label_size = 12)

VlnPlot(Gene_de, group.by = 'Cell_type' , split.by = 'type', ncol = 1,  slot = "data", split.plot = TRUE, features = c('Cd55'), pt.size = .1, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02')) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
VlnPlot(Gene_de, group.by = 'Cell_type' , split.by = 'type', ncol = 1,  slot = "data", split.plot = TRUE, features = c('Cd55'), pt.size = .1, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#d95f02')) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')

DimPlot(Gene_de, group.by = "Cell_type", label = FALSE, pt.size=1, label.size = 2)

Idents(Gene_de) <- Gene_de$Cell_type
cluster.markers <- FindAllMarkers(Gene_de, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'SCT')
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

check_heatmap <- DoHeatmap(Gene_de, features = top10$gene, group.by = 'Cell_type', assay = 'SCT',size = 2 , angle = 60) + NoLegend() + scale_fill_viridis(option="inferno")#top10$Gene
ggsave(file = 'Gene_de_top10_cluster_Heatmap.svg', plot=check_heatmap, width=10, height=15)
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg",'Cbx3','Larp1','Slc39a1','Hnf4a','Sox4','Mmp7','Dll1','Tff3',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a","Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
check_heatmap <- DoHeatmap(Gene_de, features = DotPlot_Sig, group.by = 'Cell_type', assay = 'SCT',size = 2 , angle = 60) + NoLegend() + scale_fill_viridis(option="inferno")#top10$Gene
ggsave(file = 'Gene_de_haber_cluster_Heatmap.svg', plot=check_heatmap, width=10, height=10)


write.csv(cluster.markers,'Gene_de_Sigs_Per_Clust.csv')
geneColorsDiff = brewer.pal(n=5, "RdYlBu")
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg",'Cbx3','Larp1','Slc39a1','Hnf4a','Sox4','Mmp7','Dll1','Tff3',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a","Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
DotPlot(Gene_de, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text= element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file= 'Organoid_Haber_check_dotplot.pdf', width=4.5, height=2.6, units="in")


DE_PGE2_Control <- FindMarkers(Gene_de, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0.15, assay = 'SCT')
DE_AA_Control <- FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'SCT')
DE_PGE2_AA <- FindMarkers(Gene_de, ident.1 = "Pge2", ident.2 = "AA", test.use = "MAST", logfc.threshold = log(2), min.pct = 0.15, assay = 'SCT')
Idents(Gene_de) <-  Gene_de$type
Sig_check <- FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control", test.use = "wilcox",features = gene.level.list ,logfc.threshold = 0, min.pct = 0, assay = 'SCT')

library(EnhancedVolcano)
library(extrafont)
extrafont::font_import()

EnhancedVolcano(DE_AA_Control, lab = rownames(DE_AA_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'AA versus Control', pCutoff = 10e-32, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,325), )
EnhancedVolcano(DE_PGE2_Control, lab = rownames(DE_PGE2_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'PGE2 versus Control', pCutoff = 10e-32, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,325))

write.csv(DE_AA_Control,'AA_DE.csv')
write.csv(DE_PGE2_Control,'PGE2_DE.csv')

Gene_de1 <- Gene_de
DefaultAssay(Gene_de1) <- "RNA"
Gene_de1 <- magic(Gene_de1)
DimPlot(Gene_de1)


DE_PGE2_Control <- FindMarkers(Gene_de, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0.15, assay = 'SCT')
DE_AA_Control <- FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'SCT')
DE_PGE2_AA <- FindMarkers(Gene_de, ident.1 = "Pge2", ident.2 = "AA", test.use = "MAST", logfc.threshold = log(2), min.pct = 0.15, assay = 'SCT')

library(EnhancedVolcano)
library(extrafont)
extrafont::font_import()

EnhancedVolcano(DE_AA_Control, lab = rownames(DE_AA_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'AA versus Control', pCutoff = 10e-32, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,325), )
EnhancedVolcano(DE_PGE2_Control, lab = rownames(DE_PGE2_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'PGE2 versus Control', pCutoff = 10e-32, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,325))

#AA vs PGE2 concordance

Idents(organoid) = organoid$type
DE_PGE2_Control <- FindMarkers(organoid, ident.1 = "Pge2", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.15, assay = 'SCT')
DE_AA_Control <- FindMarkers(organoid, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.15, assay = 'SCT')
write.csv(DE_PGE2_Control, 'DE_PGE2_Control.csv')
write.csv(DE_AA_Control, 'DE_AA_Control.csv')
DE_PGE2_Control <- read.csv('./data/AA_ARA_Intestine_project/DE_PGE2_Control.csv')
DE_AA_Control <- read.csv('./data/AA_ARA_Intestine_project/DE_AA_Control.csv')
rownames(DE_AA_Control) <- DE_AA_Control$X1
rownames(DE_PGE2_Control) <- DE_PGE2_Control$X1
DE_PGE2_Control$gene_name <- rownames(DE_PGE2_Control)
DE_AA_Control$gene_name <- rownames(DE_AA_Control)

matched <- intersect(DE_AA_Control$gene_name, DE_PGE2_Control$gene_name)

library(dplyr)
PGEvC <- which(DE_PGE2_Control$gene_name %in% matched)
PGEvC <- DE_PGE2_Control[PGEvC,]
PGEvC <- subset(PGEvC, select = c('gene_name','avg_logFC', 'p_val_adj'))
PGEvC <- PGEvC[order(PGEvC$gene_name),]

AAvC <- which(DE_AA_Control$gene_name %in% matched)
AAvC <- DE_AA_Control[AAvC,]
AAvC <- subset(AAvC, select = c('gene_name','avg_logFC', 'p_val_adj'))
AAvC <- AAvC[order(AAvC$gene_name),]


Concordance <- AAvC
Concordance$logFC_AA <- Concordance$avg_logFC
Concordance$p_val_AA <- Concordance$p_val_adj
Concordance$logFC_PGE <- PGEvC$avg_logFC
Concordance$p_val_PGE <- PGEvC$p_val_adj
Concordance <- subset(Concordance, select = c('gene_name','logFC_AA','logFC_PGE', 'p_val_AA','p_val_PGE'))
Concordance <- read.csv('./data/AA_ARA_Intestine_project/Concordance.csv')

AA_up <- c(nrow(Concordance[(Concordance$logFC_AA > 0 & Concordance$logFC_PGE > 0), ]), nrow(Concordance[(Concordance$logFC_AA > 0 & Concordance$logFC_PGE < 0), ]))
AA_not_UP <- c(nrow(Concordance[(Concordance$logFC_AA < 0 & Concordance$logFC_PGE > 0), ]), nrow(Concordance[(Concordance$logFC_AA < 0 & Concordance$logFC_PGE < 0), ]))
test_concor <- cbind(AA_up, AA_not_UP)
pval_concordance <- fisher.test(test_concor)

AA_not_down <- c(nrow(Concordance[(Concordance$logFC_AA > 0 & Concordance$logFC_PGE < 0), ]), nrow(Concordance[(Concordance$logFC_AA > 0 & Concordance$logFC_PGE > 0), ]))
AA_down <- c(nrow(Concordance[(Concordance$logFC_AA < 0 & Concordance$logFC_PGE < 0), ]), nrow(Concordance[(Concordance$logFC_AA < 0 & Concordance$logFC_PGE > 0), ]))
test_concor <- cbind(AA_down, AA_not_down )
pval_concordance <- fisher.test(test_concor)
pval_concordance$p.value

Concordance <- Concordance[(Concordance$p_val_AA < .05 & Concordance$p_val_PGE < .05), ]
write.csv(Concordance, 'Concordance.csv')
#plot concordance
library("tidyverse")
library("scales")
library("ggnewscale")
library("cowplot")
library("circlize")
library('EnhancedVolcano')

notable <- c('S100a6','Cldn3','Lyz1', 'Anxa3', 'Scd2', 'Ccnd1','Olfm4','Ndufa7','Agr2','Anxa2','Taldo1', 'Areg', 'Ascl2', 'Muc4', 'Bad', 'Cox20')
concordance_plot <- function(fc, x) {
  title_size <- 10
  text_size <- 8
  point_size <- 0.6
  highlights <- which(fc$gene_name %in% x)
  highlights <- fc[highlights,]
  range_fc <- range(select(fc, starts_with("log")))
  nrows <- nrow(fc) + 1
  a <-'remove this'
  for(i in fc$gene_name){
    x <- which(fc$gene_name %in% i)
    x <- fc[x,]
    if((x$logFC_AA < 0)&(x$p_val_AA < .01)&(x$logFC_PGE < 0)&(x$p_val_PGE < .01)){
      a <- c(a, "#d95f02")}
    else if ((x$p_val_AA < .01)&&(x$logFC_AA > 0)&&(x$p_val_PGE < .01)&&(x$logFC_PGE > 0)){
      a <- c(a, "#7570b3")}
    else{
      a <- c(a, "#efefef")}}
  a <- a[2:nrows]
  
  b <-'remove this'
  for(i in highlights$gene_name){
    x <- which(highlights$gene_name %in% i)
    x <- highlights[x,]
    if((x$logFC_AA < 0)&(x$p_val_AA < .01)&(x$logFC_PGE < 0)&(x$p_val_PGE < .01)){
      b <- c(b, "#d95f02")}
    else if ((x$p_val_AA < .01)&&(x$logFC_AA > 0)&&(x$p_val_PGE < .01)&&(x$logFC_PGE > 0)){
      b <- c(b, "#7570b3")}
    else{
      b <- c(b, "#efefef")}}
  b <- a[2:12]
  correlation <- cor(fc$logFC_AA,
                     fc$logFC_PGE)
  AA_up <- c(nrow(fc[(fc$logFC_AA > 0 & fc$logFC_PGE > 0), ]), nrow(fc[(fc$logFC_AA > 0 & fc$logFC_PGE < 0), ]))
  AA_down <- c(nrow(fc[(fc$logFC_AA < 0 & fc$logFC_PGE > 0), ]), nrow(fc[(fc$logFC_AA < 0 & fc$logFC_PGE < 0), ]))
  test_concor <- cbind(AA_up, AA_down)
  pval_concordance <- fisher.test(test_concor)$p.value
  figConcordance <- ggplot(fc,
                           aes(x=logFC_AA,
                               y=logFC_PGE)) +
    ggrastr::geom_point_rast(aes(color= a), show.legend = c("Downregulated", "Upregulated", ""),
                             size=point_size, alpha=0.1,  raster.dpi = 2000) +
    ggrepel::geom_text_repel(data=highlights,
                             aes(x=logFC_AA,
                                 y=logFC_PGE,
                                 label=gene_name),
                             size = 2, color='Black',
                             point.padding=0.1,
                             segment.size = 0.2,
                             min.segment.length=0.2) +
    geom_point(data=highlights,
               aes(x=logFC_AA,
                   y=logFC_PGE),
               size = point_size, color=c('Black')) +
    scale_x_continuous(limits = range_fc) +
    scale_y_continuous(limits = range_fc) +
    scale_color_manual(breaks = c("Downregulated", "Upregulated", ""), values=c( "Red", "Blue", "#b9b9b9"), guide = "legend") +
    labs(x="log2FC(AA/control)",
         y="log2FC(PGE2/control)"
    ) +
    geom_text(aes(x=range_fc[1]+0.1, y=range_fc[2]-0.1),
              label=paste0("R^2==~", round(correlation,2)), parse = TRUE,
              hjust = 0, size=3,) +
    geom_text(aes(x=range_fc[1]+1.1, y=range_fc[2]-0.11),
              label=paste0("p < .001"), parse = TRUE,
              hjust = 0, size=3,) +
    coord_equal() +
    theme_half_open() + 
    theme(axis.text = element_text(size = text_size),
          axis.title = element_text(size = title_size),
          legend.title = element_text(size=title_size),
          legend.text = element_text(size=text_size))
  
  fpath <- c('Concordance.pdf')
  ggsave(plot=figConcordance,
         fpath,
         width = 2.5,
         height = 2.5,
         units = "in")
  

}
concordance_plot(Concordance, notable)
pval_concordance$p.value

#DE analysis

#compare different clusters
Idents(organoid) <- organoid$Cell_type
Stem2v3 <- FindMarkers(organoid, ident.1 = "Stem 2", ident.2 = "Stem 3", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

Stem2v3<- Stem2v3[!grepl("mt-", rownames(Stem2v3)),]
Stem2v3<- Stem2v3[!grepl("Rpl", rownames(Stem2v3)),]
Stem2v3<- Stem2v3[!grepl("Rps", rownames(Stem2v3)),]
Stem2v3<- Stem2v3[!grepl("Atp", rownames(Stem2v3)),]
Stem2v3 <- Stem2v3[order(-Stem2v3$avg_logFC),]
volc <- EnhancedVolcano(Stem2v3, lab = rownames(Stem2v3), x = 'avg_logFC', y = 'p_val_adj', title = 'Stem 2 vs Stem 3', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-1.25, 1.25), ylim = c(0,150),  shade = c(''), shadeAlpha = 0 )
ggsave(file = 'DE_Stem2V3.svg', plot=volc, width=10, height=10)
write.csv(Stem2v3, 'DE_Stem2V3.csv')

Stem2v1 <- FindMarkers(organoid, ident.1 = "Stem 2", ident.2 = "Stem 1", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

Stem2v1<- Stem2v1[!grepl("mt-", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Rpl", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Rps", rownames(Stem2v1)),]
Stem2v1<- Stem2v1[!grepl("Atp", rownames(Stem2v1)),]
Stem2v1 <- Stem2v1[order(-Stem2v1$avg_logFC),]
volc <- EnhancedVolcano(Stem2v1, lab = rownames(Stem2v1), x = 'avg_logFC', y = 'p_val_adj', title = 'Stem 2 vs Stem 1', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-2, 1.25), ylim = c(0,330),  shade = c(''), shadeAlpha = 0 )
ggsave(file = 'DE_Stem2v1.svg', plot=volc, width=10, height=10)
write.csv(Stem2v1, 'DE_Stem2v1.csv')

Stem1v3 <- FindMarkers(organoid, ident.1 = "Stem 3", ident.2 = "Stem 1", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

Stem1v3<- Stem1v3[!grepl("mt-", rownames(Stem1v3)),]
Stem1v3<- Stem1v3[!grepl("Rpl", rownames(Stem1v3)),]
Stem1v3<- Stem1v3[!grepl("Rps", rownames(Stem1v3)),]
Stem1v3<- Stem1v3[!grepl("Atp", rownames(Stem1v3)),]
Stem1v3 <- Stem1v3[order(-Stem1v3$avg_logFC),]
volc <- EnhancedVolcano(Stem1v3, lab = rownames(Stem1v3), x = 'avg_logFC', y = 'p_val_adj', title = 'Stem 3 vs Stem 1', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-2, 2), ylim = c(0,330),  shade = c(''), shadeAlpha = 0 )
ggsave(file = 'DE_Stem1v3.svg', plot=volc, width=10, height=10)
write.csv(Stem1v3, 'DE_Stem1v3.csv')

#AA v C
AA_induced_Genes <- read_csv("AA_induced_Genes.csv")$gene_name
Gene_de <- AddModuleScore(object = Gene_de, features = 'Ptger4', name = 'Ptger4', assay = 'RNA')
VlnPlot(Gene_de, features = 'Ptger41',group.by = 'type')
ggsave('plot.svg')
Gene_de1 <-  Gene_de
Gene_de$AA_induced <- Gene_de$Cell_type
new.cluster.ids <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Idents(Gene_de) <- Gene_de$AA_induced
current.cluster.ids <- c('Not Treated','AA Induced', 'AA Induced','Not Treated', 'Not Treated','Not Treated','Not Treated','Not Treated', 'Not Treated', 'Not Treated', 'Not Treated','Not Treated')
Gene_de$AA_induced <- plyr::mapvalues(x = Gene_de$AA_induced, from = new.cluster.ids, to = current.cluster.ids)
DimPlot(Gene_de)

AA_induced_cells <-WhichCells(object = Gene_de,  expression = Ptger41 > .01)
Idents(Gene_de) <- Gene_de$AA_induced
Gene_de$AA_induced[AA_induced_cells] <- paste('AA Induced')
table(Gene_de$AA_induced, Gene_de$Type)

Idents(Gene_de) = Gene_de$AA_induced
AA_induced_1 <- subset(Gene_de,  idents =c('AA Induced'))
AA_induced_not <- subset(Gene_de,  idents =c('Not Treated'))
Idents(AA_induced_1) <- AA_induced_1$type
Idents(AA_induced_not) <- AA_induced_not$type

DE_AA_induced_vs_All <- FindMarkers(AA_induced_1, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = .05, min.pct = 0, assay = 'RNA')

DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("mt-", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Rpl", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Rps", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All<- DE_AA_induced_vs_All[!grepl("Atp", rownames(DE_AA_induced_vs_All)),]
DE_AA_induced_vs_All <- DE_AA_induced_vs_All[order(-DE_AA_induced_vs_All$avg_logFC),]
volc <- EnhancedVolcano(DE_AA_induced_vs_All, lab = rownames(DE_AA_induced_vs_All), x = 'avg_logFC', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-3, FCcutoff = .15, xlim = c(-1.25, 1.25), ylim = c(0,150), subtitle = 'Ptger4 Expressing Cells',  shade = c('Anxa3','Ly6a','Anxa2', 'Fabp1', 'S100a6', 'Areg', 'Ldha',''), shadeAlpha = 0 )
ggsave(file = 'DE_PTGER4_expressing_organoid.svg', plot=volc, width=10, height=10)
write.csv(DE_AA_induced_vs_All, 'DE_PTGER4_expressing_organoid.csv')


DE_AA_not_induced_vs_All <- FindMarkers(AA_induced_not, ident.1 = "AA", ident.2 = "Control", test.use = "MAST",  logfc.threshold = .05, min.pct = 0, assay = 'RNA')

DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("mt-", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Rpl", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Rps", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All<- DE_AA_not_induced_vs_All[!grepl("Atp", rownames(DE_AA_not_induced_vs_All)),]
DE_AA_not_induced_vs_All <- DE_AA_not_induced_vs_All[order(-DE_AA_not_induced_vs_All$avg_logFC),]

volc <- EnhancedVolcano(DE_AA_not_induced_vs_All, lab = rownames(DE_AA_not_induced_vs_All), x = 'avg_logFC', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-25, FCcutoff = .15, xlim = c(-1.25, 1.25), ylim = c(0,330), subtitle = 'Ptger4 Null Cells' )
ggsave(file = 'DE_PTGER4_NULL_organoid.svg', plot=volc, width=10, height=10)
write.csv(DE_AA_not_induced_vs_All, 'DE_PTGER4_NULL_organoid.csv')

Idents(Gene_de) <- Gene_de$type
Gene_de1 <- Gene_de
Gene_de1 <- NormalizeData(object = Gene_de1, normalization.method = "LogNormalize", assay = "RNA")
all.genes <- rownames(Gene_de1)
Gene_de1 <- ScaleData(Gene_de1, features = all.genes, assay = 'RNA')
DE_valX <- FindMarkers(Gene_de1, ident.1 = "AA", ident.2 = "Control",slot='scale.data' , test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'SCT')

DE_valX<- DE_valX[!grepl("mt-", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Rpl", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Rps", rownames(DE_valX)),]
DE_valX<- DE_valX[!grepl("Atp", rownames(DE_valX)),]
DE_valX <- DE_valX[order(-DE_valX$avg_diff),]
EnhancedVolcano(DE_valX, lab = rownames(DE_valX), x = 'avg_diff', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-25, FCcutoff = 0.25, xlim = c(-5, 5), ylim = c(0,320), subtitle = 'All Cells' )

cluster.ids <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying")
Idents(Gene_de1) <- Gene_de1$Cell_type
Sub_cluster <- subset(Gene_de1,  idents = cluster.ids)
Idents(Sub_cluster) = Sub_cluster$type
Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
DE_val1 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control",  slot='scale.data', test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'SCT')
DE_val1<- DE_val1[!grepl("mt-", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Rpl", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Rps", rownames(DE_val1)),]
DE_val1<- DE_val1[!grepl("Atp", rownames(DE_val1)),]
DE_val1 <- DE_val1[order(-DE_val1$avg_diff),]
EnhancedVolcano(DE_val1, lab = rownames(DE_val1), x = 'avg_diff', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-10, 10), ylim = c(0,330), subtitle = 'Pluripotent Cells' )


Idents(Gene_de) <- Gene_de$type
DE_val <- FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control",  test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'SCT')
EnhancedVolcano(DE_val, lab = rownames(DE_val), x = 'avg_logFC', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-8, FCcutoff = 0.25, xlim = c(-1.25, 1.25), ylim = c(0,330), subtitle = 'All Cells' )
write.csv(DE_val, 'AAVC_DE_All.csv')

cluster.ids <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying")
Idents(Gene_de) <- Gene_de$Cell_type
Sub_cluster <- subset(Gene_de,  idents = cluster.ids)
Idents(Sub_cluster) = Sub_cluster$type
Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
DE_val1 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control",  test.use = "MAST", logfc.threshold = .05, min.pct = .15, assay = 'SCT')
EnhancedVolcano(DE_val1, lab = rownames(DE_val1), x = 'avg_logFC', y = 'p_val_adj', title = 'AA vs Control', pCutoff = 10e-8, FCcutoff = 0.25, xlim = c(-1.3, 1.3), ylim = c(0,330), subtitle = 'Pluripotent Cells' )
write.csv(DE_val1, 'AAVC_DE_pluripotent.csv')

Idents(Gene_de) = Gene_de$type
DE_AA_Control <- FindMarkers(Gene_de, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = .15, assay = 'SCT')
AAVC_top_genes <- rownames(DE_AA_Control)[1:100]
Idents(Gene_de) = Gene_de$Cell_type
Sub_cluster <- subset(Gene_de,  idents = 'Stem 1')
Idents(Sub_cluster) = Sub_cluster$type
Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
DE_val1 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", features = AAVC_top_genes,  test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
DE_val1 <- DE_val1[order(rownames(DE_val1)),]
Cluster_AAvC <- subset(DE_val1, select = c('avg_logFC'))
Cluster_AAvC$gene_name <- rownames(Cluster_AAvC)

clusters <- c("Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
for(i in clusters){
  Sub_cluster <- subset(Gene_de,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  Sub_cluster <- NormalizeData(object = Sub_cluster, normalization.method = "LogNormalize", assay = "RNA")
  DE_val <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", features = AAVC_top_genes, test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
  DE_val <- DE_val[order(rownames(DE_val)),]
  DE_val <- subset(DE_val, select = c('avg_logFC'))
  Cluster_AAvC <- cbindX(Cluster_AAvC, DE_val)
  #DE_val$gene_name <- rownames(DE_val)
  #Cluster_AAvC <- merge(Cluster_AAvC, DE_val, by="gene_name", all = T)
}
Cluster_AAvC
Cluster_AAvC_genes <- Cluster_AAvC$gene_name
write.csv(Cluster_AAvC, "Cluster_AAvC.csv")
Cluster_AAvC <- read_csv("Cluster_AAvC.csv")

prominent_genes <-  c(Radiation_Genes,Fetal_Genes,Homeostatic_Genes,Regeneration_Genes,CREB_Targets_Genes,Granuloma_Genes,BCatenin_Enhanced_Genes,YAP_Targets_Genes,ECM_Induced_Genes,ASCL_Targets_Genes)

Prom_in_DE <- intersect(Cluster_AAvC_genes, prominent_genes)
YAP_in_DE <- intersect(Cluster_AAvC_genes, YAP_Targets_Genes)
WNT_in_DE <- intersect(Cluster_AAvC_genes, BCatenin_Enhanced_Genes)
CREB_in_DE <- intersect(Cluster_AAvC_genes, CREB_Targets_Genes)

library(viridis)
heatmap1 <-Cluster_AAvC
heatmap1 <- subset(heatmap1, select = -c(Gene_Name))
rownames(heatmap1) <- Cluster_AAvC$Gene_Name
dt2 <- heatmap1 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)



a <-'remove this'
for(i in genes_DE){
  if(length(intersect(i,CREB_in_DE)) == 1){
    a <- c(a, "red")}
  else{
    a <- c(a, "black")}}
a <- a[2:101]

ggplot(dt2, aes(x = rowname, y = colname, fill = value)) + geom_tile(color="white", size=0.25) + scale_fill_viridis(discrete=FALSE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Gene") + ylab("Cell Type") + labs(fill = "Log FC") + scale_y_discrete(limits = c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" ))

#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
#install.packages("Signac")
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(Gene_de) <- Gene_de$Cell_type
pseudo_stem <- WhichCells(object = Gene_de, idents = 'Stem 1')


cds <- as.cell_data_set(Gene_de)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Gene_de[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")

saveRDS(cds, file = "Organoid_pseudotime.rds")
plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="Gene_de_Pseudotime_Umap_no_trajectory.svg", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

Gene_de <- AddMetaData(
  object = Gene_de,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

Idents(Gene_de) <- Gene_de$Cell_type
Pseudotime_umap <- FeaturePlot(Gene_de, c("Pseudotime"), pt.size = 3, label = TRUE, repel = TRUE)
ggsave(file="Gene_de_Pseudotime_Umap.svg", plot=Pseudotime_umap, width=10, height=10)
Gene_de1 <- Gene_de
Idents(Gene_de) <- Gene_de$type
AA_Control <- subset(Gene_de,  idents =c('Control'))
AA_AA <- subset(Gene_de,  idents =c('AA'))

Gene_de1@meta.data$quantile <- as.factor(ntile(Gene_de1$Pseudotime, 10))
Idents(Gene_de1) <- Gene_de1$quantile

Sub_cluster <- subset(Gene_de1, idents = c('10'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_10 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_10 <- psuedo_gene_10[order(rownames(psuedo_gene_10)),]
psuedo_gene_10 <- subset(psuedo_gene_10, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('1'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_1 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_1 <- psuedo_gene_1[order(rownames(psuedo_gene_1)),]
psuedo_gene_1 <- subset(psuedo_gene_1, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('2'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_2 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_2 <- psuedo_gene_2[order(rownames(psuedo_gene_2)),]
psuedo_gene_2 <- subset(psuedo_gene_2, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('3'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_3 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_3 <- psuedo_gene_3[order(rownames(psuedo_gene_3)),]
psuedo_gene_3 <- subset(psuedo_gene_3, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('4'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_4 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_4 <- psuedo_gene_4[order(rownames(psuedo_gene_4)),]
psuedo_gene_4 <- subset(psuedo_gene_4, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('5'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_5 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_5 <- psuedo_gene_5[order(rownames(psuedo_gene_5)),]
psuedo_gene_5 <- subset(psuedo_gene_5, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('6'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_6 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_6 <- psuedo_gene_6[order(rownames(psuedo_gene_6)),]
psuedo_gene_6 <- subset(psuedo_gene_6, select = c('p_val_adj'))


Sub_cluster <- subset(Gene_de1, idents = c('7'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_7 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_7 <- psuedo_gene_7[order(rownames(psuedo_gene_7)),]
psuedo_gene_7 <- subset(psuedo_gene_7, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('8'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_8 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_8 <- psuedo_gene_8[order(rownames(psuedo_gene_8)),]
psuedo_gene_8 <- subset(psuedo_gene_8, select = c('p_val_adj'))

Sub_cluster <- subset(Gene_de1, idents = c('9'))
Idents(Sub_cluster) <- Sub_cluster$type
psuedo_gene_9 <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
psuedo_gene_9 <- psuedo_gene_9[order(rownames(psuedo_gene_9)),]
psuedo_gene_9 <- subset(psuedo_gene_9, select = c('p_val_adj'))

S100a6 <- c(psuedo_gene_1['S100a6',], psuedo_gene_2['S100a6',], psuedo_gene_3['S100a6',], psuedo_gene_4['S100a6',], psuedo_gene_5['S100a6',], psuedo_gene_6['S100a6',], psuedo_gene_7['S100a6',], psuedo_gene_8['S100a6',], psuedo_gene_9['S100a6',], psuedo_gene_10['S100a6',])
Ly6a <- c(psuedo_gene_1['Ly6a',], psuedo_gene_2['Ly6a',], psuedo_gene_3['Ly6a',], psuedo_gene_4['Ly6a',], psuedo_gene_5['Ly6a',], psuedo_gene_6['Ly6a',], psuedo_gene_7['Ly6a',], psuedo_gene_8['Ly6a',], psuedo_gene_9['Ly6a',], psuedo_gene_10['Ly6a',])
Lgr5 <- c(psuedo_gene_1['Lgr5',], psuedo_gene_2['Lgr5',], psuedo_gene_3['Lgr5',], psuedo_gene_4['Lgr5',], psuedo_gene_5['Lgr5',], psuedo_gene_6['Lgr5',], psuedo_gene_7['Lgr5',], psuedo_gene_8['Lgr5',], psuedo_gene_9['Lgr5',], psuedo_gene_10['Lgr5',])
Cd55 <- c(psuedo_gene_1['Cd55',], psuedo_gene_2['Cd55',], psuedo_gene_3['Cd55',], psuedo_gene_4['Cd55',], psuedo_gene_5['Cd55',], psuedo_gene_6['Cd55',], psuedo_gene_7['Cd55',], psuedo_gene_8['Cd55',], psuedo_gene_9['Cd55',], psuedo_gene_10['Cd55',])
Ascl2 <- c(psuedo_gene_1['Ascl2',], psuedo_gene_2['Ascl2',], psuedo_gene_3['Ascl2',], psuedo_gene_4['Ascl2',], psuedo_gene_5['Ascl2',], psuedo_gene_6['Ascl2',], psuedo_gene_7['Ascl2',], psuedo_gene_8['Ascl2',], psuedo_gene_9['Ascl2',], psuedo_gene_10['Ascl2',])
pseudo_PVal <- data.frame(S100a6, Ly6a, Lgr5, Cd55, Ascl2)
write.csv(pseudo_PVal, 'Organoid_Psuedo_gene_signficance.csv')

#make plots
pseudo_PVal <- read.csv('./data/AA_ARA_Intestine_project/Organoid_Psuedo_gene_signficance.csv')
colnames(pseudo_PVal)<- c('num' ,"S100a6", "Ly6a",  "Lgr5",  "Cd55", "Ascl2")

AA_rownames <- rownames(AA_AA@meta.data)
Control_rownames <- rownames(AA_Control@meta.data)
Cell_type_frame <- data.frame(Gene_de$type,Gene_de$orig.ident )
all <-  cbind(type = Cell_type_frame, Pseudotime = Gene_de$Pseudotime)
colnames(all) <- c("Type",'Replicate',"Pseudotime")

all <- all %>% mutate(quantile = ntile(all$Pseudotime, 10))
quart <- data.frame(var1 = tapply(all$Pseudotime, all$quantile, mean), var2 = tapply(all$Pseudotime, all$quantile, max))

features_pseudo_scores <- c("Fetal_Score", "Granuloma_Score", "Radiation_Score", "Homeostatic_Score", "Regeneration_Score")
for(i in features_pseudo_scores) {
  score_iter = i
  quart$maxi = max(Gene_de@meta.data[,i])
  AA <- subset(all, Type=='AA')
  Control <- subset(all, Type=='Control')
  
  AA <-  cbind(AA, AA_AA@meta.data[,i])
  Control <-  cbind(Control, AA_Control@meta.data[,i])
  
  colnames(AA) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(AA, quantile== i)
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
    geom_jitter(data=AA, aes(x=Pseudotime, y=score), size=.25, alpha=0.1, width=1,stroke = .3, colour = "#d95f02") +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=1) + 
    geom_smooth(data=AA, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=1) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"), position=position_dodge2(0.75), vjust=.5, inherit.aes = FALSE)
  nameoooo <-  'AA_'
  name1 <- score_iter
  Name2 <- '.svg'
  File_name <-  paste(nameoooo, name1, Name2, sep = "")
  ggsave(file= File_name, plot=Pseudotimedot, width=5, height=2.5)
}
all
features_pseudo_genes <- c( "S100a6", "Ascl2", "Lgr5", "Ly6a", "Cd55")
for(i in features_pseudo_genes) {
  score_iter = i
  quart$maxi = max(Gene_de@assays$MAGIC_RNA@data[i, ])
  AA <- subset(all, Type=='AA')
  Control <- subset(all, Type=='Control')
  
  AA <-  cbind(AA, AA_AA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, AA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(AA) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(AA, quantile== i)
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
    ggrastr::rasterise(geom_jitter(data=AA, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=AA, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + theme_cem +
    geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=5),axis.text = element_text(size = 6) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
  nameoooo <-  'AA_'
  name1 <- score_iter
  Name2 <- '.pdf'
  File_name <-  paste(nameoooo, name1, Name2, sep = "")
  ggsave(file= File_name, plot=Pseudotimedot, width = 2.2, height = 1, units="in")
}

{
  #final figure
  #Lgr5
  i = 'Lgr5'
  score_iter = i
  quart$maxi = max(Gene_de@assays$MAGIC_RNA@data[i, ])
  AA <- subset(all, Type=='AA')
  Control <- subset(all, Type=='Control')
  
  AA <-  cbind(AA, AA_AA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, AA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(AA) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(AA, quantile== i)
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
    ggrastr::rasterise(geom_jitter(data=AA, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=AA, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) + theme_cem +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  #Ascl2
  i = 'Ascl2'
  score_iter = i
  quart$maxi = max(Gene_de@assays$MAGIC_RNA@data[i, ])
  AA <- subset(all, Type=='AA')
  Control <- subset(all, Type=='Control')
  
  AA <-  cbind(AA, AA_AA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, AA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(AA) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(AA, quantile== i)
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
    ggrastr::rasterise(geom_jitter(data=AA, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=AA, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) + theme_cem +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  #Ly6a
  i = 'Ascl2'
  score_iter = i
  quart$maxi = max(Gene_de@assays$MAGIC_RNA@data[i, ])
  AA <- subset(all, Type=='AA')
  Control <- subset(all, Type=='Control')
  
  AA <-  cbind(AA, AA_AA@assays$MAGIC_RNA@data[i, ])
  Control <-  cbind(Control, AA_Control@assays$MAGIC_RNA@data[i, ])
  
  colnames(AA) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(Control) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(i in 1:10){
    AAinqunat  <- subset(AA, quantile== i)
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
    ggrastr::rasterise(geom_jitter(data=AA, aes(x=Pseudotime, y=score), size=.15, alpha=0.075, width=1,stroke = .3, colour = "#d95f02"), dpi = 300) +
    # lines
    geom_smooth(data=Control, aes(x=Pseudotime, y=score), fill="black", colour="#1b9e77", size=.5) + 
    geom_smooth(data=AA, aes(x=Pseudotime, y=score), fill="black", colour="#d95f02", size=.5) + theme_cem +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) +
    geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"),size = 1.5 , position=position_dodge2(0.75), vjust=.7, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  plot_grid(Density_diff_all, Pseudotimedot_Lgr5, Pseudotimedot_Ly6a, labels = c('Stem 3','Lgr5', 'Ascl2'), label_size = 6, ncol = 1, align = "v", hjust = -.005, label_y = .95, vjust = -.1) #'All Cells', 'Stem 1', 'Stem 2'
  ggsave(paste0("Organoid_pseudo_gene.pdf"), width = 2.2, height = 3.1, units="in")
}

#PGE2 supplement
#Gene level analysis
Idents(organoid) <- organoid$type
PGE_de <- subset(organoid,  idents =c('Control', 'Pge2'))
PGE_de <- NormalizeData(object = PGE_de ,normalization.method = "LogNormalize", assay = "RNA")
# do some abra kadabra magic

my_levels <- c('Control' ,'Pge2')
PGE_de$type <- factor(x = PGE_de$Type, levels = my_levels)
DefaultAssay(PGE_de) <- "RNA"
PGE_de <- magic(PGE_de)

Prop_table<- prop.table(x = table(PGE_de$Cell_type, PGE_de$type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
#Prop_Table1$Freq <- 1/log10(Prop_Table$Freq)*-1
Prop_Table1$Freq[2:3] <- c(0,0)
my_levels <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','Pge2'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))

ggsave(file="PGE_de_prop.svg", plot=plot, width=10, height=8)
#cluster check
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg",'Cbx3','Larp1','Slc39a1','Hnf4a','Sox4','Mmp7','Dll1','Tff3',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a","Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
DotPlot(PGE_de, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text= element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file= 'PGE2_Haber_check_dotplot.pdf', width=4.5, height=2.6, units="in")


#prop sigs
tail(PGE_de$orig.ident)
prop_sig <- table(PGE_de$Cell_type, PGE_de$orig.ident)
Prop_Sig_Table <- as.data.frame.matrix(prop_sig)
num <- length(rownames(Prop_Sig_Table))

for(i in 1:num){
  dat = data.frame(sampleID = c(colnames(Prop_Sig_Table)), current_cluster=c(Prop_Sig_Table[i,]$P5, Prop_Sig_Table[i,]$C5, Prop_Sig_Table[i,]$Control4), all_other_clusters=c( sum(Prop_Sig_Table$P5) - Prop_Sig_Table[i,]$P5, sum(Prop_Sig_Table$C5) - Prop_Sig_Table[i,]$C5, sum(Prop_Sig_Table$Control4) - Prop_Sig_Table[i,]$Control4), Treatment = c("Pge2", "Control", "Control"))
  glm_output = glm(data = dat, cbind(current_cluster, all_other_clusters) ~ Treatment, family = "binomial")
  print(coef(summary(glm_output))[,'Pr(>|z|)'])
}

gene.level.list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'Cd55')
All <- VlnPlot(PGE_de, group.by = 'Cell_type' ,  split.by = "type", features = gene.level.list, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'),log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
ggsave(file = paste0('PGE_gene_vln.pdf'), plot=All,  width=4, height=3, units="in")

Lgr5 <- VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Lgr5'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Ascl2 <- VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ascl2'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
S100a6 <- VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('S100a6'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)  + xlab('')
Ly6a <- VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Ly6a'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')
Cd55 <- VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Cd55'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')

ggsave(file = 'Lgr5_pge.svg', plot=Lgr5, width=10, height=10)
ggsave(file = 'Ascl2_pge.svg', plot=Ascl2, width=10, height=10)
ggsave(file = 'S100a6_pge.svg', plot=S100a6, width=10, height=10)
ggsave(file = 'Ly6a_pge.svg', plot=Ly6a, width=10, height=10)
ggsave(file = 'Cd55_pge.svg', plot=Cd55, width=10, height=10)

VlnPlot(PGE_de, group.by = 'Cell_type' , ncol = 1, split.by = "type", features = c('Lgr5'), pt.size = 0, assay = "MAGIC_RNA", log = TRUE, cols = c('#1b9e77' ,'#7570b3'), split.plot = TRUE, combine = TRUE) + theme(legend.position = 'none') + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0) + xlab('')


# PGE2 gene heatmaps
names <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )

Idents(Organoid.obj) <- Organoid.obj$Cell_type
for(i in names){
  Sub_cluster <- subset(Organoid.obj,  idents = i)
  Idents(Sub_cluster) <- Sub_cluster$type
  Gene_pval <- FindMarkers(Sub_cluster, ident.1 = "Organoid", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
  write.csv(Gene_pval, file = paste0(i, '_Organoid_gene_pval.csv'))
}
Enterocyte_Proximal_Organoid_gene_pval
Enterocyte_Progenitor_Organoid_gene_pval
Enteroendocrine_Organoid_gene_pval
Goblet_Organoid_gene_pval
Paneth_Organoid_gene_pval
Stem_1_Organoid_gene_pval
Transit_Amplifying_Organoid_gene_pval
Tuft_Organoid_gene_pval

rownames(Enterocyte_Organoid_gene_pval) <- Enterocyte_Organoid_gene_pval$X1
rownames(Enterocyte_Progenitor_Organoid_gene_pval) <- Enterocyte_Progenitor_Organoid_gene_pval$X1
rownames(Enteroendocrine_Organoid_gene_pval) <- Enteroendocrine_Organoid_gene_pval$X1
rownames(Goblet_Organoid_gene_pval) <- Goblet_Organoid_gene_pval$X1
rownames(Paneth_Organoid_gene_pval) <- Paneth_Organoid_gene_pval$X1
rownames(Stem_1_Organoid_gene_pval) <- Stem_1_Organoid_gene_pval$X1
rownames(Transit_Amplifying_Organoid_gene_pval) <- Transit_Amplifying_Organoid_gene_pval$X1
rownames(Tuft_Organoid_gene_pval) <- Tuft_Organoid_gene_pval$X1
rownames(Secretory_Progenitor_Organoid_gene_pval) <- Secretory_Progenitor_Organoid_gene_pval$X1
rownames(Stem_Like_Progenitor_Organoid_gene_pval) <- Stem_Like_Progenitor_Organoid_gene_pval$X1
Secretory_Progenitor_Organoid_gene_pval

Enterocyte <- ''
Enterocyte_Progenitor <- ''
Enteroendocrine <- ''
Goblet <- ''
Paneth <- ''
Stem_1 <- ''
Transit_Amplifying <- ''
Tuft <- ''
Secretory_Progenitor <- ''
Stem_Like_Progenitor <- ''
gene_list <-  c('S100a6','Lgr5','Ascl2','Ly6a','Cd55','Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')

for(i in gene_list){
  Enterocyte <- c(Enterocyte, Enterocyte_Organoid_gene_pval[i,]$p_val_adj)
  Enterocyte_Progenitor <- c(Enterocyte_Progenitor, Enterocyte_Progenitor_Organoid_gene_pval[i,]$p_val_adj)
  Enteroendocrine <- c(Enteroendocrine, Enteroendocrine_Organoid_gene_pval[i,]$p_val_adj)
  Goblet <- c(Goblet, Goblet_Organoid_gene_pval[i,]$p_val_adj)
  Paneth <- c(Paneth, Paneth_Organoid_gene_pval[i,]$p_val_adj)
  Stem_1 <- c(Stem_1, Stem_1_Organoid_gene_pval[i,]$p_val_adj)
  Transit_Amplifying <- c(Transit_Amplifying, Transit_Amplifying_Organoid_gene_pval[i,]$p_val_adj)
  Tuft <- c(Tuft, Tuft_Organoid_gene_pval[i,]$p_val_adj)
  Secretory_Progenitor <- c(Secretory_Progenitor, Secretory_Progenitor_Organoid_gene_pval[i,]$p_val_adj)
  Stem_Like_Progenitor <- c(Stem_Like_Progenitor, Stem_Like_Progenitor_Organoid_gene_pval[i,]$p_val_adj)
}

Enterocyte <- Enterocyte[2:29]
Enterocyte_Progenitor <- Enterocyte_Progenitor[2:29]
Enteroendocrine <- Enteroendocrine[2:29]
Goblet <- Goblet[2:29]
Paneth <- Paneth[2:29]
Stem_1 <- Stem_1[2:29]
Transit_Amplifying <- Transit_Amplifying[2:29]
Tuft <- Tuft[2:29]
Secretory_Progenitor <- Secretory_Progenitor[2:29]
Stem_Like_Progenitor <- Stem_Like_Progenitor[2:29]

Stem_2 <- c( 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Stem_3 <- c( 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

PGE2_heatmap <- data.frame(Stem_1, Stem_2, Stem_3, Stem_Like_Progenitor, Secretory_Progenitor, Transit_Amplifying, Enterocyte_Progenitor, Enterocyte, Enteroendocrine, Goblet, Tuft, Paneth, stringsAsFactors=FALSE)

rownames(PGE2_heatmap) <- gene_list
PGE2_heatmap[is.na(PGE2_heatmap)] <- 0
PGE2_heatmap <- data.frame(sapply(PGE2_heatmap, as.numeric))

PGE2_heatmap1 <- data.frame(PGE2_heatmap)
rownames(PGE2_heatmap1) <- gene_list
colnames(PGE2_heatmap1) <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
write.csv(PGE2_heatmap1, 'Organoid_MAST_PVALS.csv')

dt2 <- PGE2_heatmap1 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
geneColorsDiff = brewer.pal(n=5, "RdYlBu")
ggplot(dt2, aes(x = rowname, y = colname, fill = value)) + scale_fill_viridis_c(option = "magma") + geom_tile(color="white", size=0.25) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Gene") + ylab("Cell Type") + labs(fill = "Log FC") + scale_y_discrete(limits = c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" ))

#housekeeping
setwd('~/')
#saveRDS(Gene_de, file = "./data/AA_ARA_Intestine_project/Suerat_objects/Gene_de_Sobj.rds")

Gene_de <- readRDS('./data/AA_ARA_Intestine_project/Suerat_objects/Gene_de_Sobj.rds', refhook = NULL)

Idents(Gene_de) <- Gene_de$Cell_type
current.cluster.ids <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
new.cluster.ids <- c("Homeostatic Stem", "Revival Stem 1", "Revival Stem 2", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Gene_de@active.ident <- plyr::mapvalues(x = Gene_de@active.ident, from = new.cluster.ids, to = current.cluster.ids)
DimPlot(Gene_de)
Gene_de$Cell_type <- Idents(Gene_de)
DimPlot(Gene_de, group.by = "Cell_type")

Idents(Gene_de) <- Gene_de$Cell_type
current.cluster.ids <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
new.cluster.ids <- c("Homeostatic Stem", "Revival Stem 1", "Revival Stem 2", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Gene_de@active.ident <- plyr::mapvalues(x = Gene_de@active.ident, from = new.cluster.ids, to = current.cluster.ids)
DimPlot(Gene_de)
Gene_de$Cell_type <- Idents(Gene_de)
DimPlot(Gene_de, group.by = "Cell_type")

All_Genes <- Gene_de@assays$RNA@data@Dimnames[[1]]
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

writeGmtPathways(Scores, 'Beyaz_AA_final_Gene_de.gmt')

identities <- levels(organoid@ident)
my_color_palette <- hue_pal()(length(identities))


