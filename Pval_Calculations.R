#Vyom Shah
#MAST test
#organoid

names <- c('Stem 1',"Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Idents(Gene_de) <- Gene_de$Cell_type
for(i in names){
  Sub_cluster <- subset(Gene_de,  idents = i)
  Idents(Sub_cluster) <- Sub_cluster$type
  Gene_pval <- FindMarkers(Sub_cluster, ident.1 = "AA", ident.2 = "Control", test.use = "MAST", logfc.threshold = 0, min.pct = 0, assay = 'RNA')
  write.csv(Gene_pval, file = paste0(i, '_Organoid_gene_pval.csv'))
}

#wilcox test all
#Organoid Pval Calculations

names <- c( "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
scores <- c( 'Radiation_Score','Fetal_Score','Granuloma_Score','Regeneration_Score','Homeostatic_Score','YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','ASCL_Targets_Score','ECM_Induced_Score')
Idents(Gene_de) <- Gene_de$Cell_type

# scores
Sub_cluster <- subset(Gene_de,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('AA'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))
Radiation_Score <- c(wilcox.test(test_AA@meta.data$Radiation_Score, test_contr@meta.data$Radiation_Score, alternative = "two.sided")$p.value, 0 , 0)
Fetal_Score <- c(wilcox.test(test_AA@meta.data$Fetal_Score, test_contr@meta.data$Fetal_Score, alternative = "two.sided")$p.value, 0, 0)
Granuloma_Score <- c(wilcox.test(test_AA@meta.data$Granuloma_Score, test_contr@meta.data$Granuloma_Score, alternative = "two.sided")$p.value, 0 , 0)
Regeneration_Score <- c(wilcox.test(test_AA@meta.data$Regeneration_Score, test_contr@meta.data$Regeneration_Score, alternative = "two.sided")$p.value, 0 , 0)
Homeostatic_Score <- c(wilcox.test(test_AA@meta.data$Homeostatic_Score, test_contr@meta.data$Homeostatic_Score, alternative = "two.sided")$p.value, 0 , 0)
YAP_Targets_Score <- c(wilcox.test(test_AA@meta.data$YAP_Targets_Score, test_contr@meta.data$YAP_Targets_Score, alternative = "two.sided")$p.value, 0 , 0)
BCatenin_Enhanced_Score <- c(wilcox.test(test_AA@meta.data$BCatenin_Enhanced_Score, test_contr@meta.data$BCatenin_Enhanced_Score, alternative = "two.sided")$p.value, 0 , 0)
CREB_Targets_Score <- c(wilcox.test(test_AA@meta.data$CREB_Targets_Score, test_contr@meta.data$CREB_Targets_Score, alternative = "two.sided")$p.value, 0 , 0)
ASCL_Targets_Score <- c(wilcox.test(test_AA@meta.data$ASCL_Targets_Score, test_contr@meta.data$ASCL_Targets_Score, alternative = "two.sided")$p.value, 0 , 0)
ECM_Induced_Score <- c(wilcox.test(test_AA@meta.data$ECM_Induced_Score, test_contr@meta.data$ECM_Induced_Score, alternative = "two.sided")$p.value, 0 , 0)

for(i in names) {
  Sub_cluster <- subset(Gene_de,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('AA'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  Radiation_Score <- c(Radiation_Score, wilcox.test(test_AA@meta.data$Radiation_Score, test_contr@meta.data$Radiation_Score, alternative = "two.sided")$p.value)
  Fetal_Score <- c(Fetal_Score, wilcox.test(test_AA@meta.data$Fetal_Score, test_contr@meta.data$Fetal_Score, alternative = "two.sided")$p.value)
  Granuloma_Score <- c(Granuloma_Score, wilcox.test(test_AA@meta.data$Granuloma_Score, test_contr@meta.data$Granuloma_Score, alternative = "two.sided")$p.value)
  Regeneration_Score <- c(Regeneration_Score, wilcox.test(test_AA@meta.data$Regeneration_Score, test_contr@meta.data$Regeneration_Score, alternative = "two.sided")$p.value)
  Homeostatic_Score <- c(Homeostatic_Score, wilcox.test(test_AA@meta.data$Homeostatic_Score, test_contr@meta.data$Homeostatic_Score, alternative = "two.sided")$p.value)
  YAP_Targets_Score <- c(YAP_Targets_Score,wilcox.test(test_AA@meta.data$YAP_Targets_Score, test_contr@meta.data$YAP_Targets_Score, alternative = "two.sided")$p.value)
  BCatenin_Enhanced_Score <- c(BCatenin_Enhanced_Score, wilcox.test(test_AA@meta.data$BCatenin_Enhanced_Score, test_contr@meta.data$BCatenin_Enhanced_Score, alternative = "two.sided")$p.value)
  CREB_Targets_Score <- c(CREB_Targets_Score, wilcox.test(test_AA@meta.data$CREB_Targets_Score, test_contr@meta.data$CREB_Targets_Score, alternative = "two.sided")$p.value)
  ASCL_Targets_Score <- c(ASCL_Targets_Score, wilcox.test(test_AA@meta.data$ASCL_Targets_Score, test_contr@meta.data$ASCL_Targets_Score, alternative = "two.sided")$p.value)
  ECM_Induced_Score <- c(ECM_Induced_Score, wilcox.test(test_AA@meta.data$ECM_Induced_Score, test_contr@meta.data$ECM_Induced_Score, alternative = "two.sided")$p.value)
}

Organoid_PVal <- data.frame(Radiation_Score, Fetal_Score, Granuloma_Score, Regeneration_Score, Homeostatic_Score, YAP_Targets_Score, BCatenin_Enhanced_Score, CREB_Targets_Score, ASCL_Targets_Score, ECM_Induced_Score)

for(i in scores) {
  Organoid_PVal[i] <- p.adjust(as.matrix(Organoid_PVal[i]), method = 'BH')
}

write.csv(Organoid_PVal,'Organoid_Score_PVal.csv')

Cell_Type <- c('Stem 1','Stem 2','Stem 3',"Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Organoid_PVal <- data.frame(Radiation_Score, Fetal_Score, Granuloma_Score, Regeneration_Score, Homeostatic_Score, YAP_Targets_Score, BCatenin_Enhanced_Score, CREB_Targets_Score, ASCL_Targets_Score, ECM_Induced_Score)
rownames(Organoid_PVal) <- Cell_Type

write.csv(Organoid_PVal,'Organoid_Gene_lvl_PVal.csv')
Organoid_PVal <- read_csv("Organoid_Gene_lvl_PVal.csv")

Organoid_PVal <- matrix(p.adjust(as.vector(as.matrix(Organoid_PVal)), method='BH'),ncol=10)
colnames(Organoid_PVal) <-  scores
rownames(Organoid_PVal) <- Cell_Type

write.csv(Organoid_PVal,'Organoid_score_lvl_PVal.csv')


#Gene level
Pval_gene_list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'CD55')

Sub_cluster <- subset(Gene_de,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('AA'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))
DefaultAssay(test_contr) <- "MAGIC_RNA"
DefaultAssay(test_AA) <- "MAGIC_RNA"

Lgr5 <- c(wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value, 0, 0)
Ascl2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value, 0, 0)
S100a6 <- c(wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value, 0, 0)
Ly6a <- c(wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value, 0, 0)
CD55 <- c(wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value, 0, 0)

for(i in names) {
  Sub_cluster <- subset(Gene_de,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('AA'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  DefaultAssay(test_contr) <- "MAGIC_RNA"
  DefaultAssay(test_AA) <- "MAGIC_RNA"
  
  Lgr5 <- c(Lgr5, wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value)
  Ascl2 <- c(Ascl2, wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value)
  S100a6 <- c(S100a6, wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value)
  Ly6a <- c(Ly6a, wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value)
  CD55 <- c(CD55, wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value)
}
Cell_Type <- c('Stem 1','Stem 2','Stem 3',"Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Gene_Level_PVal <- data.frame(Lgr5, Ascl2, S100a6, Ly6a, CD55)
rownames(Gene_Level_PVal) <- Cell_Type
Gene_Level_PVal$Ly6a[is.nan(Gene_Level_PVal$Ly6a)] <- NA

Gene_Level_PVal1 <- matrix(p.adjust(as.vector(as.matrix(Gene_Level_PVal)), method='BH'),ncol=5)
colnames(Gene_Level_PVal1) <-  Pval_gene_list
rownames(Gene_Level_PVal1) <- Cell_Type

write.csv(Gene_Level_PVal1,'Organoid_Gene_lvl_PVal.csv')


#ARAsco Pval Calculations

names <- c('Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
scores <- c( 'Radiation_Score','Fetal_Score','Granuloma_Score','Regeneration_Score','Homeostatic_Score','YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','ASCL_Targets_Score','ECM_Induced_Score')
Idents(arasco.obj) <- arasco.obj$Cell_type

# scores
Sub_cluster <- subset(arasco.obj,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))
Radiation_Score <- c(wilcox.test(test_AA@meta.data$Radiation_Score, test_contr@meta.data$Radiation_Score, alternative = "two.sided")$p.value, 0)
Fetal_Score <- c(wilcox.test(test_AA@meta.data$Fetal_Score, test_contr@meta.data$Fetal_Score, alternative = "two.sided")$p.value, 0)
Granuloma_Score <- c(wilcox.test(test_AA@meta.data$Granuloma_Score, test_contr@meta.data$Granuloma_Score, alternative = "two.sided")$p.value, 0)
Regeneration_Score <- c(wilcox.test(test_AA@meta.data$Regeneration_Score, test_contr@meta.data$Regeneration_Score, alternative = "two.sided")$p.value, 0)
Homeostatic_Score <- c(wilcox.test(test_AA@meta.data$Homeostatic_Score, test_contr@meta.data$Homeostatic_Score, alternative = "two.sided")$p.value, 0)
YAP_Targets_Score <- c(wilcox.test(test_AA@meta.data$YAP_Targets_Score, test_contr@meta.data$YAP_Targets_Score, alternative = "two.sided")$p.value, 0)
BCatenin_Enhanced_Score <- c(wilcox.test(test_AA@meta.data$BCatenin_Enhanced_Score, test_contr@meta.data$BCatenin_Enhanced_Score, alternative = "two.sided")$p.value, 0)
CREB_Targets_Score <- c(wilcox.test(test_AA@meta.data$CREB_Targets_Score, test_contr@meta.data$CREB_Targets_Score, alternative = "two.sided")$p.value, 0)
ASCL_Targets_Score <- c(wilcox.test(test_AA@meta.data$ASCL_Targets_Score, test_contr@meta.data$ASCL_Targets_Score, alternative = "two.sided")$p.value, 0)
ECM_Induced_Score <- c(wilcox.test(test_AA@meta.data$ECM_Induced_Score, test_contr@meta.data$ECM_Induced_Score, alternative = "two.sided")$p.value, 0)

for(i in names) {
  Sub_cluster <- subset(arasco.obj,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  Radiation_Score <- c(Radiation_Score, wilcox.test(test_AA@meta.data$Radiation_Score, test_contr@meta.data$Radiation_Score, alternative = "two.sided")$p.value)
  Fetal_Score <- c(Fetal_Score, wilcox.test(test_AA@meta.data$Fetal_Score, test_contr@meta.data$Fetal_Score, alternative = "two.sided")$p.value)
  Granuloma_Score <- c(Granuloma_Score, wilcox.test(test_AA@meta.data$Granuloma_Score, test_contr@meta.data$Granuloma_Score, alternative = "two.sided")$p.value)
  Regeneration_Score <- c(Regeneration_Score, wilcox.test(test_AA@meta.data$Regeneration_Score, test_contr@meta.data$Regeneration_Score, alternative = "two.sided")$p.value)
  Homeostatic_Score <- c(Homeostatic_Score, wilcox.test(test_AA@meta.data$Homeostatic_Score, test_contr@meta.data$Homeostatic_Score, alternative = "two.sided")$p.value)
  YAP_Targets_Score <- c(YAP_Targets_Score,wilcox.test(test_AA@meta.data$YAP_Targets_Score, test_contr@meta.data$YAP_Targets_Score, alternative = "two.sided")$p.value)
  BCatenin_Enhanced_Score <- c(BCatenin_Enhanced_Score, wilcox.test(test_AA@meta.data$BCatenin_Enhanced_Score, test_contr@meta.data$BCatenin_Enhanced_Score, alternative = "two.sided")$p.value)
  CREB_Targets_Score <- c(CREB_Targets_Score, wilcox.test(test_AA@meta.data$CREB_Targets_Score, test_contr@meta.data$CREB_Targets_Score, alternative = "two.sided")$p.value)
  ASCL_Targets_Score <- c(ASCL_Targets_Score, wilcox.test(test_AA@meta.data$ASCL_Targets_Score, test_contr@meta.data$ASCL_Targets_Score, alternative = "two.sided")$p.value)
  ECM_Induced_Score <- c(ECM_Induced_Score, wilcox.test(test_AA@meta.data$ECM_Induced_Score, test_contr@meta.data$ECM_Induced_Score, alternative = "two.sided")$p.value)
}

Cell_Type <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
ARAsco_PVal <- data.frame(Radiation_Score, Fetal_Score, Granuloma_Score, Regeneration_Score, Homeostatic_Score, YAP_Targets_Score, BCatenin_Enhanced_Score, CREB_Targets_Score, ASCL_Targets_Score, ECM_Induced_Score)
rownames(ARAsco_PVal) <- Cell_Type

ARAsco_PVal <- matrix(p.adjust(as.vector(as.matrix(ARAsco_PVal)), method='BH'),ncol=10)
colnames(ARAsco_PVal) <-  scores
rownames(ARAsco_PVal) <- Cell_Type
write.csv(ARAsco_PVal,'ARAsco_Score_PVal.csv')

#Gene level
Pval_gene_list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'CD55')

Sub_cluster <- subset(arasco.obj,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))
DefaultAssay(test_contr) <- "MAGIC_RNA"
DefaultAssay(test_AA) <- "MAGIC_RNA"

Lgr5 <- c(wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value, 0)
Ascl2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value, 0)
S100a6 <- c(wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value, 0)
Ly6a <- c(wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value, 0)
CD55 <- c(wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value, 0)

for(i in names) {
  Sub_cluster <- subset(arasco.obj,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  DefaultAssay(test_contr) <- "MAGIC_RNA"
  DefaultAssay(test_AA) <- "MAGIC_RNA"
  
  Lgr5 <- c(Lgr5, wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value)
  Ascl2 <- c(Ascl2, wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value)
  S100a6 <- c(S100a6, wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value)
  Ly6a <- c(Ly6a, wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value)
  CD55 <- c(CD55, wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value)
}

Cell_Type <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
Gene_Level_PVal <- data.frame(Lgr5, Ascl2, S100a6, Ly6a, CD55)
rownames(Gene_Level_PVal) <- Cell_Type
Gene_Level_PVal$Ly6a[is.nan(Gene_Level_PVal$Ly6a)] <- NA

Gene_Level_PVal1 <- matrix(p.adjust(as.vector(as.matrix(Gene_Level_PVal)), method='BH'),ncol=5)
colnames(Gene_Level_PVal1) <-  Pval_gene_list
rownames(Gene_Level_PVal1) <- Cell_Type

write.csv(Gene_Level_PVal,'ARAsco_Gene_lvl_PVal.csv')

#PGE2 Genes
names <- c('Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
PGE2_genes <-  c('Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')
Idents(arasco.obj) = arasco.obj$Cell_type
Sub_cluster <- subset(arasco.obj,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))

Ptger1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger1")$Ptger1, FetchData(object = test_contr, vars = "Ptger1")$Ptger1, alternative = "two.sided")$p.value, 0)
Ptger2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger2")$Ptger2, FetchData(object = test_contr, vars = "Ptger2")$Ptger2, alternative = "two.sided")$p.value, 0)
Ptger3 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger3")$Ptger3, FetchData(object = test_contr, vars = "Ptger3")$Ptger3, alternative = "two.sided")$p.value, 0)
Ptger4 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger4")$Ptger4, FetchData(object = test_contr, vars = "Ptger4")$Ptger4, alternative = "two.sided")$p.value, 0)
Ptges <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptges")$Ptges, FetchData(object = test_contr, vars = "Ptges")$Ptges, alternative = "two.sided")$p.value, 0)
Ptges2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptges2")$Ptges2, FetchData(object = test_contr, vars = "Ptges2")$Ptges2, alternative = "two.sided")$p.value, 0)
Ptgs1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgs1")$Ptgs1, FetchData(object = test_contr, vars = "Ptgs1")$Ptgs1, alternative = "two.sided")$p.value, 0)
Alox12b <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox12b")$Alox12b, FetchData(object = test_contr, vars = "Alox12b")$Alox12b, alternative = "two.sided")$p.value, 0)
Ptgs2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgs2")$Ptgs2, FetchData(object = test_contr, vars = "Ptgs2")$Ptgs2, alternative = "two.sided")$p.value, 0)
Ptgis <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgis")$Ptgis, FetchData(object = test_contr, vars = "Ptgis")$Ptgis, alternative = "two.sided")$p.value, 0)
Pla2g2e <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g2e")$Pla2g2e, FetchData(object = test_contr, vars = "Pla2g2e")$Pla2g2e, alternative = "two.sided")$p.value, 0)
Ggt1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ggt1")$Ggt1, FetchData(object = test_contr, vars = "Ggt1")$Ggt1, alternative = "two.sided")$p.value, 0)
Gpx3 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx3")$Gpx3, FetchData(object = test_contr, vars = "Gpx3")$Gpx3, alternative = "two.sided")$p.value, 0)
Pla2g10 <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g10")$Pla2g10, FetchData(object = test_contr, vars = "Pla2g10")$Pla2g10, alternative = "two.sided")$p.value, 0)
Pla2g12b <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g12b")$Pla2g12b, FetchData(object = test_contr, vars = "Pla2g12b")$Pla2g12b, alternative = "two.sided")$p.value, 0)
Ephx2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ephx2")$Ephx2, FetchData(object = test_contr, vars = "Ephx2")$Ephx2, alternative = "two.sided")$p.value, 0)
Alox12 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox12")$Alox12, FetchData(object = test_contr, vars = "Alox12")$Alox12, alternative = "two.sided")$p.value, 0)
Gpx1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx1")$Gpx1, FetchData(object = test_contr, vars = "Gpx1")$Gpx1, alternative = "two.sided")$p.value, 0)
Gpx4 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx4")$Gpx4, FetchData(object = test_contr, vars = "Gpx4")$Gpx4, alternative = "two.sided")$p.value, 0)
Gpx2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx2")$Gpx2, FetchData(object = test_contr, vars = "Gpx2")$Gpx2, alternative = "two.sided")$p.value, 0)
Alox5 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox5")$Alox5, FetchData(object = test_contr, vars = "Alox5")$Alox5, alternative = "two.sided")$p.value, 0)
Alox15 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox15")$Alox15, FetchData(object = test_contr, vars = "Alox15")$Alox15, alternative = "two.sided")$p.value, 0)
Pla2g2d <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g2d")$Pla2g2d, FetchData(object = test_contr, vars = "Pla2g2d")$Pla2g2d, alternative = "two.sided")$p.value, 0)

for(i in names) {
  Sub_cluster <- subset(arasco.obj,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('Arasco'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  
  Ptger1 <- c(Ptger1, wilcox.test(FetchData(object = test_AA, vars = "Ptger1")$Ptger1, FetchData(object = test_contr, vars = "Ptger1")$Ptger1, alternative = "two.sided")$p.value)
  Ptger2 <- c(Ptger2, wilcox.test(FetchData(object = test_AA, vars = "Ptger2")$Ptger2, FetchData(object = test_contr, vars = "Ptger2")$Ptger2, alternative = "two.sided")$p.value)
  Ptger3 <- c(Ptger3, wilcox.test(FetchData(object = test_AA, vars = "Ptger3")$Ptger3, FetchData(object = test_contr, vars = "Ptger3")$Ptger3, alternative = "two.sided")$p.value)
  Ptger4 <- c(Ptger4, wilcox.test(FetchData(object = test_AA, vars = "Ptger4")$Ptger4, FetchData(object = test_contr, vars = "Ptger4")$Ptger4, alternative = "two.sided")$p.value)
  Ptges <- c(Ptges, wilcox.test(FetchData(object = test_AA, vars = "Ptges")$Ptges, FetchData(object = test_contr, vars = "Ptges")$Ptges, alternative = "two.sided")$p.value)
  Ptges2 <- c(Ptges2, wilcox.test(FetchData(object = test_AA, vars = "Ptges2")$Ptges2, FetchData(object = test_contr, vars = "Ptges2")$Ptges2, alternative = "two.sided")$p.value)
  Ptgs1 <- c(Ptgs1, wilcox.test(FetchData(object = test_AA, vars = "Ptgs1")$Ptgs1, FetchData(object = test_contr, vars = "Ptgs1")$Ptgs1, alternative = "two.sided")$p.value)
  Alox12b <- c(Alox12b, wilcox.test(FetchData(object = test_AA, vars = "Alox12b")$Alox12b, FetchData(object = test_contr, vars = "Alox12b")$Alox12b, alternative = "two.sided")$p.value)
  Ptgs2 <- c(Ptgs2, wilcox.test(FetchData(object = test_AA, vars = "Ptgs2")$Ptgs2, FetchData(object = test_contr, vars = "Ptgs2")$Ptgs2, alternative = "two.sided")$p.value)
  Ptgis <- c(Ptgis, wilcox.test(FetchData(object = test_AA, vars = "Ptgis")$Ptgis, FetchData(object = test_contr, vars = "Ptgis")$Ptgis, alternative = "two.sided")$p.value)
  Pla2g2e <- c(Pla2g2e, wilcox.test(FetchData(object = test_AA, vars = "Pla2g2e")$Pla2g2e, FetchData(object = test_contr, vars = "Pla2g2e")$Pla2g2e, alternative = "two.sided")$p.value)
  Ggt1 <- c(Ggt1, wilcox.test(FetchData(object = test_AA, vars = "Ggt1")$Ggt1, FetchData(object = test_contr, vars = "Ggt1")$Ggt1, alternative = "two.sided")$p.value)
  Gpx3 <- c(Gpx3, wilcox.test(FetchData(object = test_AA, vars = "Gpx3")$Gpx3, FetchData(object = test_contr, vars = "Gpx3")$Gpx3, alternative = "two.sided")$p.value)
  Pla2g10 <- c(Pla2g10, wilcox.test(FetchData(object = test_AA, vars = "Pla2g10")$Pla2g10, FetchData(object = test_contr, vars = "Pla2g10")$Pla2g10, alternative = "two.sided")$p.value)
  Pla2g12b <- c(Pla2g12b, wilcox.test(FetchData(object = test_AA, vars = "Pla2g12b")$Pla2g12b, FetchData(object = test_contr, vars = "Pla2g12b")$Pla2g12b, alternative = "two.sided")$p.value)
  Ephx2 <- c(Ephx2, wilcox.test(FetchData(object = test_AA, vars = "Ephx2")$Ephx2, FetchData(object = test_contr, vars = "Ephx2")$Ephx2, alternative = "two.sided")$p.value)
  Alox12 <- c(Alox12, wilcox.test(FetchData(object = test_AA, vars = "Alox12")$Alox12, FetchData(object = test_contr, vars = "Alox12")$Alox12, alternative = "two.sided")$p.value)
  Gpx1 <- c(Gpx1, wilcox.test(FetchData(object = test_AA, vars = "Gpx1")$Gpx1, FetchData(object = test_contr, vars = "Gpx1")$Gpx1, alternative = "two.sided")$p.value)
  Gpx4 <- c(Gpx4, wilcox.test(FetchData(object = test_AA, vars = "Gpx4")$Gpx4, FetchData(object = test_contr, vars = "Gpx4")$Gpx4, alternative = "two.sided")$p.value)
  Gpx2 <- c(Gpx2, wilcox.test(FetchData(object = test_AA, vars = "Gpx2")$Gpx2, FetchData(object = test_contr, vars = "Gpx2")$Gpx2, alternative = "two.sided")$p.value)
  Alox5 <- c(Alox5, wilcox.test(FetchData(object = test_AA, vars = "Alox5")$Alox5, FetchData(object = test_contr, vars = "Alox5")$Alox5, alternative = "two.sided")$p.value)
  Alox15 <- c(Alox15, wilcox.test(FetchData(object = test_AA, vars = "Alox15")$Alox15, FetchData(object = test_contr, vars = "Alox15")$Alox15, alternative = "two.sided")$p.value)
  Pla2g2d <- c(Pla2g2d, wilcox.test(FetchData(object = test_AA, vars = "Pla2g2d")$Pla2g2d, FetchData(object = test_contr, vars = "Pla2g2d")$Pla2g2d, alternative = "two.sided")$p.value)
}

Cell_Type <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
Gene_Level_PVal <- data.frame(Ptger1, Ptger2, Ptger3, Ptger4, Ptges, Ptges2, Ptgs1, Alox12b, Ptgs2, Ptgis, Pla2g2e, Ggt1, Gpx3, Pla2g10, Pla2g12b, Ephx2, Alox12, Gpx1, Gpx4, Gpx2, Alox5, Alox15, Pla2g2d)
rownames(Gene_Level_PVal) <- Cell_Type
for(i in colnames(Gene_Level_PVal)){
  Gene_Level_PVal[[i]][is.nan(Gene_Level_PVal[[i]])] <- NA
}

Gene_Level_PVal1 <- matrix(p.adjust(as.vector(as.matrix(Gene_Level_PVal)), method='BH'),ncol=23)
colnames(Gene_Level_PVal1) <-  PGE2_genes
rownames(Gene_Level_PVal1) <- Cell_Type

write.csv(Gene_Level_PVal,'ARAsco_PGE2_Gene_lvl_PVal.csv')


# Organoid
names <- c( "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
PGE2_genes <-  c('Ptger1','Ptger2','Ptger3','Ptger4','Ptges','Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d')
Idents(Gene_de) = Gene_de$Cell_type
Sub_cluster <- subset(Gene_de,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('AA'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))

Ptger1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger1")$Ptger1, FetchData(object = test_contr, vars = "Ptger1")$Ptger1, alternative = "two.sided")$p.value, 0, 0)
Ptger2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger2")$Ptger2, FetchData(object = test_contr, vars = "Ptger2")$Ptger2, alternative = "two.sided")$p.value, 0, 0)
Ptger3 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger3")$Ptger3, FetchData(object = test_contr, vars = "Ptger3")$Ptger3, alternative = "two.sided")$p.value, 0, 0)
Ptger4 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptger4")$Ptger4, FetchData(object = test_contr, vars = "Ptger4")$Ptger4, alternative = "two.sided")$p.value, 0, 0)
Ptges <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptges")$Ptges, FetchData(object = test_contr, vars = "Ptges")$Ptges, alternative = "two.sided")$p.value, 0, 0)
Ptges2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptges2")$Ptges2, FetchData(object = test_contr, vars = "Ptges2")$Ptges2, alternative = "two.sided")$p.value, 0, 0)
Ptgs1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgs1")$Ptgs1, FetchData(object = test_contr, vars = "Ptgs1")$Ptgs1, alternative = "two.sided")$p.value, 0, 0)
Alox12b <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox12b")$Alox12b, FetchData(object = test_contr, vars = "Alox12b")$Alox12b, alternative = "two.sided")$p.value, 0, 0)
Ptgs2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgs2")$Ptgs2, FetchData(object = test_contr, vars = "Ptgs2")$Ptgs2, alternative = "two.sided")$p.value, 0, 0)
Ptgis <- c(wilcox.test(FetchData(object = test_AA, vars = "Ptgis")$Ptgis, FetchData(object = test_contr, vars = "Ptgis")$Ptgis, alternative = "two.sided")$p.value, 0, 0)
Pla2g2e <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g2e")$Pla2g2e, FetchData(object = test_contr, vars = "Pla2g2e")$Pla2g2e, alternative = "two.sided")$p.value, 0, 0)
Ggt1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ggt1")$Ggt1, FetchData(object = test_contr, vars = "Ggt1")$Ggt1, alternative = "two.sided")$p.value, 0, 0)
Gpx3 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx3")$Gpx3, FetchData(object = test_contr, vars = "Gpx3")$Gpx3, alternative = "two.sided")$p.value, 0, 0)
Pla2g10 <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g10")$Pla2g10, FetchData(object = test_contr, vars = "Pla2g10")$Pla2g10, alternative = "two.sided")$p.value, 0, 0)
Pla2g12b <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g12b")$Pla2g12b, FetchData(object = test_contr, vars = "Pla2g12b")$Pla2g12b, alternative = "two.sided")$p.value, 0, 0)
Ephx2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ephx2")$Ephx2, FetchData(object = test_contr, vars = "Ephx2")$Ephx2, alternative = "two.sided")$p.value, 0, 0)
Alox12 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox12")$Alox12, FetchData(object = test_contr, vars = "Alox12")$Alox12, alternative = "two.sided")$p.value, 0, 0)
Gpx1 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx1")$Gpx1, FetchData(object = test_contr, vars = "Gpx1")$Gpx1, alternative = "two.sided")$p.value, 0, 0)
Gpx4 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx4")$Gpx4, FetchData(object = test_contr, vars = "Gpx4")$Gpx4, alternative = "two.sided")$p.value, 0, 0)
Gpx2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Gpx2")$Gpx2, FetchData(object = test_contr, vars = "Gpx2")$Gpx2, alternative = "two.sided")$p.value, 0, 0)
Alox5 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox5")$Alox5, FetchData(object = test_contr, vars = "Alox5")$Alox5, alternative = "two.sided")$p.value, 0, 0)
Alox15 <- c(wilcox.test(FetchData(object = test_AA, vars = "Alox15")$Alox15, FetchData(object = test_contr, vars = "Alox15")$Alox15, alternative = "two.sided")$p.value, 0, 0)
Pla2g2d <- c(wilcox.test(FetchData(object = test_AA, vars = "Pla2g2d")$Pla2g2d, FetchData(object = test_contr, vars = "Pla2g2d")$Pla2g2d, alternative = "two.sided")$p.value, 0, 0)

for(i in names) {
  Sub_cluster <- subset(Gene_de,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('AA'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))

  Ptger1 <- c(Ptger1, wilcox.test(FetchData(object = test_AA, vars = "Ptger1")$Ptger1, FetchData(object = test_contr, vars = "Ptger1")$Ptger1, alternative = "two.sided")$p.value)
  Ptger2 <- c(Ptger2, wilcox.test(FetchData(object = test_AA, vars = "Ptger2")$Ptger2, FetchData(object = test_contr, vars = "Ptger2")$Ptger2, alternative = "two.sided")$p.value)
  Ptger3 <- c(Ptger3, wilcox.test(FetchData(object = test_AA, vars = "Ptger3")$Ptger3, FetchData(object = test_contr, vars = "Ptger3")$Ptger3, alternative = "two.sided")$p.value)
  Ptger4 <- c(Ptger4, wilcox.test(FetchData(object = test_AA, vars = "Ptger4")$Ptger4, FetchData(object = test_contr, vars = "Ptger4")$Ptger4, alternative = "two.sided")$p.value)
  Ptges <- c(Ptges, wilcox.test(FetchData(object = test_AA, vars = "Ptges")$Ptges, FetchData(object = test_contr, vars = "Ptges")$Ptges, alternative = "two.sided")$p.value)
  Ptges2 <- c(Ptges2, wilcox.test(FetchData(object = test_AA, vars = "Ptges2")$Ptges2, FetchData(object = test_contr, vars = "Ptges2")$Ptges2, alternative = "two.sided")$p.value)
  Ptgs1 <- c(Ptgs1, wilcox.test(FetchData(object = test_AA, vars = "Ptgs1")$Ptgs1, FetchData(object = test_contr, vars = "Ptgs1")$Ptgs1, alternative = "two.sided")$p.value)
  Alox12b <- c(Alox12b, wilcox.test(FetchData(object = test_AA, vars = "Alox12b")$Alox12b, FetchData(object = test_contr, vars = "Alox12b")$Alox12b, alternative = "two.sided")$p.value)
  Ptgs2 <- c(Ptgs2, wilcox.test(FetchData(object = test_AA, vars = "Ptgs2")$Ptgs2, FetchData(object = test_contr, vars = "Ptgs2")$Ptgs2, alternative = "two.sided")$p.value)
  Ptgis <- c(Ptgis, wilcox.test(FetchData(object = test_AA, vars = "Ptgis")$Ptgis, FetchData(object = test_contr, vars = "Ptgis")$Ptgis, alternative = "two.sided")$p.value)
  Pla2g2e <- c(Pla2g2e, wilcox.test(FetchData(object = test_AA, vars = "Pla2g2e")$Pla2g2e, FetchData(object = test_contr, vars = "Pla2g2e")$Pla2g2e, alternative = "two.sided")$p.value)
  Ggt1 <- c(Ggt1, wilcox.test(FetchData(object = test_AA, vars = "Ggt1")$Ggt1, FetchData(object = test_contr, vars = "Ggt1")$Ggt1, alternative = "two.sided")$p.value)
  Gpx3 <- c(Gpx3, wilcox.test(FetchData(object = test_AA, vars = "Gpx3")$Gpx3, FetchData(object = test_contr, vars = "Gpx3")$Gpx3, alternative = "two.sided")$p.value)
  Pla2g10 <- c(Pla2g10, wilcox.test(FetchData(object = test_AA, vars = "Pla2g10")$Pla2g10, FetchData(object = test_contr, vars = "Pla2g10")$Pla2g10, alternative = "two.sided")$p.value)
  Pla2g12b <- c(Pla2g12b, wilcox.test(FetchData(object = test_AA, vars = "Pla2g12b")$Pla2g12b, FetchData(object = test_contr, vars = "Pla2g12b")$Pla2g12b, alternative = "two.sided")$p.value)
  Ephx2 <- c(Ephx2, wilcox.test(FetchData(object = test_AA, vars = "Ephx2")$Ephx2, FetchData(object = test_contr, vars = "Ephx2")$Ephx2, alternative = "two.sided")$p.value)
  Alox12 <- c(Alox12, wilcox.test(FetchData(object = test_AA, vars = "Alox12")$Alox12, FetchData(object = test_contr, vars = "Alox12")$Alox12, alternative = "two.sided")$p.value)
  Gpx1 <- c(Gpx1, wilcox.test(FetchData(object = test_AA, vars = "Gpx1")$Gpx1, FetchData(object = test_contr, vars = "Gpx1")$Gpx1, alternative = "two.sided")$p.value)
  Gpx4 <- c(Gpx4, wilcox.test(FetchData(object = test_AA, vars = "Gpx4")$Gpx4, FetchData(object = test_contr, vars = "Gpx4")$Gpx4, alternative = "two.sided")$p.value)
  Gpx2 <- c(Gpx2, wilcox.test(FetchData(object = test_AA, vars = "Gpx2")$Gpx2, FetchData(object = test_contr, vars = "Gpx2")$Gpx2, alternative = "two.sided")$p.value)
  Alox5 <- c(Alox5, wilcox.test(FetchData(object = test_AA, vars = "Alox5")$Alox5, FetchData(object = test_contr, vars = "Alox5")$Alox5, alternative = "two.sided")$p.value)
  Alox15 <- c(Alox15, wilcox.test(FetchData(object = test_AA, vars = "Alox15")$Alox15, FetchData(object = test_contr, vars = "Alox15")$Alox15, alternative = "two.sided")$p.value)
  Pla2g2d <- c(Pla2g2d, wilcox.test(FetchData(object = test_AA, vars = "Pla2g2d")$Pla2g2d, FetchData(object = test_contr, vars = "Pla2g2d")$Pla2g2d, alternative = "two.sided")$p.value)
}

Cell_Type <- c('Stem 1','Stem 2','Stem 3',"Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Gene_Level_PVal <- data.frame(Ptger1, Ptger2, Ptger3, Ptger4, Ptges, Ptges2, Ptgs1, Alox12b, Ptgs2, Ptgis, Pla2g2e, Ggt1, Gpx3, Pla2g10, Pla2g12b, Ephx2, Alox12, Gpx1, Gpx4, Gpx2, Alox5, Alox15, Pla2g2d)
rownames(Gene_Level_PVal) <- Cell_Type
for(i in colnames(Gene_Level_PVal)){
  Gene_Level_PVal[[i]][is.nan(Gene_Level_PVal[[i]])] <- NA
}
Gene_Level_PVal1 <- matrix(p.adjust(as.vector(as.matrix(Gene_Level_PVal)), method='BH'),ncol=23)
colnames(Gene_Level_PVal1) <- PGE2_genes
rownames(Gene_Level_PVal1) <- Cell_Type

write.csv(Gene_Level_PVal,'Organoid_PGE2_Gene_lvl_PVal.csv')

#PGE2
Pval_gene_list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'CD55')
Idents(PGE_de) = PGE_de$Cell_type
Sub_cluster <- subset(PGE_de,  idents = "Stem 1")
Idents(Sub_cluster) = Sub_cluster$type
test_AA <- subset(Sub_cluster,  idents =c('Pge2'))
test_contr <- subset(Sub_cluster,  idents =c('Control'))
DefaultAssay(test_contr) <- "MAGIC_RNA"
DefaultAssay(test_AA) <- "MAGIC_RNA"

Lgr5 <- c(wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value, 0, 0)
Ascl2 <- c(wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value, 0, 0)
S100a6 <- c(wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value, 0, 0)
Ly6a <- c(wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value, 0, 0)
CD55 <- c(wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value, 0, 0)

for(i in names) {
  Sub_cluster <- subset(PGE_de,  idents = i)
  Idents(Sub_cluster) = Sub_cluster$type
  test_AA <- subset(Sub_cluster,  idents =c('Pge2'))
  test_contr <- subset(Sub_cluster,  idents =c('Control'))
  DefaultAssay(test_contr) <- "MAGIC_RNA"
  DefaultAssay(test_AA) <- "MAGIC_RNA"
  
  Lgr5 <- c(Lgr5, wilcox.test(FetchData(object = test_AA, vars = "Lgr5")$Lgr5, FetchData(object = test_contr, vars = "Lgr5")$Lgr5, alternative = "two.sided")$p.value)
  Ascl2 <- c(Ascl2, wilcox.test(FetchData(object = test_AA, vars = "Ascl2")$Ascl2, FetchData(object = test_contr, vars = "Ascl2")$Ascl2, alternative = "two.sided")$p.value)
  S100a6 <- c(S100a6, wilcox.test(FetchData(object = test_AA, vars = "S100a6")$S100a6, FetchData(object = test_contr, vars = "S100a6")$S100a6, alternative = "two.sided")$p.value)
  Ly6a <- c(Ly6a, wilcox.test(FetchData(object = test_AA, vars = "Ly6a")$Ly6a, FetchData(object = test_contr, vars = "Ly6a")$Ly6a, alternative = "two.sided")$p.value)
  CD55 <- c(CD55, wilcox.test(FetchData(object = test_AA, vars = "Cd55")$Cd55, FetchData(object = test_contr, vars = "Cd55")$Cd55, alternative = "two.sided")$p.value)
}

Gene_Level_PVal <- Gene_Level_PVal[-c(4),]
Cell_Type <- c('Stem 1','Stem 2','Stem 3',"Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
Gene_Level_PVal <- data.frame(Lgr5, Ascl2, S100a6, Ly6a, CD55)
rownames(Gene_Level_PVal) <- Cell_Type
Gene_Level_PVal$Ly6a[is.nan(Gene_Level_PVal$Ly6a)] <- NA
Gene_Level_PVal$S100a6[is.nan(Gene_Level_PVal$S100a6)] <- NA

Gene_Level_PVal1 <- matrix(p.adjust(as.vector(as.matrix(Gene_Level_PVal)), method='BH'),ncol=5)
colnames(Gene_Level_PVal1) <-  Pval_gene_list
rownames(Gene_Level_PVal1) <- Cell_Type

write.csv(Gene_Level_PVal1,'Pge2_Gene_lvl_PVal.csv')

