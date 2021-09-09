enrich1 <-  ptger_human
enrich1 <-  ptger_mouse

Enrich_idents <- c('gene_name','log2FoldChange.aa_ctrl_d1','log2FoldChange.pge2_ctrl_d1','log2FoldChange.aa_ctrl_d3', 'log2FoldChange.pge2_ctrl_d3')
enrich1 <-  data.frame(enrich1)
enrich1 <- enrich1[Enrich_idents]

enrich1$V4
enrich_filter <- enrich1[enrich1$V4 %in% Enrich_idents,]
colnames(enrich_filter) <- enrich1[1,]

coef_names <- c('AAvsCtrl_D1','AAvsCtrl_D3','AAvsCtrl_D6')
enrich_filter <- enrich_filter[enrich_filter$Coefficient %in% coef_names,]
enrich_filter <-enrich1
Enrich_idents <- c('FEVR_CTNNB1_TARGETS_UP','GO_WOUND_HEALING','ZWANG_TRANSIENTLY_UP_BY_2ND_EGF_PULSE_ONLY', 'GO_LIPID_METABOLIC_PROCESS','BOQUEST_STEM_CELL_UP', 'GO_REGULATION_OF_CELL_PROLIFERATION', 'GO_REGULATION_OF_CALCIUM_ION_TRANSPORT', 'NABA_ECM_REGULATORS' )

enrich_filter$NES <- as.numeric(enrich_filter$NES)
enrich_filter$padj <- as.numeric(enrich_filter$padj)
enrich_filter$pathway <- factor(enrich_filter$pathway,levels = Enrich_idents)
pathwayColors =rev(viridis::plasma(10))
pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")

ggplot(data = enrich_filter,aes(x=pathway,y=c('log2FoldChange.aa_ctrl_d1','log2FoldChange.pge2_ctrl_d1'))) + theme_cem +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  facet_grid(cols = vars(Coefficient), scales="free", switch = 'y') +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = padj))+
  geom_point(color = 'black' ,size = 1.2 ) +
  geom_point(aes(color = padj) ,size = 1 ) +
  scale_color_gradientn(colors=pathwayColors, breaks = c(.05, .01, .001, .00001),labels =  c('.05', '.01', '.001', '<.0001'),limits = c(0,.05), trans = scales::boxcox_trans(.25)) +
  coord_flip() +
  scale_y_continuous(limits = c(0,2.5),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 

All_genes <-  c('CD55','SIK1','NR4A1','L1CAM','ANKRD1','TCF4','MYOF','AREG','CCND1','ANXA3','MSLN','LYPD6','OLFM4','CDKN2B','ACE','DCLK2','BAD','RXRA','TP53I13')
Heatmap_data <-  Heatmap_aa
Heatmap_data <- Heatmap_data[Heatmap_data$gene_name %in% All_genes,]
rownames(Heatmap_data) <-  Heatmap_data$gene_name
Heatmap_data <- Heatmap_data[,c('log2FoldChange.aa_ctrl', 'gene_name')]
Heatmap_data <- enrich1
library(reshape)
Heatmap_data <-melt(Heatmap_data)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
All_genes <-  c('CD55','SIK1','NR4A1','L1CAM','ANKRD1','TCF4','MYOF','AREG','CCND1','ANXA3','MSLN','LYPD6','OLFM4','CDKN2B','ACE','DCLK2','BAD','RXRA','TP53I13')
Heatmap_data$gene_name <- factor(Heatmap_data$gene_name,levels = All_genes)
ggplot(Heatmap_data, aes(gene_name, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-2,2)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="Ptger_Mouse_logfc_scale.pdf",  width=4, height=2, units="in")

title_size <- 10
text_size <- 8
boxplot_height_in <- 2
boxplot_width_in <- 3
point_size <- 0.8
aa_tissue$measure

aa_tissue1 <- aa_tissue
idents1 <- c('diet','metabolite','measure')
aa_tissue <- aa_tissue[idents1]


aa_tissue$measure <- log(aa_tissue$measure, 10)
aa_tissue[is.na(aa_tissue)] <- 0

aa_tissue_mean <- aa_tissue %>%
  group_by(diet, metabolite) %>%
  summarise_at(vars(measure), list(name = mean))
aa_tissue_mean

aa_tissue_split <- split(aa_tissue, aa_tissue$diet)
intestine_FC <- foldchange(aa_tissue_mean_control$arasco$name, aa_tissue_mean_control$control$name)
intestine_FC <- data.frame(intestine_FC)
intestine_FC$metabolite <- aa_tissue_mean_control$arasco$metabolite
intestine_FC1 <- intestine_FC
intestine_FC1$intestine_FC[!is.finite(intestine_FC1$intestine_FC)] <- 0
intestine_FC
aa_tissue_mean
library(reshape)
Heatmap_data <-melt(intestine_FC1)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
aa_tissue_mean
aa_tissue
length(aa_tissue_split$arasco[which(aa_tissue_split$arasco$metabolite == i),]$measure)

zscores <- zscoreT(aa_tissue$measure, df = length(aa_tissue$measure))
aa_zscore <- aa_tissue
aa_zscore$measure <- zscores

aa_zscore_mean <- aa_zscore %>%
  group_by(diet, metabolite) %>%
  summarise_at(vars(measure), list(name = mean))
aa_zscore_mean$diet

my_levels <- c('control','arasco')
aa_zscore_mean$diet <- factor(x = aa_zscore_mean$diet, levels = my_levels)

ggplot(Heatmap_data, aes(metabolite, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,3)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="Ptger_Mouse_logfc_scale.pdf",  width=4, height=2, units="in")

ggplot(aa_zscore_mean, aes(metabolite, diet, fill= name)) + geom_tile() + 
  scale_fill_gradientn(name = "Log Counts", colors = pathwayColorsDiff, limits = c(0,7.5)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="Ptger_Mouse_logfc_scale.pdf",  width=4, height=2, units="in")

#identify Bcat_intersect
signatures_fc_beta_Catenin_induced <- read.csv('~/analysis/signatures_fc_beta-Catenin-induced.csv')
signatures_fc_fetal <- read.csv('~/analysis/signatures_fc_fetal_spheroid.csv')
signatures_fc_homeostatic <- read.csv('~/analysis/signatures_fc_homeostatic-signature.csv')
signatures_fc_regeneration <- read.csv('~/analysis/signatures_fc_regeneration-induced.csv')
signatures_fc_granuloma <- read.csv('~/analysis/signatures_fc_granuloma-induced.csv')
signatures_fc_radiation <-read.csv('~/analysis/signatures_fc_radiation-induced.csv')

Fetal <- intersect(signatures_fc_beta_Catenin_induced$gene_name,signatures_fc_fetal$gene_name)
Homesostatic <- intersect(signatures_fc_beta_Catenin_induced$gene_name,signatures_fc_homeostatic$gene_name)
Regeneration <- intersect(signatures_fc_beta_Catenin_induced$gene_name,signatures_fc_regeneration$gene_name)
Granuloma <- intersect(signatures_fc_beta_Catenin_induced$gene_name,signatures_fc_granuloma$gene_name)
Radiation <- intersect(signatures_fc_beta_Catenin_induced$gene_name,signatures_fc_radiation$gene_name)
All <- list(Fetal,Homesostatic,Regeneration,Granuloma,Radiation)

write.csv(Fetal, 'Intersect_list.csv')
write.csv(Homesostatic, 'Homesostatic.csv')
write.csv(Regeneration, 'Regeneration.csv')
write.csv(Granuloma, 'Granuloma.csv')
write.csv(Radiation, 'Radiation.csv')



fc_aa_ctrl_mouse
Wnt_targets$gene_name

fc_aa_ctrl_mouse1 <- fc_aa_ctrl_mouse[fc_aa_ctrl_mouse$gene_name %in% Wnt_targets$gene_name,]
fc_aa_ctrl_mouse1[is.na(fc_aa_ctrl_mouse1)] = 0

library(reshape2)
fc_aa_ctrl_mouse1 <-melt(fc_aa_ctrl_mouse1)
pathwayColorsDiff = rev(brewer.pal(50, "RdBu"))
All_genes <-  Wnt_targets$gene_name
fc_aa_ctrl_mouse1$gene_name <- factor(fc_aa_ctrl_mouse1$gene_name,levels = All_genes)
ggplot(fc_aa_ctrl_mouse1, aes(gene_name, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, breaks = c(-4, -2, 0, 2, 4), limits = c(-4,4)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank())
ggsave(file="Bulk_logfc_bcat.pdf",  width=15, height=2, units="in")



gene.level.list <- c('Lgr5', 'Ascl2', 'S100a6', 'Ly6a', 'Cd55', 'Areg', 'Ereg','Clu')
for(i in gene.level.list) {
  All <- VlnPlot(arasco.obj, group.by = 'Cell_type' ,  split.by = "type", features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + 
    theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('arasco_Gene', i, '.pdf'), plot=All, width=2, height=2.5, units="in")
}


for(i in gene.level.list) {
  All <- VlnPlot(Gene_de, group.by = 'Cell_type' ,  split.by = "type", features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + 
    theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('AA_Gene', i, '.pdf'), plot=All, width=2, height=2.5, units="in")
}


for(i in gene.level.list) {
  All <- VlnPlot(PGE_de, group.by = 'Cell_type' ,  split.by = "type", features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + 
    theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('pge_Gene', i, '.pdf'), plot=All, width=2, height=2.5, units="in")
}

Idents(organoid) <- organoid$type
organ_Control <- subset(organoid,  idents =c('Control'))

Idents(arasco.obj) <- arasco.obj$type
ARA_Control <- subset(arasco.obj,  idents =c('Control'))
DotPlot_Sig <- c("Ptger1", "Ptger2", "Ptger3", "Ptger4") 


DotPlot(ARA_Control, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text= element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file= 'Arasco_Haber_ptger_dotplot.pdf', width=2.5, height=2.6, units="in")

DotPlot(organ_Control, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text= element_blank(), legend.title = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file= 'organoid_ptger_check_dotplot.pdf', width=2.5, height=2.6, units="in")






