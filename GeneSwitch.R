#GeneSwitch Script

library(GeneSwitches)
library(SingleCellExperiment)

#organoid
# create SCE object for GeneSwtich analysis
Idents(Gene_de) <-  Gene_de$type
#control
AA_Control <- subset(Gene_de,  idents =c('Control'))
all.genes <- rownames(AA_Control)
AA_Control <- ScaleData(AA_Control, features = all.genes, assay =  'RNA')
logexpdata <- AA_Control@assays$RNA@scale.data
rd <- Embeddings(object = AA_Control, reduction = "pca")

sce_p1 <- SingleCellExperiment(assays = List(expdata = logexpdata))
colData(sce_p1)$Pseudotime <- AA_Control$Pseudotime
colData(sce_p1)$type <- AA_Control$type
reducedDims(sce_p1) <- SimpleList(PCA = rd)
#sce_p1 <- binarize_exp(sce_p1)

h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
abline(v=0.1, col="blue")}

bn_cutoff <- 0.1
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)
sce_p1 <- find_switch_logistic_fastglm(sce_p1, show_warning = FALSE)

#AA
AA_AA <- subset(Gene_de,  idents =c('AA'))
all.genes <- rownames(AA_AA)
AA_AA <- ScaleData(AA_AA, features = all.genes, assay = 'RNA')
logexpdata <- AA_AA@assays$RNA@scale.data
rd <- Embeddings(object = AA_AA, reduction = "pca")

sce_p2 <- SingleCellExperiment(assays = List(expdata = logexpdata))
colData(sce_p2)$Pseudotime <- AA_AA$Pseudotime
colData(sce_p2)$type <- AA_AA$type
reducedDims(sce_p2) <- SimpleList(PCA = rd)

bn_cutoff <- 0.1
sce_p2 <- binarize_exp(sce_p2, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)

sce_p2 <- find_switch_logistic_fastglm(sce_p2, show_warning = FALSE)

#filter based on same cutoff
sg_p1 <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.03)
sg_p2 <- filter_switchgenes(sce_p2, allgenes = TRUE, r2cutoff = 0.03)

#plotting
sg_com <- common_genes(sg_p2, sg_p1, r2cutoff = 0.1,
                       path1name = "AA", path2name = "Control")
common_genes_plot(sg_com, timedata = sce_p1$Pseudotime)

sg_disgs <- distinct_genes(sg_p2, sg_p1, r2cutoff = 0.0,
                           path1name = "AA", path2name = "Control",
                           path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
plot_timeline_ggplot(sg_disgs, timedata = sce_p1$Pseudotime, color_by = "Paths", iffulltml = FALSE, txtsize = 2) +
  scale_color_manual(values = c("#d95f02", "#1b9e77")) +
  theme(plot.title = element_blank(), 
        text = element_text(size=6), legend.position = 'none', 
        axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Gene_Switch.pdf', width=3, height=3, units="in")

# Arasco

Idents(arasco.obj) <-  arasco.obj$type
#control
AA_Control <- subset(arasco.obj,  idents =c('Control'))
all.genes <- rownames(AA_Control)
AA_Control <- ScaleData(AA_Control, features = all.genes, assay =  'RNA')
logexpdata <- AA_Control@assays$RNA@scale.data
rd <- Embeddings(object = AA_Control, reduction = "pca")

sce_p1 <- SingleCellExperiment(assays = List(expdata = logexpdata))
colData(sce_p1)$Pseudotime <- AA_Control$Pseudotime
colData(sce_p1)$type <- AA_Control$type
reducedDims(sce_p1) <- SimpleList(PCA = rd)
#sce_p1 <- binarize_exp(sce_p1)

h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.1, col="blue")}

bn_cutoff <- 0.1
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)
sce_p1 <- find_switch_logistic_fastglm(sce_p1, show_warning = FALSE)

#AA
AA_AA <- subset(arasco.obj,  idents =c('Arasco'))
all.genes <- rownames(AA_AA)
AA_AA <- ScaleData(AA_AA, features = all.genes, assay = 'RNA')
logexpdata <- AA_AA@assays$RNA@scale.data
rd <- Embeddings(object = AA_AA, reduction = "pca")

sce_p2 <- SingleCellExperiment(assays = List(expdata = logexpdata))
colData(sce_p2)$Pseudotime <- AA_AA$Pseudotime
colData(sce_p2)$type <- AA_AA$type
reducedDims(sce_p2) <- SimpleList(PCA = rd)

bn_cutoff <- 0.1
sce_p2 <- binarize_exp(sce_p2, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)

sce_p2 <- find_switch_logistic_fastglm(sce_p2, show_warning = FALSE)

#filter based on same cutoff
sg_p1 <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.03)
sg_p2 <- filter_switchgenes(sce_p2, allgenes = TRUE, r2cutoff = 0.03)

#plotting
sg_com <- common_genes(sg_p2, sg_p1, r2cutoff = 0.2,
                       path1name = "Arasco", path2name = "Control")
common_genes_plot(sg_com, timedata = sce_p1$Pseudotime)

sg_disgs <- distinct_genes(sg_p2, sg_p1, r2cutoff = .1,
                           path1name = "Arasco", path2name = "Control",
                           path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
plot_timeline_ggplot(sg_disgs, timedata = sce_p1$Pseudotime, color_by = "Paths", iffulltml = FALSE, txtsize = 2) +
  scale_color_manual(values = c("#d95f02", "#1b9e77")) +
  theme(plot.title = element_blank(), 
        text = element_text(size=6), legend.position = 'none', 
        axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Arasco_Gene_Switch.pdf', width=3, height=3, units="in")

#specific genes
#gs_genelists <-  c('S100a6','Cd55','Ly6a','Lgr5','Areg','Ereg','Ascl2','Sapcd2','Tff3','Ccnf','Reep4','Spdl1','Anln','Aspm','Pbk','Ptger4','Stmn1')
#sg_gtypes1 <- filter_switchgenes(sce_p1, allgenes = FALSE, genelists = gs_genelists)
#sg_gtypes2 <- filter_switchgenes(sce_p1, allgenes = FALSE, genelists = gs_genelists)

## filter top 15 best fitting switching genes among all the genes
#sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, topnum = 15)

## combine switching genes and remove duplicated genes from sg_allgenes
#sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])
p#lot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)
