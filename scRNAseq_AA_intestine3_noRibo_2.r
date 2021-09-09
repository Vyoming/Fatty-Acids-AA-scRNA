library(hdf5r)
library(Seurat)
library(Seurat)
library(scales)
library(dplyr)
library(future)
library(tidyverse)
library(magrittr)
library(fgsea)
library(dplyr)
library(scales)
library(viridis)
library(biomaRt)
library(mgcv)

setwd("~/Projects/BeyazS/scRNAseq_AA/intestine3_noRibo/")


getGeneOrth = function()
{
	human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
	mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
	#genesOrth = getLDS(attributes = c("mgi_symbol"), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)
	genesOrth = getLDS(attributes = c("hgnc_symbol"), mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
	genesOrth2 = genesOrth[ genesOrth$MGI.symbol != "", ]
	return(genesOrth2)
}


getPathways = function(genesOrth = NULL)
{
	if(is.null(genesOrth)) genesOrth = getGeneOrth()
	pathwaysM = c(gmtPathways("~/GSEA/AA_sigs_mouse.gmt"))
	#pathwaysM = c(pathwaysM, gmtPathways("~/GSEA/NormalGCB_murine.gmt"))
	
	pathwaysH = c()
	pathwaysH = c(pathwaysH, gmtPathways("~/GSEA/AA_sigs_human.gmt"))
	#pathwaysH = c(pathwaysH, gmtPathways("~/GSEA/xCell_sigs.gmt"))
	#pathwaysH = c(pathwaysH, gmtPathways("~/GSEA/NormalBCell_filtered.gmt") )
	
	pathwaysH2M = lapply(pathwaysH, function(gs, map) {
	  gs = with(map, MGI.symbol[HGNC.symbol %in% gs])
	  gs
	}, genesOrth)

	pathways = c(pathwaysH2M[ !duplicated(names(pathwaysH2M))], pathwaysM[ !duplicated(names(pathwaysM))])
	return(pathways)
}


calculatePathwayScores = function(seuratObj, dimred, pathwayList, method="seurat", minCells=50,  minGenes=10)
{
	curData = seuratObj@assays$integrated@scale.data
	curData = curData[ rowSums(curData > 0) >= minCells, ]
	curData = curData[ rowSds(curData) > 0, ]
	zscore = t(curData) # scale(t(curData))
	pathwayNames = unique( gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathwayList)) )

	pathwayScores = dimred
	for(curPathway in pathwayNames)
	{
		curPathwayName = gsub("\\*", "", curPathway)
		curUpSum = rep(0, nrow(zscore))
		curDnSum = rep(0, nrow(zscore))
		anyHits = F
		curUpGenesList = c()
		curDnGenesList = c()
		
		if(grepl("\\*", curPathway) & gsub("\\*", "_UP", curPathway) %in% names(pathwayList))
		{
			curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_UP", curPathway) ]] )
		}
		if(grepl("\\*", curPathway) & gsub("\\*", "_upreg", curPathway) %in% names(pathwayList))
		{
			curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_upreg", curPathway) ]] )
		}
		if(grepl("\\*", curPathway) & gsub("\\*", "_up_", curPathway) %in% names(pathwayList))
		{
			curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_up_", curPathway) ]] )
		}

		if(grepl("\\*", curPathway) & gsub("\\*", "_DN", curPathway) %in% names(pathwayList))
		{
			curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_DN", curPathway) ]] )
		}
		if(grepl("\\*", curPathway) & gsub("\\*", "_downreg", curPathway) %in% names(pathwayList))
		{
			curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_downreg", curPathway) ]] )
		}
		if(grepl("\\*", curPathway) & gsub("\\*", "_down_", curPathway) %in% names(pathwayList))
		{
			curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_down_", curPathway) ]] )
		}
		if( length(curUpGenesList) == 0 & length(curDnGenesList) == 0 & curPathway %in% names(pathwayList))
		{
			curUpGenesList = c(curUpGenesList, pathwayList[[ curPathway ]] )
		}
		curUpGenesList = intersect(curUpGenesList, colnames(zscore) )
		curDnGenesList = intersect(curDnGenesList, colnames(zscore) )
		
		if( length(curUpGenesList) + length(curDnGenesList) <= minGenes ) # min 10 genes
		{
			next
		}
		
		if(method == "seurat")
		{
			curUpGenesList = list( curUpGenesList )
			curDnGenesList = list( curDnGenesList )
			seuratObjUp = NULL
			seuratObjDn = NULL
			tryCatch({seuratObjUp = AddModuleScore(object = seuratObj, features = curUpGenesList, ctrl = 5, name = 'Score'); curUpSum = seuratObjUp$Score1; }, error = function(e) { print(e) })
			tryCatch({seuratObjDn = AddModuleScore(object = seuratObj, features = curDnGenesList, ctrl = 5, name = 'Score'); curDnSum = seuratObjDn$Score1; }, error = function(e) { print(e) })
			curScore = curUpSum - curDnSum 
			if(is.null(names(curScore))) next
		
		}else
		{
			if( length(curUpGenesList) > 0)
			{
				curUp = zscore[, colnames(zscore) %in% curUpGenesList, drop=F]
				curUpSum = rowSums(curUp)
			}
			if( length(curDnGenesList) > 0)
			{
				curDn = zscore[, colnames(zscore) %in% curDnGenesList, drop=F]
				curDnSum = rowSums(curDn)
			}
			
			curScore = curUpSum - curDnSum
			names(curScore) = rownames(zscore)
		}
		
		curScore = data.frame(Cell=names(curScore), Score=scale(curScore))
		colnames(curScore)[2] = curPathwayName
		pathwayScores = merge(pathwayScores, curScore, by="Cell")	
	}

	pathwayScores2 = melt(pathwayScores, colnames(dimred), variable.name="Pathway", value.name="Score")
	pathwayScores2$Pathway = as.character(pathwayScores2$Pathway)
	return(pathwayScores2)
}


smoothScorePredict = function(trainDat, newGrid, knn=30)
{
	scorePredicted = data.frame(score = as.numeric(NA), 
								x = newGrid[, 1], 
								y = newGrid[, 2])
	# run KKNN
	scoreKknn = kknn::kknn(score ~ ., 
						train = trainDat, 
						test = scorePredicted, 
						kernel = "gaussian", 
						k = knn)

	scorePredicted %<>% mutate(score = fitted(scoreKknn))
	return(scorePredicted)
}

smoothScore2d = function(score, x, y, type=NULL, numGrid = 100, knn = 100, m = 2, expand=0.05, xrng=NULL, yrng=NULL)
{
	library(ash)
	curDat = data.frame(score=score, x=x, y=y)
	if(is.null(xrng)) xrng = range(curDat$x)
	if(is.null(yrng)) yrng = range(curDat$y)
	xdiff = xrng[2] - xrng[1]
	ydiff = yrng[2] - yrng[1]
	xrng[1] = xrng[1] - xdiff*expand
	xrng[2] = xrng[2] + xdiff*expand
	yrng[1] = yrng[1] - ydiff*expand
	yrng[2] = yrng[2] + ydiff*expand
	bins = bin2(cbind(curDat$x, curDat$y), ab = rbind(xrng, yrng), nbin = c(numGrid, numGrid))
	binCounts = ash2(bins, m = c(m, m))
	gridDat = data.frame( expand.grid( x = binCounts$x, y = binCounts$y), density = melt(binCounts$z)[,3] )
	gridDat2 = gridDat[ gridDat$density > 0, ]

	if(is.null(type))
	{
		return( smoothScorePredict(curDat, gridDat2, knn = knn) )
	}else
	{
		curDat$type = type
		allPredicted = ddply(curDat, "type", function(x) { smoothScorePredict( x[, colnames(x) != "type"], gridDat2, knn = knn) } )
		return(allPredicted)
	}
}


savePlot = function(outname, width = 11, height = 7, dpi = 150, pdfOutput = T, pngOutput = T)
{
	outname = gsub(".png$", "", outname)
	outname = gsub(".pdf$", "", outname)
	
	tryCatch(
	{
		if(pngOutput)
		{
			dev.copy(png, filename = paste0(outname, ".png"), width = width, height = height, units = 'in', res = dpi)
			dev.off ()
		}
		if(pdfOutput)
		{
			dev.copy2pdf(file = paste0(outname, ".pdf"), width = width, height = height)
		}
	}, error = function(e)
	{
		print("Error during output")
		stop(e)
	})
}



plotSymbols = function(dat, xName, yName, className, dotSize=3, alpha=0.5, do.legend = TRUE, do.letters = FALSE, label.size=4, colors=NULL, font.size=NULL, title = "", minCount = 25)
{
	SYMBOLS = c(LETTERS,letters,c(0:9),
              c("!","@","#","$","%","^","&","*","(",")",")","-",
                "+","_","=",";","/","|","{","}","~"))
	
	classFreq = as.data.frame(table(dat[, className]))
	classFreq2 = classFreq[ classFreq$Freq >= minCount, ]
	dat[ dat[,className] %in% classFreq[,1][classFreq[,2] < minCount] , className] = "(Other)"
	
	numClass = length(unique(dat[, className]))
	df = data.frame(x = dat[, xName], y=dat[, yName], ident=dat[, className])
	df$ident = factor(df$ident)
	p = ggplot(df, aes(x = x, y = y))
	p = p + geom_point(aes(color=ident), size=dotSize, alpha=alpha, stroke=0)
  
	if( do.letters == TRUE) 
	{
		symbols = SYMBOLS[1:numClass]
		names(symbols) = levels(df$ident)
		p = p + geom_point(aes(shape=ident), size=2*dotSize/5, color='black')
		p = p + scale_shape_manual(values=symbols)
	}
    #p = p + geom_point(aes(shape=initIdent), size=dot.size/2)
    #p = p + scale_shape_identity()
	
	if (is.null(font.size)) 
	{
      font.size = 250 * (1/numClass)
      font.size = max(font.size,5)
      font.size = min(font.size,10)
    }
	
    if (numClass > 5 & numClass < 60) 
	{
      p = p + theme(legend.position="bottom",legend.direction="vertical",
                    legend.text=element_text(size=6), legend.title = element_blank()) + guides(shape=guide_legend(ncol=5, override.aes = list(size=2.5, alpha=1)))
    }else if (numClass > 60)
	{
      p = p + theme(legend.position="bottom",legend.direction="vertical",
                    legend.text=element_text(size=6), legend.title = element_blank()) + guides(shape=guide_legend(ncol=9, override.aes = list(size=2, alpha=1)))
    }else 
	{
      p = p + theme(legend.text=element_text(size=font.size),legend.title = 
                    element_blank()) + 
					guides(color=guide_legend(ncol=1, override.aes = list(size=3, alpha=1)))
    }
  
	if(is.null(colors)) colors = iwanthue(numClass)
	lev = levels(df$ident)
	cols = colors[1:length(lev)]
	names(cols) = lev
	cols[names(cols) == 'X'] = 'black'

	p = p + scale_color_manual(values = cols)
	p = p + xlab(xName) + ylab(yName) + ggtitle(title)

	if (do.legend == FALSE) 
	{
	p = p + theme(legend.position="none")
	}

	p = p + theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"))

	return(p)
}


# parameters
minGenes = 500
minCells = 100

geneColors = brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
geneColorsDiff = brewer.pal(9, "PiYG")
dimredPlotWidth = 13
dimredPlotHeight = 10
useGraphUmap = F

colorsType = c(
	Control = "#2bbcff",
    AA = "#ffa72b",
    PGE2 = "#b82bff"
)

swatch(colorsType)


# Input files
curFilename = "AA_intestine2_SB4and5"
curTitle = "AA treatment, intestinal organoids"

inputFiles = list()
inputFiles[["SB04_C"]] = list()
inputFiles[["SB04_A"]] = list()
inputFiles[["SB05_C"]] = list()
inputFiles[["SB05_A"]] = list()
inputFiles[["SB05_P"]] = list()

inputFiles[["SB04_C"]][["Control"]] = "/home/cem2009/Projects/BeyazS/scRNAseq_AA/data/Intestinal/SB04/SB04_Ctrl_filtered_feature_bc_matrix.h5"
inputFiles[["SB04_A"]][["AA"]] = "/home/cem2009/Projects/BeyazS/scRNAseq_AA/data/Intestinal/SB04/SB04_AA_filtered_feature_bc_matrix.h5"
inputFiles[["SB05_C"]][["Control"]] = "/home/cem2009/Projects/BeyazS/scRNAseq_AA/data2/count/Beyaz_05_10x_Control/outs/filtered_feature_bc_matrix.h5"
inputFiles[["SB05_A"]][["AA"]] = "/home/cem2009/Projects/BeyazS/scRNAseq_AA/data2/count/Beyaz_05_10x_AA/outs/filtered_feature_bc_matrix.h5"
inputFiles[["SB05_P"]][["PGE2"]] = "/home/cem2009/Projects/BeyazS/scRNAseq_AA/data2/count/Beyaz_05_10x_PGE/outs/filtered_feature_bc_matrix.h5"



# Prepare Seurat object
{
	# Read and process different batches separately
	seuratObjList = list()
	allGenes = c()
	for(curDataset in names(inputFiles))
	{
		ctDat = NULL
		# Since files are in h5 read each sample separately and merge
		for(curFile in names(inputFiles[[curDataset]]))
		{
			curDat = Read10X_h5(unlist(inputFiles[[curDataset]][[curFile]]), use.names = TRUE, unique.features = TRUE)
			curDat = curDat[ !grepl("^(Rp)", rownames(curDat), ignore.case=T), ] # curDat = curDat[ !grepl("^(Rp|mt-)", rownames(curDat), ignore.case=T), ]
			colnames(curDat) = paste0(curDataset, "_", curFile, "_", colnames(curDat))
			if(is.null(ctDat)) ctDat = curDat else ctDat = cbind(ctDat, curDat)
		}
		
		# Filter and normalize 
		curSeurat = CreateSeuratObject(counts = ctDat, min.cells = minCells, min.features = minGenes)
		mitoGenes = grep(pattern = "^mt-", x = rownames(x = curSeurat@assays$RNA), value = TRUE, ignore.case = T)
		curSeurat = PercentageFeatureSet(curSeurat, features = mitoGenes, col.name = "percent.mt")
		curSeurat = SCTransform(curSeurat, vars.to.regress = "percent.mt", verbose = FALSE) 
		
		countSeurat = curSeurat@assays$RNA@counts
		
		curSeurat[["type"]] = splitGet(colnames(countSeurat), "_", 1)
		curSeurat[["batch"]] = splitGet(curDataset, "_", 1)
			
		seuratObjList[[curDataset]] = curSeurat
		allGenes = unique(c(allGenes, rownames(curSeurat@assays$RNA@counts)))	
	}

	## Integrate batches together using anchors
	options(future.globals.maxSize=Inf)
	integrationFeatures = SelectIntegrationFeatures(object.list = seuratObjList, nfeatures = 3000)
	seuratObjList = PrepSCTIntegration(object.list = seuratObjList, anchor.features = integrationFeatures, verbose = FALSE)
	seuratAnchors = FindIntegrationAnchors(object.list = seuratObjList, normalization.method = "SCT", anchor.features = integrationFeatures, verbose = FALSE)
	seuratAll = IntegrateData(anchorset = seuratAnchors, normalization.method = "SCT", verbose = FALSE)
	#DefaultAssay(seuratAll) = "integrated"
		
	countSeurat = seuratAll@assays$RNA@counts
	geneName = rownames(countSeurat)
	cellDat = data.frame(CellID=colnames(countSeurat))
	cellDat$Type = splitGet(cellDat$CellID, "_", 3)
	cellDat$Batch = splitGet(cellDat$CellID, "_", 1)
	cellDat$Sample = paste0(cellDat$Type, "_", cellDat$Batch)
	
	seuratAll[["Type"]] = cellDat$Type
	seuratAll[["Sample"]] = cellDat$Sample
	seuratAll[["Batch"]] = cellDat$Batch
	
	seuratAll$Type = factor(seuratAll$Type, levels = c("Control", "AA", "PGE2"))

	
	
	seuratAllTmp = seuratAll
	
	VlnPlot(seuratAll, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Batch") + expand_limits(y=c(0))
	ggsave("seurat_preFiltering_dist.png", width=12, height=18)
	ggsave("seurat_preFiltering_dist.pdf", width=12, height=18)
	
	seuratAll = subset(seuratAll, subset = nFeature_RNA > 500 & percent.mt < 25)
	

	#seuratAll = ScaleData(seuratAll, verbose = T)
	seuratAll = RunPCA(seuratAll, npcs = 50, verbose = FALSE)
	seuratAll = FindNeighbors(seuratAll, dims = 1:30)
	seuratAll = FindClusters(seuratAll, resolution = 0.45)
	
	if( reticulate::py_module_available("umap") & useGraphUmap == TRUE )
	{
		seuratAll = RunUMAP(seuratAll, reduction = "pca", graph = "integrated_snn", assay = NULL, 
		  n.components = 2L, metric = "euclidean", n.epochs = 500,
		  learning.rate = 1, min.dist = 0.3, spread = 2,
		  repulsion.strength = 1, negative.sample.rate = 5, seed.use = 42)
		  # n.neighbors = 15L, set.op.mix.ratio = 1, local.connectivity = 1L
	}else
	{
		seuratAll = RunUMAP(seuratAll, reduction = "pca", dims = 1:30)
	}



	DimPlot(seuratAll, reduction = "umap", label = TRUE, repel = TRUE)
	DimPlot(seuratAll, reduction = "umap", group.by = "Type")
	DimPlot(seuratAll, reduction = "umap", group.by = "Batch")
	DimPlot(seuratAll, reduction = "umap", group.by = "Sample")
	DimPlot(seuratAll, reduction = "umap", group.by = "Sample") + facet_wrap(~Sample)
	
	# if any clusters are to be removed remove them here and recluster
	# ....
	# TODO: Ask Semir to see if clusters are to be removed due to leaky sorting or unwanted cell types
	# ....




	# reorder clusters to make close clusters get similar numbers
	{
		seuratTemp = seuratAll
		dimred = data.frame(seuratAll@reductions$umap@cell.embeddings)
		dimred$Cell = rownames(dimred)
		dimred$Cluster = seuratAll$seurat_clusters
		clusterMedian = dimred %>%
							group_by(Cluster) %>%
							summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
		
		clusterMedianScaled = setRemoveRownames(as.data.frame(clusterMedian), 1)
		clusterMedianScaled = scale(clusterMedianScaled)
		clusterOrder = getTSPOrder(clusterMedianScaled)
		tempClustering = dimred$Cluster
		seuratAll$seurat_clusters = as.character(seuratAll$seurat_clusters)
		i = 1
		for(curClust in clusterOrder)
		{
			seuratAll$seurat_clusters[ seuratAll$seurat_clusters == curClust ] = paste0("C", str_pad(i, 2, pad = "0"))
			i = i + 1
		}
		clusterOrder2 = paste0("C", str_pad(as.character(as.numeric(clusterOrder) + 1), 2, pad = "0"))
		#seuratAll$seurat_clusters = gsub("C", "", seuratAll$seurat_clusters)
		seuratAll$seurat_clusters = factor(seuratAll$seurat_clusters, levels=clusterOrder2)
		names(seuratAll$seurat_clusters) = names(tempClustering)
	}

	clusterIdents = as.data.frame(table(seuratAll$seurat_clusters, Idents(seuratAll)))
	clusterIdents = clusterIdents[clusterIdents[, 3] > 0, ]
	clusterIdents2 = as.character(clusterIdents[, 1])
	names(clusterIdents2) = clusterIdents[, 2]
	seuratAll = RenameIdents(seuratAll, clusterIdents2)
	seuratAll@active.ident = factor(seuratAll@active.ident, levels = sort(levels(seuratAll@active.ident)))
	seuratAll$seurat_clusters = seuratAll@active.ident


	seuratObj = seuratAll
	dimred = data.frame(seuratObj@reductions$umap@cell.embeddings)
	dimred$Cell = rownames(dimred)
	all.equal(names(seuratObj$seurat_clusters), dimred$Cell)
	all.equal(names(seuratObj$Type), dimred$Cell)

	dimred$Cluster = factor(seuratObj$seurat_clusters, levels=sort(unique(as.character(seuratObj$seurat_clusters))))
	dimred$Sample = seuratObj$Sample
	dimred$Type = seuratObj$Type
	dimred$Batch = seuratObj$Batch


	#save.image("seurat_SCT_SB04SB05_umap.RData")
}







# Calculate cell type from SingleR
{
	if(FALSE)
	{
		# Requires R 3.6+
		library(SingleR)
		library(scRNAseq)
		library(scater)
		library(scran)
		
		

		immgen = ImmGenData()
		mousernaseq = MouseRNAseqData()
		
		haber = read.tsv2("/home/cem2009/Projects/BeyazS/scRNAseq_AA/intestine/haber/AtlasFullLength_log2TPM.txt")
		haber = setRemoveRownames(haber, 1)
		haberG = SingleCellExperiment(assays = list(logcounts = as.matrix(haber)))
		haberLabels = splitGetFromEnd(colnames(haber), "_", 1)
		haberType = c("Enterocyte", "Stem", "Goblet", "TA", "Enterocyte.Progenitor.Late", "Endocrine", "Enterocyte.Progenitor.Early", "Paneth", "Tuft.1", "Tuft.2", "Tuft")
		haberTypeName = c("Enterocyte", "Intestinal Stem cell", "Goblet", "TA", "Enterocyte Progenitor Late", "Endocrine", "Enterocyte Progenitor Early", "Paneth", "Tuft-1", "Tuft-2", "Tuft")
		haberLabels = mapvalues(haberLabels, haberType, haberTypeName)
		
		sceG = SingleCellExperiment(assays = list(counts = seuratObj@assays$RNA@counts))
		sceG = sceG[, colSums(counts(sceG)) > 0] # Remove libraries with no counts.
		sceG = logNormCounts(sceG) 
		sceG = sceG


		numcores = 1
		#bp = MulticoreParam(numcores)
		bp = SerialParam()
		singlerH = SingleR(test=sceG, ref=haberG, labels=haberLabels, de.method="wilcox", BPPARAM=bp)
		singlerM = SingleR(test=sceG, ref=immgen, labels=immgen$label.main, BPPARAM=bp)
		singlerF = SingleR(test=sceG, ref=immgen, labels=immgen$label.fine, BPPARAM=bp)
		singlerM2 = SingleR(test=sceG, ref=mousernaseq, labels=mousernaseq$label.main, BPPARAM=bp)
		singlerF2 = SingleR(test=sceG, ref=mousernaseq, labels=mousernaseq$label.fine, BPPARAM=bp)

		singler = list("ImmGen_Main"=singlerM, "ImmGen_Fine"=singlerF, "MouseRNAseq_Main"=singlerM2, "MouseRNAseq_Fine"=singlerF2, "Haber"= singlerH)
		saveRDS(singler, "singleR_result.rds")
		singlerDF = list("ImmGen_Main"=as.data.frame(singlerM), "ImmGen_Fine"=as.data.frame(singlerF), "MouseRNAseq_Main"=as.data.frame(singlerM2), "MouseRNAseq_Fine"=as.data.frame(singlerF2), "Haber"= as.data.frame(singlerH ))
		saveRDS(singlerDF, "singleR_result_df.rds")

	}else
	{
		singler = readRDS("singleR_result_df.rds")
	}
	singlerImmgenMain = addRownameColumn(singler[["ImmGen_Main"]], "Cell")[, c("Cell", "labels")]
	singlerImmgenFine = addRownameColumn(singler[["ImmGen_Fine"]], "Cell")[, c("Cell", "labels")]
	singlerMouseRNAseqMain = addRownameColumn(singler[["MouseRNAseq_Main"]], "Cell")[, c("Cell", "labels")]
	singlerMouseRNAseqFine = addRownameColumn(singler[["MouseRNAseq_Fine"]], "Cell")[, c("Cell", "labels")]
	singlerHaber = addRownameColumn(singler[["Haber"]], "Cell")[, c("Cell", "labels")]
	
	colnames(singlerImmgenMain) = c("Cell", "CellType_ImmGenMain")
	colnames(singlerImmgenFine) = c("Cell", "CellType_ImmGenFine")
	colnames(singlerMouseRNAseqMain) = c("Cell", "CellType_MouseRNAseqMain")
	colnames(singlerMouseRNAseqFine) = c("Cell", "CellType_MouseRNAseqFine")
	colnames(singlerHaber) = c("Cell", "CellType_Haber")
	
	
	dimred2 = dimred
	dimred2 = merge(dimred2, singlerImmgenMain, by="Cell", all.x=T)
	dimred2 = merge(dimred2, singlerImmgenFine, by="Cell", all.x=T)
	dimred2 = merge(dimred2, singlerMouseRNAseqMain, by="Cell", all.x=T)
	dimred2 = merge(dimred2, singlerMouseRNAseqFine, by="Cell", all.x=T)
	dimred2 = merge(dimred2, singlerHaber, by="Cell", all.x=T)


	plotSymbols(dimred2, "UMAP_1", "UMAP_2", "CellType_ImmGenMain", do.letters=T)
	ggsave("singler_umap_immgen_main.pdf", width=dimredPlotWidth, height=dimredPlotHeight)
	ggsave("singler_umap_immgen_main.png", width=dimredPlotWidth, height=dimredPlotHeight, dpi=600)
	
	plotSymbols(dimred2, "UMAP_1", "UMAP_2", "CellType_ImmGenFine", do.letters=T)
	ggsave("singler_umap_immgen_fine.pdf", width=dimredPlotWidth, height=dimredPlotHeight)
	ggsave("singler_umap_immgen_fine.png", width=dimredPlotWidth, height=dimredPlotHeight, dpi=600)
	
	plotSymbols(dimred2, "UMAP_1", "UMAP_2", "CellType_MouseRNAseqMain", do.letters=T)
	ggsave("singler_umap_mousernaseq_main.pdf", width=dimredPlotWidth, height=dimredPlotHeight)
	ggsave("singler_umap_mousernaseq_main.png", width=dimredPlotWidth, height=dimredPlotHeight, dpi=600)
	
	plotSymbols(dimred2, "UMAP_1", "UMAP_2", "CellType_MouseRNAseqFine", do.letters=T)
	ggsave("singler_umap_mousernaseq_fine.pdf", width=dimredPlotWidth, height=dimredPlotHeight)
	ggsave("singler_umap_mousernaseq_fine.png", width=dimredPlotWidth, height=dimredPlotHeight, dpi=600)

	plotSymbols(dimred2, "UMAP_1", "UMAP_2", "CellType_Haber", do.letters=T)
	ggsave("singler_umap_haber.pdf", width=dimredPlotWidth, height=dimredPlotHeight)
	ggsave("singler_umap_haber.png", width=dimredPlotWidth, height=dimredPlotHeight, dpi=600)

	
	
	for(curCol in c( "CellType_ImmGenMain",  "CellType_ImmGenFine", "CellType_MouseRNAseqMain", "CellType_MouseRNAseqFine", "CellType_Haber"))
	{
		celltype_cluster = as.data.frame(table( dimred2[, curCol], dimred2$Cluster ))
		celltype_sample = as.data.frame(table( dimred2[, curCol], dimred2$Type ))

		celltype_cluster = dcast(celltype_cluster, Var1~Var2, value.var = "Freq")
		celltype_sample = dcast(celltype_sample, Var1~Var2, value.var = "Freq")

		celltype_cluster_norm = celltype_cluster
		celltype_sample_norm = celltype_sample

		celltype_cluster_norm[2:ncol(celltype_cluster_norm)] = celltype_cluster_norm[2:ncol(celltype_cluster_norm)] / rep.row(colSums(as.matrix(celltype_cluster_norm[2:ncol(celltype_cluster_norm)])), nrow(celltype_cluster_norm))
		celltype_sample_norm[2:ncol(celltype_sample_norm)] = celltype_sample_norm[2:ncol(celltype_sample_norm)] / rep.row(colSums(as.matrix(celltype_sample_norm[2:ncol(celltype_sample_norm)])), nrow(celltype_sample_norm))


		celltype_sample_norm_filt = celltype_sample_norm[ rowMaxs(as.matrix(celltype_sample_norm[,2:ncol(celltype_sample_norm)])) > 0.005, ]

		write.tsv(celltype_cluster, paste0("SingleR_CellType_ByCluster_", curCol, ".txt"), row.names = F)
		write.tsv(celltype_sample, paste0("SingleR_CellType_BySample_", curCol, ".txt"), row.names = F)

		write.tsv(celltype_cluster_norm, paste0("SingleR_CellType_ByCluster_Norm_", curCol, ".txt"), row.names = F)
		write.tsv(celltype_sample_norm, paste0("SingleR_CellType_BySample_Norm_", curCol, ".txt"), row.names = F)

		celltype_sample_norm2 = celltype_sample_norm[ rowMeans(celltype_sample_norm[, 2:ncol(celltype_sample_norm)]) >= 0.001, ]
		celltype_sample_norm2 = melt(celltype_sample_norm2, by = "Var1")
		colnames(celltype_sample_norm2) = c("CellType", "Sample", "Proportion")
		celltype_sample_norm2$Type = splitGet(as.character(celltype_sample_norm2$Sample), "_", 1)
		celltype_sample_norm2$Type = factor(celltype_sample_norm2$Type, levels=levels(dimred$Type))

		ggplot(celltype_sample_norm2, aes(x = Type, y = Proportion, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~CellType, scales = "free") + expand_limits(y = c(0))
		ggsave(paste0("SingleR_Box_Type_BySample_", curCol, ".png"), width = 20, height = 11)
		
		ggplot(celltype_sample_norm2, aes(x = Type, y = Proportion, fill = Type)) + geom_bar(stat = "identity") + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~CellType, scales = "free") + expand_limits(y = c(0))
		ggsave(paste0("SingleR_Bar_Type_BySample_", curCol, ".png"), width = 20, height = 11)

		celltype_sample_norm3 = melt(celltype_sample_norm_filt, by = "Var1")
		colnames(celltype_sample_norm3) = c("CellType", "Sample", "Proportion")
		ggplot(celltype_sample_norm3, aes(x = Sample, y = Proportion, fill = CellType)) + geom_bar(stat = "identity") + theme_cem + scale_fill_manual(values = iwanthue(42))
		ggsave(paste0("SingleR_Bar_CellType_BySample_", curCol, ".png"), width = 15, height = 7)
	}
	
}



## Plot clusters on umap embedding
{
	library(dplyr)
	
	source("~/R/cem_polyOutline.R")
	clusterMedian = dimred %>%
						group_by(Cluster) %>%
						summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
	clusterMedian = as.data.frame(clusterMedian)
	hulls = concaveHull(dimred, "UMAP_1", "UMAP_2", "Cluster", alpha = 0.5, extend = T, minPoint = 20)
	
	numClust = length(unique(dimred$Cluster))
	numBatches = length(unique(dimred$Batch))
	numSamples = length(unique(dimred$Sample))
	numTypes = length(unique(dimred$Type))
	
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point() + theme_cem + geom_polygon(data = hulls, aes(x = x, y = y, fill = as.factor(Cluster), col = as.factor(Cluster)), alpha = 0.3) + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#444444") + theme(legend.position = "none")  + scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) + ggtitle(paste0(curTitle, ", clusters"))
	ggsave(paste0(curFilename, "dimred_cluster_hull.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_cluster_hull.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
	
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point() + theme_cem + geom_polygon(data = hulls, aes(x = x, y = y, fill = as.factor(Cluster), col = as.factor(Cluster)), alpha = 0.3) + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#444444") + theme(legend.position = "none")  + scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) + ggtitle(paste0(curTitle, ", clusters")) + facet_wrap(~Cluster)
	ggsave(paste0(curFilename, "dimred_cluster_hull_facet.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_cluster_hull_facet.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point(shape = 21, col = "black") + theme_cem + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#444444") + theme(legend.position = "none")  + scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) + ggtitle(paste0(curTitle, ", clusters"))
	ggsave(paste0(curFilename, "dimred_cluster.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_cluster.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Type), color = as.factor(Type))) + geom_point() + theme_cem  + scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + ggtitle(paste0(curTitle, ", by Type"))
	ggsave(paste0(curFilename, "dimred_type.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_type.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	nc = if(numTypes <= 3) numTypes else ceiling(sqrt(numTypes))
	nr = ceiling(numTypes/nc)
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Type), color = as.factor(Type))) + geom_point() + theme_cem  + scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + ggtitle(paste0(curTitle, ", by Type")) + facet_wrap(~Type, ncol=nc)
	ggsave(paste0(curFilename, "dimred_type_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
	ggsave(paste0(curFilename, "dimred_type_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
	
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Sample), color = as.factor(Sample))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numSamples)) + scale_color_manual(values = iwanthue(numSamples)) + ggtitle(paste0(curTitle, ", by Sample"))
	ggsave(paste0(curFilename, "dimred_sample.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_sample.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	nc = if(numSamples <= 3) numSamples else ceiling(sqrt(numSamples))
	nr = ceiling(numSamples/nc)
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Sample), color = as.factor(Sample))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numSamples)) + scale_color_manual(values = iwanthue(numSamples)) + ggtitle(paste0(curTitle, ", by Sample")) + facet_wrap(Sample~., ncol=nc)
	ggsave(paste0(curFilename, "dimred_sample_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
	ggsave(paste0(curFilename, "dimred_sample_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 

	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Batch), color = as.factor(Batch))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numBatches)) + scale_color_manual(values = iwanthue(numBatches)) + ggtitle(paste0(curTitle, ", by Batch"))
	ggsave(paste0(curFilename, "dimred_batch.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_batch.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	nc = if(numBatches <= 3) numBatches else ceiling(sqrt(numBatches))
	nr = ceiling(numBatches/nc)
	ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Batch), color = as.factor(Batch))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numBatches)) + scale_color_manual(values = iwanthue(numBatches)) + ggtitle(paste0(curTitle, ", by Batch")) + facet_wrap(~Batch, ncol=nc)
	ggsave(paste0(curFilename, "dimred_batch_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
	ggsave(paste0(curFilename, "dimred_batch_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
	

}

## Cluster frequencies
{
	clustCounts = dimred %>% group_by(Cluster, Type)  %>% summarise(Count=n())
	clustCounts = dcast(clustCounts, Cluster~Type, value.var = "Count")
	clustCounts$Cluster = factor(clustCounts$Cluster, levels=sort(unique(as.character(clustCounts$Cluster))))
	clustCounts[ is.na(clustCounts) ] = 0

	clustCountsNorm = clustCounts
	clustCountsNorm[2:ncol(clustCountsNorm)] = clustCountsNorm[2:ncol(clustCountsNorm)] / rep.row(colSums(clustCountsNorm[2:ncol(clustCountsNorm)]), nrow(clustCountsNorm))

	clustCountsNorm2 = melt(clustCountsNorm, by = "Cluster")
	clustCountsNorm2$Type = clustCountsNorm2$variable
	
	
	#clustCountsNorm2$Cluster = #factor(clustCountsNorm2$Cluster, levels = seq(0, max(as.numeric(clustCountsNorm2$Cluster))))
	ggplot(clustCountsNorm2, aes(x = Type, y = value, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
	ggsave(paste0(curFilename, "_clusterByType_box.png"), width = 13, height = 11)
	ggsave(paste0(curFilename, "_clusterByType_box.pdf"), width = 13, height = 11)
	
	ggplot(clustCountsNorm2, aes(x = Type, y = value, fill = Type)) + geom_bar(stat="identity") + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
	ggsave(paste0(curFilename, "_clusterByType_bar2.png"), width = 13, height = 11)
	ggsave(paste0(curFilename, "_clusterByType_bar2.pdf"), width = 13, height = 11)
	
	clustCountsNormSum = summarySE(clustCountsNorm2, measurevar="value", groupvars=c("Cluster", "Type"))

	ggplot(clustCountsNormSum, aes(x=Cluster, y=value, fill=Type, group=Type)) + geom_bar(stat="identity", position=position_dodge(0.9)) + theme_cem  + scale_fill_manual(values=colorsType)  + scale_color_manual(values=colorsType) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")# + geom_point(data=clustCountsNorm2, shape=21, position=position_dodge(0.9)) + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, position=position_dodge(0.9)) 
	ggsave(paste0(curFilename, "_clusterByType_bar.png"), width = 13, height = 11)
	ggsave(paste0(curFilename, "_clusterByType_bar.pdf"), width = 13, height = 11)
	
	
	
	clustCountsNormByCluster = clustCounts
	clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)] = clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)] / rep.col(rowSums(clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)]), ncol(clustCountsNormByCluster))

	clustCountsNormByCluster2 = melt(clustCountsNormByCluster, by = "Cluster")
	clustCountsNormByCluster2$Type = clustCountsNorm2$variable
	
	clustCountsNormByCluster3 = clustCountsNormByCluster2 %>% group_by(Cluster, Type) %>% summarise(value=sum(value))
	
	ggplot(clustCountsNormByCluster2, aes(x = Type, y = value, fill = Type)) + geom_bar(stat="identity") + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
	ggsave(paste0(curFilename, "_clusterByTypePerCluster_bar.png"), width = 13, height = 11)
	ggsave(paste0(curFilename, "_clusterByTypePerCluster_bar.pdf"), width = 13, height = 11)
}



## Pathway scores
{
	
	geneOrth = getGeneOrth()
	pathways = getPathways(geneOrth)
	pathwayNames = gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathways))
	pathwaysSelected = pathwayNames
	#pathwaysSelected = pathways[ pathwayNames %in% c( "SCHUHMACHER_MYC_TARGETS*", "HALLMARK_MYC_TARGETS_V1", "DECP_vs_DECN*", "GCB_GFP_MYC*") ]
	#pathwaysSelected = pathways[ grep("Haber", pathwayNames) ]
	
	pathwayScoresZscore = calculatePathwayScores(seuratObj, dimred, pathways, method="zscore")
	
	
	
	
	curPath = getwd()
	plotPath = paste0(curPath, "/PathwayPlots")
	dir.create(plotPath, showWarnings = F)
	dir.create(paste0(plotPath, "/umap_pointAll"), showWarnings = F)
	dir.create(paste0(plotPath, "/umap_pointByType"), showWarnings = F)
	dir.create(paste0(plotPath, "/umap_densityAll"), showWarnings = F)
	dir.create(paste0(plotPath, "/umap_densityByType"), showWarnings = F)
	dir.create(paste0(plotPath, "/umap_densityDifference"), showWarnings = F)
	dir.create(paste0(plotPath, "/byCluster_density"), showWarnings = F)
	setwd(plotPath)
	
	pathwayColors = viridis::viridis(10)
	pathwayColorsDiff = brewer.pal(9, "PuOr")
	
	for(curPathway in unique(pathwayScoresZscore$Pathway))
	{	
		curPathwayScore = pathwayScoresZscore[pathwayScoresZscore$Pathway %in% curPathway, ]
		if(nrow(curPathwayScore) == 0 | sum(is.nan(curPathwayScore$Score) > 0)) next
		
		rangeL = quantile(curPathwayScore$Score, 0.01, na.rm = T)
		rangeH = quantile(curPathwayScore$Score, 0.95, na.rm = T)
		if(rangeL == rangeH) { rangeL = min(curPathwayScore$Score); rangeH = max(curPathwayScore$Score); }
		
		ggplot(curPathwayScore, aes(x=UMAP_1, y=UMAP_2, fill=Score, color=Score)) + geom_point(alpha=0.7) + theme_cem + scale_color_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + scale_fill_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + ggtitle(paste0(curPathway, "\nModule score, z-score"))
		ggsave(paste0("umap_pointAll/", curFilename, "pathwayScoreZscore_umap_pointAll_", curPathway, ".png"), width=dimredPlotWidth, height=dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_pointAll/", curFilename, "pathwayScoreZscore_umap_pointAll_", curPathway, ".pdf"), width=dimredPlotWidth, height=dimredPlotHeight)
		
		ggplot(curPathwayScore, aes(x=UMAP_1, y=UMAP_2, fill=Score, color=Score)) + geom_point(alpha=0.7) + theme_cem + scale_color_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + scale_fill_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + facet_grid(~Type)
		ggsave(paste0("umap_pointByType/", curFilename, "pathwayScoreZscore_umap_pointByType_", curPathway, ".png"), width=dimredPlotWidth*3-1, height=dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_pointByType/", curFilename, "pathwayScoreZscore_umap_pointByType_", curPathway, ".pdf"), width=dimredPlotWidth*3-1, height=dimredPlotHeight)
		

		sm = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, numGrid=100, knn=50, m=2)
		ggplot(sm) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = pathwayColors) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0("umap_densityAll/", curFilename, "pathwayScoreZscore_umap_densityAll_", curPathway, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_densityAll/", curFilename, "pathwayScoreZscore_umap_densityAll_", curPathway, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
		
		
		sm2 = ddply(curPathwayScore, "Type", function(x) { smoothScore2d( x$Score, x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(curPathwayScore$UMAP_1), yrng=range(curPathwayScore$UMAP_2)) } )	
		ggplot(sm2) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = pathwayColors) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666") + facet_grid(~Type)
		ggsave(paste0("umap_densityByType/", curFilename, "pathwayScoreZscore_umap_densityByType_", curPathway, ".png"), width = dimredPlotWidth*3-1, height = dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_densityByType/", curFilename, "pathwayScoreZscore_umap_densityByType_", curPathway, ".pdf"), width = dimredPlotWidth*3-1, height = dimredPlotHeight)	
		
		
		#rangeL = quantile(curPathwayScore$Score, 0.001, na.rm = T)
		#rangeH = quantile(curPathwayScore$Score, 0.999, na.rm = T)
		#ggplot(curPathwayScore, aes_string(x = "Score", fill = "Type", color = "Type")) + geom_density(alpha=0.5) + facet_grid(Cluster~.) + theme_cem + ggtitle(paste0(curPathway, "\nModule score, z-score by cluster"))+ scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + coord_trans(limx=c(rangeL, rangeH))
		#ggsave(paste0("byCluster_density/", curFilename, "pathwayScoreZscore_byCluster_density_", curPathway, ".png"), width = dimredPlotHeight*1.5, height = dimredPlotWidth, dpi=150)
		#ggsave(paste0("byCluster_density/", curFilename, "pathwayScoreZscore_byCluster_density_", curPathway, ".pdf"), width = dimredPlotHeight*1.5, height = dimredPlotWidth)	
		
		##
		
		sm3 = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, type=curPathwayScore$Type, numGrid=100, knn=50, m=2)
		sm3W = dcast(sm3, x+y~type, value.var="score")
		sm3W$Diff = sm3W$AA - sm3W$Control
		rangeDiff = max(abs(sm3W$Diff))
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = Diff)) + theme_cem + scale_fill_gradientn(name = "AA vs Control difference", colors = pathwayColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curPathway, " AA vs Control Difference\nModule score, z-score\n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_AAvsControl.png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_AAvsControl.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
		
		##
		
		sm3 = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, type=curPathwayScore$Type, numGrid=100, knn=50, m=2)
		sm3W = dcast(sm3, x+y~type, value.var="score")
		sm3W$Diff = sm3W$PGE2 - sm3W$Control
		rangeDiff = max(abs(sm3W$Diff))
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = Diff)) + theme_cem + scale_fill_gradientn(name = "PGE2 vs Control difference", colors = pathwayColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curPathway, " PGE2 vs Control Difference\nModule score, z-score\n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_PGE2vsControl.png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_PGE2vsControl.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
		
		##
		
				
		sm3 = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, type=curPathwayScore$Type, numGrid=100, knn=50, m=2)
		sm3W = dcast(sm3, x+y~type, value.var="score")
		sm3W$Diff = sm3W$AA - sm3W$PGE2
		rangeDiff = max(abs(sm3W$Diff))
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = Diff)) + theme_cem + scale_fill_gradientn(name = "AA vs PGE2 difference", colors = pathwayColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curPathway, " AA vs PGE2 Difference\nModule score, z-score\n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_AAvsPGE2.png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_AAvsPGE2.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)	
		
		
		
	}
	setwd(curPath)


}



## Plot gene expression on dimred embedding
{
	assay = t(seuratObj@assays$integrated@scale.data)
	
	geneList = c("Ptger1", "Ptger2", "Ptger3", "Ptger4", "Ptgs1", "Ptgs2", "Ptges", "Ptges2", "Lgr5", "Olfm4", "Lyz1", "Dclk1", "Myc", "Tcf", "Lef", "Ppar", "Eph", "Il18", "Il33", "S100a6", "Pdk4", "Fabp1", "Cpt1a", "Defa21", "Defa23", "Defa35", "Defa25", "Defa30", "Defa22", "Defa31", "Defa3", "Defa5", "Defa27", "Defa29", "Defa2", "Defa33", "Defa36", "Defa32", "Defa37", "Defa28", "Defa26", "Defa17", "Defa34", "Defa24", "Muc15", "Muc1", "Muc3", "Muc3a", "Muc6", "Muc2", "Muc5ac", "Muc5b", "Muc16", "Muc19", "Mucl1", "Mucl2", "Muc4", "Muc20", "Muc13")
	
	
	curPath = getwd()
	plotPath = paste0(curPath, "/GenePlots")
	dir.create(plotPath)
	setwd(plotPath)
	
	geneList = intersect(geneList, colnames(assay))
	dimredG = merge(dimred, as.matrix(assay[, geneList, drop=F]), by.x="Cell", by.y="row.names")
	geneColors = brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
	geneColorsDiff = brewer.pal(9, "PiYG")
	
	for(curGene in geneList)
	{
		rangeL = quantile(dimredG[, curGene], 0.01, na.rm = T)
		rangeH = quantile(dimredG[, curGene], 0.95, na.rm = T)

		if(rangeL == rangeH) { rangeL = min(dimredG[, curGene]); rangeH = max(dimredG[, curGene]); }
		
		ggplot(dimredG, aes_string(x = "UMAP_1", y = "UMAP_2", fill = curGene, color = curGene)) + geom_point(alpha=0.7) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + scale_color_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + ggtitle(paste0(curGene)) + scale_alpha(range = c(0.3, 1))
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_pointAll_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_pointAll_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
		
		ggplot(dimredG, aes_string(x = "UMAP_1", y = "UMAP_2", fill = curGene, color = curGene)) + geom_point(alpha=0.7) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + scale_color_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + ggtitle(paste0(curGene)) + scale_alpha(range = c(0.3, 1)) + facet_grid(Type~.)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_pointByType_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight*2+1, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_pointByType_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight*2+1)
		
		sm = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, numGrid=100, knn=50, m=2)
		ggplot(sm) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors) + ggtitle(paste0(curGene)) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityAll_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityAll_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
		
		
		sm2 = ddply(dimredG, "Type", function(x) { smoothScore2d( x[,  curGene], x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(dimredG$UMAP_1), yrng=range(dimredG$UMAP_2)) } )	
		ggplot(sm2) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors) + ggtitle(paste0(curGene)) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666") + facet_grid(~Type)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityByType_", curGene, ".png"), width = dimredPlotWidth*3-1, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityByType_", curGene, ".pdf"), width = dimredPlotWidth*3-1, height = dimredPlotHeight)	
		
		
		sm3 = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, type=dimredG$Type, numGrid=100, knn=50, m=2)
		sm3W = dcast(sm3, x+y~type, value.var="score")
		sm3W$AA_vs_Control = sm3W$AA - sm3W$Control
		sm3W$PGE2_vs_Control = sm3W$PGE2 - sm3W$Control
		sm3W$AA_vs_PGE2 = sm3W$AA - sm3W$PGE2
		rangeDiff = max(c(abs(sm3W$AA_vs_Control), abs(sm3W$PGE2_vs_Control), abs(sm3W$AA_vs_PGE2)))
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = AA_vs_Control)) + theme_cem + scale_fill_gradientn(name = "AA vs Control", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curGene, ", AA vs Control difference \n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_AAvsControl_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_AAvsControl_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)	
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = PGE2_vs_Control)) + theme_cem + scale_fill_gradientn(name = "PGE2 vs Control", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curGene, ", PGE2 vs Control difference \n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_PGE2vsControl_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_PGE2vsControl_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)	
		
		ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = AA_vs_PGE2)) + theme_cem + scale_fill_gradientn(name = "AA vs PGE2", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) + ggtitle(paste0(curGene, ", AA vs PGE2 difference \n(calculations are experimental, confounded by population differences)")) + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#666666")
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_AAvsPGE2_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		ggsave(paste0(curFilename, "dimred_geneExpression_umapAll_densityDifference_AAvsPGE2_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)	
		
		#
		#ggplot(dimredG, aes_string(x = curGene, fill = "Type", color = "Type")) + geom_density(alpha=0.5) + facet_wrap(Cluster~.) + theme_cem + ggtitle(paste0(curGene, " expression by cluster"))+ scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType)
		#ggsave(paste0(curFilename, "dimred_geneExpression_byCluster_density_", curGene, ".png"), width = dimredPlotHeight, height = dimredPlotWidth, dpi=150)
		#ggsave(paste0(curFilename, "dimred_geneExpression_byCluster_density_", curGene, ".pdf"), width = dimredPlotHeight, height = dimredPlotWidth)	
		#
		#VlnPlot(seuratObj, features = curGene, split.by = "Type", pt.size=0.3, cols = colorsType[order(names(colorsType))])
		#ggsave(paste0(curFilename, "dimred_geneExpression_byCluster_violin_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
		#ggsave(paste0(curFilename, "dimred_geneExpression_byCluster_violin_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
	}
	
	setwd(curPath)
}


## Markers & cluster specific genes
{
	DotPlot(seuratObj, features = geneList, cols =  geneColors[c(1, length(geneColors)-1)]) + RotatedAxis()
	ggsave(paste0(curFilename, "_SeuratPlots_markerDot.png"), width=12, height=8)
	ggsave(paste0(curFilename, "_SeuratPlots_markerDot.pdf"), width=12, height=8)
	markersAll = FindAllMarkers(seuratObj, thresh.test = 3, test.use = "roc", do.print = TRUE)
	markersFiltered = subset(markersAll, avg_diff > 0 & power > 0.6)
	
	write.tsv(markersAll, "seurat_clusterMarkers_all.txt")
	write.tsv(markersFiltered, "seurat_clusterMarkers_filtered.txt")
	
	DoHeatmap(seuratObj, features = markersFiltered$gene)
	ggsave(paste0(curFilename, "_SeuratPlots_markerHeatmap.png"), width=18, height=14)
	ggsave(paste0(curFilename, "_SeuratPlots_markerHeatmap.pdf"), width=18, height=14)

	
	markersFiltered2 = markersAll %>% filter(avg_diff > 0) 
	for(curClust in unique(markersFiltered2$cluster))
	{
		DoHeatmap(seuratObj, features = markersFiltered2$gene[ markersFiltered2$cluster %in% curClust])
		ggsave(paste0(curFilename, "_SeuratPlots_markerHeatmap_cluster", curClust, ".png"), width=18, height=8)
		ggsave(paste0(curFilename, "_SeuratPlots_markerHeatmap_cluster", curClust, ".pdf"), width=18, height=8)
	}

	for(curCluster in unique(dimred$Cluster))
	{
		ggplot() +  guides(colour = FALSE) + theme_cem +
		  geom_path(data = hulls[ hulls$Cluster != curCluster, ], aes(x = x, y = y, group = Cluster), alpha = 0.3, col = "#000000") +
		  geom_path(data = hulls[ hulls$Cluster == curCluster, ], aes(x = x, y = y, group = Cluster), alpha = 1, size=2, col = "#ff0000") + 
		  geom_point(data = dimred[ dimred$Cluster == curCluster, ], aes(x = UMAP_1, y = UMAP_2), color = "#ff0000", alpha=0.6) + 
		  geom_text(data = clusterMedian[ clusterMedian$Cluster != curCluster, ], aes(label = Cluster, group = as.factor(Cluster), x=UMAP_1, y=UMAP_2), size = 6, family = "mono", fontface = "bold", color = "#888888") + 
		  geom_text(data = clusterMedian[ clusterMedian$Cluster == curCluster, ], aes(label = Cluster, group = as.factor(Cluster), x=UMAP_1, y=UMAP_2), size = 8, family = "mono", fontface = "bold", color = "#000000") 
		ggsave(paste0("ClusterOutline_dimred_", curCluster, ".png"), width=11, height=9, dpi=150)
	}

}

## Differential markers between types
{

	seuratObjAll = seuratObj
	seuratObjAll$ClusterType = paste(Idents(seuratObjAll), seuratObjAll$Type, sep = "_")
	seuratObjAll$Cluster = Idents(seuratObjAll)
	Idents(seuratObjAll) = "ClusterType"

	allMarkers = data.frame()
	for(curClust in unique(seuratObjAll$Cluster))
	{
		for(curCmp1 in c("AA", "PGE2"))
		{
			for(curCmp2 in c("PGE2", "Control"))
			{
				if(curCmp1 == curCmp2) next
				
				id1 = paste0(curClust, "_", curCmp1)
				id2 = paste0(curClust, "_", curCmp2)
				
				if( sum(seuratObjAll$ClusterType %in% id1) >= 5 & sum(seuratObjAll$ClusterType %in% id2) >= 5 )
				{
					print(paste0(id1, " ", id2))
					curMarkers = FindMarkers(seuratObjAll, ident.1 = id1, ident.2 = id2, verbose = FALSE)
					curMarkers = addRownameColumn(curMarkers, "Gene")
					allMarkers = rbind(allMarkers, data.frame(Cluster=curClust, Group1=curCmp1, Group2=curCmp2, curMarkers))
				}
			}
		}
	}
	write.tsv(allMarkers, "TypeComparison_byCluster.txt", row.names=F)
	allMarkersSig = allMarkers[ which(allMarkers$p_val_adj < 0.01), ]
	allMarkersSig = allMarkersSig[ order(allMarkersSig$p_val_adj), ]

}




## Plot the density difference
{
	library(MASS)
	library(reshape2)
	library(scales)

	# Calculate the common x and y range for geyser1 and geyser2
	xrng = range(dimred$UMAP_1)
	yrng = range(dimred$UMAP_2)

	extendRange = 0.06 # extend range by %6
	xDiff = (xrng[2] - xrng[1])
	yDiff = (yrng[2] - yrng[1])
	
	xrng[1] = xrng[1] - xDiff* extendRange
	xrng[2] = xrng[2] + xDiff * extendRange
	yrng[1] = yrng[1] - yDiff * extendRange
	yrng[2] = yrng[2] + yDiff * extendRange
	
	
	
	for(caseType in unique(dimred$Type))
	{
		for(ctrlType in unique(dimred$Type))
		{
			if(caseType == ctrlType) next
			
			caseVsCtrlName = paste0(caseType, "_vs_", ctrlType)
			
			colorLow = colorsType[ctrlType]
			colorMid = "white"
			colorHigh = colorsType[caseType]
			
			d_Case = kde2d( dimred$UMAP_1[dimred$Type == caseType], dimred$UMAP_2[dimred$Type == caseType ], lims = c(xrng, yrng), n = 500)
			d_Ctrl = kde2d( dimred$UMAP_1[dimred$Type == ctrlType  ], dimred$UMAP_2[dimred$Type == ctrlType ], lims = c(xrng, yrng), n = 500)

			# Confirm that the grid points for each density estimate are identical
			identical(d_Case$x, d_Ctrl$x) # TRUE
			identical(d_Case$y, d_Ctrl$y) # TRUE

			# Calculate the difference between the 2d density estimates
			diff_CaseVsCtrl = d_Ctrl 
			diff_CaseVsCtrl$z = d_Case$z - d_Ctrl$z
			diff_CaseVsCtrl$z = diff_CaseVsCtrl$z / max(diff_CaseVsCtrl$z)

			## Melt data into long format
			# First, add row and column names (x and y grid values) to the z-value matrix
			rownames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$x
			colnames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$y

			# Now melt it to long format
			diff_CaseVsCtrlM = melt(diff_CaseVsCtrl$z, id.var = rownames(diff_CaseVsCtrl))
			names(diff_CaseVsCtrlM) = c("UMAP_1", "UMAP_2", caseVsCtrlName)

			# Plot difference between geyser2 and geyser1 density
			ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
				geom_tile(aes_string(fill = caseVsCtrlName), alpha = 1) +
				scale_fill_gradient2(low = colorLow, mid = colorMid, high = colorHigh, midpoint = 0) +
				coord_cartesian(xlim = xrng, ylim = yrng) +
				scale_color_manual(values = colorsType) +
				guides(colour = FALSE) + theme_cem +
				geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#444444") + 
				geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.3, fill=NA, color="grey50") + 
				ggtitle(paste0("Density comparison of ", caseType, " vs ", ctrlType)) 

			ggsave(paste0(curFilename, "dimred_densityDiff_cluster_", caseVsCtrlName, ".png"), width = dimredPlotWidth, height = dimredPlotHeight)
		}
	}
}


## Trajectory analysis
{

	library(Seurat)
	library(scales)
	library(dplyr)
	library(dyno)
	library(tidyverse)
	library(Matrix)
	library(pheatmap)
	
	seuratExp = t(seuratObj@assays$integrated@scale.data)
	seuratCt = t(seuratObj@assays$SCT@counts)
	commonCells = intersect(rownames(seuratExp), rownames(seuratCt))
	commonGenes = intersect(colnames(seuratExp), colnames(seuratCt))
	seuratExp2 = seuratExp[ commonCells, commonGenes ]
	seuratCt2 = seuratCt[ commonCells, commonGenes ]

	dataset = wrap_expression(
	  counts = seuratCt2 ,
	  expression = seuratExp2
	)

	ti_slingshot = dynwrap::create_ti_method_container("dynverse/ti_slingshot:latest") # dynmethods::ti_slingshot() 
	model = infer_trajectory(dataset, ti_slingshot(), verbose = TRUE)
	#modelMonocle = infer_trajectory(dataset, ti_monocle_ddrtree(), verbose = TRUE)


	
	overall_feature_importances = dynfeature::calculate_overall_feature_importance(model, expression_source = dataset$expression)
	#branching_milestone = model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()
	#branch_feature_importance = calculate_branch_feature_importance(model, expression_source=dataset$expression)
	#branchpoint_feature_importance = calculate_branching_point_feature_importance(model, expression_source=dataset$expression, milestones_oi = branching_milestone)



	topN = 100
	features_overall = overall_feature_importances %>% 
	  top_n(topN, importance) %>% 
	  pull(feature_id)
	  
	groupingCluster = factor(dimred$Cluster)
	names(groupingCluster) = dimred$Cell
	  
	groupingType = factor(dimred$Type)
	names(groupingType) = dimred$Cell
	

	plot_heatmap(
	  model, 
	  expression_source = dataset$expression, 
	  features_oi = features_overall,
	  grouping = groupingCluster
	)
	ggsave(paste0(curFilename, "_heatmapFeatureOverall.png"), width=12, height=14, dpi=150)
	ggsave(paste0(curFilename, "_heatmapFeatureOverall.pdf"), width=12, height=14)


	
	pseudotime = dynwrap::calculate_pseudotime(model)
	features_overall = overall_feature_importances[ !grepl("^Rp", overall_feature_importances$feature_id), ] %>% top_n(topN, importance) %>% pull(feature_id)
	#features_branch = unique( branch_feature_importance[ !grepl("^Rp", branch_feature_importance$feature_id), ] %>% top_n(topN, importance) %>%  pull(feature_id) )
	seuratExpSel = t(seuratExp[ names(pseudotime)[ order(pseudotime) ], features_overall ] )
	
	#clust = hclust(as.dist(dynutils::correlation_distance(seuratExpSel)), method = "ward.D2")
	#seuratExpSel = seuratExpSel[ rownames(seuratExpSel)[clust$order], ]
	geneOrder = getTSPOrder(seuratExpSel)
	seuratExpSel = seuratExpSel[geneOrder, ]
	seuratExpSelCtrl = seuratExpSel[, colnames(seuratExpSel) %in% dimred$Cell[ dimred$Type == "Control"] ]
	seuratExpSelAA = seuratExpSel[, colnames(seuratExpSel) %in% dimred$Cell[ dimred$Type == "AA"] ]
	seuratExpSelPGE2 = seuratExpSel[, colnames(seuratExpSel) %in% dimred$Cell[ dimred$Type == "PGE2"] ]
	
	
	annColumn = data.frame(row.names=names(pseudotime), Pseudotime=pseudotime, Type=splitGet(names(pseudotime), "_", 2))
	annColumn$Type = mapvalues(as.character(annColumn$Type), c("C", "A", "P"), c("Control", "AA", "PGE2"))
	annColors = list(Pseudotime = brewer.pal(11, "Spectral"), Type = colorsType)
	
	colScale = viridis::inferno(50)
	maxVal = max(seuratExpSel)
	minVal = min(seuratExpSel)
	breaks = c(minVal, seq(quantile(seuratExpSel, 0.2), quantile(seuratExpSel, 0.99), length.out = 49), maxVal)
	
	heatmapWidth = 18
	heatmapHeight = 13
	
	pheatmap(seuratExpSel, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_all.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAll Cells")
	pheatmap(seuratExpSel, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_all.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAll Cells")
	
	pheatmap(seuratExpSelCtrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_Control.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nControl Cells")
	pheatmap(seuratExpSelCtrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_Control.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nControl Cells")
	
	pheatmap(seuratExpSelAA, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_AA.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAA Cells")
	pheatmap(seuratExpSelAA, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_AA.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAA Cells")
	
	pheatmap(seuratExpSelPGE2, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_PGE2.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nPGE2 Cells")
	pheatmap(seuratExpSelPGE2, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_PGE2.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nPGE2 Cells")
	
	
	
	#############
	
	
	
	source("~/R/cem_polyOutline.R")
	
	dimredP = merge(dimred, annColumn[, c("Pseudotime"), drop=F], by.x="Cell", by.y="row.names")
	
	clusterMedian = dimred %>%
						group_by(Cluster) %>%
						summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
	clusterMedian = as.data.frame(clusterMedian)
	hulls = concaveHull(dimred, "UMAP_1", "UMAP_2", "Cluster", alpha = 0.5, extend = T, minPoint = 3)
	

	
	ddimred = data.frame(row.names=dimred$Cell, comp_1=dimred$UMAP_1, comp_2=dimred$UMAP_2)
	cell_positions = ddimred %>% as.data.frame() %>% rownames_to_column("cell_id")
	cell_positions = as.data.frame(cell_positions)
	cell_positions = left_join(cell_positions, model$milestone_percentages %>% 
									group_by(cell_id) %>% arrange(desc(percentage)) %>% 
									filter(dplyr::row_number() == 1) %>% select(cell_id, milestone_id), "cell_id")
	color_trajectory = "none"	
	waypoints = dynwrap::select_waypoints(model)
	trajectory_projection_sd = sum(model$milestone_network$length) * 0.05
	waypoint_projection = dynplot:::project_waypoints(trajectory = model, 
				cell_positions = cell_positions, waypoints = waypoints, 
				trajectory_projection_sd = trajectory_projection_sd, 
				color_trajectory = color_trajectory)

	arrow = if (any(model$milestone_network$directed))	arrow(type = "closed", length = (unit(0.1, "inches")))
		
	edges = as.data.frame(waypoint_projection$edges)
	edgesArrow = edges[ which(edges$arrow), ]
	ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime)) + theme_cem + scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + ggtitle(paste0(curTitle, ", UMAP by Pseudotime")) + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 7, family = "mono", fontface = "bold", color = "#444444") + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.2, fill=NA, color="grey50") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edgesArrow , arrow = arrow, color = "#333333", size = 1, linejoin = "mitre", lineend = "butt") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edges, size = 1, color = "#333333")
	ggsave(paste0(curFilename, "dimred_pseudotime_hull.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_pseudotime_hull.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)

	ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime)) + theme_cem + scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + ggtitle(paste0(curTitle, ", UMAP by Pseudotime")) + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edgesArrow, arrow = arrow, color = "#333333", size = 1, linejoin = "mitre", lineend = "butt") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), edges, size = 1, color = "#333333")
	ggsave(paste0(curFilename, "dimred_pseudotime.png"), width = dimredPlotWidth, height = dimredPlotHeight)
	ggsave(paste0(curFilename, "dimred_pseudotime.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
	
	
	
	############
	
	
	
	#seuratExpSel2 = seuratExpSel - rep.col(rowMins(seuratExpSel), ncol(seuratExpSel))
	seuratExpSelMed = as.data.frame(seuratExpSel)
	seuratExpSelMed$Gene = rownames(seuratExpSelMed)
	seuratExpSelMed = melt(seuratExpSelMed, "Gene", variable.name = "Cell", value.name = "Expression")
	seuratExpSelMed = merge(seuratExpSelMed, data.frame(Cell=names(pseudotime), Pseudotime=pseudotime), by="Cell")
	seuratExpSelMed = merge(seuratExpSelMed, dimred[, c("Cell", "Type")], by="Cell")
	#seuratExpSelMed$PseudotimeFactor = cut(seuratExpSelMed$Pseudotime, 100)
	seuratExpSelMed$PseudotimeFactor = Hmisc::cut2(seuratExpSelMed$Pseudotime, g=100)
	seuratExpSelMedSum1 = seuratExpSelMed %>% group_by(Gene, PseudotimeFactor) %>% summarise(Expression=median(Expression, na.rm=T))
	seuratExpSelMedSum2 = seuratExpSelMed %>% group_by(Gene, PseudotimeFactor, Type) %>% summarise(Expression=median(Expression, na.rm=T))
	seuratExpSelMedSumW = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum1, Gene~PseudotimeFactor, value.var="Expression")))
	seuratExpSelMedSumW_Ctrl = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum2[ seuratExpSelMedSum2$Type == "Control", ], Gene~PseudotimeFactor, value.var="Expression")))
	seuratExpSelMedSumW_Case = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum2[ seuratExpSelMedSum2$Type == "AA", ], Gene~PseudotimeFactor, value.var="Expression")))
	
	cellCounts = seuratExpSelMed %>% dplyr::select(Cell, Type, PseudotimeFactor) %>% distinct()  #%>% group_by(Type, PseudotimeFactor) %>% summarise(CellCount = n())
	cellCounts$Pseudotime = as.numeric(gsub("[\\[\\(\\ ]", "", splitGet(as.character(cellCounts$PseudotimeFactor), ",", 1)))
	
	#clust = hclust(as.dist( 1 - (cor(seuratExpSelMedSumW)+1)/2 ), method = "ward.D2")
	#clust = hclust(dist(seuratExpSelMedSumW), method = "ward.D2")
	#geneOrder = clust$order
	
	#geneOrder = getTSPOrder(seuratExpSelMedSumW)
	seuratExpSelMedSumW = seuratExpSelMedSumW[geneOrder, ]
	seuratExpSelMedSumW_Ctrl = seuratExpSelMedSumW_Ctrl[geneOrder, ]
	seuratExpSelMedSumW_Case = seuratExpSelMedSumW_Case[geneOrder, ]
	seuratExpSelMedSumW_Diff = seuratExpSelMedSumW_Case - seuratExpSelMedSumW_Ctrl
	
	annColumn2 = setRemoveRownames(unique(cellCounts[, c("PseudotimeFactor", "Pseudotime")]), 1)
	
	colScale = viridis::inferno(50)
	maxVal = max(seuratExpSelMedSumW)
	minVal = min(seuratExpSelMedSumW)
	breaks = c(minVal, seq(quantile(seuratExpSelMedSumW, 0.2), quantile(seuratExpSelMedSumW, 0.99), length.out = 49), maxVal)
	
	pheatmap(seuratExpSelMedSumW, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_all.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, All Cells")
	pheatmap(seuratExpSelMedSumW, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_all.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, All Cells")
	
	pheatmap(seuratExpSelMedSumW_Case, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_AA.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, AA Cells")
	pheatmap(seuratExpSelMedSumW_Case, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_AA.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, AA")
	
	
	pheatmap(seuratExpSelMedSumW_Ctrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_WT.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, WT Cells")
	pheatmap(seuratExpSelMedSumW_Ctrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_WT.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, WT Cells")
	


	colScale = colorRampPalette(rev(c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")))(50)
	maxVal = max(abs(seuratExpSelMedSumW_Diff))
	breaks = c(-maxVal, seq(-2, 2, length.out = 49), maxVal)
	
	annColumn3 = annColumn2[, "Pseudotime", drop = F]
	pheatmap(seuratExpSelMedSumW_Diff, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn3, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpressionDifference_AAvsCtrl.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAA - Ctrl median expression difference across pseudotime")
	pheatmap(seuratExpSelMedSumW_Diff, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn3, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpressionDifference_AAvsCtrl.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAA - Ctrl median expression difference across pseudotime")
		
	ggplot(cellCounts, aes(x=Pseudotime, fill=Type, color=Type)) + geom_density(alpha=0.7, adjust=0.5) + theme_cem + scale_fill_manual(values=colorsType) + scale_color_manual(values=colorsType)
	ggsave(paste0(curFilename, "_densityByPseudotimeFactor.pdf"),  width=heatmapWidth, height=heatmapHeight)
	ggsave(paste0(curFilename, "_densityByPseudotimeFactor.png"),  width=heatmapWidth, height=heatmapHeight, dpi=150)
	
	ggplot(annColumn, aes(x=Pseudotime, fill=Type, color=Type)) + geom_density(alpha=0.7, adjust=0.5) + theme_cem + scale_fill_manual(values=colorsType) + scale_color_manual(values=colorsType)
	ggsave(paste0(curFilename, "_densityByPseudotime.pdf"),  width=heatmapWidth, height=heatmapHeight)
	ggsave(paste0(curFilename, "_densityByPseudotime.png"),  width=heatmapWidth, height=heatmapHeight, dpi=150)
	
}



sort( sapply(ls(),function(x){object.size(get(x))}), decreasing = T) 

rm(	seuratObjAll, seuratAll, seuratTemp, seuratAll2, seuratAllTmp,  seuratAnchors, seuratObjList, 
	sceG, curSeurat, dataset, seuratCt, countSeurat, assay, seuratExp, seuratExp2, seuratCt2, ctDat, curDat, haberG,
	haber, singler, immgen, singlerF, seuratExpSelMed, singlerDF, mousernaseq) 

save.image("seurat_SCT_SB04SB05_umap_2020.07.26.RData")	
