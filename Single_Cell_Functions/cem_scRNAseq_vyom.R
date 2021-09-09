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
	pathwaysM = c(gmtPathways("~/analysis/scRNA_Intestine/Beyaz_AA_final.gmt"))
	return(pathwaysM)
}

calculatePathwayScores = function(seuratObj, dimred, pathwayList, method="seurat", minCells=50,  minGenes=10)
{
	curData = seuratObj@assays$SCT@scale.data
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

