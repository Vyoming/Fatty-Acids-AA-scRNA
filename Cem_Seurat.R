`:=` <-
function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}
`?` <-
function(x, y) {
  xs <- as.list(substitute(x))
  if (xs[[1]] == as.name("<-")) x <- eval(xs[[3]])
  r <- eval(sapply(strsplit(deparse(substitute(y)), ":"), function(e) parse(text = e))[[2 - as.logical(x)]])
  if (xs[[1]] == as.name("<-")) {
    xs[[3]] <- r
        eval.parent(as.call(xs))
  } else {
    r
  }
}
addFootnote <-
function(ggplotPlot, footnoteText, fontsize=15, fontface="italic", plot=T)
{
	library(grid)
	library(gridExtra)
	ggGrob = arrangeGrob(ggplotPlot, bottom = textGrob(footnoteText, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = fontface, fontsize = fontsize)))
	if(plot)
	{
		grid.draw(ggGrob)
	}
	class(ggGrob) = c("arrange","ggplot", class(ggGrob))
	return(ggGrob)
}
addRownameColumn <-
function(df, idColName="id")
{
	ids = rownames(df)
	rownames(df) = NULL
	#df = cbind(ids, df)
	df = data.frame(ids, df)
	colnames(df)[1] = idColName
	return(df)
}
anyNA <-
function(x) { any(is.na(x)) }
calculatePathwayScores <-
function(seuratObj, dimred, pathwayList, method="seurat", minCells=50,  minGenes=10)
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
cemPalette <-
function(x) 	{ if(x<=12) brewer.pal(x, "Paired") else (c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Pastel1")))[1:x]	}
colMax2nd <-
function (x)
{
	apply(x, 2, function(row) max(row[-which.max(row)], na.rm=T))
}
colMaxNth <-
function (x, N=2, returnIndex=FALSE)
{
	if(returnIndex)
	{
		return( apply(x, 2, function(row) maxNidx(row, N)) )
	}else
	{
		return( apply(x, 2, function(row) maxN(row, N)) )
	}
}
color.bar <-
function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
concaveHull <-
function(df, X, Y, group, alpha, margin=NA, extend=F, minPoint = 3, marginMultiplier = 1)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConcaveHull, X, Y, alpha, margin, extend, minPoint)
	#hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	# hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls)
}
concaveHullSmooth <-
function(df, X, Y, group, alpha, margin=NA, extend=F, minPoint = 3, marginMultiplier = 1)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConcaveHull, X, Y, alpha, margin, extend, minPoint)
	#hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls2)
}
convexHull <-
function(df, X, Y, group, margin=NA, marginMultiplier = 3)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConvexHull, X, Y, margin)
	hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	# hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls2)
}
cp <-
function(dat) { write.tsv(dat, "~/Rcopy.txt") }
custom_fortify <-
function(x, ...) {
  df <- fortify(x, ...)
  df %>% group_by(id) %>% do(fixfeature(.))
}
dcast <-
function(...)
{
	return(reshape2::dcast(...))
}
detachAllPackages <-
function() 
{
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
eval.string.dplyr <-
function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df  
}
extendHull <-
function(df, X, Y, margin=10)
{
	library(polyclip)
	extended = polyoffset(list(x=df[, X], y=df[, Y]), margin, jointype="round", arctol=abs(margin)/40)
	extended2 = data.frame(x=extended[[1]]$x, y=extended[[1]]$y)
	return(extended2)
}
facetAdjust <-
function(x, pos = c("up", "down"))
{
	library(grid)
	library(ggplot2)
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p); dev.off()
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  n <- space - panels
  if(panels != space){
    idx <- (space - ncol - n + 1):(space - ncol)
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  class(gtable) <- c("facetAdjust", "gtable", "ggplot", "gTree","grob","gDesc" ); gtable
}
findConcaveHull <-
function(df, X, Y, alpha, hullMargin=NA, extend=F, minPoint=3)
{
	library(igraph)
	library(alphahull)
	library(polyclip)
	
	if(nrow(df) < 3)
	{
		replicationMargin = 3*hullMargin/100
		df2 = rbind( data.frame(X=df[, X] + replicationMargin, Y=df[, Y] + replicationMargin), 
					 data.frame(X=df[, X] + replicationMargin, Y=df[, Y] - replicationMargin),
					 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] + replicationMargin),
					 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] - replicationMargin) )
	}else
	{
		df2 = data.frame(X=df[, X], Y=df[, Y])
	}
	
	if(minPoint <= 3) minPoint = 3
	
	if(is.na(hullMargin)) hullMargin = min( (max(df2$X, na.rm=T) - min(df2$X, na.rm=T) ), (max(df2$Y, na.rm=T) - min(df2$Y, na.rm=T)) ) * 0.05
	
	#shape = ashape(df2$X, df2$Y, alpha)
	
	##plot(shape)
	#shapeGraph = graph.edgelist(cbind(as.character(shape$edges[, "ind1"]), as.character(shape$edges[, "ind2"])), directed = FALSE)
	##plot(shapeGraph)

	shape2 = ahull(df2$X, df2$Y, alpha)	
	shapeArcs = shape2$arcs
	shapeArcs = shapeArcs[ shapeArcs[, "end1"] != shapeArcs[, "end2"], ]
	shapeGraph = graph.edgelist(cbind(as.character(shapeArcs[, "end1"]), as.character(shapeArcs[, "end2"])), directed = FALSE)

	graphComponents = decompose.graph(shapeGraph)
	edgesList = list()
	curIdx = 1
	for(curGraph in graphComponents)
	{
		if(length(V(curGraph)) < minPoint) next
		if(length(E(curGraph)) < minPoint) next
		#if (!is.connected(curGraph)) next
		#if (any(degree(curGraph) != 2)) next
		
		# find chain end points
		for(i in 1:length(E(curGraph)))
		{
			cutg = curGraph - E(curGraph)[i]
			ends = names(which(degree(cutg) == 1))
			if(length(ends) >= 2) break
		}
		
		if(length(ends) >= 2)
		{
			path = get.shortest.paths(cutg, ends[1], ends[2])[[1]]
			# this is an index into the points
			pathX = as.numeric(V(curGraph)[path[[1]]]$name)
			# join the ends
			pathX = c(pathX, pathX[1])
			edgesList[[curIdx]] = shape2$x[pathX, ]
			curIdx = curIdx + 1
		}else
		{
			# TODO: Edge case
			# Graph is shaped weirdly, do depth-first-traversal or similar way
			# dfs
		}
	}

	if(length(edgesList)==1)
	{
		curEdges = edgesList[[1]]
		allEdges = as.data.frame(curEdges)
		colnames(allEdges) = c("x", "y")
		allEdges$group = 1
		if(extend == T)
		{
			extended = polyoffset(list(x=allEdges[, 1], y=allEdges[, 2]), hullMargin, jointype="round", arctol=abs(hullMargin)/40)
			extended2 = data.frame(x=extended[[1]]$x, y=extended[[1]]$y, group=1)
			return(extended2)
		}else
		{
			return(allEdges)
		}
	}else if(length(edgesList) > 1)
	{
		mergedEdges = edgesList[[1]]
		mergedEdges = mergedEdges[1:(nrow(mergedEdges)-1), ]
		otherEdges = list()
		for(i in 2:length(edgesList))
		{
			curEdges = edgesList[[i]]
			curEdges = curEdges[1:(nrow(curEdges)-1), ]
			otherEdges[[i-1]] = list(x=curEdges[,1], y=curEdges[,2])
		}

		mergedShape = polyclip(list(x=mergedEdges[,1], y=mergedEdges[,2]), otherEdges, op="xor")
			
		mergedShape2 = data.frame()
		for(j in 1:length(mergedShape))
		{
			mergedShape[[j]]$x = c(mergedShape[[j]]$x, mergedShape[[j]]$x[1]) ## Extend the last point, otherwise ggplot freaks out with holes
			mergedShape[[j]]$y = c(mergedShape[[j]]$y, mergedShape[[j]]$y[1])
			if(extend == F)
			{
				mergedShape2 = rbind(mergedShape2, data.frame(x=mergedShape[[j]]$x, y=mergedShape[[j]]$y, group=j))
			}else
			{			
				mergedShapeEx = polyoffset(mergedShape[[j]], hullMargin/2, jointype="round", arctol=abs(hullMargin)/40)
				mergedShapeDF = data.frame(x=mergedShapeEx[[1]]$x, y=mergedShapeEx[[1]]$y, group=j)
				#mergedShapeDF = simplifyDF(mergedShapeDF, "x", "y", hullMargin, T)
				mergedShapeDF = data.frame(x=mergedShapeDF$x, y=mergedShapeDF$y, group=j)
				mergedShape2 = rbind(mergedShape2, mergedShapeDF)
			}
		}

		allEdges = as.data.frame(mergedShape2)
		allEdges = fixfeature(allEdges)
		#ggplot(allEdges, aes(x=x, y=y, group=group)) + geom_polygon()
		
		if(extend == T)
		{
			extended = polyoffset(list(x=allEdges$x, y=allEdges$y), hullMargin/2, jointype="round", arctol=abs(hullMargin)/40)
			extended = data.frame(x=extended[[1]]$x, y=extended[[1]]$y, group=1)
			return(extended)
		}else
		{
			return(data.frame(x=allEdges$x, y=allEdges$y, group=allEdges$group))
		}
	}else
	{
		return(data.frame(x=numeric(0), y=numeric(0), group=numeric(0)))
	}
}
findConvexHull <-
function(df, X, Y, margin)
{
	if(nrow(df) < 3)
	{
		replicationMargin = margin/100
		df2 = rbind( data.frame(X=df[, X] + replicationMargin, Y=df[, Y] + replicationMargin), 
					 data.frame(X=df[, X] + replicationMargin, Y=df[, Y] - replicationMargin),
					 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] + replicationMargin),
					 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] - replicationMargin) )
	}else
	{
		df2 = data.frame(X=df[, X], Y=df[, Y])
	}
	hullIdx = chull(df2$X, df2$Y)
	hulls = df2[hullIdx, ]
	return(hulls)
}
fixfeature <-
function(df) {
  ringstarts <- which(!duplicated(df$group))
  if(length(ringstarts) < 2) {
    return(df)
  } else {
    ringstarts <- c(ringstarts, nrow(df))
    indicies <- c(1:(ringstarts[2]-1), do.call(c, lapply(2:(length(ringstarts)-1), function(x) {
      c(1, ringstarts[x]:(ringstarts[x+1]-1))
    })), nrow(df))
    return(df[indicies,])
  }
}
geom_outline <-
function(dat, xVar, yVar, groupVar, method="convex", convexHullMargin=NA, concaveHullAlpha=0.4, concaveHullExtend = T, ellipseLevel=0.90, alpha=0.3, linetype="solid")
{
	library(ggplot2)
	dat2 = dat[, c(xVar, yVar, groupVar)]
	dat2[, xVar] = as.numeric(dat2[, xVar])
	dat2[, yVar] = as.numeric(dat2[, yVar])
	dat2[, groupVar] = as.factor(dat2[, groupVar])
	if(method=="convex")
	{
		hulls = convexHull(dat2, xVar, yVar, groupVar, convexHullMargin)
		return(geom_polygon(data = hulls, alpha = alpha, linetype=linetype, aes_string(x="x", y="y", fill=groupVar, col=groupVar)))
	}else if(method == "ellipse")
	{
		
		return(stat_ellipse(data=dat2, aes_string(group=groupVar), level=ellipseLevel, geom = "polygon", alpha=alpha, linetype=linetype))
	}else if(method == "concave")
	{

		hulls = concaveHull(dat2, xVar, yVar, groupVar, concaveHullAlpha, convexHullMargin, extend = concaveHullExtend)
		return(geom_polygon(data = hulls, alpha = alpha, linetype=linetype, aes_string(x="x", y="y", fill=groupVar, col=groupVar)))
	}else if(method == "concaveSmooth")
	{

		hulls = concaveHullSmooth(dat2, xVar, yVar, groupVar, concaveHullAlpha, convexHullMargin, extend = concaveHullExtend)
		return(geom_polygon(data = hulls, alpha = alpha, linetype=linetype, aes_string(x="x", y="y", fill=groupVar, col=groupVar)))
	}
}
getGeneOrth <-
function()
{
    human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
    mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
    #genesOrth = getLDS(attributes = c("mgi_symbol"), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)
    genesOrth = getLDS(attributes = c("hgnc_symbol"), mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
    genesOrth2 = genesOrth[ genesOrth$MGI.symbol != "", ]
    return(genesOrth2)
}
getLegend <-
function(a.gplot)
{
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)	
}
getListIndex <-
function(x, n)
{
	if(n > length(x)) return("") else return(x[[n]])
}
getPathways <-
function(genesOrth = NULL)
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
getSqDist <-
function(p1, p2) 
{
    dx = p1$x - p2$x
    dy = p1$y - p2$y
    return(dx * dx + dy * dy)
}
getSqSegDist <-
function(p, p1, p2) 
{
    x = p1$x
    y = p1$y
    dx = p2$x - x
    dy = p2$y - y

    if (dx != 0 | dy != 0) 
	{
        t = ((p$x - x) * dx + (p$y - y) * dy) / (dx * dx + dy * dy)

        if (t > 1) {
            x = p2$x
            y = p2$y

        }else if (t > 0) {
            x = x + dx * t
            y = y + dy * t
        }
    }

    dx = p$x - x;
    dy = p$y - y;

    return(dx * dx + dy * dy)
}
getTSPOrder <-
function(dat)
{
	library(TSP)
	curDist = dist(dat)
	curDist = curDist / max(curDist)
	curTsp = TSP(curDist)
	curTsp = insert_dummy(curTsp, label = "cut")
	curTour = solve_TSP(curTsp, method="concorde", exe="/home/cem2009/.local/bin/concorde")
	curTour = cut_tour(curTour, "cut")
	curOrder = labels(curTour)
	return(curOrder)
}
getTSPOrderDist <-
function(curDist, concordePath="/home/cem2009/.local/bin/")
{
	library(TSP)
	curTsp = TSP(curDist)
	curTsp = insert_dummy(curTsp, label = "cut")
	#concorde_path(concordePath)
	#curTour = solve_TSP(curTsp, method="concorde")
	curTour = solve_TSP(curTsp)
	curTour = cut_tour(curTour, "cut")
	curOrder = labels(curTour)
	return(curOrder)
}
is.installed <-
function (mypkg) 
{
    is.element(mypkg, utils::installed.packages()[, 1])
}
iwanthue <-
function(n, hmin=0, hmax=360, cmin=39, cmax=80, lmin=35, lmax=80, 
                     plot=FALSE, random=FALSE) {
  library(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, 
            hmax <= 360, cmax <= 180, lmax <= 100, 
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), 
                                   seq(-100, 100, 5), 
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin & 
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    lab <- as(hcl, 'LAB')    
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  curColors = hex(LAB(clus$centers))
  return(unname(curColors))
  
}
loadRda <-
function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
	x <- env[[nm]]
    return(x)
}
maxN <-
function(x, N=2)
{
	x = x[!is.na(x)]
	len = length(x)
	if(N>len)
	{
		warning('N greater than length(x).  Setting N=length(x)')
		N = length(x)
	}
	return( sort(x, partial = len - N + 1)[len - N + 1] )
}
maxNidx <-
function(x, N=2)
{
	x = x[!is.na(x)]
	len = length(x)
	if(N>len)
	{
		warning('N greater than length(x).  Setting N=length(x)')
		N = length(x)
	}
	return( order(x, decreasing = T)[N] )
}
melt <-
function(...)
{
	return(reshape2::melt(...))
}
memUse <-
function()
{
	sort( sapply(ls(),function(x){object.size(get(x))}), decreasing = T)
}
merge.with.order <-
function(x,y, ..., sort = F, keep_order = 1)
{
	# this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
	add.id.column.to.data <- function(DATA)
	{
		data.frame(DATA, id... = seq_len(nrow(DATA)))
	}
	# add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
	order.by.id...and.remove.it <- function(DATA)
	{
		# gets in a data.frame with the "id..." column.  Orders by it and returns it
		if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")

		ss_r <- order(DATA$id...)
		ss_c <- colnames(DATA) != "id..."
		DATA[ss_r, ss_c]
	}

	# tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
	# tmp()

	if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
	if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
	
	warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
	return(merge(x=x,y=y,..., sort = sort))
}
multiplot <-
function(..., plotlist=NULL, file, cols=1, layout=NULL) 
{
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      #print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
      #                                layout.pos.col = matchidx$col))
	  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
normDataWithin <-
function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}
plotSymbols <-
function(dat, xName, yName, className, dotSize=3, alpha=0.5, do.legend = TRUE, do.letters = FALSE, label.size=4, colors=NULL, font.size=NULL, title = "", minCount = 25)
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
print.arrange <-
function(x) grid.draw(x)
print.facetAdjust <-
function(x, newpage = is.null(vp), vp = NULL) {
  if(newpage)
    grid.newpage()
	x2 = x
  if(is.null(vp)){
    grid.draw(x2)
  } else {
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(x2)
    upViewport()
  }
  invisible(x2)
}
push <-
function(listDat, element)
{
	listDat[[ length(listDat) + 1]] = element
	return(listDat)
}
read.tsv <-
function(file, header=T, ...)
{
	return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A", "NULL"), ...))
}
read.tsv2 <-
function(file, header=T, row.names=F, ...)
{
	if(row.names==F)
	{
		library(data.table)
		df = fread(file, stringsAsFactors=F, na.strings = c("NA","#N/A","N/A", "NULL"), data.table=F, ...)
		colnames(df) = make.names(colnames(df))
		return(df)
	}else
	{
		return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A", "NULL"), ...))
	}
}
rep.col <-
function(x,n)
{
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
rep.row <-
function(x,n)
{
   matrix(rep(x,each=n),nrow=n)
}
rmAll <-
function()
{
	rm(list = setdiff(ls(), lsf.str()))
}
rowMax2nd <-
function (x)
{
	apply(x, 1, function(row) max(row[-which.max(row)], na.rm=T))
}
rowMaxNth <-
function (x, N=2, returnIndex=FALSE)
{
	if(returnIndex)
	{
		return( apply(x, 1, function(row) maxNidx(row, N)) )
	}else
	{
		return( apply(x, 1, function(row) maxN(row, N)) )
	}
}
s_arrange <-
function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}
s_filter <-
function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}
s_group_by <-
function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}
s_mutate <-
function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}
s_select <-
function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}
s_summarise <-
function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
}
saveCurHistory <-
function (...) 
{
    if (interactive() & run & R.version$major >= 3) {
        try(utils::savehistory(Sys.getenv("R_HISTFILE")), silent = TRUE)
        return(TRUE)
    }
    else return(TRUE)
}
savePlot <-
function(outname, width = 11, height = 7, dpi = 150, pdfOutput = T, pngOutput = T)
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
scale_colour_Publication <-
function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
scale_fill_Publication <-
function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
setRemoveRownames <-
function(df, rownameColumn=1)
{
	rownames(df) = df[, rownameColumn]
	df2 = df[, -rownameColumn, drop=F]
	if(class(df2) != "data.frame")
	{
		df2 = data.frame(row.names=df[, rownameColumn], df[, -rownameColumn])
	}
	return(df2)
}
simplify <-
function(points, tolerance, highestQuality) 
{
    if (length(points) <= 2) return(points)

    if(!is.na(tolerance)) sqTolerance = tolerance*tolerance else sqTolerance = 1
	
	if(!highestQuality) points = simplifyRadialDist(points, sqTolerance)
    points = simplifyDouglasPeucker(points, sqTolerance)

    return(points)
}
simplifyDF <-
function(pointsDat, xVal, yVal, tolerance, highestQuality) 
{
	library(data.table)
    if (nrow(pointsDat) <= 2) return(pointsDat)
	
	points = list()
	for(i in 1:nrow(pointsDat)) { points[[i]] = data.frame(x=pointsDat[i, xVal], y=pointsDat[i, yVal])}

    pointsSimp = simplify(points, tolerance, highestQuality) 
	pointsSimp = rbindlist(pointsSimp)
    return(pointsSimp)
}
simplifyDouglasPeucker <-
function(points, sqTolerance) 
{
    last = length(points)

    simplified = list(points[[1]])
    simplified = simplifyDPStep(points, 1, last, sqTolerance, simplified)
    simplified = push(simplified, points[[last]])

    return(simplified)
}
simplifyDPStep <-
function(points, first, last, sqTolerance, simplified) 
{
    maxSqDist = sqTolerance
	if(first == last) return(simplified)
    for (i in (first + 1):last) 
	{
        sqDist = getSqSegDist(points[[i]], points[[first]], points[[last]])

        if (sqDist > maxSqDist) 
		{
            index = i
            maxSqDist = sqDist
        }
    }

    if (maxSqDist > sqTolerance) 
	{
        if (index - first > 1) simplified = simplifyDPStep(points, first, index, sqTolerance, simplified)
        simplified = push(simplified, points[[index]])
        if (last - index > 1) simplified = simplifyDPStep(points, index, last, sqTolerance, simplified)
    }
	return(simplified)
}
simplifyRadialDist <-
function(points, sqTolerance) 
{

    prevPoint = points[[1]]
    newPoints = c(prevPoint)
    
	len = length(points)
    for (i in 2:len) 
	{
        point = points[[i]];

        if (getSqDist(point, prevPoint) > sqTolerance) {
            newPoints = c(newPoints, point)
            prevPoint = point;
        }
    }

    if (prevPoint$x != point$x | prevPoint$y != point$y) newPoints = c(newPoints, point)

    return(newPoints)
}
smoothHull <-
function(df, X, Y)
{
	#theta = runif(n <- 300, 0, 2 * pi)
	#r = sqrt(runif(n, 0.25^2, 0.5^2))
	#x = cbind(0.5 + r * cos(theta), 0.5 + r * sin(theta))
	#df2 = as.data.frame(x)
	#colnames(df2) = c("X", "Y")

	n = nrow(df)
	df2 = cbind(1:(n+1), df[c(1:n, 1), X], df[c(1:n, 1), Y])
	t.coarse = seq(1, n+1, 0.05)
	x.coarse = approx(df2[, 1], df2[, 2], xout=t.coarse)$y
	y.coarse = approx(df2[, 1], df2[, 3], xout=t.coarse)$y
	## create a Bezier curve                                           
	library(Hmisc)
	bz = as.data.frame(bezier(x.coarse, y.coarse))
	return(bz)
}
smoothScore2d <-
function(score, x, y, type=NULL, numGrid = 100, knn = 100, m = 2, expand=0.05, xrng=NULL, yrng=NULL)
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
smoothScorePredict <-
function(trainDat, newGrid, knn=30)
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
splitGet <-
function (strings, sep, n)
{
	splitVals = strsplit(strings, sep)
	return ( sapply(splitVals, getListIndex, n) )
}
splitGetFromEnd <-
function (strings, sep, n)
{
	splitVals = strsplit(strings, sep)
	strLength = sapply(splitVals,length) - n + 1
	return(mapply("[[", splitVals, strLength))
}
subtree <-
function(x, i, ..., error_value, warn, simplify=TRUE){
    if(missing(error_value)) error_value <- NULL
    if(missing(warn)) warn <- is.null(error_value)
    if(is.null(error_value) && warn < 1){
        warning("With no `error_value` `warn` is ignored and all errors break the execution")
        warn <- 1
    }
    ret <- if(is.function(i)){
        if(missing(...)){
            i(x)
        } else {
            i(subtree(x, ..., error_value=error_value, warn=warn, simplify=simplify))
        }
    } else if(missing(...)){
        if(is.null(error_value)){
            x[i]
        } else {
            coerce_class <- function(x){
                x <- as(x, class(error_value))
                if(length(x) != length(error_value))
                    stop(sprintf("values must be length %i, but result is length %i",
                                 length(error_value), length(x)))
                x
            }
            lapply(x[i], function(xi){
                if(warn < 1) tryCatch({
                    coerce_class(xi)
                }, error = function(err){
                    if(warn == 0)
                        warning(err$message)
                    error_value
                }) else {
                    coerce_class(xi)
                }
            })
        }
    } else {
        lapply(x[i], subtree, ..., error_value=error_value, warn=warn, simplify=simplify)
    }
    if(simplify){
        if(length(ret) == 1){
            ret <- ret[[1]]
        } else if(all(sapply(ret, length) == 1)){
            ret <- unlist(ret, recursive=FALSE)
        } else if(all(sapply(lapply(ret, dim), is.null)) &&
                  all(sapply(ret, length) == length(ret[[1]]))){
            ret.class <- sapply(ret, class)
            ret.na <- sapply(ret, function(x) all(is.na(x)))
            i <- head(which(!ret.na), 1)
            if((is.numeric(ret[[i]]) || is.character(ret[[i]]) || is.logical(ret[[i]])) &&
                    length(unique(ret.class[!ret.na])) == 1){
                ret <- do.call(cbind, ret)
            }
        }
        if(is.null(dim(ret)) && !is.null(dim(x)) && length(ret) == length(x)){
            ret <- array(ret, dim=dim(x), dimnames=dimnames(x))
        }
    }
    ret
}
summarySE <-
function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
summarySEwithin <-
function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
swatch <-
function(x) {
  par(mai=c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.2, 0.4))
  barplot(rep(1, length(x)), col=rev(x), space = 0.1, axes=FALSE, 
          names.arg=rev(x), cex.names=0.8, horiz=T, las=1)
  return(invisible(NULL))
}
theme_Publication <-
function(base_size=14) {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90, vjust = 2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(size=12), 
			   axis.text.x = element_text(angle=-90, hjust=0, vjust=0.5),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}
ti_slingshot <-
function (cluster_method = "pam", ndim = 20L, shrink = 1L, reweight = TRUE, 
    reassign = TRUE, thresh = 0.001, maxit = 10L, stretch = 2L, 
    smoother = "smooth.spline", shrink.method = "cosine") 
{
    new_defaults <- as.list(environment())[formalArgs(param_overrider_fun)]
    param_names <- map_chr(definition$parameters$parameters, 
        "id")
    for (param_name in names(new_defaults)) {
        if (param_name %in% param_names) {
            definition$parameters$parameters[[which(param_name == 
                param_names)]]$default <- new_defaults[[param_name]]
        }
        else {
            stop("Unknown parameter: ", param_name)
        }
    }
    definition
}
toNum <-
function(x)
{
	return(as.numeric(as.character(x)))
}
wait <-
new("CFunc", .Data = function () 
.Primitive(".C")(<pointer: 0x0>), code = "#include <R.h>\n\n#include <sys/wait.h>\n\nextern \"C\" {\n  void file2d5207ddbda2e (  );\n}\n\nvoid file2d5207ddbda2e (  ) {\nint wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};\n}\n")
write.tsv <-
function(object, file, row.names=T, quote=F, col.names=NA)
{
	if(row.names==F & is.na(col.names)) col.names=T
	write.table(object, file, sep="\t", quote=quote, row.names=row.names, col.names=col.names)
}
