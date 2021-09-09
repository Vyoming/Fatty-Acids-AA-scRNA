anyNA = function(x) { any(is.na(x)) }

fixfeature <- function(df) {
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

custom_fortify <- function(x, ...) {
  df <- fortify(x, ...)
  df %>% group_by(id) %>% do(fixfeature(.))
}

findConvexHull = function(df, X, Y, margin)
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

##


findConcaveHull = function(df, X, Y, alpha, hullMargin=NA, extend=F, minPoint=3)
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

extendHull = function(df, X, Y, margin=10)
{
	library(polyclip)
	extended = polyoffset(list(x=df[, X], y=df[, Y]), margin, jointype="round", arctol=abs(margin)/40)
	extended2 = data.frame(x=extended[[1]]$x, y=extended[[1]]$y)
	return(extended2)
}

smoothHull = function(df, X, Y)
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

convexHull = function(df, X, Y, group, margin=NA, marginMultiplier = 3)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConvexHull, X, Y, margin)
	hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	# hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls2)
}

concaveHull = function(df, X, Y, group, alpha, margin=NA, extend=F, minPoint = 3, marginMultiplier = 1)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConcaveHull, X, Y, alpha, margin, extend, minPoint)
	#hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	# hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls)
}

concaveHullSmooth = function(df, X, Y, group, alpha, margin=NA, extend=F, minPoint = 3, marginMultiplier = 1)
{
	if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
	hulls = ddply(df, group, findConcaveHull, X, Y, alpha, margin, extend, minPoint)
	#hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
	hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
	
	return(hulls2)
}

geom_outline = function(dat, xVar, yVar, groupVar, method="convex", convexHullMargin=NA, concaveHullAlpha=0.4, concaveHullExtend = T, ellipseLevel=0.90, alpha=0.3, linetype="solid")
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


###

push = function(listDat, element)
{
	listDat[[ length(listDat) + 1]] = element
	return(listDat)
}

getSqDist = function(p1, p2) 
{
    dx = p1$x - p2$x
    dy = p1$y - p2$y
    return(dx * dx + dy * dy)
}

# square distance from a point to a segment
getSqSegDist = function(p, p1, p2) 
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

#basic distance-based simplification
simplifyRadialDist = function(points, sqTolerance) 
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

simplifyDPStep = function(points, first, last, sqTolerance, simplified) 
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

# simplification using Ramer-Douglas-Peucker algorithm
simplifyDouglasPeucker = function(points, sqTolerance) 
{
    last = length(points)

    simplified = list(points[[1]])
    simplified = simplifyDPStep(points, 1, last, sqTolerance, simplified)
    simplified = push(simplified, points[[last]])

    return(simplified)
}

# both algorithms combined
simplify = function(points, tolerance, highestQuality) 
{
    if (length(points) <= 2) return(points)

    if(!is.na(tolerance)) sqTolerance = tolerance*tolerance else sqTolerance = 1
	
	if(!highestQuality) points = simplifyRadialDist(points, sqTolerance)
    points = simplifyDouglasPeucker(points, sqTolerance)

    return(points)
}

simplifyDF = function(pointsDat, xVal, yVal, tolerance, highestQuality) 
{
	library(data.table)
    if (nrow(pointsDat) <= 2) return(pointsDat)
	
	points = list()
	for(i in 1:nrow(pointsDat)) { points[[i]] = data.frame(x=pointsDat[i, xVal], y=pointsDat[i, yVal])}

    pointsSimp = simplify(points, tolerance, highestQuality) 
	pointsSimp = rbindlist(pointsSimp)
    return(pointsSimp)
}

