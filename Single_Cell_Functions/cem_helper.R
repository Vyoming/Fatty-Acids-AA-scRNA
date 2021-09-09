library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(matrixStats)


options(stringsAsFactors=FALSE)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_cem = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())
 
cemPalette = function(x) 	{ if(x<=12) brewer.pal(x, "Paired") else (c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Pastel1")))[1:x]	}


theme_Publication <- function(base_size=14) {
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

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


cp = function(dat) { write.tsv(dat, "~/Rcopy.txt") }

loadRda = function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
	x <- env[[nm]]
    return(x)
}

rmAll = function()
{
	rm(list = setdiff(ls(), lsf.str()))
}
anyNA = function(x) { any(is.na(x)) }

rep.row = function(x,n)
{
   matrix(rep(x,each=n),nrow=n)
}
rep.col = function(x,n)
{
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

toNum = function(x)
{
	return(as.numeric(as.character(x)))
}

maxN = function(x, N=2)
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

maxNidx = function(x, N=2)
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


rowMax2nd = function (x)
{
	apply(x, 1, function(row) max(row[-which.max(row)], na.rm=T))
}

rowMaxNth = function (x, N=2, returnIndex=FALSE)
{
	if(returnIndex)
	{
		return( apply(x, 1, function(row) maxNidx(row, N)) )
	}else
	{
		return( apply(x, 1, function(row) maxN(row, N)) )
	}
}

colMax2nd = function (x)
{
	apply(x, 2, function(row) max(row[-which.max(row)], na.rm=T))
}

colMaxNth = function (x, N=2, returnIndex=FALSE)
{
	if(returnIndex)
	{
		return( apply(x, 2, function(row) maxNidx(row, N)) )
	}else
	{
		return( apply(x, 2, function(row) maxN(row, N)) )
	}
}

memUse = function()
{
	sort( sapply(ls(name=.GlobalEnv),function(x){object.size(get(x))}), decreasing = T) 
}

addFootnote = function(ggplotPlot, footnoteText, fontsize=15, fontface="italic", plot=T)
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

print.arrange <- function(x) grid.draw(x)

`?` <- function(x, y) {
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

':=' <- function(lhs, rhs) {
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

merge.with.order <- function(x,y, ..., sort = F, keep_order = 1)
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


#rowSds <- function(x, center=NULL, ...) {
#	n <- !is.na(x);
#	n <- rowSums(n);
#	n[n <= 1] <- NA;
#	
#	
#	if (is.null(center)) {
#	center <- rowMeans(x, ...);
#	}
#	
#	x <- x - center;
#	x <- x*x;
#	x <- rowSums(x, ...);
#	x <- x/(n-1);
#	sqrt(x);
#}


read.tsv = function(file, header=T, ...)
{
	return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A", "NULL"), ...))
}

read.tsv2 = function(file, header=T, row.names=F, ...)
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

write.tsv = function(object, file, row.names=T, quote=F, col.names=NA)
{
	if(row.names==F & is.na(col.names)) col.names=T
	write.table(object, file, sep="\t", quote=quote, row.names=row.names, col.names=col.names)
}

#setRemoveRownames = function(df, rownameColumn=1)
#{
#	rownames(df) = df[, rownameColumn]
#	df = df[, -rownameColumn, drop=F]
#	return(df)
#}

setRemoveRownames = function(df, rownameColumn=1)
{
	rownames(df) = df[, rownameColumn]
	df2 = df[, -rownameColumn, drop=F]
	if(class(df2) != "data.frame")
	{
		df2 = data.frame(row.names=df[, rownameColumn], df[, -rownameColumn])
	}
	return(df2)
}

addRownameColumn = function(df, idColName="id")
{
	ids = rownames(df)
	rownames(df) = NULL
	#df = cbind(ids, df)
	df = data.frame(ids, df)
	colnames(df)[1] = idColName
	return(df)
}



#
#splitGet=function (strings, sep, n)
#{
#	return ( sapply(strsplit(strings, sep), "[[", n) )
#}

getListIndex = function(x, n)
{
	if(n > length(x)) return("") else return(x[[n]])
}
splitGet=function (strings, sep, n)
{
	splitVals = strsplit(strings, sep)
	return ( sapply(splitVals, getListIndex, n) )
}

splitGetFromEnd=function (strings, sep, n)
{
	splitVals = strsplit(strings, sep)
	strLength = sapply(splitVals,length) - n + 1
	return(mapply("[[", splitVals, strLength))
}
		
GetColRowLayout = function(n)
{
	if(n == 3) return(c(3, 1))
	if(n == 7) return(c(4, 2))
	if(n == 8) return(c(4, 2))
	if(n == 10) return(c(5, 2))
	if(n == 14) return(c(5, 3))
	if(n == 15) return(c(5, 3))
	if(n == 17) return(c(6, 3))
	if(n == 18) return(c(6, 3))
	
	return( n2mfrow(n))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
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

getLegend<-function(a.gplot)
{
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)	
}

detachAllPackages <- function() 
{
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
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
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
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

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
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



#' Extract a subset of a tree of nested lists
#' 
#' Modeling results produced by \code{\link{evaluate}} comes in the
#' form of nested lists. This function can be used to subset or rearrange parts
#' of the results into vectors, matrices or data frames.
#' Also note the \code{\link[emil]{select}} function that provides an extension
#' to the \pkg{dplyr} package for data manipulation.
#' 
#' This function can only be used to extract data, not to assign.
#' 
#' @param x List of lists.
#' @param i Indexes to extract on the first level of the tree. Can also be a
#'   function that will be applied to the downstream result of the function.
#' @param ... Indexes to extract on subsequent levels.
#' @param error_value A template for the return value in case it is missing or
#'   invalid. Note that \code{NA} is a \code{\link{logical}} by default,
#'   causing \code{subtree} to also convert existing results to logicals.
#'   To get around this, please specify it as \code{as.numeric(NA)},
#'   \code{as.character(NA)}, or similar (see the example below).
#' @param warn Specifies whether warnings should be displayed (\code{0}),
#'   ignored (\code{-1}), or break execution (\code{1}). Works like the
#'   \code{\link{options}} parameter \code{warn}.
#' @param simplify Whether to collapse results into vectors or matrices when
#'   possible (\code{TRUE}) or to preserve the original tree structure as a
#'   list (\code{FALSE}).
#' @return A subset of the list tree.
#' @example examples/subtree.r
#' @seealso \code{\link{select}}, \code{\link{get_prediction}},
#'   \code{\link{get_importance}}, \code{\link{get_tuning}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
subtree <- function(x, i, ..., error_value, warn, simplify=TRUE){
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


# # pos - where to add new labels
# # newpage, vp - see ?print.ggplot
# facetAdjust2 <- function(x, pos = c("up", "down"), 
                        # newpage = is.null(vp), vp = NULL)
# {
	# library(grid)
	# library(ggplot2)
  # # part of print.ggplot
  # ggplot2:::set_last_plot(x)
  # if(newpage)
    # grid.newpage()
  # pos <- match.arg(pos)
  # p <- ggplot_build(x)
  # gtable <- ggplot_gtable(p)
  # # finding dimensions
  # dims <- apply(p$panel$layout[2:3], 2, max)
  # nrow <- dims[1]
  # ncol <- dims[2]
  # # number of panels in the plot
  # panels <- sum(grepl("panel", names(gtable$grobs)))
  # space <- ncol * nrow
  # # missing panels
  # n <- space - panels
  # # checking whether modifications are needed
  # if(panels != space){
    # # indices of panels to fix
    # idx <- (space - ncol - n + 1):(space - ncol)
    # # copying x-axis of the last existing panel to the chosen panels 
    # # in the row above
    # gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    # if(pos == "down"){
      # # if pos == down then shifting labels down to the same level as 
      # # the x-axis of last panel
      # rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   # gtable$layout$name)
      # lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      # gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    # }
  # }
  # # again part of print.ggplot, plotting adjusted version
  # if(is.null(vp)){
    # grid.draw(gtable)
  # }
  # else{
    # if (is.character(vp)) 
      # seekViewport(vp)
    # else pushViewport(vp)
    # grid.draw(gtable)
    # upViewport()
  # }
  # invisible(p)
# }

facetAdjust <- function(x, pos = c("up", "down"))
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

print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
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

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

# Helper functions that allow string arguments for  dplyr's data modification functions like arrange, select etc. 
# Author: Sebastian Kranz

# Examples are below

#' Modified version of dplyr's filter that uses string arguments
#' @export
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}

#' Modified version of dplyr's select that uses string arguments
#' @export
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_arrange = function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_mutate = function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}

#' Modified version of dplyr's summarise that uses string arguments
#' @export
s_summarise = function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @export
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}

#' Internal function used by s_filter, s_select etc.
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df  
}

#' Plot colour swatches for a vector of colours
#' 
#' Plot named colour swatches for a vector of colours.
#' 
#' @param x a vector of colours, specified as: colour names (i.e.
#' colour names returned by \code{colors()}); numeric indices into 
#' \code{palette()}, or hexadecimal strings in the form \code{"#RRGGBB"}, where 
#' \code{RR}, \code{GG}, and \code{BB} are pairs of hexadecimal digits 
#' representing red, green, and blue components, in the range \code{00} to 
#' \code{FF}.
#' @return \code{NULL}. The colour swatch is plotted to the active plotting 
#'   device.
#' @seealso \code{\link{iwanthue}}
#' @export
#' @examples 
#' swatch(colours()[1:10])
#' swatch(1:4)
#' swatch(iwanthue(5))
swatch <- function(x) {
  curNames = rev(x)
  if(!is.null(names(x))) curNames = paste0(rev(names(x)), " - ", curNames)
  par(mai=c(0.4, max(strwidth(x, "inch") + 2, na.rm = TRUE), 0.2, 0.4))
  barplot(rep(1, length(x)), col=rev(x), space = 0.1, axes=FALSE, 
          names.arg=curNames, cex.names=0.8, horiz=T, las=1)
  return(invisible(NULL))
}

#' Generate a colour palette by k-means clustering of LAB colour space.
#' 
#' Generate a palette of distinct colours through k-means clustering of LAB 
#' colour space.
#' 
#' @param n Numeric. The number of colours to generate.
#' @param hmin Numeric, in the range [0, 360]. The lower limit of the hue range 
#'   to be clustered.
#' @param hmax Numeric, in the range [0, 360]. The upper limit of the hue range 
#'   to be clustered.
#' @param cmin Numeric, in the range [0, 180]. The lower limit of the chroma 
#'   range to be clustered.
#' @param cmax Numeric, in the range [0, 180]. The upper limit of the chroma 
#'   range to be clustered.
#' @param lmin Numeric, in the range [0, 100]. The lower limit of the luminance 
#'   range to be clustered.
#' @param lmax Numeric, in the range [0, 100]. The upper limit of the luminance 
#'   range to be clustered.
#' @param plot Logical. Should the colour swatches be plotted (using 
#'   \code{\link{swatch}})?
#' @param random Logical. If \code{TRUE}, clustering will be determined by the 
#'   existing RNG state. If \code{FALSE}, the seed will be set to \code{1} for 
#'   clustering, and on exit, the function will restore the pre-existing RNG 
#'   state.
#' @return A vector of \code{n} colours (as hexadecimal strings), representing 
#'   centers of clusters determined through k-means clustering of the LAB colour
#'   space delimited by \code{hmin}, \code{hmax}, \code{cmin}, \code{cmax}, 
#'   \code{lmin} and \code{lmax}.
#' @details Note that \code{iwanthue} currently doesn't support \code{hmin} 
#'   greater than \code{hmax} (which should be allowed, since hue is circular).
#' @references 
#' \itemize{
#'   \item \href{http://tools.medialab.sciences-po.fr/iwanthue/}{iwanthue - colors for data scientists}
#'   \item \href{https://github.com/medialab/iwanthue}{iwanthue on
#'   GitHub}   
#' }
#' @seealso \code{\link{swatch}}
#' @export
#' @importFrom colorspace LAB hex coords
#' @examples 
#' iwanthue(5)
#' iwanthue(5, plot=TRUE)
#' iwanthue(5, 0, 240, 0, 24, 0, 100, plot=TRUE) # shades
#' iwanthue(5, 0, 360, 0, 54, 67, 100, plot=TRUE) # pastel
#' iwanthue(5, 0, 360, 54, 180, 27, 67, plot=TRUE) # pimp
#' iwanthue(5, 0, 360, 36, 180, 13, 73, plot=TRUE) #intense
#' iwanthue(3, 0, 300, 60, 180, 73, 100, plot=TRUE) # fluoro
#' iwanthue(3, 220, 260, 12, 150, 0, 53, plot=TRUE) # blue ocean
#iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100, 
iwanthue <- function(n, hmin=0, hmax=360, cmin=39, cmax=80, lmin=35, lmax=80, 
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

dcast = function(...)
{
	return(reshape2::dcast(...))
}

melt = function(...)
{
	return(reshape2::melt(...))
}

getTSPOrder = function(dat, concordePath="/home/cem2009/.local/bin/")
{
	library(TSP)
	curDist = dist(dat)
	curDist = curDist / max(curDist)
	curTsp = TSP(curDist)
	curTsp = insert_dummy(curTsp, label = "cut")
	#concorde_path(concordePath)
	#curTour = solve_TSP(curTsp, method="concorde")
	curTour = solve_TSP(curTsp)
	curTour = cut_tour(curTour, "cut")
	curOrder = labels(curTour)
	return(curOrder)
}

getTSPOrderDist = function(curDist, concordePath="/home/cem2009/.local/bin/")
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

median.quartile = function(x)
{
	out = quantile(x, probs = c(0.25, 0.5, 0.75))
	names(out) = c("ymin", "y", "ymax")
	return(out) 
}


library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')


GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
        )
}