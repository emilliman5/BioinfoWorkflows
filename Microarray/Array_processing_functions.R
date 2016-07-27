##
## Array_processing_functions.R
##
## Some functions which may prove helpful in analyzing microarray data in R

plotSmoothScatter <- function
(x, y=NULL, bandwidthN=300,
 transformation=function(x)x^0.25,
 xlim=NULL, ylim=NULL, nbin=256,
 colramp=colorRampPalette(c("white", "lightblue", "blue", "orange", "orangered2")),
 doTest=FALSE,
 naAction=c("remove", "floor0", "floor1"),
 applyRangeCeiling=TRUE,
 ...)
{
   ## Purpose is to wrapper the smoothScatter() function in order to increase
   ## the apparent detail by adjusting the bandwidth parameters in a somewhat
   ## more automated/intuitive way than the default parameters.
   ##
   ## The smoothScatter() function uses some default number, which results
   ## in a low amount of detail, and in my opinion loses important features.
   ##
   ## To see an example of the difference, and hopefully the visual improvement,
   ## run this function plotSmoothScatter(doTest=TRUE)
   ##
   ## If provided with xlim and/or ylim values, and applyRangeCeiling=TRUE,
   ## this function will "cap" values outside the range so they fit on the
   ## edge of the visible plot.  The reason is to indicate visually whether
   ## there are many points outside the viewing range, while still allowing
   ## one to zoom to a reasonable numerical range for the majority of data points.
   ##
   ## naAction defines how NA values are handled.  By default, smoothScatter() chokes,
   ## but this function can either remove the offending points, or floor them at 1 or 0.
   ##
   naAction <- match.arg(naAction);
   if (is.null(colramp)) {
      colramp <- colorRampPalette(c("white", "lightblue", "blue", "orange", "orangered2"));
   }
   if (doTest) {
      testN <- 20000;
      x <- matrix(ncol=2, data=rnorm(testN*2));
      testN1 <- as.integer(testN*0.9);
      testN2 <- testN - testN1;
      x[1:testN1,2] <- x[1:testN1,1] + rnorm(testN1)*0.1;
      x[(testN1+1):(testN1+testN2),2] <- x[(testN1+1):(testN1+testN2),1] + rnorm(testN2)*0.9;
      oPar <- par();
      par("mfrow"=c(3,3));
      x[(testN1+1):(testN1+testN2),2] <- x[(testN1+1):(testN1+testN2),1] + rnorm(testN2)*0.3;
      smoothScatter(x, colramp=colramp,
                    main="smoothScatter\n10% noise at rnorm(x)*0.3\ndefault bandwidth",
                    xlim=c(-4,4), ylim=c(-4,4));
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.3\ncustom bandwidthN=", bandwidthN/2),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN/2);
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.3\ncustom bandwidthN=", bandwidthN),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN);
      x[(testN1+1):(testN1+testN2),2] <- x[(testN1+1):(testN1+testN2),1] + rnorm(testN2)*0.6;
      smoothScatter(x, colramp=colramp,
                    main="smoothScatter\n10% noise at rnorm(x)*0.6\ndefault bandwidth",
                    xlim=c(-4,4), ylim=c(-4,4));
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.6\ncustom bandwidthN=", bandwidthN/2),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN/2);
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.6\ncustom bandwidthN=", bandwidthN),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN);
      x[(testN1+1):(testN1+testN2),2] <- x[(testN1+1):(testN1+testN2),1] + rnorm(testN2)*0.9;
      smoothScatter(x, colramp=colramp,
                    main="smoothScatter\n10% noise at rnorm(x)*0.9\ndefault bandwidth",
                    xlim=c(-4,4), ylim=c(-4,4));
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.9\ncustom bandwidthN=", bandwidthN/2),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN/2);
      plotSmoothScatter(x, doTest=FALSE, colramp=colramp,
                        main=paste(sep="", "smoothScatter\n10% noise at rnorm(x)*0.9\ncustom bandwidthN=", bandwidthN),
                        xlim=c(-4,4), ylim=c(-4,4), bandwidthN=bandwidthN);
      par(oPar);
      return(NULL);
   }

   if (!is.null(y) && is.numeric(y) && is.numeric(x)) {
      x <- matrix(ncol=2, c(x, y));
   }
   if ((is.matrix(x) || is.data.frame(x)) && ncol(x) == 2) {
      y <- x[,2];
      x <- x[,1];
   } else if (is.null(Y) && (is.matrix(x) || is.data.frame(x))) {
      stop("Cannot handle matrix input x, when Y is NULL.");
   }
   ## Deal with NA values
   if (naAction == "remove") {
      naValues <- is.na(x) | is.na(y);
      x <- x[!naValues];
      y <- y[!naValues];
   } else if (naAction == "floor0") {
      naValuesX <- is.na(x);
      x[naValuesX] <- 0;
      naValuesY <- is.na(y);
      y[naValuesY] <- 0;
   } else if (naAction == "floor1") {
      naValuesX <- is.na(x);
      x[naValuesX] <- 1;
      naValuesY <- is.na(y);
      y[naValuesY] <- 1;
   }

   bandwidthN <- rep(bandwidthN, length.out=2);
   if (is.null(xlim)) {
      xlim <- range(x);
   }
   if (is.null(ylim)) {
      ylim <- range(y);
   }
   ## Apply a ceiling to values outside the range
   if (applyRangeCeiling) {
      tooHighX <- x > max(xlim);
      tooLowX <- x < min(xlim);
      x[tooHighX] <- max(xlim);
      x[tooLowX] <- min(xlim);
      tooHighY <- y > max(ylim);
      tooLowY <- y < min(ylim);
      y[tooHighY] <- max(ylim);
      y[tooLowY] <- min(ylim);
   }
   xlim4 <- (c(-1,1) * diff(xlim)*0.02) + xlim;
   ylim4 <- (c(-1,1) * diff(ylim)*0.02) + ylim;
   bandwidthXY <- c(diff(xlim4)/bandwidthN[1], diff(ylim4)/bandwidthN[2]);
   smoothScatter(x=x, y=y,
                 #xaxt="n", yaxt="n",
                 transformation=transformation,
                 bandwidth=bandwidthXY, nbin=nbin,
                 xlim=xlim4, ylim=ylim4, xaxs="i", yaxs="i",
                 colramp=colramp, ...);
   invisible(list(x=x, y=y,
                 transformation=transformation,
                 bandwidth=bandwidthXY, nbin=nbin,
                 xlim=xlim4, ylim=ylim4, xaxs="i", yaxs="i",
                 colramp=colramp));
}

centerGeneData <- function
(indata, floor=NA,
 controlSamples=NA, needsLog=NULL,
 mean=FALSE, returnGroupedValues=FALSE,
 ...)
{
   ## Purpose is to provide quick data-centering (subtracting average from each row, usually
   ## in log2 space).
   ##
   ## If mean=FALSE, then median is used for centering, sometimes better than
   ## using mean because it is less sensitive to outliers
   ##
   ## New: if data contains text and numeric columns, center the numeric columns only
   ##
   hasCharColumns <- FALSE;
   if (!class(indata) %in% c("matrix", "numeric")) {
      colClass <- sapply(1:ncol(indata), function(i){class(indata[,i])});
      colClassChar <- which(!colClass %in% c("numeric", "integer"));
      if (length(colClassChar) > 0) {
         hasCharColumns <- TRUE;
         indataChar <- indata[,colClassChar,drop=FALSE];
         indata <- as.matrix(indata[,-colClassChar]);
      }
   }
   controls <- colnames(indata) %in% controlSamples;#colnames(indata[,controlSamples]);

   if (length(controls[controls]) == 0) {
      controls <- rep(TRUE, dim(indata)[2]);
   }
   if (is.null(needsLog)) {
      ## We will try to auto-detect whether to log2 transform the data
      needsLog <- max(indata, na.rm=TRUE) > 40;
   }
   if (needsLog) {
      indata <- log2(indata);
   }
   ## Note: Switched to using sweep() because it is much, much faster than apply()
   if (mean) {
      indataMeans <- rowMeans(indata[,controls,drop=FALSE]);
      centeredData <- sweep(indata, 1, indataMeans);
      if (returnGroupedValues) {
         centeredData <- cbind(centeredData, "mean"=indataMeans);
      }
   } else {
      if (require(matrixStats, quietly=TRUE)) {
         indataMedians <- rowMedians(indata[,controls,drop=FALSE], na.rm=TRUE);
      } else {
         printDebug("Note: Install matrixStats package for markedly faster centerGeneData() operations.", fgText="yellow");
         indataMedians <- apply(indata[,controls,drop=FALSE], 1, function(x){
            median(x, na.rm=TRUE);
         });
      }
      centeredData <- sweep(indata, 1, indataMedians);
      if (returnGroupedValues) {
         centeredData <- cbind(centeredData, "median"=indataMedians);
      }
   }

   if (hasCharColumns) {
      printDebug("Keeping non-numeric columns as-is.");
      centeredData <- cbind(indataChar, centeredData)
   }
   return(centeredData);
}

plotSAsmooth <- function
(fit, xlab="Average log-expression", ylab="log2(sigma)",
 zero.weights=FALSE, pch=16, cex=0.2, useRaster=TRUE, ...)
{
   if (!is(fit, "MArrayLM"))
      stop("fit must be a MArrayLM object")
   x <- fit$Amean
   y <- log2(fit$sigma)
   if (!is.null(fit$weights) && !zero.weights) {
      w <- fit$weights
      w[is.na(w)] <- 0
      w[w < 0] <- 0
      allzero <- apply(w == 0, 1, all)
      y[allzero] <- NA
   }
   plotSmoothScatter(x, y, xlab=xlab, ylab=ylab, pch=pch, cex=cex,
      useRaster=useRaster, ...)
   lines(lowess(x, y, f=0.4), col="red")
   if (!is.null(fit$s2.prior)) {
      if (length(fit$s2.prior) == 1) {
         abline(h=log2(fit$s2.prior)/2, col="blue")
      } else {
         o <- order(x)
         lines(x[o], log2(fit$s2.prior[o])/2, col="blue")
         legend("topright", legend=c("lowess", "prior"),
         col=c("red", "blue"), lty=1)
      }
   }
   invisible()
}

eBayes2TopTables <- function
(lmFit3, defineHits=TRUE, cutoffAveExpr=1, cutoffAdjPVal=1, cutoffPVal=0.05, cutoffFold=1.5,
 returnFold=TRUE, mergeDF=FALSE, renameHeaders=TRUE, includeAveExpr=FALSE,
 ...)
{
   ## Purpose is to convert the lmFit3 results of eBayes() into a list of data.frames
   logratio2foldchange <- function(logratio, base=2) {
      retval <- base^(logratio)
      retval <- ifelse(retval < 1, -1/retval, retval);
      retval;
   }

   ## Pull out topTables for each coefficient (contrast)
   coefNames <- colnames(coef(lmFit3));
   #coefLabels <- nameVector(gsub("([0-9]+)", "\\1hr", coefNames), coefNames);
   coefLabels <- nameVector(coefNames, coefNames);
   lmTopTables <- lapply(nameVector(coefNames), function(i){
      iLabel <- coefLabels[i];
      iTopTable <- topTable(lmFit3, coef=i, sort.by="none", number=nrow(lmFit3));
      if (defineHits) {
         iTopTable[,"hit"] <- (iTopTable[,"AveExpr"] >= cutoffAveExpr & 
            iTopTable[,"adj.P.Val"] <= cutoffAdjPVal & 
            iTopTable[,"P.Value"] <= cutoffPVal & 
            abs(iTopTable[,"logFC"]) >= log2(cutoffFold))*
            sign(iTopTable[,"logFC"]);
      }
      if (returnFold) {
         iTopTable[,"fold"] <- logratio2foldchange(iTopTable[,"logFC"]);
      }
      if (renameHeaders) {
         colnames(iTopTable) <- paste(colnames(iTopTable), iLabel);
      }
      if (includeAveExpr) {
         iTopTableDF <- data.frame(check.names=FALSE, "Gene"=rownames(iTopTable),
             iTopTable[,provigrep(c("hit", "p.val", "logFC", "fold", "AveExpr"), colnames(iTopTable)),drop=FALSE]);
      } else {
         iTopTableDF <- data.frame(check.names=FALSE, "Gene"=rownames(iTopTable),
             iTopTable[,provigrep(c("hit", "p.val", "logFC", "fold"), colnames(iTopTable)),drop=FALSE]);
      }
      rownames(iTopTableDF) <- iTopTableDF[,"Gene"];
      iTopTableDF;
   });
   if (defineHits) {
      attr(lmTopTables, "cutoffAveExpr") <- cutoffAveExpr;
      attr(lmTopTables, "cutoffAdjPVal") <- cutoffAdjPVal;
      attr(lmTopTables, "cutoffPVal") <- cutoffPVal;
      attr(lmTopTables, "cutoffFold") <- cutoffFold;
   }
   if (mergeDF) {
      lmTopTablesAll <- mergeDFs(lmTopTables);
      lmTopTablesAll[,"AveExpr"] <- topTable(lmFit3, coef=1, sort.by="none", number=nrow(lmFit3))[,"AveExpr"];
      lmTopTablesAll <- lmTopTablesAll[,provigrep(c("Gene", "^hit", "."), colnames(lmTopTablesAll))];
      rownames(lmTopTablesAll) <- lmTopTablesAll[,"Gene"];
      if (defineHits) {
         attr(lmTopTablesAll, "cutoffAveExpr") <- cutoffAveExpr;
         attr(lmTopTablesAll, "cutoffAdjPVal") <- cutoffAdjPVal;
         attr(lmTopTablesAll, "cutoffPVal") <- cutoffPVal;
         attr(lmTopTablesAll, "cutoffFold") <- cutoffFold;
      }
      return(lmTopTablesAll);
   } else {
      return(lmTopTables);
   }
}

mergeDFs <- function
(inList)
{
   x <- inList[[1]];
   for (i in 2:length(inList)) {
      y <- inList[[i]];
      ## If we have merge.data.frame() arguments to pass,
      ## use do.call() so we can separately pass those function arguments
      x <- merge(x, y, all.x=TRUE, all.y=TRUE);
   }
   x;
}

provigrep <- function
(patterns, x, maxValues=NULL, sortFunc=c,
 rev=FALSE,
 ...)
{
   ## Purpose is to provide "progressive vigrep()" (which is value-returning, case-insensitive)
   ## mainly to allow for prioritized ordering of matching entries. For example,
   ## one might grep for a few patterns, but want certain pattern hits to come back first.
   ##
   ## rev will return the entries in reverse order, which is effective if you use a
   ## set of patterns to exclude. The entries at the end will have made it through the
   ## exclude patterns. However, entries within each subset are still sorted in the proper
   ## order.
   ##
   ## sortFunc is intended to allow for sorting each set of matched entries along the way,
   ## examples are sort and mirnaSort.
   valueSetL <- lapply(patterns, function(i){
      x <- grep(pattern=i, x=x, ignore.case=TRUE, value=TRUE, ...);
      if (!is.null(sortFunc)) {
         x <- sortFunc(x);
      }
      x;
   });
   if (rev) {
      valueSetL <- rev(valueSetL);
   }
   valueSet <- unique(unlist(valueSetL));
   if (!is.null(maxValues) && maxValues > 0 && length(valueSet) > 0) {
      #valueSet <- valueSet[1:maxValues];
      valueSet <- head(valueSet, maxValues);
   }
   return(valueSet);
}

getDate <- function()
{
   ## Purpose is to define a data in the format
   ## 2015may05 (YYYYmmDD)
   tolower(format(Sys.time(), "%Y%b%d"));
}

logratio2foldchange <- function (logratio, base=2)
{
    retval <- base^(logratio);
    retval <- ifelse(retval < 1, -1/retval, retval);
    retval;
}

######
##Eric's super cool functions
######

threeWayVenn<-function(x, filename="venn", categoryNames=c("Group 1", "Group 2", "Group 3"))
{
    venn.diagram(imagetype="png",cat.dist=c(0.1,0.1, 0.05) ,x = list(rownames(x[x[,1]!=0,]), 
                rownames(x[x[,2]!=0,]),rownames(x[x[,3]!=0,])), category.names = categoryNames,
                margin=c(0.2,0.2),filename=paste0("../results/",getDate(),"/",filename), 
                fill=c("firebrick4","darkblue", "darkorange"),alpha=c(0.5,0.5, 0.5))   
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

panel.density <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- density(x)
    polygon(h, col="skyblue", border="skyblue4", ...)
}


