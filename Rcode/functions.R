## readTables
##
##    This function reads eight tables, one per each file of serbgo database, merge them
##    and generates an output list with two components. First component is a data.frame
##    with the merged table, and second component is is a table with the id and the associated
##    name when year is not 2005, other wise is NULL
##    a table with the 
##
## Argument(s)
##
##     path :    
##   labels : 
##     year : 


readTables <- function(path, labels, year)
{
    filenames <- file.path(datDir, paste0("20", year), paste0(labels, year, ".txt"))

    main <- read.table(file = filenames[1], header = TRUE, sep = "\t")
    type <- read.table(file = filenames[2], header = TRUE, sep = "\t")
    spec <- read.table(file = filenames[3], header = TRUE, sep = "\t")
    data <- read.table(file = filenames[4], header = TRUE, sep = "\t")
    anot <- read.table(file = filenames[5], header = TRUE, sep = "\t")
    stat <- read.table(file = filenames[6], header = TRUE, sep = "\t")
    outp <- read.table(file = filenames[7], header = TRUE, sep = "\t")
    tags <- read.table(file = filenames[8], header = TRUE, sep = "\t")

    if(year=="05")
    {
        out.toolnames <- main[c("mgoid", "nameT")]
        out.toolnames <- as.matrix(out.toolnames)
        out.toolnames[10, 2] <- "GOstat"
        out.toolnames <- as.data.frame(out.toolnames)
    }else{
        out.toolnames <- NULL
    }
    
    fields <- c("nameT", "promoter", "link", "refer")

    main.mat <- as.matrix(main[-match(fields, names(main))])
    type.mat <- as.matrix(type)
    spec.mat <- as.matrix(spec)
    data.mat <- as.matrix(data)
    anot.mat <- as.matrix(anot)
    stat.mat <- as.matrix(stat)
    outp.mat <- as.matrix(outp)

    out.df <- merge(main.mat, type.mat, by.x = "mgoid", by.y = "tgoid")
    out.df <- merge(out.df, spec.mat, by.x = "mgoid", by.y = "sgoid")
    out.df <- merge(out.df, data.mat, by.x = "mgoid", by.y = "dgoid")
    out.df <- merge(out.df, anot.mat, by.x = "mgoid", by.y = "agoid")
    out.df <- merge(out.df, stat.mat, by.x = "mgoid", by.y = "xgoid")
    out.df <- merge(out.df, outp.mat, by.x = "mgoid", by.y = "ogoid")

    out <- list(capabilities = out.df, toolnames = out.toolnames)

    return(out)
}


## funcBySection
##
##     This functions provides a list of matrices, one per section of SerbGO functionalities.
##
## Argument(s)
##
##     sec : list of the sections where each slot consist of the names of the fields in SerbGO
##     mat : matrix with the capabilities associated with each GOTool


funcBySection <- function(sec, mat)
{
    out <- lapply(sec, function(x, mat)
                       {
                           new.mat <- mat[, x]
                           return(new.mat)
                       },
                       mat = mat)

    return(out)
}


## absFreqBySection
##
## Argument(s)
##
## mat.sec

absFreqBySection <- function(mat.sec)
{

    out <- lapply(mat.sec, function(mat)
                           {
                               abs.freq <- apply(mat, 1, sum)
                               return(abs.freq)
                           })
    return(out)
}


## absFreq.df
##
##   From a list of lists where each slot is a section of functionalities of  serbgo describing GO Tools,
##   generates a data.frame where first two columns are "section", "GO tool", and the following columns
##   are the absolute frequencies per each year reviewed and summarized.
##
## Argument(s)
##
##   lst : list of lists of sections of serbgo describing capabilities per each tool

absFreq.df <- function(lst)
{
   gotool <- names(lst[[1]][[1]])
   year <- names(lst)
   
   gotool.sec <- lapply(lst, as.data.frame)
   out.ls <- lapply(gotool.sec, melt)
   out.df <- as.data.frame(out.ls)
   out.df <- out.df[, -agrep(".variable", colnames(out.df))[-1]]
   out <- data.frame(out.df[, 1], rep(gotool, length(year)), out.df[, -1])
   colnames(out) <- c("Section", "GOTool", year)
                     
   return(out)
}

## relFreq.df
##
##   Given a data.frame with absolute frequencies per GO tool nd section and the number of available
##   functionalites, calulates relative frequencies associated with each position of this data.frame
## Argument(s)
##
##   dtf :
##   lst :

relFreq.df <- function(dtf, lst)
{

    dtf.ls <- split(dtf[,-c(1,2)], f = dtf$Section)
    out.ls <- lapply(names(lst), function(i, x, y)
                                 {
                                     rel.freq <- (x[[i]] / y[[i]])  * 100
                                     return(rel.freq)
                                 },
                                 x = dtf.ls,
                                 y = lst)
    out.df <- do.call("rbind", out.ls)
    colnames(out.df) <- gsub("Abs", "Rel", colnames(out.df))
    out <- cbind(dtf[,1:2],out.df)
    
    return(out)
}



## relFreqBySection
##
##
## Argument(s)
##
##   abs.freq :
##      n.sec :

relFreqBySection <- function(abs.freq, n.sec)
{
    sec.names <- names(n.sec)
    
    out <- lapply(sec.names, function(x, abs.freq, n.sec)
                             {
                                 n <- n.sec[[x]]
                                 rel.freq <- (abs.freq[[x]] / n) * 100
                                 return(rel.freq)
                             },
                             abs.freq = abs.freq,
                             n.sec = n.sec)
    names(out) <- sec.names
    
    return(out)
}


## splitBySection
##
##   Splits a data.frame by a factor such that it is one of the first two columns of the data.frame
##
## Argument(s)
##
##   df : data.frame where two first columns are factors
##    j : columns number of the data.frame indicating the factor used to split the data.frame

splitBySection <- function(df, j)
{
    f.split <- df[, j]
    f.names <-unique(as.character( df[, ifelse(j==1, 2, 1)]))
        
    out.ls <- split(df[, -c(1,2)], f = f.split)
    out <- lapply(out.ls, function(x, i)
                          {
                            mat <- as.matrix(x)
                            rownames(mat) <- i
                            return(mat)
                          },
                          i = f.names)
           
    return(out)
}


## barplot.2
##
##    Creates a bar plot with vertical for each matrix, associated with a section of SerbGO capabilities, and
##    corresponding to the number of functionalities of GO tools (rows) at three different years (columns).
##    Frequency matrices are provided in slots of a list.
##
## Argument(s)
##
##        lst : list of matrices containing absolute frequencies of GO at three different years
##        sec : character vector with the names associated with each section, that are going to be used as a
##              titles of the bar plots.
##     outDir : character vector with the name of the output path. NULL by default.
##       file : character vector with the name of the output file. NULL by default
##  save.file : if TRUE, then saves each plot in a .png file. TRUE by default

barplot.2 <- function(lst, sec, outDir = NULL, file = NULL, save.file = TRUE)
{
    ptm <- proc.time()
    
    i <- names(lst)
    names(sec) <- i
    
    lapply(i, function(i, lst, sec, outDir, file, save.file)
              {
                  if(!(save.file))
                  {
                      dev.new()
                      out.plot <- paste("Barplot of", sec[i])
                         
                  }else{
                      if(is.null(outDir))
                      {
                          pdf(file = paste0(file, ".pdf"))
                          out.plot <- paste("Barplot of", sec[i], "saved in", paste0(file, ".pdf"))
                      }else{
                          pdf(file = file.path(outDir, paste(file, sec[i], "pdf", sep = ".")))
                          out.plot <- paste("Barplot of", sec[i],
                                            "saved in", file.path(outDir, paste(file, sec[i], "pdf", sep = ".")))
                      }
                  }

                  barplot(as.matrix(t(lst[[i]])), names.arg = rownames(lst[[i]]), horiz = FALSE , beside = TRUE, 
                          las = 3,
                          col = c("tomato", "lightyellow",  "royalblue2"),
                          ylab = "% of functionalities",
                          main = sec[i],
                          cex.names = 0.7,                                
                          off = 0)
                  legend("topleft",
                         legend = c("2005", "2007", "2009"),
                         fill = c("tomato", "lightyellow",  "royalblue2"))
                  
                  if(save.file) dev.off()
                  
                  return(out.plot)
              },
              lst = lst,
              sec = sec,
              outDir = outDir,
              file = file,
              save.file)

    out <- proc.time() - ptm

    return(out)
}




## matching.coeff
##
##
##
## Argument(s)
##
##   x :

matching.coeff <- function(x)
{

  n <- dim(x)[1]
  m <- dim(x)[2]

  dist.mat <- matrix(0,n,n)
  colnames(dist.mat) <- rownames(x)
  rownames(dist.mat) <- rownames(x)

  for(i in 1:n)
  {
    g1 <- x[i,]
    for(j in i:n)
    {
      g2 <- x[j,]
      partial.sum <- apply(rbind(g1, g2), 2, sum)
      a <- table(partial.sum)["2"]
      if(is.na(a)) a <- 0
      d <- table(partial.sum)["0"]
      if(is.na(d)) d <- 0
      num <- a + d
      dist.mat[i,j] <- 1 - num / m
      dist.mat[j,i] <- dist.mat[i,j]
    }
  }
  return(as.dist(dist.mat))
}

## avgSilVal
##
##   Given a distance object, yields a matrix with the numbers of clusters and the associated average silhouette
##   width for each cluster
##
## Argument(s)
##
##   d : distance object

avgSilVal <- function(d)
{
    n <- nrow(as.matrix(d))
    out <- cbind(1:n, numeric(n))
              
    for (k in 2:(n-1))
    {
        out[k, 2] <- pam(d, k)$silinfo$avg.width
    }
                 
    colnames(out) <- c("Num.Clusters", "Avg.Silh.Width")

    return(out)
}

## plotAvgSilVal
##
##    Plots the average silhouette width for each number of clusters. 
##
## Argument(s)
##
##       d : distance object
##    main : character with the title of the plot
##    xlab : character with the label of X axis
##    ylab : character with the label of Y axis

plotAvgSilVal <- function(d, main = "Average Silhouette Values per Number of Clusters",
                          xlab = "Num. of clusters", ylab = "Avg. Silhouette Value")
{
    ptm <- proc.time()

    asw <- avgSilVal(d)
    k <- nrow(asw)
    
    plot(asw[-1, ],
         main = main,
         sub = paste("\nOptimal Num. of Clusters: ", which.max(asw[,2]),
                     ", Silhouette: ", round(max(asw[, 2]), 3)),
         xlab = xlab, ylab = ylab, lwd = 3,
         type = "h", pch = 19, col = "blue", cex.axis = 0.6, cex.sub = 0.6, cex.lab = 0.6) #, xaxt = "n")
    #  axis(side = 1, at = 0:k, labels = 0:k, cex.axis = 0.7)
    
    axis(side = 1, at = which.max(asw[, 2]),
         paste("optimum", which.max(asw[, 2]), sep = "\n"),
         col = "red", font = 1.5, col.axis = "red")
    points(which.max(asw[, 2]), max(asw[, 2]), pch = 19, col = "red", cex = 1)

    out <- proc.time() - ptm

    return(out)
}

## plotSil
##
##   Given a distance object, creates a Silhouette plot of Partitioning Around Medoids (PAM) with the
##   optimal number of clusters
##
## Argument(s)
##
##           d : distance object
##        main : character with the title of the plot
##   si.colors : character vector with the colors associated with each cluster

plotSil <- function(d, main = NULL, si.colors = NULL, cex = 0.6, do.n.k = FALSE)
{
    ptm <- proc.time()

    asw <- avgSilVal(d)
    
    if(is.null(si.colors))
    {
        k <- which.max(asw[, 2])
        si.colors <- rainbow(k)
    }else{
        k <- length(si.colors)
    }
    
    pam.k <- pam(d, k)
    si <- silhouette(pam.k)
    opt <- par(cex = cex)
    plot(si, main = main, col = si.colors, do.n.k = do.n.k)
    par(opt)

    out <- proc.time() - ptm

    return(out)
}



## leafColors
##
##    Given a dendrogram, looks for leaves, and if so gets the color among colors provided, which are
##    associated with clusters.
##
## Argument(s)
##
##              n : node from the dendrogram
##    labelColors : character vector with the colors associated with each cluster
##     clusMember : numeric vector with the number of the cluster associated with each objects of the
##                  dendrogram
##         branch : TRUE if branches must be colored. FALSE by default

leafColors <- function(n, labelColors, clusMember, branch = FALSE)
{
    if (is.leaf(n))
    {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        if(branch)
        {
            attr(n, "edgePar") <- c(a$nodePar, list(col = labCol))
        }
    }
   
    return(n)
}


## plotDendro
##
##    Plots amore customized hierarchical cluster base on the function 'dendrapply' that can
##    be used to apply a function to all nodes of a dendrgoram. This comes very handy if one wants
##    to add some color to the labels.
##
## Argument(s)
##
##             hc : hclust object
##           main : character with the title of the plot
##           ylab : character with the label of Y axis
##           type : type of edges of the dendrogram. By default takes "rectangle", but it can take "triangle"
##                  shape
##    labelColors : character vector with the colors associated with each cluster
##     clusMember : numeric vector with the number of the cluster associated with each objects of the
##                  dendrogram
##         branch : TRUE if branches must be colored. FALSE by default
##            cex :

plotDendro <- function(hc, main = NULL, ylab = "Distance",
                       type = "rectangle", labelColors, clusMember, branch = FALSE, cex = 0.8)                  
{
    ptm <- proc.time()
    
    hcd <- as.dendrogram(hc)

    clusDendro <- dendrapply(hcd, leafColors, labelColors, clusMember, branch)
    op <- par(mar = par("mar") + c(4, 0, 0, 0), cex = cex)
    plot(clusDendro, main = main, type = type, ylab = ylab)
    par(op)

    out <- proc.time() - ptm

    return(out)
}

## adequacy
##
##    Calculates the agreement measure of P_{m}^2 proposed by Mardia et al
##
## Argument(s)
##
##    v : numeric vector with the eigenvalues of a classical MDS, yielded with the function 'cmdscale' from


adequacy <- function(v)
{
  adq <- numeric()
  n <- which(v < 0)[1]
  if (is.na(n)) n <- length(v)
  for(i in 1:n)
  {
    adq[i] <- sum((v[1:i])^2) / sum(v^2)
  }
  return(adq)
}


## adqCMDS
##
##   Given a list of list resulting classical multidimensinal scaling results performed with 'cmdscale'
##   function from package 'stats', yields a matrix where for each component (rows) of the MDS (columns)
##   gives the variability
##
## Argument(s)
##
##         lst : a list containing the lists resulting from the classical MDS computed with cmdscale.
##   col.names : character vector with the names associated with each MDS

adqCMDS <- function(lst, col.names)
{
   out.ls <- lapply(lst, function(x) {adequacy(x$eig)})
   
   out <- matrix(NA, nrow = max(unlist(lapply(out.ls, length))), ncol = length(out.ls))

   for(j in 1:ncol(out))
   {
       for(i in 1:length(out.ls[[j]]))
       {
           out[i, j] <- out.ls[[j]][i]
       }
   }
   rownames(out) <- paste("Component", 1:nrow(out))
   colnames(out) <- col.names
   
   return(out)   
}

## plotAdq
##
##    Plots the accumulated percentage of variability explained for each number of components m of the
##    representation of a classical MDS 
##
## Argument(s)
##
##     adq : matrix with the accumulated relative frequencies (adequacies) generated wit 'adqCMDS' function
##    main : character with the title of the plot

plotAdq <- function(adq, main = "Adequacy of the m-dimensional Representation")
{
   ptm <- proc.time()

   plot(1:length(adq), adq * 100, main = main, 
        xlab = "Num. of components", ylab = "% of Agreement", lwd = 3, las = 1,
        type = "h", col = "blue", cex.axis = 0.6, cex.lab = 0.6,cex.main = 0.8)

   axis(side = 2, at = adq[2] * 100, paste(round(adq[2] * 100, 2), "%"),
        col.axis = "red", las = 1, font = 1.5, cex.axis = 0.6)
   axis(side = 1, at = 2, "2", col.axis = "red", las = 1, font = 1.5, cex.axis = 0.6)
   points(2, adq[2] * 100, pch = 19, col = "red", cex = 1)

   axis(side = 2, at = adq[3] * 100, paste(round(adq[3] * 100, 2), "%"),
        col.axis = "orange", las = 1, font = 1.5, cex.axis = 0.6)
   axis(side = 1, at = 3, "3", col.axis = "orange", las = 1, font = 1.5, cex.axis = 0.6)
   points(3, adq[3] * 100, pch = 19, col = "orange", cex = 1)

   out <- proc.time() - ptm

   return(out)
}

## plotMDS
##
##    Plots a 2-dimensional configuration point of a classical metric multidimensional scaling
##
## Argument(s)
##
##      mds : matrix with the coordinates, where rows are objects and columns are components. If argument 'adq'
##            is NULL, then coord must be a list object, resulting from isoMDS, with two slots; one, the points,
##            and two, the stress
##      adq : vector with the percentage of variability associated with each component to be ploted in axes
##            labels. 
##    Xlims : vector of tw elemnts with the limits of the X axis
##    Ylims : vector of tw elemnts with the limits of the Y axis
##        g : vector with the numbers of each group at which each object belogns to
##     gcol : color names associated with each group
##   pnames : character vector with the names of ach object to be plotted
##     main : character with the title of the plot
##
## Examples
##
##   swiss.x <- as.matrix(swiss[, -1])
##   swiss.dist <- dist(swiss.x)
##
##   swiss.cmd <- cmdscale(swiss.dist, eig = TRUE)
##   swiss.adq <- adequacy(swiss.cmd$eig)
##   Xlims.cmd <- c(min(swiss.cmd$points[, 1]), max(swiss.cmd$points[, 1]))
##   Ylims.cmd <- c(min(swiss.cmd$points[, 2]), max(swiss.cmd$points[, 2]))
##   plotMDS(mds = swiss.cmd$points, adq = swiss.adq, Xlims = Xlims.cmd, Ylims = Ylims.cmd, g = 1,
##           gcol="blue", pnames = rownames(swiss.x), main = "Classical MDS", cex = 0.8)
##
##   swiss.nmd <- isoMDS(swiss.dist)
##   Xlims.nmd <- c(min(swiss.nmd$points[, 1]), max(swiss.nmd$points[, 1]))
##   Ylims.nmd <- c(min(swiss.nmd$points[, 2]), max(swiss.nmd$points[, 2]))
##   plotMDS(mds = swiss.nmd, adq = NULL, Xlims = Xlims.nmd, Ylims = Ylims.nmd, g = 1,
##           gcol="blue", pnames = rownames(swiss.x), main = "Non-metric MDS", cex = 0.8)

plotMDS <- function(mds, adq, Ylims, Xlims, g, gcol, pnames, main = NULL, cex = 0.8)
{
    ptm <- proc.time()

    if(is.null(adq))
    {
        coord <- mds$points
        
        Xlab <- "Component 1"
        Ylab <- "Component 2"
        main <- paste(main, "(Stress = ", round(mds$stress, 2), "%)")
    }else{
        coord <- mds
        
        Xlab <-  paste0("Component 1 (", round(adq[1] * 100, 2), "%)")
        Ylab <-  paste0("Component 2 (", round((adq[2] - adq[1]) * 100, 2), "%)")
    }
    
    opt <- par(cex = cex)
    plot(coord, xlab = Xlab, ylab = Ylab, type = "n",  main = main, ylim = Ylims, xlim = Xlims)
    par(opt)
    
    for(i in 1:max(g))
    {
        points(coord[g==i, 1],coord[g==i, 2], pch = 19, col = gcol[i])
    }
    
    text(coord[, 1], coord[, 2], pnames, cex = cex, pos = 4)

    out <- proc.time() - ptm

    return(out)
}

## stress
##
##    Yields a numeric vector with the stress percentages of a non-metric MDS calculated with 'isoMDS' function
##
## Argument(s)
##
##   d : distance object
##
## Example(s)
##
##   m <- matrix(rnorm(100, 15, 20), nrow =20, ncol=5)
##   rownames(m) <- letters[1:20]
##   colnames(m) <- paste0("s",1:5)
##
##   d <- dist(m)
##   stress(d)

stress <- function(d)
{
    out <- NULL
    
    n.rows <- nrow(as.matrix(d))
    
    for(i in 1:(n.rows-1))
    {
        if (!any(dim(cmdscale(d, k = i)) != c(n.rows, i)))
        {
            out.i <- isoMDS(d, k = i, trace = FALSE)
            out <- c(out, out.i$stress)
        }
    }
        
    return(out)

}

## screePlot
##
##    Plots the stress percentage of non-metric MDS solutioncalculated with 'isoMDS' function per each dimension
##
## Argument(s)
##
##       d : distance object
##    main : character with the title of the plot
##
##
## Example(s)
##
##   m <- matrix(rnorm(100, 15, 20), nrow =20, ncol=5)
##   rownames(m) <- letters[1:20]
##   colnames(m) <- paste0("s",1:5)
##
##   d <- dist(m)
##   plotStress(d, main = "Stress Plot")

screePlot <- function(d, main = "Scree Plot")
{
    ptm <- proc.time()

    st <- stress(d)

    plot(x = 1:length(st), y = st,
         xlab = "Dimension", ylab = "% of stress", main = main,
         type = "b", pch = 19, col = "blue", cex.axis = 0.6, cex.lab = 0.6,cex.main = 0.8)
    
    axis(side = 2, at = st[2], paste(round(st[2], 2) , "%"),
         col.axis = "red", las = 1, font = 1.5, cex.axis = 0.6)
    axis(side = 1, at = 2, "2", col.axis = "red", las = 1, font = 1.5, cex.axis = 0.6)
    points(2, st[2], pch = 19, col = "red", cex = 1, lwd = 3)
    
    axis(side = 2, at = st[which(st<5)[1]], paste(round(st[which(st<5)[1]], 2), "%"),
         col.axis = "orange", las = 1, font = 1.5, cex.axis = 0.6)
    axis(side = 1, at = which(st<5)[1], which(st<5)[1], col.axis = "orange", las = 1, font = 1.5, cex.axis = 0.6)
    points(which(st<5)[1], st[which(st<5)[1]], pch = 19, col = "orange", cex = 1, lwd = 3)
    
    out <- proc.time() - ptm

    return(out)
}

## shepardPlot
##
##    Plots a Shepard diagram given a the (original) distance matrix and the non-metric MDS solution.
## 
## Argument(s)
##
##       d : distance object
##    nmds : list resulting from the Kruskal's non-metric MDS performed with 'isoMDS' function
##    main : character with the title of the plot
##    xlim : limits of the x axis
##    ylim : limits of the y axis
##
## Example(s)
##
##

shepardPlot <- function(d, nmds, main = "Shepard Diagram", xlim, ylim)
{
    ptm <- proc.time()
    
    sh <- Shepard(d, nmds$points)

    R2.stress <- 1 - ((nmds$stress) / 100) ^ 2

    lm.sh <- lm(sh$yf ~ sh$x)
    R2.linear <- summary(lm.sh)$r.squared
    
    plot(sh, pch = 19, col = "blue", main = main,
         xlim = xlim, ylim = ylim,
         xlab = "Observed Proximities", ylab = "Distance and Disparity")
    abline(lm.sh, col = "black")
    lines(sh$x, sh$yf, type = "S", col = "red", lwd = 2)

    legend("topleft",
       legend = c(paste("R-square (Stress) =", round(R2.stress, 2)),
                  paste("R-square (Linear fit) =", round(R2.linear, 2))),
       lty = 1, pch = 19,
       col = c("red", "black"),
#       ncol = 2,
       bty = "n",
       cex = 0.5,
       text.col = c("red", "black"),
       inset = 0.01)    
    
    out <- proc.time() - ptm

    return(out)
}
