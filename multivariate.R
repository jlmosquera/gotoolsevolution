##   Thesis Title : Methods and Models for the Analysis of Biological Significance on High-Throughput Data
## PhD Candidate  : Jose Luis Mosquera
## Thesis Advisor : Alex Sanchez
##
## Description : This script is devoted to perform a multivariate data analysis based on the data
##               downloaded from SerbGO database and already preprocessed for the analysis. It also
##               allows studying the correlation between the distance matrices.
##     License : Use and distribution subject to the terms of the GPL-2 license. See LICENSE for details.


###########################################
## Chunk 1: Load packages
###########################################

library("vegan")  # vegdist
library("MASS")   # isoMDS
library("fpc")    # pam
library("rgl")

###########################################
## Chunk 2: Read R source functions
###########################################

source(file.path(codeDir, "functions.R"))

###########################################
## Chunk 3: Dissimilarity matrices
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Computation of dissimilarities
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jaccard coefficient

Dist05.jac <- vegdist(M05, method = "jaccard")
Dist07.jac <- vegdist(M07, method = "jaccard")
Dist09.jac <- vegdist(M09, method = "jaccard")

Dist.jac <- list(Dist05.jac, Dist07.jac, Dist09.jac)

## Matching coefficient

Dist05.match <- matching.coeff(M05)
Dist07.match <- matching.coeff(M07)
Dist09.match <- matching.coeff(M09)

Dist.match <- list(Dist05.match, Dist07.match, Dist09.match)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save distance matrices in a single csv
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dist.jac.m <- lapply(Dist.jac, as.matrix)
Dist.match.m <- lapply(Dist.match, as.matrix)

sapply(Dist.jac.m, write.table, file = file.path(resDir, "Dist.Jac.csv"), append = TRUE, sep = ",")
sapply(Dist.match.m, write.table, file = file.path(resDir, "Dist.Math.csv"), append = TRUE, sep = ",")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Ckecking dissimilarities matrices
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jaccard coefficient

symnum(as.matrix(Dist05.jac))
symnum(as.matrix(Dist07.jac))
symnum(as.matrix(Dist09.jac))

## Matching coefficient

symnum(as.matrix(Dist05.match))
symnum(as.matrix(Dist07.match))
symnum(as.matrix(Dist09.match))


###########################################
## Chunk 4: Hierarchical Clustering
###########################################

##-----------------------------------------
## Compute HCs
##-----------------------------------------

## Jaccard coefficient

hc05.jac.avg <- hclust(Dist05.jac, method = "average")
hc07.jac.avg <- hclust(Dist07.jac, method = "average")
hc09.jac.avg <- hclust(Dist09.jac, method = "average")

hc.jac <- list(hc05.jac.avg, hc07.jac.avg, hc09.jac.avg)

## Matching coefficient

hc05.match.avg <- hclust(Dist05.match, method = "average")
hc07.match.avg <- hclust(Dist07.match, method = "average")
hc09.match.avg <- hclust(Dist09.match, method = "average")

hc.match <- list(hc05.match.avg, hc07.match.avg, hc09.match.avg)

##-----------------------------------------
## Optimal number of clusters
##-----------------------------------------

## Silhouettes: A graphical aid to the interpretation and validation of cluster analysis

pdf(file = file.path(resDir, "OptimalClustersJac.pdf"))
   par(mfrow = c(2, 3), cex = 0.8, las = 1, oma = c( 0, 0, 2, 0 ))
   for(i in 1:3)
   {
       plotAvgSilVal(Dist.jac[[i]], main = paste0("20", years[i]))    # k = 2, 9, 2
   }
   for(i in 1:3)
   {
       plotSil(Dist.jac[[i]], main = "Silhouette Plot", cex = 0.5, do.n.k = FALSE)
   }
   title( "Optimal Number of Clusters (Jaccard Coefficient)", outer = TRUE )
dev.off()

pdf(file = file.path(resDir, "OptimalClustersMatch.pdf"))
   par(mfrow = c(2, 3), cex = 0.8, las = 1, oma = c( 0, 0, 2, 0 ))
   for(i in 1:3)
   {
        plotAvgSilVal(Dist.match[[i]], main = paste0("20", years[i])) # k = 3, 16, 10
    }
   for(i in 1:3)
   {
        plotSil(Dist.match[[i]], main = "Silhouette Plot", cex = 0.5, do.n.k = FALSE)
   }
   title( "Optimal Number of Clusters (Matching Coefficient)", outer = TRUE )
dev.off()

##-----------------------------------------
## Cut trees into 3 groups
##-----------------------------------------

## Jaccard coefficient

g05.jac <- cutree(hc05.jac.avg, ngroups)
g07.jac <- cutree(hc07.jac.avg, ngroups)
g09.jac <- cutree(hc09.jac.avg, ngroups)

## Matching coefficient

g05.match <- cutree(hc05.match.avg, ngroups)
g07.match <- cutree(hc07.match.avg, ngroups)
g09.match <- cutree(hc09.match.avg, ngroups)

##-----------------------------------------
## Comparison of clusters
##-----------------------------------------

g <- cbind(g05.jac, g07.jac, g09.jac, g05.match, g07.match, g09.match)
write.csv2(x = g, file = file.path(resDir, "clusters.csv"), row.names = TRUE)


###########################################
## Chunk 5: Multidimensional Scaling
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Coordinates of configuration points
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##-----------------------------------------
## Classical MDS
##-----------------------------------------

## Jaccard coefficient

cm05 <- cmdscale(Dist05.jac, k = 2, eig = TRUE)
cm07 <- cmdscale(Dist07.jac, k = 2, eig = TRUE)
cm09 <- cmdscale(Dist09.jac, k = 2, eig = TRUE)

cm.jac <- list(cm05, cm07, cm09)

## Matching coefficient

cm05.match <- cmdscale(Dist05.match, k = 2, eig = TRUE)
cm07.match <- cmdscale(Dist07.match, k = 2, eig = TRUE)
cm09.match <- cmdscale(Dist09.match, k = 2, eig = TRUE)

cm.match <- list(cm05.match, cm07.match, cm09.match)

##-----------------------------------------
## Kruskal Non-metric MDS
##-----------------------------------------

## Jaccard coefficient

nm05 <- isoMDS(Dist05.jac, k = 2, trace = FALSE)
nm07 <- isoMDS(Dist07.jac, k = 2, trace = FALSE)
nm09 <- isoMDS(Dist09.jac, k = 2, trace = FALSE)

nm.jac <- list(nm05, nm07, nm09)

## Matching coefficient

nm05.match <- isoMDS(Dist05.match, k = 2, trace = FALSE)
nm07.match <- isoMDS(Dist07.match, k = 2, trace = FALSE)
nm09.match <- isoMDS(Dist09.match, k = 2, trace = FALSE)

nm.match <- list(nm05.match, nm07.match, nm09.match)


###########################################
## Chunk 6: Evaluation of MDS solutions
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Adequacy of classical MDS solution
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jaccard coefficient

cm.jac.ls <- list(cm05, cm07, cm09)
adq.jac <- adqCMDS(cm.jac.ls, paste0("20", years))
write.csv2(x = adq.jac, file = file.path(resDir, "adquacy.MDS.jac.csv"), row.names = TRUE)

## Matching coefficient

cm.match.ls <- list(cm05.match, cm07.match, cm09.match)
adq.match <- adqCMDS(cm.match.ls, paste0("20", years))
write.csv2(x = adq.match, file = file.path(resDir, "adquacy.MDS.match.csv"), row.names = TRUE)

##-----------------------------------------
## Adequacy Plots
##-----------------------------------------

pdf(file = file.path(resDir, "adqequacyPlots.pdf"))
   par(mfrow = c(2, 3), oma = c( 0, 0, 2, 0 ), cex = 0.8, las = 1)
   for(i in 1:3)
   {
      plotAdq(adq = adq.jac[, i], main = paste0("20", years[i], "\nJaccard Coefficient"))
   }
   for(i in 1:3)
   {
      plotAdq(adq = adq.match[, i], main = paste0("20", years[i], "\nMatching Coefficient"))
   }
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Scree Plots and Shepard Diagrams
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jaccard Coefficient

pdf(file = file.path(resDir, "ScreeAndShepardPlotsJac.pdf"))
   par(mfrow = c(2, 3), oma = c( 0, 0, 2, 0 ), cex = 0.8, las = 1)
   for(i in 1:3)
   {
      screePlot(d = Dist.jac[[i]], main = paste0("20", years[i]))
   }
   for(i in 1:3)
   {
      shepardPlot(d = Dist.jac[[i]], nmds = nm.jac[[i]],
                  xlim = c(0.3, 1), ylim = c(0, 1.6),
                  main = NULL)
   }
dev.off()

## Matching Coefficient

pdf(file = file.path(resDir, "ScreeAndShepardPlotsMatch.pdf"))
   par(mfrow = c(2, 3), oma = c( 0, 0, 2, 0 ), cex = 0.8, las = 1)
   for(i in 1:3)
   {
      screePlot(d = Dist.match[[i]], main = paste0("20", years[i]))
   }
   for(i in 1:3)
   {
      shepardPlot(d = Dist.match[[i]], nm.match[[i]],
                  xlim = c(0, 0.6), ylim = c(0, 1.3),
                  main = NULL)
   }
dev.off()


###########################################
## Chunk 7: Plots of Clustering Analysis
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Limits of X and Y axes (MDS)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jaccard coefficient

minX <- min(cm05$points[, 1], cm07$points[, 1], cm09$points[, 1],
            nm05$points[, 1], nm07$points[, 1], nm09$points[, 1])
maxX <- max(cm05$points[, 1], cm07$points[, 1], cm09$points[, 1],
            nm05$points[, 1], nm07$points[, 1], nm09$points[, 1])
minY <- min(cm05$points[, 2], cm07$points[, 2], cm09$points[, 2],
            nm05$points[, 2], nm07$points[, 2], nm09$points[, 2])
maxY <- max(cm05$points[, 2], cm07$points[, 2], cm09$points[, 2],
            nm05$points[, 2], nm07$points[, 2], nm09$points[, 2])

Xlims <- c(minX, maxX)     
Ylims <- c(minY, maxY)

## Matching coefficient

minX.match <- min(cm05.match$points[, 1], cm07.match$points[, 1], cm09.match$points[, 1],
                  nm05.match$points[, 1], nm07.match$points[, 1], nm09.match$points[, 1])

maxX.match <- max(cm05.match$points[, 1], cm07.match$points[, 1], cm09.match$points[, 1],
                  nm05.match$points[, 1], nm07.match$points[, 1], nm09.match$points[, 1])
minY.match <- min(cm05.match$points[, 2], cm07.match$points[, 2], cm09.match$points[, 2],
                  nm05.match$points[, 2], nm07.match$points[, 2], nm09.match$points[, 2])
maxY.match <- max(cm05.match$points[, 2], cm07.match$points[, 2], cm09.match$points[, 2],
                  nm05.match$points[, 2], nm07.match$points[, 2], nm09.match$points[, 2])

Xlims.match <- c(minX.match, maxX.match)     
Ylims.match <- c(minY.match, maxY.match)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## HC and MDS plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in 1:3)
{
    pdf(file = file.path(resDir, paste0("HCandMDS.20", years[i], ".pdf")))
        par(mfrow = c(2, 3), oma = c( 0, 0, 2, 0 ), cex = 0.8, las = 1)
    
        ## Jaccard coefficient
    
        plotDendro(hc.jac[[i]], ylab = "Jaccard Coefficient", main = "Hierarchical Clustering",
                   type = "rectangle", labelColors = g.colors, clusMember =  g[,i], cex = 0.4)
        plotMDS(mds = cm.jac[[i]]$points, adq = adq.jac[1:2, i], Ylims = Ylims, Xlims = Xlims,
                main = "Classical MDS", g = g[, i], gcol = g.colors, pnames = GONames, cex = 0.4)
 
        plotMDS(mds = nm.jac[[i]], adq = NULL, Ylims = Ylims, Xlims = Xlims,
                main = "Non-metric MDS", g = g[, i], gcol = g.colors, pnames = GONames, cex = 0.4)
 
        ## Matching coefficient
    
        plotDendro(hc.match[[i]], ylab = "Matching Coefficient", main = "Hierarchical Clustering",
                   type = "rectangle", labelColors = g.colors, clusMember =  g[, i + 3], cex = 0.4)
        plotMDS(mds = cm.match[[i]]$points, adq = adq.match[1:2, i], Ylims = Ylims.match, Xlims = Xlims.match,
                main = "Classical MDS", g = g[, i + 3], gcol = g.colors, pnames = GONames, cex = 0.4)

        plotMDS(mds = nm.match[[i]], adq = NULL, Ylims = Ylims.match, Xlims = Xlims.match,
                main = "Non-metric MDS", g = g[, i + 3], gcol = g.colors, pnames = GONames, cex = 0.4)
    
        title(paste0("20", years[i]), outer = TRUE )
    dev.off()
}


###########################################
## Chunk 6: Mantel Tests
###########################################

##-----------------------------------------
## Single Matel Tests
##-----------------------------------------

## Jaccard coefficient

m05vs07.jac <- mantel(Dist05.jac, Dist07.jac, permutations = 9999)
m07vs09.jac <- mantel(Dist07.jac, Dist09.jac, permutations = 9999)
m05vs09.jac <- mantel(Dist05.jac, Dist09.jac, permutations = 9999)

m.jac <- list(m05vs07.jac, m07vs09.jac, m05vs09.jac)

## Matching coefficient

m05vs07.match <-mantel(Dist05.match, Dist07.match, permutations = 9999)
m07vs09.match <-mantel(Dist07.match, Dist09.match, permutations = 9999)
m05vs09.match <-mantel(Dist05.match, Dist09.match, permutations = 9999)

m.match <- list(m05vs07.match, m07vs09.match, m05vs09.match)

##-----------------------------------------
## Partial Mantel Test
##-----------------------------------------

pm.jac <- mantel.partial(Dist05.jac, Dist09.jac, Dist07.jac, permutations = 9999)
pm.match <- mantel.partial(Dist05.match, Dist09.match, Dist07.match, permutations = 9999)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save Matel Tests in a single csv
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.jac.lst <- lapply(m.jac, function(x){out <- c(r = round(x$statistic, 3), PValue = x$signif)})
m.jac.df <- do.call("rbind", m.jac.lst)
m.jac.df <- rbind(m.jac.df, c(round(pm.jac$statistic, 3), pm.jac$signif))
m.jac.df <- data.frame(Mantel.Test = c("$r_M(D_{2005}^{J},D_{2007}^{J})$",
                                 "$r_M(D_{2007}^{J},D_{2009}^{J})$",
                                 "$r_M(D_{2005}^{J},D_{2009}^{J})$",
                                 "$r_M(D_{2005}^{J},D_{2009}^{J}|D_{2007}^{J})$"),
                 m.jac.df)

m.match.lst <- lapply(m.match, function(x){out <- c(r = round(x$statistic, 3), PValue = x$signif)})
m.match.df <- do.call("rbind", m.match.lst)
m.match.df <- rbind(m.match.df, c(round(pm.match$statistic, 3), pm.match$signif))
m.match.df <- data.frame(Mantel.Test = c("$r_M(D_{2005}^{M},D_{2007}^{M})$",
                                   "$r_M(D_{2007}^{M},D_{2009}^{M})$",
                                   "$r_M(D_{2005}^{M},D_{2009}^{M})$",
                                   "$r_M(D_{2005}^{M},D_{2009}^{M}|D_{2007}^{M})$"),
                 m.match.df)

m.df <- rbind(m.jac.df, m.match.df)

write.csv2(x = m.df, file = file.path(resDir, "mantel.tests.csv"), row.names = FALSE)
