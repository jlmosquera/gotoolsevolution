## Title of Thesis : Methods And Models for The Analysis Of Biological Significance On
##                   High-Throughput Data
## PhD Candidate : Jose Luis Mosquera
##    Supervisor : Alex Sanchez
##
##      Script : DATA ANALYSIS for studying the evolution and clustering of GO Tools available
##               in SerbGO database
## Description : In order to study the evolution of the 30 GO tools originally included in
##               SerbGO database along five years, namely, 2005, 2007 and 2009, three transversal
##               sections over a this period of time have been considered.

######################################
## Chunk 1: Global parameters
######################################

  preprocess <- FALSE
 descriptive <- TRUE
multivariate <- TRUE

######################################
## Chunk 2: Directories
######################################

## Working directory
           
wrkDir <- "~/gotoolsevolution"
setwd(wrkDir)

## Input directories

datDir <- file.path(wrkDir, "data")

## Output directories

resDir <- file.path(wrkDir, "results")
tmpDir <- file.path(wrkDir, "tmp")

## R Coding directory

codeDir <- file.path (wrkDir, "Rcode")

######################################
## Chunk 3: Pre-processing data
######################################

years <- c("05", "07", "09")

if(preprocess) source(file.path(wrkDir, "preprocess.R"))

#######################################
## Chunk 4: Data Analysis
#######################################

## Descriptive and Inference Statistics

d.colors <- c("tomato", "lightyellow",  "royalblue2")
colores <- colorRampPalette(d.colors)

if(descriptive) source(file.path(wrkDir, "descriptive.R"))

## Multivariate Statistics

ngroups <- 3
g.colors <-  c("tomato", "gold", "royalblue2")

if(multivariate) source(file.path(wrkDir, "multivariate.R"))




