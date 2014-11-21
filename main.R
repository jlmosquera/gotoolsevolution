##   Thesis Title : Methods and Models for The Analysis Of Biological Significance Based o High-Throughput Data
##  PhD Candidate : Jose Luis Mosquera
## Thesis Advisor : Alex Sanchez
##
## Description : It is the main file to perform the data analysis for studying the evolution and clustering of GO Tools
##               available in SerbGO database. This script is intended to provide main parameters required during the
##               analysis.
##    License : See LICENSE.md for details.

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




