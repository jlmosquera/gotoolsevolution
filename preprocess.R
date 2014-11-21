##   Thesis Title : Methods and Models for the Analysis of Biological Significance on High-Throughput Data
## PhD Candidate  : Jose Luis Mosquera
## Thesis Advisor : Alex Sanchez
##
## Description : This script is devoted to preprocess data files in order to study the evolution of GO tools
##               classified and stored in SerbGO database.
##     License : Use and distribution subject to the terms of the GPL-2 license. See LICENSE for details.

######################################
## Chunk 1: Read R source functions
######################################

source(file.path(codeDir, "functions.R"))


######################################
## Chunk 2: Read data
######################################

labels <- c("main", "type", "species", "data", "annotation", "statistics", "outputs", "tags")

## Capabilities

m05 <- readTables(path = datDir, labels = labels, year = "05")
M05 <- m05$capabilities
M07 <- readTables(path = datDir, labels = labels, year = "07")$capabilities
M09 <- readTables(path = datDir, labels = labels, year = "09")$capabilities

print(t(M05))
print(t(M07))
print(t(M09))

## GO tool names

noms <- m05$toolnames
GOtools <- as.vector(t(noms["nameT"]))

print(GOtools)


######################################
## Chunk 3: Normalization of matrices
######################################

dim(M05) # 26 x 215
dim(M07) # 26 x 215
dim(M09) # 26 x 215

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Homogeneize field names
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Assign M09 field names of M09 to M07

colnames(M07) <- colnames(M09)

## Update field names of M05 to M07 (and M09)

old.field.pos <- grep(FALSE , (colnames(M07) %in% colnames(M05)))
old.field.names <- colnames(M05)[old.field.pos]

colnames(M05)[5] <- "tspotted"

colnames(M05)[19] <- "tfree4fau"
colnames(M05)[56] <- "devc"
colnames(M05)[72] <- "dentrez"
colnames(M05)[81] <- "dclusid"
colnames(M05)[82] <- "did"
colnames(M05)[93] <- "aembl"
colnames(M05)[94] <- "aensembl"
colnames(M05)[107] <- "aentrez"
colnames(M05)[137] <- "adisease"
colnames(M05)[167] <- "xsngset"
colnames(M05)[169] <- "xmltset"
colnames(M05)[180] <- "xtestother"

colnames(M05)[194] <- "oscsv"
colnames(M05)[195] <- "oxls"
colnames(M05)[199] <- "oother"

M05[197] <- apply(M05[196:197],1,sum)

M05[196] <- M07[195]
colnames(M05)[196] <- "ocsv"
aux <- M05[194]
M05[194] <- M05[195]
colnames(M05)[194] <- "oxls"

M05[195] <- M05[196]
colnames(M05)[195] <- "ocsv"
M05[196] <- aux
colnames(M05)[196] <- "oscsv"

## Remove 1st column that contains GO tool id (mgoid)

M05 <- M05[,-1]
M07 <- M07[,-1]
M09 <- M09[,-1]

## Change 9 values corresponding to unvalidated capabilities by 0

M05[M05==9] <- 0 
M07[M07==9] <- 0
M09[M09==9] <- 0

## Assign GO tool names to the rows

rownames(M05)<-noms[,2]
rownames(M07)<-noms[,2]
rownames(M09)<-noms[,2]
M05 <- as.matrix(M05)
M07 <- as.matrix(M07)
M09 <- as.matrix(M09)

## Aggregate " DNA Microarray" functionality

var2erase <- c("tspotted", "tinsitua")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "Automated Updating Source" functionality

var2erase <- c("dweekly", "dmonthly", "dquarter")
M05 <- M05[,-match(var2erase, colnames(M05))]
M07 <- M07[,-match(var2erase, colnames(M07))]
M09 <- M09[,-match(var2erase, colnames(M09))]

## Aggregate "UniGene" functionality

var2erase <- c("dclusid", "did", "dnames", "dugsymbl")
M05[,"dclusid"] <- apply(M05[, var2erase],1,sum)
M05 <- M05[,-match(c("did", "dnames", "dugsymbl"), colnames(M05))]
M05[which(M05[,"dclusid"]>1),"dclusid"] <- 1
colnames(M05)[which(colnames(M05)=="dclusid")] <- "dunigene"

M07[,"dclusid"] <- apply(M07[, var2erase],1,sum)
M07 <- M07[,-match(c("did", "dnames", "dugsymbl"), colnames(M07))]
M07[which(M07[, "dclusid"] > 1), "dclusid"] <- 1
colnames(M07)[which(colnames(M07)=="dclusid")] <- "dunigene"

M09[,"dclusid"] <- apply(M09[, var2erase],1,sum)
M09 <- M09[,-match(c("did", "dnames", "dugsymbl"), colnames(M09))]
M09[which(M09[, "dclusid"] > 1), "dclusid"] <- 1
colnames(M09)[which(colnames(M09)=="dclusid")] <- "dunigene"

## Aggregate "Gene Ontology" functionality

var2erase <- c("amf", "abp", "acc")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[,-match(var2erase, colnames(M07))]
M09 <- M09[,-match(var2erase, colnames(M09))]

## Aggregate "Keyword searching" functionality

var2erase <- c("ablastsc", "aoverlap", "afindids")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "Enrichment GO Terms" functionality

var2erase <- c("xunder", "xover")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "Define cutoff for" functionality

var2erase <- c("xpvalue", "xqvalue")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "Correction for Multiple Testing" functionality

var2erase <- c("xbonferr", "xfdr", "xindep",	"xnoindep", "xfwer", "xholm", "xwestfal", "xothers")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "List limited" functionality

var2erase <- c("onterms",	"othresh")
M05 <- M05[, -match(var2erase, colnames(M05))]
M07 <- M07[, -match(var2erase, colnames(M07))]
M09 <- M09[, -match(var2erase, colnames(M09))]

## Aggregate "Visualization" functionality

var2erase <- c("ovisuali", "oviewgot", "odag", "otree", "obarchar", "oviewpath", "ogobrows", 
               "oamigo", "oownbrow")
M05[,"ovisuali"] <- apply(M05[, var2erase], 1, sum)
M07[,"ovisuali"] <- apply(M07[, var2erase], 1, sum)
M09[,"ovisuali"] <- apply(M09[, var2erase], 1, sum)

var2erase <- var2erase[which(var2erase!="ovisuali")]
M05 <- M05[, -match(var2erase, colnames(M05))]
M05[which(M05[,"ovisuali"]>1),"ovisuali"] <- 1

M07 <- M07[, -match(var2erase, colnames(M07))]
M07[which(M07[,"ovisuali"]>1),"ovisuali"] <- 1

M09 <- M09[, -match(var2erase, colnames(M09))]
M09[which(M09[,"ovisuali"]>1),"ovisuali"] <- 1

######################################
## Chunk 4: Save data
######################################

save(noms, GOtools, M05, M07, M09, file = paste(datDir, "matrices.Rdata", sep="/"))
