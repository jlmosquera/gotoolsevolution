##   Thesis Title : Methods and Models for the Analysis of Biological Significance Based on High-Throughput Data
## PhD Candidate  : Jose Luis Mosquera
## Thesis Advisor : Alex Sanchez
##
## Description : This script is devoted to perform a descriptive statistic and an inferential analysis
##               based on the data downloaded from SerbGO database and already preprocessed for the analysis.
##     License : See LICENSE.md for details.


###########################################
## Chunk 1: Load packages
###########################################

library("reshape")
library("coin")
library("multcomp")
library("colorspace")

###########################################
## Chunk 2: Read R source functions
###########################################

source(file.path(codeDir, "functions.R"))


###########################################
## Chunk 3: Read data
###########################################

load(file.path(datDir, "matrices.Rdata"))

## Names to facilitate visualization

GONames <- substr(GOtools, 1, 7)
GONames[which(GONames=="Onto-To")] <- "OntoTools"
GONames[which(GONames=="ontolog")] <- "Traverser"

## Names of functionalities fields

func.names <- colnames(M05)
original.table <- as.factor(substr(func.names, 1, 1))
func.names.bytable <- split(x = func.names, f = original.table)

section <- list(tool.type = as.character(unlist(func.names.bytable[c("m", "t")])),
                supp.spec = func.names.bytable$s,
                inputdata = func.names.bytable$d,
                annotatio = func.names.bytable$a,
                stat.anal = func.names.bytable$x,
                out.resul = func.names.bytable$o)
## Matrix dimension

n.tools <- nrow(M05)                # number of tools :  26
m.func <- ncol(M05)       # number of functionalities : 178


###########################################
## Chunk 4: Global Analysis
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Frequencies
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

M05.rows <- apply(M05, 1, sum) 
M07.rows <- apply(M07, 1, sum)  
M09.rows <- apply(M09, 1, sum)
M05.rows.freq <- (M05.rows / m.func)*100
M07.rows.freq <- (M07.rows / m.func)*100
M09.rows.freq <- (M09.rows / m.func)*100

eines <- data.frame(AbsFreq.05 = M05.rows, AbsFreq.07 = M07.rows, AbsFreq.09 = M09.rows)
eines.freq <- data.frame(RelFreq.05 = M05.rows.freq,
                         RelFreq.07 = M07.rows.freq,
                         RelFreq.09 = M09.rows.freq)


GOTools.desc <- data.frame(AbsFreq.05 = eines[, 1], RelFreq.05 = round(eines.freq[, 1], 3),
                           AbsFreq.07 = eines[, 2], RelFreq.07 = round(eines.freq[, 2], 3),
                           AbfFreq.09 = eines[, 3], RelFreq.09 = round(eines.freq[, 3], 3))
rownames(GOTools.desc) <- rownames(eines)

write.csv2(x = GOTools.desc, file = file.path(resDir, "desc.global.csv"), row.names = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = file.path(resDir, "desc.barplot.global.pdf"))
   barplot(as.matrix(t(eines.freq)), names.arg = GONames, horiz = FALSE , beside = TRUE, 
           las = 3, ylim = c(0, 100), col = d.colors, ylab = "% of functionalities",
           cex.names = 0.7, off = 0)
   legend("topleft", legend = c("2005", "2007", "2009"), fill = d.colors)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Chi-square Test for Homogeneity
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eines
test.tools <- t(eines)
rownames(test.tools) <- paste0("20", years)
tools.df <- melt(test.tools)
colnames(tools.df) <- c("year", "GOTool", "N")
tools.df$year <- factor(tools.df$year)


## H0: The distribution of functionalities is the same for the year x and y
## H1: The distribution of functionalities is NOT the same for the year x and y

cs.05vs07 <- chisq.test(test.tools[c("2005", "2007"), ], correct = FALSE)
cs.07vs09 <- chisq.test(test.tools[c("2007", "2009"), ], correct = FALSE)
cs.05vs09 <- chisq.test(test.tools[c("2005", "2009"), ], correct = FALSE)

cs.X2 <- c(cs.05vs07$statistic, cs.07vs09$statistic, cs.05vs09$statistic)
cs.df <- c(cs.05vs07$parameter, cs.07vs09$parameter, cs.05vs09$parameter)
cs.p <- c(cs.05vs07$p.value, cs.07vs09$p.value, cs.05vs09$p.value)
cs.a <- p.adjust(p = cs.p, method = "fdr")
cs <- data.frame(Test = c("2005 vs 2007", "2007 vs 2009", "2005 vs 2009"),
                 Chi.Square = cs.X2, df = cs.df, PValue = cs.p, Adj.Pvalue = cs.a)
         
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit a non-parametric curve (loess)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2005 vs 2007

loess.05vs07 <- loess(eines$AbsFreq.07 ~  eines$AbsFreq.05, data = eines, model = TRUE)
pred.05vs07 <- predict(object = loess.05vs07, newdata = data.frame(year.05 = eines$AbsFreq.05), se = TRUE)

## 2005 vs 2009

loess.05vs09 <- loess(eines$AbsFreq.09 ~  eines$AbsFreq.05, data = eines, model = TRUE)
pred.05vs09 <- predict(object = loess.05vs09, newdata = data.frame(year.05 = eines$AbsFreq.05), se = TRUE)

## 2007 vs 2009

loess.07vs09 <- loess(eines$AbsFreq.09 ~  eines$AbsFreq.07, data = eines, model = TRUE)
pred.07vs09 <- predict(object = loess.07vs09, newdata = data.frame(year.07 = eines$AbsFreq.07), se = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Write results
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink(file = file.path(resDir, "fitting.models.txt"))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Chi-square Test for Homogeneity\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2007\n")
cat("-----------------------------------------\n")
cat("\n")
print(cs.05vs07)
cat("\n")
cat("-----------------------------------------\n")
cat("2007 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
print(cs.07vs09)
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
print(cs.05vs09)
cat("\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("loess\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.05vs09)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.05vs09)
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2007\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.05vs07)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.05vs07)
cat("\n")
cat("-----------------------------------------\n")
cat("2007 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.07vs09)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.07vs09)
sink()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots associated with fitting models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Xmax <- Ymax <- max(eines) + (10 - (max(eines) %% 10))

pdf(file = file.path(resDir, "desc.fit.global.pdf"))

par(mfrow = c(2, 2))

## Global boxplots

boxplot(N ~ year, data=tools.df, ylim = c(0,Ymax), col = d.colors, ylab = "Num. of functinalities")
legend("topleft", legend = c("2005", "2007", "2009"), fill = d.colors)

## Loess 2005 vs 2007

plot(eines$AbsFreq.07 ~ eines$AbsFreq.05, xlab = "2005", ylab = "2007", xlim = c(0, Xmax), ylim = c(0, Ymax),
     pch = 19, col = "royalblue2")

points(eines$AbsFreq.05[order(eines$AbsFreq.05)], pred.05vs07$fit[order(eines$AbsFreq.05)],
       type = "l", col = "tomato")
points(eines$AbsFreq.05[order(eines$AbsFreq.05)],
       pred.05vs07$fit[order(eines$AbsFreq.05)] -
       qt(0.975, pred.05vs07$df) * pred.05vs07$se[order(eines$AbsFreq.05)],
       type= "l", lty = 2, col = "tomato")
points(eines$AbsFreq.05[order(eines$AbsFreq.05)], 
       pred.05vs07$fit[order(eines$AbsFreq.05)] +
       qt(0.975, pred.05vs07$df) * pred.05vs07$se[order(eines$AbsFreq.05)],
       type= "l", lty = 2, col = "tomato")

text(eines$AbsFreq.05, eines$AbsFreq.07, GONames, cex = 0.6, pos = 4)

## Loess 2005 vs 2009

plot(eines$AbsFreq.09 ~ eines$AbsFreq.05, xlab = "2005", ylab = "2009", xlim = c(0, Xmax), ylim = c(0, Ymax),
     pch = 19, col = "royalblue2")

points(eines$AbsFreq.05[order(eines$AbsFreq.05)], pred.05vs09$fit[order(eines$AbsFreq.05)],
       type = "l", col = "tomato")
points(eines$AbsFreq.05[order(eines$AbsFreq.05)], 
       pred.05vs09$fit[order(eines$AbsFreq.05)] -
       qt(0.975, pred.05vs09$df) * pred.05vs09$se[order(eines$AbsFreq.05)],
       type= "l", lty = 2, col = "tomato")
points(eines$AbsFreq.05[order(eines$AbsFreq.05)], 
       pred.05vs09$fit[order(eines$AbsFreq.05)] +
       qt(0.975, pred.05vs09$df) * pred.05vs09$se[order(eines$AbsFreq.05)],
       type= "l", lty = 2, col = "tomato")

text(eines$AbsFreq.05, eines$AbsFreq.09, GONames, cex = 0.6, pos = 4)

## Loess 2007 vs 2009

plot(eines$AbsFreq.09 ~ eines$AbsFreq.07, xlab = "2007", ylab = "2009", xlim = c(0, Xmax), ylim = c(0, Ymax),
     pch = 19, col = "royalblue2")

points(eines$AbsFreq.07[order(eines$AbsFreq.07)], pred.07vs09$fit[order(eines$AbsFreq.07)],
       type = "l", col = "tomato")
points(eines$AbsFreq.07[order(eines$AbsFreq.07)], 
       pred.07vs09$fit[order(eines$AbsFreq.07)] -
       qt(0.975, pred.07vs09$df) * pred.07vs09$se[order(eines$AbsFreq.07)],
       type= "l", lty = 2, col = "tomato")
points(eines$AbsFreq.07[order(eines$AbsFreq.07)], 
       pred.07vs09$fit[order(eines$AbsFreq.07)] +
       qt(0.975, pred.07vs09$df) * pred.07vs09$se[order(eines$AbsFreq.07)],
       type= "l", lty = 2, col = "tomato")

text(eines$AbsFreq.07, eines$AbsFreq.09, GONames, cex = 0.6, pos = 4)

dev.off()

pdf(file = file.path(resDir, "desc.lines.global.pdf"))
   matplot(test.tools, col = rainbow(n.tools), type = 'l', lty = 1, ylim = range(test.tools),
           ylab = "Num. of functinalities", axes = FALSE)
   axis(1, 1:nrow(t(as.matrix(eines))), paste0("20", years))
   axis(2)
   box() #- to make it look "as usual"

   par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
   legend("topleft", colnames(t(as.matrix(eines))), col = rainbow(n.tools), cex = .65, pch = 19)
dev.off()


###########################################
## Chunk 5: Analysis by Sections
###########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Frequencies
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

M05.section <-funcBySection(sec = section, mat = M05)
M07.section <-funcBySection(sec = section, mat = M07)
M09.section <-funcBySection(sec = section, mat = M09)

m.func.section <- lapply(M05.section, ncol)

M05.section.rows <- absFreqBySection(M05.section)
M07.section.rows <- absFreqBySection(M07.section)
M09.section.rows <- absFreqBySection(M09.section)
eines.section.ls <- list(AbsFreq.05 = M05.section.rows,
                         AbsFreq.07 = M07.section.rows,
                         AbsFreq.09 = M09.section.rows)
eines.section.abs <- absFreq.df(eines.section.ls)
eines.section.rel <- relFreq.df(dtf = eines.section.abs, lst = m.func.section)

GOTools.section.desc <- data.frame(eines.section.abs[, 1:3], RefFreq.05 = eines.section.rel[, 3],
                                    AbsFreq.07 = eines.section.abs[, 4], RefFreq.07 =  eines.section.rel[, 4],
                                    AbsFreq.09 = eines.section.abs[, 5], RefFreq.09 =  eines.section.rel[, 5])


write.csv2(x = GOTools.section.desc, file = file.path(resDir, "desc.section.csv"), row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eines.section.ls <- splitBySection(df = eines.section.rel, j = 1)

sec.names <- c("Type Of Tool", "Supported Species", "Input Data", "Annotation Functionalities",
               "Statistical Methods", "Outputs")
names(sec.names) <- names(eines.section.ls)

## barplot.2(lst = eines.section.ls,  sec = sec.names, outDir = resDir, file = "desc.barplot", save.file = TRUE)

pdf(file = file.path(resDir, "desc.barplot.sections.pdf"))
   par(mfrow = c(2, 3), oma = c( 0, 0, 2, 0 ), cex = 0.8, las = 1)
   for(i in names(eines.section.ls))
   {
       barplot(as.matrix(t(eines.section.ls[[i]])), names.arg = GONames, horiz = FALSE , beside = TRUE, 
               las = 3, ylim = c(0, 100), col = c("tomato", "gold",  "royalblue2"),
               ylab = "% of functionalities", border = NA,
               cex.names = 0.4, cex.axis = 0.6, cex.main = 0.8, cex.lab = 0.6,
               main = sec.names[i], off = 0)
       
       legend("topleft", legend = c("2005", "2007", "2009"),
              fill = c("tomato", "lightyellow",  "royalblue2"), cex= 0.6)
   }
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Chi-square Test for Homogeneity
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eines.section.abs.ls <- splitBySection(df = eines.section.abs, j = 1)

test.tools.ls <- lapply(eines.section.abs.ls, function(x, y)
                                          {
                                              out.tools <- t(x)
                                              rownames(out.tools) <- paste0("20", y)
                                              out.df <-  melt(out.tools)
                                              colnames(out.df) <- c("year", "GOTool", "N")
                                              return(out.df)
                                          },
                                          y = years)

## Is there any significant difference due to year?

cs.ls <- lapply(eines.section.abs.ls,
                function(df)
                {
                    idx.all.zero <- which(apply(df==0, 1, sum)>1)
                    if(sum(idx.all.zero)>0)
                    {
                        test.tools <- df[-idx.all.zero, ]
                    }else{
                        test.tools <- df
                    }
                    
                    colnames(test.tools) <- paste0("20", years)                          
                    
                    out.05vs07 <- chisq.test(test.tools[, c("2005", "2007")], correct = FALSE)
                    out.07vs09 <- chisq.test(test.tools[, c("2007", "2009")], correct = FALSE)
                    out.05vs09 <- chisq.test(test.tools[, c("2005", "2009")], correct = FALSE)
                    
                    out.X2 <- c(out.05vs07$statistic, out.07vs09$statistic, out.05vs09$statistic)
                    out.df <- c(out.05vs07$parameter, out.07vs09$parameter, out.05vs09$parameter)
                    out.p <- c(out.05vs07$p.value, out.07vs09$p.value, out.05vs09$p.value)
                    out.a <- p.adjust(p = out.p, method = "fdr")
                    out <- data.frame(Test = c("2005 vs 2007", "2007 vs 2009", "2005 vs 2009"),
                                      Chi.Square = out.X2, df = out.df, PValue = out.p, Adj.Pvalue = out.a)
                    
                    return(out)
                })

cs.df <- do.call("rbind", cs.ls)

cs.sec.df <- data.frame(Section = rep(sec.names, each = 3, times = 1), cs.df)

write.csv2(x = cs.sec.df, file = file.path(resDir, "chisq.sec.csv"), row.names = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit a non-parametric curves (loess)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2005 vs 2007

loess.05vs07.ls <- lapply(eines.section.abs.ls, function(x, cols)
                                            {
                                                df <- as.data.frame(x)
                                                out <- loess(df[, cols[2]] ~  df[, cols[1]], data = df,
                                                             model = TRUE)
                                                return(out)
                                            },
                                            cols = c("AbsFreq.05", "AbsFreq.07"))

pred.05vs07.ls <- lapply(names(loess.05vs07.ls),
                         function(sec, loess.ls, eines.ls, col)
                         {
                             predict(object = loess.ls[[sec]],
                                     newdata = data.frame(year.05 = eines.ls[[sec]][, col]),
                                     se = TRUE)
                         },
                         loess.ls = loess.05vs07.ls,
                         eines.ls = eines.section.abs.ls,
                         col = "AbsFreq.05")
names(pred.05vs07.ls) <- names(loess.05vs07.ls)

## 2005 vs 2009

loess.05vs09.ls <- lapply(eines.section.abs.ls, function(x, cols)
                                            {
                                                df <- as.data.frame(x)
                                                out <- loess(df[, cols[2]] ~  df[, cols[1]], data = df,
                                                             model = TRUE)
                                                return(out)
                                            },
                                            cols = c("AbsFreq.05", "AbsFreq.09"))

pred.05vs09.ls <- lapply(names(loess.05vs09.ls),
                         function(sec, loess.ls, eines.ls, col)
                         {
                             predict(object = loess.ls[[sec]],
                                     newdata = data.frame(year.05 = eines.ls[[sec]][, col]),
                                     se = TRUE)
                         },
                         loess.ls = loess.05vs09.ls,
                         eines.ls = eines.section.abs.ls,
                         col = "AbsFreq.05")
names(pred.05vs09.ls) <- names(loess.05vs09.ls)



## 2007 vs 2009

#loess.07vs09 <- loess(eines$AbsFreq.09 ~  eines$AbsFreq.07, data = eines, model = TRUE)
#pred.07vs09 <- predict(object = loess.07vs09, newdata = data.frame(year.07 = eines$AbsFreq.07), se = TRUE)



## 2007 vs 2009

loess.07vs09.ls <- lapply(eines.section.abs.ls, function(x, cols)
                                            {
                                                df <- as.data.frame(x)
                                                out <- loess(df[, cols[2]] ~  df[, cols[1]], data = df,
                                                             model = TRUE)
                                                return(out)
                                            },
                                            cols = c("AbsFreq.07", "AbsFreq.09"))

pred.07vs09.ls <- lapply(names(loess.07vs09.ls),
                         function(sec, loess.ls, eines.ls, col)
                         {
                             predict(object = loess.ls[[sec]],
                                     newdata = data.frame(year.07 = eines.ls[[sec]][, col]),
                                     se = TRUE)
                         },
                         loess.ls = loess.07vs09.ls,
                         eines.ls = eines.section.abs.ls,
                         col = "AbsFreq.07")
names(pred.07vs09.ls) <- names(loess.07vs09.ls)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Write results
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink(file = file.path(resDir, "fitting.models.section.txt"))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("loess\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.05vs09.ls)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.05vs09.ls)
cat("\n")
cat("-----------------------------------------\n")
cat("2005 vs 2007\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.05vs07.ls)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.05vs07.ls)
cat("\n")
cat("-----------------------------------------\n")
cat("2007 vs 2009\n")
cat("-----------------------------------------\n")
cat("\n")
cat("Local Polynomial Regression Fitting .....\n")
cat("\n")
print(loess.07vs09.ls)
cat("\n")
cat("Model Predictions .......................\n")
cat("\n")
print(pred.07vs09.ls)

sink()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots associated with fitting models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lapply(names(test.tools.ls),
       function(sec, lst.test.tools, lst.eines, lst.pred, lst.pred.5vs7, lst.pred.5vs9, lst.pred.7vs9,
                sec.name, outDir, file, save.file)
       {
           test.tools.df <- lst.test.tools[[sec]]
           eines.df <- lst.eines[[sec]]
           pred.5vs7   <- lst.pred.5vs7[[sec]]
           pred.5vs9   <- lst.pred.5vs9[[sec]]
           pred.7vs9   <- lst.pred.7vs9[[sec]]

           if(!(save.file))
           {
               dev.new()
           }else{
               if(is.null(outDir))
               {
                   pdf(file = paste0(file, ".pdf"))
                   out.plot <- paste("Plots associated with", sec.name[sec], "saved in", paste0(file, ".pdf"))
               }else{
                   pdf(file = file.path(outDir, paste(file, sec.name[sec], "pdf", sep = ".")))
                   out.plot <- paste("Plots associated with", sec.name[sec],
                                     "saved in", file.path(outDir, paste(file, sec.name[sec], "pdf", sep = ".")))
               }
           }
           
           Xmax <- Ymax <- max(eines.df) + (10 - (max(eines.df) %% 10))
           
           par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
          
           # boxplot

           boxplot(N ~ year, data = test.tools.df,ylim = c(0, Ymax), col = d.colors,
                   ylab = "Num. of functinalities")
           legend("topleft", legend = paste0("20", years), fill = d.colors)

           # plot loess 2005 vs 2007

           plot(eines.df[, "AbsFreq.07"] ~ eines.df[, "AbsFreq.05"], xlab = "2005", ylab = "2007",
                xlim = c(0, Xmax), ylim = c(0, Ymax), pch = 19, col = "royalblue2")
           
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs7$fit[order(eines.df[, "AbsFreq.05"])],
                  type = "l", col = "tomato")
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs7$fit[order(eines.df[, "AbsFreq.05"])] -
                  qt(0.975, pred.5vs7$df) * pred.5vs7$se[order(eines.df[, "AbsFreq.05"])],
                  type= "l", lty = 2, col = "tomato")
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs7$fit[order(eines.df[, "AbsFreq.05"])] +
                  qt(0.975, pred.5vs7$df) * pred.5vs7$se[order(eines.df[, "AbsFreq.05"])],
                  type= "l", lty = 2, col = "tomato")
           
           text(eines.df[, "AbsFreq.05"], eines.df[, "AbsFreq.07"], rownames(eines.df), cex = 0.6, pos = 4)


           # plot loess 2005 vs 2009

           plot(eines.df[, "AbsFreq.09"] ~ eines.df[, "AbsFreq.05"], xlab = "2005", ylab = "2009",
                xlim = c(0, Xmax), ylim = c(0, Ymax), pch = 19, col = "royalblue2")
           
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs9$fit[order(eines.df[, "AbsFreq.05"])],
                  type = "l", col = "tomato")
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs9$fit[order(eines.df[, "AbsFreq.05"])] -
                  qt(0.975, pred.5vs9$df) * pred.5vs9$se[order(eines.df[, "AbsFreq.05"])],
                  type= "l", lty = 2, col = "tomato")
           points(eines.df[, "AbsFreq.05"][order(eines.df[, "AbsFreq.05"])], 
                  pred.5vs9$fit[order(eines.df[, "AbsFreq.05"])] +
                  qt(0.975, pred.5vs9$df) * pred.5vs9$se[order(eines.df[, "AbsFreq.05"])],
                  type= "l", lty = 2, col = "tomato")
           
           text(eines.df[, "AbsFreq.05"], eines.df[, "AbsFreq.09"], rownames(eines.df), cex = 0.6, pos = 4)
 
           # plot loess 2007 vs 2009

           plot(eines.df[, "AbsFreq.09"] ~ eines.df[, "AbsFreq.07"], xlab = "2007", ylab = "2009", 
                xlim = c(0, Xmax), ylim = c(0, Ymax), pch = 19, col = "royalblue2")
           
           points(eines.df[, "AbsFreq.07"][order(eines.df[, "AbsFreq.07"])], 
                  pred.7vs9$fit[order(eines.df[, "AbsFreq.07"])], type = "l", col = "tomato")
           points(eines.df[, "AbsFreq.07"][order(eines.df[, "AbsFreq.07"])], 
                  pred.7vs9$fit[order(eines.df[, "AbsFreq.07"])] -
                  qt(0.975, pred.7vs9$df) * pred.7vs9$se[order(eines.df[, "AbsFreq.07"])],
                  type= "l", lty = 2, col = "tomato")
           points(eines.df[, "AbsFreq.07"][order(eines.df[, "AbsFreq.07"])], 
                  pred.7vs9$fit[order(eines.df[, "AbsFreq.07"])] +
                  qt(0.975, pred.7vs9$df) * pred.7vs9$se[order(eines.df[, "AbsFreq.07"])],
                  type= "l", lty = 2, col = "tomato")
           
           text(eines.df[, "AbsFreq.07"], eines.df[, "AbsFreq.09"], rownames(eines.df), cex = 0.6, pos = 4)
         
           mtext(sec.name[sec], outer = TRUE, cex = 1.5)

            if(save.file) dev.off()
       },
       lst.test.tools = test.tools.ls,
       lst.eines = eines.section.abs.ls,
       lst.pred.5vs7 = pred.05vs07.ls,
       lst.pred.5vs9 = pred.05vs09.ls,
       lst.pred.7vs9 = pred.07vs09.ls,
       sec.name = sec.names,
       outDir = resDir,
       file = "desc.fit",
       save.file = TRUE)
