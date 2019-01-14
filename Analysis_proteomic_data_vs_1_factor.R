setwd("your/repository")
datafr = read.csv2("example_data_set.csv", sep=";", dec=".")

library(broom) ## glance used in plotCorr function
library(factoextra) ## fviz used in plotCorr function
library(seriation) ## seriate used in plotCorr function
library(superheat) ## superheat used in plotCorr function
library(Hmisc) ## rcorr to get the correlation and the pvalue of the correlation



varnames = colnames(datafr)[-(1 : 3)] ## vector of variable names
data.prot = datafr[, varnames]
var.factor = datafr[, "group"] 
var.lim = datafr[, "adiposity"] #Need to be changed according to your "discriminant value" i.e. you need to refer to the name of the 3rd column
nvars = ncol(data.prot)
alpha = 0.05 
df.res = data.frame(matrix(ncol=0, nrow=ncol(data.prot)))

#Shapiro test if the data follows a normale distribution
df.res$pval.shapiro = sapply(1 : nvars, function(i) shapiro.test(data.prot[, i])$p.value)
df.res$normality = sapply(df.res$pval.shapiro, function(x) ifelse(x < alpha, F, T))

#Then, according to True or False, T test (for TRUE value, meaning they are following a normal distribution) or Kruskal-Wallis (for FALSE results) are performed
df.res$pval.testTorKW = sapply(1 : nvars, function(i) {
  if (df.res$normality[i]) glance(lm(data.prot[, i] ~ var.factor))$p.value
  else glance(kruskal.test(data.prot[, i] ~ var.factor))$p.value
})
df.res$signif.testTorKW = sapply(df.res$pval.testTorKW, function(x) ifelse(x < alpha, T, F)) #tell you if it's significative or not

#Give  results of T test and Kruskal-Wallis in any case 
df.res$pval.Ttest = sapply(1 : nvars, function(i) glance(lm(data.prot[, i] ~ var.factor))$p.value)
df.res$signif.Ttest = sapply(df.res$pval.testT, function(x) ifelse(x < alpha, T, F))
df.res$pval.KW = sapply(1 : nvars, function(i) glance(kruskal.test(data.prot[, i] ~ var.factor))$p.value)
df.res$signif.KW = sapply(df.res$pval.KW, function(x) ifelse(x < alpha, T, F))

#Creation in the user working directory of a csv file with all the results (p-value and if it is or not considered as significative)
rownames(df.res)=varnames
write.csv2(df.res,file="results_normality_significativity.csv")

source("plotCorr.R")
plotCorr(var.lim, data.prot)

data.prot.group=datafr[,c(-1,-2)]
data.prot.mat=as.matrix(data.prot.group)

matricecorre=rcorr(data.prot.mat, type="pearson")

pvalcorr=matricecorre$P
correlation=matricecorre$r
l=ncol(pvalcorr)

write.csv2(correlation,file="correlation_matrix.csv")
write.csv2(pvalcorr,file="p-value_correlation_matrix.csv")
