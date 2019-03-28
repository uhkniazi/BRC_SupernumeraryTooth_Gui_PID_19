# File: 10_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 13/03/2018

## load the data
source('header.R')
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
g_did = 39
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 39) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 39')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

fReplicates = sapply(strsplit(dfSample$description, ';'), function(x) return(x[2]))
dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))

## first normalise the data
# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 5703 18413 
mData = mData[!(i< 3),]
dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))

## delete sample section after testing
mData.norm = round(mData.norm, 0)

# set.seed(123)
# i = sample(1:nrow(mData.norm), 600, replace = F)
# dfData = data.frame(t(mData.norm[i,]))

dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
str(dfSample.2)
dfData$fTreatment = factor(dfSample.2$group1, levels = c('WT', 'Mut'))
dfData$fTime = factor(dfSample.2$group3)
dfData$fBatch = dfData$fTreatment:dfData$fTime

dfData$Coef = factor(dfData$fBatch:dfData$ind)

dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef), ]
str(dfData)

## load the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load(file='temp/fit.stan.nb_13Mar.rds')
print(fit.stan, c('sigmaRan1[1]', 'iSize'), digits=3)
print(fit.stan, c('rGroupsJitter1'))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

getDifferenceVector = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  return(d)
}



## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
d$`1` = d$`1`:d$`2`
d = d[,-4]
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base1='Mut:13.5', deflection1='Mut:14.5', base2='WT:13.5', deflection2='WT:14.5') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif1 = getDifferenceVector(ivData = mCoef[,c[deflection1]], ivBaseline = mCoef[,c[base1]])
  dif2 = getDifferenceVector(ivData = mCoef[,c[deflection2]], ivBaseline = mCoef[,c[base2]])
  dif = getDifference(ivData = dif1, ivBaseline = dif2)
  r = data.frame(ind= as.character(d$ind[c[base1]]), coef.base1=mean(mCoef[,c[base1]]), 
                 coef.deflection1=mean(mCoef[,c[deflection1]]),
                 coef.base2=mean(mCoef[,c[base2]]), coef.deflection2=mean(mCoef[,c[deflection2]]),
                 zscore=dif$z, pvalue=dif$p)
  r$difference = (r$coef.deflection1 - r$coef.base1) - (r$coef.deflection2 - r$coef.base2)
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
library(org.Mm.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))
df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'X', fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'X', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'X')
table(dfResults$adj.P.Val < 0.1)
## save the results 
setwd(gcswd)
setwd('interactionContrasts/')
dir.create('results')
write.csv(dfResults, file='results/DEAnalysisInterContMut14.13VsWT14.13.xls')
setwd(gcswd)
