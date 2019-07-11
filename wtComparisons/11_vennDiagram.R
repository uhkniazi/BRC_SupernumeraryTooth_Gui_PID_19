# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 10/07/2019

source('header.R')

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('results//DEAnalysis', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

for (i in seq_along(cvTitle)){
  df = ldfData[[i]]
  hist(df$logFC, xlab='Log Fold Change', ylab='', main=paste('Fold Changes ', cvTitle[i], sep=''))
  f_plotVolcano(df, cvTitle[i], 0.01, fc.lim=range(df$logFC))
}

names(ldfData)
i = grep('Mut', names(ldfData))
ldfData[i] = NULL
cvTitle = cvTitle[-i]
names(ldfData)
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.1,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.1,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.1,]

library(VennDiagram)

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub))
names(ldfData)
names(lVenn) = cvTitle
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
setwd('wtComparisons/')
venn.diagram(lVenn, filename = 'results/venn_all_contrasts_WTTime.tif', margin=0.1)

###########################################################
############ load the count matrix and normalise
###########################################################
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
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

mData.orig = mData
## first normalise the data
# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))

###########################################################

# subset the data only to wt
dfSample = dfSample.2[dfSample.2$group1 == 'WT',]
mData = mData.norm[, as.character(dfSample$fReplicates)]
identical(colnames(mData), as.character(dfSample$fReplicates))

# grouping factor, time
fGroups = factor(dfSample$group3)
levels(fGroups)

names(lVenn)

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = cvCommonGenes %in% lVenn[[i]]
}
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

# create groups in the data based on ncol^2-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.1)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## put these results together
dfCommonGenes = data.frame(mCommonGenes, sig.pvals=rowSums(mCommonGenes), groups=cp, Symbol=ldfData[[1]][cvCommonGenes, 'SYMBOL'])
## gene names have an X before them remove those if present
rownames(dfCommonGenes) = (gsub(pattern = '^X', replacement = '', x = rownames(dfCommonGenes)))
head(dfCommonGenes)

write.csv(dfCommonGenes, file='results/commonDEGenes.xls')

################################################################################
###### Gui's questions
################################################################################
###1) Is it possible to obtain a figure 
### for trends in gene regulation for WT x WT comparisons between
### all the stages (E13.5, E14.5 and E15.5)? Like a Heatmap figure?
library(NMF)
library(RColorBrewer)
mMat = mData[rownames(dfCommonGenes), ]
colnames(mMat) = as.character(fGroups)
mMat = mMat[,order(fGroups)]
#mMat = t(scale(t(log(mMat+0.5))))
levels(fGroups)
rownames(mMat) = dfCommonGenes$Symbol

pdf('results/heatmap_all_1009.pdf')
aheatmap(log(mMat+0.5), labRow=NA, annRow = NA, scale = 'row', Rowv = T, Colv=NA, cexRow=0.7, cexCol = 0.7, 
         #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))
dev.off(dev.cur())

#### 2) Could we obtain graphs like the interaction contrast ones for the DE genes
####    throughout dental development (E13.5 to E15.5 – WT)?
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load(file='../temp/fit.stan.nb_13Mar.rds')

## format data to extract
dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
str(dfSample.2)
identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))
dfData$fTreatment = factor(dfSample.2$group1, levels = c('WT', 'Mut'))
dfData$fTime = factor(dfSample.2$group3)
dfData$fBatch = dfData$fTreatment:dfData$fTime
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef), ]
str(dfData)

mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)

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

### extract the genes that need plotting
temp = gsub('X', '', as.character(d$split))
head(temp)
d$split = factor(temp)
head(d)
str(d)
dfGenes = dfCommonGenes[dfCommonGenes$groups == 3,]
head(dfGenes)
d.bk = d[as.character(d$split) %in% rownames(dfGenes),]
## drop the mut samples
i = grep('Mut', as.character(d.bk$fBatch))
d.bk = d.bk[-i,]
d.bk = droplevels.data.frame(d.bk)
library(org.Mm.eg.db)
df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(d.bk$split), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(as.character(d.bk$split), df$ENTREZID)
df = df[i,]
d.bk$SYMBOL = df$SYMBOL
identical(as.character(d.bk$split), df$ENTREZID)
head(d.bk)
d.bk$coef = colMeans(mCoef[,d.bk$cols])
temp = gsub('WT:|Mut:', '', as.character(d.bk$fBatch))
d.bk$time = factor(temp)
library(lattice)
xyplot(coef ~ time | SYMBOL, data=d.bk, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='Model Estimated log Deflections from Intercept', main=list(label='40 Genes DE expressed at 3 time points in WT', cex=0.8))

#### 3) To translate the data into a more meaningful biological context and to 
# characterize more thoroughly sets of functionally related genes, is it possible
# to organize the differentially expressed datasets into gene ontology groupings (figures)?
library(GOstats)

goTest = function(cvSeed, univ = keys(org.Mm.eg.db, 'ENTREZID')){
  ## set up universe background
  dfUniv = AnnotationDbi::select(org.Mm.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
  dfUniv = na.omit(dfUniv)
  univ = unique(dfUniv$ENTREZID)
  
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(cvSeed),
               annotation='org.Mm.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = hyperGTest(params) 
  return(oGOStat)
}

lGO.results = lapply(unique(dfCommonGenes$groups), function(group){
  return(goTest(rownames(dfCommonGenes)[dfCommonGenes$groups == group]))
})

oFile.go = file('results/GO_groups.csv', 'wt')
temp = sapply(unique(dfCommonGenes$groups), function(group){
  p1 = paste('Contrast Comparison ', group)
  df = summary(lGO.results[[group]])
  p2 = paste(colnames(df), collapse = ',')
  writeLines(p1, oFile.go)
  writeLines(p2, oFile.go)
  sapply(1:10, function(x){
    p3 = gsub(',', replacement = '-', df[x,])
    p3 = paste(p3, collapse=',')
    writeLines(p3, oFile.go)
  })
})

close(oFile.go)

################################################################################


