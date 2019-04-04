# File: scratch.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 12/03/2019

source('header.R')

## load the data
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

g = scan(what=character())
mData.norm = mData.norm[g,]
dim(mData.norm)

library(lattice)
library(org.Mm.eg.db)
## remove X from annotation names
rownames(mData.norm)
df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(rownames(mData.norm)), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(rownames(mData.norm), df$ENTREZID)
df = df[i,]
identical(rownames(mData.norm), df$ENTREZID)
rownames(mData.norm) = df$SYMBOL
df = data.frame(t(log(mData.norm+0.5)))
df = stack(df)
f1 = factor(dfSample.2$group1, levels = c('WT', 'Mut'))
f2 = factor(dfSample.2$group3)
df$fBatch = factor(f2:f1)
bwplot(values ~ fBatch | ind, data=df, 
       scales=list(x=list(cex=0.5, rot=45), y=list(cex=1)))

