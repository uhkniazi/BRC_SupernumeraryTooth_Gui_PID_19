# File: scratch.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: random collection of tasks and code snippets
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

## make plots of genes of choice
g = scan(what=character())
g = AnnotationDbi::select(org.Mm.eg.db, keys = 'Gas1', columns='ENTREZID', keytype = 'SYMBOL')

mData.norm = mData.norm[g$ENTREZID,]
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

bwplot(values ~ fBatch , data=df, 
       scales=list(x=list(cex=0.5, rot=45), y=list(cex=1)))


########## print stan coefficients of certain genes
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

library(org.Mm.eg.db)
d$ind2 = gsub('X', '', as.character(d$ind))
str(d)
d$SYMBOL = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(d$ind2), columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
str(d)
head(d)

## 
d.bk = d
i = which(d$SYMBOL == 'Gas1')

# set outliers to NA
df = mCoef[,d$cols[i]]
colnames(df) = as.character(d$fBatch[i])
df = apply(df, 2, function(x){
  x[x < quantile(x, 0.01) | x > quantile(x, 0.99)] = NA
  return(x)
})

df = data.frame(df)
df = stack(df)
df$f1 = factor(dfSample.2$group1, levels = c('WT', 'Mut'))
df$f2 = factor(dfSample.2$group3)


histogram( ~ exp(values) | ind, data=df, scales=list(relation='free'))
histogram( ~ values | ind, data=df, scales=list(relation='free'))

## format data for plotting
# set outliers to NA
df = mCoef[,d$cols[i]]
colnames(df) = as.character(d$fBatch[i])
# df = apply(df, 2, function(x){
#   x[x < quantile(x, 0.01) | x > quantile(x, 0.99)] = NA
#   return(x)
# })
mModules = df
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = rownames(d)
## split this factor into sub factors
f = strsplit(d$mods, ':')
d = cbind(d, do.call(rbind, f))
colnames(d) = c(colnames(d)[1:5], c('treatment', 'time'))
d$time = factor(d$time)

dotplot(m+s1+s2 ~ treatment | time, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=1), scales=list(relation='free'), 
        main='', ylab='Gas1 - Model estimated Coefficients (log scale)', xlab='')

dotplot(m+s1+s2 ~ treatment | time, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=1), #scales=list(relation='free'), 
        main='', ylab='Gas1 - Model estimated Coefficients (log scale)', xlab='')
## example of xyplot with confidence interval bars
xyplot(m ~ treatment | time, data=d,
       panel=function(x,y,ulim,llim, ...) {
         lsegments(x, s1, x, s2, lwd=1)
         panel.xyplot(x,y, ...)
       }
       , type=c('p'), pch=19, groups=treatment, cex=0.5,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5)), main='41 Modules - baseline',
       auto.key = list(columns=3))


