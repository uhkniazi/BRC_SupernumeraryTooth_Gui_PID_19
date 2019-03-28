# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 14/03/2019

source('header.R')
setwd(gcswd)
setwd('interactionContrasts/')
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
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.1,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.1,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.1,]


library(VennDiagram)

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), 
             rownames(dfContrast2.sub), 
             rownames(dfContrast3.sub)
             )
names(ldfData)
names(lVenn) = cvTitle
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
venn.diagram(lVenn, filename = 'results/venn_all_contrasts.tif', margin=0.1)

### repeat the analysis but separate up and down regulated genes
## select significant genes
dfContrast1.up = dfContrast1.sub[dfContrast1.sub$logFC > 0, ]
dfContrast1.down = dfContrast1.sub[dfContrast1.sub$logFC < 0, ]

dfContrast2.up = dfContrast2.sub[dfContrast2.sub$logFC > 0, ]
dfContrast2.down = dfContrast2.sub[dfContrast2.sub$logFC < 0, ]

dfContrast3.up = dfContrast3.sub[dfContrast3.sub$logFC > 0, ]
dfContrast3.down = dfContrast3.sub[dfContrast3.sub$logFC < 0, ]

# create a list for overlaps
lVenn = list(rownames(dfContrast1.up), rownames(dfContrast1.down),
             rownames(dfContrast2.up), rownames(dfContrast2.down), 
             rownames(dfContrast3.up), rownames(dfContrast3.down)
             )


#cvTitle = gsub('results//DEAnalysis(\\w+)VsControl.xls', '\\1', names(ldfData))
cvTitle.up = paste(cvTitle, 'up', sep='-')
cvTitle.down = paste(cvTitle, 'down', sep='-')
o = c(1, 4, 2, 5, 3, 6)
cvTitle = c(cvTitle.up, cvTitle.down)[o]
names(lVenn) = cvTitle

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
setwd(gcswd)
