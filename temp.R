genes = scan(what=character())
mMat = mData.norm[genes,]
colnames(mMat) = as.character(fGroups)
mMat = mMat[,order(fGroups)]
#mMat = t(scale(t(log(mMat+0.5))))
levels(fGroups)
genes.sym = scan(what=character())
rownames(mMat) = genes.sym


aheatmap(log(mMat[,colnames(mMat) %in% c('WT:15.5', 'Mut:15.5')]+0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.7, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
col=c(brewer.pal(9, 'BrBG')))

aheatmap(log(mMat[,colnames(mMat) %in% c('WT:13.5', 'Mut:13.5')]+0.5), annRow = NA, scale = 'none', Rowv = NA, Colv=NA, cexRow=0.7, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))


genes.sym = scan(what=character())
rownames(mMat) = genes.sym
aheatmap((mMat[,colnames(mMat) %in% c('WT:14.5', 'Mut:14.5')]), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=1, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
col=c('white', brewer.pal(9, 'YlOrRd')))

aheatmap(log(mMat[,colnames(mMat) %in% c('WT:14.5', 'Mut:14.5')]+0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=1, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
col=c(brewer.pal(9, 'BrBG')))


genes = scan(what=character())
mMat = mData.norm[rownames(mData.norm) %in% genes,]
dim(mMat)
colnames(mMat) = as.character(fGroups)
mMat = mMat[,order(fGroups)]
#mMat = t(scale(t(log(mMat+0.5))))
levels(fGroups)
df = AnnotationDbi::select(org.Mm.eg.db, rownames(mMat), columns = 'SYMBOL', keytype = 'ENTREZID')
genes.sym = df$SYMBOL
#genes.sym = scan(what=character())
rownames(mMat) = genes.sym


aheatmap(log(mMat+0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.7, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))

aheatmap(log(mMat[,colnames(mMat) %in% c('WT:13.5', 'WT:14.5', 'WT:15.5')] +0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.8, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))

aheatmap(log(mMat[,colnames(mMat) %in% c('WT:15.5', 'Mut:15.5')] +0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.8, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))

aheatmap(log(mMat[,colnames(mMat) %in% c('WT:15.5', 'Mut:15.5')] +0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.8, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))


aheatmap(log(mMat[,colnames(mMat) %in% c('Mut:13.5', 'Mut:14.5', 'Mut:15.5')] +0.5), annRow = NA, scale = 'row', Rowv = NA, Colv=NA, cexRow=0.8, cexCol = 0.7, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'),
         col=c(brewer.pal(9, 'BrBG')))

# figures
temp = gsub('X', '', as.character(d$split))
head(temp)
d$split = factor(temp)
head(d)
str(d)

d.bk = d[d$split %in% genes,]

df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(d.bk$split), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(as.character(d.bk$split), df$ENTREZID)
df = df[i,]
d.bk$SYMBOL = df$SYMBOL
identical(as.character(d.bk$split), df$ENTREZID)

d.bk$coef = colMeans(mCoef[,d.bk$cols])
temp = gsub('(WT|Mut):.+', '\\1', as.character(d.bk$fBatch))
d.bk$groups = factor(temp)
temp = gsub('WT:|Mut:', '', as.character(d.bk$fBatch))
d.bk$time = factor(temp)
xyplot(coef ~ time | SYMBOL, groups=groups, data=d.bk, type=c('l', 'p'), auto.key=list(columns=2), scales=list(relation='free'),
       ylab='Model Estimated log Deflections from Intercept')

l = tapply(d$cols, d$split, FUN = function(x) {
  ms = colMeans(mCoef[,d.bk$cols])
  r = data.frame(ind= as.character(d$ind[c[base1]]), coef.base1=mean(mCoef[,c[base1]]), 
                 coef.deflection1=mean(mCoef[,c[deflection1]]),
                 coef.base2=mean(mCoef[,c[base2]]), coef.deflection2=mean(mCoef[,c[deflection2]]),
                 zscore=dif$z, pvalue=dif$p)
  r$difference = (r$coef.deflection1 - r$coef.base1) - (r$coef.deflection2 - r$coef.base2)
  #return(format(r, digi=3))
  return(r)
})
