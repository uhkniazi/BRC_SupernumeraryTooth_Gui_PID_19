# File: 16_cgraphClust.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the common DE genes to create a network
# Date: 15/07/2019


source('header.R')
library(downloader)
library(org.Mm.eg.db)

url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/experimental/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

dfData = read.csv('wtComparisons/results/commonDEGenes.xls', header=T, row.names=1)

url = 'https://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt'
dir.create('dataExternal', showWarnings = F)
csReactomeFile = 'dataExternal/NCBI2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')

dfGraph = data.frame('Gene'=dfReactome$V1, 'Pathway'=dfReactome$V2, 'Description'=dfReactome$V4)
dim(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% rownames(dfData), ]
dim(dfGraph)
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))

oCGbp.reactome = CGraph.bipartite2(dfGraph[,1:2], ivWeights = c(2, 1, 0), mix.prior = c(5, 3, 1))
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c(''), bDropOrphans = T)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red'), bDropOrphans = T)
set.seed(123)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

# #### use the GO database
# dfGraph = AnnotationDbi::select(org.Mm.eg.db,
#                                 rownames(dfData)
#                                 , 'GO', 'ENTREZID')
# dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
# dfGraph = dfGraph[,c('ENTREZID', 'GO')]
# dfGraph = na.omit(dfGraph)
# str(dfGraph)
# length(unique(dfGraph$ENTREZID))
# 
# oCGbp.gobp = CGraph.bipartite(dfGraph, ivWeights = c(2, 1, 0))
# table(E(getProjectedGraph(oCGbp.gobp))$weight)
# plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)
# plot.projected.graph(oCGbp.gobp, cDropEdges = c('red'), bDropOrphans = T)

## plot the graph with respective colours for tooth development genes
cvGenes = scan(what=character())

ig = getProjectedGraph(oCGbp.reactome)
ig = delete.edges(ig, which(E(ig)$weight < 2))
ig = delete.vertices(ig, which(degree(ig) == 0))

cvVertices = V(ig)$name
V(ig)$color = 'lightgrey'
V(ig)[cvGenes[cvGenes %in% cvVertices]]$color = 'lightgreen'

df = select(org.Mm.eg.db, keys = V(ig)$name, columns = 'SYMBOL', keytype = 'ENTREZID')
dim(df)
vcount(ig)
identical(V(ig)$name, df$ENTREZID)

V(ig)$label = df$SYMBOL

pdf('temp/network.pdf')
set.seed(123)
par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label.cex=0.2, layout=layout_with_fr(ig, weights=E(ig)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', vertex.size=5)
#legend('topright', legend = c('TB', 'Sepsis', 'Common'), fill = c('lightgreen', 'moccasin', 'lightskyblue1'))

cvGenes.overlap = cvGenes[cvGenes %in% cvVertices]

p = get.all.shortest.paths(ig, V(ig)[cvGenes.overlap], V(ig)[cvGenes.overlap], weights=E(ig)$green)

ig.sub = induced_subgraph(ig, V(ig)[unique(unlist(p[[1]]))])

plot(ig.sub, vertex.label.cex=1, layout=layout_with_fr(ig.sub, weights=E(ig.sub)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', vertex.size=5)
dev.off(dev.cur())
## extract the neighbourhood, i.e. names of vertices in one order
# v = unique(names(unlist(neighborhood(ig, 1, V(ig)[cvGenes.overlap]))))
# 
# ig.sub = ig
# ig.sub = induced_subgraph(ig.sub, V(ig.sub)[v])
# p = get.all.shortest.paths(ig.sub, V(ig.sub)[cvGenes.overlap], V(ig.sub)[cvGenes.overlap], weights=E(ig.sub)$green)
# ## extract edges of shortest paths to colour them differently
# ivEdges = numeric()
# for (o in 1:length(p$res)){
#   iFirst = p$res[[o]][1]
#   iLast = p$res[[o]][length(p$res[[o]])]
#   # if the vertex shortest path is to itself so array will have 
#   # only one dimension
#   if (iFirst == iLast) next;
#   # make an edge list
#   for (i in 1:(length(p$res[[o]])-1) ) {
#     ivE = c(p$res[[o]][i], p$res[[o]][i+1])
#     ivEdges = c(ivEdges, ivE)
#   } # inner for
# } # outer for
# # get the ids for these edges
# e = get.edge.ids(ig.sub, ivEdges)
# E(ig.sub)$color = 'grey'
# E(ig.sub)[e]$color = 'red'


