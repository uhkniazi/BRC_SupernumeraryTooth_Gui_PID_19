# File: 16_cgraphClust.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the common DE genes to create a network
# Date: 20/03/2019


source('header.R')
library(downloader)
library(org.Mm.eg.db)

url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/experimental/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

dfData = read.csv('results/commonDEGenes.xls', header=T, row.names=1)

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

oCGbp.reactome = CGraph.bipartite(dfGraph[,1:2], ivWeights = c(2, 1, 0))
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c(''), bDropOrphans = T)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red'), bDropOrphans = T)
set.seed(123)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

#### use the GO database

dfGraph = AnnotationDbi::select(org.Mm.eg.db,
                                rownames(dfData)
                                , 'GO', 'ENTREZID')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
dfGraph = dfGraph[,c('ENTREZID', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$ENTREZID))

oCGbp.gobp = CGraph.bipartite(dfGraph, ivWeights = c(2, 1, 0))
table(E(getProjectedGraph(oCGbp.gobp))$weight)
plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)
plot.projected.graph(oCGbp.gobp, cDropEdges = c('red'), bDropOrphans = T)

# create a template graph for making correlation graphs
ig = CGraph.union(getProjectedGraph(oCGbp.reactome),
                  getProjectedGraph(oCGbp.gobp))
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 1))
vcount(ig)
ecount(ig)
ig.p = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig.p)
plot(ig.p, vertex.label=NA, vertex.size=2, layout=layout_with_fr, vertex.frame.color=NA)

### load the count matrix from previous script 10_deAnalysis.R
# create 2 correlation graphs
oCGcor = CGraph.cor(ig.template = ig.p, mCor = cor(mCounts), ivWeights = c(2, 1, 0)) 
table(E(getProjectedGraph(oCGcor))$weight)
plot.projected.graph(oCGcor)

ig = CGraph.union(getProjectedGraph(oCGbp.reactome),
                  getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGcor))
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 1))
vcount(ig)
ecount(ig)
ig = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig)

set.seed(123)
plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$weight))


## save the graph object in graphml format to use in cytoscape
write.graph(ig, file= 'temp/gui_graph_multiWeights.graphml', format='graphml')


