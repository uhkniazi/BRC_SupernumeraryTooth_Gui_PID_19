# Name: createSpatialGmtFile.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 18/03/2019
# Desc: Create a gmt file - after binning each genomic location with the genes annotated in those.


#### logic
## 1 - get list of annotated genes as GRanges objects
## 2 - overlap this with a set of GRanges that are binned genomic coordinates
## 3 - 

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
oGRgenes = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Function: f_bin_vector
# DESC: Takes the start and end values of the vector; and the number of
#       bins, and returns a data.frame with the start and ends of the 
#       bins
# ARGS: start coordinate, end coordinate, and number of bins
# RETS: a data.frame object with starts and ends of each bin
f_bin_vector = function(start, end, bins){
  s = floor(seq(start, end, length.out=bins+1))
  e = s-1
  e[length(e)] = s[length(s)]
  length(s) = length(s)-1
  e = e[2:length(e)]
  return(data.frame(start=s, end=e))
}# f_bin_vector

# Function: f_dfGetMatchingStrandsFromGRanges
# Desc: Match 2 GRanges objects and return a dataframe with matching indices of the 2 ranges
#       including their strand symbols
# Args: 2 GRanges objects i.e. query and subject
# Rets: a data frame with matching indices and their strands
f_dfGetMatchingStrandsFromGRanges = function(oGRquery, oGRsubject){
  # get the ranges without considering strands where the 2 match
  df = findOverlaps(oGRquery, oGRsubject, ignore.strand=T)
  # convert the hits object to a data frame
  df = as.data.frame(df)
  s.query = as.vector(strand(oGRquery[df$queryHits]))
  s.subject=as.vector(strand(oGRsubject[df$subjectHits]))
  df = cbind(df, s.query, s.subject)
  return(df)
} # function


# use the seqlengths information from seqinfo object to create
# binning ranges
ivEnds = seqlengths(oGRgenes)

lGmt = vector(mode = 'list', length(ivEnds))
names(lGmt) = names(ivEnds)

for (i in seq_along(ivEnds)){
  gr = oGRgenes[seqnames(oGRgenes) == names(ivEnds)[i]]
  if (length(gr) == 0) next;
  gr = sort(gr)
  grb = f_bin_vector(1, ivEnds[i], 100) 
  grb = GRanges(names(ivEnds)[i], ranges = IRanges(grb$start, grb$end), strand='*')
  df = f_dfGetMatchingStrandsFromGRanges(gr, grb)
  # add the gene names
  df$names = names(gr)[df$queryHits]
  # create names for each coordinate
  n = sapply(df$subjectHits, function(x) paste(range(grb[x])))
  df$location = n
  l = tapply(df$names, factor(df$location), function(x) paste0(x))
  lGmt[[i]] = l
}

table(sapply(lGmt, is.null))
lGmt.sub = lGmt[!sapply(lGmt, is.null)]

oFile = file('results/mm10GeneChromosome.gmt', 'wt')

for (i in 1:length(lGmt.sub)){
  x = names(lGmt.sub[[i]])
  for(y in seq_along(x)){
    # remove single gene bins
    if(length(lGmt.sub[[i]][[x[y]]]) < 2) next;
    p = paste(x[y], names(lGmt.sub)[i], paste(lGmt.sub[[i]][[x[y]]], collapse='\t'), sep='\t')
    writeLines(p, oFile)
  }
}
close(oFile)
