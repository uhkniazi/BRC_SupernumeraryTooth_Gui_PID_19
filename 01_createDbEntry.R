# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 4/3/2019


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
setwd('dataExternal/')
setwd('raw/')
setwd('fastq/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz')

# each sample has 2 files 
fSplit = gsub('_[1|2].fastq.gz', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metaData.csv', header=T, stringsAsFactors = F)
dfMeta$Title = gsub(' ', '', dfMeta$Title)
dfMeta$Sample = gsub(' ', '', dfMeta$Sample)
dfMeta$Genotype = gsub(' ', '', dfMeta$Genotype)
str(dfMeta)
# sanity check
table(as.character(dfMeta$Title) %in% unique(fSplit))
## order the table in the same sequence as file names
i = match(names(lFiles), as.character(dfMeta$Title))
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$Title), names(lFiles))
identical(as.character(dfMeta$Title), unique(fSplit))

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=unique(fSplit), 
                       description= paste('sample name', as.character(dfMeta$Sample),  
                                          'group1 is Genotype',
                                          'group2 is technical replicate',
                                          'group3 is developmental stage', sep=';'),
                       group1 = dfMeta$Genotype, group2= dfMeta$Replicate, group3=dfMeta$Stage)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# # write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 39;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
# get the names of the samples
temp = lapply(dfSamples$title, function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[dfSamples$title == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
