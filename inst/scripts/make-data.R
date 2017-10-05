### =========================================================================
### restfulSEData 
### -------------------------------------------------------------------------
###

#banoSEMeta : metadata RangedSummarizedExperiment shell for banovichSE

library("yriMulti")
data("banovichSE")
banovichSE@assays = Assays(SimpleList())
save(banovichSE, file = "banovichSEMeta.rda")

# st100k.rda : metadata RangedSummarizedExperiment shell for 100k cells from 10x genomics 1.3 million neuron dataset
## Download the "1M_neurons_filtered_gene_bc_matrices_h5.h5" file 
## from url : "https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276"

library(EnsDb.Mmusculus.v79)
library(TENxGenomics)
tx = TENxGenomics("1M_neurons_filtered_gene_bc_matrices_h5.h5") #load the 1M_neurons_filtered_gene_bc_matrices_h5.h5 file
t100k = as.matrixSummarizedExperiment(tx[,1:100000]) # create a SummarizedExperiment with the required number of cells(100000,400000,1306127) )
save(st100k, file="st100k.rda") #raw 
allg = genes(EnsDb.Mmusculus.v79)
dum28k = GRanges(rep("1", 27998), IRanges(rep(2,27998), width=0))
mcols(dum28k) = DataFrame(gene_id=NA_character_, gene_name=NA_character_,
                          entrezid=NA_character_, gene_biotype=NA_character_, seq_coord_system=
                            NA_character_, symbol=NA_character_)
seqlevels(dum28k) = seqlevels(allg)
seqinfo(dum28k) = seqinfo(allg)
names(dum28k) = rownames(t10k)
okn = intersect(names(dum28k), names(allg))
ranges(dum28k[okn]) = ranges(allg[names(dum28k[okn])])
seqnames(dum28k[okn]) = seqnames(allg[names(dum28k[okn])])
mcols(dum28k[okn]) = mcols(allg[okn])
nna = which(!is.na(dum28k$gene_id))
all.equal(names(dum28k)[nna], dum28k[nna]$gene_id)
rowRanges(t100k) = dum28k 
st100k = sort(t100k)
st100k@assays = Assays(SimpleList()) # to retrieve the metadata alone
save(st100k, file="st100k.rda")


# st400k.rda : metadata RangedSummarizedExperiment shell for 400k cells from 10x genomics 1.3 million neuron dataset
## Following the same procedure as above, to subset 400k cells, add rowRanges and retreive metadata alone.
t400k = as.matrixSummarizedExperiment(tx[,1:400000])
st400k@assays = Assays(SimpleList()) # to retrieve the metatdata alone
save(st1400k, file="st400k.rda")

# full_1Mneurons full 1.3 million cell dataset from 10x genomics 
## Following the same procedure as above, to subset 400k cells, add rowRanges, no sorting, retreive metadata alone.
full_1Mneurons = as.matrixSummarizedExperiment(tx[1:27998,1:1306127])
save(full_1Mneurons, file="full_1Mneurons.rda")

# gr450k.rda : GRanges with metadata for illumina 450k methylation assay}
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
## Subsetting three metadata columns we are interested in.
gr450k <- hm450[, c("channel","percentGC","probeType")]
save(gr450k, file="gr450k.rda")

# gtexRecount.rda : metadata RangedSummarizedExperiment shell for RECOUNT gtex rse_gene
## Download the data from the source url :https://jhubiostatistics.shinyapps.io/recount
rse_gene@assays = Assays(SimpleList())
save(rse_gene, file="gtexRecount.rda")

# tasicST6.rda  : Supplemental table from Tasic et al. 2016
#source url : http://www.nature.com/doifinder/10.1038/nn.4216, supplementary table 6. 
library(xlsx)
library(plyr)
df <- read.xlsx("nn.4216-S8 (1).xlsx","final_summary_12") 
data_st6 <- subset(df, select=c("Final.Cluster.ID","Transcriptomic.type","NA.","Present.Markers"))
tasicST6 <- rename(data_st6, c("Final.Cluster.ID"="clid", "Transcriptomic.type"="txtype1", "NA."="txtype2", "Present.Markers"="genes"))
save(tasicST6, file = "tasicST6.rda")