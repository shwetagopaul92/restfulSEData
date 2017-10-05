### =========================================================================
### restfulSEData metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("banoSEMeta", "st100k", "st400k", "full_1Mneurons", "gr450k", "gtexRecount", "tasicST6"),
            
  Description = c(paste0("Metadata RangedSummarizedExperiment shell for banovichSE"),
                  paste0("Metadata RangedSummarizedExperiment shell for 100k cells from 10x genomics 1.3 million neuron dataset"),
                  paste0("Metadata RangedSummarizedExperiment shell for 400k cells range-sorted from 10xgenomics 1.3 million neuron dataset"),
                  paste0("Metadata RangedSummarizedExperiment shell for the full 1.3 million neuron dataset from 10x genomics"),
                  paste0("GRanges with metadata for illumina 450k methylation assay"),
                  paste0("Metadata RangedSummarizedExperiment shell for RECOUNT gtex rse_gene"),
                  paste0("Supplemental table from Tasic et al. 2016")),
  BiocVersion = c(rep(3.5, 7)),
  Genome = c("hg19", "mm10", "mm10", "mm10", "hg19", "hg19", "mm10"), 
  SourceType = c(rep(".RData", 6), ".xlsx" ),
  SourceUrl = c("10.1371/journal.pgen.1004663", "https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276", "https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276", "https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276", "https://bioconductor.org/packages/release/data/annotation/html/FDb.InfiniumMethylation.hg19.html", "https://jhubiostatistics.shinyapps.io/recount", "http://www.nature.com/doifinder/10.1038/nn.4216"), 
  SourceVersion = "October 5 2017",
  Species = c("Homo Sapiens", "Mus musculus (E18 mice)", "Mus Musculus (E18 mice)", "Mus Musculus (E18 mice)", "Homo sapiens", "Homo Sapiens", "Mus Musculus"), 
  TaxonomyId = c(9606,10090,10090,10090,9606,9606,10090),
  Coordinate_1_based = TRUE,
  DataProvider = c("yriMulti", "10x Genomics", "10x Genomics", "10x Genomics", "Illumina 450 methylation assay", "GTex", "GEO"), 
  Maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
  RDataClass = c(rep("RangedSummarizedExperiment", 4), "GRanges", "RangedSummarizedExperiment", "DataFrame"),
  DispatchClass = rep("rda", 7),
  ResourceName = c(paste0("banoSEMeta.rda"),
                   paste0("st100k.rda"),
                   paste0("st400k.rda"),
                   paste0("full_1Mneurons.rda"),
                   paste0("gr450k.rda"), 
                   paste0("gtexRecount.rda"),
                   paste0("tasicST6.rda"))
  
)

write.csv(meta, file="~/Documents/GITHUB/restfulSEData/inst/extdata/metadata.csv", row.names=FALSE)
