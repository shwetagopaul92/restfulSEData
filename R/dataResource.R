### =========================================================================
### restfulSEData 
### A data.frame to help users understand the datasets
### -------------------------------------------------------------------------
###
datR <- data.frame(
  Datasets = c("banoSEMeta : Metadata RangedSummarizedExperiment shell for banovichSE ", 
            "st100k : Metadata RangedSummarizedExperiment shell for 100k cells from 10x genomics 1.3 million neuron dataset", 
            "st400k : Metadata RangedSummarizedExperiment shell for 400k cells range-sorted from 10xgenomics 1.3 million neuron dataset", 
            "full_1Mneurons : Metadata SummarizedExperiment shell for the full 1.3 million neuron dataset from 10x genomics", 
            "gr450k : GRanges with metadata for illumina 450k methylation assay", 
            "gtexRecount : Metadata RangedSummarizedExperiment shell for RECOUNT gtex rse_gene", 
            "tasicST6 : Supplemental table from Tasic et al. 2016")
  )
#' Convenience functions to explore the datasets
#' @examples  
#' dataResource()
#' @export
dataResource = function(){
  datR$Datasets
}

