#' onenesstwoness - predictor of HRD phenotype
#'
#' This package provides tools for predicting HRD (Homologous Recombination Deficiency) phenotype 
#' using various genomic and statistical methods. It integrates multiple libraries and utilities 
#' to facilitate genomic data analysis and modeling.
#'
#' @name onenesstwoness-package
#' @docType package
#' @importFrom gUtils grl.in parse.gr gr.fix grl.pivot
#' @import khtools
#' @import randomForest
#' @importFrom dplyr mutate
#' @importFrom data.table fread
#' @importFrom reshape2 dcast
#' @importFrom signature.tools.lib HRDetect_pipeline
#' @importFrom rtracklayer TwoBitFile
#' @import GxG
#' @import BiocGenerics 
#' @import S4Vectors
#' @import GenomicRanges
#' @import data.table
#' @import gGnome
#' @importFrom magrittr `%>%`
#' @import VariantAnnotation
#' @importFrom Rsamtools TabixFile
registerS3method("merge", "data.table", data.table::merge.data.table)
 "_PACKAGE"
