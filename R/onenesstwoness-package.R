#' onenesstwoness - predictor of HRD phenotype
#'
#' @importFrom gUtils grl.in parse.gr gr.fix grl.pivot
#' @import khtools
#' @import randomForest
#' @import data.table
#' @import GenomicRanges
#' @importFrom dplyr mutate
#' @importFrom data.table fread
#' @importFrom reshape2 dcast
#' @importFrom signature.tools.lib HRDetect_pipeline
#' @importFrom GxG homeology.wrapper
#' @importFrom rtracklayer TwoBitFile
#' @import BiocGenerics S4Vectors
#' @import gGnome
#' @importFrom magrittr `%>%`
registerS3method("merge", "data.table", data.table::merge.data.table)
 "_PACKAGE"
