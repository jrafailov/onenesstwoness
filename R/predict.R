#' @title predict_hrd
#' @name predict_hrd
#'
#' @description function to predict the HRD score
#'
#' @param model Path to the model file
#' @param complex Path to the complex file
#' @param homeology Path to the homeology file
#' @param homeology_stats Path to the homeology stats file
#' @param hrdetect_results Path to the HRDetect results file
#' @param libdir Path to the library directory
#' @return
#' @examples
#' @export
predict_hrd <- function(model, complex, homeology, homeology_stats, hrdetect_results, libdir) {
    model <- readRDS(model)
    gg <- readRDS(complex)
    gg <- gGnome::refresh(gg)
    all.events <- gg$meta$events
    ev.types <- c("qrppos", "qrpmin", "qrpmix", "qrdup", "qrdel", "tib")
    all.events$type <- factor(all.events$type, ev.types)
    expl_variables <- all.events %>% reshape2::dcast("sample" ~ type, fun.aggregate = length, drop = FALSE)

    if (NROW(all.events) > 0) {
        .new <- c(expl_variables$qrppos, expl_variables$qrpmin, expl_variables$qrpmix)
        .old <- c(expl_variables$qrdup, expl_variables$qrdel, expl_variables$tib)
        if (!identical(.new, .old)) {
            if (all(.new == 0) && any(.old) > 0) {
                expl_variables$qrppos <- expl_variables$qrdup
                expl_variables$qrpmin <- expl_variables$qrdel
                expl_variables$qrpmix <- expl_variables$tib
            } else if (all(.old == 0) && any(.new) > 0) {
                expl_variables$qrdup <- expl_variables$qrppos
                expl_variables$qrdel <- expl_variables$qrpmin
                expl_variables$tib <- expl_variables$qrpmix
            }
        }
    } else {
        expl_variables$qrppos <- 0
        expl_variables$qrpmin <- 0
        expl_variables$qrpmix <- 0
        expl_variables$qrdup <- 0
        expl_variables$qrdel <- 0
        expl_variables$tib <- 0
    }

    message("Processing homeologous dels")
    jhom <- fread(homeology)
    if (NROW(jhom) > 0) {
        bp1 <- parse.gr(jhom$bp1)
        bp2 <- parse.gr(jhom$bp2)
        bp1 <- gr.fix(bp1, bp2)
        bp2 <- gr.fix(bp2, bp1)
        jhom$jspan <- jJ(grl.pivot(GRangesList(bp1, bp2)))$span
    }
    jhom_stats <- fread(homeology_stats)
    dels <- jhom[!is.na(jhom$del), colnames(jhom), drop = F, with = F]
    if (NROW(jhom_stats)) {
        dels <- merge.repl(
            dels,
            jhom_stats[, .(
                hlen = max(max(ifelse(na2false(.SD$r > 0.9), .SD$minpx, 0L)), 0L)
            ), by = edge.id],
            by = "edge.id"
        )
    }
    num_ihdels <- NROW(dels[dels$hlen >= 10 & dels$jspan > 1000, ])
    expl_variables$ihdel <- num_ihdels

    message("Processing HRDetect inputs: del.mh.prop, RS3, RS5, hrd-LOH score, SNV3, SNV8")
    res <- readRDS(hrdetect_results)
    hrd <- res$data_matrix
    if (!identical(type(hrd), "double")) hrd <- data.matrix(hrd)
    hrd <- as.data.table(hrd)
    hrd <- setcols(hrd, c("SV3", "SV5"), c("RS3", "RS5"))
    expl_variables <- expl_variables %>% mutate(hrd)

    if (!grepl("\\/", model)) {
        mod_path <- paste0(libdir, "/", model)
    } else {
        mod_path <- model
    }
    message("Path to model:\n", mod_path)
    mod <- readRDS(mod_path)

    tmp <- expl_variables[, c("tib", "ihdel", "qrdup", "qrdel", "RS3", "RS5", "del.mh.prop", "hrd", "SNV8", "SNV3"), drop = F]

    ### begin glmnet_utils.R
    ### begin .Rprofile
    options(stringsAsFactors = FALSE)
    options(bitmapType = "cairo")
    options(device = grDevices::png)
    options(scipen = 0)

    #######################
    #######################
    #######################

    Sys.setenv(R_DATATABLE_NUM_THREADS = 1)
    Sys.setenv(R_REMOTES_UPGRADE = "never")
    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
    Sys.setenv("GENCODE_DIR" = "~/DB/GENCODE")

    Sys.setenv("BASH_FUNC_blip()" = "() { echo \"hoohah\"; }")

    Sys.setenv(DEFAULT_GENOME = "~/DB/references/hg19/human_g1k_v37_decoy.chrom.sizes")
    Sys.setenv(DEFAULT_BSGENOME = "~/DB/references/hg19/human_g1k_v37_decoy.chrom.sizes")

    ww <- with
    wn <- within

    # go.R <- function() {
    #     eval(
    #         quote(
    #             .libPaths(
    #                 unique(
    #                     c( ## "/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2_KH",
    #                         .libPaths()
    #                     )
    #                 )
    #             )
    #         ), globalenv()
    #     )
    #     evalq(
    #         {
    #             suppressWarnings({
    #                 suppressPackageStartupMessages({
    #                     ### begin startup.R
    #                     .libPaths(unique(c("/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2_KH", .libPaths())))

    #                     Sys.setenv(DEFAULT_BSGENOME = "/gpfs/commons/groups/imielinski_lab/DB/GATK/human_g1k_v37_decoy.chrom.sizes")
    #                     Sys.setenv(DEFAULT_GENOME = "/gpfs/commons/groups/imielinski_lab/DB/GATK/human_g1k_v37_decoy.chrom.sizes")

    #                     Sys.setenv(PATH = paste("~/software/bcftools-1.9", Sys.getenv("PATH"), sep = ":"))
    #                     Sys.setenv(BCFTOOLS_PLUGINS = "/gpfs/commons/groups/imielinski_lab/Software/bcftools-1.9/plugins")

    #                     startup()
    #                     library(withr)
    #                     library(skitools)
    #                     library(gGnome)
    #                     library(Flow)
    #                     library(dplyr)
    #                     library(bamUtils)
    #                     library(skidb)
    #                     library(naturalsort)
    #                     library(magrittr)
    #                     library(signature.tools.lib)
    #                     library(tidyr)
    #                     library(MASS)
    #                     library(khtools)
    #                     library(wesanderson)
    #                     forceload(.force = T)

    #                     ## require3(
    #                     ##     withr,
    #                     ##     skitools,
    #                     ##     gGnome,
    #                     ##     Flow,
    #                     ##     dplyr,
    #                     ##     bamUtils,
    #                     ##     skidb,
    #                     ##     naturalsort,
    #                     ##     magrittr,
    #                     ##     signature.tools.lib,
    #                     ##     tidyr,
    #                     ##     MASS,
    #                     ##     khtools,
    #                     ##     wesanderson
    #                     ## )

    #                     ## source("~/lab/home/khadi/git/khscripts/patch_gtrack.R")

    #                     exprs <- expression({
    #                         brewer.master <- skitools::brewer.master
    #                         `%$%` <- gUtils::`%$%`
    #                         `%&%` <- gUtils::`%&%`
    #                         ## `%Q%` <- gUtils::`%Q%`
    #                         `%Q%` <- khtools::`%Q%`
    #                         `%+%` <- gUtils::`%+%`
    #                         `%-%` <- gUtils::`%-%`
    #                         `%^%` <- gUtils::`%^%`
    #                         width <- GenomicRanges::width
    #                         reduce <- GenomicRanges::reduce
    #                         between <- data.table::between
    #                         setnames <- data.table::setnames
    #                         select <- dplyr::select
    #                         update <- Flow::update
    #                         cache <- Flow::cache
    #                         set <- data.table::set
    #                         set_names <- rlang::set_names
    #                         matches <- dplyr::matches
    #                         n <- dplyr::n
    #                         last <- dplyr::last
    #                         first <- dplyr::first
    #                         tailf <- khtools::tailf
    #                         coalesce <- khtools::coalesce
    #                         ppdf <- khtools::ppdf
    #                         ppng <- khtools::ppng
    #                         with_libpaths <- withr::with_libpaths
    #                         with_options <- withr::with_options
    #                         melt <- data.table::melt
    #                         overwritefun(khtools::gr.flipstrand, gUtils::gr.flipstrand)
    #                         gr.flipstrand <- khtools::gr.flipstrand
    #                         .S3method("merge", "Junction", merge.Junction)
    #                         registerS3method("merge", "Junction", merge.Junction, envir = globalenv())
    #                         .S3method("dcast", "data.table", dcast.data.table)
    #                         registerS3method("dcast", "data.table", dcast.data.table, envir = globalenv())
    #                     })

    #                     eix <- seq(2, length(exprs[[1]]))
    #                     for (i in eix) {
    #                         print(exprs[[1]][[i]])
    #                         tryCatch(
    #                             {
    #                                 eval(exprs[[1]][[i]], envir = globalenv())
    #                             },
    #                             error = function(e) NULL
    #                         )
    #                     }

    #                     ## tryCatch({
    #                     ##     .S3method("merge", "Junction", merge.Junction)
    #                     ##     registerS3method("merge", "Junction", merge.Junction, envir = parent.frame())
    #                     ## }, error = function(e) NULL)



    #                     ## tryCatch({
    #                     ##     .S3method("dcast", "data.table", dcast.data.table)
    #                     ##     registerS3method("dcast", "data.table", dcast.data.table, envir = parent.frame())
    #                     ## }, error = function(e) NULL)
    #                     ### end startup.R
    #                 })
    #             })
    #         },
    #         globalenv()
    #     )
    # }
}
