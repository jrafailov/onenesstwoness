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
#' @return A list containing the following elements:
#' \itemize{
#' \item{expl_variables}{A data frame containing the explanatory variables}
#' \item{ot_scores}{A data frame containing the Oneness and Twoness scores}
#' }
#' @export
#' @importFrom gUtils grl.in parse.gr gr.fix grl.pivot
#' @import gGnome
#' @importFrom data.table fread
predict_hrd <- function(complex,
    homeology,  
    homeology_stats, 
    hrdetect_results,
    model = system.file("data/model", "stash.retrained.model.rds", package = "onenesstwoness"), 
    libdir) {
    
    if (is.null(complex) | is.null(homeology) | is.null(homeology_stats) | is.null(hrdetect_results) | is.null(model)) {
        stop("One or more required parameters are missing.")
    }

    withr::with_package(c("skitools", "gGnome", "dplyr", "magrittr", "reshape", "khtools"), {

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
            jhom <- jhom %>% merge(gg$edges[type == "ALT"]$dt, by = "sedge.id", all= TRUE, suffixes = c(".x", "")) %>% as.data.table()
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

        mod <- readRDS(model)

        expl_variables$DUP_1kb_100kb <- 0
        if(NROW(jhom[class == "DUP-like"]) > 0) {
            expl_variables$DUP_1kb_100kb <- 
                jhom[class == "DUP-like"][jspan >= 1e3 & jspan <= 1e5, .N]
        }

        ot_scores <- predict(mod, expl_variables, type = "prob")

        message("Predicting Oneness Twoness scores")

        outputs <- list(
            expl_variables = expl_variables,
            ot_scores = ot_scores
        )

        message("Caching Oneness/Twoness results")
        saveRDS(outputs, "onenesstwoness_results.rds")

        return(outputs)
    })
}
