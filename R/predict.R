#' @title run_hrdetect
#' @name run_hrdetect
#' @description function to run HRDetect
#' @param snv Path to the SNV file
#' @param indel Path to the indel file if not present in SNV
#' @param jabba Jabba rds file with $asegstats
#' @param sv Path to the SV file (vcf/rds)
#' @param mask Path to the mask file
#' @param hets Path to the hets file
#' @param genome enum of hg19/hg38/mm10/canFam3
#' @param ref Path to the reference genome
#' @param outdir Path to the output directory
#' @export
#' @importFrom signature.tools.lib HRDetect_pipeline
run_hrdetect <- function(snv,
    indel,
    jabba,
    sv,
    mask,
    hets,
    genome,
    ref,
    outdir) {
    
    message("Running HRDetect")

    setDTthreads(1)

    system(paste('mkdir -p', outdir))

    v <- VariantAnnotation::scanVcfHeader(snv)
    sl <- seqlengths(v)

    wg <- si2gr(sl[which(!grepl("GL|MT|M|NC", names(sl)))]) ## modified line

    wg <- keepStandardChromosomes(wg, pruning.mode = "coarse")
    wg <- sort(sortSeqlevels(wg), ignore.strand = TRUE)

    gr.mask <- withv(mask, {
        if (!is.null(x) && file.exists(x)) {
            message("mask provided")
            if (grepl("bed(.gz)?$", x))
                import(x)
            else if (grepl(".rds$", x))
                readRDS(x)
            else {
                message("mask file in incorrect format")
                message("no region filtering will be done")
                wg[0]
            }
        } else {
            wg[0]
        }
    })

    if (!inherits(gr.mask, "GRanges")) {
        warning("mask specified incorrectly\nproceeding without filtering")
        gr.mask <- wg[0]
    }

    if (!identical(gr.mask, GRanges())) {
        message("parsing mask file, and saving to regions file for subsetting VCF\n")
    } else {
        message("No mask regions provided, using whole SNV / InDel VCFs")
    }

    good.wg <- gr.setdiff(wg, gr.mask) %>% keepStandardChromosomes(pruning.mode = "coarse")

    regions.file <- "./good_rfile.txt"
    regions.bed <- "./good_rfile.bed"
    mask.file <- "./mask_rfile.txt"

    fwrite(gr2dt(good.wg)[, 1:3, with = FALSE], regions.file, sep = "\t", col.names = FALSE)

    fwrite(setnames(gr2dt(good.wg)[, 1:3, with = FALSE][, start := start - 1],
                c("#chrom", "chromStart", "chromEnd")), regions.bed, sep = "\t", col.names = TRUE)

    message("Processing SNVs/indels")

    snv.tmp <- process_snv(snv, regions.bed)
    indel.tmp <- process_indel(indel, regions.bed, ref)
    sv.tmp <- process_sv(jabba, sv)
    cnv.tmp <- process_cnv(jabba, hets)

    ########## run pipeline ##########
    res <- signature.tools.lib::HRDetect_pipeline(SNV_vcf_files = snv.tmp,
                                                Indels_vcf_files = indel.tmp,
                                                SV_bedpe_files = sv.tmp,
                                                CNV_tab_files = cnv.tmp,
                                                genome.v = genome, SNV_signature_version = 'COSMICv3.2')

    saveRDS(res, file.path(outdir, 'hrdetect_results.rds'))

    if (NROW(res$hrdetect_output) > 0)
        fwrite(res$hrdetect_output, file.path(outdir, 'hrdetect_output.txt'))
    else {
        message("HRDetect score was not calculated!!")
        message("Inputs may not be available")
    }
}



#' @title predict_hrd
#' @name predict_hrd
#'
#' @description function to predict the B1+2 score
#'
#' @param model Path to the model file
#' @param complex Path to the complex file
#' @param homeology Path to the homeology file
#' @param homeology_stats Path to the homeology stats file
#' @param hrdetect_results Path to the HRDetect results file
#' @param save Logical indicating whether to save the results
#' @return A list containing the following elements:
#' \itemize{
#' \item{expl_variables}{A data frame containing the explanatory variables}
#' \item{ot_scores}{A data frame containing the Oneness and Twoness scores}
#' }
#' @export
#' @importFrom gUtils grl.in parse.gr gr.fix grl.pivot
#' @import khtools gGnome randomForest
#' @importFrom dplyr mutate
#' @importFrom data.table fread
#' @importFrom reshape2 dcast
#' @importFrom GxG homeology.wrapper
predict_hrd <- function(complex,
    homeology,  
    homeology_stats, 
    hrdetect_results,
    model = system.file("data/model", "stash.retrained.model.rds", package = "onenesstwoness"),
    cores = 4,
    save = TRUE) {
    
    if (is.null(complex) | is.null(homeology) | is.null(homeology_stats) | is.null(hrdetect_results) | is.null(model)) {
        stop("One or more required parameters are missing.")
    }

    gg <- readRDS(complex)
    gg <- gGnome::refresh(gg)
    all.events <- gg$meta$events
    ev.types <- c("qrppos", "qrpmin", "qrpmix", "qrdup", "qrdel", "tib")
    all.events$type <- factor(all.events$type, ev.types)
    expl_variables <- all.events %>% dcast("sample" ~ type, fun.aggregate = length, drop = FALSE)

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

    if(!is.null(homeology) & !is.null(homeology_stats)){
        jhom <- fread(homeology)
        jhom_stats <- fread(homeology_stats)
    } else{
        message("Running homeology")
        hom.run <- homeology.wrapper(
            gg,
            width = 200,
            pad = 20,
            thresh = 1,
            stride = 8,
            savegMatrix = TRUE,
            annotate = FALSE,
            bidirectional = TRUE,
            flip = FALSE,
            cores = cores
        )
        jhom <- hom.run[[3]]
        jhom_stats <- hom.run[[2]]
    }

    message("Processing homeologous dels")
    if (NROW(jhom) > 0) {
        bp1 <- parse.gr(jhom$bp1)
        bp2 <- parse.gr(jhom$bp2)
        bp1 <- gr.fix(bp1, bp2)
        bp2 <- gr.fix(bp2, bp1)
        jhom$jspan <- jJ(grl.pivot(GRangesList(bp1, bp2)))$span
        jhom <- jhom %>% merge(gg$edges[type == "ALT"]$dt, by = "sedge.id", all= TRUE, suffixes = c(".x", "")) %>% as.data.table()
    }
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

    if (save) {
        message("Saving Oneness/Twoness results")
        saveRDS(outputs, "onenesstwoness_results.rds")
    }

    return(outputs)
     
}
