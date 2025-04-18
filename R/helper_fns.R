#' @title bcfindex
#' @name bcfindex
#' @description function to index bcf/vcf file
bcfindex <- function (vcf, force = TRUE) {
    if (!grepl(".[bv]cf(.gz)?$", vcf) & !grepl(".[bv]cf(.bgz)?$", vcf)) {
        stop("check if you have a valid bcf/vcf file")
    }
    if (!file.exists(paste0(vcf, ".tbi")) & !file.exists(paste0(vcf,
        ".csi")) || isTRUE(force)) {
        system(sprintf("bcftools index --tbi %s", vcf))
    }
    vcf
}

#' @title hrdetect_process_snv
#' @name hrdetect_process_snv
#' @description function to process SNV file
hrdetect_process_snv <- function(snv, regions.bed) {
    if (file.exists(snv) && file.info(snv)$size) {
        snv.tmp = "./snv.vcf.gz"

        system2("rm", c(paste0(snv.tmp, "*")), stderr = "/dev/null")

        writeLines(system(sprintf("bcftools query -l %s", snv), intern = T),
                "./excls.txt")

        message("Processing SNV VCF")
        if (grepl("gz$", snv)) {
            cmd = sprintf("{ vcftools --gzvcf %s --remove-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/){ gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools view -v snps | bcftools norm -Oz -m-any; } > %s", snv, normalizePath(regions.bed), snv.tmp)
        } else {
            cmd = sprintf("{ vcftools --vcf %s --remove-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -v snps | bcftools norm -Oz -m-any; } > %s", snv, normalizePath(regions.bed), snv.tmp)
        }
        message("Performing VCFtools and bcftools operations...")
        out = system(cmd)
        if (out != 0)  {
            message("trying snv processing again")
            if (grepl("gz$", snv)) {
                cmd = sprintf("(vcftools --gzvcf %s --remove-indels --remove-filtered-all --recode --stdout | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools norm -Oz -m-any) > %s", snv, snv.tmp)
            } else {
                cmd = sprintf("(vcftools --vcf %s --remove-indels --remove-filtered-all --recode --stdout | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools norm -Oz -m-any) > %s", snv, snv.tmp)
            }
            out = system(cmd)
            if (out != 0) warning("could not be intersected with non-mask regions bed")
        }
        bcfindex(snv.tmp)
        names(snv.tmp) = "sample_1"
    } else {
        snv.tmp = NULL
    }
    return(snv.tmp)
}

#' @title hrdetect_process_indel
#' @name hrdetect_process_indel
#' @description function to process indel file
hrdetect_process_indel <- function(indel, regions.bed, ref, fasta) {
    if (is.null(indel) || !file.exists(indel) || identical(indel, "/dev/null")) {
        message("no proper indel file provided... assuming indels are in snv vcf file")
        indel = snv
    }

    if (!is.null(indel) && file.exists(indel) && file.info(indel)$size) {
        indel.tmp = "./indel.vcf.bgz"
        system2("rm", c(paste0(indel.tmp, "*")), stderr = "/dev/null")

        writeLines(system(sprintf("bcftools query -l %s", indel), intern = T),
                "./excls.txt")

        message("Processing indel input")
        if (grepl("gz$", indel)) {
            cmd = sprintf("{ vcftools --gzvcf %s --keep-only-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b ./good_rfile.bed -header | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools view -v indels | bcftools norm -Ov -m-any | bcftools norm -f %s --check-ref s | bgzip -c; } > ./indel.vcf.bgz", indel, fasta)
        } else {
            cmd = sprintf("{ vcftools --vcf %s --keep-only-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -v indels | bcftools norm -Ov -m-any | bcftools norm -f %s --check-ref s | bgzip -c; } > %s", indel, regions.bed, ref, indel.tmp)
        }

        message("Performing VCFtools and bcftools operations...")
        system(cmd)

        bcfindex(indel.tmp)

        indel_fix = readVcf(TabixFile(indel.tmp))
        rr = rowRanges(indel_fix)
        ins = which(width(rr) == 0 & nchar(as.character(rr$REF)) == 0)
        del = which(width(rr) > 0 & nchar(as.character(unlist(rr$ALT))) == 0)

        if (length(c(ins, del))) {
            message("fixing indels")
        }

        if (length(ins)) {
            rowRanges(indel_fix)[ins] = GenomicRanges::shift(
                gr.resize(rr[ins],
                        1,
                        pad = FALSE,
                        fix = "start"), 1L)
            insref = unname(subseq(unlist(VariantAnnotation::alt(indel_fix)[ins]), 1, 1))
            VariantAnnotation::ref(indel_fix)[ins] = insref
        }

        if (length(del)) {
            rowRanges(indel_fix)[del] = gr.resize(rr[del],
                                                width(rr[del]) + 1,
                                                pad = FALSE, fix = "end")
            delref = unname(get_seq(Rsamtools::FaFile(fasta),
                                    rowRanges(indel_fix)[del]))
            delalt = unname(split(subseq(delref, 1, 1), seq_along(delref)))
            VariantAnnotation::alt(indel_fix)[del] = delalt
            VariantAnnotation::ref(indel_fix)[del] = delref
        }

        suppressWarnings(rm("rr"))

        if (length(c(ins, del))) {
            message("indels fixed... writing to vcf")
            writeVcf(indel_fix, file_path_sans_ext(indel.tmp), index = TRUE)
        }

        suppressWarnings(rm("indel_fix"))
        names(indel.tmp) = "sample_1"
    } else {
        indel.tmp = NULL
    }
    return(indel.tmp)
}

#' @importFrom dplyr select left_join
#' @title hrdetect_process_sv
#' @name hrdetect_process_sv
#' @description function to process SV file
hrdetect_process_sv <- function(jabba, sv) {
    message("Processing SV input")
    if (!(identical(jabba, "/dev/null") && identical(sv, "/dev/null"))) {
        fe.jabba <- !file.not.exists(jabba)
        fe.sv <- !file.not.exists(sv)

        strand.conv <- setnames(rbind(
            data.table(cbind("+", "+"), cbind("-", "+")),
            data.table(cbind("-", "-"), cbind("+", "-")),
            data.table(cbind("+", "-"), cbind("-", "-")),
            data.table(cbind("-", "+"), cbind("+", "+"))
        ), c("new.strand1", "new.strand2", "strand1", "strand2"))

        if (fe.jabba & fe.sv) {
            warning("both jabba and sv inputs provided.")
        }

        if (fe.jabba | fe.sv) {
            gg <- if (fe.jabba) {
                gG(jabba = jabba)
            } else {
                gG(junctions = gGnome:::read.juncs(sv))
            }
        } else {
            stop("must provide at least one of jabba or sv vcf inputs")
        }

        altgrl <- gg$edges[type == "ALT"]$grl
        sv.tmp <- "./sv.bedpe"

        if (length(altgrl)) {
            bdpe <- grl2bedpe(sort.GRangesList(gr.noval(altgrl), ignore.strand = TRUE))
            these.cols <- colnames(bdpe)
            fwrite(dplyr::select(dplyr::left_join(bdpe, strand.conv), !!these.cols, everything()) %>%
                mutate(strand1 = new.strand1,
                        strand2 = new.strand2,
                        new.strand1 = NULL,
                        new.strand2 = NULL,
                        sample = "sample_1"),
                sv.tmp, sep = "\t")
        } else {
            empt <- data.table(chrom1 = "", start1 = "", end1 = "",
                                chrom2 = "", start2 = "", end2 = "",
                                name = "", score = "", strand1 = "", strand2 = "",
                                sample = "")[0,]
            fwrite(empt, sv.tmp, sep = "\t")
        }
        names(sv.tmp) <- "sample_1"
    } else {
        sv.tmp <- NULL
    }
    return(sv.tmp)
}

#' @importFrom dplyr rename_at vars
#' @import MatrixGenerics
#' @title hrdetect_process_cnv
#' @name hrdetect_process_cnv
#' @description function to process CNV file
hrdetect_process_cnv <- function(jabba, hets) {
    if (file.exists(jabba) && !identical(jabba, "/dev/null")) {
        message("Processing CNV input")
        jabd.simple <- readRDS(jabba)
        if(inherits(jabd.simple, "gGraph")) {
            jabd.simple <- gg2jab(jabd.simple)
        }
        jabba.has.hets <- all(c("asegstats", "aadj", "agtrack") %in% names(jabd.simple))
        hets.exist <- !file.not.exists(hets)

        if (isTRUE(hets.exist) && isFALSE(jabba.has.hets)) {
            message("Pulling het_pileups data")
            if (!file.exists(paste(dirname(jabba), 'hets.gr.rds', sep = '/'))) {
                if (!is.null(hets)) {
                    if (is.character(hets)) {
                        if (grepl(".rds$", hets)) {
                            hets <- readRDS(hets)
                        } else {
                            hets <- fread(hets)
                        }
                    } else {
                        hets <- hets
                    }
                    if (!is.data.table(hets)) {
                        hets <- as.data.table(hets)
                    }
                    if (inherits(hets, "data.frame")) {
                        if (!is.null(hets$alt.count.n) & !is.null(hets$ref.count.n)) {
                            hets$ref.frac.n <- hets$alt.count.n / (hets$alt.count.n + hets$ref.count.n)
                            hets.gr <- dt2gr(hets[pmin(ref.frac.n, 1 - ref.frac.n) > 0.2 & (ref.count.n + alt.count.n) >= 2, ])
                            hets.gr$alt <- hets.gr$alt.count.t
                            hets.gr$ref <- hets.gr$ref.count.t
                        } else {
                            hets.gr <- dt2gr(hets)
                            if (all(c("alt", "ref") %in% colnames(hets))) {
                                message("Valid hets already")
                            } else if (all(c("alt.count.t", "ref.count.t") %in% colnames(hets))) {
                                hets.gr$alt <- hets.gr$alt.count.t
                                hets.gr$ref <- hets.gr$ref.count.t
                                hets.gr <- hets.gr %Q% (!is.na(alt)) %Q% (!is.na(ref))
                            } else {
                                message("hets is not in valid format, ignore")
                                hets.gr <- NULL
                            }
                        }
                    } else if (inherits(hets, "GRanges")) {
                        if (all(c("alt.count.t", "ref.count.t") %in% colnames(values(hets)))) {
                            hets.gr <- hets
                            hets.gr$alt <- hets$alt.count.t
                            hets.gr$ref <- hets$ref.count.t
                            hets.gr <- hets.gr %Q% (!is.na(alt) & !is.na(ref))
                        }
                    } else {
                        message("hets is neither data.table nor GRanges, ignore.")
                    }
                    if (!is.null(hets.gr)) {
                        saveRDS(hets.gr, paste(dirname(jabba), '/hets.gr.rds', sep = '/'))
                    }
                }
            } else {
                hets.gr <- readRDS(paste(dirname(jabba), "hets.gr.rds", sep = "/"))
            }
            if (!all(c('asegstats', 'aadj', 'agtrack') %in% names(jabd.simple))) {
                message("Running allelic JaBbA...")
                jaba <- tryCatch(JaBbA:::jabba.alleles(jabd.simple, hets.gr, verbose = TRUE, uncoupled = TRUE)[c('asegstats', 'aadj', 'agtrack')], error = function(e) NULL)
                if (!is.null(jaba)) {
                    jabd.simple <- c(jabd.simple, jaba)
                } else {
                        message("could not estimate allelic CN due to no aberrant edges in graph")
                    if (sum(mcols(jabd.simple$junctions)$cn > 0, na.rm = T) == 0) {
                    } else {
                        message("could not estimate allelic CN... check inputs")
                    }
                    message("assuming high and low are 1/2 of total cn")
                    jabd.simple$asegstats <- grbind(within(jabd.simple$segstats, {cn <- ceiling(cn / 2); type <- "high"}),
                                                    within(jabd.simple$segstats, {cn <- floor(cn / 2); type <- "low"}))
                }
            }
        }
        cnv <- as.data.table(keepStandardChromosomes(jabd.simple$asegstats, pruning.mode = "coarse") %Q%
                                (strand == "+" & !grepl("M|MT", seqnames) & width != 1) %>%
                                sortSeqlevels() %>%
                                sort(ignore.strand = TRUE) %>%
                                gr.nochr)[!duplicated(paste(seqnames, start, end, width, strand, type))]
        cnv <- setDT(rename_at(dcast.data.table(cnv, seqnames + start + end + width + strand ~ type, value.var = "cn"), dplyr::vars(high, low), ~paste0(., "_cn")))
        cnv <- cnv[, .(seg_no = seq_len(.N),
                        Chromosome = seqnames,
                        chromStart = start,
                        chromEnd = end,
                        total.copy.number.inNormal = 2,
                        minor.copy.number.inNormal = 1,
                        total.copy.number.inTumour = high_cn + low_cn,
                        minor.copy.number.inTumour = low_cn)] %>%
            subset2(complete.cases(x))
    }

    if (!exists('cnv') || is.null(cnv)) {
        stop("no allelic CN available...")
    }

    if (is.character(cnv) && cnv == "/dev/null") {
        cnv.tmp <- NULL
    } else {
        cnv.tmp <- "./cna.txt"
        fwrite(cnv, cnv.tmp, sep = "\t")
        names(cnv.tmp) <- "sample_1"
    }
    return(cnv.tmp)
}

file.not.exists = function (x, nullfile = "/dev/null", bad = c(NA, "NA", "NULL", 
    "")) {
    isnul = (is.null(x))
    if(isnul) {
        return(TRUE)
    }
    isbadfile = (x %in% bad | x == nullfile)
    isgoodfile = which(!isbadfile)
    isbadfile[isgoodfile] = !file.exists(as.character(x[isgoodfile]))
    isnolength = len(x) == 0
    return(isnul | isnolength | isbadfile)
}

gg2jab = function (gg, purity = NA, ploidy = NA) {
    if (!inherits(gg, "gGraph")) {
        stop("Must supply gGraph")
    }
    ab.edges = array(NA, dim = c(length(gg$junctions[type == 
        "ALT"]), 3, 2), dimnames = list(NULL, c("from", "to", 
        "edge.ix"), c("+", "-")))
    ab.edges[, 1, 1] = gg$sedgesdt[sedge.id > 0][match(gg$junctions[type == 
        "ALT"]$dt$edge.id, edge.id), from]
    ab.edges[, 2, 1] = gg$sedgesdt[sedge.id > 0][match(gg$junctions[type == 
        "ALT"]$dt$edge.id, edge.id), to]
    ab.edges[, 3, 1] = gg$junctions[type == "ALT"]$dt$edge.id
    ab.edges[, 1, 2] = gg$sedgesdt[sedge.id < 0][match(gg$junctions[type == 
        "ALT"]$dt$edge.id, edge.id), from]
    ab.edges[, 2, 2] = gg$sedgesdt[sedge.id < 0][match(gg$junctions[type == 
        "ALT"]$dt$edge.id, edge.id), to]
    ab.edges[, 3, 2] = gg$junctions[type == "ALT"]$dt$edge.id
    adj.dims = dim(gg$adj)
    if (!is.null(gg$sedgesdt$cn)) {
        adj = Matrix::sparseMatrix(i = gg$sedgesdt[, from], j = gg$sedgesdt[, 
            to], x = gg$sedgesdt[, cn], dims = adj.dims)
    }
    else {
        adj = Matrix::sparseMatrix(i = gg$sedgesdt[, from], j = gg$sedgesdt[, 
            to], x = 1, dims = adj.dims)
    }
    res = list(segstats = gg$gr, adj = adj, junctions = gg$junctions[type == 
        "ALT"]$grl, ab.edges = ab.edges)
    if (!is.na(ploidy)) {
        res$ploidy = ploidy
    }
    else {
        res$ploidy = weighted.mean(gg$nodes$dt[, cn], gg$nodes$dt[, 
            width], na.rm = TRUE)
    }
    if (!is.na(purity)) {
        res$purity = purity
    }
    else {
        res$purity = 1
    }
    return(res)
}
