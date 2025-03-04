# Function to preprocess data
my_preprocess <- function(dat) {
  lapply(dat, function(x) {
    if (inherits(x, c("numeric", "integer"))) {
      unlist(mystats(x))
    }
  })
}

# Function to predict data
my_predict <- function(mypre, newdat, method = "robust") {
  goods <- sapply(mypre, Negate(is.null))
  nm <- names(goods[goods])
  qm <- qmat(newdat, , nm)
  if (identical(method, "robust")) {
    mp <- mapply(function(x, y) {
      (x - y["x.median"]) / (y["x.q90"] - y["x.q10"])
    }, qm, mypre[nm], SIMPLIFY = F)
  }
  do.assign(newdat, setColnames(as.data.frame(mp), nm))
}

# Function to fit random forest model
fit_rforest <- function(mod, data) {
  data <- copy(data)
  classnames <- mod$classes
  scores <- as.data.frame.matrix(predict(mod, data, type = "prob"))
  out <- copy(data)
  for (i in 1:NCOL(scores)) {
    out[[colnames(scores)[[i]]]] <- scores[[i]]
  }
  out$classes <- predict(mod, data, type = "class")
  out$NRSUM <- rowSums(qmat(out, , classnames[-1]))
  return(out)
}


#' @title predict.pp
#' @description This function predicts using preprocessed data.
#'
#' @param pp.res A list containing the preprocessed data.
#' @param newdat A data frame containing the new data to be predicted.
#' @param apply.prefun A logical value indicating whether the pre-processing function should be applied.
#' @return A list containing the preprocessed data, the new data, the pre-processing function, the variables, the method, and the original new data.
#' @export 
# Function to predict using preprocessed data
predict.pp <- function(pp.res, newdat, apply.prefun = TRUE) {
  require(caret)
  require(dplyr)
  odat <- copy3(newdat)
  if (is.function(apply.prefun) || (!is.null(pp.res$prefun) && isTRUE(apply.prefun) && is.function(pp.res$prefun))) {
    if (is.function(apply.prefun)) pp.res$prefun <- apply.prefun
    for (v in pp.res$vars) {
      newdat[[v]] <- pp.res$prefun(newdat[[v]])
    }
  } else {
    pp.res$prefun <- NULL
  }
  if (!is.null(pp.res$pp) && !identical(pp.res$method, "robust")) {
    newdat <- predict(pp.res$pp, newdat)
  } else if (identical(pp.res$method, "robust")) {
    newdat <- my_predict(pp.res$pp, newdat)
  }
  return(structure(list(pp = pp.res$pp, newdat = newdat, prefun = pp.res$prefun, vars = pp.res$vars, method = pp.res$method, orig_newdat = odat), class = "carpost"))
}

#' @title mystats
#' @description This function calculates various statistics for a given numeric vector.
#'
#' @param x A numeric vector for which the statistics are to be calculated.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return A list of statistics including mean, median, standard deviation, median absolute deviation, variance, minimum, maximum, and various quantiles.
#' @examples
#' mystats(c(1, 2, 3, 4, 5)) # returns a list of statistics
#' mystats(c(NA, 2, 3, 4, 5), na.rm = TRUE) # returns a list of statistics with NA removed
#' @export
mystats <- function(x, na.rm = TRUE) {
  if (is.null(dim(x))) {
    x <- as.data.frame(x)
  }
  fun <- function(x) {
    if (is.numeric(x)) {
      c(
        mean = mean(x, na.rm = na.rm),
        median = median(x, na.rm = na.rm),
        sd = sd(x, na.rm = na.rm),
        mad = mad(x, na.rm = na.rm),
        var = stats::var(x, na.rm = na.rm),
        min = min(x, na.rm = na.rm),
        max = max(x, na.rm = na.rm),
        q25 = unname(quantile(x, 0.25, na.rm = na.rm)),
        q75 = unname(quantile(x, 0.75, na.rm = na.rm)),
        q10 = unname(quantile(x, 0.10, na.rm = na.rm)),
        q90 = unname(quantile(x, 0.90, na.rm = na.rm))
      )
    } else {
      c(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        mad = NA_real_,
        var = NA_real_,
        min = NA_real_,
        max = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_,
        q10 = NA_real_,
        q90 = NA_real_
      )
    }
  }
  cn <- colnames(x)
  out.lst <- lapply(x, fun)
  out.lst
}




forceload = function(envir = globalenv(), .force = FALSE) {
    if (!exists("envload_34507213048974")) {
        envload_34507213048974 = new.env(parent = globalenv())
        globasn(envload_34507213048974)
        sesh =  sessionInfo()
        pkgs = c(sesh$basePkgs,
                 names(sesh$otherPkgs),
                 names(sesh$loadedOnly))
        pkvec = rep(FALSE, length(pkgs))
        names(pkvec) = pkgs
        envload_34507213048974$pkvec = pkvec
    }
    force = function(x) x
    ## pkgs = gsub("package:", "", grep('package:', search(), value = TRUE))
    ## pkgs = c(pkgs, names(sessionInfo()$loadedOnly))
    if (!exists("sesh")) sesh =  sessionInfo()
    if (!exists("pkgs")) {
        pkgs = c(sesh$basePkgs,
                 names(sesh$otherPkgs),
                 names(sesh$loadedOnly))
    }
    pkvec = envload_34507213048974$pkvec
    notloaded_firsttime = setdiff(pkgs, names(pkvec))
    pkvec = c(pkvec, setNames(rep_len(FALSE, length(notloaded_firsttime)), notloaded_firsttime))
    if (.force)
        notloaded = names(pkvec)
    else
        notloaded = names(which(pkvec == FALSE))
    if (length(notloaded)) {
        for (pkg in notloaded) {
            tryCatch( {
                message("force loading ", pkg)
                ## invisible(eval(as.list((asNamespace(pkg))), envir = envir))
                ## invisible(eval(eapply(asNamespace(pkg), force, all.names = TRUE), envir = envir))
                invisible(eval(parse(text = sprintf("as.list((asNamespace(\"%s\")))", pkg)), envir = envir))
                invisible(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = envir))
                pkvec[pkg] = TRUE
            }, error = function(e) message("could not force load ", pkg))
        }
        envload_34507213048974$pkvec = pkvec
    } else {
        message("nothing to forceload")
    }
    invisible()
}


require3 = function (...)
{
    names2 = function(x) {
        nm = names(x)
        if (is.null(nm))
            return(rep_len("", length(x)))
        else
            return(nm)
    }
    suppressMessages(forceload(.force = T))
    lst.arg = as.list(match.call(expand.dots = F))$`...`
    nm = names2(lst.arg)
    otherarg = lst.arg[nzchar(nm)]
    pkgarg = lst.arg[!nzchar(nm)]
    pkgarg = pkgarg[sapply(pkgarg, function(x) is.call(x) || is.character(x))]
    charvec = as.character(all.vars(match.call()))
    if (length(charvec)) {
        notfound= { set.seed(10); paste0("notfound_", round(runif(1) * 1e9)); }
        vars = mget(charvec, ifnotfound=notfound, mode = "character", inherits = T)
        ## charvec = unlist(strsplit(toString(vars[[1]]), ", "))
        charvec = unique(c(names2(vars[vars == notfound]), unlist(vars[vars != notfound])))
    }
    charvec = c(charvec, unlist(as.vector(sapply(pkgarg,
                                                 function(x) tryCatch(eval(x), error = function(e) NULL)))))
    for (lib in charvec) {
        pev = packageEvent(lib, "onLoad")
        gh = getHook(pev)
        if (length(gh) == 0 || is.null(gh$forceall12340987)) {
            setHook(pev,
                    list("forceall12340987" = function(...) forceall(envir = asNamespace(lib))))
        }
        ## do.call(require, c(alist(package = lib, character.only = T),
        ##                    otherarg))
        if (NROW(otherarg)) {
            is.char = sapply(otherarg, is.character)
            otherarg[is.char] = paste0("\"", otherarg[is.char], "\"")
            otherargs = paste(paste(names(otherarg), "=", unlist(otherarg)), collapse = ",")
            eval(parse(text = sprintf("require(%s,%s)", lib, otherargs)), globalenv())
        } else {
            eval(parse(text = sprintf("require(%s)", lib)), globalenv())
        }
    }
    suppressMessages(forceload(.force = T))
    invisible()
}


globasn = function (obj, var = NULL, return_obj = TRUE, envir = .GlobalEnv,
                     verbose = TRUE, vareval = F)
{
    var = as.list(match.call())$var
    if (is.null(var)) {
        globx = as.character(substitute(obj))
    }
    else {
        if (is.name(var)) {
            if (isFALSE(vareval))
                var = as.character(var)
            else var = eval(var, parent.frame())
        }
        else if (!is.character(var)) {
            stop("var must be coercible to a character")
        }
        if (inherits(var, "character")) {
            globx = var
        }
        else {
            globx = as.character(substitute(var))
        }
    }
    if (verbose)
        message("variable being assigned to ", globx)
    assign(globx, value = obj, envir = envir)
    if (return_obj) {
        invisible(obj)
    }
    else {
        invisible()
    }
}


forceall = function(invisible = TRUE, envir = parent.frame(), evalenvir = parent.frame()) {
    if (!exists("envload_34507213048974")) {
        envload_34507213048974 = new.env(parent = globalenv())
        globasn(envload_34507213048974)
        sesh =  sessionInfo()
        pkgs = c(sesh$basePkgs,
                 names(sesh$otherPkgs),
                 names(sesh$loadedOnly))
        pkvec = rep(FALSE, length(pkgs))
        names(pkvec) = pkgs
        envload_34507213048974$pkvec = pkvec
    }
    pkg = environmentName(envir)
    pkvec = envload_34507213048974$pkvec
    if ( { pkg %in% names(pkvec) && isFALSE(pkvec[pkg]); } ||
         { ! pkg %in% names(pkvec); } ) {
        if (invisible == TRUE)  {
            ## invisible(eval(as.list(envir), envir = evalenvir))
            ## invisible(eval(eapply(envir, force, all.names = TRUE), envir = evalenvir))
            invisible(eval(parse(text = sprintf("as.list(asNamespace(\"%s\"))", pkg)), evalenvir))
            invisible(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = evalenvir))
        } else {
            ## print(eval(as.list(envir), envir = evalenvir))
            ## print(eval(eapply(envir, force, all.names = TRUE), envir = evalenvir))
            print(eval(parse(text = sprintf("as.list(asNamespace(\"%s\"))", pkg)), evalenvir))
            print(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = evalenvir))
        }
        addon = TRUE
        names(addon) = pkg
        envload_34507213048974$pkvec = c(pkvec, addon)
    } else {
        message("nothing to load")
    }
}

