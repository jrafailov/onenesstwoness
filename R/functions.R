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
