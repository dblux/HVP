#' @include utils.R
#'
NULL

#' Recursive variance partitioning (RVP)
#'
#' `RVP` calculates the proportion of variance associated with batch effects
#' in a data set (the "RVP" value of a data set). To determine whether batch
#' effects are statistically significant in a data set, a permutation test can
#' be performed by setting `nperm` to a number above 100. `RVP` is
#' an S3 generic function; methods can be added for new classes. S3 methods
#' for class: array-like objects (default), `SummarizedExperiment`, 
#' `SingleCellExperiment` and `Seurat` are provided.
#'
#' @param x object to calculate RVP for.
#' @param ... additional arguments to pass to S3 methods.
#'
#' @returns List containing the following components:
#'   \describe{
#'     \item{`RVP`}{the proportion of variance associated with batch effects.}
#'     \item{`sum.squares`}{matrix of sum of squares between batch and total
#'       sum of squares for all features.}
#'     \item{`p.value`}{p-value of permutation test}
#'     \item{`null.distribution`}{numeric, null distribution of RVP values.}
#'   }
#'   Last two components are only present if permuation test is performed.
#'
#' @author Wei Xin Chan
#'
#' @rdname RVP
#' @export
#'
RVP <- function(x, ...) UseMethod("RVP", x)


#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#' @param nperm numeric indicating number of permutations to simulate in the
#'   Monte Carlo permutation test. We recommend a value no less than 1000.
#'   By default, no permutation test is performed.
#' @param use.sparse logical indicating whether to use sparse matrices when
#'   computing RVP. N.B. Using sparse matrices may lead to slight increase
#'  in run time. 
#'
#' @details Default S3 method is for class data frame or matrix with
#'   dimensions (nfeatures, nsamples).
#'
#' @importFrom progress progress_bar
#' @importFrom utils capture.output
#'
#' @rdname RVP
#' @export
#'
RVP.default <- function(
  x, batch, cls = NULL,
  nperm = 0, use.sparse = FALSE,
  ...
) {
  if (!use.sparse) {
    res <- .RVP(x, batch, cls)
  } else {
    res <- .RVP_sparseMatrix(x, batch, cls)
  }

  # Permutation test
  if (nperm > 0L) {
    # TODO: Multiprocessing for permtest
    if (nperm < 100)
      stop("nperm has to be above 100!")

    pb <- progress_bar$new(
      format = "Permutations: [:bar] :current/:total in :elapsed.",
      total = nperm, clear = FALSE, width = 75
    )
    capture.output(pb$tick(0), file = nullfile())

    null_distr <- numeric()
    for (i in seq_len(nperm)) {
      shuffled_batch <- sample(batch)
      if (!use.sparse) {
        null_rvp <- .RVP(x, shuffled_batch, cls)$RVP
      } else {
        null_rvp <- .RVP_sparseMatrix(x, shuffled_batch, cls)$RVP
      }
      null_distr <- c(null_distr, null_rvp)
      pb$tick()
    }
    res$p.value <- sum(null_distr > res$RVP) / nperm
    res$null.distribution <- null_distr
  }
  res
}


#' @param batchname character, name of column in metadata indicating batch.
#' @param classname character, name of column/s in metadata indicating class. 
#'
#' @rdname RVP
#' @export
#'
RVP.Seurat <- function(
  x, batchname, classname = NULL,
  nperm = 0, use.sparse = FALSE,
  ...
) {
  # Suggests: SeuratObject
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("Please install SeuratObject package!")

  # GetAssayData is for Seurat assay v3/v4
  assay_data <- SeuratObject::GetAssayData(x)
  batch <- x@meta.data[[batchname]]
  cls <- if (is.null(classname) || is.na(classname)) {
    NULL
  } else if (length(classname) == 1L) {
    # x[[name]] returns dataframe instead of vector!
    x@meta.data[[classname]]
  } else {
    lapply(classname, function(name) x@meta.data[[name]])
  }
  RVP.default(assay_data, batch, cls, nperm, use.sparse, ...)
}


#' @param assayname character, name of assay to use. By default the first
#'   assay is used.
#'
#' @details S3 method for `SummarizedExperiment` is applicable for the
#'   `SingleExperiment` class as well, as it inherits from the
#'   `SummarizedExperiment` class.
#'
#' @rdname RVP
#' @export
#'
RVP.SummarizedExperiment <- function(
  x, batchname, classname = NULL,
  assayname = NULL,
  nperm = 0, use.sparse = FALSE,
  ...
) {
  # Suggests: SummarizedExperiment 
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Please install SummarizedExperiment package!")

  assay_data <- if (is.null(assayname) || is.na(assayname)) {
    SummarizedExperiment::assay(x)
  } else {
    SummarizedExperiment::assay(x, assayname)
  }
  batch <- x[[batchname]]
  cls <- if (is.null(classname) || is.na(classname)) {
    NULL
  } else if (length(classname) == 1L) {
    x[[classname]]
  } else {
    lapply(classname, function(name) x[[name]])
  }
  RVP.default(assay_data, batch, cls, nperm, use.sparse, ...)
}


#' Helper function that calculates RVP
#'
#' @param X data frame or matrix with dim (nfeatures, nsamples).
#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#'
#' @keywords internal
#' @noRd
#'
.RVP <- function(X, batch, cls) {
  ### CHECK ARGS ###
  if (ncol(X) != length(batch))
    stop("Length of batch does not match number of columns in X!")

  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0

  ### COMPUTE RVP ###
  if (is.null(cls)) {
    feature_names <- rownames(X)
    feature_means <- rowMeans(X)
    ss_total <- rowSums((X - feature_means) ^ 2)
    X_batches <- split_cols(X, batch, drop = TRUE)
    batch_means <- lapply(X_batches, rowMeans)
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(
      (batch_means - feature_means) ^ 2 %*% diag(nperbatches)
    )
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- feature_names 
    colnames(SS) <- c("ss_batch", "ss_total")

    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_names <- rownames(X)
    ss_total <- rowSums((X - rowMeans(X)) ^ 2)
    X_classes <- split_cols(X, cls, drop = TRUE)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: Recursive call
    SS_classes <- lapply(
      mapply(
        .RVP, X_classes, batch_classes,
        MoreArgs = list(cls = NULL),
        SIMPLIFY = FALSE
      ),
      function(obj) obj$sum.squares
    )
    # Filters out obj$sum.squares == NA
    SS_classes <- SS_classes[!is.na(SS_classes)]
    if (length(SS_classes) == 0L) {
      confound_message <- paste(
        "RVP is unable to quantify batch effects as batch and class",
        "are completely confounded!"
      )
      stop(confound_message)
    }
    ss_batch_classes <- lapply(SS_classes, function(X) X[, "ss_batch"])
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- feature_names
    colnames(SS) <- c("ss_batch", "ss_total")

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}


#' Helper function that calculates RVP using sparse matrices
#'
#' @param X data frame or matrix with dim (nfeatures, nsamples).
#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @keywords internal
#' @noRd
#'
.RVP_sparseMatrix <- function(X, batch, cls = NULL) {
  ### CHECK ARGS ###
  if (ncol(X) != length(batch))
    stop("Length of batch does not match number of columns in X!")
  
  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0
  
  ### COMPUTE RVP ###
  if (is.null(cls)) {
    feature_names <- rownames(X)
    feature_means <- rowMeans(X, sparseResult = TRUE)
    ss_total <- rowSums((X - feature_means) ^ 2, sparseResult = TRUE)
    X_batches <- split_cols(X, batch, drop = TRUE)
    batch_means <- lapply(
      X_batches,
      function(X) as(rowMeans(X, sparseResult = TRUE), "sparseMatrix")
    )
    # N.B. cbind only works on sparseMatrix and not sparseVector
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(
      (batch_means - feature_means) ^ 2 %*% .sparseDiagonal(x = nperbatches),
      sparseResult = TRUE
    )
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(feature_names, c("ss_batch", "ss_total")),
      sparse = TRUE
    )
    # print("RVP.sparseMatrix: vars")
    # print(sapply(ls(), function(x) object_size(mget(x, inherits = TRUE))))
    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_names <- rownames(X)
    ss_total <- rowSums(
      (X - rowMeans(X, sparseResult = TRUE)) ^ 2,
      sparseResult = TRUE
    )
    X_classes <- split_cols(X, cls, drop = TRUE)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: Recursive call
    SS_classes <- lapply(
      mapply(
        .RVP_sparseMatrix, X_classes, batch_classes,
        MoreArgs = list(cls = NULL),
        SIMPLIFY = FALSE
      ),
      function(obj) obj$sum.squares
    )
    # Filters out obj$sum.squares == NA
    SS_classes <- SS_classes[!is.na(SS_classes)]
    if (length(SS_classes) == 0L) {
      confound_message <- paste(
        "RVP is unable to quantify batch effects as batch and class",
        "are completely confounded!"
      )
      stop(confound_message)
    }
    ss_batch_classes <- lapply(
      SS_classes,
      # N.B. drop = False to return sparseMatrix instead of dense vector
      function(X) X[, "ss_batch", drop = FALSE]
    )
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    # TODO: Check if there is a problem with ss_batch being a sparseMatrix 
    # and ss_total being a sparseVector
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(feature_names, c("ss_batch", "ss_total")),
      sparse = TRUE
    )

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}
