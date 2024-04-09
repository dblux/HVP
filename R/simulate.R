#' Simulate log-transformed gene expression microarray data
#'
#' @param m number of features.
#' @param crosstab matrix of contingency table specifying number of samples in
#'   each class-batch condition, with classes as rows and batches as columns.
#' @param delta magnitude of additive batch effects.
#' @param gamma magnitude of multiplicative batch effects. Variance parameter
#'   of Normal distribution modelling batch effects across samples in a batch.
#' @param phi percentage of differentially expressed features.
#' @param c shape parameter of Gamma distribution modelling log fold-change.
#' @param d rate parameter of Gamma distribution modelling log fold-change.
#' @param epsilon magnitude of feature-wise variation.
#' @param kappa magnitude of sample-specific variation (scaling factor).
#' @param a shape parameter of Gamma distribution modelling basal expression.
#' @param b rate parameter of Gamma distribution modelling basal expression.
#' @param dropout logical indicating whether to perform dropout
#' @param r inverse scale parameter of the sigmoid function.
#' @param s midpoint parameter of the sigmoid function.
#' @param seed numeric specifying random seed. Defaults to no seed.
#'
#' @returns A list containing the following components:
#'   \describe{
#'     \item{`X`}{matrix with dimensions `(m, n)` of log expression values.}
#'     \item{`metadata`}{data frame with `n` rows of sample metadata.}
#'     \item{`diff.features`}{character vector, names of differentially
#'       expressed features.}
#'     \item{`W`}{matrix with dimensions `(m, n)` of log expression values
#'       without addition of batch effect term.}
#'     \item{`batch.terms`}{matrix with dimensions (m, n) of batch effect terms.}
#'     \item{`class.logfc`}{matrix of mean log fold change for each feature
#'       for classes.}
#'     \item{`batch.logfc`}{matrix of mean log fold change for each features
#'       for batches.}
#'     \item{`params`}{list of parameters supplied.}
#'   }
#'
#' @author Wei Xin Chan
#'
#' @export
#'
simulate_microarray <- function(
  m,
  crosstab,
  delta = 1,
  gamma = 0.5,
  phi = 0.2,
  c = 10,
  d = 6,
  epsilon = 0.5, # limit = (, 1)
  kappa = 0.2, # limit = (, 0.3)
  a = 40,
  b = 5,
  dropout = FALSE,
  r = 2,
  s = -6,
  seed = NA
) {
  params <- list(
    crosstab = crosstab,
    delta = delta, gamma = gamma,
    phi = phi, c = c, d = d,
    epsilon = epsilon, kappa = kappa,
    a = a, b = b, 
    dropout = dropout, r = r, s = s,
    seed = seed
  )

  if (!is.na(seed))
    set.seed(seed)

  n <- sum(crosstab)
  n_class <- nrow(crosstab)
  n_batch <- ncol(crosstab)
  gs <- rep(rep(seq_len(n_class), n_batch), crosstab) # class encoding
  ks <- rep(rep(seq_len(n_batch), each = n_class), crosstab) # batch encoding

  # Metadata
  gs_alphabet <- LETTERS[gs]
  sid <- paste(paste0("ID", seq_len(n)), gs_alphabet, ks, sep = "_")
  metadata <- data.frame(
    class = gs_alphabet,
    batch = as.factor(ks),
    row.names = sid
  )

  log_psi <- rgamma(m, a, rate = b)
  # Log fold-change factors for each class
  # Class A has zero log fold change w.r.t. itself
  log_rho <- matrix(0, m, n_class)
  colnames(log_rho) <- LETTERS[seq_len(n_class)]
  diff.features <- NULL
  if (n_class > 1) {
    n_upreg <- n_downreg <- round(phi * m / 2, 0)
    n_diffexpr <- n_upreg + n_downreg
    for (g in seq(2, n_class)) {
      diff.features <- sort(sample(seq_len(m), n_diffexpr))
      upreg.features <- sort(sample(diff.features, n_upreg))
      downreg.features <- setdiff(diff.features, upreg.features)
      for (i in upreg.features) {
        log_rho[i, g] <- rgamma(1, c, rate = d)
      }
      for (i in downreg.features) {
        log_rho[i, g] <- -rgamma(1, c, rate = d)
      }
    }
  }

  # Base expression values with class effects
  Z <- matrix(0, m, n)
  colnames(Z) <- sid
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      g <- gs[j]
      Z[i, j] <- rnorm(1, log_psi[i] + log_rho[i, g], epsilon)
    }
  }

  # Sample specific scaling term (in log space)
  log_alpha <- rnorm(n, 0, kappa)
  W <- sweep(Z, 2, log_alpha, `+`)

  # Batch effects
  log_beta <- matrix(rnorm(m * n_batch, 0, delta), m, n_batch)
  omega <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      omega[i, j] <- rnorm(1, log_beta[i, k], gamma)
    }
  }

  X <- W + omega
  X[X < 0] <- 0 # set negative values to zero

  if (dropout) {
    P <- sigmoid(X, c, d)
    indicator <- matrix(0, m, n)
    for (i in seq_len(m)) {
      for (j in seq_len(n)) {
        indicator[i, j] <- rbinom(1, 1, P[i, j])
      }
    }
    X <- X * indicator
  }

  list(
    X = X, metadata = metadata,
    diff.features = diff.features,
    W = W, batch.terms = omega,
    class.logfc = log_rho, batch.logfc = log_beta,
    params = params
  )
}


#' Sigmoid function
#'
#' @param numeric scalar/vector/matrix
#' @param r inverse scale parameter of the sigmoid function
#' @param s midpoint parameter of the sigmoid function
#'
#' @returns A numeric scalar/vector/matrix of the same dimensions containing
#'   the transformed values.
#'
#' @export
#'
sigmoid <- function(x, r = 1, s = 0) 1 / (1 + exp(-(r * (x + s))))
