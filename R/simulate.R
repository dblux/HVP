#' Simulate log-transformed bulk gene expression data
#'
#' @param crosstab matrix of contingency table specifying number of samples in
#'   each class-batch condition, with classes as rows and batches as columns.
#' @param m number of genes.
#' @param delta magnitude of additive batch effects (i.e. standard deviation
#'   of normal distribution modelling batch log fold change means of all genes).
#' @param gamma magnitude of multiplicative batch effects (i.e. standard
#'   deviation of normal distribution modelling log batch effect terms of all
#'   samples in a batch).
#' @param phi percentage of differentially expressed genes.
#' @param c shape parameter of Gamma distribution modelling class log fold
#'   change means of all genes.
#' @param d rate parameter of Gamma distribution modelling class log fold
#'   change means of all genes.
#' @param epsilon magnitude of random noise across samples (i.e. standard
#'   deviation of normal distribution modelling log expression values with
#'   class effects only).
#' @param kappa standard deviation of normal distribution modelling log scaling
#'   factors of all samples.
#' @param a shape parameter of Gamma distribution modelling basal log mean
#'   expression of all genes.
#' @param b rate parameter of Gamma distribution modelling basal log mean
#'   expression of all genes.
#' @param dropout logical indicating whether to perform dropout
#' @param r inverse scale parameter of the sigmoid function used to calculate
#'    probability of dropout for each value.
#' @param s midpoint parameter of the sigmoid function used to calculate
#'    probability of dropout for each value.
#'
#' @returns A list containing the following components:
#'   \describe{
#'     \item{`X`}{matrix with dimensions `(m, n)` of log expression values.}
#'     \item{`metadata`}{data frame with `n` rows of sample metadata.}
#'     \item{`diff.genes`}{
#'       character vector, names of differentially expressed genes.}
#'     \item{`Y`}{
#'       matrix with dimensions `(m, n)` of log expression values with class
#'       effects only.}
#'     \item{`batch.terms`}{
#'       matrix with dimensions `(m, n)` of log batch effect terms.}
#'     \item{`class.logfc`}{
#'       matrix of class log fold change means for each gene in each class}
#'     \item{`batch.logfc`}{
#'       matrix of batch log fold change means for each gene in each batch.}
#'     \item{`params`}{list of parameters supplied.}
#'   }
#'
#' @author Wei Xin Chan
#'
#' @importFrom stats rbinom rgamma rnorm
#' @importFrom utils capture.output
#'
#' @examples
#'
#' crosstab <- matrix(10, 3, 2)
#' data <- simulate_bulkexpr(crosstab, 100)
#'
#' @export
#'
simulate_bulkexpr <- function(
  crosstab,
  m,
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
  s = -6
) {
  params <- list(
    crosstab = crosstab,
    delta = delta, gamma = gamma,
    phi = phi, c = c, d = d,
    epsilon = epsilon, kappa = kappa,
    a = a, b = b, 
    dropout = dropout, r = r, s = s
  )
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

  psi <- rgamma(m, a, rate = b)
  # Log fold-change factors for each class
  # Class A has zero log fold change w.r.t. itself
  rho <- matrix(0, m, n_class)
  colnames(rho) <- LETTERS[seq_len(n_class)]
  diff.genes <- NULL
  if (n_class > 1) {
    n_upreg <- n_downreg <- round(phi * m / 2, 0)
    n_diffexpr <- n_upreg + n_downreg
    for (g in seq(2, n_class)) {
      diff.genes <- sort(sample(seq_len(m), n_diffexpr))
      upreg.genes <- sort(sample(diff.genes, n_upreg))
      downreg.genes <- setdiff(diff.genes, upreg.genes)
      for (i in upreg.genes) {
        rho[i, g] <- rgamma(1, c, rate = d)
      }
      for (i in downreg.genes) {
        rho[i, g] <- -rgamma(1, c, rate = d)
      }
    }
  }

  # Base expression values with class effects
  Z <- matrix(0, m, n)
  colnames(Z) <- sid
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      g <- gs[j]
      Z[i, j] <- rnorm(1, psi[i] + rho[i, g], epsilon)
    }
  }

  # Sample specific scaling term (in log space)
  alpha <- rnorm(n, 0, kappa)
  Y <- sweep(Z, 2, alpha, `+`)

  # Batch effects
  beta <- matrix(rnorm(m * n_batch, 0, delta), m, n_batch)
  omega <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      omega[i, j] <- rnorm(1, beta[i, k], gamma)
    }
  }

  X <- Y + omega
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
    X = X,
    metadata = metadata,
    diff.genes = diff.genes,
    Y = Y,
    batch.terms = omega,
    class.logfc = rho,
    batch.logfc = beta,
    params = params
  )
}


#' Sigmoid function
#'
#' @param x numeric scalar/vector/matrix
#' @param r inverse scale parameter of the sigmoid function
#' @param s midpoint parameter of the sigmoid function
#'
#' @returns A numeric scalar/vector/matrix of the same dimensions containing
#'   the transformed values.
#'
#' @examples
#'
#' p <- sigmoid(0.5)
#'
#' @export
#'
sigmoid <- function(x, r = 1, s = 0) 1 / (1 + exp(-(r * (x + s))))
