set.seed(0)

########## SETUP ########## 

# Data
n <- 200
m <- 1000
sids <-  sprintf("SID%03d", seq_len(n)) # n is 3 digits
x <- rgamma(m * n, 40, 5)
missing <- rbinom(m * n, 1, 0.8) == 1
x[missing] <- 0
dense <- matrix(
  x, m, n,
  dimnames = list(seq_len(m), sids)
)
sparse <- as(dense, "sparseMatrix")
data <- data.frame(dense)

# Metadata
ncond <- n / 4
batch <- rep(1:2, each = ncond * 2)
class1 <- rep(rep(LETTERS[1:2], each = ncond), 2)
class1_int <- ifelse(class1 == "A", 1, 2)
class2 <- rep(LETTERS[3:4], ncond * 2)
class3 <- rep(LETTERS[3:4], each = ncond * 2)
metadata <- data.frame(batch, class1, class2, row.names = colnames(dense))

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(logcounts = sparse),
  colData = metadata
)
seu <- SeuratObject::CreateSeuratObject(
  counts = sparse,
  meta.data = metadata
)

########## UNIT TESTS ##########

# Test: Utility functions
test_that("split_cols works when arg: f is a list of vectors", {
  list_dense <- split_cols(dense, list(class1, class2), drop = TRUE)

  expect_named(list_dense, c("A.C", "B.C", "A.D", "B.D"))
  expect_identical(
    unname(sapply(list_dense, nrow)),
    rep(nrow(dense), 4)
  )
})

# Test: HVP
HVP_TRUE <- 0.01011387
test_that("HVP works on a dense matrix", {
  # test that arg: batch, cls can be of type {numeric, character, factor}
  expect_equal(HVP(dense, batch, class1)$HVP, HVP_TRUE)
  expect_equal(HVP(dense, as.factor(batch), class1)$HVP, HVP_TRUE)
  expect_equal(HVP(dense, as.character(batch), class1)$HVP, HVP_TRUE)
  expect_equal(HVP(dense, batch, as.factor(class1))$HVP, HVP_TRUE)
  expect_equal(HVP(dense, batch, class1_int)$HVP, HVP_TRUE)
})

test_that("HVP works when arg: f is a list of vectors", {
  expect_no_error(HVP(dense, batch, list(class1, class2)))
})

test_that("HVP throws error when batch and class are completely confounded", {
  expect_error(
    HVP(dense, batch, list(class1, class3)),
    "batch and class are completely confounded"
  )
})

test_that("HVP works on a sparse matrix", {
  expect_equal(HVP(sparse, batch, class1)$HVP, HVP_TRUE)
})

test_that("HVP works on a data frame", {
  expect_equal(HVP(data, batch, class1)$HVP, HVP_TRUE)
})

test_that("HVP works on a SingleCellExperiment object", {
  expect_equal(HVP(sce, "batch", "class1")$HVP, HVP_TRUE)
})

test_that("HVP works on a Seurat object", {
  expect_equal(HVP(seu, "batch", "class1")$HVP, HVP_TRUE)
})

PVALUE <- 0.404
cat("Running permutation test...", fill = TRUE)
test_that("Permutation test returns p.value and null.distribution", {
  expect_no_error({
    obj <- HVP(dense, batch, class1, nperm = 1000)
    expect_equal(obj$p.value, PVALUE)
  })
})

test_that("Helper HVP function that uses sparse matrices works", {
  expect_equal(HVP(sparse, batch, class1, use.sparse = TRUE)$HVP, HVP_TRUE)
})
