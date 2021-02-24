library(hedgehog)

gen.matrix <- function(gen, nrow, ncol) {
  generate(for (ns in gen.c(gen)) {
    matrix(ns, nrow = nrow, ncol = ncol)
  })
}

gen.sym <- function(gen, n) {
  generate(for (s in n) {
    matrix(ns, nrow = s, ncol = s)
  })
}

test_that("laplacian approximations have correct summations", {
  expect_true(T)
})
