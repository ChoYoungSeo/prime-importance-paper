test_that("tied projected neighbors exclude self without changing valid entries", {
  neighbors <- getFromNamespace(".prime_query_neighbors", "prime")
  cases <- list(matrix(rep(c(0,1,3,7),each=3),ncol=1),
                matrix(rep(rep(c(0,1,3,7),each=3),5),ncol=5),
                matrix(0,20,1))
  for (X in cases) for (k in c(1L,4L)) {
    n <- nrow(X)
    algorithm <- if(n<=2000 && ncol(X)>=5) "brute" else "kd_tree"
    original <- FNN::get.knn(X,k=k,algorithm=algorithm)$nn.index
    fixed <- neighbors(X,seq_len(n),k)
    self <- original==row(original)
    expect_false(any(fixed==row(fixed)))
    expect_identical(fixed[!self],original[!self])
    expect_true(all(apply(fixed,1,function(r)!anyDuplicated(r))))
    for (i in seq_len(n)) {
      distances <- rowSums(sweep(X,2,X[i,],"-")^2)
      distances[i] <- Inf
      expect_lte(max(distances[fixed[i,]]),sort(distances)[k]+1e-12)
    }
  }
})

test_that("untied full-query and existing partial-query behavior are preserved", {
  neighbors <- getFromNamespace(".prime_query_neighbors", "prime")
  set.seed(918)
  X <- matrix(runif(200),40,5)
  original <- FNN::get.knn(X,k=8,algorithm="brute")$nn.index
  expect_identical(neighbors(X,seq_len(40),8L),original)
  q <- c(2L,8L,19L)
  partial <- neighbors(X,q,8L)
  expect_false(any(partial==matrix(q,nrow(partial),ncol(partial))))
  expect_identical(partial,original[q,,drop=FALSE])
})
