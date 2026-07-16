# Random technology with the materials-balance identity satisfied by
# construction: b = share * (u'x - v y), share in (0, 1).
make_random_tech <- function(L = 30, N = 3, groups = c("A", "B"),
                             v = 0.01, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(L * N, 10, 100), L, N,
              dimnames = list(NULL, paste0("in", seq_len(N))))
  u <- rep(1, N)
  potential <- as.vector(x %*% u)
  y <- runif(L, 1, 0.5 * min(potential) / max(v, 1e-12))
  y <- pmin(y, 0.5 * potential / max(v, 1e-12))
  b <- (potential - v * y) * runif(L, 0.6, 0.98)
  group <- factor(sample(groups, L, replace = TRUE))
  pgt_tech(x = x, y = y, b = b, u = u, v = v, group = group)
}
