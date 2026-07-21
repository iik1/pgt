# Internal linear-programming kernels. All solvers evaluate one DMU at a
# time against a reference (peer) set and return the minimum feasible bad
# output b* (or, for the directional model, the maximal joint expansion)
# together with the peer weights and selected constraint duals.
#
# Formulations follow Rodseth (2025, JPA): Eq. 6 in reduced form for the
# weak-G-disposability model and Eq. 13 (abatement fixed at the evaluated
# DMU's own level) for the factorially determined multi-output /
# directional representation. The unit tests replicate the paper's
# numerical example (Tables 2-3) through these kernels.
#
# Material flow coefficients are DMU-specific (one row of the L x N x P
# array per DMU). The objective acts on a single pollutant p per solve;
# its per-DMU potential enters through mb_rhs = u_i'x_i - v_i y_i. In the
# weak-G-disposability model the caps of the remaining pollutants
# constrain the peer mix as well (one row per extra pollutant).
#
# The kernels rebuild the LP per DMU. For a fixed peer set the constraint
# matrices are DMU-invariant (only RHS values and single coefficients
# change), so a set.rhs-based reuse of one LP object is a possible future
# optimisation; DMU-invariant vectors are already hoisted via .solve_ctx().

# Pollutant potential u_i'x_i and materials-balance cap
# u_i'x_i - v_i'y_i for pollutant p (one value per DMU). The cap equals
# the uncontrolled emission z_i when the account closes.
.mb_potential <- function(tech, p = 1L) {
  rowSums(matrix(tech$u[, , p], tech$L, tech$N) * tech$x)
}
.mb_cap <- function(tech, p = 1L) {
  .mb_potential(tech, p) - .retained(tech, p)
}

# Weak-G-disposability estimator (Rodseth 2025, Eq. 6/7) in reduced
# form, for pollutant p. Equation 6 anchors only the intended outputs
# at the evaluated unit (inputs are decision variables), so
# substituting the summing-up condition u_i'eps_x + v_i'eps_y = eps_b
# and setting the input slacks to zero collapses the programme to
#
#   b*_i = min_lambda sum_l lambda_l b_l + v_i'(sum_l lambda_l y_l - y_i)
#   s.t.  sum_l lambda_l y_ml >= y_mi   (m = 1..M)
#         sum_l lambda_l = 1  [VRS only],  lambda >= 0
#
# with the evaluated unit's own coefficients v_i (Eq. 9 form under
# producer-specific coefficients). Self-reference is always feasible,
# so scores b*/b lie in (0, 1] with no feasibility screen. The optimal
# peer mix also delivers the implied uncontrolled emission
# z* = sum lambda z_l, abatement a* = sum lambda a_l (when observed)
# and input point x* = sum lambda x_l, the quantities in the paper's
# Table 2. The reported output duals are the total derivatives
# d b*/d y_mi = mu_m - v_mi, where mu_m is the output-row dual; they
# can be negative when the retained content v is large (producing more
# output binds more pollutant into the product).
.lp_wgd_one <- function(i, Y, b, v_i, peers, vrs = TRUE,
                        z = NULL, a = NULL) {
  L <- length(peers)
  M <- ncol(Y)

  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  penalty <- as.vector(Y[peers, , drop = FALSE] %*% v_i)
  lpSolveAPI::set.objfn(lp, b[peers] + penalty)
  for (m in seq_len(M)) {
    lpSolveAPI::add.constraint(lp, Y[peers, m], ">=", Y[i, m])
  }
  if (vrs) {
    lpSolveAPI::add.constraint(lp, rep(1, L), "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, L))

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, lambda = NULL,
                dual_output = rep(NA_real_, M),
                z_star = NA_real_, a_star = NA_real_))
  }
  lambda <- lpSolveAPI::get.variables(lp)
  duals <- lpSolveAPI::get.sensitivity.rhs(lp)$duals
  list(
    status = status,
    b_star = lpSolveAPI::get.objective(lp) - sum(v_i * Y[i, ]),
    lambda = lambda,
    dual_output = duals[seq_len(M)] - v_i,
    z_star = if (is.null(z)) NA_real_ else sum(lambda * z[peers]),
    a_star = if (is.null(a)) NA_real_ else sum(lambda * a[peers])
  )
}

# The package's input-anchored benchmark (previously model = "wgd"):
# minimal emissions holding the evaluated unit's inputs fixed, with the
# peer emission envelope and the unit's own materials-balance cap. This
# is NOT Eq. 6 of Rodseth (2025), whose inputs are free; it is the
# current-input companion the decomposition's early stages build on.
#
#   min_{lambda, bq} bq
#   s.t. sum_l lambda_l y_ml         >= y_mi  (m = 1..M)
#        sum_l lambda_l x_nl         <= x_ni  (n = 1..N, optional)
#        sum_l lambda_l b_l  -  bq   <= 0
#        bq                          <= u_i'x_i - v_i'y_i   [cap]
#        sum_l lambda_l b_ql         <= u_qi'x_i - v_qi'y_i [extra pollutants]
#        sum_l lambda_l               = 1                   [VRS only]
#        lambda >= 0, bq >= 0
#
# A DMU violating any pollutant's identity loses self-reference
# feasibility here, so its LP solves only if some peer mix meets every
# row within the caps; infeasibility (status != 0) is confined to
# violating DMUs.
.lp_wgd_anchored_one <- function(i, X, Y, b, mb_rhs, peers, vrs = TRUE,
                                 input_constraints = TRUE,
                                 b_other = NULL, cap_other = NULL) {
  L <- length(peers)
  N <- ncol(X)
  M <- ncol(Y)

  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  lpSolveAPI::set.objfn(lp, c(rep(0, L), 1))

  for (m in seq_len(M)) {
    lpSolveAPI::add.constraint(lp, c(Y[peers, m], 0), ">=", Y[i, m])
  }
  if (input_constraints) {
    for (n in seq_len(N)) {
      lpSolveAPI::add.constraint(lp, c(X[peers, n], 0), "<=", X[i, n])
    }
  }
  lpSolveAPI::add.constraint(lp, c(b[peers], -1), "<=", 0)
  lpSolveAPI::add.constraint(lp, c(rep(0, L), 1), "<=", mb_rhs[i])
  if (!is.null(b_other)) {
    for (q in seq_len(ncol(b_other))) {
      lpSolveAPI::add.constraint(lp, c(b_other[peers, q], 0), "<=",
                                 cap_other[i, q])
    }
  }
  if (vrs) {
    lpSolveAPI::add.constraint(lp, c(rep(1, L), 0), "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, L + 1))

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, lambda = NULL,
                dual_output = rep(NA_real_, M), mb_rhs = mb_rhs[i]))
  }
  sol <- lpSolveAPI::get.variables(lp)
  duals <- lpSolveAPI::get.sensitivity.rhs(lp)$duals
  # The bad-output (emission-envelope) dual is structurally -1 whenever
  # the LP solves, so it carries no information and is not returned; the
  # informative MB diagnostic is the cap headroom mb_rhs - b_star.
  list(
    status = status,
    b_star = sol[L + 1],
    lambda = sol[seq_len(L)],
    dual_output = duals[seq_len(M)],
    mb_rhs = mb_rhs[i]
  )
}

# Convex lower (y, b)-envelope LP: Rodseth Eq. 6 with inputs free.
#
#   min_lambda sum_l lambda_l b_l
#   s.t. sum_l lambda_l y_l >= y_i,  sum lambda = 1 [VRS],  lambda >= 0
#
# With inputs free the weak-G slack equality u'eps_x + v eps_y = eps_b is
# absorbed into the reallocation (inputs adjust with the peer mix and
# emissions follow the pollutant they carry), so the program reduces to
# the convex lower envelope of the (y, b) scatter. Self-reference is
# always feasible, so scores b*/b lie in (0, 1] with no screen.
.lp_envelope_one <- function(i, Y, b, peers, vrs = TRUE) {
  L <- length(peers)
  M <- ncol(Y)
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  lpSolveAPI::set.objfn(lp, b[peers])
  for (m in seq_len(M)) {
    lpSolveAPI::add.constraint(lp, Y[peers, m], ">=", Y[i, m])
  }
  if (vrs) {
    lpSolveAPI::add.constraint(lp, rep(1, L), "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, L))

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, lambda = NULL,
                dual_output = rep(NA_real_, M), mb_rhs = NA_real_))
  }
  duals <- lpSolveAPI::get.sensitivity.rhs(lp)$duals
  list(
    status = status,
    b_star = lpSolveAPI::get.objective(lp),
    lambda = lpSolveAPI::get.variables(lp),
    dual_output = duals[seq_len(M)],
    mb_rhs = NA_real_
  )
}

# Factorially determined multi-output / directional model (Rodseth 2025,
# Eq. 13 with abatement fixed at the evaluated DMU's own level a_i, the
# specification that yields the paper's Table 3). Pollutant p.
#
#   max_{lambda, thy, thb} thy + thb
#   s.t. sum_l lambda_l y_l  - thy        >= y_i           (good output)
#        sum_l lambda_l a_l                >= a_i           (abatement held
#                                                            at own level)
#        sum_l lambda_l x_nl              <= x_ni  (n=1..N)
#        v_i thy - thb                     = z_i - b_i - a_i (materials
#                                                            balance)
#        sum_l lambda_l                    = 1             (VRS)
#        lambda, thy, thb >= 0
#
# where z_i = u_i'x_i - v_i y_i is the evaluated DMU's uncontrolled
# emission. Holding abatement at the DMU's own level (its pollution
# control absorbs resources that cannot expand the good output) is what
# keeps high-abatement DMUs on the frontier. thy is the good-output
# expansion and thb the bad-output contraction; gross inefficiency
# thy + thb decomposes into good-output efficiency thy and bad-output
# efficiency thb, and the maximal good output is y_i + thy. When
# z_i = b_i + a_i (uncontrolled = controlled + abatement) the materials
# balance forces thb = v_i thy.
.lp_fdmo_one <- function(i, X, y, a, v, z, b, peers, vrs = TRUE,
                         input_constraints = TRUE) {
  L <- length(peers)
  N <- ncol(X)
  # column order: lambda_1..L, thy, thb
  ncol_lp <- L + 2L
  ithy <- L + 1L
  ithb <- L + 2L

  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = ncol_lp)
  invisible(lpSolveAPI::lp.control(lp, sense = "max"))
  obj <- numeric(ncol_lp)
  obj[ithy] <- 1
  obj[ithb] <- 1
  lpSolveAPI::set.objfn(lp, obj)

  # good-output frontier: sum lambda y - thy >= y_i
  r <- numeric(ncol_lp); r[seq_len(L)] <- y[peers]; r[ithy] <- -1
  lpSolveAPI::add.constraint(lp, r, ">=", y[i])
  # abatement fixed at own level: sum lambda a >= a_i
  r <- numeric(ncol_lp); r[seq_len(L)] <- a[peers]
  lpSolveAPI::add.constraint(lp, r, ">=", a[i])
  # inputs
  if (input_constraints) {
    for (n in seq_len(N)) {
      r <- numeric(ncol_lp); r[seq_len(L)] <- X[peers, n]
      lpSolveAPI::add.constraint(lp, r, "<=", X[i, n])
    }
  }
  # materials balance: v_i thy - thb = z_i - b_i - a_i
  r <- numeric(ncol_lp); r[ithy] <- v[i]; r[ithb] <- -1
  lpSolveAPI::add.constraint(lp, r, "=", z[i] - b[i] - a[i])
  if (vrs) {
    r <- numeric(ncol_lp); r[seq_len(L)] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, ncol_lp))

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, gross = NA_real_, theta_y = NA_real_,
                theta_b = NA_real_, lambda = NULL))
  }
  sol <- lpSolveAPI::get.variables(lp)
  list(
    status = status,
    gross = sol[ithy] + sol[ithb],
    theta_y = sol[ithy],
    theta_b = sol[ithb],
    lambda = sol[seq_len(L)]
  )
}

# By-production intersection technology (Murty, Russell and Levkoff
# 2012). The technology is the intersection of two independent
# sub-technologies, solved as two LPs:
#
#   T1 (intended production): the good output is producible from the
#   inputs with free disposal. Output efficiency is the reciprocal of
#   the maximal radial expansion of y,
#     E1 = 1 / max{ phi : sum_l lambda_l y_l >= phi y_i,
#                         sum_l lambda_l x_nl <= x_ni,
#                         [sum lambda = 1], lambda >= 0 }.
#
#   T2 (residual generation): the bad output is caused by the
#   emission-generating inputs and is costly to dispose (bounded below).
#   Emission efficiency is the maximal radial contraction of b,
#     E2 = min{ psi : sum_l mu_l b_l <= psi b_i,
#                     sum_l mu_l x_nl >= x_ni  (n in polluting),
#                     [sum mu = 1], mu >= 0 }.
#
# The arithmetic mean of the two sub-efficiencies (in the spirit of the
# Fare-Grosskopf-Lovell graph measure) is
# E_FGL = (E1 + E2) / 2. E2 = b*/b is the environmental efficiency,
# reported as the principal score; E1 and E_FGL are returned alongside.
.lp_byprod_one <- function(i, X, Y, b, pol, peers, vrs = TRUE) {
  L <- length(peers)
  N <- ncol(X)

  # T1: output efficiency (max phi).
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp1, sense = "max"))
  o <- numeric(L + 1); o[L + 1] <- 1; lpSolveAPI::set.objfn(lp1, o)
  for (m in seq_len(ncol(Y))) {
    r <- c(Y[peers, m], -Y[i, m])
    lpSolveAPI::add.constraint(lp1, r, ">=", 0)
  }
  for (n in seq_len(N)) {
    lpSolveAPI::add.constraint(lp1, c(X[peers, n], 0), "<=", X[i, n])
  }
  if (vrs) lpSolveAPI::add.constraint(lp1, c(rep(1, L), 0), "=", 1)
  lpSolveAPI::set.bounds(lp1, lower = rep(0, L + 1))
  st1 <- lpSolveAPI::solve.lpExtPtr(lp1)
  phi <- if (st1 == 0) lpSolveAPI::get.variables(lp1)[L + 1] else NA_real_
  lam <- if (st1 == 0) lpSolveAPI::get.variables(lp1)[seq_len(L)] else NULL

  # T2: emission efficiency (min psi).
  lp2 <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp2, sense = "min"))
  o <- numeric(L + 1); o[L + 1] <- 1; lpSolveAPI::set.objfn(lp2, o)
  r <- c(b[peers], -b[i]); lpSolveAPI::add.constraint(lp2, r, "<=", 0)
  for (n in pol) {
    lpSolveAPI::add.constraint(lp2, c(X[peers, n], 0), ">=", X[i, n])
  }
  if (vrs) lpSolveAPI::add.constraint(lp2, c(rep(1, L), 0), "=", 1)
  lpSolveAPI::set.bounds(lp2, lower = rep(0, L + 1))
  st2 <- lpSolveAPI::solve.lpExtPtr(lp2)
  psi <- if (st2 == 0) lpSolveAPI::get.variables(lp2)[L + 1] else NA_real_
  mu <- if (st2 == 0) lpSolveAPI::get.variables(lp2)[seq_len(L)] else NULL

  # A single failed sub-LP invalidates the whole intersection measure:
  # return every score as NA so status != 0 always means "scores NA".
  status <- if (st1 == 0 && st2 == 0) 0L else max(st1, st2)
  if (status != 0L) {
    return(list(status = status, output_eff = NA_real_,
                emission_eff = NA_real_, fgl = NA_real_,
                b_star = NA_real_, lambda = NULL, mu = NULL))
  }
  e1 <- if (phi <= 0) NA_real_ else 1 / phi
  list(status = status, output_eff = e1, emission_eff = psi,
       fgl = (e1 + psi) / 2, b_star = psi * b[i],
       lambda = lam, mu = mu)
}

# Materials-balance cost model (Coelli, Lauwers and Van Huylenbroeck
# 2007). Environmental efficiency is the ratio of minimal to observed
# aggregate material inflow u'x, decomposed into technical and
# environmental-allocative efficiency, EE = TE x EAE, with the material
# flow coefficients u_i playing the role of prices:
#
#   TE = min{ theta : sum_l lambda_l x_nl <= theta x_ni,
#                     sum_l lambda_l y_l >= y_i, [sum lambda = 1] }
#   EE = min{ u_i'xproj : sum_l lambda_l x_nl <= xproj_n,
#                         sum_l lambda_l y_l >= y_i, [sum lambda = 1] }
#            / (u_i'x_i)
#   EAE = EE / TE.
#
# TE holds the input mix fixed (radial); EE lets the input mix move
# toward low-pollutant inputs, so EE <= TE and EAE lies in (0, 1]. The
# implied minimal controlled emission is u_i'xproj* - v_i y_i.
.lp_mbcost_one <- function(i, X, Y, u_row, ret_i, peers, vrs = TRUE) {
  L <- length(peers)
  N <- ncol(X)
  pot_i <- sum(u_row * X[i, ])
  # No pollutant potential (the DMU uses none of the pollutant-bearing
  # inputs): the material-inflow ratio is undefined.
  if (pot_i <= .Machine$double.eps) {
    return(list(status = 2L, mbe = NA_real_, te = NA_real_,
                eae = NA_real_, b_star = NA_real_, lambda = NULL))
  }

  # radial technical efficiency (min theta)
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp1, sense = "min"))
  o <- numeric(L + 1); o[L + 1] <- 1; lpSolveAPI::set.objfn(lp1, o)
  for (n in seq_len(N)) {
    r <- c(X[peers, n], -X[i, n]); lpSolveAPI::add.constraint(lp1, r, "<=", 0)
  }
  for (m in seq_len(ncol(Y))) {
    lpSolveAPI::add.constraint(lp1, c(Y[peers, m], 0), ">=", Y[i, m])
  }
  if (vrs) lpSolveAPI::add.constraint(lp1, c(rep(1, L), 0), "=", 1)
  lpSolveAPI::set.bounds(lp1, lower = c(rep(0, L), 0))
  st1 <- lpSolveAPI::solve.lpExtPtr(lp1)
  te <- if (st1 == 0) lpSolveAPI::get.variables(lp1)[L + 1] else NA_real_

  # material-inflow minimisation (min u_i'xproj)
  nc <- L + N
  lp2 <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp2, sense = "min"))
  o <- numeric(nc); o[(L + 1):nc] <- u_row; lpSolveAPI::set.objfn(lp2, o)
  for (m in seq_len(ncol(Y))) {
    lpSolveAPI::add.constraint(lp2, c(Y[peers, m], rep(0, N)), ">=",
                               Y[i, m])
  }
  for (n in seq_len(N)) {
    r <- numeric(nc); r[seq_len(L)] <- X[peers, n]; r[L + n] <- -1
    lpSolveAPI::add.constraint(lp2, r, "<=", 0)
  }
  if (vrs) {
    lpSolveAPI::add.constraint(lp2, c(rep(1, L), rep(0, N)), "=", 1)
  }
  lpSolveAPI::set.bounds(lp2, lower = rep(0, nc))
  st2 <- lpSolveAPI::solve.lpExtPtr(lp2)
  pot_star <- if (st2 == 0) lpSolveAPI::get.objective(lp2) else NA_real_
  lam <- if (st2 == 0) lpSolveAPI::get.variables(lp2)[seq_len(L)] else NULL

  # As in .lp_byprod_one: one failed sub-LP invalidates the EE = TE x EAE
  # decomposition, so status != 0 always means every score is NA.
  status <- if (st1 == 0 && st2 == 0) 0L else max(st1, st2)
  if (status != 0L) {
    return(list(status = status, mbe = NA_real_, te = NA_real_,
                eae = NA_real_, b_star = NA_real_, lambda = NULL))
  }
  ee <- pot_star / pot_i
  list(status = status, mbe = ee, te = te, eae = ee / te,
       b_star = pot_star - ret_i, lambda = lam)
}

# Weak-disposability model (Kuosmanen 2005 correct VRS formulation),
# a reference axiom system for cross-checking. The technology splits the
# intensity weights into an active part z (which scales the good and bad
# outputs together) and an abatement part w (which scales inputs only):
#
#   min phi
#   s.t. sum_l z_l y_l           >= y_i
#        sum_l z_l b_l            =  phi b_i        (weak disposability)
#        sum_l (z_l + w_l) x_nl   <= x_ni  (n=1..N)
#        sum_l (z_l + w_l)        =  1             (VRS)
#        z, w >= 0
#
# phi = b*/b is the emission efficiency: the bad output cannot be
# reduced without proportionally scaling down the good output. Weak
# disposability, unlike weak-G-disposability, imposes no materials
# balance, so this model is included only as a reference cross-check.
.lp_wd_one <- function(i, X, Y, b, peers, vrs = TRUE) {
  L <- length(peers)
  N <- ncol(X)
  # columns: z_1..L, w_1..L, phi
  nc <- 2L * L + 1L
  iphi <- nc
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  o <- numeric(nc); o[iphi] <- 1; lpSolveAPI::set.objfn(lp, o)

  # good outputs: sum z y_m >= y_mi
  for (m in seq_len(ncol(Y))) {
    r <- numeric(nc); r[seq_len(L)] <- Y[peers, m]
    lpSolveAPI::add.constraint(lp, r, ">=", Y[i, m])
  }
  # bad output: sum z b - phi b_i = 0
  r <- numeric(nc); r[seq_len(L)] <- b[peers]; r[iphi] <- -b[i]
  lpSolveAPI::add.constraint(lp, r, "=", 0)
  # inputs: sum (z + w) x <= x_i
  for (n in seq_len(N)) {
    r <- numeric(nc); r[seq_len(L)] <- X[peers, n]
    r[(L + 1):(2L * L)] <- X[peers, n]
    lpSolveAPI::add.constraint(lp, r, "<=", X[i, n])
  }
  if (vrs) {
    r <- numeric(nc); r[seq_len(2L * L)] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, nc))
  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, lambda = NULL,
                dual_output = NA_real_, mb_rhs = NA_real_))
  }
  phi <- lpSolveAPI::get.variables(lp)[iphi]
  z <- lpSolveAPI::get.variables(lp)[seq_len(L)]
  list(status = status, b_star = phi * b[i], lambda = z,
       dual_output = NA_real_, mb_rhs = NA_real_)
}

# DMU-invariant quantities of a per-DMU solve loop, computed once per
# fit and passed to every .lp_solve_one() call (they would otherwise be
# recomputed L times per fit and L x B times in boot_pgt()).
.solve_ctx <- function(tech, model, p) {
  ctx <- list(b_p = tech$b[, p])
  if (model == "wgd") {
    ctx$v_p <- matrix(tech$v[, , p], tech$L, tech$M)
    ctx$z <- .mb_cap(tech, p)
    ctx$a_p <- if (is.null(tech$a)) NULL else tech$a[, p]
  }
  if (model %in% c("wgd_anchored", "fdmo")) {
    ctx$mb_cap <- .mb_cap(tech, p)
  }
  if (model == "wgd_anchored" && tech$P > 1L) {
    q <- setdiff(seq_len(tech$P), p)
    ctx$b_other <- tech$b[, q, drop = FALSE]
    ctx$cap_other <- vapply(q, function(pp) .mb_cap(tech, pp),
                            numeric(tech$L))
    if (!is.matrix(ctx$cap_other)) {
      ctx$cap_other <- matrix(ctx$cap_other, ncol = length(q))
    }
  }
  if (model == "mb_cost") {
    ctx$ret <- .retained(tech, p)
  }
  if (model == "byprod") {
    ctx$pol <- .polluting_inputs(tech, p)
  }
  ctx
}

# Dispatch a single-DMU solve for the requested model and pollutant.
.lp_solve_one <- function(model, i, tech, peers, vrs, p = 1L,
                          input_constraints = TRUE, ctx = NULL) {
  if (is.null(ctx)) ctx <- .solve_ctx(tech, model, p)
  switch(
    model,
    wgd = .lp_wgd_one(i, tech$y, ctx$b_p, ctx$v_p[i, ], peers,
                      vrs = vrs, z = ctx$z, a = ctx$a_p),
    wgd_anchored = .lp_wgd_anchored_one(
      i, tech$x, tech$y, ctx$b_p, ctx$mb_cap, peers, vrs = vrs,
      input_constraints = input_constraints,
      b_other = ctx$b_other, cap_other = ctx$cap_other),
    envelope = .lp_envelope_one(i, tech$y, ctx$b_p, peers, vrs = vrs),
    fdmo = {
      if (is.null(tech$a)) {
        stop("model = \"fdmo\" requires an abatement output 'a' in ",
             "pgt_tech().", call. = FALSE)
      }
      .lp_fdmo_one(i, tech$x, tech$y[, 1L], tech$a[, p],
                   as.vector(tech$v[, 1L, p]),
                   ctx$mb_cap, ctx$b_p,
                   peers, vrs = vrs, input_constraints = input_constraints)
    },
    byprod = .lp_byprod_one(i, tech$x, tech$y, ctx$b_p,
                            ctx$pol, peers, vrs = vrs),
    mb_cost = .lp_mbcost_one(i, tech$x, tech$y, tech$u[i, , p],
                             ctx$ret[i], peers, vrs = vrs),
    wd = .lp_wd_one(i, tech$x, tech$y, ctx$b_p, peers, vrs = vrs),
    stop("unknown model '", model, "'", call. = FALSE)
  )
}

# Stage programme of the Rodseth (2025) Eq. 11 decomposition, built on
# the extended representation of Eq. 9. Minimises controlled emissions
# z - a for pollutant p while holding, per the stage, the production
# inputs (hold_xp), the pollution-control inputs (hold_xa), the
# abatement output at its own level (hold_a) and the quality variable
# rho at zero (hold_quality). Variables: lambda (L), z, a, eps_z,
# eps_a, rho, one output slack per good output, and one input slack per
# held input. The quality row defines rho as the peer mix's aggregate
# coefficient gap relative to the evaluated unit (Eq. 9); holding
# quality pins rho to zero. The summing-up condition runs over the held
# inputs and the output slacks with the evaluated unit's coefficients.
.lp_wgd_stage <- function(i, tech, peers, vrs, p = 1L,
                          hold_xp = TRUE, hold_xa = TRUE,
                          hold_a = TRUE, hold_quality = TRUE) {
  L <- length(peers)
  M <- tech$M
  x_abate <- tech$x_abate
  x_prod <- setdiff(seq_len(tech$N), x_abate)
  held <- c(if (hold_xp) x_prod, if (hold_xa) x_abate)
  H <- length(held)
  u_i <- matrix(tech$u[, , p], tech$L, tech$N)[i, ]
  v_i <- matrix(tech$v[, , p], tech$L, tech$M)[i, ]
  # peer values under the evaluated unit's coefficients, their true
  # uncontrolled emissions, and the quality gaps q_l = z_l - w_l
  w <- as.vector(tech$x %*% u_i) - as.vector(tech$y %*% v_i)
  z_acct <- .mb_cap(tech, p)
  q <- z_acct - w
  # Uncontrolled emissions close by definition (z = u'x - v'y). When
  # abatement is unobserved, the implied abatement is the closure gap
  # z_l - b_l. When abatement IS observed but an account does not close
  # exactly, the accounting z and the measured b + a disagree; the
  # stage programmes use the measurement-consistent z_l = b_l + a_l so
  # that peers' controlled emissions equal their measured b_l and the
  # fully freed stage collapses to the reduced form of Eq. 6. The
  # quality row keeps the coefficient-gap data q_l; both constructions
  # coincide when every account closes, the case of Rodseth (2025).
  a_l <- if (is.null(tech$a)) z_acct - tech$b[, p] else tech$a[, p]
  a_i <- a_l[i]
  z_data <- tech$b[, p] + a_l
  w <- z_data - q

  # columns: lambda (L), z, a, eps_z, eps_a, rho, eps_y (M), eps_x (H)
  iz <- L + 1L; ia <- L + 2L; iez <- L + 3L; iea <- L + 4L
  irho <- L + 5L
  iey <- L + 5L + seq_len(M)
  iex <- L + 5L + M + seq_len(H)
  nc <- L + 5L + M + H

  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  obj <- numeric(nc); obj[iz] <- 1; obj[ia] <- -1
  lpSolveAPI::set.objfn(lp, obj)

  for (m in seq_len(M)) {
    r <- numeric(nc); r[seq_len(L)] <- tech$y[peers, m]; r[iey[m]] <- -1
    lpSolveAPI::add.constraint(lp, r, "=", tech$y[i, m])
  }
  r <- numeric(nc); r[seq_len(L)] <- q[peers]; r[irho] <- -1
  lpSolveAPI::add.constraint(lp, r, "=", 0)
  r <- numeric(nc); r[seq_len(L)] <- w[peers]
  r[iez] <- 1; r[irho] <- 1; r[iz] <- -1
  lpSolveAPI::add.constraint(lp, r, "=", 0)
  r <- numeric(nc); r[seq_len(L)] <- a_l[peers]; r[iea] <- -1
  if (hold_a) {
    lpSolveAPI::add.constraint(lp, r, "=", a_i)
  } else {
    r[ia] <- -1
    lpSolveAPI::add.constraint(lp, r, "=", 0)
  }
  for (h in seq_len(H)) {
    n <- held[h]
    r <- numeric(nc); r[seq_len(L)] <- tech$x[peers, n]; r[iex[h]] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", tech$x[i, n])
  }
  # Freed inputs remain decision variables in Eq. 9: their disposal
  # slacks can absorb any amount, so with a freed positive-content
  # input the summing-up condition binds only from below
  # (eps_z + eps_a >= held terms); it stays an equality when every
  # input is held or the freed inputs carry no material content.
  freed <- setdiff(seq_len(tech$N), held)
  sum_op <- if (length(freed) && any(u_i[freed] > 0)) "<=" else "="
  r <- numeric(nc)
  if (H > 0) r[iex] <- u_i[held]
  r[iey] <- v_i; r[iez] <- -1; r[iea] <- -1
  lpSolveAPI::add.constraint(lp, r, sum_op, 0)
  if (vrs) {
    r <- numeric(nc); r[seq_len(L)] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", 1)
  }

  lower <- rep(0, nc)
  lower[c(iz, ia, irho)] <- -Inf
  upper <- rep(Inf, nc)
  if (hold_quality) {
    lower[irho] <- 0; upper[irho] <- 0
  }
  if (hold_a) {
    lower[ia] <- 0; upper[ia] <- 0
  }
  lpSolveAPI::set.bounds(lp, lower = lower, upper = upper)

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, z_star = NA_real_,
                a_star = NA_real_, rho = NA_real_, lambda = NULL))
  }
  sol <- lpSolveAPI::get.variables(lp)
  a_star <- if (hold_a) a_i else sol[ia]
  list(
    status = status,
    b_star = sol[iz] - a_star,
    z_star = sol[iz],
    a_star = a_star,
    rho = sol[irho],
    lambda = sol[seq_len(L)]
  )
}

# Principal environmental-efficiency score of one solution under `model`:
# the single definition shared by pgt(), boot_pgt() and compare_models().
.model_score <- function(model, sol, b_i) {
  if (!is.null(sol$status) && sol$status != 0) return(NA_real_)
  switch(model,
    byprod = sol$emission_eff,
    mb_cost = sol$mbe,
    sol$b_star / b_i)
}

# Emission-generating input columns for the by-production T2
# sub-technology: the user-supplied partition, else the inputs with a
# positive material flow coefficient for pollutant p, else all inputs.
.polluting_inputs <- function(tech, p = 1L) {
  if (length(tech$polluting)) {
    return(tech$polluting)
  }
  pos <- which(colSums(tech$u[, , p, drop = FALSE] > 0) > 0)
  if (length(pos)) pos else seq_len(tech$N)
}

# Peer index sets: pooled or own-group.
.peer_sets <- function(tech, peers) {
  if (peers == "all") {
    return(lapply(seq_len(tech$L), function(i) seq_len(tech$L)))
  }
  if (is.null(tech$group)) {
    stop("peers = \"group\" requires a 'group' in pgt_tech().",
         call. = FALSE)
  }
  idx <- split(seq_len(tech$L), tech$group)
  lapply(seq_len(tech$L), function(i) idx[[as.character(tech$group[i])]])
}
