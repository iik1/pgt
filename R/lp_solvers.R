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

# Pollutant potential u_i'x_i and materials-balance cap u_i'x_i - v_i y_i
# for pollutant p (one value per DMU).
.mb_potential <- function(tech, p = 1L) {
  rowSums(tech$u[, , p, drop = TRUE] * tech$x)
}
.mb_cap <- function(tech, p = 1L) {
  .mb_potential(tech, p) - tech$v[, p] * tech$y
}

# Weak-G-disposability LP (Rodseth Eq. 6, reduced form) for pollutant p.
#
#   min_{lambda, bq} bq
#   s.t. sum_l lambda_l y_l          >= y_i
#        sum_l lambda_l x_nl         <= x_ni   (n = 1..N, optional)
#        sum_l lambda_l b_l  -  bq   <= 0
#        bq                          <= u_i'x_i - v_i y_i    [materials balance]
#        sum_l lambda_l b_ql         <= u_qi'x_i - v_qi y_i  [per extra pollutant q]
#        sum_l lambda_l               = 1                    [VRS only]
#        lambda >= 0, bq >= 0
#
# The materials-balance row keeps the projected bad output consistent
# with the evaluated DMU's own pollutant potential. With several
# pollutants the peer mix's emission of every other pollutant q must
# respect that pollutant's cap for the evaluated DMU as well (the
# envelope and cap rows for q collapse to one row because q carries no
# objective term). A DMU violating any pollutant's identity loses
# self-reference feasibility, so its LP solves only if some peer mix
# meets every row within the caps; infeasibility (status != 0) is
# therefore confined to violating DMUs, but most violators still solve.
.lp_wgd_one <- function(i, X, y, b, mb_rhs, peers, vrs = TRUE,
                        input_constraints = TRUE,
                        b_other = NULL, cap_other = NULL) {
  L <- length(peers)
  N <- ncol(X)

  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  lpSolveAPI::set.objfn(lp, c(rep(0, L), 1))

  lpSolveAPI::add.constraint(lp, c(y[peers], 0), ">=", y[i])
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
                dual_output = NA_real_, mb_rhs = mb_rhs[i]))
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
    dual_output = duals[1],
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
.lp_envelope_one <- function(i, y, b, peers, vrs = TRUE) {
  L <- length(peers)
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  lpSolveAPI::set.objfn(lp, b[peers])
  lpSolveAPI::add.constraint(lp, y[peers], ">=", y[i])
  if (vrs) {
    lpSolveAPI::add.constraint(lp, rep(1, L), "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = rep(0, L))

  status <- lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(status = status, b_star = NA_real_, lambda = NULL,
                dual_output = NA_real_, mb_rhs = NA_real_))
  }
  duals <- lpSolveAPI::get.sensitivity.rhs(lp)$duals
  list(
    status = status,
    b_star = lpSolveAPI::get.objective(lp),
    lambda = lpSolveAPI::get.variables(lp),
    dual_output = duals[1],
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
# reported as the headline; E1 and E_FGL are returned alongside.
.lp_byprod_one <- function(i, X, y, b, pol, peers, vrs = TRUE) {
  L <- length(peers)
  N <- ncol(X)

  # T1: output efficiency (max phi).
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
  invisible(lpSolveAPI::lp.control(lp1, sense = "max"))
  o <- numeric(L + 1); o[L + 1] <- 1; lpSolveAPI::set.objfn(lp1, o)
  r <- c(y[peers], -y[i]); lpSolveAPI::add.constraint(lp1, r, ">=", 0)
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
.lp_mbcost_one <- function(i, X, y, u_row, v_i, peers, vrs = TRUE) {
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
  lpSolveAPI::add.constraint(lp1, c(y[peers], 0), ">=", y[i])
  if (vrs) lpSolveAPI::add.constraint(lp1, c(rep(1, L), 0), "=", 1)
  lpSolveAPI::set.bounds(lp1, lower = c(rep(0, L), 0))
  st1 <- lpSolveAPI::solve.lpExtPtr(lp1)
  te <- if (st1 == 0) lpSolveAPI::get.variables(lp1)[L + 1] else NA_real_

  # material-inflow minimisation (min u_i'xproj)
  nc <- L + N
  lp2 <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp2, sense = "min"))
  o <- numeric(nc); o[(L + 1):nc] <- u_row; lpSolveAPI::set.objfn(lp2, o)
  lpSolveAPI::add.constraint(lp2, c(y[peers], rep(0, N)), ">=", y[i])
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
       b_star = pot_star - v_i * y[i], lambda = lam)
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
.lp_wd_one <- function(i, X, y, b, peers, vrs = TRUE) {
  L <- length(peers)
  N <- ncol(X)
  # columns: z_1..L, w_1..L, phi
  nc <- 2L * L + 1L
  iphi <- nc
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp, sense = "min"))
  o <- numeric(nc); o[iphi] <- 1; lpSolveAPI::set.objfn(lp, o)

  # good output: sum z y >= y_i
  r <- numeric(nc); r[seq_len(L)] <- y[peers]
  lpSolveAPI::add.constraint(lp, r, ">=", y[i])
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
  if (model %in% c("wgd", "fdmo")) {
    ctx$mb_cap <- .mb_cap(tech, p)
  }
  if (model == "wgd" && tech$P > 1L) {
    q <- setdiff(seq_len(tech$P), p)
    ctx$b_other <- tech$b[, q, drop = FALSE]
    ctx$cap_other <- vapply(q, function(pp) .mb_cap(tech, pp),
                            numeric(tech$L))
    if (!is.matrix(ctx$cap_other)) {
      ctx$cap_other <- matrix(ctx$cap_other, ncol = length(q))
    }
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
    wgd = .lp_wgd_one(i, tech$x, tech$y, ctx$b_p, ctx$mb_cap,
                      peers, vrs = vrs,
                      input_constraints = input_constraints,
                      b_other = ctx$b_other, cap_other = ctx$cap_other),
    envelope = .lp_envelope_one(i, tech$y, ctx$b_p, peers, vrs = vrs),
    fdmo = {
      if (is.null(tech$a)) {
        stop("model = \"fdmo\" requires an abatement output 'a' in ",
             "pgt_tech().", call. = FALSE)
      }
      .lp_fdmo_one(i, tech$x, tech$y, tech$a[, p], tech$v[, p],
                   ctx$mb_cap, ctx$b_p,
                   peers, vrs = vrs, input_constraints = input_constraints)
    },
    byprod = .lp_byprod_one(i, tech$x, tech$y, ctx$b_p,
                            ctx$pol, peers, vrs = vrs),
    mb_cost = .lp_mbcost_one(i, tech$x, tech$y, tech$u[i, , p],
                             tech$v[i, p], peers, vrs = vrs),
    wd = .lp_wd_one(i, tech$x, tech$y, ctx$b_p, peers, vrs = vrs),
    stop("unknown model '", model, "'", call. = FALSE)
  )
}

# Headline environmental-efficiency score of one solution under `model`:
# the single definition shared by pgt(), boot_pgt() and compare_models().
.headline_score <- function(model, sol, b_i) {
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
