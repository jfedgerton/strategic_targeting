// ============================================================================
// m2_ame.stan  -- v2 strategic targeting, model family M2 (bespoke AME)
//
// Identical to m1_latent_additive.stan plus a multiplicative initiator x target
// term u_i' v_j with K latent dimensions in the stage-2 utility:
//
//   U[r] = -1/2 (x[r] - mu[n])' Sigma^{-1} (x[r] - mu[n]) + b_j + Z[r] gamma + u_i' v_j
//
// U and V are transformed parameters and are saved; u_i'v_j is recomputed in R
// (an S x S generated quantity would be too large for real-data S).
// Identification of U (S x K) and V (S x K):
//   * V is column-centred (a column constant in V adds u_i' c, constant across
//     alternatives for a given initiator, so it cancels in the softmax).
//   * U is column-centred too: adding h to U[, k] and replacing b_j by b_j - h V[j, k]
//     leaves the likelihood unchanged, so without it u_i'v_j would be identified only up
//     to a target-specific shift absorbed by b_j. The centring is a hard constraint: the
//     entry of the LAST non-anchor state in each column equals minus the sum of the other
//     entries, so the anchor block (zeros above the diagonal, positive diagonal) is intact.
//     With both U and V centred and the anchor block, u_i'v_j is identified up to the
//     rotation handled by the anchors and the Procrustes step.
//   * The K rows of U belonging to the `anchor` states form a lower-triangular
//     block with a positive diagonal (fixes rotation and reflection). Anchors must
//     be supplied by the data: v2 uses the K states with the MOST initiations in
//     the estimation sample (audit F6: the old code anchored on whichever ccodes
//     sorted first as strings). Non-anchor rows of U are free.
//   * Scale between U and V (u_k -> c u_k, v_k -> v_k / c) is fixed only by the
//     N(0,1) priors on both; this is standard for AME-type models.
//   * Post-hoc: 03_extract.R applies an orthogonal Procrustes alignment of each
//     chain's U/V draws to a reference chain before pooling (safety net; with the
//     anchor block the rotation should already be identified).
// K must be >= 1 and <= S - 1. Works for K = 1, 2, 3.
// ============================================================================
functions {
  // D = 2: quadratic form in (dM, dR); D = 1: squared distance on column dcol only.
  vector util_block(int n, int s, int e, matrix x, row_vector mu_n, vector sigma_x,
                    real rho, vector b_j, array[] int alt_target, vector zg,
                    row_vector u_i, matrix V, int D, int dcol) {
    int L = e - s + 1;
    vector[L] q;
    if (D == 2) {
      vector[L] d1 = (segment(col(x, 1), s, L) - mu_n[1]) / sigma_x[1];
      vector[L] d2 = (segment(col(x, 2), s, L) - mu_n[2]) / sigma_x[2];
      q = (d1 .* d1 - 2 * rho * (d1 .* d2) + d2 .* d2) / (1 - rho * rho);
    } else {
      vector[L] d1 = (segment(col(x, dcol), s, L) - mu_n[1]) / sigma_x[1];
      q = d1 .* d1;
    }
    array[L] int tg = segment(alt_target, s, L);
    return -0.5 * q + b_j[tg] + segment(zg, s, L) + V[tg] * u_i';
  }
}
data {
  int<lower=1> S;
  int<lower=1> N1;
  int<lower=1> N2;
  int<lower=1> E;
  int<lower=0> P1;
  int<lower=1> P;
  int<lower=0> Pz;
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=S> sy_state;
  matrix[N1, P1] W;
  matrix[N1, P] X;
  array[N1] int<lower=1, upper=N2> cs_start;
  array[N1] int<lower=1, upper=N2> cs_end;
  array[N1] int<lower=1> n_full;
  array[N2] int<lower=1, upper=S> alt_target;
  matrix[N2, 2] x;
  matrix[N2, Pz] Z;
  array[E] int<lower=1, upper=N1> ev_sy;
  array[E] int<lower=1, upper=N2> ev_row;      // chosen alternative's row
  array[E] int<lower=1, upper=N2> ev_start;    // the event's OWN choice-set rows
  array[E] int<lower=1, upper=N2> ev_end;
  int<lower=0, upper=1> joint;
  int<lower=0, upper=1> sampled;
  int<lower=0, upper=1> prior_only;     // 1 = skip both likelihoods (prior checks)
  int<lower=1, upper=3> D_mode;         // 1 = material only (dM), 2 = relational only (dR), 3 = both
  vector[P] x_ref;
  real<lower=0> prior_beta_sd;
  // ---- AME
  int<lower=1> K;
  array[K] int<lower=1, upper=S> anchor;   // distinct state indices
}
transformed data {
  int D = (D_mode == 3) ? 2 : 1;        // targeting dimensions used in stage 2
  int dcol = (D_mode == 2) ? 2 : 1;     // which column of x when D = 1
  int Dre = 1 + D;                      // random-effect vector (c_i, a_i) length
  vector[N1] iv_corr = rep_vector(0, N1);
  array[S - K] int nonanchor;              // state indices that are not anchors
  for (n in 1:N1) {
    if (cs_end[n] < cs_start[n]) reject("empty choice set for state-year ", n);
    if (sampled == 1) iv_corr[n] = log(n_full[n]) - log(cs_end[n] - cs_start[n] + 1);
  }
  for (e in 1:E) {
    if (ev_row[e] < ev_start[e] || ev_row[e] > ev_end[e])
      reject("event ", e, " chosen row outside its own choice set");
    if (y1[ev_sy[e]] != 1) reject("event ", e, " belongs to a state-year with y1 == 0");
    if (sampled == 0 && (ev_start[e] != cs_start[ev_sy[e]] || ev_end[e] != cs_end[ev_sy[e]]))
      reject("event ", e, ": with full choice sets the event set must be the state-year set");
    if (sampled == 1 && ev_start[e] <= cs_end[N1])
      reject("event ", e, ": sampled event rows must follow the state-year (IV) rows");
  }
  if (K > S - 1) reject("K must be <= S - 1");
  {
    int idx = 1;
    for (s in 1:S) {
      int is_anchor = 0;
      for (k in 1:K) if (anchor[k] == s) is_anchor += 1;
      if (is_anchor > 1) reject("duplicate anchor state ", s);
      if (is_anchor == 0) {
        nonanchor[idx] = s;
        idx += 1;
      }
    }
    if (idx != S - K + 1) reject("anchor bookkeeping failed");
  }
}
parameters {
  real alpha0;
  vector[P1] alpha;
  real tau;
  matrix[P, D] beta;
  vector<lower=0>[D] sigma_x;
  array[joint * (D - 1)] real<lower=-1, upper=1> rho_x;          // D = 2, joint only
  vector<lower=0>[Dre] sd_re;                                    // sd of (c, a_1[, a_2])
  array[joint * (D - 1)] cholesky_factor_corr[3] L_Om_joint;     // D = 2, joint only
  array[(1 - joint) * (D - 1)] real<lower=-1, upper=1> r_ca1;    // D = 2, independent only
  array[(1 - joint) * (D - 1)] real<lower=-1, upper=1> q_ca2;
  array[2 - D] real<lower=-1, upper=1> r_ca_1d;                  // D = 1: corr(c, a)
  matrix[Dre, S] z_re;
  real<lower=0> sigma_b;
  vector[S] z_b;
  vector[Pz] gamma;
  // ---- AME
  vector<lower=0>[K] u_diag;               // positive diagonal of the anchor block
  vector[K * (K - 1) / 2] u_lower;         // strictly lower-triangular anchor entries
  matrix[S - K - 1, K] U_free;             // rows for non-anchor states except the last (determined by centring)
  matrix[S, K] V_raw;
}
transformed parameters {
  matrix[Dre, Dre] L_Om;
  real rho;
  real r_c_a1;                 // corr(c, a1), shared by the "a" and "b" models (corr(c, a) when D = 1)
  real r_c_a2;                 // corr(c, a2) (0 when D = 1)
  if (D == 2) {
    if (joint == 1) {
      L_Om = L_Om_joint[1];
      rho = rho_x[1];
      r_c_a1 = L_Om[2, 1] * L_Om[1, 1];
      r_c_a2 = L_Om[3, 1] * L_Om[1, 1];
    } else {
      real r1 = r_ca1[1];
      real r2 = q_ca2[1] * sqrt(1 - r1 * r1);        // guarantees r1^2 + r2^2 < 1 (PD)
      matrix[3, 3] R_a = [[1, r1, r2], [r1, 1, 0], [r2, 0, 1]];
      L_Om = cholesky_decompose(R_a);
      rho = 0;
      r_c_a1 = r1;
      r_c_a2 = r2;
    }
  } else {
    real r1 = r_ca_1d[1];
    matrix[2, 2] R_1 = [[1, r1], [r1, 1]];
    L_Om = cholesky_decompose(R_1);
    rho = 0;
    r_c_a1 = r1;
    r_c_a2 = 0;
  }
  matrix[S, Dre] re = (diag_pre_multiply(sd_re, L_Om) * z_re)';
  vector[S] c_i = col(re, 1);
  matrix[S, D] a_i = re[, 2:Dre];
  vector[S] b_j = sigma_b * (z_b - mean(z_b));
  // ---- assemble U with the anchor block lower-triangular, positive diagonal
  matrix[S, K] U;
  matrix[S, K] V;
  {
    int idx = 1;
    for (k in 1:K) {
      for (l in 1:K) {
        if (l == k) U[anchor[k], l] = u_diag[k];
        else if (l < k) { U[anchor[k], l] = u_lower[idx]; idx += 1; }
        else U[anchor[k], l] = 0;
      }
    }
    if (S - K > 1) U[nonanchor[1:(S - K - 1)]] = U_free;
    U[nonanchor[S - K]] = rep_row_vector(0, K);
    for (k in 1:K) U[nonanchor[S - K], k] = -sum(col(U, k));   // hard column sum-to-zero
    for (k in 1:K) V[, k] = V_raw[, k] - mean(V_raw[, k]);
  }
}
model {
  alpha0 ~ normal(0, 2.5);
  alpha ~ normal(0, 1);
  tau ~ normal(0, 1);
  to_vector(beta) ~ normal(0, prior_beta_sd);
  sigma_x ~ lognormal(0, 1);
  sd_re ~ normal(0, 1);
  if (D == 1) {
    // one targeting dimension: corr(c, a) gets the same Beta(2.5, 2.5) marginal (Var 1/6) as the
    // shared correlations in the a / b models, so m / r / a / b differ only in their restrictions
    target += beta_lpdf((r_c_a1 + 1) / 2 | 2.5, 2.5);
  } else if (joint == 1) {
    L_Om_joint[1] ~ lkj_corr_cholesky(2);
  } else {
    // Same marginal prior on the shared correlations as LKJ(2) implies for a 3x3 matrix:
    // (r + 1) / 2 ~ Beta(2.5, 2.5), i.e. Var(r) = 1/6, applied symmetrically to corr(c, a1)
    // and corr(c, a2), with corr(a1, a2) = 0 and PD enforced through r2 = q sqrt(1 - r1^2)
    // (Jacobian 0.5 * log(1 - r1^2)). So "a" differs from "b" only by the cross-dimensional
    // restrictions (Sigma off-diagonal and corr(a1, a2)), not by the priors they share.
    target += beta_lpdf((r_c_a1 + 1) / 2 | 2.5, 2.5) + beta_lpdf((r_c_a2 + 1) / 2 | 2.5, 2.5);
    target += 0.5 * log1m(square(r_c_a1));
  }
  to_vector(z_re) ~ std_normal();
  sigma_b ~ normal(0, 1);
  z_b ~ std_normal();
  gamma ~ normal(0, 1);
  u_diag ~ normal(0, 1);
  u_lower ~ normal(0, 1);
  to_vector(U[nonanchor]) ~ normal(0, 1);  // all non-anchor rows, including the determined one
  to_vector(V_raw) ~ normal(0, 1);
  {
    matrix[N1, D] Mu = X * beta + a_i[sy_state];
    vector[N2] zg = (Pz > 0) ? Z * gamma : rep_vector(0, N2);
    vector[N2] util;
    vector[N1] lse;
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block(n, s, e, x, Mu[n], sigma_x, rho, b_j, alt_target, zg,
                             U[sy_state[n]], V, D, dcol);
      lse[n] = log_sum_exp(util[s:e]);
    }
    if (prior_only == 0) y1 ~ bernoulli_logit(alpha0 + W * alpha + c_i[sy_state] + tau * (lse + iv_corr));
    for (ev in 1:E) {
      int s = ev_start[ev];
      int e = ev_end[ev];
      if (sampled == 1)
        util[s:e] = util_block(ev_sy[ev], s, e, x, Mu[ev_sy[ev]], sigma_x, rho, b_j, alt_target, zg,
                               U[sy_state[ev_sy[ev]]], V, D, dcol);
      if (prior_only == 0) target += util[ev_row[ev]] - log_sum_exp(util[s:e]);
    }
  }
}
generated quantities {
  vector[E] log_lik_stage2;
  vector[N1] log_lik_stage1;
  matrix[S, D] mu_state;      // mu_i at x_ref: x_ref' beta + a_i   (D columns)
  matrix[E, D] mu_event;      // mu at the event's own X: X[ev_sy] beta + a_i
  matrix[D, D] Sigma;
  matrix[Dre, Dre] Omega;     // covariance of (c, a_1[, a_2])
  corr_matrix[Dre] R_Om;      // correlation of (c, a_1[, a_2])
  {
    matrix[N1, D] Mu = X * beta + a_i[sy_state];
    vector[N2] zg = (Pz > 0) ? Z * gamma : rep_vector(0, N2);
    vector[N2] util;
    vector[N1] lse;
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block(n, s, e, x, Mu[n], sigma_x, rho, b_j, alt_target, zg,
                             U[sy_state[n]], V, D, dcol);
      lse[n] = log_sum_exp(util[s:e]);
    }
    vector[N1] eta1 = alpha0 + W * alpha + c_i[sy_state] + tau * (lse + iv_corr);
    for (n in 1:N1) log_lik_stage1[n] = bernoulli_logit_lpmf(y1[n] | eta1[n]);
    for (ev in 1:E) {
      int s = ev_start[ev];
      int e = ev_end[ev];
      if (sampled == 1)
        util[s:e] = util_block(ev_sy[ev], s, e, x, Mu[ev_sy[ev]], sigma_x, rho, b_j, alt_target, zg,
                               U[sy_state[ev_sy[ev]]], V, D, dcol);
      log_lik_stage2[ev] = util[ev_row[ev]] - log_sum_exp(util[s:e]);
      mu_event[ev] = Mu[ev_sy[ev]];
    }
    for (s in 1:S) mu_state[s] = x_ref' * beta + a_i[s];
    if (D == 2) Sigma = [[sigma_x[1]^2, rho * sigma_x[1] * sigma_x[2]],
                         [rho * sigma_x[1] * sigma_x[2], sigma_x[2]^2]];
    else Sigma = [[sigma_x[1]^2]];
    R_Om = multiply_lower_tri_self_transpose(L_Om);
    Omega = quad_form_diag(R_Om, sd_re);
  }
}
