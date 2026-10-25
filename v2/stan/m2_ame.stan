// ============================================================================
// m2_ame.stan  -- v2 strategic targeting, model family M2 (bespoke AME)
//
// Identical to m1_latent_additive.stan plus a multiplicative initiator x target
// term u_i' v_j with K latent dimensions in the stage-2 utility:
//
//   U[r] = -1/2 (x[r] - mu[n])' Sigma^{-1} (x[r] - mu[n]) + b_j + Z[r] gamma + u_i' v_j
//   mu[n] = mu0 + X[n] beta + a_i   (mu0 is the population-average ideal gap, added 2026-10-23;
//                                    see m1_latent_additive.stan; prior normal(0, prior_mu0_sd))
//
// Event weights ev_w (D8, 2026-10-25) and the "a" correlation prior (B4): identical to m1_latent_additive.stan.
// Choice-set flags in_iv / in_s2 (Codex review P1-1, 2026-10-24): identical to m1_latent_additive.stan.
// A state-year's row block is the union of its stage-1 IV set (January rule, in_iv = 1) and its
// stage-2 choice set (reference-date rule, in_s2 = 1); IV[n] is over the in_iv rows, the stage-2
// softmax over the in_s2 rows; blocks with every row in both sets (cs_mixed = 0) take the
// pre-review code path. With sampled = 1 every row must have in_iv = in_s2 = 1 (resolved in R).
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
  // Utility of the alternative rows s..e for one state-year, given its mu (row of length D),
  // the dyadic-control contribution zg_local (length e - s + 1) and the initiator's u_i.
  // D = 2: quadratic form in (dM, dR); D = 1: squared distance on column dcol only.
  vector util_block(int s, int e, matrix x, row_vector mu_n, vector sigma_x,
                    real rho, vector b_j, array[] int alt_target, vector zg_local, int D, int dcol,
                    row_vector u_i, matrix V) {
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
    return -0.5 * q + b_j[tg] + zg_local + V[tg] * u_i';
  }
  // Dyadic-control contribution for rows s..e (per block, so partial sums are independent).
  vector zg_block(matrix Z, vector gamma, int s, int e, int Pz) {
    int L = e - s + 1;
    if (Pz > 0) return block(Z, s, 1, L, Pz) * gamma;
    return rep_vector(0, L);
  }
  // Partial log-likelihood (stage 1 + stage 2) for a slice of state-years; used by reduce_sum.
  // sy_slice holds state-year indices; events of a state-year are ev_sorted[ev_first[n] .. +ev_count[n]-1].
  // cs_mixed[n] = 0: every row of the block is in both sets (IV = log_sum_exp(u), stage-2 softmax over u, the
  // pre-review path). cs_mixed[n] = 1: IV over the block-local IV rows iv_loc[iv_first[n] .. +iv_count[n]-1],
  // the stage-2 softmax over s2_loc[s2_first[n] .. +s2_count[n]-1] (sampled = 0 only).
  real sy_partial_sum(array[] int sy_slice, int start, int end,
                      array[] int y1, array[] int sy_state, matrix W, matrix X,
                      array[] int cs_start, array[] int cs_end, vector iv_corr,
                      array[] int cs_mixed, array[] int iv_first, array[] int iv_count, array[] int iv_loc,
                      array[] int s2_first, array[] int s2_count, array[] int s2_loc,
                      matrix x, array[] int alt_target, matrix Z, int Pz,
                      array[] int ev_first, array[] int ev_count, array[] int ev_sorted,
                      array[] int ev_row, array[] int ev_start, array[] int ev_end, vector ev_w, int sampled,
                      int D, int dcol,
                      real alpha0, vector alpha, real tau, row_vector mu0, matrix beta, matrix a_i, vector c_i,
                      vector sigma_x, real rho, vector b_j, vector gamma,
                      matrix U, matrix V) {
    int L = end - start + 1;
    vector[L] eta;
    real lp = 0;
    for (k in 1:L) {
      int n = sy_slice[k];
      int s = cs_start[n];
      int e = cs_end[n];
      row_vector[D] mu_n = mu0 + X[n] * beta + a_i[sy_state[n]];
      vector[e - s + 1] u;
      real lse;
      real lse2;
      u = util_block(s, e, x, mu_n, sigma_x, rho, b_j, alt_target, zg_block(Z, gamma, s, e, Pz), D, dcol, U[sy_state[n]], V);
      if (cs_mixed[n] == 0) {
        lse = log_sum_exp(u);                                    // both sets == the whole block
        lse2 = lse;
      } else {
        lse = log_sum_exp(u[segment(iv_loc, iv_first[n], iv_count[n])]);
        lse2 = (ev_count[n] > 0 && sampled == 0) ? log_sum_exp(u[segment(s2_loc, s2_first[n], s2_count[n])]) : lse;
      }
      eta[k] = alpha0 + W[n] * alpha + c_i[sy_state[n]] + tau * (lse + iv_corr[n]);
      for (j in 1:ev_count[n]) {
        int ev = ev_sorted[ev_first[n] + j - 1];
        if (sampled == 0) {
          lp += ev_w[ev] * (u[ev_row[ev] - s + 1] - lse2);       // event set == the state-year's in_s2 rows; weighted (D8)
        } else {
          int s2 = ev_start[ev];
          int e2 = ev_end[ev];
          vector[e2 - s2 + 1] u2 = util_block(s2, e2, x, mu_n, sigma_x, rho, b_j, alt_target,
                                              zg_block(Z, gamma, s2, e2, Pz), D, dcol, U[sy_state[n]], V);
          lp += ev_w[ev] * (u2[ev_row[ev] - s2 + 1] - log_sum_exp(u2));
        }
      }
    }
    lp += bernoulli_logit_lpmf(y1[sy_slice] | eta);
    return lp;
  }
}
data {
  int<lower=1> S;
  int<lower=1> N1;
  int<lower=1> N2;
  int<lower=0> E;
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
  array[N2] int<lower=0, upper=1> in_iv;   // row is in its state-year's stage-1 IV set (January rule)
  array[N2] int<lower=0, upper=1> in_s2;   // row is in its state-year's stage-2 choice set (reference-date rule)
  array[E] int<lower=1, upper=N1> ev_sy;
  array[E] int<lower=1, upper=N2> ev_row;      // chosen alternative's row
  array[E] int<lower=1, upper=N2> ev_start;    // the event's OWN choice-set rows
  array[E] int<lower=1, upper=N2> ev_end;
  vector<lower=0, upper=1>[E] ev_w;            // stage-2 event weights (D8, 2026-10-25): 1 / the initiator's targets in the same
                                               // dispute and year; multiply the stage-2 log-lik terms (composite likelihood);
                                               // all 1 = the unweighted program; log_lik_stage2 is unweighted per event
  int<lower=0, upper=1> joint;
  int<lower=0, upper=1> sampled;
  int<lower=0, upper=1> prior_only;     // 1 = skip both likelihoods (prior checks)
  int<lower=1, upper=3> D_mode;         // 1 = material only (dM), 2 = relational only (dR), 3 = both
  vector[P] x_ref;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_mu0_sd;           // sd of the normal prior on the mu intercept mu0
  int<lower=1> grainsize;               // reduce_sum_static slice size over state-years (>= N1 = one serial slice)
  // ---- AME
  int<lower=1> K;
  array[K] int<lower=1, upper=S> anchor;   // distinct state indices
}
transformed data {
  int D = (D_mode == 3) ? 2 : 1;        // targeting dimensions used in stage 2
  int dcol = (D_mode == 2) ? 2 : 1;     // which column of x when D = 1
  int Dre = 1 + D;                      // random-effect vector (c_i, a_i) length
  vector[N1] iv_corr = rep_vector(0, N1);
  array[N1] int sy_idx;                 // 1..N1, the array reduce_sum slices over
  array[N1] int ev_count = rep_array(0, N1);
  array[N1] int ev_first;
  array[E] int ev_sorted;               // event ids grouped by state-year (counting sort, original order kept within)
  for (n in 1:N1) sy_idx[n] = n;
  for (e in 1:E) ev_count[ev_sy[e]] += 1;
  {
    int pos = 1;
    for (n in 1:N1) { ev_first[n] = pos; pos += ev_count[n]; }
    array[N1] int fill = ev_first;
    for (e in 1:E) { ev_sorted[fill[ev_sy[e]]] = e; fill[ev_sy[e]] += 1; }
  }
  array[S - K] int nonanchor;              // state indices that are not anchors
  // ---- choice-set flags: block-local index lists of the IV rows and the stage-2 rows per state-year
  array[N1] int cs_mixed = rep_array(0, N1);   // 1 when some row of the block is outside the IV set or the stage-2 set
  array[N1] int iv_count = rep_array(0, N1);
  array[N1] int s2_count = rep_array(0, N1);
  array[N1] int iv_first;
  array[N1] int s2_first;
  for (n in 1:N1) {
    if (cs_end[n] < cs_start[n]) reject("empty choice set for state-year ", n);
    for (r in cs_start[n]:cs_end[n]) {
      if (in_iv[r] == 0 && in_s2[r] == 0) reject("row ", r, " of state-year ", n, " is in neither choice set");
      iv_count[n] += in_iv[r];
      s2_count[n] += in_s2[r];
    }
    if (iv_count[n] == 0) reject("state-year ", n, " has an EMPTY stage-1 IV set (no row with in_iv = 1)");
    if (iv_count[n] < cs_end[n] - cs_start[n] + 1 || s2_count[n] < cs_end[n] - cs_start[n] + 1) cs_mixed[n] = 1;
    if (sampled == 1 && cs_mixed[n] == 1) reject("state-year ", n, ": with sampled = 1 every row must have in_iv = in_s2 = 1");
    if (sampled == 1) iv_corr[n] = log(n_full[n]) - log(cs_end[n] - cs_start[n] + 1);
  }
  // sizes from the per-block counts: with sampled = 1 the event-block rows after the state-year blocks also carry
  // in_iv = in_s2 = 1 and must not be counted here
  int n_iv_rows = sum(iv_count);
  int n_s2_rows = sum(s2_count);
  array[n_iv_rows] int iv_loc;
  array[n_s2_rows] int s2_loc;
  {
    int p_iv = 1;
    int p_s2 = 1;
    for (n in 1:N1) {
      iv_first[n] = p_iv;
      s2_first[n] = p_s2;
      for (r in cs_start[n]:cs_end[n]) {
        if (in_iv[r] == 1) { iv_loc[p_iv] = r - cs_start[n] + 1; p_iv += 1; }
        if (in_s2[r] == 1) { s2_loc[p_s2] = r - cs_start[n] + 1; p_s2 += 1; }
      }
    }
    if (p_iv != n_iv_rows + 1 || p_s2 != n_s2_rows + 1) reject("choice-set flag bookkeeping failed");
  }
  for (e in 1:E) {
    if (ev_row[e] < ev_start[e] || ev_row[e] > ev_end[e])
      reject("event ", e, " chosen row outside its own choice set");
    if (in_s2[ev_row[e]] != 1) reject("event ", e, ": chosen row is not in the stage-2 set (in_s2 = 0)");
    if (y1[ev_sy[e]] != 1) reject("event ", e, " belongs to a state-year with y1 == 0");
    if (sampled == 0 && (ev_start[e] != cs_start[ev_sy[e]] || ev_end[e] != cs_end[ev_sy[e]]))
      reject("event ", e, ": with full choice sets the event set must be the state-year set");
    if (sampled == 0 && s2_count[ev_sy[e]] < 2)
      reject("event ", e, ": its state-year's stage-2 set is empty or holds only the chosen row");
    if (sampled == 1 && ev_start[e] <= cs_end[N1])
      reject("event ", e, ": sampled event rows must follow the state-year (IV) rows");
    if (sampled == 1 && ev_end[e] - ev_start[e] + 1 < 2)
      reject("event ", e, ": sampled event set holds only the chosen row");
    if (sampled == 1) for (r in ev_start[e]:ev_end[e]) if (in_iv[r] == 0 || in_s2[r] == 0)
      reject("event ", e, ": with sampled = 1 every event-block row must have in_iv = in_s2 = 1");
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
  row_vector[D] mu0;                                             // mu intercept (population-average ideal gap)
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
  profile("priors") {
    alpha0 ~ normal(0, 2.5);
    alpha ~ normal(0, 1);
    tau ~ normal(0, 1);
    mu0 ~ normal(0, prior_mu0_sd);
    to_vector(beta) ~ normal(0, prior_beta_sd);
    sigma_x ~ lognormal(0, 1);
    sd_re ~ normal(0, 1);
    if (D == 1) {
      // one targeting dimension: corr(c, a) gets the same Beta(2.5, 2.5) marginal (Var 1/6) as the
      // shared correlations in the a / b models, so m / r / a / b differ only in their restrictions
      target += beta_lpdf((r_c_a1 + 1) / 2 | 2.5, 2.5);
    } else if (joint == 1) {
      L_Om_joint[1] ~ lkj_corr_cholesky(2);
      // Sigma's correlation (b only): explicit (rho + 1) / 2 ~ Beta(2.5, 2.5), the same marginal as the shared
      // correlations (2026-10-25 P3; before this the prior was the implicit uniform(-1, 1) of the bounds)
      target += beta_lpdf((rho_x[1] + 1) / 2 | 2.5, 2.5);
    } else {
      // "a" prior on (r1, r2) = (corr(c, a1), corr(c, a2)) with corr(a1, a2) = 0 (review fix 2026-10-25, B4). LKJ(2) on the
      // 3x3 matrix is proportional to det(R) = 1 - r1^2 - r2^2 on the slice corr(a1, a2) = 0, so the density on the unit
      // disk is p(r1, r2) propto (1 - r1^2 - r2^2); its marginals are exactly (r + 1) / 2 ~ Beta(2.5, 2.5) (Var 1/6), the
      // LKJ(2) 3x3 marginals of the "b" model. (Two independent Beta(2.5, 2.5) marginals truncated to the disk, the
      // previous code, give Var about 0.158 and are not that density.) Parameterisation r2 = q sqrt(1 - r1^2), q in (-1, 1):
      // Jacobian |dr2 / dq| = sqrt(1 - r1^2), i.e. 0.5 * log1m(r1^2).
      target += log1m(square(r_c_a1) + square(r_c_a2));
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
  }
  // ---- likelihood: stage 1 + stage 2, sliced over state-years
  if (prior_only == 0) {
    target += reduce_sum_static(sy_partial_sum, sy_idx, grainsize,
                                y1, sy_state, W, X, cs_start, cs_end, iv_corr,
                                cs_mixed, iv_first, iv_count, iv_loc, s2_first, s2_count, s2_loc,
                                x, alt_target, Z, Pz,
                                ev_first, ev_count, ev_sorted, ev_row, ev_start, ev_end, ev_w, sampled, D, dcol,
                                alpha0, alpha, tau, mu0, beta, a_i, c_i, sigma_x, rho, b_j, gamma, U, V);
  }
}
generated quantities {
  vector[E] log_lik_stage2;
  vector[N1] log_lik_stage1;
  matrix[S, D] mu_state;      // mu_i at x_ref: mu0 + x_ref' beta + a_i   (D columns)
  matrix[E, D] mu_event;      // mu at the event's own X: mu0 + X[ev_sy] beta + a_i
  matrix[D, D] Sigma;
  matrix[Dre, Dre] Omega;     // covariance of (c, a_1[, a_2])
  corr_matrix[Dre] R_Om;      // correlation of (c, a_1[, a_2])
  {
    matrix[N1, D] Mu = rep_matrix(mu0, N1) + X * beta + a_i[sy_state];
    vector[N2] util;
    vector[N1] lse;       // IV: over the in_iv rows
    vector[N1] lse_s2;    // stage-2 normaliser (sampled = 0): over the in_s2 rows
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block(s, e, x, Mu[n], sigma_x, rho, b_j, alt_target, zg_block(Z, gamma, s, e, Pz), D, dcol,
                             U[sy_state[n]], V);
      if (cs_mixed[n] == 0) {
        lse[n] = log_sum_exp(util[s:e]);
        lse_s2[n] = lse[n];
      } else {
        vector[e - s + 1] u = util[s:e];
        lse[n] = log_sum_exp(u[segment(iv_loc, iv_first[n], iv_count[n])]);
        lse_s2[n] = (s2_count[n] > 0) ? log_sum_exp(u[segment(s2_loc, s2_first[n], s2_count[n])]) : not_a_number();
      }
    }
    vector[N1] eta1 = alpha0 + W * alpha + c_i[sy_state] + tau * (lse + iv_corr);
    for (n in 1:N1) log_lik_stage1[n] = bernoulli_logit_lpmf(y1[n] | eta1[n]);
    for (ev in 1:E) {
      int s = ev_start[ev];
      int e = ev_end[ev];
      if (sampled == 1) {
        util[s:e] = util_block(s, e, x, Mu[ev_sy[ev]], sigma_x, rho, b_j, alt_target, zg_block(Z, gamma, s, e, Pz), D, dcol,
                               U[sy_state[ev_sy[ev]]], V);
        log_lik_stage2[ev] = util[ev_row[ev]] - log_sum_exp(util[s:e]);
      } else {
        log_lik_stage2[ev] = util[ev_row[ev]] - lse_s2[ev_sy[ev]];
      }
      mu_event[ev] = Mu[ev_sy[ev]];
    }
    for (s in 1:S) mu_state[s] = mu0 + x_ref' * beta + a_i[s];
    if (D == 2) Sigma = [[sigma_x[1]^2, rho * sigma_x[1] * sigma_x[2]],
                         [rho * sigma_x[1] * sigma_x[2], sigma_x[2]^2]];
    else Sigma = [[sigma_x[1]^2]];
    R_Om = multiply_lower_tri_self_transpose(L_Om);
    Omega = quad_form_diag(R_Om, sd_re);
  }
}
