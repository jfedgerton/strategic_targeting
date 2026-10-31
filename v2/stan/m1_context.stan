// ============================================================================
// m1_context.stan  -- v2 strategic targeting, model M1c (dual typology stage b, PI-approved spec 2026-10-27)
//
// Derived from the joint ("b") variant of m1_latent_additive.stan: the same two-stage spatial choice model, estimated
// jointly, with ONE change in stage 2 and ONE in stage 1.
//
// Stage 1 (monadic, state-year n):  y1[n] ~ Bernoulli(inv_logit(eta1[n]))
//   eta1[n] = alpha0 + W[n] alpha + c_i + tau * IV[n]
//   W carries one extra outcome-blind control (multi_mid_jan, appended by assemble_cell(CONTEXT_MODEL = TRUE)); the
//   program does not know which column it is (P1 is data).
//
// Stage 2 (dyadic choice, event e made by state-year n = ev_sy[e], in group g = ev_grp[e]):
//   U_g[r] = -1/2 (x[r] - mu_e)' Sigma_g^{-1} (x[r] - mu_e) + b_j + Z[r] gamma
//   mu_e   = mu0 + delta_lead * [g == 1] + delta_join * [g == 2] + X[n] beta + a_i          (row_vector[2])
//   Sigma_g = diag(sigma_g[g]) R diag(sigma_g[g]),  R = [[1, rho], [rho, 1]]  shared across the three groups
//   P(chosen = r | e) = softmax of U_g over the state-year's stage-2 set (in_s2 rows)
// Groups (data ev_grp, from assemble_cell()): 0 = bilateral (reference: no shift, sigma_g[1]), 1 = coalition leader
// (delta_lead, sigma_g[2]), 2 = joiner (delta_join, sigma_g[3]). A state-year can hold events of several groups (a
// state leading one coalition and joining another in the same year), so the utilities are computed per group present.
//
// The inclusive value IV[n] = log_sum_exp over the in_iv rows of the BILATERAL utility U_0 (mu_n = mu0 + X[n] beta + a_i,
// Sigma_1): the opportunity value of acting alone, defined for every state-year whether or not it initiates (a
// non-initiating state-year has no group). Implementation decision, flagged for the PI in the README.
//
// Priors: as the b variant (mu0 ~ normal(0, prior_mu0_sd); beta ~ normal(0, prior_beta_sd); each sigma_g[g] ~ lognormal(0, 1),
// the prior of sigma_x; sd_re ~ half-normal(0, 1); LKJ(2) on the 3x3 correlation of (c, a_1, a_2); (rho + 1)/2 ~ Beta(2.5, 2.5);
// gamma, alpha, tau ~ normal(0, 1); alpha0 ~ normal(0, 2.5)), plus delta_lead, delta_join ~ normal(0, prior_mu0_sd) per element,
// the same family and scale as mu0 (a shift of the population ideal point is on the scale of the ideal point itself).
//
// Event weights (D8), the choice-set flags (in_iv / in_s2), the reduce_sum_static slicing over state-years and the generated
// quantities follow m1_latent_additive.stan. Restrictions: D_mode = 3 (both axes), joint = 1 and sampled = 0 are REQUIRED
// (rejected otherwise); the data list of build_stan_data() plus ev_grp is the input.
// Generated quantities: log_lik_stage1[N1], log_lik_stage2[E] (unweighted, with the event's own group), mu_state[S, 2] (the
// BILATERAL ideal point at x_ref: mu0 + x_ref' beta + a_i; the leader / joiner ideal points are mu_state + delta_g),
// mu_event[E, 2] (at the event's own X, WITH its group shift), Sigma (bilateral), Sigma_g[3], mu0_lead, mu0_join, Omega, R_Om.
// ============================================================================
functions {
  // Utility of the alternative rows s..e for one state-year and one group: quadratic form in (dM, dR) around mu_n with scales
  // sig (length 2) and the shared correlation rho, plus the target intercept and the dyadic controls.
  vector util_block_c(int s, int e, matrix x, row_vector mu_n, vector sig, real rho, vector b_j, array[] int alt_target, vector zg_local) {
    int L = e - s + 1;
    vector[L] d1 = (segment(col(x, 1), s, L) - mu_n[1]) / sig[1];
    vector[L] d2 = (segment(col(x, 2), s, L) - mu_n[2]) / sig[2];
    vector[L] q = (d1 .* d1 - 2 * rho * (d1 .* d2) + d2 .* d2) / (1 - rho * rho);
    return -0.5 * q + b_j[segment(alt_target, s, L)] + zg_local;
  }
  vector zg_block(matrix Z, vector gamma, int s, int e, int Pz) {
    int L = e - s + 1;
    if (Pz > 0) return block(Z, s, 1, L, Pz) * gamma;
    return rep_vector(0, L);
  }
  // Partial log-likelihood (stage 1 + stage 2) for a slice of state-years; used by reduce_sum.
  real sy_partial_sum_c(array[] int sy_slice, int start, int end,
                        array[] int y1, array[] int sy_state, matrix W, matrix X,
                        array[] int cs_start, array[] int cs_end,
                        array[] int cs_mixed, array[] int iv_first, array[] int iv_count, array[] int iv_loc,
                        array[] int s2_first, array[] int s2_count, array[] int s2_loc,
                        matrix x, array[] int alt_target, matrix Z, int Pz,
                        array[] int ev_first, array[] int ev_count, array[] int ev_sorted,
                        array[] int ev_row, array[] int ev_grp, vector ev_w, array[] int has_lead, array[] int has_join,
                        real alpha0, vector alpha, real tau, row_vector mu0, row_vector delta_lead, row_vector delta_join,
                        matrix beta, matrix a_i, vector c_i, array[] vector sigma_g, real rho, vector b_j, vector gamma) {
    int L = end - start + 1;
    vector[L] eta;
    real lp = 0;
    for (k in 1:L) {
      int n = sy_slice[k];
      int s = cs_start[n];
      int e = cs_end[n];
      row_vector[2] mu_n = mu0 + X[n] * beta + a_i[sy_state[n]];
      vector[e - s + 1] zg = zg_block(Z, gamma, s, e, Pz);
      vector[e - s + 1] u = util_block_c(s, e, x, mu_n, sigma_g[1], rho, b_j, alt_target, zg);   // bilateral: IV and group-0 events
      vector[e - s + 1] u_ld;
      vector[e - s + 1] u_jn;
      real lse;
      real lse2 = 0;
      real lse_ld = 0;
      real lse_jn = 0;
      if (cs_mixed[n] == 0) {
        lse = log_sum_exp(u);
        lse2 = lse;
      } else {
        lse = log_sum_exp(u[segment(iv_loc, iv_first[n], iv_count[n])]);
        lse2 = (ev_count[n] > 0) ? log_sum_exp(u[segment(s2_loc, s2_first[n], s2_count[n])]) : lse;
      }
      eta[k] = alpha0 + W[n] * alpha + c_i[sy_state[n]] + tau * lse;
      if (has_lead[n] == 1) {
        u_ld = util_block_c(s, e, x, mu_n + delta_lead, sigma_g[2], rho, b_j, alt_target, zg);
        lse_ld = (cs_mixed[n] == 0) ? log_sum_exp(u_ld) : log_sum_exp(u_ld[segment(s2_loc, s2_first[n], s2_count[n])]);
      }
      if (has_join[n] == 1) {
        u_jn = util_block_c(s, e, x, mu_n + delta_join, sigma_g[3], rho, b_j, alt_target, zg);
        lse_jn = (cs_mixed[n] == 0) ? log_sum_exp(u_jn) : log_sum_exp(u_jn[segment(s2_loc, s2_first[n], s2_count[n])]);
      }
      for (j in 1:ev_count[n]) {
        int ev = ev_sorted[ev_first[n] + j - 1];
        int r = ev_row[ev] - s + 1;
        if (ev_grp[ev] == 0) lp += ev_w[ev] * (u[r] - lse2);
        else if (ev_grp[ev] == 1) lp += ev_w[ev] * (u_ld[r] - lse_ld);
        else lp += ev_w[ev] * (u_jn[r] - lse_jn);
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
  int<lower=0> P;                              // 0 = no X (conventional comparison, 2026-10-31: beta is 0 x 2, mu = mu0 + shifts + a_i)
  int<lower=0> Pz;
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=S> sy_state;
  matrix[N1, P1] W;
  matrix[N1, P] X;
  array[N1] int<lower=1, upper=N2> cs_start;
  array[N1] int<lower=1, upper=N2> cs_end;
  array[N1] int<lower=1> n_full;               // unused here (sampled = 0 required); kept for the common data list
  array[N2] int<lower=1, upper=S> alt_target;
  matrix[N2, 2] x;
  matrix[N2, Pz] Z;
  array[N2] int<lower=0, upper=1> in_iv;
  array[N2] int<lower=0, upper=1> in_s2;
  array[E] int<lower=1, upper=N1> ev_sy;
  array[E] int<lower=1, upper=N2> ev_row;
  array[E] int<lower=1, upper=N2> ev_start;
  array[E] int<lower=1, upper=N2> ev_end;
  vector<lower=0, upper=1>[E] ev_w;
  array[E] int<lower=0, upper=2> ev_grp;       // stage-b group of the event: 0 bilateral, 1 coalition leader, 2 joiner
  int<lower=0, upper=1> joint;                 // must be 1
  int<lower=0, upper=1> sampled;               // must be 0
  int<lower=0, upper=1> prior_only;
  int<lower=1, upper=3> D_mode;                // must be 3
  vector[P] x_ref;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_mu0_sd;                  // sd of the normal prior on mu0 AND on delta_lead / delta_join
  int<lower=1> grainsize;
}
transformed data {
  array[N1] int sy_idx;
  array[N1] int ev_count = rep_array(0, N1);
  array[N1] int ev_first;
  array[E] int ev_sorted;
  array[N1] int has_lead = rep_array(0, N1);
  array[N1] int has_join = rep_array(0, N1);
  if (D_mode != 3) reject("m1_context: D_mode must be 3 (both axes), got ", D_mode);
  if (joint != 1) reject("m1_context: joint must be 1, got ", joint);
  if (sampled != 0) reject("m1_context: sampled alternatives are not supported (sampled must be 0)");
  for (n in 1:N1) sy_idx[n] = n;
  for (e in 1:E) {
    ev_count[ev_sy[e]] += 1;
    if (ev_grp[e] == 1) has_lead[ev_sy[e]] = 1;
    if (ev_grp[e] == 2) has_join[ev_sy[e]] = 1;
  }
  {
    int pos = 1;
    for (n in 1:N1) { ev_first[n] = pos; pos += ev_count[n]; }
    array[N1] int fill = ev_first;
    for (e in 1:E) { ev_sorted[fill[ev_sy[e]]] = e; fill[ev_sy[e]] += 1; }
  }
  array[N1] int cs_mixed = rep_array(0, N1);
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
  }
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
    if (ev_row[e] < ev_start[e] || ev_row[e] > ev_end[e]) reject("event ", e, " chosen row outside its own choice set");
    if (in_s2[ev_row[e]] != 1) reject("event ", e, ": chosen row is not in the stage-2 set (in_s2 = 0)");
    if (y1[ev_sy[e]] != 1) reject("event ", e, " belongs to a state-year with y1 == 0");
    if (ev_start[e] != cs_start[ev_sy[e]] || ev_end[e] != cs_end[ev_sy[e]]) reject("event ", e, ": with full choice sets the event set must be the state-year set");
    if (s2_count[ev_sy[e]] < 2) reject("event ", e, ": its state-year's stage-2 set is empty or holds only the chosen row");
  }
}
parameters {
  real alpha0;
  vector[P1] alpha;
  real tau;
  row_vector[2] mu0;                           // bilateral (reference) population ideal gap
  row_vector[2] delta_lead;                    // coalition-leader shift of the ideal gap (material, relational)
  row_vector[2] delta_join;                    // joiner shift
  matrix[P, 2] beta;
  array[3] vector<lower=0>[2] sigma_g;         // tolerance scales per group: 1 bilateral, 2 leader, 3 joiner
  real<lower=-1, upper=1> rho_x;               // the shared correlation of Sigma_g
  vector<lower=0>[3] sd_re;
  cholesky_factor_corr[3] L_Om;
  matrix[3, S] z_re;
  real<lower=0> sigma_b;
  vector[S] z_b;
  vector[Pz] gamma;
}
transformed parameters {
  matrix[S, 3] re = (diag_pre_multiply(sd_re, L_Om) * z_re)';
  vector[S] c_i = col(re, 1);
  matrix[S, 2] a_i = re[, 2:3];
  vector[S] b_j = sigma_b * (z_b - mean(z_b));
}
model {
  profile("priors") {
    alpha0 ~ normal(0, 2.5);
    alpha ~ normal(0, 1);
    tau ~ normal(0, 1);
    mu0 ~ normal(0, prior_mu0_sd);
    delta_lead ~ normal(0, prior_mu0_sd);
    delta_join ~ normal(0, prior_mu0_sd);
    to_vector(beta) ~ normal(0, prior_beta_sd);
    for (g in 1:3) sigma_g[g] ~ lognormal(0, 1);
    sd_re ~ normal(0, 1);
    L_Om ~ lkj_corr_cholesky(2);
    target += beta_lpdf((rho_x + 1) / 2 | 2.5, 2.5);
    to_vector(z_re) ~ std_normal();
    sigma_b ~ normal(0, 1);
    z_b ~ std_normal();
    gamma ~ normal(0, 1);
  }
  if (prior_only == 0) {
    target += reduce_sum_static(sy_partial_sum_c, sy_idx, grainsize,
                                y1, sy_state, W, X, cs_start, cs_end,
                                cs_mixed, iv_first, iv_count, iv_loc, s2_first, s2_count, s2_loc,
                                x, alt_target, Z, Pz,
                                ev_first, ev_count, ev_sorted, ev_row, ev_grp, ev_w, has_lead, has_join,
                                alpha0, alpha, tau, mu0, delta_lead, delta_join, beta, a_i, c_i, sigma_g, rho_x, b_j, gamma);
  }
}
generated quantities {
  vector[E] log_lik_stage2;
  vector[N1] log_lik_stage1;
  matrix[S, 2] mu_state;                       // bilateral ideal point at x_ref
  matrix[E, 2] mu_event;                       // the event's ideal point at its own X, with its group shift
  matrix[2, 2] Sigma;                          // bilateral Sigma_g[1]
  array[3] matrix[2, 2] Sigma_g;
  row_vector[2] mu0_lead = mu0 + delta_lead;
  row_vector[2] mu0_join = mu0 + delta_join;
  matrix[3, 3] Omega;
  corr_matrix[3] R_Om;
  {
    matrix[N1, 2] Mu = rep_matrix(mu0, N1) + X * beta + a_i[sy_state];
    vector[N2] util;
    vector[N1] lse;
    vector[N1] lse_s2;
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block_c(s, e, x, Mu[n], sigma_g[1], rho_x, b_j, alt_target, zg_block(Z, gamma, s, e, Pz));
      if (cs_mixed[n] == 0) {
        lse[n] = log_sum_exp(util[s:e]);
        lse_s2[n] = lse[n];
      } else {
        vector[e - s + 1] u = util[s:e];
        lse[n] = log_sum_exp(u[segment(iv_loc, iv_first[n], iv_count[n])]);
        lse_s2[n] = (s2_count[n] > 0) ? log_sum_exp(u[segment(s2_loc, s2_first[n], s2_count[n])]) : not_a_number();
      }
    }
    vector[N1] eta1 = alpha0 + W * alpha + c_i[sy_state] + tau * lse;
    for (n in 1:N1) log_lik_stage1[n] = bernoulli_logit_lpmf(y1[n] | eta1[n]);
    for (ev in 1:E) {
      int n = ev_sy[ev];
      int s = cs_start[n];
      int e = cs_end[n];
      if (ev_grp[ev] == 0) {
        log_lik_stage2[ev] = util[ev_row[ev]] - lse_s2[n];
        mu_event[ev] = Mu[n];
      } else {
        row_vector[2] mu_e = Mu[n] + ((ev_grp[ev] == 1) ? delta_lead : delta_join);
        vector[e - s + 1] u_g = util_block_c(s, e, x, mu_e, sigma_g[ev_grp[ev] + 1], rho_x, b_j, alt_target, zg_block(Z, gamma, s, e, Pz));
        real lse_g = (cs_mixed[n] == 0) ? log_sum_exp(u_g) : log_sum_exp(u_g[segment(s2_loc, s2_first[n], s2_count[n])]);
        log_lik_stage2[ev] = u_g[ev_row[ev] - s + 1] - lse_g;
        mu_event[ev] = mu_e;
      }
    }
    for (s in 1:S) mu_state[s] = mu0 + x_ref' * beta + a_i[s];
    for (g in 1:3) Sigma_g[g] = [[sigma_g[g][1]^2, rho_x * sigma_g[g][1] * sigma_g[g][2]],
                                 [rho_x * sigma_g[g][1] * sigma_g[g][2], sigma_g[g][2]^2]];
    Sigma = Sigma_g[1];
    R_Om = multiply_lower_tri_self_transpose(L_Om);
    Omega = quad_form_diag(R_Om, sd_re);
  }
}
