// ============================================================================
// m2_context.stan  -- v2 strategic targeting, model M2c: AME with initiation types (PI-approved 2026-11-04)
//
// m2_ame.stan (M2b, the joint "b" variant: D = 2, full Sigma, full Omega) PLUS the role structure of m1_context.stan (M1c):
//   U_g[r] = -1/2 (x[r] - mu_e)' Sigma_g^{-1} (x[r] - mu_e) + b_j + Z[r] gamma + u_i' v_j
//   mu_e   = mu0 + delta[g] * [g >= 1] + X[n] beta + a_i,   Sigma_g = diag(sigma_g[g + 1]) R diag(sigma_g[g + 1]),  R = [[1, rho], [rho, 1]] shared
//   group g = ev_grp[e] in 0..G-1 (data G, ev_grp; group 0 = the reference: no shift, sigma_g[1]); the stage-2 softmax over the event's allowed rows
//   IV[n]  = log_sum_exp over the in_iv rows of the REFERENCE-group utility U_0 (mu0 + X[n] beta + a_i, sigma_g[1], rho, AND u_i' v_j), exactly as
//            M1c uses the reference utility and M2b puts the latent term in the IV
// It is derived from m1_context.stan (same data, same transformed data, same parameter names, so an M1c fit initialises it) with the AME block of
// m2_ame.stan added: data K and anchor; parameters u_diag, u_lower, U_free, V_raw; U and V transformed parameters with the M2b identification
// (V column-centred; U column-centred by a hard constraint on the last non-anchor row; the K anchor rows of U lower-triangular with a positive
// diagonal; scale fixed by the N(0, 1) priors; Procrustes alignment across chains in 03_extract.R). The AME term enters every group's utility
// identically (V[tg] * u_i'). With V = 0 (V_raw constant) the program is M1c term for term (test_m2_context.R).
// Priors: M1c's for the role parameters (delta[g] ~ normal(0, prior_mu0_sd), sigma_g[g] ~ lognormal(0, 1), (rho + 1)/2 ~ Beta(2.5, 2.5)) and
// M2b's for everything else (identical to M1c's for the shared ones; u_diag, u_lower, U[nonanchor], V_raw ~ normal(0, 1)).
// Event-level choice sets (2026-11-04, switch joiner_set): as m1_context.stan (n_ev_alt, ev_alt_first, ev_alt_n, ev_alt_loc).
// Restrictions: D_mode = 3, joint = 1, sampled = 0 (rejected otherwise); K in 1..S-1 with distinct anchors.
// Generated quantities: log_lik_stage1[N1], log_lik_stage2[E] (unweighted, the event's own group and allowed set), mu_state[S, 2] (reference
// group at x_ref), mu_event[E, 2] (with the group shift), Sigma (reference group), Sigma_g[G], mu0_g[G], Omega, R_Om; U, V (transformed parameters).
// ============================================================================
functions {
  // Utility of the alternative rows s..e for one state-year and one group: quadratic form in (dM, dR) around mu_n with scales
  // sig (length 2) and the shared correlation rho, plus the target intercept and the dyadic controls.
  // ... plus the AME term V[tg] u_i' (the order of the terms is M2b's)
  vector util_block_c(int s, int e, matrix x, row_vector mu_n, vector sig, real rho, vector b_j, array[] int alt_target, vector zg_local,
                      row_vector u_i, matrix V) {
    int L = e - s + 1;
    vector[L] d1 = (segment(col(x, 1), s, L) - mu_n[1]) / sig[1];
    vector[L] d2 = (segment(col(x, 2), s, L) - mu_n[2]) / sig[2];
    vector[L] q = (d1 .* d1 - 2 * rho * (d1 .* d2) + d2 .* d2) / (1 - rho * rho);
    array[L] int tg = segment(alt_target, s, L);
    return -0.5 * q + b_j[tg] + zg_local + V[tg] * u_i';
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
                        array[] int ev_row, array[] int ev_grp, vector ev_w, int G, array[,] int has_g,
                        array[] int ev_restr, array[] int ev_alt_first, array[] int ev_alt_n, array[] int ev_alt_loc,
                        real alpha0, vector alpha, real tau, row_vector mu0, array[] row_vector delta,
                        matrix beta, matrix a_i, vector c_i, array[] vector sigma_g, real rho, vector b_j, vector gamma,
                        matrix U, matrix V) {
    int L = end - start + 1;
    vector[L] eta;
    real lp = 0;
    for (k in 1:L) {
      int n = sy_slice[k];
      int s = cs_start[n];
      int e = cs_end[n];
      row_vector[2] mu_n = mu0 + X[n] * beta + a_i[sy_state[n]];
      vector[e - s + 1] zg = zg_block(Z, gamma, s, e, Pz);
      array[G] vector[e - s + 1] u_g;          // u_g[1]: reference group (IV and group-0 events); u_g[g + 1]: group g, filled when present
      vector[G] lse_g = rep_vector(0, G);      // lse_g[1] = the stage-2 lse of the reference group (lse2 before 2026-10-31)
      real lse;
      u_g[1] = util_block_c(s, e, x, mu_n, sigma_g[1], rho, b_j, alt_target, zg, U[sy_state[n]], V);
      if (cs_mixed[n] == 0) {
        lse = log_sum_exp(u_g[1]);
        lse_g[1] = lse;
      } else {
        lse = log_sum_exp(u_g[1][segment(iv_loc, iv_first[n], iv_count[n])]);
        lse_g[1] = (ev_count[n] > 0) ? log_sum_exp(u_g[1][segment(s2_loc, s2_first[n], s2_count[n])]) : lse;
      }
      eta[k] = alpha0 + W[n] * alpha + c_i[sy_state[n]] + tau * lse;
      for (g in 2:G) {
        if (has_g[n, g] == 1) {
          u_g[g] = util_block_c(s, e, x, mu_n + delta[g - 1], sigma_g[g], rho, b_j, alt_target, zg, U[sy_state[n]], V);
          lse_g[g] = (cs_mixed[n] == 0) ? log_sum_exp(u_g[g]) : log_sum_exp(u_g[g][segment(s2_loc, s2_first[n], s2_count[n])]);
        }
      }
      for (j in 1:ev_count[n]) {
        int ev = ev_sorted[ev_first[n] + j - 1];
        int r = ev_row[ev] - s + 1;
        int gg = ev_grp[ev] + 1;
        if (ev_restr[ev] == 0) {
          lp += ev_w[ev] * (u_g[gg][r] - lse_g[gg]);
        } else {
          lp += ev_w[ev] * (u_g[gg][r] - log_sum_exp(u_g[gg][segment(ev_alt_loc, ev_alt_first[ev], ev_alt_n[ev])]));   // restricted set (2026-11-04)
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
  int<lower=2> G;                              // number of stage-2 groups (data since 2026-10-31; 3 = bilateral / coalition leader / joiner)
  array[E] int<lower=0, upper=G - 1> ev_grp;   // stage-b group of the event, 0..G-1 (0 = the reference group)
  int<lower=0> n_ev_alt;                       // 2026-11-04: event-level allowed stage-2 rows (ragged): total length
  array[E] int<lower=1> ev_alt_first;          //   first position of event e's rows in ev_alt_loc
  array[E] int<lower=1> ev_alt_n;              //   number of allowed rows of event e (its whole stage-2 set unless restricted)
  array[n_ev_alt] int<lower=1> ev_alt_loc;     //   block-local positions (r - cs_start + 1), strictly increasing within an event
  int<lower=0, upper=1> joint;                 // must be 1
  int<lower=0, upper=1> sampled;               // must be 0
  int<lower=0, upper=1> prior_only;
  int<lower=1, upper=3> D_mode;                // must be 3
  vector[P] x_ref;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_mu0_sd;                  // sd of the normal prior on mu0 AND on delta[g]
  int<lower=1> grainsize;
  // ---- AME (as m2_ame.stan)
  int<lower=1> K;
  array[K] int<lower=1, upper=S> anchor;      // distinct state indices (the K states with most initiations)
}
transformed data {
  array[N1] int sy_idx;
  array[N1] int ev_count = rep_array(0, N1);
  array[N1] int ev_first;
  array[E] int ev_sorted;
  array[N1, G] int has_g = rep_array(0, N1, G);   // has_g[n, g + 1] = 1 when state-year n holds an event of group g
  if (D_mode != 3) reject("m2_context: D_mode must be 3 (both axes), got ", D_mode);
  if (joint != 1) reject("m2_context: joint must be 1, got ", joint);
  if (sampled != 0) reject("m2_context: sampled alternatives are not supported (sampled must be 0)");
  for (n in 1:N1) sy_idx[n] = n;
  for (e in 1:E) {
    ev_count[ev_sy[e]] += 1;
    has_g[ev_sy[e], ev_grp[e] + 1] = 1;
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
  // event-level allowed sets (2026-11-04): within the state-year's stage-2 set, strictly increasing, holding the chosen row; restricted when smaller
  array[E] int ev_restr = rep_array(0, E);
  for (e in 1:E) {
    int n = ev_sy[e];
    int L = cs_end[n] - cs_start[n] + 1;
    int found = 0;
    if (ev_alt_first[e] + ev_alt_n[e] - 1 > n_ev_alt) reject("event ", e, ": its allowed-row index runs past ev_alt_loc");
    if (ev_alt_n[e] > s2_count[n]) reject("event ", e, ": more allowed rows than its state-year's stage-2 set");
    for (k in 1:ev_alt_n[e]) {
      int l = ev_alt_loc[ev_alt_first[e] + k - 1];
      if (l < 1 || l > L) reject("event ", e, ": allowed row outside its state-year's block");
      if (in_s2[cs_start[n] + l - 1] != 1) reject("event ", e, ": allowed row not in the stage-2 set");
      if (k > 1 && l <= ev_alt_loc[ev_alt_first[e] + k - 2]) reject("event ", e, ": allowed rows must be strictly increasing");
      if (cs_start[n] + l - 1 == ev_row[e]) found = 1;
    }
    if (found == 0) reject("event ", e, ": the chosen row is not among its allowed rows");
    if (ev_alt_n[e] < s2_count[n]) ev_restr[e] = 1;
  }
  // ---- AME anchors (as m2_ame.stan)
  array[S - K] int nonanchor;
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
  row_vector[2] mu0;                           // bilateral (reference) population ideal gap
  array[G - 1] row_vector[2] delta;            // shift of the ideal gap of group g = 1..G-1 vs the reference (material, relational)
  matrix[P, 2] beta;
  array[G] vector<lower=0>[2] sigma_g;         // tolerance scales per group: sigma_g[g + 1] for group g (1 = the reference)
  real<lower=-1, upper=1> rho_x;               // the shared correlation of Sigma_g
  vector<lower=0>[3] sd_re;
  cholesky_factor_corr[3] L_Om;
  matrix[3, S] z_re;
  real<lower=0> sigma_b;
  vector[S] z_b;
  vector[Pz] gamma;
  // ---- AME (as m2_ame.stan)
  vector<lower=0>[K] u_diag;                   // positive diagonal of the anchor block
  vector[K * (K - 1) / 2] u_lower;             // strictly lower-triangular anchor entries
  matrix[S - K - 1, K] U_free;                 // rows for non-anchor states except the last (determined by centring)
  matrix[S, K] V_raw;
}
transformed parameters {
  matrix[S, 3] re = (diag_pre_multiply(sd_re, L_Om) * z_re)';
  vector[S] c_i = col(re, 1);
  matrix[S, 2] a_i = re[, 2:3];
  vector[S] b_j = sigma_b * (z_b - mean(z_b));
  // ---- assemble U with the anchor block lower-triangular, positive diagonal; U and V column-centred (as m2_ame.stan)
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
    for (g in 1:(G - 1)) delta[g] ~ normal(0, prior_mu0_sd);
    to_vector(beta) ~ normal(0, prior_beta_sd);
    for (g in 1:G) sigma_g[g] ~ lognormal(0, 1);
    sd_re ~ normal(0, 1);
    L_Om ~ lkj_corr_cholesky(2);
    target += beta_lpdf((rho_x + 1) / 2 | 2.5, 2.5);
    to_vector(z_re) ~ std_normal();
    sigma_b ~ normal(0, 1);
    z_b ~ std_normal();
    gamma ~ normal(0, 1);
    u_diag ~ normal(0, 1);
    u_lower ~ normal(0, 1);
    to_vector(U[nonanchor]) ~ normal(0, 1);  // all non-anchor rows, including the determined one
    to_vector(V_raw) ~ normal(0, 1);
  }
  if (prior_only == 0) {
    target += reduce_sum_static(sy_partial_sum_c, sy_idx, grainsize,
                                y1, sy_state, W, X, cs_start, cs_end,
                                cs_mixed, iv_first, iv_count, iv_loc, s2_first, s2_count, s2_loc,
                                x, alt_target, Z, Pz,
                                ev_first, ev_count, ev_sorted, ev_row, ev_grp, ev_w, G, has_g,
                                ev_restr, ev_alt_first, ev_alt_n, ev_alt_loc,
                                alpha0, alpha, tau, mu0, delta, beta, a_i, c_i, sigma_g, rho_x, b_j, gamma, U, V);
  }
}
generated quantities {
  vector[E] log_lik_stage2;
  vector[N1] log_lik_stage1;
  matrix[S, 2] mu_state;                       // reference-group ideal point at x_ref
  matrix[E, 2] mu_event;                       // the event's ideal point at its own X, with its group shift
  matrix[2, 2] Sigma;                          // reference-group Sigma_g[1]
  array[G] matrix[2, 2] Sigma_g;
  array[G] row_vector[2] mu0_g;                // group population ideal points: mu0_g[1] = mu0, mu0_g[g + 1] = mu0 + delta[g]
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
      util[s:e] = util_block_c(s, e, x, Mu[n], sigma_g[1], rho_x, b_j, alt_target, zg_block(Z, gamma, s, e, Pz), U[sy_state[n]], V);
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
        if (ev_restr[ev] == 0) {
          log_lik_stage2[ev] = util[ev_row[ev]] - lse_s2[n];
        } else {
          vector[e - s + 1] u0 = util[s:e];
          log_lik_stage2[ev] = util[ev_row[ev]] - log_sum_exp(u0[segment(ev_alt_loc, ev_alt_first[ev], ev_alt_n[ev])]);
        }
        mu_event[ev] = Mu[n];
      } else {
        row_vector[2] mu_e = Mu[n] + delta[ev_grp[ev]];
        vector[e - s + 1] u_g = util_block_c(s, e, x, mu_e, sigma_g[ev_grp[ev] + 1], rho_x, b_j, alt_target, zg_block(Z, gamma, s, e, Pz), U[sy_state[n]], V);
        real lse_g = (ev_restr[ev] == 1) ? log_sum_exp(u_g[segment(ev_alt_loc, ev_alt_first[ev], ev_alt_n[ev])])
                     : ((cs_mixed[n] == 0) ? log_sum_exp(u_g) : log_sum_exp(u_g[segment(s2_loc, s2_first[n], s2_count[n])]));
        log_lik_stage2[ev] = u_g[ev_row[ev] - s + 1] - lse_g;
        mu_event[ev] = mu_e;
      }
    }
    for (s in 1:S) mu_state[s] = mu0 + x_ref' * beta + a_i[s];
    mu0_g[1] = mu0;
    for (g in 2:G) mu0_g[g] = mu0 + delta[g - 1];
    for (g in 1:G) Sigma_g[g] = [[sigma_g[g][1]^2, rho_x * sigma_g[g][1] * sigma_g[g][2]],
                                 [rho_x * sigma_g[g][1] * sigma_g[g][2], sigma_g[g][2]^2]];
    Sigma = Sigma_g[1];
    R_Om = multiply_lower_tri_self_transpose(L_Om);
    Omega = quad_form_diag(R_Om, sd_re);
  }
}
