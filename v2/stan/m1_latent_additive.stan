// ============================================================================
// m1_latent_additive.stan  -- v2 strategic targeting, model family M1
//
// Two-stage spatial choice model of MID initiation, estimated JOINTLY.
//
// Stage 1 (monadic, state-year n):  y1[n] ~ Bernoulli(inv_logit(eta1[n]))
//   eta1[n] = alpha0 + W[n] alpha + c_i + tau * IV[n]
//
// Stage 2 (dyadic choice, event e made by state-year n = ev_sy[e]):
//   U[r] = -1/2 (x[r] - mu[n])' Sigma^{-1} (x[r] - mu[n]) + b_j + Z[r] gamma
//   P(chosen = r) = exp(U[r]) / sum_{r' in choice set of n} exp(U[r'])
//   mu[n] = X[n] beta + a_i                         (beta is P x 2)
//   IV[n] = log_sum_exp_{r in choice set of n} U[r]  (+ sampling correction)
//
// Random effects: (c_i, a_i1, a_i2) ~ MVN(0, Omega), Omega = diag(sd_re) R diag(sd_re)
//                 b_j = sigma_b * (z_b - mean(z_b))   (soft sum-to-zero; a constant
//                 added to every b_j cancels in the softmax so it is not identified)
//
// Flag `joint`:  1 = full Sigma (rho_x free) and full 3x3 correlation R ("b")
//                0 = Sigma diagonal AND corr(a_i1, a_i2) = 0 ("a"); corr(c, a1) and
//                    corr(c, a2) remain free, parameterised so R is always PD:
//                    r1 = corr(c,a1) in (-1,1), r2 = q * sqrt(1 - r1^2), q in (-1,1)
//                    => det(R) = 1 - r1^2 - r2^2 > 0.
// Flag `sampled`: 0 = FULL choice sets (main specification; the event's set IS its
//                 state-year's set, rows cs_start..cs_end).
//                 1 = two kinds of sampled sets, both built in R (sample_alternatives):
//                 (i) IV rows: for EVERY state-year, m alternatives drawn uniformly without
//                     replacement from the full set, independently of outcomes (the same rule
//                     for initiating and non-initiating rows, chosen targets NOT forced in);
//                     IV = log_sum_exp over them + log(n_full / n_sub). This estimate of the
//                     log partition sum is downward biased by Jensen's inequality
//                     (E log S_hat <= log E S_hat = log S); the bias shrinks with m and is
//                     identical in expectation for initiating and non-initiating rows.
//                 (ii) event rows: EACH event gets its own set = chosen target + m alternatives
//                     drawn uniformly from the rest, independently per event, stored after the
//                     IV rows (ev_start..ev_end). Under uniform sampling the conditional
//                     logit on the event's own set has the positive conditioning property
//                     (McFadden 1978) so its score is unbiased.
//
// The normal-density constant -log(2 pi |Sigma|^{1/2}) is constant across
// alternatives and cancels in the softmax, so only the quadratic form enters.
// The term -1/2 mu' Sigma^{-1} mu is also constant across alternatives; it is kept
// here (cheap) so that IV has the full systematic-utility level.
//
// Ragged choice sets are flattened: rows cs_start[n]..cs_end[n] belong to state-year n.
// ============================================================================
functions {
  // Utility of every alternative row for state-year n, given mu[n] (length-2 row).
  vector util_block(int n, int s, int e, matrix x, row_vector mu_n, vector sigma_x,
                    real rho, vector b_j, array[] int alt_target, vector zg) {
    int L = e - s + 1;
    vector[L] d1 = (segment(col(x, 1), s, L) - mu_n[1]) / sigma_x[1];
    vector[L] d2 = (segment(col(x, 2), s, L) - mu_n[2]) / sigma_x[2];
    vector[L] q = (d1 .* d1 - 2 * rho * (d1 .* d2) + d2 .* d2) / (1 - rho * rho);
    return -0.5 * q + b_j[segment(alt_target, s, L)] + segment(zg, s, L);
  }
}
data {
  // ---- dimensions
  int<lower=1> S;             // number of states (index space for c_i, a_i, b_j)
  int<lower=1> N1;            // state-years (stage-1 rows; each owns a choice set)
  int<lower=1> N2;            // flattened alternative rows
  int<lower=1> E;             // initiation events
  int<lower=0> P1;            // stage-1 controls W (no intercept column)
  int<lower=1> P;             // predictors of mu (initiator power changes + controls)
  int<lower=0> Pz;            // dyadic controls Z
  // ---- stage 1 / state-year level
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=S> sy_state;
  matrix[N1, P1] W;
  matrix[N1, P] X;
  array[N1] int<lower=1, upper=N2> cs_start;
  array[N1] int<lower=1, upper=N2> cs_end;
  array[N1] int<lower=1> n_full;        // full choice-set size (for sampling correction)
  // ---- alternative rows
  array[N2] int<lower=1, upper=S> alt_target;
  matrix[N2, 2] x;                      // (dM_ij, dR_ij) = (M_i - M_j, R_i - R_j) at t-1
  matrix[N2, Pz] Z;
  // ---- events
  array[E] int<lower=1, upper=N1> ev_sy;
  array[E] int<lower=1, upper=N2> ev_row;      // chosen alternative's row
  array[E] int<lower=1, upper=N2> ev_start;    // the event's OWN choice-set rows
  array[E] int<lower=1, upper=N2> ev_end;
  // ---- flags and reference point
  int<lower=0, upper=1> joint;
  int<lower=0, upper=1> sampled;
  vector[P] x_ref;                      // reference X at which mu_state is reported
  real<lower=0> prior_beta_sd;
}
transformed data {
  vector[N1] iv_corr = rep_vector(0, N1);
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
}
parameters {
  real alpha0;
  vector[P1] alpha;
  real tau;
  matrix[P, 2] beta;
  vector<lower=0>[2] sigma_x;
  array[joint] real<lower=-1, upper=1> rho_x;
  vector<lower=0>[3] sd_re;                      // sd of (c, a1, a2)
  array[joint] cholesky_factor_corr[3] L_Om_joint;
  array[1 - joint] real<lower=-1, upper=1> r_ca1;
  array[1 - joint] real<lower=-1, upper=1> q_ca2;
  matrix[3, S] z_re;
  real<lower=0> sigma_b;
  vector[S] z_b;
  vector[Pz] gamma;
}
transformed parameters {
  matrix[3, 3] L_Om;
  real rho;
  if (joint == 1) {
    L_Om = L_Om_joint[1];
    rho = rho_x[1];
  } else {
    real r1 = r_ca1[1];
    real r2 = q_ca2[1] * sqrt(1 - r1 * r1);
    matrix[3, 3] R_a = [[1, r1, r2], [r1, 1, 0], [r2, 0, 1]];
    L_Om = cholesky_decompose(R_a);
    rho = 0;
  }
  matrix[S, 3] re = (diag_pre_multiply(sd_re, L_Om) * z_re)';
  vector[S] c_i = col(re, 1);
  matrix[S, 2] a_i = re[, 2:3];
  vector[S] b_j = sigma_b * (z_b - mean(z_b));
}
model {
  // ---- priors
  alpha0 ~ normal(0, 2.5);
  alpha ~ normal(0, 1);
  tau ~ normal(0, 1);
  to_vector(beta) ~ normal(0, prior_beta_sd);
  sigma_x ~ lognormal(0, 1);
  sd_re ~ normal(0, 1);
  if (joint == 1) L_Om_joint[1] ~ lkj_corr_cholesky(2);
  // r_ca1, q_ca2: implicit uniform(-1, 1)
  to_vector(z_re) ~ std_normal();
  sigma_b ~ normal(0, 1);
  z_b ~ std_normal();
  gamma ~ normal(0, 1);

  // ---- likelihood
  {
    matrix[N1, 2] Mu = X * beta + a_i[sy_state];
    vector[N2] zg = (Pz > 0) ? Z * gamma : rep_vector(0, N2);
    vector[N2] util;
    vector[N1] lse;
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block(n, s, e, x, Mu[n], sigma_x, rho, b_j, alt_target, zg);
      lse[n] = log_sum_exp(util[s:e]);
    }
    y1 ~ bernoulli_logit(alpha0 + W * alpha + c_i[sy_state] + tau * (lse + iv_corr));
    for (ev in 1:E) {
      int s = ev_start[ev];
      int e = ev_end[ev];
      if (sampled == 1)   // event rows are separate from the IV rows: compute their utilities
        util[s:e] = util_block(ev_sy[ev], s, e, x, Mu[ev_sy[ev]], sigma_x, rho, b_j, alt_target, zg);
      target += util[ev_row[ev]] - log_sum_exp(util[s:e]);
    }
  }
}
generated quantities {
  vector[E] log_lik_stage2;
  vector[N1] log_lik_stage1;
  matrix[S, 2] mu_state;      // mu_i at x_ref: x_ref' beta + a_i
  matrix[E, 2] mu_event;      // mu at the event's own X: X[ev_sy] beta + a_i
  matrix[2, 2] Sigma;
  matrix[3, 3] Omega;         // covariance of (c, a1, a2)
  corr_matrix[3] R_Om;        // correlation of (c, a1, a2)
  {
    matrix[N1, 2] Mu = X * beta + a_i[sy_state];
    vector[N2] zg = (Pz > 0) ? Z * gamma : rep_vector(0, N2);
    vector[N2] util;
    vector[N1] lse;
    for (n in 1:N1) {
      int s = cs_start[n];
      int e = cs_end[n];
      util[s:e] = util_block(n, s, e, x, Mu[n], sigma_x, rho, b_j, alt_target, zg);
      lse[n] = log_sum_exp(util[s:e]);
    }
    vector[N1] eta1 = alpha0 + W * alpha + c_i[sy_state] + tau * (lse + iv_corr);
    for (n in 1:N1) log_lik_stage1[n] = bernoulli_logit_lpmf(y1[n] | eta1[n]);
    for (ev in 1:E) {
      int s = ev_start[ev];
      int e = ev_end[ev];
      if (sampled == 1)
        util[s:e] = util_block(ev_sy[ev], s, e, x, Mu[ev_sy[ev]], sigma_x, rho, b_j, alt_target, zg);
      log_lik_stage2[ev] = util[ev_row[ev]] - log_sum_exp(util[s:e]);
      mu_event[ev] = Mu[ev_sy[ev]];
    }
    for (s in 1:S) mu_state[s] = x_ref' * beta + a_i[s];
    Sigma = [[sigma_x[1]^2, rho * sigma_x[1] * sigma_x[2]],
             [rho * sigma_x[1] * sigma_x[2], sigma_x[2]^2]];
    R_Om = multiply_lower_tri_self_transpose(L_Om);
    Omega = quad_form_diag(R_Om, sd_re);
  }
}
