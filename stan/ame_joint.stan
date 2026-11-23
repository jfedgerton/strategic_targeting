// stan/ame_joint.stan -- joint two-stage AME (docs/ame_design.md rev 2), cells 1-3 (PI 2026-11-23: K = 0 main, cut posterior, K = 1).
// Stage 1: conditional logit over each event's alternatives; utility = X1 * beta1 + b_j + eta_j,era + coal * delta_j [+ lam1 * u_i * v_j].
// Stage 2: masked categorical {out (ref), A, B, O}; utility = alpha_c + X_c * beta_c + theta_c * bbar_d + eps_dc [+ lam2 * u_k * vbar_opp].
// bbar_d = mean over the dispute's onset target set of (b_j + eta_j,era(d)), from the draws (or data for the cut posterior).
// Threading: each reduce_sum slice computes its own utilities from data rows and SMALL shared parameter vectors (no large autodiff
// vectors are passed as shared arguments, which reduce_sum would deep-copy per slice).
functions {
  // stage 1 over a slice of events (events are contiguous blocks of rows)
  real s1_partial(array[] int ev_slice, int start, int end, matrix X1, vector beta1, array[] int tgt1, array[] int era1, vector coal1,
                  array[] int ini1, vector b, matrix eta, vector del, int K, vector uu, vector vv, real lam1,
                  array[] int e_start, array[] int e_end, array[] int e_chosen, vector w_e) {
    real lp = 0;
    for (k in 1:size(ev_slice)) { int e = ev_slice[k]; int s = e_start[e]; int t = e_end[e]; int n = t - s + 1;
      vector[n] u = X1[s:t] * beta1;
      for (r in 1:n) { int row = s + r - 1; u[r] += b[tgt1[row]] + eta[tgt1[row], era1[row]] + coal1[row] * del[tgt1[row]];
        if (K) u[r] += lam1 * uu[ini1[row]] * vv[tgt1[row]]; }
      lp += w_e[e] * (u[e_chosen[e] - s + 1] - log_sum_exp(u)); }
    return lp;
  }
  // stage 2 over a slice of candidate rows
  real s2_partial(array[] int r_slice, int start, int end, matrix XA, matrix XB, matrix XO, matrix alpha, array[] int era_d, vector betaA, vector betaB, vector betaO,
                  vector theta, vector bbar, matrix eps, array[] int d2, array[] int y2, array[,] int av2,
                  int K, vector uu, vector vA, vector vB, real lam2, array[] int k2) {
    int n = end - start + 1;
    vector[n] UA = XA[start:end] * betaA; vector[n] UB = XB[start:end] * betaB; vector[n] UO = XO[start:end] * betaO;
    real lp = 0;
    for (r in 1:n) { int row = start + r - 1; int d = d2[row]; vector[4] Um; Um[1] = 0;
      real la = 0; real lb = 0;
      if (K) { la = lam2 * uu[k2[row]] * vA[d]; lb = lam2 * uu[k2[row]] * vB[d]; }
      int ed = era_d[d];
      Um[2] = av2[row, 1] ? alpha[1, ed] + UA[r] + theta[1] * bbar[d] + eps[d, 1] + la : negative_infinity();
      Um[3] = av2[row, 2] ? alpha[2, ed] + UB[r] + theta[2] * bbar[d] + eps[d, 2] + lb : negative_infinity();
      Um[4] = av2[row, 3] ? alpha[3, ed] + UO[r] + theta[3] * bbar[d] + eps[d, 3] + la : negative_infinity();
      lp += Um[y2[row]] - log_sum_exp(Um); }
    return lp;
  }
}
data {
  int<lower=0, upper=1> use_s1; int<lower=0, upper=1> use_s2; int<lower=0, upper=1> fix_bbar; int<lower=0, upper=1> K;
  int<lower=1> J; int<lower=1> NE_ERA;
  int<lower=0> N1; int<lower=0> P1; matrix[N1, P1] X1; array[N1] int<lower=1, upper=J> tgt1; array[N1] int<lower=1, upper=NE_ERA> era1;
  vector<lower=0, upper=1>[N1] coal1; array[N1] int<lower=1, upper=J> ini1;
  int<lower=0> E; array[E] int e_start; array[E] int e_end; array[E] int e_chosen; vector<lower=0>[E] w_e;
  int<lower=0> N2; int<lower=0> PA; int<lower=0> PB; int<lower=0> PO; matrix[N2, PA] XA; matrix[N2, PB] XB; matrix[N2, PO] XO;
  array[N2] int<lower=1, upper=4> y2; array[N2, 3] int<lower=0, upper=1> av2; array[N2] int<lower=1, upper=J> k2;
  int<lower=0> ND; array[N2] int<lower=1, upper=ND> d2;
  int<lower=0> NT; array[NT] int<lower=1, upper=J> ts_j; array[NT] int<lower=1, upper=ND> ts_d; array[ND] int<lower=1, upper=NE_ERA> era_d;
  vector<lower=0>[ND] ts_n; vector[ND] bbar_data;
  matrix[ND, J] W_oppA; matrix[ND, J] W_oppB; matrix[J, J] Pv; vector<lower=0, upper=1>[J] u_set;   // u_set: states appearing as initiators or stage-2 candidates (u centred over these)
  int<lower=1> grainsize; int<lower=1> grainsize2;
}
transformed data { array[E] int ev_idx; array[N2] int r_idx; for (e in 1:E) ev_idx[e] = e; for (n in 1:N2) r_idx[n] = n; }
parameters {
  vector[P1] beta1;
  real<lower=0> sigma_b; real<lower=0> sigma_del; real<lower=0> sigma_eta;
  vector[J] b_raw; vector[J] del_raw; matrix[J, NE_ERA] eta_raw;
  matrix[3, NE_ERA] alpha; vector[PA] betaA; vector[PB] betaB; vector[PO] betaO; vector[3] theta;
  vector<lower=0>[3] sigma_eps; matrix[ND, 3] eps_raw;
  vector[K ? J : 0] u_free; vector[K ? J : 0] v_free; array[K] real<lower=0> lam1_raw; array[K] real lam2_raw;   // lambda1 > 0 fixes the sign; lambda2 free
}
transformed parameters {
  // centred (audit P1-1): a constant added to b, to an era column of eta, or to del cancels inside every stage-1 event, so only the prior
  // would pin it; uncentred it would pass into bbar and let theta act as era intercepts in stage 2. Stage 2 has its own era intercepts.
  vector[J] b = sigma_b * (b_raw - mean(b_raw));
  vector[J] del = sigma_del * (del_raw - mean(del_raw));
  matrix[J, NE_ERA] eta;
  for (e in 1:NE_ERA) eta[, e] = sigma_eta * (eta_raw[, e] - mean(eta_raw[, e]));
  vector[ND] bbar;
  if (fix_bbar) bbar = bbar_data;
  else { bbar = rep_vector(0, ND); for (t in 1:NT) bbar[ts_d[t]] += b[ts_j[t]] + eta[ts_j[t], era_d[ts_d[t]]]; bbar = bbar ./ ts_n; }
}
model {
  vector[J] uu = rep_vector(0, J); vector[J] vv = rep_vector(0, J); real lam1 = 0; real lam2 = 0;
  beta1 ~ normal(0, 2.5); to_vector(alpha) ~ normal(0, 2.5); betaA ~ normal(0, 2.5); betaB ~ normal(0, 2.5); betaO ~ normal(0, 2.5); theta ~ normal(0, 2.5);
  sigma_b ~ normal(0, 1); sigma_del ~ normal(0, 1); sigma_eta ~ normal(0, 1); sigma_eps ~ normal(0, 1);
  b_raw ~ student_t(4, 0, 1); del_raw ~ std_normal(); to_vector(eta_raw) ~ std_normal(); to_vector(eps_raw) ~ std_normal();
  if (K) { u_free ~ std_normal(); v_free ~ std_normal(); lam1_raw ~ normal(0, 1); lam2_raw ~ normal(0, 1);
    uu = (u_free - dot_product(u_free, u_set) / sum(u_set)) .* u_set; vv = Pv * v_free; lam1 = lam1_raw[1]; lam2 = lam2_raw[1]; }
  if (use_s1) target += reduce_sum(s1_partial, ev_idx, grainsize, X1, beta1, tgt1, era1, coal1, ini1, b, eta, del, K, uu, vv, lam1, e_start, e_end, e_chosen, w_e);
  if (use_s2) {
    matrix[ND, 3] eps = eps_raw * diag_matrix(sigma_eps);
    vector[ND] vA = rep_vector(0, ND); vector[ND] vB = rep_vector(0, ND);
    if (K) { vA = W_oppA * vv; vB = W_oppB * vv; }
    target += reduce_sum(s2_partial, r_idx, grainsize2, XA, XB, XO, alpha, era_d, betaA, betaB, betaO, theta, bbar, eps, d2, y2, av2, K, uu, vA, vB, lam2, k2);
  }
}
generated quantities {
  vector[K ? J : 0] u_out; vector[K ? J : 0] v_out;
  if (K) { u_out = (u_free - dot_product(u_free, u_set) / sum(u_set)) .* u_set; v_out = Pv * v_free; }
}
