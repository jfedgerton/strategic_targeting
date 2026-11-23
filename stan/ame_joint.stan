// stan/ame_joint.stan -- joint two-stage AME (docs/ame_design.md rev 2), cells 1-3 (PI 2026-11-23: K = 0 main, cut posterior, K = 1).
// Stage 1: conditional logit over each event's alternatives; utility = X1 * beta1 + b_j + era_j + coal * delta_j [+ lam1 * u_i * v_j].
// Stage 2: masked categorical {out (ref), A, B, O}; utility = alpha_c + X_c * beta_c + theta_c * bbar_d + eps_dc [+ lam2 * u_k * vbar_opp].
// bbar_d = mean over the dispute's onset target set of (b_j + era_j(d)), computed from the draws (or supplied as data for the cut posterior).
// Flags: use_s1, use_s2 (which likelihoods enter), fix_bbar (cut posterior: bbar is data), K in {0, 1}.
functions {
  real s1_partial(array[] int ev_slice, int start, int end, vector eta1, array[] int e_start, array[] int e_end, array[] int e_chosen, vector w_e) {
    real lp = 0;
    for (k in 1:size(ev_slice)) { int e = ev_slice[k];
      lp += w_e[e] * (eta1[e_chosen[e]] - log_sum_exp(eta1[e_start[e]:e_end[e]])); }
    return lp;
  }
  // stage 2 over a slice of candidate rows: masked 4-category log likelihood (U2 columns: A, B, O; out = 0)
  real s2_partial(array[] int r_slice, int start, int end, matrix U2, array[] int y2, array[,] int av2) {
    real lp = 0;
    for (k in 1:size(r_slice)) { int n = r_slice[k]; vector[4] Um; Um[1] = 0;
      for (c in 1:3) Um[c + 1] = av2[n, c] ? U2[n, c] : negative_infinity();
      lp += Um[y2[n]] - log_sum_exp(Um); }
    return lp;
  }
}
data {
  int<lower=0, upper=1> use_s1; int<lower=0, upper=1> use_s2; int<lower=0, upper=1> fix_bbar; int<lower=0, upper=1> K;
  int<lower=1> J; int<lower=1> NE_ERA;                       // states, eras
  // stage 1
  int<lower=0> N1; int<lower=0> P1; matrix[N1, P1] X1; array[N1] int<lower=1, upper=J> tgt1; array[N1] int<lower=1, upper=NE_ERA> era1;
  vector<lower=0, upper=1>[N1] coal1; array[N1] int<lower=1, upper=J> ini1;
  int<lower=0> E; array[E] int e_start; array[E] int e_end; array[E] int e_chosen; vector<lower=0>[E] w_e;
  // stage 2
  int<lower=0> N2; int<lower=0> PA; int<lower=0> PB; int<lower=0> PO; matrix[N2, PA] XA; matrix[N2, PB] XB; matrix[N2, PO] XO;
  array[N2] int<lower=1, upper=4> y2;                       // 1 out, 2 A, 3 B, 4 O
  array[N2, 3] int<lower=0, upper=1> av2;                   // availability of A, B, O
  array[N2] int<lower=1, upper=J> k2;                        // candidate state
  int<lower=0> ND; array[N2] int<lower=1, upper=ND> d2;      // dispute index
  int<lower=0> NT; array[NT] int<lower=1, upper=J> ts_j; array[NT] int<lower=1, upper=ND> ts_d; array[ND] int<lower=1, upper=NE_ERA> era_d;
  vector<lower=0>[ND] ts_n;                                  // targets per dispute
  vector[ND] bbar_data;                                      // used only when fix_bbar = 1
  // latent (K = 1): side-mean opponent positions as weight matrices over states, and the projection that orthogonalizes v
  matrix[ND, J] W_oppA; matrix[ND, J] W_oppB;                // rows: dispute; weights 1/|side| on the B side (for A, O) and A side (for B)
  matrix[J, J] Pv;                                           // I - H, H = projection on [1, Rbar_s, Mbar_s]
  int<lower=1> grainsize; int<lower=1> grainsize2;
}
transformed data { array[E] int ev_idx; array[N2] int r_idx; for (e in 1:E) ev_idx[e] = e; for (n in 1:N2) r_idx[n] = n; }
parameters {
  vector[P1] beta1;
  real<lower=0> sigma_b; real<lower=0> sigma_del; real<lower=0> sigma_eta;
  vector[J] b_raw; vector[J] del_raw; matrix[J, NE_ERA] eta_raw;
  vector[3] alpha; vector[PA] betaA; vector[PB] betaB; vector[PO] betaO; vector[3] theta;
  vector<lower=0>[3] sigma_eps; matrix[ND, 3] eps_raw;
  vector[K ? J : 0] u_free; vector[K ? J : 0] v_free; vector<lower=0>[K ? 2 : 0] lam;
}
transformed parameters {
  vector[J] b = sigma_b * b_raw;                             // b_raw ~ student_t(4): heavy left tail for never-chosen states
  vector[J] del = sigma_del * del_raw;
  matrix[J, NE_ERA] eta = sigma_eta * eta_raw;
  vector[ND] bbar;
  if (fix_bbar) bbar = bbar_data;
  else { bbar = rep_vector(0, ND); for (t in 1:NT) bbar[ts_d[t]] += b[ts_j[t]] + eta[ts_j[t], era_d[ts_d[t]]]; bbar = bbar ./ ts_n; }
}
model {
  beta1 ~ normal(0, 2.5); alpha ~ normal(0, 2.5); betaA ~ normal(0, 2.5); betaB ~ normal(0, 2.5); betaO ~ normal(0, 2.5); theta ~ normal(0, 2.5);
  sigma_b ~ normal(0, 1); sigma_del ~ normal(0, 1); sigma_eta ~ normal(0, 1); sigma_eps ~ normal(0, 1);
  b_raw ~ student_t(4, 0, 1); del_raw ~ std_normal(); to_vector(eta_raw) ~ std_normal(); to_vector(eps_raw) ~ std_normal();
  vector[K ? J : 0] u; vector[K ? J : 0] v;
  if (K) { u_free ~ std_normal(); v_free ~ std_normal(); lam ~ normal(0, 1); u = u_free - mean(u_free); v = Pv * v_free; }
  if (use_s1) {
    vector[N1] eta1 = X1 * beta1;
    for (n in 1:N1) eta1[n] += b[tgt1[n]] + eta[tgt1[n], era1[n]] + coal1[n] * del[tgt1[n]];
    if (K) for (n in 1:N1) eta1[n] += lam[1] * u[ini1[n]] * v[tgt1[n]];
    target += reduce_sum(s1_partial, ev_idx, grainsize, eta1, e_start, e_end, e_chosen, w_e);
  }
  if (use_s2) {
    matrix[ND, 3] eps = eps_raw * diag_matrix(sigma_eps);
    vector[N2] UA = alpha[1] + XA * betaA; vector[N2] UB = alpha[2] + XB * betaB; vector[N2] UO = alpha[3] + XO * betaO;
    vector[ND] vA; vector[ND] vB;
    if (K) { vA = W_oppA * v; vB = W_oppB * v; }
    matrix[N2, 3] U2;
    U2[, 1] = UA + theta[1] * bbar[d2] + eps[d2, 1]; U2[, 2] = UB + theta[2] * bbar[d2] + eps[d2, 2]; U2[, 3] = UO + theta[3] * bbar[d2] + eps[d2, 3];
    if (K) { vector[N2] uk = lam[2] * u[k2]; U2[, 1] += uk .* vA[d2]; U2[, 3] += uk .* vA[d2]; U2[, 2] += uk .* vB[d2]; }
    target += reduce_sum(s2_partial, r_idx, grainsize2, U2, y2, av2);
  }
}
generated quantities {
  vector[K ? J : 0] u_out; vector[K ? J : 0] v_out;
  if (K) { u_out = u_free - mean(u_free); v_out = Pv * v_free; }
}
