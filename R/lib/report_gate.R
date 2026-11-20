## report_gate -- the resampling reporting rule (delta review N4; GPT review D6): a quantity gets an interval / positivity share / diagnostic p
## only if B >= 200, >= 90% of draws are usable for it, AND its observed full-sample estimate is finite (a valid observed fit).
report_gate <- function(B, okq, est) B >= 200L & okq / B >= 0.90 & is.finite(est)
