#!/usr/bin/env bash
## ============================================================================
## config.sh -- every setting of the multiverse loop in one place. run_all.sh sources this file.
## ============================================================================

## ---- inputs (edit these paths) ----
## Built event data from the main repo (01_build_data.R --entry=later). Linked into this repo as ./data_entry
DATA_ENTRY="${DATA_ENTRY:-"$HOME/strategic_targeting_relsweep/v2/data_entry"}"
## Raw network files (dipl_dyad.csv, destaRaw.rda, dca_dyad.csv)
NETDATA="${NETDATA:-"$HOME/Dropbox/state_targeting/relsweep_data"}"
## Reuse already-built positions? The core and IGO arms do not contain the diplomacy layer, so the 2026-11-18 files are unchanged by the
## no-future-diplomacy fix; the diplomacy arm must come from the 2026-11-19 no-future build. Set REUSE_POSITIONS=0 to rebuild everything (~4 h).
REUSE_POSITIONS="${REUSE_POSITIONS:-0}"   # 2026-11-22 (GPT review SRY4): full rebuild required; run_all.sh refuses 1
REUSE_CORE_IGO="${REUSE_CORE_IGO:-"$HOME/strategic_targeting_relsweep/v2/output/rel_sweep/run_20261118_124722/positions"}"
REUSE_DIPL="${REUSE_DIPL:-"$HOME/strategic_targeting_relsweep/v2/output/rel_sweep/nofut_20261119_104425/positions"}"

## ---- compute ----
CORES="${CORES:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"

## ---- main specification (PI 2026-11-19) ----
TYPES="${TYPES:-bcd}"                 # three types: bilateral / initiating-side coalition (originators within the window + attack-side joiners) / target-side joiners
WINDOW_MAIN=1             # originator window in months (definite bound); windows 0-6 run for the appendix
WINDOWS="0 1 2 3 4 5 6"
SAMPLE_MAIN=fmid          # fatal MIDs; "mid" = all MID onsets (failed-bargain step A)
INTENSITY="jf fd2 fd3 fd4 fd5 fd6"   # intensity samples run through the full multiverse (jointly fatal; dispute fatality level >= 2..6)
MAIN_SPEC='BOOKc12_ptaany|igo+def+off+pta+dca|prpos'   # the main relational measure (PI 2026-11-23): the book / JOP 2024 measure built exactly
# (R/book_measure_c12.R: contig12 + sender/receiver/year FE, logit / PPML, observed - expected, mean of 5, symmetric), positive-part PageRank.
# The signed version (|prsigned) goes in the appendix. Built by run_stage2.sh step 1 (not part of the stage-1 grid run).
# Was 'W2contig_ptaall|igo+def+off+neu+nonagg+consul+dca+un|pr' (2026-11-18) for the stage-1 run in output_init; that run holds every measure, but its main-measure bootstrap is the old one

## ---- resampling ----
B_BOOT=2000               # pigeonhole bootstrap draws
B_PERM=1000               # permutation draws
N_RANDOM_BOOT=25          # random plausible measures bootstrapped at the main window (seed 123)

## ---- robustness grid (run on the 201 Phase A measures at the main window) ----
MERROR_K="0.25 0.5 1"     # measurement-error strengths (x within-year SD of the position); def/contig run with both signs
N_NOISE_DRAWS=5           # independent random-noise draws per strength (seeds 123..127)
