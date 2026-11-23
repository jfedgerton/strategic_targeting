#!/usr/bin/env bash
## tests/test_s2_runner.sh -- resume and provenance behaviour of R/s2_05_sweep.R and R/s2_07_boot_cells.R (GPT review P1-2, P2-4, P2-8).
## Usage: tests/test_s2_runner.sh BUILD_DIR BUILDS_ROOT DUMPS_DIR EDGES_DIR POSITIONS_DIR (needs w1_fmid under BUILDS_ROOT and
## positions_core_W2contig_ptaall.rds with its edge file). About 3 min.
set -uo pipefail
BD=$1; BR=$2; DD=$3; ED=$4; PD=$5; T=$(mktemp -d); ok=1
chk() { if eval "$2"; then printf "%-90s PASS\n" "$1"; else printf "%-90s FAIL\n" "$1"; ok=0; fi; }
cp $PD/positions_core_W2contig_ptaall.rds $T/positions_core_W2contig_ptaall.rds
A="--build=$BD --dump=$DD/long_init_bon_w1d_relyear_fmid.rds --edges_dir=$ED --out=$T/sw --positions=$T/positions_core_W2contig_ptaall.rds --max_cells=2 --chunk=2 --cores=2"
Rscript R/s2_05_sweep.R $A > $T/a.log 2>&1; chk "1 sweep runs and reports DONE with a census" "grep -q S2_SWEEP_DONE $T/a.log && grep -q census $T/a.log"
Rscript R/s2_05_sweep.R $A > $T/b.log 2>&1; chk "1 rerun with identical inputs skips finished chunks" "grep -q 'S2_SWEEP_DONE (0 measures' $T/b.log"
Rscript -e 'p <- readRDS(commandArgs(TRUE)[1]); p$R[, 1] <- p$R[, 1] + 1e-3; saveRDS(p, commandArgs(TRUE)[1])' $T/positions_core_W2contig_ptaall.rds
Rscript R/s2_05_sweep.R $A > $T/c.log 2>&1; rc=$?
chk "2 changed positions CONTENTS with the same spec ids: the sweep refuses to reuse the folder" "[ $rc -ne 0 ] && grep -q 'different code or inputs' $T/c.log"
cp $PD/positions_core_W2contig_ptaall.rds $T/positions_core_W2contig_ptaall.rds
printf 'file,spec_id,window,sample,rule,n_boot,cell_id\npositions_core_W2contig_ptaall.rds,W2contig_ptaall|igo|pr,1,fmid,test,4,b001\n' > $T/cells.csv
cp $PD/positions_core_W2contig_ptaall.rds $T/pd_copy.rds; mkdir -p $T/pd; mv $T/pd_copy.rds $T/pd/positions_core_W2contig_ptaall.rds
B="--cells=$T/cells.csv --posdir=$T/pd --builds_root=$BR --dumps_dir=$DD --edges_dir=$ED --out=$T/bt --cores=2 --block=2"
Rscript R/s2_07_boot_cells.R $B > $T/d.log 2>&1; chk "3 bootstrap runs; summary has per-quantity finite counts" "grep -q S2_BOOT_DONE $T/d.log && head -1 $T/bt/b001/summary.csv | grep -q n_finite"
## a worker-lost draw (try-error) is retried on rerun; numerical failures are kept
Rscript -e 'f <- commandArgs(TRUE)[1]; x <- readRDS(f); x$draws[[1]] <- structure("Error : killed", class = "try-error"); saveRDS(x, f)' $T/bt/b001/draws_001.rds; rm $T/bt/b001/summary.csv
Rscript R/s2_07_boot_cells.R $B > $T/e.log 2>&1
chk "4 a draw lost to a worker is refitted on rerun" "Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]); quit(status = !is.numeric(x\$draws[[1]]))' $T/bt/b001/draws_001.rds"
Rscript R/s2_07_boot_cells.R ${B/--block=2/--block=3} > $T/f.log 2>&1; rc=$?
chk "5 a changed --block on an existing cell folder stops" "[ $rc -ne 0 ] && grep -q 'holds a different cell' $T/f.log"
rm -rf $T; [ $ok = 1 ] && echo "ALL PASS" || { echo "SOME TESTS FAILED"; exit 1; }
