#!/usr/bin/env bash
## tests/test_hash_cache.sh -- GPT review E1 regression (portable to bash 3.2). A file is hashed in one invocation, replaced by DIFFERENT bytes of
## the SAME size with its mtime preserved, and hashed again in a new invocation: the hash must change. Within one invocation, hinv after a
## same-second rewrite must also force a re-hash. Uses the h()/hinv() definitions from run_all.sh verbatim.
set -euo pipefail
cd "$(dirname "$0")/.."
DEFS=$(sed -n '/^HC=\$(mktemp/,/^hinv()/p' run_all.sh)
T=$(mktemp -d); printf 'a,b\n1,2\n' > $T/f.csv; touch -t 202001010000 $T/f.csv
inv() { bash -c "set -euo pipefail; $DEFS; $1" ; }
H1=$(inv "h $T/f.csv")
cp -p $T/f.csv $T/ref; printf 'a,b\n1,3\n' > $T/f.csv; touch -r $T/ref $T/f.csv          # same size, same mtime, different bytes
[ "$(wc -c < $T/f.csv)" = "$(wc -c < $T/ref)" ] && [ "$(date -r $T/f.csv +%s)" = "$(date -r $T/ref +%s)" ] || { echo "fixture failed"; exit 1; }
H2=$(inv "h $T/f.csv")
[ "$H1" != "$H2" ] || { echo "FAIL: new invocation reused the old hash ($H1)"; exit 1; }
## within one invocation: hash, same-size same-mtime rewrite, hinv, re-hash -> must differ
H3=$(inv "a=\$(h $T/f.csv); printf 'a,b\n1,4\n' > $T/f.csv; touch -r $T/ref $T/f.csv; hinv $T/f.csv; b=\$(h $T/f.csv); [ \"\$a\" != \"\$b\" ] && echo ok")
[ "$H3" = "ok" ] || { echo "FAIL: hinv did not force a re-hash"; exit 1; }
rm -rf $T; echo "test_hash_cache: all checks passed"
