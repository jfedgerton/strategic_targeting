# Trade-agreement sources (copied 2026-10-11 from ottoman_byzantine/data/raw/replication_sweep)
| File | Source | Unit | Coverage | Use |
|---|---|---|---|---|
| desta_list_of_treaties_02_02.txt | DESTA list of treaties v02_02 (Dataverse doi:10.7910/DVN/Q9SX4V, via Isgandarov 2026 replication) | treaty; members as ISO numeric (c1..cN) | entry into force 1948+ | PTA layer (any agreement in force) |
| PTAMIDTables.dta | Mansfield & Pevehouse 2000 IO, PTAMID (people.reed.edu/~ahm/Projects/PTAMID) | dyad-year; dyadid = ccode1*1000 + ccode2; PTASame = shared PTAs | 1947-2000 | Alternative PTA layer |
| prd_data.csv | Isgandarov 2026 II replication | politically relevant dyad-years; joint_pta | 1960-2010, PRD only | Cross-check only (not a full panel) |
Limitation: DESTA list has no withdrawal/termination dates, so an agreement counts as in force from entry into force through 2014.
