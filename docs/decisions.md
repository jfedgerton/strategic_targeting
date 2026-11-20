# Design decisions (settled; with dates)

## Specification
- **Three initiator types** (2026-11-19): bilateral / initiating-side coalition / target-side joiner. Pooling coalition originators with attack-side joiners is a substantive choice (both act on the attacking side of a coalition dispute) and an assumption. Failing to reject equal surfaces (p 0.55–0.98 on the main measure at windows of 1 month or more) is consistent with it but is not evidence of equivalence, and with 141 originator events the test has limited power. The four-type model, which drops the assumption, is reported in the appendix at every window.
- **Originator window = 1 month** (2026-11-19), definite date bound. Rationale: entry dates are recorded at the month, often with bounds; one month is the smallest window that treats "entered the same or next month" as entering at onset, measured with date error. Windows 0–6 are in the appendix. Under three types the window only moves 25–31 bilateral initiators whose partner arrived within it.
- **Window promotions are side-specific** (2026-11-19): only joiners on the initiating side can become coalition originators; an early entrant on the target side stays a target-side joiner.
- **Coalition status at onset** (`--bilateral=onset`): only partners present at initiation make an initiator a coalition initiator.
- **Fatal MIDs, all years, WWII included, full choice set** (standing decisions).

## Joiner side (2026-11-19, after the GPT full audit)
COW-MID Side A is the side of the first codable action, not the aggressor. A joiner is on the *initiating* side when its initial COW side equals the dispute's initiating side:
1. the COW side of events between the dispute's own originators (participant-level originator flags, not the pair-level orig1/orig2, which mark some joiners as originators);
2. otherwise the project's direction rule: the side of the COW war initiator for war-linked disputes, Side A for the rest.

Case review of the 13 disputes rule 1 could not settle (encoded as regression checks in `R/05_dump_events.R`):

| Dispute | Initiating side | Initiating-side joiners | Target-side joiners |
|---|---|---|---|
| 112 Italy–Papal States 1860 | Italy | — | France |
| 362 Korea 1976 | North Korea | — | US |
| 1293 Cyprus 1974 | Turkey | — | Greece |
| 1363 Laos 1960–61 | USSR | North Vietnam | US, Thailand |
| 3098 Libya–Sudan 1981 | Libya | — | Somalia, US, Egypt, Oman |
| 1035 Six-Day War 1967 | Israel | — | Iraq, Egypt, Jordan, Saudi Arabia, Kuwait |
| 200 Suez 1956 (**hard**) | Israel (COW war 155) | UK, France, US | USSR |
| 611 Vietnam 1964 | US (COW war 163) | South Vietnam, Australia, South Korea, NZ, Philippines, Thailand, Cambodia | USSR, China |
| 1046 Yom Kippur 1971–73 | Egypt/Syria (COW war 181) | Egypt, Jordan, Iraq, Saudi Arabia | — |
| 3826 Romania–Hungary 1944 | Germany (WWII, COW war 139) | Germany | — |
| 2150 Shaba I 1977 | Angola/Cuba (Side A) | Congo-Brazzaville, Uganda | France, Morocco, Egypt, CAR |
| 3131 Rhodesia–Zambia 1977–79 | Rhodesia (Side A) | — | Tanzania |
| 1379 Korean DMZ 1966–69 (**hard**) | US/South Korea (Side A) | USSR (1969) | — |

Hard cases: 200 (COW codes the US on Israel's side although the US opposed the invasion) and 1379 (USSR on the US/SK side in 1969; North Korea arguably initiated the DMZ conflict). Both are dropped in the `dropside_hard` check; all 13 in `dropside_reviewed`.

## Measurement
- **No future information in the diplomacy layer** (2026-11-19): each dyad-year uses the latest diplomatic-exchange wave at or before its year; zero before the first wave; a state that entered after its source wave is zero until observed.
- **Plausibility rules** (2026-11-18, set before any coefficients): US and USSR above the median 1950–89; China rises from 1950–75 to 1990–2014; mean of UK/France/Germany above the median 1816–1913; North Korea below the 90th percentile after 1991. Measures that fail are still run and reported.
- **Main measure** (2026-11-18): contiguity-residualized PageRank on IGO + the five ATOP types + DCA + UN agreement, PTA coded "all", no diplomacy.

## Other fixes from the GPT full audit (2026-11-19)
- Aggression window anchored on the GML dispute start; defense partners counted in the year before the dispute start (entry-year version kept as a separate check).
- M0 pools every control (one model under every scheme); each model has its own fit validity.
- Lag comparisons can use one pooled scale so only information differs; the decomposition uses an earlier non-overlapping baseline and can add the omitted quadratic terms (tested jointly).
- Drops remove an event if any of its source disputes is listed.
