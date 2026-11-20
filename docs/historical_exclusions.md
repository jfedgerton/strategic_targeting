# History-chosen dispute exclusions — coding rule, ledger and lists (PI approved 2026-11-20)

These lists are chosen on historical grounds, separately from the influence-chosen top-5 / top-10 drops (which are adversarial). They were written **after** the influential disputes were known, so they are retrospective robustness sets, not a preregistered test. An event is dropped if **any** of its source disputes is listed, so a listed MID is removed whole, including any non-qualifying episodes inside it (overcoverage noted per case).

## Coding rule (fixed 2026-11-20, before any run with these lists)
1. **Closed list of triggers.** The coalition's main public justification at the time is one of: (a) aggression by the target against another state; (b) mass atrocities; (c) attacks on foreign nationals or diplomats; (d) defiance of UN Security Council resolutions; (e) an attack by the target state or by forces it harbours. Restoring an elected government is **not** a trigger (expanded list only).
2. **Main justification.** A qualifying rationale counts only if it was the main public justification at the time. Territorial, colonial, succession, debt-collection and alliance-war aims as the main aim exclude the case.
3. **Three states using force.** At least three states on the coalition side have a use-of-force action (COW-MID highest action ≥ 14) in the dispute, by MID participant dates. Threats, displays and support do not count.
4. **Both sides screened.** Every fatal-sample dispute with three or more distinct initiators on the initiating side, or three or more target-side joiners, was screened (fmid dump, 1-month window).
5. **UN-authorized subset.** The Security Council authorized or recommended the use of force for that operation (Korea: Res. 83/84 recommend assistance; Gulf: Res. 678 authorizes force). Recognition of self-defense (Res. 1368) or demands without a force mandate do not count.

## Ledger (screened coalition disputes)
Force = states on the coalition side with a use-of-force action (GML participant table). Role = the coalition's type in the bcd frame.

| Dispute | Case | Role | Force | Trigger / main justification | Mandate | List | Source |
|---|---|---|---|---|---|---|---|
| 51 | Korea 1950–53 | defender | 16 | (a) North Korean invasion of the South | UNSC 83/84 | UN, strict, expanded | [Res. 83](https://mandates.un.org/mandate/S/RES/83%281950%29), [Res. 84](https://www.refworld.org/legal/resolution/unsc/1950/en/80514) |
| 3957 | Gulf War 1990–91 | defender | 14 | (a) Iraqi invasion of Kuwait | UNSC 678 | UN, strict, expanded | [Res. 678](https://digitallibrary.un.org/record/102245?ln=en) |
| 31 | Boxer intervention 1900 | chooser | 8 | (c) attacks on legations and nationals; indemnity and garrisons followed (mixed aims; relief was the stated aim) | — | strict, expanded | [FRUS 1900](https://history.state.gov/historicaldocuments/frus1900/d183), [US Army](https://history.army.mil/Research/Reference-Topics/Army-Campaigns/Brief-Summaries/China-Relief-Expedition/) |
| 145 | Nanking incident 1927 (UK, US, France, Japan) | chooser | 4 | (c) attacks on foreign consulates and nationals; shelling to cover evacuation | — | strict, expanded | [FRUS 1927](https://history.state.gov/historicaldocuments/frus1927v02/d201), [Hansard 28 Mar 1927](https://api.parliament.uk/historic-hansard/commons/1927/mar/28/sir-a-chamberlains-official-account) |
| 146 | Navarino 1827 | chooser | 3 | (b) Ottoman–Egyptian conduct in Greece, enforcing the 1827 treaty; **borderline** (also intervention in an independence war; battle not ordered as such) | — | strict, expanded (dropped in strict_noborder) | [Hansard 1833](https://hansard.parliament.uk/commons/1833-03-26/debates/e972a6e2-8efd-41f0-9aa1-ca7c028085f9/CommonsChamber), [Hansard 1877](https://hansard.parliament.uk/Commons/1877-05-14/debates/cf79aab8-8890-4eba-9518-f57ae9c6268c/TheEasternQuestion%E2%80%94TheResolutions%28MrGladstone%29) |
| 3552 | Iraq no-fly zones 1992–93 | chooser | 4 (US, UK, France, Saudi Arabia) | (d)/(b) enforcing UNSC 688 against repression; no force mandate | claimed, not authorized | strict, expanded | [UK MOD briefing](https://webarchive.nationalarchives.gov.uk/ukgwa/20171123122521mp_/http://www.iraqinquiry.org.uk/media/223045/2009-11-06-briefing-mod-iraq-no-fly-zones.pdf), [AFHRA](https://www.afhra.af.mil/Portals/16/documents/Airmen-at-War/Jamieson_SouthernIraq30Sep15.pdf) |
| 4137 | Kosovo 1998–99 | chooser | 14 | (b) humanitarian catastrophe in Kosovo | none | strict, expanded | [NATO 23 Mar 1999](https://nato.int/docu/pr/1999/p99-040e.htm) |
| 4283 | Afghanistan 2001 | chooser | 9 | (e) harbouring the 9/11 attackers; self-defense | Res. 1368 recognizes self-defense (not authorization) | strict, expanded | [Res. 1368](https://documents.un.org/api/symbol/access?l=en&s=s%2Fres%2F1368%282001%29&t=pdf) |
| 3058 | Grenada 1983 | chooser | 7 (US + Caribbean states) | (c) protecting US nationals, OECS request; **borderline** (contested legality; GA Res. 38/7 condemned it) | none | strict, expanded (dropped in strict_noborder) | [Reagan 25 Oct 1983](https://www.reaganlibrary.gov/research/speeches/102583e), [GA 38/7](https://www.refworld.org/legal/resolution/unga/1983/en/6579) |
| 4273 | Iraq 1997–2003 | chooser | 5 (US, UK, Saudi Arabia 1997–99; Australia, Poland 2003) | (d) non-compliance with UNSC resolutions; Desert Fox (1998) was US–UK; the 2003 invasion is a separate episode inside the same MID | none | expanded | [UK Hansard 2003](https://publications.parliament.uk/pa/cm200203/cmhansrd/vo030324/text/30324w05.htm) |
| 4251 | Sierra Leone 1997 (ECOMOG) | chooser | 3 (Nigeria, Guinea, Ghana) | restoring the elected government after the coup — not a listed trigger | Res. 1132 (embargo enforcement only) | expanded | [Res. 1132](https://digitallibrary.un.org/record/244598?ln=en) |
| 2150 | Shaba I 1977 (Morocco, Egypt, Uganda, CAR with Zaire; France is not a MID participant) | defender | 5 | assistance to Zaire against an invasion from Angola, at Zaire's request; alliance assistance rather than collective punishment | — | expanded | source to add at coder review |
| 135 | Mexico 1861 | chooser | 3 | debt collection (main aim); later French imperial project | — | none | [US Office of the Historian](https://history.state.gov/milestones/1861-1865/french-intervention) |
| 89 | Spain 1823 | chooser | 1 | monarchical restoration (Holy Alliance) | — | none | — |
| 2208, 1129 | China 1926; Hodeida 1934 | chooser; defender | 0; 0 | protection of nationals, but no coalition state used force | — | none | — |
| 200 | Suez 1956 | chooser/defender | 3 | canal and strategic aims | — | none | — |
| 194, 261, 375, 88, 1482, 1250 | 19th-century unification, succession and regional wars | both | — | territorial / succession / alliance | — | none | — |
| 197, 257, 258, 324, 1265, 57 | WWI, WWII, Russia 1919, Crimean War | both | — | great-power / alliance wars | — | none | — |
| 611, 1035, 1046, 1039, 607, 4022, 4339, 1286, 1703, 210 | Vietnam, Arab–Israeli wars, Syria 1957, Korea 1994, Second Congo War, Cominform–Yugoslavia, Spanish Civil War | both | — | alliance, regional or civil-war intervention | — | none | — |
| Libya 2011 | — | — | — | — | UNSC 1973 | not in the fatal sample (no dispute in the build) | [Res. 1973](https://mandates.un.org/mandate/S/RES/1973%282011%29) |

## Lists (inputs/)
| File | Disputes | Purpose |
|---|---|---|
| `exclusions_un.txt` | 51, 3957 | UN-authorized coalitions (both defender-side) |
| `exclusions_strict.txt` | 51, 3957, 31, 145, 146, 3552, 4137, 4283, 3058 | Main history-chosen check |
| `exclusions_strict_noborder.txt` | strict minus 146, 3058 | Borderline cases kept in |
| `exclusions_expanded.txt` | strict + 4273, 4251, 2150 | Open-trigger / mixed-episode check |

**Before treating these as final:** a second coder should review this ledger without seeing any estimates (GPT delta review §3). Changing a list later reruns only its variant and sweep (the variant is stamped with the list).
