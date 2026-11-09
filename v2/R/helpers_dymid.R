## ============================================================================
## helpers_dymid.R -- the "interacting dyads" event foundation (switch events = dymid, PI-approved 2026-11-04; README section 6.9)
##
## Functions only (sourced by 01_build_data.R under --events=dymid and by tests/testthat/test_dymid.R). Sources, both in v2/data_raw/dymid
## (provenance in v2/data_raw/SOURCE.md, sha256 checked by the build against the constants below):
##   dyadic_mid_4.03.csv     Maoz et al., Dyadic MID 4.03 (1816-2014): one row per DIRECTED dyad-year of a dyadic dispute, both directions
##                           listed (A, B) and (B, A); only directly interacting pairs. In 4.03 the column dyindex is an integer equal to disno
##                           on every row (the .001 dyad suffix of 4.01 is lost), so it cannot identify a dyad: the dyadic dispute is keyed by
##                           (disno, unordered pair) here, and a pair re-entering the same dispute (separate dyindices in 4.01; 25 pairs, 18 in
##                           WWI / WWII) is ONE event dated by its first spell; every spell is its own activity period in the pair-dispute
##                           spells (2026-11-08, dymid_override_pdi()). The two directed copies of a dyad-dispute do not always agree
##                           (2026-11-08): dymid_dyads() validates them field by field and resolves every disagreement through the checked-in
##                           adjudication table data_raw/dymid/adjudication.csv (rules in DYMID_ADJ_RULES); an uncovered one stops the build.
##   directed_dyadic_war.csv Maoz, Directed Dyadic Interstate War data (1816-2010): one row per directed war dyad-year (warnum = COW war number,
##                           disno, dyindex with the decimal suffix). A dyadic war = direct battles with >= 1,000 battle deaths across both states.
##                           batdtha / batdthb are sometimes national totals: never used as dyadic magnitudes.
##   adjudication.csv        the copy disagreements of dyadic_mid_4.03.csv the build resolves (one row per dyad, field and sample it applies to:
##                           disno, a, b, field, value_copy1, value_copy2, resolved_value, applies_to, rule, note). copy1 = the row listing the
##                           smaller ccode as state A, copy2 = the other row mirrored into the same (a, b) frame. Checked against the data on
##                           every build: a row whose copy values differ from the data, a disagreement without a row, a row without a
##                           disagreement, or a resolved value that does not follow from its rule stops the build.
## ============================================================================
DYMID_CSV_SHA256 <- "9d26eedd792a0067f590f038c3b7b240bd373cd70200ee729241e1078c6957c3"   # data_raw/dymid/dyadic_mid_4.03.csv
DDIW_CSV_SHA256  <- "aa2eb5244a0cb0b579867e50eb02e4297fb29250e67120fb132e3c0618dceda0"   # data_raw/dymid/directed_dyadic_war.csv
DYMID_DIR_LABELS <- c(a_joiner_orig = "(a) exactly one joiner: the joiner initiates", a_joiner_later = "(a) two joiners: the later entrant initiates",
                      b_dyadic = "(b) both originators (or a joiner tie): the dataset's dyadic initiator (dyad_rolea / wardyadrolea = 1)",
                      c_cowwar = "(c) war dyad that (a) and (b) do not decide (no dyadic initiator): COW war-initiating side (D1)",
                      fallback_sideA = "no rule decides (no dyadic initiator: simultaneous first action or the two copies contradict each other; no D1 decision): COW Side A kept")
## the codebook values of DYMID recip (1 = reciprocated, initiator not identifiable, roles assigned randomly; 0 = not). 4.03 also carries recip = 2 on
## 4 rows (dispute 4721, 2012: 731 - 732 and 731 - 740, both directions), a value the codebook does not define: accepted explicitly and treated as
## NOT reciprocated (recip flag = recip == 1); any other value stops read_dymid()
DYMID_RECIP_CODEBOOK <- 0:1; DYMID_RECIP_UNDOCUMENTED <- 2L
## ---- copy disagreements (2026-11-08) -------------------------------------------------------------------------------------------------------
## The dyad-level fields of dymid_dyads() that the build uses and that are validated between the two directed copies of every dyad-dispute, with
## the rule the adjudication table may apply to each (DYMID_ADJ_RULES; a disagreement on a field with no rule, or on a field not listed here, has
## no resolution and stops the build). Values as written in adjudication.csv: start_ym / end_ym as "YYYY-MM", fatality and hostility codes as
## integers with -9 for missing, init_dyadic as the initiator's ccode or NA, war_dymid as 0 / 1, side_a_dy / side_b_dy as the distinct per-year
## codes of the copy joined by ";" (NA = not used).
DYMID_COPY_FIELDS <- c("start_ym", "end_ym", "fatlev_max", "fatleva_max", "fatlevb_max", "war_dymid", "init_dyadic", "orig_a_dy", "orig_b_dy", "recip",
                       "hihost", "hihosta", "hihostb", "side_a_dy", "side_b_dy")
DYMID_ADJ_RULES <- c(
  direction_uninformative = "the copies name different dyadic initiators: rule (b) is uninformative for this dyad (init_dyadic = NA), so the direction falls to the next rule, (c) D1, else the Side A fallback (PI-approved 2026-11-08)",
  war_defer_to_war_data   = "the copies disagree on the war flag: the dyad is a war dyad iff the Directed Dyadic War data has a war dyad for (disno, pair) (PI-approved 2026-11-08)",
  fatal_strict_lower      = "the copies disagree on a fatality code: the strict sample fmid takes the LOWER code (PI-approved 2026-11-08); the order is 0 < missing (-9) < 1 < ... < 6, so a missing code never makes a dyad strict-fatal",
  fatal_ub_higher         = "the copies disagree on a fatality code: fmidub (upper bound) and fmidany ('any') take the HIGHER code (PI-approved 2026-11-08); the order is 0 < missing (-9) < 1 < ... < 6, so a missing code is the upper bound of a 0 and a positive code the upper bound of a missing one",
  earliest_start          = "the copies disagree on the dyad's start month: the earliest start, as the dyad's start is the minimum over its rows within a copy (the convention of the build before 2026-11-08, when the A < B copy was used and the mismatch counted)",
  latest_end              = "the copies disagree on the dyad's end month: the latest end, as the dyad's end is the maximum over its rows within a copy",
  hostility_higher        = "the copies disagree on a hostility level: the higher level, as the dyad's hostility is the maximum over its rows within a copy (hostlev1 / hostlev2 / hostlev_ev are not sample inputs under dymid)",
  sides_uninformative     = "the copies disagree on a state's side: DYMID's sides are not used for this dyad (NA in every year); the participation table takes COW's sides, else the side opposite the partner (DYMID sides enter only for states COW lacks, none in this data)")
## the fatality order of the two fatal rules: 0 < missing < 1 < ... < 6 (missing codes sit between 0 and the positive codes, see DYMID_ADJ_RULES)
fat_rank <- function(v) ifelse(is.na(v), 0.5, as.numeric(v))

## na9(v): the DYMID / COW missing codes (negative) as NA
na9 <- function(v) ifelse(v < 0L, NA_integer_, as.integer(v))
## max_na(v): the maximum of the non-missing values, NA when all are missing
max_na <- function(v) { v <- v[!is.na(v)]; if (length(v)) max(v) else NA_integer_ }

## read_dymid(path): the columns of Dyadic MID 4.03 the build uses, as integers; asserts the structure the code relies on (both directions of every
## dyad-year present; dyindex == disno, the 4.03 truncation documented above)
read_dymid <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("disno", "dyindex", "statea", "stateb", "year", "strtday", "strtmnth", "strtyr", "endday", "endmnth", "endyear", "fatlev", "hihost", "recip",
            "sideaa", "sideab", "fatleva", "fatlevb", "hihosta", "hihostb", "orignata", "orignatb", "rolea", "roleb", "dyad_rolea", "dyad_roleb", "war")
  miss <- setdiff(need, names(x)); if (length(miss)) stop("read_dymid: columns absent: ", paste(miss, collapse = ", "))
  x <- as.data.frame(lapply(x[, need], as.integer))
  stopifnot(!anyNA(x[, c("disno", "statea", "stateb", "year", "strtyr", "strtmnth", "endyear", "endmnth", "orignata", "orignatb", "dyad_rolea", "war")]),
            all(x$statea != x$stateb), all(x$strtmnth %in% 1:12), all(x$endmnth %in% 1:12), all(x$war %in% 0:1), all(x$orignata %in% 0:1),
            all(x$dyad_rolea %in% 1:3), all(x$strtyr <= x$year))
  bad_r <- is.na(x$recip) | !x$recip %in% c(DYMID_RECIP_CODEBOOK, DYMID_RECIP_UNDOCUMENTED)
  if (any(bad_r)) stop("read_dymid: recip takes values outside {0, 1} and the documented undocumented code 2 on ", sum(bad_r), " rows: ", paste(unique(x$recip[bad_r]), collapse = ", "))
  attr_r2 <- x[x$recip == DYMID_RECIP_UNDOCUMENTED, c("disno", "statea", "stateb", "year", "recip")]
  if (!all(x$dyindex == x$disno)) stop("read_dymid: dyindex differs from disno on some rows; the (disno, pair) key of helpers_dymid.R assumes the 4.03 truncation, review")
  ## both directions of every dyad-year: the multiset of (disno, a, b, year) keys of the A < B rows equals that of the A > B rows. The start month of
  ## the two copies differs on a few rows (3250 Russia - China 1900: July / August); the A < B copy is used, the count is in attr "n_start_mismatch"
  k <- function(d, s = FALSE) sort(paste(d$disno, pmin(d$statea, d$stateb), pmax(d$statea, d$stateb), d$year, if (s) paste(d$strtyr, d$strtmnth) else ""))
  if (!identical(k(x[x$statea < x$stateb, ]), k(x[x$statea > x$stateb, ]))) stop("read_dymid: the two directions of the dyad-years do not mirror each other")
  attr(x, "n_start_mismatch") <- length(setdiff(k(x[x$statea < x$stateb, ], TRUE), k(x[x$statea > x$stateb, ], TRUE)))
  attr(x, "recip2_rows") <- attr_r2          # the recip = 2 rows (not in the codebook; treated as recip = 0), for the build log
  x
}
## dymid_year_conflicts(x): the rows of read_dymid() whose listed year falls AFTER their own end date (year > endyear; year < strtyr already stops
## read_dymid()), one row per (disno, a < b, spell start, spell end) with the conflicting years of each direction (";"-joined, sorted). Codex second
## review P1 (2026-11-09): the aggregation used to take y_last from such years and dymid_override_pdi() clamped the explicit end date forward into
## them, labelled exact. 4.03 has one case (4182 Syria - Israel, rows 1997-1999 dated 1995-05 to 1996-10, both directions).
dymid_year_conflicts <- function(x) {
  b <- x[x$year > x$endyear, c("disno", "statea", "stateb", "year", "strtyr", "strtmnth", "endyear", "endmnth")]
  if (!nrow(b)) return(data.frame(disno = integer(0), a = integer(0), b = integer(0), spell_start = character(0), spell_end = character(0),
                                  years_ab = character(0), years_ba = character(0), stringsAsFactors = FALSE))
  b |> mutate(a = pmin(statea, stateb), b = pmax(statea, stateb), ab = statea < stateb,
              spell_start = ym_label(ym_index(strtyr, strtmnth)), spell_end = ym_label(ym_index(endyear, endmnth))) |>
    group_by(disno, a, b, spell_start, spell_end) |>
    summarise(years_ab = paste(sort(unique(year[ab])), collapse = ";"), years_ba = paste(sort(unique(year[!ab])), collapse = ";"), .groups = "drop") |>
    arrange(disno, a, b, spell_start) |> as.data.frame()
}
## dymid_resolve_year_conflicts(x, tbl): resolves every dymid_year_conflicts() case through the checked-in table data_raw/dymid/year_date_conflicts.csv
## (columns disno, a, b, spell_start, spell_end, years_outside, resolution, note; PI decision 2026-11-09). The table is checked against the data both
## ways: a conflict without a row, a row without a conflict, or listed years that differ from the data (either direction) stops the build. The only
## resolution is trust_dates: the explicit start / end dates stand and the rows whose year falls outside them are dropped (both directions), so no
## spell is extended past its stated end. Returns x without those rows (attr "year_conflicts" = the applied table, for the build log).
DYMID_YEAR_RESOLUTIONS <- c(trust_dates = "the explicit dates stand; rows listed in years outside them are dropped (both directions)")
dymid_resolve_year_conflicts <- function(x, tbl) {
  need <- c("disno", "a", "b", "spell_start", "spell_end", "years_outside", "resolution", "note")
  miss <- setdiff(need, names(tbl)); if (length(miss)) stop("dymid_resolve_year_conflicts: columns absent: ", paste(miss, collapse = ", "))
  tbl <- tbl[, need]; for (cc in c("disno", "a", "b")) tbl[[cc]] <- as.integer(tbl[[cc]])
  stopifnot(!anyNA(tbl$disno), !anyNA(tbl$a), !anyNA(tbl$b), all(tbl$a < tbl$b))
  bad <- setdiff(tbl$resolution, names(DYMID_YEAR_RESOLUTIONS)); if (length(bad)) stop("dymid_resolve_year_conflicts: unknown resolution(s): ", paste(bad, collapse = ", "))
  key <- function(z) paste(z$disno, z$a, z$b, z$spell_start, z$spell_end)
  if (anyDuplicated(key(tbl))) stop("dymid_resolve_year_conflicts: duplicated (disno, a, b, spell_start, spell_end) rows")
  cf <- dymid_year_conflicts(x)
  unc <- cf[!key(cf) %in% key(tbl), ]
  if (nrow(unc)) stop("dymid_resolve_year_conflicts: ", nrow(unc), " spell(s) of dyadic_mid_4.03.csv list years after their own end date and are not in ",
                      "data_raw/dymid/year_date_conflicts.csv (historical adjudication needed; the build does not extend a stated end date): ",
                      paste(sprintf("dispute %d %d-%d %s..%s years %s / %s", unc$disno, unc$a, unc$b, unc$spell_start, unc$spell_end, unc$years_ab, unc$years_ba), collapse = "; "))
  stale <- tbl[!key(tbl) %in% key(cf), ]
  if (nrow(stale)) stop("dymid_resolve_year_conflicts: ", nrow(stale), " table row(s) match no year / date conflict of the data (stale table): ",
                        paste(sprintf("dispute %d %d-%d %s..%s", stale$disno, stale$a, stale$b, stale$spell_start, stale$spell_end), collapse = "; "))
  m <- match(key(tbl), key(cf))
  badv <- tbl$years_outside != cf$years_ab[m] | tbl$years_outside != cf$years_ba[m]
  if (any(badv)) stop("dymid_resolve_year_conflicts: listed years differ from the data (both directions must list exactly years_outside): ",
                      paste(sprintf("dispute %d %d-%d table %s, data %s / %s", tbl$disno[badv], tbl$a[badv], tbl$b[badv], tbl$years_outside[badv], cf$years_ab[m][badv], cf$years_ba[m][badv]), collapse = "; "))
  drop <- x$year > x$endyear
  out <- x[!drop, , drop = FALSE]; rownames(out) <- NULL
  stopifnot(nrow(out) == nrow(x) - sum(vapply(strsplit(tbl$years_outside, ";"), length, 1L)) * 2L, nrow(dymid_year_conflicts(out)) == 0L)
  k <- function(d) sort(paste(d$disno, pmin(d$statea, d$stateb), pmax(d$statea, d$stateb), d$year))
  if (!identical(k(out[out$statea < out$stateb, ]), k(out[out$statea > out$stateb, ]))) stop("dymid_resolve_year_conflicts: the two directions no longer mirror each other")
  for (a_ in c("n_start_mismatch", "recip2_rows")) attr(out, a_) <- attr(x, a_)
  attr(out, "year_conflicts") <- tbl
  out
}
## read_ddiw(path): the columns of the Directed Dyadic Interstate War data the build uses
read_ddiw <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("warnum", "disno", "dyindex", "statea", "stateb", "warstrtmnth", "warstrtday", "warstrtyr", "warendmnth", "warenday", "warendyr", "year",
            "warolea", "waroleb", "wardyadrolea", "wardyadroleb", "batdtha", "batdthb", "batdths")
  miss <- setdiff(need, names(x)); if (length(miss)) stop("read_ddiw: columns absent: ", paste(miss, collapse = ", "))
  x <- x[, need]
  for (cc in setdiff(need, "dyindex")) x[[cc]] <- as.integer(x[[cc]])
  ## one data error: USA - Austria-Hungary 1918 (war 106) has start month 24 / day 10 and end month 3 / day 11, i.e. day and month swapped in both
  ## dates (24 October and 3 November 1918). Rows with a start month above 12 and a start day of at most 12 get both dates swapped (attr "n_swapped").
  sw <- x$warstrtmnth > 12L & x$warstrtday %in% 1:12
  if (any(sw)) {
    m <- x$warstrtmnth[sw]; x$warstrtmnth[sw] <- x$warstrtday[sw]; x$warstrtday[sw] <- m
    m <- x$warendmnth[sw]; x$warendmnth[sw] <- x$warenday[sw]; x$warenday[sw] <- m
  }
  attr_sw <- sum(sw)
  stopifnot(!anyNA(x[, c("warnum", "disno", "statea", "stateb", "warstrtmnth", "warstrtyr", "warendmnth", "warendyr", "year", "warolea", "wardyadrolea")]),
            all(x$statea != x$stateb), all(x$warstrtmnth %in% 1:12), all(x$warendmnth %in% 1:12), all(x$warolea %in% 1:4), all(x$wardyadrolea %in% c(1L, 3L)),
            all(x$warstrtyr <= x$year), all(x$year <= x$warendyr))
  k <- function(d) sort(paste(d$warnum, d$disno, pmin(d$statea, d$stateb), pmax(d$statea, d$stateb), d$year))
  if (!identical(k(x[x$statea < x$stateb, ]), k(x[x$statea > x$stateb, ]))) stop("read_ddiw: the two directions of the war dyad-years do not mirror each other")
  stopifnot(all(ym_index(x$warstrtyr, x$warstrtmnth) <= ym_index(x$warendyr, x$warendmnth)))
  attr(x, "n_swapped") <- attr_sw
  x
}

## dymid_copies(dy): the two directed copies of every dyad-year of read_dymid() in ONE frame (a < b): copy1 = the rows with statea < stateb as they
## are; copy2 = the rows with statea > stateb with the per-state columns swapped (statea / stateb, sideaa / sideab, fatleva / fatlevb, hihosta /
## hihostb, orignata / orignatb, rolea / roleb, dyad_rolea / dyad_roleb), so that both describe the same (a, b) and can be compared field by field.
dymid_copies <- function(dy) {
  c2 <- dy[dy$statea > dy$stateb, ]
  for (p in list(c("statea", "stateb"), c("sideaa", "sideab"), c("fatleva", "fatlevb"), c("hihosta", "hihostb"), c("orignata", "orignatb"), c("rolea", "roleb"), c("dyad_rolea", "dyad_roleb"))) {
    tmp <- c2[[p[1]]]; c2[[p[1]]] <- c2[[p[2]]]; c2[[p[2]]] <- tmp
  }
  list(copy1 = dy[dy$statea < dy$stateb, ], copy2 = c2)
}
## dymid_dyads_copy(r): the dyad-level aggregation of ONE copy (rows with statea < stateb, or a mirrored copy); the columns of dymid_dyads() below
## without the adjudication columns; attr "years" and attr "spells".
dymid_dyads_copy <- function(r) {
  r <- r |> rename(a = statea, b = stateb) |>
    mutate(rs = ym_index(strtyr, strtmnth), re = ym_index(endyear, endmnth), fatlev = na9(fatlev), fatleva = na9(fatleva), fatlevb = na9(fatlevb))
  ## row-level dates are rs / re so that no summarise() expression below sees an already-reduced column (recip_code2 before recip, for the same reason)
  out <- r |> group_by(disno, a, b) |> arrange(rs, year, .by_group = TRUE) |>
    summarise(onset_year = min(strtyr), start_ym = min(rs), end_ym = max(re), y_first = min(year), y_last = max(year),
              n_spells = n_distinct(rs), n_rows = n(),
              fatlev_max = max_na(fatlev), fatleva_max = max_na(fatleva), fatlevb_max = max_na(fatlevb),
              fatlev_last = max_na(fatlev[year == max(year)]), fatleva_last = max_na(fatleva[year == max(year)]), fatlevb_last = max_na(fatlevb[year == max(year)]),
              war_dymid = any(war == 1L),
              role_first = list(unique(dyad_rolea[rs == min(rs)])), role_varies = n_distinct(dyad_rolea) > 1L,
              orig_a_dy = max(orignata), orig_b_dy = max(orignatb), orig_varies = n_distinct(orignata) > 1L | n_distinct(orignatb) > 1L,
              recip_code2 = any(recip == DYMID_RECIP_UNDOCUMENTED), recip = any(recip == 1L), hihosta = max_na(hihosta), hihostb = max_na(hihostb), hihost = max_na(hihost), .groups = "drop")
  out <- out |> mutate(role_varies_first = lengths(role_first) > 1L,
                       role1 = vapply(role_first, function(z) if (length(z) == 1L) as.integer(z) else NA_integer_, 1L),
                       init_dyadic = case_when(role1 == 1L ~ a, role1 == 3L ~ b, TRUE ~ NA_integer_)) |> select(-role_first, -role1)
  stopifnot(nrow(out) == nrow(distinct(r, disno, a, b)), all(out$onset_year == out$y_first), all(ym_year(out$start_ym) == out$onset_year), !anyNA(out$war_dymid))
  attr(out, "years") <- r |> group_by(disno, a, b, year) |>
    summarise(start_ym = min(rs), end_ym = max(re), side_a_dy = as.integer(sideaa[1]), side_b_dy = as.integer(sideab[1]),
              n_side = n_distinct(sideaa), .groups = "drop")
  ## one row per spell (distinct start month within the dyad-dispute): its end = the latest end over its rows, its years = the rows' years
  attr(out, "spells") <- r |> group_by(disno, a, b, spell_start_ym = rs) |>
    summarise(spell_end_ym = max(re), y_first = min(year), y_last = max(year), .groups = "drop") |> arrange(disno, a, b, spell_start_ym) |>
    group_by(disno, a, b) |> mutate(spell = row_number()) |> ungroup()
  out
}
## dymid_adjudication_read(path): the adjudication table (see the header), validated: the columns, every rule known, applies_to in the allowed set,
## no duplicated (disno, a, b, field, applies_to). An absent file is an error (the build needs the checked-in table).
dymid_adjudication_read <- function(path) {
  if (!file.exists(path)) stop("dymid_adjudication_read: the adjudication table ", path, " does not exist")
  adj <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0))
  need <- c("disno", "a", "b", "field", "value_copy1", "value_copy2", "resolved_value", "applies_to", "rule", "note")
  miss <- setdiff(need, names(adj)); if (length(miss)) stop("dymid_adjudication_read: columns absent: ", paste(miss, collapse = ", "))
  adj <- adj[, need]; for (cc in c("disno", "a", "b")) adj[[cc]] <- as.integer(adj[[cc]])
  stopifnot(!anyNA(adj$disno), !anyNA(adj$a), !anyNA(adj$b), all(adj$a < adj$b))
  bad <- setdiff(adj$rule, names(DYMID_ADJ_RULES)); if (length(bad)) stop("dymid_adjudication_read: unknown rule(s): ", paste(bad, collapse = ", "))
  bad <- setdiff(adj$field, DYMID_COPY_FIELDS); if (length(bad)) stop("dymid_adjudication_read: field(s) the build does not validate: ", paste(bad, collapse = ", "))
  stopifnot(all(adj$applies_to %in% c("all", "fmid", "fmidany,fmidub")))
  if (anyDuplicated(adj[, c("disno", "a", "b", "field", "applies_to")])) stop("dymid_adjudication_read: duplicated (disno, a, b, field, applies_to) rows")
  adj
}
## dymid_adj_format(field, v): a field value as adjudication.csv writes it (see DYMID_COPY_FIELDS)
dymid_adj_format <- function(field, v) {
  if (field %in% c("start_ym", "end_ym")) return(ifelse(is.na(v), "NA", ym_label(v)))
  if (field %in% c("fatlev_max", "fatleva_max", "fatlevb_max", "hihost", "hihosta", "hihostb")) return(ifelse(is.na(v), "-9", as.character(as.integer(v))))
  if (field == "war_dymid") return(as.character(as.integer(v)))
  ifelse(is.na(v), "NA", as.character(v))
}
## dymid_copy_disagreements(d1, d2): the dyad-level disagreements between two dymid_dyads_copy() results (same dyads, asserted) on DYMID_COPY_FIELDS:
## one row per (disno, a, b, field) with value_copy1 / value_copy2 as formatted strings and the raw values (v1 / v2 lists). The sides are compared
## over the dyad's years (attr "years"): a disagreement in any year lists the distinct codes of each copy.
dymid_copy_disagreements <- function(d1, d2) {
  stopifnot(nrow(d1) == nrow(d2), all(d1$disno == d2$disno), all(d1$a == d2$a), all(d1$b == d2$b))
  y1 <- attr(d1, "years"); y2 <- attr(d2, "years")
  stopifnot(nrow(y1) == nrow(y2), all(y1$disno == y2$disno), all(y1$a == y2$a), all(y1$b == y2$b), all(y1$year == y2$year))
  rows <- list()
  for (f in setdiff(DYMID_COPY_FIELDS, c("side_a_dy", "side_b_dy"))) {
    v1 <- d1[[f]]; v2 <- d2[[f]]
    d <- !((is.na(v1) & is.na(v2)) | (!is.na(v1) & !is.na(v2) & v1 == v2))
    if (any(d)) rows[[f]] <- data.frame(disno = d1$disno[d], a = d1$a[d], b = d1$b[d], field = f, value_copy1 = dymid_adj_format(f, v1[d]), value_copy2 = dymid_adj_format(f, v2[d]), stringsAsFactors = FALSE)
  }
  for (f in c("side_a_dy", "side_b_dy")) {
    v1 <- y1[[f]]; v2 <- y2[[f]]
    d <- !((is.na(v1) & is.na(v2)) | (!is.na(v1) & !is.na(v2) & v1 == v2))
    if (any(d)) {
      k <- distinct(y1[d, c("disno", "a", "b")])
      key_y <- paste(y1$disno, y1$a, y1$b); key_k <- paste(k$disno, k$a, k$b)
      pv <- function(v) vapply(key_k, function(z) paste(unique(dymid_adj_format(f, v[key_y == z])), collapse = ";"), "")
      rows[[f]] <- data.frame(disno = k$disno, a = k$a, b = k$b, field = f, value_copy1 = unname(pv(v1)), value_copy2 = unname(pv(v2)), stringsAsFactors = FALSE)
    }
  }
  out <- bind_rows(rows)
  if (!nrow(out)) out <- data.frame(disno = integer(0), a = integer(0), b = integer(0), field = character(0), value_copy1 = character(0), value_copy2 = character(0), stringsAsFactors = FALSE)
  out |> arrange(disno, a, b, field)
}
## dymid_dyads(dy, adj, war_pairs): one row per dyadic dispute (disno, a < b) of read_dymid(), built from BOTH directed copies: the aggregation of
## dymid_dyads_copy() on each copy, validated field by field (DYMID_COPY_FIELDS), every disagreement resolved by the adjudication table adj
## (dymid_adjudication_read(); NULL = an empty table, so any disagreement stops) under its rule (DYMID_ADJ_RULES), the table itself checked against
## the data (copy values, coverage both ways, resolved value recomputed from the rule). war_pairs = distinct (disno, a, b) of the war dyads of the
## Directed Dyadic War data (ddiw_dyads()), needed by the rule war_defer_to_war_data. Columns:
##   onset_year / start_ym       the dyad's own first start (dyad-level dates: strtyr / strtmnth of its earliest spell; month index ym_index)
##   end_ym                      the latest dyad end (endyear / endmnth) over its rows
##   n_spells                    distinct start months of the dyad within the dispute (> 1: the pair re-entered; one event, first spell; every
##                               spell its own activity period in the pair-dispute spells, attr "spells")
##   fatlev / fatleva / fatlevb  the MAXIMUM non-missing value over the dyad's rows (the codebook lists them for the last year of the dyadic dispute
##                               only; -9 = missing -> NA); *_last = the value in the dyad's last year (max over that year's rows), for the log.
##                               Under a copy disagreement *_max carries the STRICT (lower) code (fmid) and *_ub the HIGHER code (fmidany, fmidub);
##                               without one *_ub == *_max
##   war_dymid                   any row with war == 1 (a copy disagreement: the war data decides)
##   init_dyadic                 the dyadic initiator of the FIRST spell (dyad_rolea == 1 -> a, == 3 -> b, == 2 simultaneous -> NA); NA as well when
##                               the first-spell rows disagree (counted: role_varies_first) or the copies disagree (copy_adjudicated, field
##                               init_dyadic); role_varies = dyad_rolea differs over all the dyad's rows
##   orig_a_dy / orig_b_dy       DYMID orignata / orignatb (max over the rows; orig_varies counts disagreement)
##   recip                       any row with recip == 1 (the codebook: initiator unidentifiable, roles assigned randomly); recip = 2 (undocumented,
##                               read_dymid()) counts as NOT reciprocated, explicitly; recip_code2 = any row with recip == 2 (for the log)
##   hihosta / hihostb / hihost  max over the rows (a copy disagreement: the higher level)
##   copy_adjudicated            TRUE when any field of the dyad was resolved by the table; adj_fields lists them (";")
## attr "years": one row per (disno, a, b, year) with the row's spell start / end month and the DYMID sides of a / b (sideaa / sideab; NA for a dyad
## whose sides were adjudicated as uninformative). attr "spells": one row per (disno, a, b, spell) with spell_start_ym, spell_end_ym, y_first,
## y_last (the resolved dates for the first / last spell under earliest_start / latest_end). attr "adjudication": the applied table (the
## disagreements with rule, resolved value and applies_to), for the build log.
dymid_dyads <- function(dy, adj = NULL, war_pairs = NULL) {
  cp <- dymid_copies(dy)
  d1 <- dymid_dyads_copy(cp$copy1); d2 <- dymid_dyads_copy(cp$copy2)
  dis <- dymid_copy_disagreements(d1, d2)
  if (is.null(adj)) adj <- data.frame(disno = integer(0), a = integer(0), b = integer(0), field = character(0), value_copy1 = character(0), value_copy2 = character(0),
                                      resolved_value = character(0), applies_to = character(0), rule = character(0), note = character(0), stringsAsFactors = FALSE)
  adj <- adj[order(adj$disno, adj$a, adj$b, adj$field, adj$applies_to), , drop = FALSE]; rownames(adj) <- NULL
  key <- function(z) paste(z$disno, z$a, z$b, z$field)
  ## coverage both ways and the copy values
  unc <- dis[!key(dis) %in% key(adj), ]
  if (nrow(unc)) stop("dymid_dyads: ", nrow(unc), " copy disagreement(s) of dyadic_mid_4.03.csv are not covered by the adjudication table (add a row per field, with a rule from DYMID_ADJ_RULES): ",
                      paste(sprintf("dispute %d %d-%d %s (copy1 %s, copy2 %s)", unc$disno, unc$a, unc$b, unc$field, unc$value_copy1, unc$value_copy2), collapse = "; "))
  stale <- adj[!key(adj) %in% key(dis), ]
  if (nrow(stale)) stop("dymid_dyads: ", nrow(stale), " adjudication row(s) match no copy disagreement of the data (stale table): ",
                        paste(sprintf("dispute %d %d-%d %s", stale$disno, stale$a, stale$b, stale$field), collapse = "; "))
  m <- match(key(adj), key(dis))
  badv <- adj$value_copy1 != dis$value_copy1[m] | adj$value_copy2 != dis$value_copy2[m]
  if (any(badv)) stop("dymid_dyads: adjudication row(s) whose copy values differ from the data: ",
                      paste(sprintf("dispute %d %d-%d %s (table %s / %s, data %s / %s)", adj$disno[badv], adj$a[badv], adj$b[badv], adj$field[badv], adj$value_copy1[badv], adj$value_copy2[badv],
                                    dis$value_copy1[m][badv], dis$value_copy2[m][badv]), collapse = "; "))
  ## the rule allowed per field, and the applies_to each rule needs
  fld_rule <- list(init_dyadic = "direction_uninformative", war_dymid = "war_defer_to_war_data", fatlev_max = c("fatal_strict_lower", "fatal_ub_higher"),
                   fatleva_max = c("fatal_strict_lower", "fatal_ub_higher"), fatlevb_max = c("fatal_strict_lower", "fatal_ub_higher"), start_ym = "earliest_start", end_ym = "latest_end",
                   hihost = "hostility_higher", hihosta = "hostility_higher", hihostb = "hostility_higher", side_a_dy = "sides_uninformative", side_b_dy = "sides_uninformative")
  rule_applies <- c(fatal_strict_lower = "fmid", fatal_ub_higher = "fmidany,fmidub")
  for (i in seq_len(nrow(adj))) {
    ok <- adj$field[i] %in% names(fld_rule) && adj$rule[i] %in% fld_rule[[adj$field[i]]]
    if (!ok) stop(sprintf("dymid_dyads: rule %s is not a rule for field %s (dispute %d %d-%d); no resolution", adj$rule[i], adj$field[i], adj$disno[i], adj$a[i], adj$b[i]))
    want <- if (adj$rule[i] %in% names(rule_applies)) rule_applies[[adj$rule[i]]] else "all"
    if (adj$applies_to[i] != want) stop(sprintf("dymid_dyads: rule %s needs applies_to = %s, the table has %s (dispute %d %d-%d %s)", adj$rule[i], want, adj$applies_to[i], adj$disno[i], adj$a[i], adj$b[i], adj$field[i]))
  }
  ## a fatality disagreement needs BOTH rules (fmid and fmidany / fmidub)
  fk <- dis[dis$field %in% c("fatlev_max", "fatleva_max", "fatlevb_max"), ]
  for (i in seq_len(nrow(fk))) {
    r <- adj[key(adj) == key(fk)[i], ]
    if (!setequal(r$rule, c("fatal_strict_lower", "fatal_ub_higher"))) stop(sprintf("dymid_dyads: fatality disagreement dispute %d %d-%d %s needs one fatal_strict_lower row (fmid) and one fatal_ub_higher row (fmidany,fmidub)", fk$disno[i], fk$a[i], fk$b[i], fk$field[i]))
  }
  ## resolve: start from copy1, apply each row; recompute the resolved value from the rule and compare with the table
  out <- d1 |> mutate(fatlev_ub = fatlev_max, fatleva_ub = fatleva_max, fatlevb_ub = fatlevb_max, copy_adjudicated = FALSE, adj_fields = "")
  years <- attr(d1, "years"); spells1 <- attr(d1, "spells"); spells2 <- attr(d2, "spells")
  stopifnot(nrow(spells1) == nrow(spells2), all(spells1$disno == spells2$disno), all(spells1$a == spells2$a), all(spells1$b == spells2$b), all(spells1$spell == spells2$spell),
            all(spells1$y_first == spells2$y_first), all(spells1$y_last == spells2$y_last))
  ## a spell-level date disagreement must be the dyad-level one (first spell's start / the spell carrying the latest end), else there is no table row for it
  sk <- paste(spells1$disno, spells1$a, spells1$b)
  dk_s <- paste(dis$disno, dis$a, dis$b)[dis$field == "start_ym"]; dk_e <- paste(dis$disno, dis$a, dis$b)[dis$field == "end_ym"]
  s_dis <- spells1$spell_start_ym != spells2$spell_start_ym; e_dis <- spells1$spell_end_ym != spells2$spell_end_ym
  ## the end rule resolves ONE spell per dyad, the one carrying the latest end over both copies (which.max, the same choice as below); an end
  ## disagreement on any other spell has no row and stops (audit P2-2, 2026-11-08)
  emax <- pmax(spells1$spell_end_ym, spells2$spell_end_ym)
  kmax <- ave(seq_along(emax), sk, FUN = function(i) rep(i[which.max(emax[i])], length(i)))
  inner <- (s_dis & !(spells1$spell == 1L & sk %in% dk_s)) | (e_dis & !(seq_along(sk) == kmax & sk %in% dk_e))
  if (any(inner)) stop("dymid_dyads: a copy disagreement on the dates of a spell that is not the dyad's start or end (no adjudication field for it): ",
                       paste(sprintf("dispute %d %d-%d spell %d", spells1$disno[inner], spells1$a[inner], spells1$b[inner], spells1$spell[inner]), collapse = "; "))
  spells <- spells1
  adj$resolved_check <- rep(NA_character_, nrow(adj))
  for (i in seq_len(nrow(adj))) {
    j <- which(out$disno == adj$disno[i] & out$a == adj$a[i] & out$b == adj$b[i]); stopifnot(length(j) == 1L)
    f <- adj$field[i]; rl <- adj$rule[i]; v1 <- d1[[f]][j]; v2 <- d2[[f]][j]
    res <- switch(rl,
      direction_uninformative = NA_integer_,
      war_defer_to_war_data = { if (is.null(war_pairs)) stop("dymid_dyads: the rule war_defer_to_war_data needs war_pairs (the war data's dyads)")
        any(war_pairs$disno == adj$disno[i] & war_pairs$a == adj$a[i] & war_pairs$b == adj$b[i]) },
      fatal_strict_lower = c(v1, v2)[which.min(fat_rank(c(v1, v2)))],
      fatal_ub_higher = c(v1, v2)[which.max(fat_rank(c(v1, v2)))],
      earliest_start = min(v1, v2), latest_end = max(v1, v2), hostility_higher = max(v1, v2, na.rm = TRUE),
      sides_uninformative = NA_integer_)
    if (f %in% c("side_a_dy", "side_b_dy")) {
      adj$resolved_check[i] <- "NA"
      years[[f]][years$disno == adj$disno[i] & years$a == adj$a[i] & years$b == adj$b[i]] <- NA_integer_
    } else {
      adj$resolved_check[i] <- dymid_adj_format(f, res)
      if (rl == "fatal_ub_higher") out[[sub("_max$", "_ub", f)]][j] <- res else out[[f]][j] <- res
      if (f == "start_ym") { out$onset_year[j] <- ym_year(res); k <- spells$disno == adj$disno[i] & spells$a == adj$a[i] & spells$b == adj$b[i] & spells$spell == 1L; spells$spell_start_ym[k] <- res }
      if (f == "end_ym") { k <- spells$disno == adj$disno[i] & spells$a == adj$a[i] & spells$b == adj$b[i]; kk <- which(k)[which.max(emax[k])]; stopifnot(kmax[kk] == kk); spells$spell_end_ym[kk] <- res }
    }
    out$copy_adjudicated[j] <- TRUE; out$adj_fields[j] <- paste(unique(c(strsplit(out$adj_fields[j], ";")[[1]], f)), collapse = ";")
  }
  badr <- adj$resolved_value != adj$resolved_check
  if (any(badr)) stop("dymid_dyads: adjudication row(s) whose resolved_value does not follow from the rule: ",
                      paste(sprintf("dispute %d %d-%d %s rule %s (table %s, rule gives %s)", adj$disno[badr], adj$a[badr], adj$b[badr], adj$field[badr], adj$rule[badr], adj$resolved_value[badr], adj$resolved_check[badr]), collapse = "; "))
  stopifnot(all(ym_year(out$start_ym) == out$onset_year), all(out$start_ym <= out$end_ym), all(fat_rank(out$fatlev_ub) >= fat_rank(out$fatlev_max)),
            all(fat_rank(out$fatleva_ub) >= fat_rank(out$fatleva_max)), all(fat_rank(out$fatlevb_ub) >= fat_rank(out$fatlevb_max)),
            all(spells$spell_start_ym <= spells$spell_end_ym), sum(out$copy_adjudicated) == nrow(distinct(dis, disno, a, b)))
  attr(out, "years") <- years; attr(out, "spells") <- spells
  attr(out, "adjudication") <- adj |> select(-resolved_check)
  out
}
## ddiw_dyads(wr): one row per war dyad (warnum, disno, a < b) of the directed war data: war_start_ym / war_end_ym (the dyadic war's dates, first and
## last battle), war_onset_year, init_wardyadic (wardyadrolea == 1 -> a, == 3 -> b, from the earliest row), warole_a / warole_b (1 primary initiator,
## 2 joiner on the initiator's side, 3 primary target, 4 joiner on the target's side), dyindex (as recorded), batdths (max; logging only).
## attr "years": one row per (warnum, disno, a, b, year).
ddiw_dyads <- function(wr) {
  r <- wr[wr$statea < wr$stateb, ] |> rename(a = statea, b = stateb) |>
    mutate(ws = ym_index(warstrtyr, warstrtmnth), we = ym_index(warendyr, warendmnth))
  out <- r |> group_by(warnum, disno, a, b) |> arrange(ws, year, .by_group = TRUE) |>
    summarise(war_start_ym = min(ws), war_end_ym = max(we), war_onset_year = min(warstrtyr), y_first = min(year), y_last = max(year),
              init_wardyadic = ifelse(wardyadrolea[1] == 1L, a[1], b[1]), wrole_varies = n_distinct(wardyadrolea) > 1L,
              warole_a = warolea[1], warole_b = waroleb[1], dyindex = paste(unique(dyindex), collapse = ";"), batdths = max_na(batdths), .groups = "drop")
  stopifnot(all(out$y_first == out$war_onset_year), all(out$war_start_ym <= out$war_end_ym))
  attr(out, "years") <- r |> distinct(warnum, disno, a, b, year, ws, we) |> group_by(warnum, disno, a, b, year) |>
    summarise(start_ym = min(ws), end_ym = max(we), .groups = "drop")
  out
}
## dymid_link(dyads, wdy): ONE event dyad per (disno, a, b).
##   merged      : a war dyad of the war data whose (disno, a, b) is a DYMID dyadic dispute (the 4.03 dyindex is not a dyad key, so the link is
##                 disno + pair; the war data's decimal dyindex is kept for the log). One event: the DYMID dyadic dispute (its dates and onset, its
##                 dyadic initiator), flagged war, war_ddiw = TRUE, warnum_ddiw = the COW war number.
##   ddiw_only   : a war dyad with no DYMID dyadic dispute: an event of its own dated by the war start (onset = war start year), its dyadic war
##                 initiator for rule (b), its war roles for the originator / side fallback.
##   dymid_only  : every other DYMID dyadic dispute (war = war_dymid; a DYMID war == 1 dyad absent from the war data is still a war dyad).
## A (disno, pair) with two war dyads (two COW wars in one dispute) stops the build. Returns the event dyads with source, war, war_ddiw, war_dymid,
## warnum_ddiw, onset_year, start_ym, end_ym, init_b (the rule-(b) dyadic initiator: DYMID's for DYMID dyads, the war data's for ddiw_only),
## init_wardyadic (for the merged-dyad comparison), the DYMID fatality / role columns (NA for ddiw_only); attr "years" (disno, a, b, year,
## start_ym, end_ym, side_a_dy, side_b_dy; sides NA for war-only years) and attr "counts".
dymid_link <- function(dyads, wdy) {
  stopifnot(!anyDuplicated(dyads[, c("disno", "a", "b")]))
  if (anyDuplicated(wdy[, c("disno", "a", "b")])) stop("dymid_link: a (disno, pair) carries two war dyads of the war data; review")
  m <- dyads |> left_join(wdy |> select(disno, a, b, warnum_ddiw = warnum, war_start_ym, war_end_ym, war_onset_year, init_wardyadic, warole_a, warole_b, dyindex_ddiw = dyindex),
                          by = c("disno", "a", "b")) |>
    mutate(war_ddiw = !is.na(warnum_ddiw), source = ifelse(war_ddiw, "merged", "dymid_only"), init_b = init_dyadic)
  wo <- wdy |> anti_join(dyads, by = c("disno", "a", "b")) |>
    transmute(disno, a, b, onset_year = war_onset_year, start_ym = war_start_ym, end_ym = war_end_ym, y_first, y_last, n_spells = 1L, n_rows = NA_integer_,
              war_dymid = FALSE, war_ddiw = TRUE, warnum_ddiw = warnum, war_start_ym, war_end_ym, war_onset_year, init_wardyadic, warole_a, warole_b,
              dyindex_ddiw = dyindex, source = "ddiw_only", init_b = init_wardyadic, init_dyadic = NA_integer_, recip = FALSE, recip_code2 = FALSE, role_varies = FALSE,
              role_varies_first = FALSE, orig_varies = FALSE, copy_adjudicated = FALSE, adj_fields = "")
  out <- bind_rows(m, wo) |> mutate(war = war_dymid | war_ddiw) |> arrange(disno, a, b)
  stopifnot(!anyDuplicated(out[, c("disno", "a", "b")]), nrow(out) == nrow(dyads) + nrow(wo), !anyNA(out$war), !anyNA(out$onset_year),
            all(c("fatlev_ub", "fatleva_ub", "fatlevb_ub", "copy_adjudicated") %in% names(out)), !anyNA(out$copy_adjudicated))
  wy <- attr(wdy, "years") |> semi_join(wo, by = c("disno", "a", "b")) |> transmute(disno, a, b, year, start_ym, end_ym, side_a_dy = NA_integer_, side_b_dy = NA_integer_, n_side = 1L)
  attr(out, "years") <- bind_rows(attr(dyads, "years"), wy) |> arrange(disno, a, b, year)
  ## the spells of the event dyads (dymid_override_pdi()): DYMID's per dyad; a war-only dyad is one spell dated by the dyadic war
  ws <- wo |> transmute(disno, a, b, spell_start_ym = start_ym, spell_end_ym = end_ym, y_first, y_last, spell = 1L)
  attr(out, "spells") <- bind_rows(attr(dyads, "spells"), ws) |> arrange(disno, a, b, spell)
  stopifnot(nrow(distinct(attr(out, "spells"), disno, a, b)) == nrow(out))
  attr(out, "counts") <- list(n_dymid = nrow(dyads), n_war_dyads = nrow(wdy), n_merged = sum(out$source == "merged"), n_ddiw_only = nrow(wo),
                              n_dymid_only_war = sum(out$source == "dymid_only" & out$war_dymid), n_merged_dymid_war0 = sum(out$source == "merged" & !out$war_dymid),
                              n_merged_role_disagree = sum(out$source == "merged" & !is.na(out$init_dyadic) & out$init_dyadic != out$init_wardyadic))
  out
}

## dymid_participation(dd, evd, years): the participation table the event machinery of the dymid build reads (escalation and choice-set spells,
## participant entry dates, roles, joiner sets), in the schema of COW's directed dispute-years: dispnum, ccode1, ccode2, year, sidea1, sidea2, orig1,
## orig2 (both directions). It is COW's table (dd) with the (dispnum, pair) of every event dyad REPLACED by the event dyad's own years, and the event
## dyads COW lacks ADDED (DYMID-only dispute numbers, pairs absent from COW's dispute). COW pairs that are not event dyads (e.g. USA - Finland in
## dispute 258) stay as COW has them: they are co-participants of the dispute, so the already-fighting rule, the participant spells and the joiner
## sets keep their current definition. Per participant (dispnum, ccode):
##   originator flag : COW's orig1 when the state is a COW participant of that dispute; else DYMID orignat; else (a war-only dyad's state) its war
##                     role (1 / 3 primary = originator, 2 / 4 joiner)
##   side in a year  : COW's sidea1 in that year (a state with no COW row in that year: its side in its first COW year); a state that is not a COW
##                     participant: DYMID's side in the dispute's COW orientation (DYMID's Side A is inverted against COW's in some disputes; the
##                     orientation is the majority over the states both code), else the side opposite its partner; a war-only state: its war role
## Every event-dyad row must have opposite sides. A state with two COW sides in a year (dispute 258) takes the one opposite its partner; when the
## rule still yields the same side for both (DYMID and COW disagree on who opposed whom), the partner of the COW participant (or of the originator,
## or state b) takes the opposite side for that dyad-year (counted: n_side_conflict; it becomes a two-sided participant-year, which the role and
## joiner-set code already handle). Returns the table; attr "diagnostics".
dymid_participation <- function(dd, evd, years) {
  stopifnot(all(c("dispnum", "ccode1", "ccode2", "year", "sidea1", "sidea2", "orig1", "orig2") %in% names(dd)), all(c("disno", "a", "b", "warole_a", "warole_b", "orig_a_dy", "orig_b_dy") %in% names(evd)),
            all(c("disno", "a", "b", "year", "side_a_dy", "side_b_dy") %in% names(years)))
  cow <- dd |> transmute(dispnum = as.integer(dispnum), ccode1 = as.integer(ccode1), ccode2 = as.integer(ccode2), year = as.integer(year), sidea1 = as.integer(sidea1),
                         sidea2 = as.integer(sidea2), orig1 = as.integer(orig1), orig2 = as.integer(orig2))
  cp <- cow |> distinct(dispnum, ccode = ccode1, year, side = sidea1, orig = orig1)
  co <- cp |> group_by(dispnum, ccode) |> summarise(n_o = n_distinct(orig), orig_cow = first(orig), first_year = min(year), .groups = "drop")
  if (any(co$n_o != 1L)) stop("dymid_participation: COW originator status differs across a participant's dispute-years")
  cs_y <- cp |> group_by(dispnum, ccode, year) |> summarise(s1 = any(side == 1L), s0 = any(side == 0L), .groups = "drop")
  cs_f <- cp |> inner_join(co |> select(dispnum, ccode, first_year), by = c("dispnum", "ccode")) |> filter(year == first_year) |>
    group_by(dispnum, ccode) |> summarise(f1 = any(side == 1L), f0 = any(side == 0L), .groups = "drop")
  ## orientation of DYMID's sides against COW's, per dispute, over the states both code (first DYMID year of each such state)
  ylong <- bind_rows(years |> transmute(disno, ccode = a, year, side_dy = side_a_dy), years |> transmute(disno, ccode = b, year, side_dy = side_b_dy)) |> filter(!is.na(side_dy))
  orient <- ylong |> group_by(disno, ccode) |> summarise(side_dy = side_dy[which.min(year)], .groups = "drop") |>
    inner_join(cs_f, by = c("disno" = "dispnum", "ccode")) |> filter(xor(f1, f0)) |>
    group_by(disno) |> summarise(n_both = n(), n_agree = sum(side_dy == as.integer(f1)), .groups = "drop") |>
    mutate(inverted = 2L * n_agree < n_both, tie = 2L * n_agree == n_both)
  ## originator flags per participant
  dyo <- bind_rows(evd |> transmute(disno, ccode = a, o_dy = orig_a_dy, wr = warole_a), evd |> transmute(disno, ccode = b, o_dy = orig_b_dy, wr = warole_b)) |>
    group_by(disno, ccode) |> summarise(o_dy = max_na(o_dy), wr = max_na(wr), .groups = "drop")
  po <- dyo |> left_join(co |> select(dispnum, ccode, orig_cow), by = c("disno" = "dispnum", "ccode")) |>
    mutate(orig = case_when(!is.na(orig_cow) ~ orig_cow, !is.na(o_dy) ~ as.integer(o_dy), wr %in% c(1L, 3L) ~ 1L, wr %in% c(2L, 4L) ~ 0L, TRUE ~ NA_integer_),
           orig_src = case_when(!is.na(orig_cow) ~ "cow", !is.na(o_dy) ~ "dymid", TRUE ~ "war_role"))
  if (anyNA(po$orig)) stop("dymid_participation: a participant has no originator flag in any source")
  ## sides per event-dyad year
  y <- years |> left_join(evd |> select(disno, a, b, warole_a, warole_b), by = c("disno", "a", "b"))
  look <- function(cc) y |> select(disno, ccode = all_of(cc), year) |> left_join(cs_y, by = c("disno" = "dispnum", "ccode", "year")) |>
    left_join(cs_f, by = c("disno" = "dispnum", "ccode")) |> mutate(cow = !is.na(f1),
                                                                    c1 = ifelse(!is.na(s1), s1, f1), c0 = ifelse(!is.na(s0), s0, f0))
  la <- look("a"); lb <- look("b")
  inv <- coalesce(orient$inverted[match(y$disno, orient$disno)], FALSE)
  dy_a <- ifelse(is.na(y$side_a_dy), NA_integer_, ifelse(inv, 1L - y$side_a_dy, y$side_a_dy))
  dy_b <- ifelse(is.na(y$side_b_dy), NA_integer_, ifelse(inv, 1L - y$side_b_dy, y$side_b_dy))
  wr_a <- ifelse(y$warole_a %in% 1:2, 1L, ifelse(y$warole_a %in% 3:4, 0L, NA_integer_)); wr_b <- ifelse(y$warole_b %in% 1:2, 1L, ifelse(y$warole_b %in% 3:4, 0L, NA_integer_))
  ## a unique COW side -> that side; two COW sides -> NA here, resolved against the partner below; not a COW participant -> DYMID (oriented), else war role
  uni <- function(l) ifelse(l$cow & xor(l$c1, l$c0), as.integer(l$c1), NA_integer_)
  sa <- uni(la); sb <- uni(lb)
  two_a <- la$cow & la$c1 & la$c0; two_b <- lb$cow & lb$c1 & lb$c0
  sa <- ifelse(is.na(sa) & !la$cow, coalesce(dy_a, wr_a), sa); sb <- ifelse(is.na(sb) & !lb$cow, coalesce(dy_b, wr_b), sb)
  sa <- ifelse(is.na(sa) & two_a & !is.na(sb), 1L - sb, sa); sb <- ifelse(is.na(sb) & two_b & !is.na(sa), 1L - sa, sb)
  ## both two-sided (or still unknown): DYMID's oriented sides, else the opposite of the partner, else a = side A
  sa <- ifelse(is.na(sa), coalesce(dy_a, ifelse(is.na(sb), 1L, 1L - sb)), sa); sb <- ifelse(is.na(sb), coalesce(dy_b, 1L - sa), sb)
  stopifnot(!anyNA(sa), !anyNA(sb))
  conflict <- sa == sb
  ## anchor: the COW participant (one of the two), else the originator, else a
  oa <- po$orig[match(paste(y$disno, y$a), paste(po$disno, po$ccode))]; ob <- po$orig[match(paste(y$disno, y$b), paste(po$disno, po$ccode))]
  anchor_b <- (lb$cow & !la$cow) | (la$cow == lb$cow & ob == 1L & oa == 0L)
  sa_fix <- ifelse(conflict & anchor_b, 1L - sb, sa); sb_fix <- ifelse(conflict & !anchor_b, 1L - sa, sb)
  stopifnot(all(sa_fix != sb_fix))
  ev_rows <- y |> transmute(dispnum = disno, a, b, year, side_a = sa_fix, side_b = sb_fix) |>
    left_join(po |> select(dispnum = disno, a = ccode, orig_a = orig), by = c("dispnum", "a")) |> left_join(po |> select(dispnum = disno, b = ccode, orig_b = orig), by = c("dispnum", "b"))
  stopifnot(!anyNA(ev_rows$orig_a), !anyNA(ev_rows$orig_b))
  ev_dd <- bind_rows(ev_rows |> transmute(dispnum, ccode1 = a, ccode2 = b, year, sidea1 = side_a, sidea2 = side_b, orig1 = orig_a, orig2 = orig_b),
                     ev_rows |> transmute(dispnum, ccode1 = b, ccode2 = a, year, sidea1 = side_b, sidea2 = side_a, orig1 = orig_b, orig2 = orig_a))
  ## COW pairs that are event dyads are replaced (both directions); everything else of COW stays
  evk <- evd |> distinct(disno, a, b)
  keep_cow <- cow |> mutate(.a = pmin(ccode1, ccode2), .b = pmax(ccode1, ccode2)) |> anti_join(evk, by = c("dispnum" = "disno", ".a" = "a", ".b" = "b")) |> select(-.a, -.b)
  out <- bind_rows(keep_cow, ev_dd) |> arrange(dispnum, ccode1, ccode2, year)
  assert_no_dup_keys(out, c("dispnum", "ccode1", "ccode2", "year"), "dymid participation table")
  stopifnot(all(out$sidea1 + out$sidea2 == 1L), all(out$orig1 %in% 0:1), all(out$orig2 %in% 0:1))
  chk <- out |> group_by(dispnum, ccode1) |> summarise(n = n_distinct(orig1), .groups = "drop")
  if (any(chk$n != 1L)) stop("dymid_participation: originator status not constant per participant in the merged table")
  nk <- cow |> mutate(.a = pmin(ccode1, ccode2), .b = pmax(ccode1, ccode2)) |> distinct(dispnum, .a, .b)
  attr(out, "diagnostics") <- list(
    n_cow_rows = nrow(cow), n_cow_rows_kept = nrow(keep_cow), n_cow_pairs = nrow(nk), n_cow_pairs_replaced = nrow(semi_join(nk, evk, by = c("dispnum" = "disno", ".a" = "a", ".b" = "b"))),
    n_cow_pairs_not_event = nrow(anti_join(nk, evk, by = c("dispnum" = "disno", ".a" = "a", ".b" = "b"))),
    n_event_pairs_added = nrow(anti_join(evk, nk, by = c("disno" = "dispnum", "a" = ".a", "b" = ".b"))),
    n_disputes_added = length(setdiff(unique(evd$disno), unique(cow$dispnum))), n_event_rows = nrow(ev_rows),
    n_orig_src = table(po$orig_src), n_orig_disagree = sum(!is.na(po$orig_cow) & !is.na(po$o_dy) & po$orig_cow != po$o_dy),
    orig_disagree = po |> filter(!is.na(orig_cow), !is.na(o_dy), orig_cow != o_dy) |> select(disno, ccode, orig_cow, orig_dy = o_dy),
    n_inverted_disputes = sum(orient$inverted), n_orient_tie = sum(orient$tie), n_two_sided_resolved = sum(two_a | two_b),
    n_side_conflict = sum(conflict), side_conflict = y[conflict, c("disno", "a", "b", "year")])
  out
}

## dymid_disputes(disp_cow, dy, evd, dd_ev): the dispute-level table of every dispute in dd_ev, in the schema of cow_mid_disps that
## dispute_date_bounds() reads (dispnum, styear, stmon, mindur, maxdur, hiact, hostlev, fatality). COW's rows for the COW disputes; for a DYMID-only
## dispute number (DYMID split some COW disputes, e.g. 2581, 3587) a constructed row: start = its earliest dyad start (year, month), duration = days from
## the earliest start to the latest end (unknown days: the first / last of the month), hostlev = max hihost (5 when it holds a war dyad), fatality = max
## fatlev, hiact NA. attr "constructed" lists them. hostlev_ev = the event machinery's dispute hostility: max(COW hostlev, DYMID max hihost, 5 when the
## dispute holds a war dyad), used for the war-level dispute filter of D1 and war_disp.
dymid_disputes <- function(disp_cow, dy, evd, dd_ev) {
  stopifnot(all(c("dispnum", "styear", "stmon", "mindur", "maxdur", "hiact", "hostlev", "fatality") %in% names(disp_cow)))
  ids <- sort(unique(dd_ev$dispnum))
  cw <- disp_cow |> transmute(dispnum = as.integer(dispnum), styear = as.integer(styear), stmon = as.integer(stmon), mindur = as.numeric(mindur), maxdur = as.numeric(maxdur),
                              hiact = as.integer(hiact), hostlev = as.integer(hostlev), fatality = as.integer(fatality)) |> filter(dispnum %in% ids)
  dd_d <- function(y, m, d, last) { d0 <- as.Date(sprintf("%d-%02d-01", y, m)); if (!last) d0 + pmax(coalesce(d, 1L), 1L) - 1L else {
    nm <- as.Date(sprintf("%d-%02d-01", y + (m == 12L), ifelse(m == 12L, 1L, m + 1L))) - 1L; ifelse(!is.na(d) & d > 0L, d0 + d - 1L, nm) } }
  r <- dy |> mutate(sd = dd_d(strtyr, strtmnth, strtday, FALSE), ed = as.Date(dd_d(endyear, endmnth, endday, TRUE), origin = "1970-01-01"))
  dy_disp <- r |> group_by(dispnum = disno) |> summarise(sd = min(sd), ed = max(ed), hihost = max_na(hihost), fatlev = max_na(na9(fatlev)), .groups = "drop")
  wd <- evd |> filter(war) |> distinct(dispnum = disno) |> mutate(has_war = TRUE)
  new_ids <- setdiff(ids, cw$dispnum)
  ## a war-only dyad whose dispute number is in neither COW nor DYMID (1617, Cambodia (811) - North Vietnam / DRV (816), 1970): dated by the dyadic war
  wo <- evd |> filter(source == "ddiw_only", disno %in% setdiff(new_ids, dy_disp$dispnum))
  if (nrow(wo)) wo <- wo |> group_by(dispnum = disno) |>
    summarise(sd = as.Date(sprintf("%d-%02d-01", ym_year(min(start_ym)), min(start_ym) %% 12L + 1L)),
              ed = as.Date(sprintf("%d-%02d-01", ym_year(max(end_ym)), max(end_ym) %% 12L + 1L)), hihost = 5L, fatlev = NA_integer_, .groups = "drop")
  if (nrow(wo)) dy_disp <- bind_rows(dy_disp, wo)
  cons <- dy_disp |> filter(dispnum %in% new_ids) |>
    transmute(dispnum, styear = as.integer(format(sd, "%Y")), stmon = as.integer(format(sd, "%m")), mindur = as.numeric(ed - sd) + 1, maxdur = as.numeric(ed - sd) + 1,
              hiact = NA_integer_, hostlev = as.integer(hihost), fatality = coalesce(as.integer(fatlev), -9L)) |>
    left_join(wd, by = "dispnum") |> mutate(hostlev = ifelse(coalesce(has_war, FALSE), 5L, hostlev)) |> select(-has_war)
  if (length(setdiff(new_ids, cons$dispnum))) stop("dymid_disputes: disputes without any dispute-level source: ", paste(setdiff(new_ids, cons$dispnum), collapse = ", "))
  out <- bind_rows(cw, cons) |> arrange(dispnum)
  out <- out |> left_join(dy_disp |> select(dispnum, hihost), by = "dispnum") |> left_join(wd, by = "dispnum") |>
    mutate(hostlev_ev = pmax(hostlev, coalesce(as.integer(hihost), 0L), ifelse(coalesce(has_war, FALSE), 5L, 0L))) |> select(-hihost, -has_war)
  stopifnot(nrow(out) == length(ids), !anyNA(out$styear), all(out$stmon %in% 1:12), !anyNA(out$hostlev_ev), all(out$maxdur >= out$mindur))
  attr(out, "constructed") <- cons
  out
}

## dymid_override_pdi(pdi, spells): the pair-dispute spells of the event dyads (pair_dispute_intervals() on the participation table) replaced by the
## dyads' OWN spells (attr "spells" of dymid_link(): one per distinct DYMID start month within the dyad-dispute, dated by its start and its latest
## end; a war-only dyad's spell is the dyadic war), instead of the dispute / participant bounds. A pair_dispute_intervals() spell is a run of
## consecutive participation years, which BRIDGED a dyad that left the dispute and re-entered it in a later year (USA - Bulgaria in 258: Dec 1941 -
## Jun 1942 and Nov 1943 - Apr 1944 were one spell Dec 1941 - Apr 1944) until 2026-11-08: now each DYMID spell is its own row (spell_first /
## spell_last = the spell's years, exact months clamped into them, counted: the source dates sometimes disagree with the listed years), so the
## already-fighting rule, the choice-set reference dates and the escalation comparisons see the gap. Spells of one dyad may share a year (a pair
## leaving and re-entering within one year, 8 dyad-years): the year rule and the already-fighting keys take any() over them, as they do over the
## spells of a pair. `spell` is renumbered per directed pair in date order (spell 1 = the dyad's first spell = the event's date, pair_onset
## unchanged). date_source = "dymid", start_exact = TRUE. Every other spell (COW pairs that are not event dyads) is unchanged, as is the row
## order of the other pairs. attr: n_dymid (event-dyad rows after the split), n_clamped, n_split_pairs (directed pairs with more rows than before),
## n_split_rows (rows added).
dymid_override_pdi <- function(pdi, spells) {
  stopifnot(all(c("dispnum", "ccode1", "ccode2", "spell", "spell_first", "spell_last", "pair_onset", "start_lo", "start_hi", "end_lo", "end_hi", "start_exact", "date_source", "part_reason") %in% names(pdi)),
            all(c("disno", "a", "b", "spell", "spell_start_ym", "spell_end_ym", "y_first", "y_last") %in% names(spells)))
  sb <- bind_rows(spells |> transmute(dispnum = disno, ccode1 = a, ccode2 = b, sp = spell, sp_s = spell_start_ym, sp_e = spell_end_ym, y_first, y_last),
                  spells |> transmute(dispnum = disno, ccode1 = b, ccode2 = a, sp = spell, sp_s = spell_start_ym, sp_e = spell_end_ym, y_first, y_last))
  p <- pdi |> mutate(.r = row_number())
  ## the dyad spells overlapping each pair_dispute_intervals() run of an event dyad; a spell must lie inside one run, and the runs of a pair must
  ## be covered exactly by its spells (the participation years ARE the dyads' years)
  j <- p |> select(.r, dispnum, ccode1, ccode2, spell_first, spell_last) |> inner_join(sb, by = c("dispnum", "ccode1", "ccode2"), relationship = "many-to-many") |>
    filter(y_last >= spell_first, y_first <= spell_last)
  stopifnot(all(j$y_first >= j$spell_first), all(j$y_last <= j$spell_last))
  cov <- j |> group_by(.r, spell_first, spell_last) |> summarise(n_y = n_distinct(unlist(mapply(seq.int, y_first, y_last, SIMPLIFY = FALSE))), .groups = "drop")
  stopifnot(all(cov$n_y == cov$spell_last - cov$spell_first + 1L))
  ev_rows <- p |> semi_join(sb, by = c("dispnum", "ccode1", "ccode2"))
  if (!setequal(ev_rows$.r, unique(j$.r))) stop("dymid_override_pdi: an event dyad's pair-dispute run has no dyad spell")
  if (nrow(anti_join(sb, p, by = c("dispnum", "ccode1", "ccode2")))) stop("dymid_override_pdi: a dyad spell belongs to no pair-dispute of the participation table")
  new <- p[j$.r, ] |> mutate(spell_first = j$y_first, spell_last = j$y_last,
                             s_c = pmin(pmax(j$sp_s, ym_index(spell_first, 1L)), ym_index(spell_first, 12L)), e_c = pmin(pmax(j$sp_e, ym_index(spell_last, 1L)), ym_index(spell_last, 12L)),
                             e_c = pmax(e_c, s_c), clamped = s_c != j$sp_s | e_c != j$sp_e,
                             start_lo = s_c, start_hi = s_c, end_lo = e_c, end_hi = e_c, start_exact = TRUE, date_source = "dymid",
                             part_reason = ifelse(clamped, "dymid dyad dates (clamped into the spell's years)", "dymid dyad dates"), .sp = j$sp)
  if ("start_clipped" %in% names(new)) new$start_clipped <- FALSE
  n_clamped <- sum(new$clamped); new <- new |> select(-s_c, -e_c, -clamped)
  stopifnot(nrow(distinct(new, dispnum, ccode1, ccode2, .sp)) == nrow(new))
  ## renumber the event dyads' spells in date order; the other pairs keep their rows and numbering; the original row order (dispnum, ccode1, ccode2, spell)
  new <- new |> arrange(dispnum, ccode1, ccode2, .sp) |> group_by(dispnum, ccode1, ccode2) |> mutate(spell = row_number()) |> ungroup() |> select(-.sp)
  out <- bind_rows(p[!p$.r %in% ev_rows$.r, ], new) |> arrange(dispnum, ccode1, ccode2, spell) |> select(-.r)
  stopifnot(all(out$start_lo <= out$start_hi), all(out$start_hi <= out$end_hi), all(out$end_lo <= out$end_hi),
            all(out$spell_first[out$spell == 1L] == out$pair_onset[out$spell == 1L]), nrow(out) == nrow(pdi) - nrow(ev_rows) + nrow(new),
            nrow(distinct(out, dispnum, ccode1, ccode2)) == nrow(distinct(pdi, dispnum, ccode1, ccode2)))
  assert_no_dup_keys(out, c("dispnum", "ccode1", "ccode2", "spell"), "dymid pair dispute intervals")
  ## the union of years per pair is unchanged (every year of the participation table stays a spell year)
  yrs <- function(d) d |> group_by(dispnum, ccode1, ccode2) |> summarise(y = list(sort(unique(unlist(mapply(seq.int, spell_first, spell_last, SIMPLIFY = FALSE))))), .groups = "drop")
  stopifnot(identical(yrs(out), yrs(pdi)))
  n_before <- ev_rows |> count(dispnum, ccode1, ccode2, name = "n0"); n_after <- new |> count(dispnum, ccode1, ccode2, name = "n1")
  split <- inner_join(n_before, n_after, by = c("dispnum", "ccode1", "ccode2")) |> filter(n1 > n0)
  attr(out, "n_dymid") <- nrow(new); attr(out, "n_clamped") <- n_clamped; attr(out, "n_split_pairs") <- nrow(split); attr(out, "n_split_rows") <- nrow(new) - nrow(ev_rows)
  attr(out, "split_pairs") <- split
  out
}

## dymid_fatal_codes(fatleva, fatlevb, fatlev, war): the fatality codings of the dymid samples (fatality levels already max-over-years, NA = missing)
##   strict (fmid, MAIN) : fatleva >= 1 AND fatlevb >= 1 (both states suffered deaths in this dyad's fighting), OR a war dyad; missing -> not fatal
##   any    (fmidany)    : fatlev >= 1 (the dyadic dispute's fatality, either side) or fatleva >= 1 or fatlevb >= 1, OR a war dyad
##   ub     (fmidub)     : upper bound for the missing per-state fields (PI decision 2026-11-05): the strict rule, OR a war dyad, OR one side >= 1
##                         and the other side missing, OR both sides missing and the dyadic fatlev >= 1 or missing. A dyad whose dyadic fatlev is
##                         0 (coded non-fatal) is never made fatal by a missing side (it is fatal under ub only through the strict rule or war).
## Returns a data.frame of three logicals, never NA, plus ub_branch (the first branch that makes a dyad ub-fatal: "strict", "one_side_missing",
## "both_missing", or "none"); strict implies any and ub (asserted).
dymid_fatal_codes <- function(fatleva, fatlevb, fatlev, war) {
  n <- length(war); stopifnot(length(fatleva) == n, length(fatlevb) == n, length(fatlev) == n, !anyNA(war))
  f1 <- function(v) v %in% 1:6
  strict <- war | (f1(fatleva) & f1(fatlevb))
  any_ <- war | f1(fatlev) | f1(fatleva) | f1(fatlevb)
  dy_not0 <- !(fatlev %in% 0L)                                                    # the dyadic fatlev is not coded 0 (>= 1 or missing)
  one_miss <- dy_not0 & ((f1(fatleva) & is.na(fatlevb)) | (is.na(fatleva) & f1(fatlevb)))
  both_miss <- is.na(fatleva) & is.na(fatlevb) & (f1(fatlev) | is.na(fatlev))
  ub <- strict | one_miss | both_miss
  stopifnot(!anyNA(strict), !anyNA(any_), !anyNA(ub), all(any_[strict]), all(ub[strict]), !any(ub & !strict & fatlev %in% 0L))
  data.frame(strict = strict, any = any_, ub = ub,
             ub_branch = ifelse(strict, "strict", ifelse(one_miss, "one_side_missing", ifelse(both_miss, "both_missing", "none"))))
}

## dymid_direction(ev, psd): the initiator of every event dyad. ev = one row per event dyad oriented COW Side A -> Side B in its onset year (dispnum,
## initiator, target, year, sidea1 = 1, sidea2 = 0, orig1, orig2, fatality1 / 2, hostlev1 / 2) with war (logical), d1_init (the D1 initiator from
## cow_war_direction(): the COW war-initiating side's state, NA when D1 does not decide) and init_b (the rule-(b) dyadic initiator, NA when simultaneous).
## Rules, in the PI's order (PI decision 2026-11-05: (a) -> (b) -> (c)):
##   (a) the entry rule of section 6.8 adapted to dyads (entry_later_direction(), the existing participant-date logic): exactly one joiner -> the joiner
##       initiates; two joiners -> the later entrant initiates; an exact tie falls through
##   (b) both originators (or a joiner tie): the dataset's dyadic initiator (DYMID dyad_rolea == 1 / the war data's wardyadrolea == 1 for war-only dyads)
##   (c) only when (a) and (b) cannot decide (no dyadic initiator) and the dyad is a war dyad: the COW war-initiating side (D1, as in the pairs build)
##   none decides (no dyadic initiator: simultaneous first action, or the two directed copies contradict each other (adjudication rule
##       direction_uninformative, init_b = NA); and no D1 decision): COW Side A kept, listed as fallback_sideA
## Until 2026-11-05 (c) ran before (b) for war dyads; dir_changed marks the events whose initiator differs from what that order gave (a war dyad
## where D1 and (b) disagree and (a) does not decide). No half-weight / both-directions option (PI declined). Returns ev reoriented, with dir_source
## (one of names(DYMID_DIR_LABELS)), dir_rule (the same plus " (joiner tie)" for ties), dir_flipped (= the initiator is not on COW Side A),
## entry_reversed / entry_rule (entry_later_direction()), b_c_disagree (D1 decides and differs from (b), any dyad), dir_changed.
dymid_direction <- function(ev, psd) {
  stopifnot(all(c("dispnum", "initiator", "target", "year", "sidea1", "sidea2", "orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2", "war", "d1_init", "init_b") %in% names(ev)),
            all(ev$sidea1 == 1L), all(ev$sidea2 == 0L), !anyNA(ev$war))
  stopifnot(all(is.na(ev$d1_init) | ev$d1_init == ev$initiator | ev$d1_init == ev$target), all(is.na(ev$init_b) | ev$init_b == ev$initiator | ev$init_b == ev$target))
  pre <- case_when(!is.na(ev$init_b) ~ "b_dyadic", ev$war & !is.na(ev$d1_init) ~ "c_cowwar", TRUE ~ "fallback_sideA")
  init0 <- ifelse(pre == "c_cowwar", ev$d1_init, ifelse(pre == "b_dyadic", ev$init_b, ev$initiator))
  out <- ev; f <- init0 == ev$target
  swap <- function(x, y) { tmp <- out[[x]][f]; out[[x]][f] <<- out[[y]][f]; out[[y]][f] <<- tmp }
  swap("initiator", "target"); swap("sidea1", "sidea2"); swap("orig1", "orig2"); swap("fatality1", "fatality2"); swap("hostlev1", "hostlev2")
  stopifnot(all(out$initiator == init0))
  out <- entry_later_direction(out, psd)
  out$dir_source <- ifelse(out$entry_rule == "joiner_orig", "a_joiner_orig", ifelse(out$entry_rule == "joiner_joiner_later", "a_joiner_later", pre))
  out$dir_rule <- ifelse(out$entry_rule == "joiner_joiner_tie", paste0(pre, " (joiner tie)"), out$dir_source)
  out$dir_flipped <- out$sidea1 == 0L
  out$b_c_disagree <- !is.na(ev$d1_init) & !is.na(ev$init_b) & ev$d1_init != ev$init_b
  ## the previous order ((c) before (b) for war dyads) gave the D1 initiator wherever (a) did not decide; the entry rule (a) does not depend on the
  ## pre-orientation except through ties, which fall through to the pre-rule, so the change is exactly: war, disagreement, (a) not deciding
  out$dir_changed <- ev$war & out$b_c_disagree & !out$dir_source %in% c("a_joiner_orig", "a_joiner_later")
  stopifnot(all(out$initiator[out$dir_changed] == ev$init_b[out$dir_changed]))
  stopifnot(all(out$dir_source %in% names(DYMID_DIR_LABELS)), !any(out$orig1 == 1L & out$orig2 == 0L), all(out$sidea1 + out$sidea2 == 1L))
  out
}

## dymid_ww2(evd, years, war_inter, ddiw): the DYMID-only dispute numbers (not in COW, so not in WW2_DISPNUMS) that belong to the Second World War:
## a dyad-year in COW inter-state war 139 or a war-data dyad of war 139, and the dispute starts in 1939-09..1945-09 (the window of ww2_disputes()).
## Returns a data.frame (disno, start, n_cow139, n_ddiw139, is_ww2).
dymid_ww2 <- function(evd, years, war_inter, cow_ids) {
  w139 <- war_inter |> filter(warnum == 139) |> transmute(a = pmin(as.integer(ccode1), as.integer(ccode2)), b = pmax(as.integer(ccode1), as.integer(ccode2)), year = as.integer(year)) |> distinct()
  new <- evd |> filter(!disno %in% cow_ids) |> group_by(disno) |> summarise(start = min(start_ym), n_ddiw139 = sum(warnum_ddiw %in% 139L), .groups = "drop")
  n139 <- years |> filter(!disno %in% cow_ids) |> semi_join(w139, by = c("a", "b", "year")) |> count(disno, name = "n_cow139")
  new |> left_join(n139, by = "disno") |> mutate(n_cow139 = coalesce(n_cow139, 0L),
                                                  is_ww2 = start >= ym_index(1939L, 9L) & start <= ym_index(1945L, 9L) & (n_cow139 > 0L | n_ddiw139 > 0L))
}
