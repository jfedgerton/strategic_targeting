## ============================================================================
## helpers_dymid.R -- the "interacting dyads" event foundation (switch events = dymid, PI-approved 2026-11-04; README section 6.9)
##
## Functions only (sourced by 01_build_data.R under --events=dymid and by tests/testthat/test_dymid.R). Sources, both in v2/data_raw/dymid
## (provenance in v2/data_raw/SOURCE.md, sha256 checked by the build against the constants below):
##   dyadic_mid_4.03.csv     Maoz et al., Dyadic MID 4.03 (1816-2014): one row per DIRECTED dyad-year of a dyadic dispute, both directions
##                           listed (A, B) and (B, A); only directly interacting pairs. In 4.03 the column dyindex is an integer equal to disno
##                           on every row (the .001 dyad suffix of 4.01 is lost), so it cannot identify a dyad: the dyadic dispute is keyed by
##                           (disno, unordered pair) here, and a pair re-entering the same dispute (separate dyindices in 4.01; 27 pairs, mostly
##                           WWI / WWII) is ONE event dated by its first spell (the later spells stay in the already-fighting spells).
##   directed_dyadic_war.csv Maoz, Directed Dyadic Interstate War data (1816-2010): one row per directed war dyad-year (warnum = COW war number,
##                           disno, dyindex with the decimal suffix). A dyadic war = direct battles with >= 1,000 battle deaths across both states.
##                           batdtha / batdthb are sometimes national totals: never used as dyadic magnitudes.
## ============================================================================
DYMID_CSV_SHA256 <- "9d26eedd792a0067f590f038c3b7b240bd373cd70200ee729241e1078c6957c3"   # data_raw/dymid/dyadic_mid_4.03.csv
DDIW_CSV_SHA256  <- "aa2eb5244a0cb0b579867e50eb02e4297fb29250e67120fb132e3c0618dceda0"   # data_raw/dymid/directed_dyadic_war.csv
DYMID_DIR_LABELS <- c(a_joiner_orig = "(a) exactly one joiner: the joiner initiates", a_joiner_later = "(a) two joiners: the later entrant initiates",
                      c_cowwar = "(c) war dyad, both originators (or a joiner tie): COW war-initiating side (D1)",
                      b_dyadic = "(b) both originators (or a joiner tie): the dataset's dyadic initiator (dyad_rolea / wardyadrolea = 1)",
                      fallback_sideA = "no rule decides (simultaneous first action, no D1 decision): COW Side A kept")

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
  if (!all(x$dyindex == x$disno)) stop("read_dymid: dyindex differs from disno on some rows; the (disno, pair) key of helpers_dymid.R assumes the 4.03 truncation, review")
  ## both directions of every dyad-year: the multiset of (disno, a, b, year) keys of the A < B rows equals that of the A > B rows. The start month of
  ## the two copies differs on a few rows (3250 Russia - China 1900: July / August); the A < B copy is used, the count is in attr "n_start_mismatch"
  k <- function(d, s = FALSE) sort(paste(d$disno, pmin(d$statea, d$stateb), pmax(d$statea, d$stateb), d$year, if (s) paste(d$strtyr, d$strtmnth) else ""))
  if (!identical(k(x[x$statea < x$stateb, ]), k(x[x$statea > x$stateb, ]))) stop("read_dymid: the two directions of the dyad-years do not mirror each other")
  attr(x, "n_start_mismatch") <- length(setdiff(k(x[x$statea < x$stateb, ], TRUE), k(x[x$statea > x$stateb, ], TRUE)))
  x
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

## dymid_dyads(dy): one row per dyadic dispute (disno, a < b) from the A < B rows of read_dymid() (the A > B rows are their mirror). Columns:
##   onset_year / start_ym       the dyad's own first start (dyad-level dates: strtyr / strtmnth of its earliest spell; month index ym_index)
##   end_ym                      the latest dyad end (endyear / endmnth) over its rows
##   n_spells                    distinct start months of the dyad within the dispute (> 1: the pair re-entered; one event, first spell)
##   fatlev / fatleva / fatlevb  the MAXIMUM non-missing value over the dyad's rows (the codebook lists them for the last year of the dyadic dispute
##                               only; -9 = missing -> NA); *_last = the value in the dyad's last year (max over that year's rows), for the log
##   war_dymid                   any row with war == 1
##   init_dyadic                 the dyadic initiator of the FIRST spell (dyad_rolea == 1 -> a, == 3 -> b, == 2 simultaneous -> NA); NA as well when
##                               the first-spell rows disagree (counted: role_varies_first); role_varies = dyad_rolea differs over all the dyad's rows
##   orig_a_dy / orig_b_dy       DYMID orignata / orignatb (max over the rows; orig_varies counts disagreement)
##   recip                       any row with recip == 1 (the codebook: initiator unidentifiable, roles assigned randomly)
##   hihosta / hihostb / hihost  max over the rows
## attr "years": one row per (disno, a, b, year) with the row's spell start / end month and the DYMID sides of a / b (sideaa / sideab).
dymid_dyads <- function(dy) {
  r <- dy[dy$statea < dy$stateb, ] |> rename(a = statea, b = stateb) |>
    mutate(rs = ym_index(strtyr, strtmnth), re = ym_index(endyear, endmnth), fatlev = na9(fatlev), fatleva = na9(fatleva), fatlevb = na9(fatlevb))
  ## row-level dates are rs / re so that no summarise() expression below sees an already-reduced column
  out <- r |> group_by(disno, a, b) |> arrange(rs, year, .by_group = TRUE) |>
    summarise(onset_year = min(strtyr), start_ym = min(rs), end_ym = max(re), y_first = min(year), y_last = max(year),
              n_spells = n_distinct(rs), n_rows = n(),
              fatlev_max = max_na(fatlev), fatleva_max = max_na(fatleva), fatlevb_max = max_na(fatlevb),
              fatlev_last = max_na(fatlev[year == max(year)]), fatleva_last = max_na(fatleva[year == max(year)]), fatlevb_last = max_na(fatlevb[year == max(year)]),
              war_dymid = any(war == 1L),
              role_first = list(unique(dyad_rolea[rs == min(rs)])), role_varies = n_distinct(dyad_rolea) > 1L,
              orig_a_dy = max(orignata), orig_b_dy = max(orignatb), orig_varies = n_distinct(orignata) > 1L | n_distinct(orignatb) > 1L,
              recip = any(recip == 1L), hihosta = max_na(hihosta), hihostb = max_na(hihostb), hihost = max_na(hihost), .groups = "drop")
  out <- out |> mutate(role_varies_first = lengths(role_first) > 1L,
                       role1 = vapply(role_first, function(z) if (length(z) == 1L) as.integer(z) else NA_integer_, 1L),
                       init_dyadic = case_when(role1 == 1L ~ a, role1 == 3L ~ b, TRUE ~ NA_integer_)) |> select(-role_first, -role1)
  stopifnot(nrow(out) == nrow(distinct(r, disno, a, b)), all(out$onset_year == out$y_first), all(ym_year(out$start_ym) == out$onset_year), !anyNA(out$war_dymid))
  attr(out, "years") <- r |> group_by(disno, a, b, year) |>
    summarise(start_ym = min(rs), end_ym = max(re), side_a_dy = as.integer(sideaa[1]), side_b_dy = as.integer(sideab[1]),
              n_side = n_distinct(sideaa), .groups = "drop")
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
              dyindex_ddiw = dyindex, source = "ddiw_only", init_b = init_wardyadic, init_dyadic = NA_integer_, recip = FALSE, role_varies = FALSE,
              role_varies_first = FALSE, orig_varies = FALSE)
  out <- bind_rows(m, wo) |> mutate(war = war_dymid | war_ddiw) |> arrange(disno, a, b)
  stopifnot(!anyDuplicated(out[, c("disno", "a", "b")]), nrow(out) == nrow(dyads) + nrow(wo), !anyNA(out$war), !anyNA(out$onset_year))
  wy <- attr(wdy, "years") |> semi_join(wo, by = c("disno", "a", "b")) |> transmute(disno, a, b, year, start_ym, end_ym, side_a_dy = NA_integer_, side_b_dy = NA_integer_, n_side = 1L)
  attr(out, "years") <- bind_rows(attr(dyads, "years"), wy) |> arrange(disno, a, b, year)
  attr(out, "counts") <- list(n_dymid = nrow(dyads), n_war_dyads = nrow(wdy), n_merged = sum(out$source == "merged"), n_ddiw_only = nrow(wo),
                              n_dymid_only_war = sum(out$source == "dymid_only" & out$war_dymid), n_merged_dymid_war0 = sum(out$source == "merged" & !out$war_dymid),
                              n_merged_role_disagree = sum(out$source == "merged" & !is.na(out$init_dyadic) & out$init_dyadic != out$init_wardyadic))
  out
}

## dymid_participation(dd, evd, years): the participation table the event machinery of the dymid build reads (escalation and choice-set spells,
## participant entry dates, roles, joiner sets), in the schema of COW's directed dispute-years: dispnum, ccode1, ccode2, year, sidea1, sidea2, orig1,
## orig2 (both directions). It is COW's table (dd) with the (dispnum, pair) of every event dyad REPLACED by the event dyad's own years, and the event
## dyads COW lacks ADDED (DYMID-only dispute numbers, pairs absent from COW's dispute). COW pairs that are not event dyads (e.g. Japan - Poland in
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
  ## a war-only dyad whose dispute number is in neither COW nor DYMID (1617, Cambodia - South Vietnam 1970): dated by the dyadic war
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

## dymid_override_pdi(pdi, years): the pair-dispute spells of the event dyads (pair_dispute_intervals() on the participation table) dated by the
## dyads' own dates instead of the dispute / participant bounds. A spell (run of consecutive years) of an event dyad starts at the earliest DYMID
## spell start among its years and ends at the latest end (war-only dyads: the dyadic war's dates), exact months, clamped into the spell's first /
## last year (counted: the source dates sometimes disagree with the listed years). date_source = "dymid", start_exact = TRUE. Every other spell
## (COW pairs that are not event dyads) is unchanged.
dymid_override_pdi <- function(pdi, years) {
  stopifnot(all(c("dispnum", "ccode1", "ccode2", "spell", "spell_first", "spell_last", "start_lo", "start_hi", "end_lo", "end_hi", "start_exact", "date_source", "part_reason") %in% names(pdi)))
  yb <- bind_rows(years |> transmute(dispnum = disno, ccode1 = a, ccode2 = b, year, start_ym, end_ym), years |> transmute(dispnum = disno, ccode1 = b, ccode2 = a, year, start_ym, end_ym))
  p <- pdi |> mutate(.r = row_number())
  j <- p |> select(.r, dispnum, ccode1, ccode2, spell_first, spell_last) |> inner_join(yb, by = c("dispnum", "ccode1", "ccode2"), relationship = "many-to-many") |>
    filter(year >= spell_first, year <= spell_last) |> group_by(.r, spell_first, spell_last) |>
    summarise(s = min(start_ym), e = max(end_ym), n_y = n_distinct(year), .groups = "drop") |>
    mutate(s_c = pmin(pmax(s, ym_index(spell_first, 1L)), ym_index(spell_first, 12L)), e_c = pmin(pmax(e, ym_index(spell_last, 1L)), ym_index(spell_last, 12L)),
           e_c = pmax(e_c, s_c), clamped = s_c != s | e_c != e)
  stopifnot(all(j$n_y == j$spell_last - j$spell_first + 1L))
  p$start_lo[j$.r] <- j$s_c; p$start_hi[j$.r] <- j$s_c; p$end_lo[j$.r] <- j$e_c; p$end_hi[j$.r] <- j$e_c
  p$start_exact[j$.r] <- TRUE; p$date_source[j$.r] <- "dymid"; p$part_reason[j$.r] <- ifelse(j$clamped, "dymid dyad dates (clamped into the spell's years)", "dymid dyad dates")
  if ("start_clipped" %in% names(p)) p$start_clipped[j$.r] <- FALSE
  stopifnot(all(p$start_lo <= p$start_hi), all(p$start_hi <= p$end_hi), all(p$end_lo <= p$end_hi))
  attr(p, "n_dymid") <- nrow(j); attr(p, "n_clamped") <- sum(j$clamped)
  p |> select(-.r)
}

## dymid_fatal_codes(fatleva, fatlevb, fatlev, war): the fatality codings of the dymid samples (fatality levels already max-over-years, NA = missing)
##   strict (fmid, MAIN) : fatleva >= 1 AND fatlevb >= 1 (both states suffered deaths in this dyad's fighting), OR a war dyad; missing -> not fatal
##   any    (fmidany)    : fatlev >= 1 (the dyadic dispute's fatality, either side) or fatleva >= 1 or fatlevb >= 1, OR a war dyad
##   ub     (fmidub)     : the strict rule with a missing fatleva / fatlevb counted as fatal (bounds check), OR a war dyad
## Returns a data.frame of three logicals, never NA; strict implies any and ub.
dymid_fatal_codes <- function(fatleva, fatlevb, fatlev, war) {
  n <- length(war); stopifnot(length(fatleva) == n, length(fatlevb) == n, length(fatlev) == n, !anyNA(war))
  strict <- war | (fatleva %in% 1:6 & fatlevb %in% 1:6)
  any_ <- war | fatlev %in% 1:6 | fatleva %in% 1:6 | fatlevb %in% 1:6
  ub <- war | ((fatleva %in% 1:6 | is.na(fatleva)) & (fatlevb %in% 1:6 | is.na(fatlevb)))
  stopifnot(all(any_[strict]), all(ub[strict]))
  data.frame(strict = strict, any = any_, ub = ub)
}

## dymid_direction(ev, psd): the initiator of every event dyad. ev = one row per event dyad oriented COW Side A -> Side B in its onset year (dispnum,
## initiator, target, year, sidea1 = 1, sidea2 = 0, orig1, orig2, fatality1 / 2, hostlev1 / 2) with war (logical), d1_init (the D1 initiator from
## cow_war_direction(): the COW war-initiating side's state, NA when D1 does not decide) and init_b (the rule-(b) dyadic initiator, NA when simultaneous).
## Rules, in the PI's order:
##   (a) the entry rule of section 6.8 adapted to dyads (entry_later_direction(), the existing participant-date logic): exactly one joiner -> the joiner
##       initiates; two joiners -> the later entrant initiates; an exact tie falls through
##   (c) war dyads with both originators (or a joiner tie): the COW war-initiating side (D1, as in the pairs build) where D1 decides
##   (b) otherwise: the dataset's dyadic initiator (DYMID dyad_rolea == 1 / the war data's wardyadrolea == 1 for war-only dyads)
##   neither decides (simultaneous first action and no D1 decision): COW Side A kept, listed as fallback_sideA
## (c) is placed before (b) for war dyads so a war dyad keeps the direction the pairs and entry builds give it ("as now"); the build logs how often (b)
## and (c) disagree. Returns ev reoriented, with dir_source (one of names(DYMID_DIR_LABELS)), dir_rule (the same plus " (joiner tie)" for ties),
## dir_flipped (= the initiator is not on COW Side A), entry_reversed / entry_rule (entry_later_direction()), init_c / init_b kept for the log.
dymid_direction <- function(ev, psd) {
  stopifnot(all(c("dispnum", "initiator", "target", "year", "sidea1", "sidea2", "orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2", "war", "d1_init", "init_b") %in% names(ev)),
            all(ev$sidea1 == 1L), all(ev$sidea2 == 0L), !anyNA(ev$war))
  stopifnot(all(is.na(ev$d1_init) | ev$d1_init == ev$initiator | ev$d1_init == ev$target), all(is.na(ev$init_b) | ev$init_b == ev$initiator | ev$init_b == ev$target))
  pre <- case_when(ev$war & !is.na(ev$d1_init) ~ "c_cowwar", !is.na(ev$init_b) ~ "b_dyadic", TRUE ~ "fallback_sideA")
  init0 <- ifelse(pre == "c_cowwar", ev$d1_init, ifelse(pre == "b_dyadic", ev$init_b, ev$initiator))
  out <- ev; f <- init0 == ev$target
  swap <- function(x, y) { tmp <- out[[x]][f]; out[[x]][f] <<- out[[y]][f]; out[[y]][f] <<- tmp }
  swap("initiator", "target"); swap("sidea1", "sidea2"); swap("orig1", "orig2"); swap("fatality1", "fatality2"); swap("hostlev1", "hostlev2")
  stopifnot(all(out$initiator == init0))
  out <- entry_later_direction(out, psd)
  out$dir_source <- ifelse(out$entry_rule == "joiner_orig", "a_joiner_orig", ifelse(out$entry_rule == "joiner_joiner_later", "a_joiner_later", pre))
  out$dir_rule <- ifelse(out$entry_rule == "joiner_joiner_tie", paste0(pre, " (joiner tie)"), out$dir_source)
  out$dir_flipped <- out$sidea1 == 0L
  out$b_c_disagree <- ev$war & !is.na(ev$d1_init) & !is.na(ev$init_b) & ev$d1_init != ev$init_b
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
