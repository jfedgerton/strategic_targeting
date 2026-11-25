## tests/test_grid_provenance.R -- grid positions carry the hashes of grid.R and grid_layers.rds (GPT P2-5); a mismatch stops before fitting
source("R/lib/grid.R"); td <- tempfile(); dir.create(td); gl <- file.path(td, "grid_layers.rds")
saveRDS(data.frame(a = 1:3), gl); pz <- list(provenance = grid_provenance(gl))
stopifnot(isTRUE(grid_check_provenance(pz, gl)))
saveRDS(data.frame(a = 1:4), gl)                                                  # layers rebuilt after the positions
r <- tryCatch(grid_check_provenance(pz, gl), error = function(e) conditionMessage(e)); stopifnot(grepl("different grid.R or grid_layers", r))
r <- tryCatch(grid_check_provenance(list(), gl), error = function(e) conditionMessage(e)); stopifnot(grepl("no provenance", r))
cat("TEST_GRID_PROVENANCE_OK\n")
