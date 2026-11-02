# Single source of truth for the FROZEN aligned-search grid (revision step 1.2, frozen
# 11 Oct 2026 after the two pre-specified expansions; STEP_1.2_INTEGRATION_SPEC.md, run 3).
# Every alignment run after that date — the regenerated DB1 rerun, FVC2002 DB3 — uses this.
FROZEN_GRID <- list(theta_max = 50, t_max = 40)
apply_frozen_grid <- function(config) {
  config$aligned$theta_max <- FROZEN_GRID$theta_max
  config$aligned$t_max     <- FROZEN_GRID$t_max
  config
}
