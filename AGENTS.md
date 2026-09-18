# Codex guide for recurSurvTests

Before modifying statistical code in this repository:

1. Read `docs/chat-context.md` for the statistical motivation, method distinctions,
   source papers, and known implementation caveats.
2. Preserve the distinction between:
   - Wang-Chang marginal recurrent **gap-time survival**;
   - mean cumulative recurrent-event functions over calendar time; and
   - recurrent cause-specific cumulative-incidence functions.
3. Do not call the METRON RCIF test a log-rank test.
4. Do not call `wc_logrank()` an exact Zhao et al. (2020) implementation. It is
   the Luo-Huang `G_rho*` implementation based on a Wang-Chang-equivalent
   weighted risk set.
5. The Zhao et al. supplement available during the original chat contains
   simulations/residual diagnostics but no source code and not enough formulas
   by itself to reconstruct the exact test. If the full paper becomes available,
   add a separate implementation and tests.
6. Preserve the `ss_rcif()` `cause_weighting` sensitivity option. The printed
   Equation 15 in the 2023 paper is potentially non-additive under literal use.
7. Run `roxygen2::roxygenise()` and `devtools::check()` after changes whenever R
   is available.

Public functions should retain roxygen2 documentation and references.
