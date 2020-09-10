## Null simulations for the backcross simulations

This directory contains code (and results) to do simulations under the
null hypothesis of no QTL, to be used in the backcross simulations to
assess power.

The main script [`1_null_sims_backcross_scantwo.R`](1_null_sims_backcross_scantwo.R)
was run in 96 batches, using `SUB` replaced successively by the values
1, 2, ..., 96.

[`2_combine_results.R`](2_combine_results.R) combines the results into
one object, and [`3_calc_penalties.R`](3_calc_penalties.R) turns those
results into estimated penalties.
