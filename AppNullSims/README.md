## Null simulations for the application

This directory contains code (and results) to do simulations under the
null hypothesis of no QTL, for the B6&times;BTBR application.

The main script [`1_get_threshold_batch.R`](1_get_threshold_batch.R)
was run in 96 batches, using `SUB` replaced successively by the values
1, 2, ..., 96.

[`2_combine_results.R`](2_combine_results.R) combines the results into
one object, and [`3_calc_penalties.R`](3_calc_penalties.R) turns those
results into estimated penalties.
