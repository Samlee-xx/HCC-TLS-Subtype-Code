CN model inputs and final analysis cohort
=======================================

The original K9 CN model was fitted on 696 regions containing 1,415,572 cells.
The manuscript analyses use 667 TLSs containing 1,414,248 cells: 325 BF-TLSs and
342 TNC-TLSs.

The upstream_696 directory retains the original fitting cohort, neighborhood
coordinates and fixed-K fitting parameters. These inputs are needed by
Figure2_CN_clustering.py to reproduce the original CN model. Refitting on only
the final 667 TLSs would change the fitting dataset and could change CN labels.
The coordinate table contains only the 696 regions actually used by that fit.

cn_cell_assignment_k9_final.csv contains the frozen labels for the 696-region
fit. The neighborhood mean, observed/expected and Z-score matrices describe
that model. They are retained as the definitions of its nine CN groups.

Final analysis membership is defined by
../TLS/tls_subtype_assignment_final.csv.
Figure2_CN.py joins the frozen assignments to the final 667-TLS layer table.
The cell-type counts and downstream CN ROE tables use this final cohort.

TLS_cohort_transition_696_to_667.csv records which historical regions belong
to the final cohort and the exclusion rules for the other 29 regions.

The larger fitting cohort is retained to reproduce the model. It must not be
used as the denominator for the manuscript's 667-TLS subtype comparisons.

Figure2_CN_clustering.py uses K=9.
