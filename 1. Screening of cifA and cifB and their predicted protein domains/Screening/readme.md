The initial search was based on a [dataset](https://github.com/Nykfel/ARTHROCIF/tree/main/1.%20Screening%20of%20cifA%20and%20cifB%20and%20their%20predicted%20protein%20domains/Screening/Support/RefDomains) comprising 14 pairs of PD-(D/E)XK nuclease representing the 11 referenced motif types (3 pairs for Type X, 2 for Type IX).

A first search was conducted using the PSI-BLAST algorithm from `blastp` against the Non-redundant protein sequences database (`nr`) with two iterations. Two additional searches were performed using `tblastn` and `tblastx` against the Core nucleotide database (`core_nt`).

For each method, the search was restricted to eukaryotic organisms. For each species of interest identified, only the two BLAST hits with the lowest E-values were retained, and this was done separately for each method in order to avoid method-specific bias (filtering script: [Filtered.R](https://github.com/Nykfel/ARTHROCIF/blob/main/1.%20Screening%20of%20cifA%20and%20cifB%20and%20their%20predicted%20protein%20domains/Screening/Support/Filtered.R)).

Sequences with an E-value greater than 0.005 were considered non-significant and were therefore excluded.
