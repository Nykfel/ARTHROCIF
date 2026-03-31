## Screening of predicted cifA and cifB protein domains

The domain characterization of all sequences of interest was performed using the HHpred platform, with the following reference databases: SCOPe70_2.08, Pfam-A_v37.4, COG_KOG_v1.0, and SMART_v6.

These sequences were then compared using BLASTp against a custom database composed of isolated domain sequences from cifB genes reported by Amoros et al. (2025). This approach allowed the detection of domains or alignments that may have been missed by HHpred, particularly short or weakly conserved regions.

This workflow was complemented by a manual alignment step using the UGENE software, in order to distinguish true homologous regions from potential false positives.

## Data structure and processing

Information related to the presence of domains for each sequence is provided in the table Domains_cifAB-like.xlsx.

In a broad sens, columns labeled "STOP" indicate an ORF disruption event.

The script plot_domain_architecture.R used to visualize proteins and their associated domains requires the data to be in a long format rather than a wide format. Therefore, a data transformation step was necessary, which was performed using the script pivot_domains_long.R.
