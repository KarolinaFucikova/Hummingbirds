# Hummingbirds
This R project contains files used to calculate Pagel's lambda using the hummingbird phylogeny.
The R script in the src directory is annotated to indicate how the analyses were done and how specific lines can be modified to adapt to different data.
The underlying data are located in the data directory and contain a Bayesian consensus tree in .tre format as well as a csv table of measurements of body parts (bill and foot) in different hummingbird species.
The tree and table do not correspond 100% - names match but there are taxa that are in the tree but not in the table (or do not have data for all variables in the table), and vice versa. Missing data are dropped from the individual analyses in R.