# Hummingbirds
Citation: Robert K. Colwell, Thiago F. Rangel, Karolina Fučíková, Diego Sustaita, Gregor M.
Yanega, and Alejandro Rico-Guevara: Repeated evolution of unorthodox feeding styles drives a negative correlation between foot size and bill length in hummingbirds. American Naturalist.

Corresponding authors; email: colwell@uconn.edu; thiago.rangel@ufg.br; colibri@uw.edu
Author of this script: karolina.fucikova@gmail.com
ORCIDs: Colwell, https://https://orcid.org/0000-0002-1384-0354; Rangel, http://orcid.org/ 0000-
0002-2001-7382; Fucíková, http://orcid.org/0000-0002-2177-4692; Sustaita,
https://orcid.org/0000-0001-9932-909X; Yanega, https://orcid.org/0000-0001-7896-7738; Rico Guevara https://orcid.org/0000-0003-4067-5312

Abstract:
Differences among hummingbird species in bill length and shape have rightly been viewed as adaptive, in relation to the morphology of the flowers they visit for nectar. In this study we examine functional variation in a behaviorally related, but neglected feature: hummingbird feet. We gathered records of hummingbirds clinging by their feet to feed legitimately as pollinators or illegitimately as nectar-robbers—“unorthodox” feeding behaviors. We measured key features of bills and feet for 220 species of hummingbirds and compared the 66 known “clinger” species (covering virtually the entire scope of hummingbird body size) to the 144 presumed “non-clinger” species. Once the effects of body size, elevation, and phylogenetic signal are accounted for statistically, hummingbirds display a surprising, but functionally interpretable negative correlation. Short-billed species have evolved exceptionally long hallux (hind-toe) claws—independently, more than 20 times, and in every major clade. Their biomechanically enhanced feet allow them to save energy by clinging to feed legitimately on short-corolla flowers and by stealing nectar from long-corolla flowers, especially at high elevations, where hovering is exceptionally costly. In contrast, long-billed species have smaller hallux claws, as plant species with long-corolla flowers enforce hovering to feed, simply by the way they present their flowers.

This R project contains files used to calculate Pagel's lambda and infer the evolution of feeding styles using the hummingbird phylogeny.

src

The R scripts are the src directory - Pagel_Rscript.R for Pagel's lambda calculations and Character_evolution.R for character evolution mapping analyses. Each script is annotated to indicate how the analyses were done and how specific lines can be modified to adapt to different data. The annotations also include links to relevant tutorials and blogs.

data

The underlying data are located in the data directory and contain a Bayesian consensus tree in .tre format (HumtreeCons.tre) as well as a csv table (HumData2.csv) of measurements of body parts (bill and foot, and the relevant standard error - SE - values) and feeding styles in different hummingbird species. Feeding styles are recorded in two ways in the table - FeedStyle column contains the specific feeding styles 1-7 (e.g., feeds legitimately on the wing, feeds through openings while clinging...), whereas the Clinging column contains the same data coded as binary: clinging (1) vs. non-clinging (0). The Clinging column was used to map the evolution of feeding on the phylogeny.

The tree and table do not correspond 100% - names match but there are taxa that are in the tree but not in the table (or do not have data for all variables in the table), and vice versa. Missing data are dropped from the individual analyses in R - the first part of each script handles the 'cleanup' of the data, then each script proceeds with analyses.

figures

The figures folder includes simmap results under the ER model (the better-fitting model selected based on AIC) as well as under the more complex ARD model (for illustration only). The results are plotted in two ways for each model: tree with pie charts at nodes showing odds of a particular state (feeding style) being present in an ancestor, and tree with color gradients on branches representing the density of character state inference across the 100 simulation performed.