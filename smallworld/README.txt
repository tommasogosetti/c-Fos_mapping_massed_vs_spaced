TOPOLOGY + SMALL-WORLDNESS ANALYSIS
====================================

CREATE ONE NEW FOLDER, e.g.:

C:\Users\tomgo\Desktop\mapping\1Script paper\topology_smallworld_reviewer

PUT INSIDE IT:

MANDATORY
---------
1. topology_smallworld_permutation.py
2. ss net.csv
3. ms net.csv
4. sc net.csv
5. mc net.csv

OPTIONAL
--------
6. hc net.csv
Only if you want the script to calculate the observed HC small-world value too.
HC is NOT used in the Spaced-vs-Massed permutation test.

YOU DO NOT NEED THE GRAPHML FILES.

The script validates the four reconstructed experimental networks using the
already established edge counts and giant-component sizes:
SS 374 edges / GC 121
MS 237 edges / GC 114
SC 253 edges / GC 121
MC 430 edges / GC 113

FIRST: QUICK OBSERVED CHECK
---------------------------
Double click:
TEST_OBSERVED_ONLY.bat

This reconstructs the four networks and calculates:
- Global efficiency
- Characteristic path length
- Clustering coefficient
- Transitivity
- Small-world sigma

No permutations.

SECOND: SMOKE TEST
------------------
Double click:
TEST_100_PERMUTATIONS.bat

This runs 100 permutations and only checks that the complete pipeline works.

THIRD: FINAL ANALYSIS
---------------------
Double click:
FINAL_5000_PERMUTATIONS.bat

This runs:
- 5000 stratified subject-label permutations
- 20 connected degree-preserving random reference networks for each reconstructed network

This is computationally intensive and may take hours.

OUTPUT
------
A new folder is created automatically:

topology_smallworld_results

Most important files:
- observed_topology_metrics.csv
- network_validation.csv
- spaced_vs_massed_permutation_results.csv
- permutation_distributions.csv
- analysis_settings.json

PRIMARY STATISTIC
-----------------
For each metric:

mean[(SS - MS), (SC - MC)]

Labels are permuted independently within the spatial and cue strata.
No subject bootstrap with replacement is used.

SMALL-WORLD DEFINITION
----------------------
sigma = (C / C_rand) / (L / L_rand)

C = observed average clustering coefficient
L = observed characteristic path length
C_rand and L_rand = means of connected degree-preserving randomized networks

sigma > 1 is consistent with small-world organization relative to the
degree-preserving random reference.

DEPENDENCIES
------------
Python packages:
- numpy
- pandas
- scipy
- networkx

The BAT files install/update these packages automatically.
