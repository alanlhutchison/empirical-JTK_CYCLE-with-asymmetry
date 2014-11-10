Hi!

This is the read-me for running the analysis described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2014 (in review).

The commands should be run in this order:

jtk7.py
make_permutations.sh
get_empP_do_BH.py


use "head script_name" or "./script_name -h" to see usage.


It is important that the same time series file and period/phase/asymmetry parameters (determined by the files found in ref_files) that are used in jtk7.py are used in make_permutations.sh. 

It is important that the time series data you are analyzing have the same period/phase/asymmetry search parameters and the same time point distribution (header) as the permutation file that you are referencing. For example, if you analyzed data under one condition that had the header:

# ZT0 ZT4 ZT8 ZT12 ZT16 ZT20 ZT0 ZT4 ZT8 ZT12 ZT16 ZT20

then you could not use those permutations on this header:

# ZT0 ZT4 ZT8 ZT12 ZT16 ZT20 ZT0 ZT4 ZT8 ZT16 ZT20

instead, you would need to recalculate the permutations.










