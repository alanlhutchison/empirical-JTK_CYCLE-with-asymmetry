Hi!

This is the read-me for running the analysis described in Hutchison AL, Maienscein-Cline M, Chiang AH, Tabei SMA, Gudjonson H, Bahroos N, Allada R, Dinner AR. “Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data.” PLoS Computational Biology 2015 Mar. Vol. 11, No. 3, pp. e1004094, DOI: 10.1371/journal.pcbi.1004094

An application of the method to 16S gut microbiome data can be found at: Leone VA, Gibbons SM, Martinez K, Hutchison AL, Huang EY, Cham CM, Pierre JF, Heneghan AF, Nadimpalli A, Hubert N, Zale E, Wang Y, Huang Y, Theriault B, Dinner AR, Musch MW, Kudsk KA, Prendergast BJ, Gilbert JA, Chang EB. “Effects of diurnal variation of gut microbes and high fat feeding on host circadian clock function and metabolism” Cell Host-Microbe 2015. doi:10.1016/j.chom.2015.03.006 (in press)

It is based on the original JTK_CYCLE code from Hughes ME, Hogenesch JB, Kornacker K. JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets. J Biol Rhythms. 2010 Oct;25(5):372-80. doi: 10.1177/0748730410379711. 

This implementation is in Python and shell script.

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

Currently, making permutations is the time-limiting step of the empirical calculation. To facilitate the use of empirical JTK, we have provided the null distributions for some time point sets suggested by users. If you are interested the null distribution for a particular set of time points, please contact the first author at alanlhutchison ...at... uchicago ...dot... edu. Please note that the headers are labeled with ZT, but will work equally well for CT time points. 

This code is released with the MIT License. See the License.txt file for more information.
