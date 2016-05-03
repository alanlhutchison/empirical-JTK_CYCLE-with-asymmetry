#Empirical JTK_CYCLE

Hi!

This is the README for running the rhythm detection method described in [Hutchison AL, Maienscein-Cline M, Chiang AH, Tabei SMA, Gudjonson H, Bahroos N, Allada R, Dinner AR. “Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data.” PLoS Computational Biology 2015 Mar. Vol. 11, No. 3, pp. e1004094](doi:10.1371/journal.pcbi.1004094)

It is based on the original JTK_CYCLE code from Hughes ME, Hogenesch JB, Kornacker K. JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets. J Biol Rhythms. 2010 Oct;25(5):372-80. doi: 10.1177/0748730410379711.

A previous version of this method empirically calculated the p-values using many permutations. We have recently sped up this calculation by approximating the null tau distribution using a Gamma distribution based on 1000 permutations. P-values can be derived from this model, allowing for p-values independent of the number of permutations. A pre-print demonstrating the accuracy of this process is forthcoming.


Applications:

An application of the method to mouse and human pancreas and liver can be found in Perelis et al. "Pancreatic β cell enhancers regulate rhythmic transcription of genes controlling insulin secretion" (2015). 350(6261). aac4250. doi:10.1126/science.aac4250

An application of the method to Drosophila neurons can be found in Flourakis M et al. "A conserved bicycle model for circadian clock control of membrane excitability" Cell (2015). 162(4). 836-848. doi:10.1016/j.cell.2015.07.036

An application of the method to 16S gut microbiome data can be found at: Leone VA et al. “Effects of diurnal variation of gut microbes and high fat feeding on host circadian clock function and metabolism” Cell Host-Microbe (2015). 17(5). 681-689. 13 May doi:10.1016/j.chom.2015.03.006

License:
This code is released with the MIT License. See the License.txt file for more information.



Running this command:

./eJTK-CalcP.py -f example/TestInput4.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT


Will produce 3 files

1) TestInput4_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout.txt
   This is the output of eJTK.py, it contains the best reference waveform matching each time series. Best is defined as having the highest Tau value. This becomes input for CalcP.py.


2) TestInput4_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtknull1000.txt
   This is the output of eJTK.py unless otherwise specified by the -n flag (see eJTK-CalcP.py -h for more information). It similar to *jtkout.txt only it contains the results of 1000 runs of Gaussian noise. It is also an input for CalcP.py


3) TestInput4_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt
   This is the output of CalcP.py. It is the equivalent of *jtkout.txt, only now with correct p-values as estimated by fitting the time series to a Gamma distribution. It also contains a column of these p-values adjusted with the Benjamini-Hochberg correction.
