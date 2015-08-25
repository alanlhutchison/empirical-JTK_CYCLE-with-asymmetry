#!/bin/bash
# Created by Alan L. Hutchison, alanlhutchison@uchicago.edu, 
# Aaron Dinner Group, University of Chicago
# Use the -h flag to see the usage
#

usage ()
{
echo "  -f ) Filename of List of series filenames"
echo "  -w ) Waveform_filename"
echo "	-p ) Period filename"
echo "	-s ) Phase filename"
echo "	-a ) Asymmetry filename"
echo "  -x ) Prefix"
echo "  -t ) paTh (directory from where script is run)"
echo "  -h ) usage"
}

FPATH=""

prefix=""
fn=""
fn_waveform=""
fn_period=""
fn_phase=""
fn_width=""
fn_output=""

while getopts "f:x:w:p:s:a:t:oh" opt; do
    case $opt in
        f ) fn=$OPTARG;;
	x ) prefix=$OPTARG;;
        w ) fn_waveform=$OPTARG;;
	p ) fn_period=$OPTARG;;
        s ) fn_phase=$OPTARG;;
	a ) fn_width=$OPTARG;;
	t ) FPATH=$OPTARG;;
        h ) usage
        exit 0;;
        *) usage
        exit 1;;
    esac
done




#for num in 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
for num in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    fn_out=${prefix}_permutation_${num}.txt
    fn_out_2=${fn_out%.txt}_jtkout.txt
    make_permutations/make_permutations_off_file_jtk_permute_3_extract_best_wrap.sh ${FPATH} $fn ${fn_out} $fn_waveform $fn_period $fn_phase $fn_width $fn_out_2 >> output_file
done


cat ${prefix}_permutation_*_jtkout.best_ps.txt > ${prefix}_permutation_01-20_jtkout.best_ps.txt

