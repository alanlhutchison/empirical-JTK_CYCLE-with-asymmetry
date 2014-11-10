#!/bin/bash
# Author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron Dinner Group, University of Chicago

FPATH=$1
fn=$2
fn_out=$3
fn_waveform=$4
fn_period=$5
fn_phase=$6
fn_width=$7
fn_out_2=$8



./make_permutations_off_file.py ${FPATH}/$fn ${FPATH}/$fn_out
./jtk_permute_3.py -f ${FPATH}/$fn_out --waveform ${FPATH}/$fn_waveform -p ${FPATH}/$fn_period --phase ${FPATH}/$fn_phase -w ${FPATH}/$fn_width -o ${FPATH}/$fn_out_2
./extract_best_p_values_3.py ${FPATH}/$fn_out_2
