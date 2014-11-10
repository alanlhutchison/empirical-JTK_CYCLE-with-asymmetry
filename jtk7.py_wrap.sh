#!/bin/bash
#
# Use -h to see the usage
#

usage ()
{
echo "  -f ) Filename of List of series filenames"
echo "  -w ) Waveform_filename"
echo "	-o ) Output filename (optional)"
echo "	-p ) Period filename"
echo "	-s ) Phase filename"
echo "	-a ) Asymmetry filename"
echo "  -x ) Prefix"
echo "  -h ) usage"
}

FPATH=$(pwd)

prefix=""
fn_list=""
fn_waveform=""
fn_period=""
fn_phase=""
fn_width=""
fn_output=""

while getopts "f:x:w:p:s:a:oh" opt; do
    case $opt in
        f ) fn_list=$OPTARG;;
	x ) prefix=$OPTARG;;
        w ) fn_waveform=$OPTARG;;
	o ) fn_output=$OPTARG;;
	p ) fn_period=$OPTARG;;
        s ) fn_phase=$OPTARG;;
	a ) fn_width=$OPTARG;;
        h ) usage
        exit 0;;
        *) usage
        exit 1;;
    esac
done


list=$(cat ${FPATH}/${fn_list})
for fn in $list
do
    ffn=${fn%%.txt}
    fn_output=${ffn}_jtkout.txt

    #sbatch jtk7.py_sbatch -f ${FPATH}/$fn -x $prefix -w ${FPATH}/$fn_waveform -o ${FPATH}/$fn_output -p ${FPATH}/$fn_period -s ${FPATH}/$fn_phase -a ${FPATH}/$fn_width



    # Strip widths, width, and .txt off of the width file name
    w_prefix=${fn_width#widths_}
    ww_prefix=${w_prefix#width_}
    www_prefix=${ww_prefix%.txt}


    # Strip phases and .txt off of the width file name
    ph_prefix=${fn_phase#phases_}
    pph_prefix=${ph_prefix%.txt}


    # Strip period and .txt off of the width file name
    er_prefix=${fn_period#periods_}
    per_prefix=${er_prefix#period_}
    pper_prefix=${per_prefix#period}
    ppper_prefix=${pper_prefix%.txt}


    prefix="per${ppper_prefix}_ph${pph_prefix}_a${www_prefix}"
    echo $prefix
    echo "/project/dinner/alanlhutchison/Fly/run/JTK/jtk7.py_sbatch -f $fn -x $prefix -w $fn_waveform -p $fn_period -s $fn_phase -a $fn_width"
    sbatch_i=$(sbatch /project/dinner/alanlhutchison/Fly/run/JTK/jtk7.py_sbatch -f $fn -x $prefix -w $fn_waveform -p $fn_period -s $fn_phase -a $fn_width)
    sbatch_id=$(echo $sbatch_i | awk '{print $4}')

    /project/dinner/alanlhutchison/Fly/run/JTK/sample_data/make_permutations/make_permutations_off_file_jtk_permute_3_extract_best.py_sbatch_wrap.sh ${FPATH} $fn $fn_waveform $fn_period $fn_phase $fn_width $sbatch_id

    

done
