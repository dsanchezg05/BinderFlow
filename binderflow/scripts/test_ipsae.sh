#!/bin/bash

script_dir=$1
dir_dir=$2
loc=$3

python3 $script_dir/ipsae.py $loc/AF3_run/AF3_results/$dir_dir/$dir_dir"_confidences.json" $loc/AF3_run/AF3_results/$dir_dir/$dir_dir"_model.cif" 15 15

wait

ipsae=$(cat $loc/AF3_run/AF3_results/$dir_dir/*_15.txt | awk '{print $6}' | tail -3 | sort -n | tail -1)

echo $ipsae