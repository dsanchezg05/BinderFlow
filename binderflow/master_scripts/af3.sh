#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition AF3
#SBATCH --exclusive
#SBATCH --gres=gpu:1

#conda activate base
script_dir=$1
env_dir=$2

# collect the data from hits and make .json1 
function json_1 {
    sc_dir=$script_dir
    envd=$env_dir

    env sc_dir="$sc_dir" envd="$envd" python3 - <<'END'
import sys
import os
from os import environ
from pathlib import Path
import warnings
sys.path.append(f"{environ.get('sc_dir')}/binderflow/scripts")
import af3_json_generator

warnings.filterwarnings("ignore")

pdb_list = []
loc = environ.get("envd")

hit_dir = os.listdir(f"{loc}/hits/")
if "pdbs" in hit_dir:
    pdb_path = os.path.abspath(f"{loc}/hits/pdbs/")
    for i in Path(pdb_path).iterdir():
        if str(i).endswith(".pdb"):
            print(i)
            pdb_list.append(i)

af3_json_generator.json_1(pdb_list, loc)
END
}

json_1

#execute af3_cryo two times, one for json1 to collect data adn second for json2 (loop) to make prediction for every project

#json1


template_0_ref=$(find $env_dir/AF3_run/AF3_json/Structure* -name "Template*" | head -1)

if [[ ! -f "$(find $env_dir/AF3_run/AF3_results/msa_* -name "*.json")" ]]; then
    echo "Running AF3 template data on $(basename $template_0_ref)"
    bash $script_dir/binderflow/scripts/submit_af3_cryo.sh --input "$template_0_ref" --data True --directory "$script_dir" --local_dir "$env_dir"
else
    template_msa=$(find $env_dir/AF3_run/AF3_results/msa_* -name "*.json")
    echo
    echo "MSA already performed in $(basename $template_msa), skipping..."
    echo

fi


wait


source $script_dir/config.sh
conda activate $AF3_env
#make a3m files
msa_dir=$(find $env_dir/AF3_run/AF3_results/msa_* -name "msa_*.json")
python3 $script_dir/binderflow/scripts/AF3_get_msa.py --input_json_path "$msa_dir" --output_dir "$env_dir/AF3_run/AF3_json/"

wait


function name {
    n_dir=$pdb_file python3 - <<END
import sys
import os
from os import environ
import json
with open('{}'.format(environ.get('n_dir'))) as f:
    pdb_dic= json.load(f)
pdb_name=pdb_dic["name"]
print(pdb_name)
END
}
#json2

project_counter=0;total_name="ini"
pdbs_inference=$(find $env_dir/AF3_run/AF3_json/Structure_* -name "Inference*")
for pdb_file in $pdbs_inference; do
    pdb_name=$(name)
    pdb_basename=$(basename $pdb_file)
    if [[ ! -d "$env_dir/AF3_run/AF3_results/$pdb_name" ]]; then
        total_name=$total_name/"$pdb_name"
        echo "Running AF3 inference on $pdb_name"
        bash $script_dir/binderflow/scripts/submit_af3_cryo.sh --input "$pdb_file" --directory "$script_dir" --project "$project_counter" --local_dir "$env_dir"
        wait
        let project_counter=$project_counter+1
    else
        echo
        echo "Inference already performed on $pdb_name, skipping..."
        echo
    fi

done



#organise results and update stats.csv from af2ig

source $script_dir/config.sh
conda activate $BINDERFLOW_ENV
echo "$total_name"
python3 $script_dir/binderflow/scripts/AF3_summary_info.py --directory "$script_dir" --names "$total_name" --local_dir "$env_dir"
wait

echo 
echo
echo "Done"
echo "$project_counter PDBs infered succesfully"
if [[ $project_counter -gt 0 ]]; then
    echo "PDBs aviable in $env_dir/AF3_run/AF3_results/run_name/seed_"
fi