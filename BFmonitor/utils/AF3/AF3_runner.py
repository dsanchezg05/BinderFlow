'''
Script to run AF3 crossvalidation for the selected binders.

THINGS TO DO:
- Add the option to get sequences from the silent file, so you do not have to extract the PDBs before running AF3
- 
'''
import sys
sys.stdout.reconfigure(line_buffering=True) # to update af3 logfile while dash is open

import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import subprocess
import argparse
import glob
import time
SCRIPT_DIR = os.path.dirname(__file__)
sys.path.append(f'{SCRIPT_DIR}/BFmonitor/utils/AF3/')
from af3_utils import *
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--hits_dir", "-p", help = "Hits directory")
parser.add_argument("--directory", "-d", help='Working directory')
args = parser.parse_args()

hits_dir = args.hits_dir
dir = args.directory



# load variables from config file
## Get script directory
script_dir = Path(__file__).parent
print(script_dir)
config_file = str(script_dir) + '/../../../config.sh'
print(config_file)
variables = load_bash_vars(config_file)

# Define the template directory
template_dir = f"{dir}/AF3_run/AF3_json/Target_data/"
output_dir = f"{dir}/AF3_run/AF3_output/"

# Generate the JSON file for data
AF3_target_json_generator(hits_dir, dir)

# Run AF3 in data_mode
## First, check if AF3 has already been run in data mode
if AF3_check_done(output_dir, num_inferences=-1):
    print("AF3 data mode has already been run. Skipping this step.")
else:  
    print('running AF3 in data mode...')
    print(variables)
    command = f'sbatch --nodes="{variables["NODES_AF3"]}" -p "{variables["PARTITION_AF3"]}" --open-mode=append --gres="{variables["GRES_AF3"]}" --exclusive --cpus-per-gpu="{variables["CPUS_PER_GPU_AF3"]}" -o {dir}/AF3_run/slurm_logs/%j.out '
    command += f'{script_dir}/af3.sh -i {template_dir}/Template_data.json -o {output_dir} -d True -c {config_file}'
    print(command)
    subprocess.run(command, shell=True)

while not AF3_check_done(output_dir, num_inferences=-1):
    print("AF3 is still running...")
    time.sleep(60)  # Check every minute

# Generate the A3M files that you will need for the next steps
data_json_result = glob.glob(f'{output_dir}/msa_template*/*.json')[0]
data_dir = os.path.dirname(os.path.abspath(data_json_result))

command = f'conda run -n {variables["AF3_ENV"]} --no-capture-output python {script_dir}/AF3_get_msa.py --input_json_path {data_json_result} --output_dir {data_dir}'
# use `conda run -n <env>` so the command works in non-interactive shells without sourcing conda
subprocess.run(command, shell=True)

# Generate the JSON files for the input pdbs
AF3_binder_json_generator(hits_dir, data_dir,dir)


# Send the AF3 runs for the binder structures
num_inferences=0
for json_file in glob.glob(f'{dir}/AF3_run/AF3_json/Structure_*/Inference_*.json'):
    print(json_file)
    run_dir = output_dir + 'inference_' + json_file.split('/')[-1].split('Inference_')[-1].split('.')[0]
    print(run_dir)
    num_inferences=num_inferences+1
    if os.path.exists(run_dir):
        print(f"AF3 has already been run for {json_file}. Skipping this job submission.")
        continue
    else:
        print("Submitting AF3 job for:", json_file)
        command = f'sbatch --nodes="{variables["NODES_AF3"]}" -p "{variables["PARTITION_AF3"]}" --open-mode=append --gres="{variables["GRES_AF3"]}" --exclusive --cpus-per-gpu="{variables["CPUS_PER_GPU_AF3"]}" -o {dir}/AF3_run/slurm_logs/%j.out'
        command += f' {script_dir}/af3.sh -i {json_file} -o {output_dir} -c {config_file}'
        subprocess.run(command, shell=True)
        
print(f"All AF3 jobs for binder structures have been submitted, {num_inferences}")

while not AF3_check_done(output_dir, num_inferences, data = False):
    print("AF3 is still running...")
    time.sleep(60)  # Check every minute

# Once all AF3 jobs are done, extract the summary information
print("Running Summary info")
AF3_summary_info(dir)


time.sleep(30)
print("Loading RMSD script")
SCRIPT_DIR= os.path.dirname(__file__)
path_rmsd=SCRIPT_DIR.split("/")[:-3]
p="/".join(path_rmsd)
rmsd=f"{p}/binderflow/scripts"

path=Path(hits_dir)
hits=path.parent.absolute()
#print(rmsd)
#print(hits)

if num_inferences != 0:
    print(num_inferences)
    clstr=int(round(np.log10(num_inferences)*5,1))
    print(clstr)
    subprocess.run([f"python3 {rmsd}/RMSD_rosetta.py --directory {hits} --t {clstr}"],shell=True)
else:
    print("RMSD script actually run")
    pass
