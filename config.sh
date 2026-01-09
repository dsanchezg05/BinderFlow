#Manually set the paths for every other software you have installed and environment name

RFD_PATH="/apps/rosetta/RFDifussion"
PMPNN_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as SILENT and AF2IG
AF2IG_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and SILENT
BINDERFLOW_PATH="$(dirname "$(realpath "${BASH_SOURCE[0]}")")" #Do not change this
SILENT_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and AF2IG

# AF3 paths
AF3_PATH="/apps/alphafold3" #Path to AlphaFold3 repo
HMMER_PATH="/apps/alphafold3/hmmer/bin"
AF3_DB_DIR="/emdata_fast/cryoemadmin/models/AF3/"
AF3_MODEL_DIR="/emdata_fast/cryoemadmin/models/AF3/"
# MANUALLY SET ENVIRONMENTS NAMES (example names set)

RFD_ENV="SE3nv4090"
PMPNN_ENV="dl_binder_design"
AF2_ENV="af2_binder_design"
BINDERFLOW_ENV="watcher"
AF3_ENV="AF3"


#MANUALLY SET SBATCH CONFIGURATIONS (examples in place)

NODES=1
PARTITION=RFD
CPUS_PER_GPU=12
GRES=gpu:1

# SLURM CONFIGURATIONS FOR AF3

NODES_AF3=1
PARTITION_AF3=AF3
CPUS_PER_GPU_AF3=12
GRES_AF3=gpu:1

source /apps/profile.d/load_all.sh
