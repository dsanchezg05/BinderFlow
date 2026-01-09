#!/bin/bash


unified="False"
data="None"
# Add all things needed to the path
# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input) input="$2" ; shift ;;
        -p|--project) project="$2"; shift;;
        -u|--unified_memory) unified="$2" ; shift ;;
        -i|--directory) directory="$2"; shift;;
        -l|--local_dir) local_dir="$2"; shift;;
        -d|--data) data="$2" ; shift ;; # Run only in data mode, no inference. When running in data mode, the JSON file should not contain anything else than the sequence (no MSApath)
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift # Shift past the current argument
done

source /apps/profile.d/load_all.sh
#echo $directory
source $directory/config.sh
conda activate $AF3_env
#conda activate AF3
export PATH="/apps/alphafold3/hmmer/bin:/apps/miniconda3/envs/AlphaFold3:$PATH"

# Add Linear Algebra Accelerator mods
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

if [ $unified = "False" ]; then 
    export XLA_PYTHON_CLIENT_PREALLOCATE=true
    export XLA_CLIENT_MEM_FRACTION=0.95
else
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2
fi
if [ "$data" == "None" ]; then
    echo ""Runnning in inference mode with input: $input"" 
    file_basename=$(basename $input)
    python3 $AF3_PATH/run_alphafold.py --json_path "$input" --output_dir $local_dir/AF3_run/AF3_results/ --db_dir /emdata_fast/cryoemadmin/models/AF3/ --model_dir /emdata_fast/cryoemadmin/models/AF3/ 
fi
if [ "$data" == "True" ]; then
    echo "Running in data mode, no inference will be performed."
    python3 $AF3_PATH/run_alphafold.py --json_path "$input" --output_dir $local_dir/AF3_run/AF3_results/ --db_dir /emdata_fast/cryoemadmin/models/AF3/ --model_dir /emdata_fast/cryoemadmin/models/AF3/ --norun_inference
fi

