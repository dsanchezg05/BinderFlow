#!/bin/bash

# Get the input arguments
unified="False"
data="None"
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input) input="$2" ; shift ;;
        -c|--config) config="$2" ; shift ;;
        -o|--output) output="$2"; shift;;
        -u|--unified_memory) unified="$2" ; shift ;;
        -d|--data) data="$2" ; shift ;; # Run only in data mode, no inference. When running in data mode, the JSON file should not contain anything else than the sequence (no MSApath)
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift # Shift past the current argument
done

# Load the config file
source $config

# Run AF3 

conda activate $AF3_ENV
#conda activate AF3
export PATH="$HMMER_PATH:$AF3_PATH:$PATH"

# Add Linear Algebra Accelerator mods
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

if [ "$unified" = "False" ]; then
    export XLA_PYTHON_CLIENT_PREALLOCATE=true
    export XLA_CLIENT_MEM_FRACTION=0.95
else
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2
fi
if [ "$data" == "None" ]; then
    echo ""Runnning in inference mode with input: $input"" 
    python3 $AF3_PATH/run_alphafold.py --json_path "$input" --output_dir $output --db_dir $AF3_DB_DIR --model_dir $AF3_MODEL_DIR
fi
if [ "$data" == "True" ]; then
    echo "Running in data mode, no inference will be performed."
    python3 $AF3_PATH/run_alphafold.py --json_path "$input" --output_dir $output --db_dir $AF3_DB_DIR --model_dir $AF3_MODEL_DIR --norun_inference
fi

