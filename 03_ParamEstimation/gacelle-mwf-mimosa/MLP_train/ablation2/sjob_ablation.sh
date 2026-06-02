#!/bin/bash
#SBATCH --job-name sjob_ablation
#SBATCH --account=inverse
#SBATCH --partition=rtx6000,rtx8000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=128G
#SBATCH --time=6-12:00:00
#SBATCH --gpus=1
#SBATCH --array=1-42%20
#SBATCH --output=slurm-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

cd /autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation2

CFG_FILE=ablation_list.txt
CFG_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CFG_FILE})

echo "Running ablation: ${CFG_NAME}"
echo "JobID: ${SLURM_JOB_ID}  TaskID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"

# English comment: you must set MATLAB Runtime root for compiled app.
# Option A: if your cluster has a module for MCR, load it here (recommended).
# module load mcr/<version>

# Option B: set MCRROOT manually.
# Example (placeholder): MCRROOT=/path/to/MATLAB_Runtime/vXX
export MCRROOT=/autofs/cluster/matlab/24.2
./run_run_one_ablation.sh "${MCRROOT}" "${CFG_NAME}"

if [ -f "./run_run_one_ablation.sh" ]; then
  if [ -z "${MCRROOT}" ]; then
    echo "ERROR: MCRROOT is not set, but run_run_one_ablation.sh exists."
    echo "Set MCRROOT to your MATLAB Runtime root (e.g., export MCRROOT=/path/to/MATLAB_Runtime/vXX)."
    exit 1
  fi
  ./run_run_one_ablation.sh "${MCRROOT}" "${CFG_NAME}"
elif [ -x "./run_one_ablation" ]; then
  # English comment: some builds can run directly if runtime is discoverable.
  ./run_one_ablation "${CFG_NAME}"
else
  echo "ERROR: Cannot find compiled launcher (run_run_one_ablation.sh) or executable (run_one_ablation)."
  ls -lah
  exit 1
fi