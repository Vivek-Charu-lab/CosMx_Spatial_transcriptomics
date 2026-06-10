#!/bin/bash

# code-server launcher via Apptainer
# Image: linuxserver/code-server (4.123.0)
# Pull (if missing):
#   mkdir -p /scratch/users/franzake/codeserver
#   APPTAINER_CACHEDIR=/scratch/users/franzake/SINGULARITY_CACHE \
#     apptainer pull --arch amd64 /scratch/users/franzake/codeserver/code-server.sif \
#     docker://linuxserver/code-server

# Create Job session (example — adjust resources as needed)
# =========================================================
# srun --partition=gpu --time=08:00:00 --cpus-per-task=8 --mem=128G --pty --gres=gpu:1 bash
source ~/.bashrc
micromamba activate /oak/stanford/groups/bhowitt/conda_envs/workspace


# Set Remote connection params
# ============================
SIF=/scratch/users/franzake/codeserver/code-server.sif
PORT=8888
CODESERVER_TMP="${SCRATCH}/codeserver_tmp"
WORKSPACE="${SCRATCH}"   # default workspace opened in the browser

mkdir -p "${CODESERVER_TMP}/user-data" "${CODESERVER_TMP}/extensions"

# Use a simple password (set PASSWORD env var to override, default: codeserver)
export PASSWORD="${PASSWORD:-codeserver}"

echo "Starting code-server on port ${PORT}..."
echo "Connect via: http://localhost:${PORT}"
echo "Password: ${PASSWORD}"
echo ""

CONDA_ENV=/oak/stanford/groups/bhowitt/conda_envs/workspace
apptainer exec \
    --env PASSWORD="${PASSWORD}" \
    --env PATH="/app/code-server/bin:${CONDA_ENV}/bin:${PATH}" \
    --env LD_LIBRARY_PATH="${CONDA_ENV}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    --bind /share/software:/share/software \
    --bind "${SCRATCH}:${SCRATCH}" \
    --bind /oak:/oak \
    "${SIF}" \
    /app/code-server/bin/code-server \
        --bind-addr "0.0.0.0:${PORT}" \
        --auth password \
        --user-data-dir "${CODESERVER_TMP}/user-data" \
        --extensions-dir "${CODESERVER_TMP}/extensions" \
        --disable-telemetry \
        --disable-update-check \
        "${WORKSPACE}"
