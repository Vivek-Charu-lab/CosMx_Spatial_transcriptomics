#!/bin/bash

# RStudio Server launcher via Apptainer
# Image: rocker/rstudio:4.6.0

# Create Job session
# ==================
# srun --partition=gpu --time=08:00:00 --cpus-per-task=8 --mem=128G --pty --gres=gpu:1 bash 
source ~/.bashrc
micromamba activate /oak/stanford/groups/bhowitt/conda_envs/workspace


# Set Remote connection params
# ============================
SIF=/scratch/users/franzake/rstudio_server/rstudio_4.6.0.sif
PORT=8787
RSTUDIO_TMP="${SCRATCH}/rstudio_tmp"

mkdir -p "${RSTUDIO_TMP}/run" "${RSTUDIO_TMP}/var-lib-rstudio-server" "${RSTUDIO_TMP}/var-run"

# - Write rserver config
cat > "${RSTUDIO_TMP}/rserver.conf" <<EOF
www-port=${PORT}
server-daemonize=0
server-data-dir=${RSTUDIO_TMP}/var-run
EOF

# - Use a simple password (set PASSWORD env var to override, default: rstudio)
export PASSWORD="${PASSWORD:-rstudio}"

echo "Starting RStudio Server on port ${PORT}..."
echo "Connect via: http://localhost:${PORT}"
echo "Username: $(whoami) | Password: ${PASSWORD}"
echo ""

CONDA_ENV=/oak/stanford/groups/bhowitt/conda_envs/workspace
apptainer exec \
    --env PASSWORD="${PASSWORD}" \
    --env RSTUDIO_WHICH_R="${CONDA_ENV}/bin/R" \
    --env LD_LIBRARY_PATH="${CONDA_ENV}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    --bind "${RSTUDIO_TMP}/run:/run" \
    --bind "${RSTUDIO_TMP}/var-lib-rstudio-server:/var/lib/rstudio-server" \
    --bind "${RSTUDIO_TMP}/rserver.conf:/etc/rstudio/rserver.conf" \
    --bind /share/software:/share/software \
    "${SIF}" \
    rserver --server-pid-file="${RSTUDIO_TMP}/rserver.pid" \
            --server-user="$(whoami)" \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --www-port=${PORT}
