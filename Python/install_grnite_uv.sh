#!/usr/bin/env bash
# =============================================================================
# GRNITE — install from GitHub using uv + CUDA 12.1 PyTorch wheels
# =============================================================================
# Upstream repository: https://github.com/aliaaz99/GRNITE
#
#
# Prerequisites: git, uv (https://docs.astral.sh/uv/), NVIDIA driver compatible
# with CUDA 12.1 runtime shipped by the PyTorch wheels.
#
# Usage (from anywhere):
#   bash /path/to/GastroDeconv2FateMap/scripts/install_grnite_uv.sh
#
# Optional environment variables:
#   GRNITE_DIR   Target directory for the clone (default: ~/src/GRNITE)
#   UV_PYTHON    Python version for the venv (default: 3.11)
#   SKIP_CLONE   If set to 1, skip git clone/pull and only install into
#                GRNITE_DIR (must already contain Main.py).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_URL="${REPO_URL:-https://github.com/aliaaz99/GRNITE.git}"
GRNITE_DIR="${GRNITE_DIR:-"${HOME}/src/GRNITE"}"
UV_PYTHON="${UV_PYTHON:-3.11}"
SKIP_CLONE="${SKIP_CLONE:-0}"

TORCH_INDEX_URL="https://download.pytorch.org/whl/cu121"
# PyG binary wheels track torch minor; 2.5.1+cu121 matches this wheel line.
PYG_FIND_LINKS="${PYG_FIND_LINKS:-https://data.pyg.org/whl/torch-2.5.0+cu121.html}"

REQUIREMENTS_FILE="${SCRIPT_DIR}/grnite_uv_requirements.txt"

if [[ ! -f "${REQUIREMENTS_FILE}" ]]; then
  echo "Missing ${REQUIREMENTS_FILE}" >&2
  exit 1
fi

if [[ "${SKIP_CLONE}" != "1" ]]; then
  if [[ ! -d "${GRNITE_DIR}/.git" ]]; then
    mkdir -p "$(dirname "${GRNITE_DIR}")"
    git clone "${REPO_URL}" "${GRNITE_DIR}"
  else
    git -C "${GRNITE_DIR}" pull --ff-only
  fi
else
  if [[ ! -f "${GRNITE_DIR}/Main.py" ]]; then
    echo "SKIP_CLONE=1 but Main.py not found in ${GRNITE_DIR}" >&2
    exit 1
  fi
fi

cd "${GRNITE_DIR}"

uv venv --python "${UV_PYTHON}" .venv
# shellcheck source=/dev/null
source .venv/bin/activate

echo "Installing PyTorch stack (CUDA 12.1 wheels)..."
uv pip install \
  --extra-index-url "${TORCH_INDEX_URL}" \
  torch==2.5.1+cu121 \
  torchvision==0.20.1+cu121 \
  torchaudio==2.5.1+cu121

echo "Installing PyG extension wheels (safe to skip failures on some platforms)..."
set +e
uv pip install \
  --extra-index-url "${TORCH_INDEX_URL}" \
  --find-links "${PYG_FIND_LINKS}" \
  pyg-lib \
  torch-scatter \
  torch-sparse \
  torch-cluster \
  torch-spline-conv
set -e

echo "Installing torch-geometric + pinned scientific stack..."
uv pip install \
  --extra-index-url "${TORCH_INDEX_URL}" \
  -r "${REQUIREMENTS_FILE}"

echo ""
echo "Done. Activate and run (example):"
echo "  source \"${GRNITE_DIR}/.venv/bin/activate\""
echo "  cd \"${GRNITE_DIR}\""
echo "  python Main.py --step 1"
echo ""
echo "Upstream docs and data layout: https://github.com/aliaaz99/GRNITE"
