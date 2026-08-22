#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

CONFIG="config/01a_clustree/vbct/small_recluster_qc_pass_only.json"
if [[ $# -gt 0 && "${1}" != --* ]]; then
  CONFIG="${1}"
  shift
fi

python3 01a_run_clustree_from_config.py --config "${CONFIG}" "$@"
