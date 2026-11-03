#!/usr/bin/env bash
set -euo pipefail

matlab_root="${MATLAB_ROOT:-/workspace/MATLAB/R2026b}"
matlab_system_libs="${MATLAB_SYSTEM_LIBS:-/workspace/MATLAB/system-libs}"
matlab_bootstrap_work="${MATLAB_BOOTSTRAP_WORK:-/workspace/work/matlab-bootstrap}"
matlab_license_endpoint_file="${MATLAB_LICENSE_ENDPOINT_FILE:-$matlab_bootstrap_work/license-server.txt}"
if [[ -z "${MLM_LICENSE_FILE:-}" && -r "$matlab_license_endpoint_file" ]]; then
    IFS= read -r MLM_LICENSE_FILE < "$matlab_license_endpoint_file" || [[ -n "${MLM_LICENSE_FILE:-}" ]]
    export MLM_LICENSE_FILE
fi
if [[ ! -x "$matlab_root/bin/matlab" ]]; then
    printf 'MATLAB is not installed at %s. Run cloud-matlab/bootstrap.sh first.\n' "$matlab_root" >&2
    exit 69
fi
if [[ -z "${MLM_LICENSE_FILE:-}" ]]; then
    printf 'Set MLM_LICENSE_FILE to the network license endpoint (port@host) in the environment.\n' >&2
    exit 78
fi

export MATLAB_PREFDIR="${MATLAB_PREFDIR:-/workspace/work/matlab-preferences/$(basename -- "$matlab_root")}"
mkdir -p "$MATLAB_PREFDIR"

# Only add extracted libraries; never replace the inherited proxy or CA settings.
for library_dir in "$matlab_system_libs/usr/lib/x86_64-linux-gnu" "$matlab_system_libs/lib/x86_64-linux-gnu"; do
    if [[ -d "$library_dir" ]]; then
        export LD_LIBRARY_PATH="$library_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
done
exec "$matlab_root/bin/matlab" "$@"
