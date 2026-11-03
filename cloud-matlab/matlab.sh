#!/usr/bin/env bash
set -euo pipefail

matlab_release="${MATLAB_RELEASE:-R2026a}"
matlab_root="${MATLAB_ROOT:-/workspace/MATLAB/$matlab_release}"
matlab_license_mode="${MATLAB_LICENSE_MODE:-online}"
case "$matlab_license_mode" in
    network|online|onlinelicensing) ;;
    *) printf 'MATLAB_LICENSE_MODE must be network, online or onlinelicensing.\n' >&2; exit 64 ;;
esac
matlab_system_libs="${MATLAB_SYSTEM_LIBS:-/workspace/MATLAB/system-libs}"
matlab_bootstrap_work="${MATLAB_BOOTSTRAP_WORK:-/workspace/work/matlab-bootstrap}"
matlab_license_endpoint_file="${MATLAB_LICENSE_ENDPOINT_FILE:-$matlab_bootstrap_work/license-server.txt}"
if [[ "$matlab_license_mode" == network && -z "${MLM_LICENSE_FILE:-}" && -r "$matlab_license_endpoint_file" ]]; then
    IFS= read -r MLM_LICENSE_FILE < "$matlab_license_endpoint_file" || [[ -n "${MLM_LICENSE_FILE:-}" ]]
    export MLM_LICENSE_FILE
fi
if [[ ! -x "$matlab_root/bin/matlab" ]]; then
    printf 'MATLAB is not installed at %s. Run cloud-matlab/bootstrap.sh first.\n' "$matlab_root" >&2
    exit 69
fi
if [[ "$matlab_license_mode" == network && -z "${MLM_LICENSE_FILE:-}" ]]; then
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
for executable_dir in "$matlab_system_libs/usr/bin" "$matlab_system_libs/usr/sbin" "$matlab_system_libs/bin" "$matlab_system_libs/sbin"; do
    if [[ -d "$executable_dir" ]]; then
        export PATH="$executable_dir:$PATH"
    fi
done
if [[ "$matlab_license_mode" != network ]]; then
    # Select the documented account mode explicitly; do not reuse a network endpoint.
    exec env -u MLM_LICENSE_FILE "$matlab_root/bin/matlab" -licmode "$matlab_license_mode" "$@"
fi
exec "$matlab_root/bin/matlab" "$@"
