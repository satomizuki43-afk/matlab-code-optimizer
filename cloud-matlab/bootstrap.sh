#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
matlab_release="${MATLAB_RELEASE:-R2026a}"
matlab_license_mode="${MATLAB_LICENSE_MODE:-online}"
case "$matlab_license_mode" in
    network|online|onlinelicensing) ;;
    *) printf 'MATLAB_LICENSE_MODE must be network, online or onlinelicensing.\n' >&2; exit 64 ;;
esac
export MATLAB_ROOT="${MATLAB_ROOT:-/workspace/MATLAB/$matlab_release}"
matlab_bootstrap_work="${MATLAB_BOOTSTRAP_WORK:-/workspace/work/matlab-bootstrap}"
matlab_tmp_base="${MATLAB_INSTALL_TMP_BASE:-/tmp}"
matlab_license_endpoint_file="${MATLAB_LICENSE_ENDPOINT_FILE:-$matlab_bootstrap_work/license-server.txt}"
if [[ "$matlab_license_mode" == network && -z "${MLM_LICENSE_FILE:-}" && -r "$matlab_license_endpoint_file" ]]; then
    IFS= read -r MLM_LICENSE_FILE < "$matlab_license_endpoint_file" || [[ -n "${MLM_LICENSE_FILE:-}" ]]
    export MLM_LICENSE_FILE
fi
check_only=false
install_only=false
if [[ "${1:-}" == --check && $# == 1 ]]; then
    check_only=true
elif [[ "${1:-}" == --install-only && $# == 1 ]]; then
    install_only=true
elif [[ $# != 0 ]]; then
    printf 'Usage: bash cloud-matlab/bootstrap.sh [--check|--install-only]\n' >&2
    exit 64
fi

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    printf 'This bootstrap requires Linux x86_64.\n' >&2
    exit 69
fi
# This dependency snapshot targets Debian 13, as used by this cloud environment.
. /etc/os-release
if [[ "$ID" != debian || "${VERSION_ID%%.*}" != 13 || "$matlab_release" != R2026a ]]; then
    printf 'This bootstrap targets Debian 13 and MATLAB R2026a. Update the OS/dependency plan before choosing another release.\n' >&2
    exit 69
fi

printf 'Target: MATLAB %s, Linux x86_64, Debian %s\n' "$matlab_release" "$VERSION_ID"
printf 'Install directory: %s\n' "$MATLAB_ROOT"
printf 'License mode: %s (authorization not yet verified).\n' "$matlab_license_mode"
df -h /workspace "$matlab_tmp_base"
if [[ "$matlab_license_mode" != network ]]; then
    printf 'Account sign-in is required on this cloud machine; saved network endpoints are ignored.\n'
elif [[ -n "${MLM_LICENSE_FILE:-}" ]]; then
    printf 'Network license endpoint configured: yes (value hidden).\n'
else
    printf 'Network license endpoint configured: no.\n'
fi

if "$check_only"; then
    printf 'Read-only checks; network and license checkout have not been verified.\n'
    bash "$script_dir/runtime-deps.sh" --check
    exit 0
fi

bash "$script_dir/runtime-deps.sh"
if [[ ! -x "$MATLAB_ROOT/bin/matlab" ]]; then
    workspace_available="$(df -PB1 /workspace | awk 'NR==2 {print $4}')"
    temporary_available="$(df -PB1 "$matlab_tmp_base" | awk 'NR==2 {print $4}')"
    if (( workspace_available < 6*1024*1024*1024 || temporary_available < 3*1024*1024*1024 )); then
        printf 'Insufficient install headroom: need at least 6 GiB in /workspace and 3 GiB temporary space. Installation peak usage may require more.\n' >&2
        exit 73
    fi
    mkdir -p "$matlab_bootstrap_work" "$(dirname -- "$MATLAB_ROOT")"
    # No credentials are needed to download/install; MATLAB checks the license at runtime.
    curl --fail --location --proto '=https' --tlsv1.2 \
        --retry 2 --connect-timeout 20 --max-time 300 \
        https://www.mathworks.com/mpm/glnxa64/mpm \
        --output "$matlab_bootstrap_work/mpm.part"
    chmod 700 "$matlab_bootstrap_work/mpm.part"
    mv -- "$matlab_bootstrap_work/mpm.part" "$matlab_bootstrap_work/mpm"

    matlab_install_tmp="$(mktemp -d "$matlab_tmp_base/matlab-mpm.XXXXXX")"
    trap 'rm -rf -- "$matlab_install_tmp"' EXIT
    TMPDIR="$matlab_install_tmp" "$matlab_bootstrap_work/mpm" install \
        --release="$matlab_release" \
        --destination="$MATLAB_ROOT" \
        --products MATLAB 2>&1 | tee "$matlab_bootstrap_work/install.log"
fi

if "$install_only" || [[ "$matlab_license_mode" != network ]]; then
    printf 'MATLAB files installed. Runtime and license checkout have not been verified.\n'
    if [[ "$matlab_license_mode" != network ]]; then
        printf 'Next: authenticate interactively with MATLAB_LICENSE_MODE=%s bash cloud-matlab/matlab.sh -nodisplay\n' "$matlab_license_mode"
    fi
    exit 0
fi
if [[ -z "${MLM_LICENSE_FILE:-}" ]]; then
    printf 'MATLAB files installed. Configure MLM_LICENSE_FILE before the runtime/license check.\n' >&2
    exit 78
fi
bash "$script_dir/matlab.sh" -batch "disp(version); assert(license('test','MATLAB')); disp('MATLAB startup and license check passed.');"
