#!/usr/bin/env bash
set -euo pipefail

# Install missing runtime libraries into a writable prefix, without sudo.
# The package list is the official MathWorks Debian 13 / R2026b list.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
matlab_system_libs="${MATLAB_SYSTEM_LIBS:-/workspace/MATLAB/system-libs}"
matlab_bootstrap_work="${MATLAB_BOOTSTRAP_WORK:-/workspace/work/matlab-bootstrap}"
check_only=false
if [[ "${1:-}" == --check && $# == 1 ]]; then
    check_only=true
elif [[ $# != 0 ]]; then
    printf 'Usage: bash runtime-deps.sh [--check]\n' >&2
    exit 64
fi

installed() {
    [[ "$(dpkg-query -W -f='${db:Status-Status}' "$1" 2>/dev/null || true)" == installed ]]
}

missing=()
while IFS= read -r package || [[ -n "$package" ]]; do
    [[ -z "$package" || "$package" == \#* ]] && continue
    # sudo is unnecessary for a user-owned installation; locales-all includes locales.
    [[ "$package" == sudo ]] && continue
    if [[ "$package" == locales ]] && installed locales-all; then
        continue
    fi
    if ! installed "$package" && [[ ! -f "$matlab_system_libs/.packages/$package" ]]; then
        missing+=("$package")
    fi
done < "$script_dir/debian13-packages.txt"

if [[ ${#missing[@]} == 0 ]]; then
    printf 'MATLAB runtime dependencies: available.\n'
    exit 0
fi
printf 'Missing runtime packages: %s\n' "${missing[*]}"
if "$check_only"; then
    exit 2
fi

for package in "${missing[@]}"; do
    case "$package" in
        lib*|ibverbs-providers) ;;
        *) printf 'Install the missing system package in the base environment: %s\n' "$package" >&2; exit 69 ;;
    esac
done

mkdir -p "$matlab_bootstrap_work" "$matlab_system_libs/.packages"
apt_work="$matlab_bootstrap_work/apt"
mkdir -p "$apt_work/state/lists/partial" "$apt_work/cache/archives/partial" \
    "$apt_work/conf.d" "$apt_work/sources.list.d"
cat > "$apt_work/apt.conf" <<EOF
Dir::Etc::main "/dev/null";
Dir::Etc::parts "$apt_work/conf.d";
Dir::Etc::sourcelist "$apt_work/sources.list";
Dir::Etc::sourceparts "$apt_work/sources.list.d";
Dir::State "$apt_work/state";
Dir::State::status "/var/lib/dpkg/status";
Dir::Cache "$apt_work/cache";
Dir::Log "$apt_work";
APT::Sandbox::User "$(id -un)";
Acquire::Languages "none";
EOF
cat > "$apt_work/sources.list" <<'EOF'
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian trixie main
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian-security trixie-security main
EOF
# Refresh private indexes; the cloud's system apt configuration can be read-only.
http_proxy="${HTTP_PROXY:-${http_proxy:-}}" \
https_proxy="${HTTPS_PROXY:-${https_proxy:-}}" \
    APT_CONFIG="$apt_work/apt.conf" apt-get update
download_dir="$(mktemp -d "$matlab_bootstrap_work/debs.XXXXXX")"
trap 'rm -rf -- "$download_dir"' EXIT
(
    cd -- "$download_dir"
    # Keep repository signatures and TLS verification enabled.
    http_proxy="${HTTP_PROXY:-${http_proxy:-}}" \
    https_proxy="${HTTPS_PROXY:-${https_proxy:-}}" \
        APT_CONFIG="$apt_work/apt.conf" apt-get download "${missing[@]}"
)
shopt -s nullglob
archives=("$download_dir"/*.deb)
if [[ ${#archives[@]} == 0 ]]; then
    printf 'No dependency archives were downloaded.\n' >&2
    exit 69
fi
for archive in "${archives[@]}"; do
    package="$(dpkg-deb -f "$archive" Package)"
    dpkg-deb --extract "$archive" "$matlab_system_libs"
    dpkg-deb -f "$archive" Version > "$matlab_system_libs/.packages/$package"
done
bash "$0" --check
