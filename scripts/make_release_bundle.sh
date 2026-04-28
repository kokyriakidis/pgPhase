#!/usr/bin/env bash
set -euo pipefail

# One-command release flow:
# 1) build
# 2) (optional) run validation checks
# 3) create portable bundle
# 4) create compressed tarball of that bundle
#
# Usage:
#   scripts/make_release_bundle.sh [bundle_dir]
#
# Default bundle_dir:
#   dist/pgphase-<os>-<arch>

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
case "$OS_NAME" in
  Linux)  PLATFORM="linux-$ARCH_NAME" ;;
  Darwin) PLATFORM="macos-$ARCH_NAME" ;;
  *)
    echo "[release] unsupported platform: $OS_NAME" >&2
    exit 1
    ;;
esac

BUNDLE_DIR="${1:-$ROOT_DIR/dist/pgphase-$PLATFORM}"
BUNDLE_BASENAME="$(basename "$BUNDLE_DIR")"
TARBALL_PATH="$ROOT_DIR/dist/${BUNDLE_BASENAME}.tar.gz"

RUN_CHECKS="${RUN_CHECKS:-0}"

if [[ "$RUN_CHECKS" == "1" ]]; then
  echo "[release] build + checks (strict)"
  make -C "$ROOT_DIR" check
else
  echo "[release] build (checks skipped; set RUN_CHECKS=1 for strict mode)"
  make -C "$ROOT_DIR" pgphase
fi

echo "[release] create portable bundle"
bash "$ROOT_DIR/scripts/make_portable_bundle.sh" "$BUNDLE_DIR"

echo "[release] create tarball: $TARBALL_PATH"
mkdir -p "$ROOT_DIR/dist"
tar -C "$ROOT_DIR/dist" -czf "$TARBALL_PATH" "$BUNDLE_BASENAME"

echo "[release] done"
echo "[release] bundle : $BUNDLE_DIR"
echo "[release] tarball: $TARBALL_PATH"
