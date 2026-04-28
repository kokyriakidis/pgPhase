#!/usr/bin/env bash
set -euo pipefail

# Build a portable runtime bundle for pgphase (Linux/macOS) by copying runtime
# shared libraries next to the binary and generating a launcher script that
# sets platform-appropriate library search paths.
#
# Usage:
#   scripts/make_portable_bundle.sh [output_dir]
#
# Default output:
#   dist/pgphase-<os>-<arch>

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
case "$OS_NAME" in
  Linux)   PLATFORM="linux-$ARCH_NAME" ;;
  Darwin)  PLATFORM="macos-$ARCH_NAME" ;;
  *)
    echo "[portable] unsupported platform: $OS_NAME" >&2
    exit 1
    ;;
esac

OUT_DIR="${1:-$ROOT_DIR/dist/pgphase-$PLATFORM}"
BIN_DIR="$OUT_DIR/bin"
LIB_DIR="$OUT_DIR/lib"
BIN_SRC="$ROOT_DIR/pgphase"

echo "[portable] building pgphase"
make -C "$ROOT_DIR" pgphase

echo "[portable] preparing bundle at: $OUT_DIR"
rm -rf "$OUT_DIR"
mkdir -p "$BIN_DIR" "$LIB_DIR"

cp "$BIN_SRC" "$BIN_DIR/pgphase"

collect_linux_deps() {
  ldd "$BIN_SRC" | awk '
    /=> \// { print $3 }
    /^\// { print $1 }
  ' | sort -u
}

collect_macos_deps() {
  # Keep non-system libs (e.g. Homebrew/local installs). System libs under
  # /usr/lib and /System are expected on macOS and are not bundled.
  otool -L "$BIN_SRC" | awk 'NR>1 { gsub(/^[[:space:]]+/, "", $0); print $1 }' \
    | awk '
      /^\// {
        if ($0 ~ "^/usr/lib/" || $0 ~ "^/System/") next;
        print $0;
      }
    ' | sort -u
}

if [[ "$OS_NAME" == "Linux" ]]; then
  mapfile -t DEPS < <(collect_linux_deps)
else
  mapfile -t DEPS < <(collect_macos_deps)
fi

echo "[portable] copying ${#DEPS[@]} shared libraries for $OS_NAME"
for lib in "${DEPS[@]}"; do
  if [[ -f "$lib" ]]; then
    cp -L "$lib" "$LIB_DIR/"
  else
    echo "[portable] warning: dependency not found on disk: $lib" >&2
  fi
done

# Create launcher that pins runtime library search to bundled libs first.
if [[ "$OS_NAME" == "Linux" ]]; then
cat > "$BIN_DIR/pgphase-portable" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_DIR="$SELF_DIR/../lib"
exec env LD_LIBRARY_PATH="$LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$SELF_DIR/pgphase" "$@"
EOF
else
cat > "$BIN_DIR/pgphase-portable" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_DIR="$SELF_DIR/../lib"
exec env DYLD_LIBRARY_PATH="$LIB_DIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" \
  "$SELF_DIR/pgphase" "$@"
EOF
fi
chmod +x "$BIN_DIR/pgphase-portable"

# Record build/runtime metadata for users.
{
  echo "pgPhase portable bundle"
  echo
  echo "Build host: $(uname -srvmo)"
  echo "Build date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo
  echo "Run with:"
  echo "  ./bin/pgphase-portable <args>"
  echo
  echo "Bundle contents:"
  echo "  bin/pgphase"
  echo "  bin/pgphase-portable"
  if [[ "$OS_NAME" == "Linux" ]]; then
    echo "  lib/*.so* (copied from ldd)"
  else
    echo "  lib/*.dylib* (copied from otool -L, excluding system libs)"
  fi
  echo
  echo "Compatibility note:"
  if [[ "$OS_NAME" == "Linux" ]]; then
    echo "  This bundle is portable across Linux systems with compatible kernel/CPU and glibc ABI."
    echo "  For maximum portability across older distributions, build on the oldest target distro."
  else
    echo "  This bundle is portable across macOS systems with compatible CPU architecture and runtime ABI."
    echo "  Build separate bundles for arm64 (Apple Silicon) and x86_64 (Intel) when needed."
  fi
} > "$OUT_DIR/README.portable.txt"

echo "[portable] done"
echo "[portable] launcher: $BIN_DIR/pgphase-portable"
