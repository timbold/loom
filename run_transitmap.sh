#!/usr/bin/env bash
# run_transitmap.sh — wrapper for topo → loom → transitmap with progress + pretty SVG
#
# Usage:
#   ./run_transitmap.sh -i ../examples/5shar.json -o ../output/5shar-full.svg -m "5шарурд"
#
# Env overrides:
#   TOOLS_DIR         Directory with binaries (default: <script_dir>/build, falls back to <script_dir>)
#   TOPO_ARGS         Args for topo        (default: --smooth 1)
#   LOOM_ARGS         Args for loom        (defaults set below)
#   TRANSITMAP_ARGS   Args for transitmap  (defaults set below; excludes --me-station)

set -euo pipefail

# ---------- pretty console ----------
bold()  { printf "\033[1m%s\033[0m\n" "$*"; }
green() { printf "\033[32m%s\033[0m\n" "$*"; }
blue()  { printf "\033[34m%s\033[0m\n" "$*"; }
red()   { printf "\033[31m%s\033[0m\n" "$*"; }

# ---------- defaults ----------
INPUT="../examples/5shar.json"
OUTPUT="../output/5shar-full.svg"
ME_STATION="5шарурд"

# Ensure UTF-8 for non-ASCII station names
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

# Resolve script dir and default tools dir to "<script_dir>/build"
SCRIPT_PATH="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd -- "$(dirname -- "$SCRIPT_PATH")" && pwd)"
: "${TOOLS_DIR:="${SCRIPT_DIR}/build"}"

# Default args (match your one-liner)
: "${TOPO_ARGS:="--smooth 1"}"
: "${LOOM_ARGS:="--in-stat-cross-pen-diff-seg=4 --in-stat-cross-pen-same-seg=2 --in-stat-sep-pen=9 --same-seg-cross-pen=6 --diff-seg-cross-pen=1 --sep-pen=2"}"
: "${TRANSITMAP_ARGS:="-l -r --render-dir-markers --render-markers-tail --highlight-terminal --tl-ratio=2.1 --compact-route-label --tail-ignore-sharp-angle --crowded-line-thresh=1 --dir-marker-spacing=2"}"

# ---------- help ----------
show_help() {
  cat <<EOF
$(bold "run_transitmap.sh") — wrapper for topo → loom → transitmap

Usage:
  $0 -i <input.json> -o <output.svg> -m "<me-station>"

Options:
  -i   Input JSON (default: ${INPUT})
  -o   Output SVG (default: ${OUTPUT})
  -m   Station for --me-station (default: ${ME_STATION})
  -h   Show this help

Env overrides:
  TOOLS_DIR         Directory with binaries (topo, loom, transitmap)
  TOPO_ARGS         Args for topo
  LOOM_ARGS         Args for loom
  TRANSITMAP_ARGS   Args for transitmap (excluding --me-station)

Example:
  TOOLS_DIR=/opt/projects/loom/build \\
  $0 -i ../examples/5shar.json -o ../output/5shar-full.svg -m "5шарурд"
EOF
}

# ---------- parse args ----------
while getopts ":i:o:m:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    m) ME_STATION="$OPTARG" ;;
    h) show_help; exit 0 ;;
    \?) red "Error: Unknown option -$OPTARG"; exit 2 ;;
    :)  red "Error: Option -$OPTARG requires an argument."; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

# ---------- binaries ----------
TOP_BIN="${TOOLS_DIR}/topo"
LOOM_BIN="${TOOLS_DIR}/loom"
TM_BIN="${TOOLS_DIR}/transitmap"

# Fallback: if not in build, try the script dir itself
if [[ ! -x "$TOP_BIN" || ! -x "$LOOM_BIN" || ! -x "$TM_BIN" ]]; then
  TOP_BIN="${SCRIPT_DIR}/topo"
  LOOM_BIN="${SCRIPT_DIR}/loom"
  TM_BIN="${SCRIPT_DIR}/transitmap"
fi

for bin in "$TOP_BIN" "$LOOM_BIN" "$TM_BIN"; do
  [[ -x "$bin" ]] || { red "Error: missing executable $bin"; exit 3; }
done

# ---------- inputs/outputs ----------
[[ -f "$INPUT" ]] || { red "Error: input not found: $INPUT"; exit 4; }
mkdir -p "$(dirname "$OUTPUT")"

TMP_SVG="$(mktemp)"
trap 'rm -f "$TMP_SVG"' EXIT

echo
bold "▶ LOOM pipeline"
blue "input      : $INPUT"
blue "output     : $OUTPUT"
blue "me-station : $ME_STATION"
blue "tools      : $(dirname "$TOP_BIN")"
echo

bold "• topo        $TOPO_ARGS"       >&2
bold "• loom        $LOOM_ARGS"       >&2
bold "• transitmap  $TRANSITMAP_ARGS --me-station=\"$ME_STATION\"" >&2

# ---------- heartbeat (used when pv is unavailable) ----------
heartbeat() {
  local parent_pid=$1
  local start_ts=$(date +%s)
  tput civis 2>/dev/null || true
  while kill -0 "$parent_pid" 2>/dev/null; do
    local bytes cpu elapsed
    bytes=$(stat -c %s "$TMP_SVG" 2>/dev/null || echo 0)
    cpu=0
    # Sum CPU of children of the subshell
    mapfile -t kids < <(pgrep -P "$parent_pid" || true)
    if (( ${#kids[@]:-0} > 0 )); then
      for k in "${kids[@]}"; do
        local c
        c=$(ps -o pcpu= -p "$k" 2>/dev/null | tr -d ' ' || true)
        [[ -n "$c" ]] && cpu=$(awk -v a="$cpu" -v b="$c" 'BEGIN{printf("%.1f", a+b)}')
      done
    fi
    elapsed=$(( $(date +%s) - start_ts ))
    printf "\r⏳ cpu %5s%% | out %12s bytes | %4ss " "${cpu:-0.0}" "$bytes" "$elapsed" >&2
    sleep 1
  done
  printf "\r\033[K" >&2
  tput cnorm 2>/dev/null || true
}

# ---------- run pipeline with multi-stage progress ----------
set -o pipefail
rc=0

if command -v pv >/dev/null 2>&1; then
  in_size=$(wc -c < "$INPUT" 2>/dev/null || echo 0)
  bold "• multi-stage pv (input, after topo, after loom)" >&2

  STDBUF=""
  if command -v stdbuf >/dev/null 2>&1; then
    STDBUF="stdbuf -oL -eL"
  fi

  {
    cat "$INPUT" \
    | pv -pte -r -b -s "$in_size" -N "read input" \
    | $STDBUF "$TOP_BIN" $TOPO_ARGS \
    | pv -pte -r -b            -N "after topo" \
    | $STDBUF "$LOOM_BIN" $LOOM_ARGS \
    | pv -pte -r -b            -N "after loom" \
    | $STDBUF "$TM_BIN" $TRANSITMAP_ARGS --me-station="$ME_STATION" \
    > "$TMP_SVG"
  } || rc=$?
else
  bold "• pv not found → heartbeat mode (install 'pv' for bars)" >&2
  (
    cat "$INPUT" \
    | "$TOP_BIN" $TOPO_ARGS \
    | "$LOOM_BIN" $LOOM_ARGS \
    | "$TM_BIN" $TRANSITMAP_ARGS --me-station="$ME_STATION" \
    > "$TMP_SVG"
  ) &
  PIPE_PID=$!
  heartbeat "$PIPE_PID"
  wait "$PIPE_PID" || rc=$?
fi

[[ $rc -eq 0 ]] || { red "✗ Pipeline failed"; exit "$rc"; }

# ---------- pretty format ----------
if command -v xmllint >/dev/null 2>&1; then
  xmllint --format "$TMP_SVG" > "$OUTPUT"
else
  mv "$TMP_SVG" "$OUTPUT"
fi

green "✓ Done: $OUTPUT"
