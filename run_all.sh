#!/usr/bin/env bash
# Usage:  ./run_all.sh          -> redraw all figures from the shipped data/ (about 1 min)
#         ./run_all.sh full     -> recompute every simulation from scratch, then redraw (slower)
# Results are written to ./output/ (the shipped figures/ folder is never overwritten).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-quick}"
OUT="$HERE/output"; mkdir -p "$OUT"; cd "$OUT"
export PYTHONPATH="$HERE/code"
PY="${PYTHON:-python3}"

if [ "$MODE" = "full" ]; then
  echo "[full] folding / melting curves (Figs 2-3)";  "$PY" "$HERE/code/fig23.py"
  echo "[full] R3C parameter scan (Fig 4)";           "$PY" "$HERE/code/r3c_scan.py" | tee r3c_scan_summary.txt
  echo "[full] stochastic-corrector model (Figs 5-6)"; "$PY" "$HERE/code/fig56.py"
  echo "[full] Tm tables and Donnan/pH solver";       "$PY" "$HERE/code/run_tm.py" | tee run_tm.txt
  "$PY" "$HERE/code/duplex_tm.py" | tee duplex_tm.txt
  "$PY" "$HERE/code/ph_gap.py"    | tee ph_gap.txt
else
  cp "$HERE"/data/*.json .
fi

echo "[plot] figures"
"$PY" "$HERE/code/plot_fig1b.py"
"$PY" "$HERE/code/plot_fig2.py"
"$PY" "$HERE/code/plot_fig3.py" | tee fig3_values.txt
"$PY" "$HERE/code/plot_fig4.py"
"$PY" "$HERE/code/plot56.py"
echo "Done. Figures in $OUT"
