#!/bin/bash
set -euo pipefail

mkdir -p temp_dir

echo "=== Compiling solvers (without compiler optimizations) ==="
# ------------------------------------
# Version 0 — FORTRAN 77 Legacy Solver
# ------------------------------------
gfortran solver_0.for -o temp_dir/sol0
echo "solver_0.for Compiled"

# ------------------------------------
# Version 1 — Direct C++ Translation
# ------------------------------------
g++ solver_1.cpp -o temp_dir/sol1
echo "solver_1.cpp Compiled"

# ------------------------------------
# Version 2 — Loop Reordered
# ------------------------------------
g++ solver_2.cpp -o temp_dir/sol2
echo "solver_2.cpp Compiled"

# ------------------------------------
# Version 3 — Flattened Memory
# ------------------------------------
g++ solver_3.cpp -o temp_dir/sol3
echo "solver_3.cpp Compiled"

# ------------------------------------
# Version 4 — Blitz++ Arrays
# ------------------------------------
g++ solver_4.cpp -o temp_dir/sol4
echo "solver_4.cpp Compiled"

# ------------------------------------
# Version 5 — Blitz++ Linear + Pointer
# ------------------------------------
g++ solver_5.cpp -o temp_dir/sol5
echo "solver_5.cpp Compiled"

# ------------------------------------
# Copy Input Data to tmp Folder
# ------------------------------------
cp INP.DAT temp_dir/
echo "INP.DAT Copied to ./temp_dir"

echo
echo "=== Simulating Solver with Perf ==="
# events: cache counters + instructions + branches (core + atom)
EVENTS=\
"cpu_core/cache-references/,cpu_atom/cache-references/,"\
"cpu_core/cache-misses/,cpu_atom/cache-misses/,"\
"cpu_core/L1-dcache-loads/,cpu_atom/L1-dcache-loads/,"\
"cpu_core/L1-dcache-load-misses/,"\
"cpu_core/LLC-loads/,cpu_atom/LLC-loads/,"\
"cpu_core/LLC-load-misses/,cpu_atom/LLC-load-misses/,"\
"cpu_core/LLC-stores/,cpu_core/LLC-store-misses/,"\
"cpu_core/instructions/,cpu_atom/instructions/,"\
"cpu_core/branch-instructions/,cpu_atom/branch-instructions/,"\
"cpu_core/branch-misses/,cpu_atom/branch-misses/"

# If perf is not in PATH, user must adjust.
for i in 0 1 2 3 4 5
do
    echo " Running perf for solver_${i} ..."
    perf stat -e $EVENTS -o temp_dir/sol${i}_perf.txt ./temp_dir/sol${i} || {
        echo "perf failed for sol${i}. Check perf installation and permissions (may need sudo)."
        exit 1
    }
    echo "  saved: sol${i}_perf.txt"
done

echo
echo "=== All done ==="
