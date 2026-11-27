#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
run_dir="$(mktemp -d)"
trap 'rm -rf "$run_dir"' EXIT
run_notebook() {
  jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=-1 \
    --ExecutePreprocessor.kernel_name=python3 \
    --output-dir="$run_dir" "$1"
}
run_notebook graph_checks.ipynb
g++ -O3 -std=c++17 random_checks.cpp -o random_checks
./random_checks
g++ -O3 -std=c++17 all_pair_checks.cpp -o all_pair_checks
./all_pair_checks
run_notebook sparse_checks.ipynb
run_notebook verify_paulis.ipynb
OPENBLAS_NUM_THREADS=1 run_notebook benchmark_generator.ipynb
run_notebook noise_budget.ipynb
run_notebook make_figures.ipynb
(cd figures && pdflatex -interaction=nonstopmode -halt-on-error protocol.tex)
