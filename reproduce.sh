#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python graph_checks.py
g++ -O3 -std=c++17 random_checks.cpp -o random_checks
./random_checks
g++ -O3 -std=c++17 all_pair_checks.cpp -o all_pair_checks
./all_pair_checks
python sparse_checks.py
python verify_paulis.py
OPENBLAS_NUM_THREADS=1 python benchmark_generator.py
python noise_budget.py
python make_figures.py
(cd figures && pdflatex -interaction=nonstopmode -halt-on-error protocol.tex)
