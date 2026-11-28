# Fundamental Limits and Robust Certification of Genuine Multipartite Steering

Pritam Roy  
S. N. Bose National Centre for Basic Sciences, Kolkata, India

This repository contains the numerical calculations and figure sources accompanying the manuscript of the same title. The calculations concern genuine multipartite steering (GMS) with prescribed untrusted qubit measurements, graph-state stabilizer certificates, and the exact visibility threshold of a globally depolarized five-cycle with two adjacent untrusted vertices.

## Reproducing the calculations

The Python notebooks require the packages in `requirements.txt`. The two rank-enumeration programs require a C++17 compiler (GCC or compatible). Figure generation uses Matplotlib with LaTeX text rendering (including type1cm, cm-super and dvipng); the circuit diagram requires pdfLaTeX and TikZ.

From the repository root:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
bash reproduce.sh
```

The script runs the notebooks and C++ programs in dependency order and writes the generated CSV files and PDF figures locally. The notebooks can also be run individually in JupyterLab. Run the graph enumeration and rank tests before the reduced-check and plotting notebooks. The all-pair tests for larger graphs can be computationally expensive.

## Calculations

| File | Calculation | Manuscript |
|---|---|---|
| `verify_paulis.ipynb` | Pauli identities, fixed-frame non-GMS construction and local noise thresholds | Secs. III–V; Appendices A, C and E |
| `verify_exact.ipynb` | Qutrit counterexample, equality of the five-cycle assemblages at the critical visibility, local dephasing identities, generator spectra, calibration and dense-graph noise checks | Secs. III A and V–VI; Appendices G and H |
| `graph_checks.ipynb` | Exhaustive census of connected nonisomorphic graphs with 3–7 vertices | Appendix F 2, Table VI |
| `random_checks.cpp` | Fixed-pair and all-pair rank tests for random graphs | Fig. 6(b) |
| `all_pair_checks.cpp` | All-pair rank tests for graphs with 17–24 vertices | Fig. 6(b), Appendix F 2 |
| `sparse_checks.ipynb` | Reduced sets of trusted stabilizer checks | Appendix F 3, Table VII |
| `benchmark_generator.ipynb` | Numerical search for feasible generator-sum strategies | Sec. V; Appendix H |
| `noise_budget.ipynb` | Gate and readout noise models and sufficient finite-shot budgets | Fig. 5, Table IV |
| `make_figures.ipynb` | Vector plots for Figs. 1–3, 5 and 6 | Main text |
| `figures/protocol.tex` | Preparation circuit, measurement settings and nine-correlator estimator | Fig. 4 |

The graph-atlas census is exhaustive. The random-graph calculations use the following sampling parameters.

| Calculation | Samples | Seed |
|---|---|---|
| Fixed pair | 100 per size, 6–20 vertices | 20261127 |
| All pairs (initial range) | 100 per size, 6–14 vertices; 20 per size, 15–16 vertices | 20261127 |
| All pairs (larger graphs) | 20 per size, 17–20 vertices; 10 per size, 21–24 vertices | 20261128 |
| Reduced trusted-check sets | Up to 3000 candidate subsets, first rank-certified graph at each size | 20261127 |
| Generator-sum search | 12 starts, at most 80 iterations, 15 bipartitions | 20261127 |

The larger-graph program includes a self-test on 140 graphs (seed 203). The illustrative dense graph in Fig. 2 uses seed 15.

## Outputs

- `atlas_results.csv`, `atlas_summary.csv`: graph-atlas census.
- `random_graphs.csv`, `random_summary.csv`, `all_pair_extension.csv`: random-graph rank tests.
- `sparse_checks.csv`: reduced stabilizer-check sets.
- `generator_benchmark.csv`, `noise_budget.csv`: generator search and statistical estimates.
- `figures/`: seven plot PDFs and the circuit PDF. Fig. 3 uses `noise.pdf` and `visibility_geometry.pdf`; Fig. 6 uses `random.pdf` and `finite.pdf`.

Adjacency masks and stabilizer labels are stored as decimal representations of binary strings, with bit `j` corresponding to vertex `j`. An all-pair flag of `-1` in `random_graphs.csv` denotes a test that was not performed.

The rank and stabilizer criteria are sufficient certificates: failure of a test does not imply absence of GMS. The exact five-cycle visibility threshold and generator-sum bound are proved analytically in the manuscript; the notebooks check the corresponding identities. Local noise and finite-shot estimates refer to the noise and independent-trial models specified there. No experimental data are used.

## Citation

Please cite *Fundamental Limits and Robust Certification of Genuine Multipartite Steering* (Pritam Roy) and the repository commit used for the calculation.
