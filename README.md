# Blockchain Parallel Validation & Construction

This repository accompanies the paper **“Exploiting Multi-Core Parallelism in Blockchain Validation and Construction”** by Arivarasan Karmegam, Lucianna Kiffer, and Antonio Fernández Anta. The code reproduces the paper’s optimisation baselines  and heuristic schedulers for both parallel block execution and parallel block construction on multicore hardware.

## Core ideas

- **Sequential execution wastes cores.** Ethereum-style block validation traditionally runs transactions one by one even when their read/write sets do not conflict. The provided solvers and heuristics expose safe parallelism while maintaining determinism.
- **Two validator problems.**
  - **Problem 1 – Parallel block execution.** Given an ordered list of transactions, assign them to `p` cores to minimize the makespan while enforcing all read/write dependencies.
  - **Problem 2 – Parallel block construction.** Select and schedule a subset of mempool transactions that maximizes reward under conflict, gas/runtime, and processor constraints.
- **Exact + heuristic methods.** The repository includes MILP baselines (solved via MATLAB's `intlinprog` with the HiGHS backend) and fast greedy heuristics that achieve near-optimal schedules at production-ready speeds.

## Repository layout

```
README.md                          Project overview (this file)

simple/                            MILP + heuristic scripts for "simple" transaction models
simple/automated_solver_p1.m       Batch MILP launcher for Problem 1
simple/automated_solver_p2.m       Batch MILP launcher for Problem 2
simple/Heuristics_P1_simple.m      Deterministic greedy scheduler for Problem 1
simple/Heuristics_P2_simple.m      Deterministic greedy selector/scheduler for Problem 2
simple/longest_chain_length.m      Utility for computing DAG heights inside heuristics
simple/solver_p1.m                 Preprocessing for Problem 1 MILP
simple/solver_p2.m                 Preprocessing for Problem 2 MILP
simple/solver_min_rounds_p1.m      Core MILP formulation for Problem 1
simple/solver_min_rounds_p2.m      Core MILP formulation for Problem 2
simple/ranges.csv                  Sample window definitions for automated runs

simple and complex/                Mixed-workload (simple + complex gas models) scripts
simple and complex/automated_solver_p1.m     Batch MILP launcher for Problem 1
simple and complex/automated_solver_p2.m     Batch MILP launcher for Problem 2
simple and complex/Heuristics_P1_complex.m   Complex-weighted P1 heuristic
simple and complex/Heuristics_P3_complex.m   Complex-weighted P2 heuristic
simple and complex/longest_chain_length.m    Shared DAG utility for heterogeneous heuristics
simple and complex/solve_min_rounds_p1.m     Core MILP formulation for Problem 1
simple and complex/solve_min_rounds_p2.m     Core MILP formulation for Problem 2
simple and complex/solver_p1.m               Preprocessing for Problem 1 MILP
simple and complex/solver_p2.m               Preprocessing for Problem 2 MILP
simple and complex/p1_batch_params.csv       Example batch configuration input for Problem 1
simple and complex/batch_params.csv          Example batch configuration input for Problem 2
```

## Requirements

- MATLAB R2025b (tested) with Optimization Toolbox.
- HiGHS mixed-integer solver, accessed through `intlinprog` (ships with MATLAB R2025b).
- Ethereum trace CSV files containing transaction hashes, read/write sets, gas usage, and reward weights (gasUsed x tipPerGas).

## How to Run

1. Place the required CSV traces in the corresponding folder (simple or simple and complex)
2. Open MATLAB, add the relevant folder (`simple/` or `simple and complex/`) to the path, and edit the configuration block at the top of the script you wish to run (e.g., row ranges, number of cores `p`, gas/runtime budgets `R`, and reward columns).
3. Run the script (F5 or `run`) to generate logs and summary CSVs under the configured `logs` directory. Summary tables capture metrics such as number of transactions, detected conflicts, MILP objective value, rounds used, and solver status.

## Data expectations

CSV inputs should expose, at a minimum:

- `tx_hash` – unique identifier for display/logging.
- `access_read`, `access_write` – semi-colon-separated storage keys (strings) used to construct conflict graphs.
- `gas_used` – execution cost per transaction.
- `weight_eth` or `weight_wei` - rewards in different units. We use eth in the paper.

## Reproducing paper experiments

To replicate the paper’s figures:

1. Acquire Ethereum mainnet traces that match the block ranges studied in the paper.
2. Configure the automated MILP scripts with the appropriate windows (`ranges.csv` / `batch_params.csv`) and run them to obtain optimal baselines.
3. Execute the heuristic scripts on the same slices to gather makespan/reward metrics.
4. Compare heuristic outputs against MILP summaries to calculate approximation ratios and speedups.

The repository does not include Ethereum transaction traces due to their large size.
Users should supply their own CSV exports matching the schema described above (tx_hash, access_read, access_write, gas_used, weight_eth).

## License

All MATLAB scripts are distributed under the GNU General Public License v3.0, as noted in the header of each file.
