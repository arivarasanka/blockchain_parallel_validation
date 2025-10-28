# Blockchain Parallel Validation

Tools in this repository analyse how blockchain transactions can be validated in parallel when conflicts are constrained by read/write access lists.  The MATLAB scripts focus on two mixed-integer optimisation problems:

- **Problem 1 (P1 - minimum rounds):** assign each transaction in a window to one of *R* rounds while observing a per-round capacity *p* and the precedence edges implied by conflicts.  The objective is to minimise the number of rounds that are ultimately used.
- **Problem 2 (P2 - maximum weight):** given the same conflict graph, select up to *p* transactions per round (across *R* rounds) in order to maximise the total transaction weight.

Both formulations are solved with MATLAB's `intlinprog`, and the repository supplements the MILPs with heuristics, batch runners, and log processing helpers that make it practical to study large traces of Ethereum transactions.

## Key capabilities

- Build conflict graphs from CSV access-list traces and compute precedence bounds automatically.
- Run single-window experiments for interactive exploration or batch sweeps over many transaction ranges.
- Generate logs and summary tables that capture solver statistics, objectives, and feasibility checks.
- Compare optimal schedules against fast greedy heuristics for back-of-the-envelope benchmarking.

## Repository layout

```
.
├── README.md              ← project documentation (this file)
└── simple/
    ├── automated_solver_p1.m   ← batch runner for P1 windows defined in `ranges.csv`
    ├── automated_solver_p2.m   ← batch runner for P2 windows defined in `ranges.csv`
    ├── Heuristics_P1_simple.m  ← greedy/heuristic schedulers for P1
    ├── Heuristics_P2_simple.m  ← heuristic selectors for P2
    ├── longest_chain_length.m  ← helper for computing precedence lower bounds
    ├── ranges.csv              ← example list of `[startId,endId,p,slackR/R]` windows
    ├── solve_min_rounds_p1.m   ← MILP model for P1 (minimise used rounds)
    ├── solver_min_rounds_p2.m  ← MILP model for P2 (maximise scheduled weight)
    ├── solver_p1.m             ← example interactive driver for a single P1 window
    └── solver_p2.m             ← example interactive driver for a single P2 window
```

All MATLAB code lives under `simple/` and only depends on core MATLAB plus the Optimization Toolbox.

## Requirements

- MATLAB R2022b or newer (earlier releases may work but are not validated).
- Optimization Toolbox (for `intlinprog`).
- Enough RAM to hold the selected slice of the access-list CSV in memory.

## Preparing input data

The scripts consume a CSV describing the access list for each transaction.  For best results, include at least the following columns:

| Column | Required | Notes |
| --- | --- | --- |
| `tx_hash` | Optional | Used for logging human-readable identifiers. |
| `blockNumber` | Optional | Helpful when cross-referencing schedules with chain height. |
| `access_read` | **Yes** | Semicolon-separated list of storage keys that are read. |
| `access_write` | **Yes** | Semicolon-separated list of storage keys that are written. |
| `weight_theoretical_eth` | Optional (P2) | Objective weight used by P2 MILP and heuristics. |

Batch experiments additionally require a `ranges.csv` with one row per window:

| Column | Description |
| --- | --- |
| `startId`, `endId` | 1-indexed (inclusive) row IDs in the main access-list CSV. |
| `p` | Per-round capacity. |
| `slackR` | (P1, optional) Extra rounds beyond the longest precedence chain. Falls back to `defaultSlackR` when omitted. |
| `R` | (P2) Number of rounds available for scheduling. |

The example file at `simple/ranges.csv` shows the required headers and formatting.  When preparing your own datasets, keep the main CSV in the repository root or adjust the `csvFile` parameter inside the scripts.

## Quick start

1. **Add the toolbox to your MATLAB path**

   ```matlab
   addpath('simple');
   ```

2. **Point the configuration to your data** by editing the "Config" block at the top of `solver_p1.m`, `solver_p2.m`, `automated_solver_p1.m`, or `automated_solver_p2.m`.  Update `csvFile`, the target window, capacity `p`, and (for batches) `defaultSlackR`/`summaryOut` as needed.【F:simple/automated_solver_p1.m†L8-L15】【F:simple/automated_solver_p2.m†L8-L14】
3. **Load and run a window** using either the single-window drivers or the automated batch scripts.

### Single-window workflows

- **P1 interactive solver:** `solver_p1.m` constructs the conflict graph for `[startId..endId]`, invokes `solve_min_rounds_p1.m`, checks feasibility, and prints the resulting schedule summary.【F:simple/solver_p1.m†L1-L125】
- **P2 interactive solver:** `solver_p2.m` solves the weighted selection problem and reports the chosen transactions, achieved weight, and per-round utilisation.【F:simple/solver_p2.m†L1-L132】

### Batch automation

- **P1 sweeps:** `automated_solver_p1.m` iterates over every row in `ranges.csv`, logging the run to `simple/logs/` and collecting metrics such as objective value, precedence bounds, exit flag, and runtime.【F:simple/automated_solver_p1.m†L1-L86】
- **P2 sweeps:** `automated_solver_p2.m` mirrors the P1 workflow for the weighted formulation and produces analogous log and summary outputs.【F:simple/automated_solver_p2.m†L1-L85】

Batch runners create `simple/logs/` on demand.  Clear the folder between large experiments to keep disk usage manageable.【F:simple/automated_solver_p1.m†L16-L18】

## Outputs

- **Log files:** Each batch iteration streams MATLAB output to `simple/logs/p?_*.log`, making it easy to audit solver progress and configuration.
- **Summary CSVs:** After a batch completes, `summary_p1_ranges.csv` or `summary_p2_ranges.csv` records statistics per window, including solver exit flags, objective values, bounds, and runtimes.【F:simple/automated_solver_p1.m†L72-L86】【F:simple/automated_solver_p2.m†L70-L83】
- **Schedule matrices:** The interactive drivers return MATLAB matrices (rows = rounds, columns = transaction IDs) that you can post-process for visualisation or downstream simulation.

## Heuristics and helper utilities

- `Heuristics_P1_simple.m` and `Heuristics_P2_simple.m` provide greedy baselines that yield feasible schedules quickly, allowing you to compare optimal vs. heuristic performance.【F:simple/Heuristics_P1_simple.m†L1-L120】【F:simple/Heuristics_P2_simple.m†L1-L118】
- `longest_chain_length.m` computes the length of the longest precedence chain, which acts as a natural lower bound on the number of rounds needed before launching the MILP.【F:simple/longest_chain_length.m†L1-L92】
- `solve_min_rounds_p1.m` and `solver_min_rounds_p2.m` expose the underlying MILP model definitions, making it straightforward to tweak solver tolerances or add custom constraints.【F:simple/solve_min_rounds_p1.m†L1-L110】【F:simple/solver_min_rounds_p2.m†L1-L111】

## Troubleshooting

- **Empty windows:** Ensure `startId` ≤ `endId` and the range falls within the bounds of the main CSV; otherwise the batch scripts will raise an error when slicing the table.【F:simple/automated_solver_p1.m†L40-L66】
- **Missing columns:** If your CSV lacks `access_read` or `access_write`, MATLAB will issue an error before solving.  Regenerate the dataset with the required headers.【F:simple/automated_solver_p1.m†L61-L66】
- **Solver convergence:** Tune `optimoptions` inside the MILP scripts to adjust time limits, tolerances, or branching heuristics for especially large windows.【F:simple/solve_min_rounds_p1.m†L72-L93】【F:simple/solver_min_rounds_p2.m†L79-L110】

## Committing your updates

When you are ready to save your local changes to this repository, run the following commands from the project root:

1. **Inspect the working tree** to review modified files.

   ```bash
   git status
   ```

2. **Stage the updates** you want to include in the commit.

   ```bash
   git add <files>
   ```

   You can stage everything at once with `git add .`, but staging individual files keeps the commit history focused.

3. **Create the commit** with a descriptive message.

   ```bash
   git commit -m "Describe what changed"
   ```

4. **Push the commit** to the remote repository.

   ```bash
   git push origin <branch-name>
   ```

Replace `<branch-name>` with the branch you are working on (for example, `main` or `work`).  Once the push succeeds, your commit will be available in the remote repository and ready for a pull request or deployment.

## Contributing

Contributions are welcome!  Please open an issue or submit a pull request describing the change you have in mind.  When adding new scripts or modifying solver behaviour, update this README so others can quickly discover the new capabilities.

