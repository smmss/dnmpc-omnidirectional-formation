# Distributed NMPC for formation control of transformable omnidirectional robots

MATLAB code for the simulations of the paper

> A. Alaei, S. M. M. SeyedSajadi, H. Atrianfar, *Nonlinear Model Predictive Formation Control with Obstacle Avoidance in Transformable Omnidirectional Mobile Robots*.

The code simulates semi-wheel agents with three omni-wheels under a distributed nonlinear model predictive controller (DNMPC). Each agent solves a local optimal control problem with IPOPT (through CasADi), subject to input bounds, obstacle constraints, and inter-agent clearance and compatibility constraints. The repository reproduces the three nominal scenarios of Section 5 and the Monte Carlo study of Section 5.8.

## Requirements

- MATLAB R2022b or later (the paper used R2022b) with the Optimization Toolbox and the Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox, only for running the Monte Carlo study in parallel.
- [CasADi](https://web.casadi.org/get/) for MATLAB (version 3.8.2 was used). It is not included in this repository: unzip the CasADi MATLAB release into a folder `casadi/` in the repository root, or add it to the MATLAB path yourself.

## Quick start

```matlab
setup_paths                 % adds src/ and casadi/ to the path
run_scenarios               % three nominal scenarios (about 10 min in total)
plot_results('nominal')     % figures in results/figures/
```

## Nominal scenarios (Sections 5.2-5.4)

`run_scenarios` simulates

| Scenario | Name in the code | Section |
|---|---|---|
| Leaderless consensus on a directed cycle | `leaderless` | 5.2 |
| Leader-follower consensus with two static obstacles | `leader_consensus` | 5.3 |
| Square formation (side 1.2 m) following a virtual leader on a circle | `formation` | 5.4 |

It saves the closed-loop trajectories in `results/nominal.mat` and prints the settling time, final error, smallest inter-agent distance and mean solve time of each scenario (in the two consensus scenarios, the agents converge to a common point, so that the smallest distance is zero). Expected output for the formation: settling time 3.93 s and smallest distance 0.897 m. A single scenario can be run with `run_scenarios({'formation'})`. All parameters (weights, horizon, graph, initial poses, obstacles) are defined in `src/scenario.m`.

## Monte Carlo study (Section 5.8, Table 6)

```matlab
run_montecarlo                                   % 7 conditions x 1000 runs, 8 workers
run_montecarlo('runs',10,'workers',0, ...
               'conditions',{'nominal'})         % quick serial test
plot_results('montecarlo')                       % box plots of error and solve time
```

The formation scenario is shortened to 20 s and started from 1000 random initial conditions (the same for all conditions, fixed seed) under seven conditions: nominal; measurement noise; wheel slip and model mismatch; message delay with packet loss; failure of the links 3→2 and 1→4 at t = 8 s, with and without the sensing-based clearance; and all imperfections combined. At the end, the script prints the columns of Table 6. A run counts as converged if its final error is below 2 cm. Under slip and model mismatch, the steady error is about 1 cm, since the controller uses the nominal model, so no settling time is given for these conditions.

Every finished run is saved in `results/montecarlo_full_runs/`, and an interrupted study continues where it stopped. The full study took about 13.5 h with 8 parallel workers (Intel Core Ultra 7 265K). One worker needs up to about 3.5 GB of memory; reduce `'workers'` on machines with less memory.

The option `'graph'` repeats the study with the alternative formation graphs compared in the response to the reviewers: `'chain'` (leader → 1 → 2 → 3 → 4) and `'chain14'` (the chain plus the link 1 → 4).

## Files

| File | Content |
|---|---|
| `setup_paths.m` | adds `src/` and CasADi to the path |
| `run_scenarios.m` | nominal scenarios of Sections 5.2-5.4 |
| `run_montecarlo.m` | Monte Carlo study of Section 5.8 |
| `plot_results.m` | figures of the nominal scenarios and of the Monte Carlo study |
| `src/scenario.m` | parameters of all scenarios |
| `src/dnmpc_sim.m` | closed-loop simulation: plant, communication, local problems, fallback |
| `src/casadi_local.m` | local optimal control problem of one agent (CasADi/IPOPT) |
| `src/casadi_central.m` | centralized NMPC used for comparison |
| `src/run_metrics.m` | performance indices of a run |

Solve times depend on the computer and on the number of simulations running in parallel.

## Citation

If you use this code, please cite the paper above. The reference will be updated after publication.

## License

MIT License, see `LICENSE`.
