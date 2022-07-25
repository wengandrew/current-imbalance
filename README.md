# current-imbalance
Current imbalance dynamics in parallel-connected cells.

May 14, 2022

MATLAB 2021b was used to generate the code and run the simulations. To run, set the path to the root. Do not run from `src` since the code assumes you are starting from one level up from `src.)


## Next steps

- [x] Understand why current imbalance does not depend on resistance
- [x] Implement CV hold dynamics
- [ ] Make plots to separate `I_Ohmic` from `I_Rebalance`
- [ ] Develop a simple model to link current imbalance and SOC imbalance to an incremental degradation per cycle
- [ ] Run `z_start --> `z_end` sweeps for both the linear and non-linear case. For the non-linear case, consider both LFP and NMC chemistries.
- [ ] Explore the impact of SOC-dependent resistances.
