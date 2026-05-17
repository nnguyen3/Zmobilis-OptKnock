# Zmobilis-OptKnock

This repository contains OptKnock workflows for ethanol and isobutanol production using anaerobic Zymomonas mobilis genome-scale metabolic models (GEMs).
## Repository Structure

### ethanol/

Contains the ethanol OptKnock benchmark workflow.

Main script:

```text
optknock_ethanol_nhi_may_16.m
```

Workflow summary:

* Runs OptKnock using `EX_etoh_e` as the target reaction
* Tests broad and filtered knockout candidate lists
* Checks whether `LDH_D` exists in the model and candidate list
* Performs manual `LDH_D` knockout analysis

### isobutanol/

Contains the isobutanol OptKnock workflow.

Main script:

```text
optknock_isobutanol_nhi_0516.m
```

Workflow summary:

* Runs OptKnock using `EX_ibtol_e` as the target reaction
* Tests broad and filtered knockout candidate lists
* Removes exchange, sink, transport, and other non-ideal reactions
* Manually validates knockout candidates
* Tests whether isobutanol production is growth-coupled
* Identifies `ALCD2x` as the strongest knockout candidate under the forced-isobutanol condition

*Key Isobutanol Result
Wild-type growth:
WT growth = 0.059677
Wild-type maximum isobutanol production:
WT max isobutanol flux = 10.030751
During normal growth, the WT model produced zero isobutanol:
WT isobutanol flux during growth = 0.000000
After forcing minimum isobutanol production, OptKnock identified ALCD2x as the strongest knockout candidate.
Manual validation showed:
ALCD2x KO growth = 0.036737
ALCD2x KO isobutanol flux during growth = 9.671502
This result suggests that ALCD2x may improve growth-associated isobutanol production in the anaerobic Z. mobilis model.

## Requirements

* MATLAB
* COBRA Toolbox
* Gurobi solver

## Notes

Large `.mat` model files are included for reproducibility.
