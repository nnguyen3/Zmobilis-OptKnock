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

## Key Isobutanol Results

Wild-type growth:

```text
WT growth = 0.059677
```

Wild-type maximum isobutanol production:

```text
WT max isobutanol flux = 10.030751
```

During normal growth, the WT model produced zero isobutanol:

```text
WT isobutanol flux during growth = 0.000000
```

This showed that isobutanol production was feasible, but not naturally growth-coupled.

After forcing minimum isobutanol production, OptKnock identified `ALCD2x` as the strongest knockout candidate.

Manual validation showed:

```text
ALCD2x KO growth = 0.036737
ALCD2x KO isobutanol flux during growth = 9.671502
```

Compared to the forced WT condition:

```text
Forced WT isobutanol flux during growth = 1.000000
ALCD2x KO isobutanol flux during growth = 9.671502
Increase = 8.671502
```

These results suggest that `ALCD2x` may improve growth-associated isobutanol production in the anaerobic *Z. mobilis* model.

### FluxRETAP-inspired Analysis

A FluxRETAP-inspired flux variability analysis (FVA) workflow
was performed to compare low- and high-isobutanol production states.

Workflow:
- Gradually forced increasing isobutanol production
- Ran Flux Variability Analysis (FVA) at each production level
- Compared reaction flux distributions between conditions
- Ranked reactions based on flux changes

Key finding:
`ALCD2x` showed one of the largest flux changes between
low- and high-isobutanol conditions. As isobutanol production
increased, flux through the ethanol synthesis pathway decreased,
suggesting competition between ethanol and isobutanol production pathways.

This result reinforced the previous OptKnock and manual validation results,
supporting `ALCD2x` as the strongest engineering target identified in this study.

## Requirements

* MATLAB
* COBRA Toolbox
* Gurobi solver

## Notes

Large `.mat` model files are included for reproducibility.
