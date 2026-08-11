# iZM_SDSU: Genome-Scale Metabolic Modeling of Zymomonas mobilis
## Citation

Nguyen, N. (2026). *nnguyen3/Zmobilis-OptKnock: iZM_SDSU Thesis Release v1.0.1* (Version v1.0.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21881065

This repository contains the final **iZM_SDSU genome-scale metabolic model (GEM)** of *Zymomonas mobilis* and computational workflows used to investigate ethanol and isobutanol production.

The repository includes the final curated metabolic model, OptKnock analyses, computational validation of predicted knockout strategies, and a FluxRETAP-inspired analysis of metabolic flux changes associated with isobutanol production.

## Repository Structure

```text
Zmobilis-OptKnock/
├── Models/
│   ├── Zm_model_april_27_anaerobic_GF_2026.mat
│   ├── Zm_model_may_06_anaerobic_GF_2026_isobutanol.mat
│   ├── iZM_SDSU.mat
│   └── iZM_SDSU.xml
│
├── ethanol/
│   └── OK_etoh_0714.mlx
│
├── isobutanol/
│   ├── OK_iso_0714.mlx
│   └── OK_iso_double_validate_0714.mlx
│
├── FluxRETAP_Isobutanol_Zmobilis/
└── README.md
```

## Final iZM_SDSU Model

The `models/` directory contains the final curated **iZM_SDSU** genome-scale metabolic model.

Two formats are provided:

- `iZM_SDSU.xml` — SBML version for model sharing and interoperability
- `iZM_SDSU.mat` — MATLAB version for use with the COBRA Toolbox

The final model contains:

```text
Reactions:   2,286
Metabolites: 1,880
Genes:       718
```

The biomass objective reaction is:

```text
BIOMASS_core
```

The SBML model was validated by re-importing the exported model and performing flux balance analysis (FBA).

The re-imported model retained:

```text
Reactions:   2,286
Metabolites: 1,880
Genes:       718
FBA status:  optimal
Growth:      0.059677 h^-1
```

## Anaerobic Simulation Conditions

Biofuel production analyses were performed under anaerobic conditions.

Key simulation constraints included:

```text
Glucose uptake (EX_glc__D_e) = -10 mmol gDW^-1 h^-1
Oxygen uptake  (EX_o2_e)     = 0 mmol gDW^-1 h^-1
```

The biomass reaction was used as the growth objective during computational validation.

For the isobutanol OptKnock analysis, a minimum biomass requirement of **50% of the wild-type growth rate** was used to maintain viable growth while searching for production-associated knockout strategies.

## Ethanol OptKnock Analysis

The `ethanol/` directory contains the final OptKnock workflow used to investigate knockout strategies for ethanol production under anaerobic conditions.

Main analysis:

```text
OK_etoh_0714.mlx
```

The workflow includes:

- Applying anaerobic simulation constraints
- Using ethanol exchange (`EX_etoh_e`) as the target reaction
- Running OptKnock to identify knockout candidates
- Computationally validating predicted knockout strategies using FBA
- Evaluating whether ethanol production is growth-coupled

Multiple knockout strategies were predicted during the OptKnock analyses.

However, subsequent computational validation showed that **none of the tested knockout strategies achieved validated growth-coupled ethanol production under the examined anaerobic conditions**.

Therefore, no final knockout strategy was selected for ethanol production.

## Isobutanol OptKnock Analysis

The `isobutanol/` directory contains the final OptKnock and computational validation workflows for isobutanol production.

Main analyses:

```text
OK_iso_0714.mlx
OK_iso_double_validate_0714.mlx
```

`OK_iso_0714.mlx` contains the final OptKnock analysis.

`OK_iso_double_validate_0714.mlx` contains the computational validation of the selected double-knockout strategy.

### Wild-Type Model

Under the simulated anaerobic condition, wild-type growth was:

```text
WT growth = 0.059677 h^-1
```

The maximum theoretical isobutanol production was:

```text
WT maximum isobutanol flux = 10.030751 mmol gDW^-1 h^-1
```

During biomass maximization, the wild-type model produced no isobutanol:

```text
WT isobutanol flux during growth = 0 mmol gDW^-1 h^-1
```

This indicates that isobutanol production is feasible in the expanded model but is not naturally coupled to biomass production.

### Final Validated Knockout Strategy

Following exclusion of `ALCD2x` from the OptKnock search space, an additional OptKnock analysis was performed to identify alternative knockout candidates.

The **PDH–PYRDC double knockout** was the most promising computationally validated strategy.

Computational validation showed:

```text
Growth = 0.036720 h^-1
Growth relative to WT = 61.53%
Isobutanol production = 9.669147 mmol gDW^-1 h^-1
```

The mutant therefore maintained approximately **61.5% of wild-type growth** while supporting substantial isobutanol production.

These results identify the **PDH–PYRDC double knockout** as the most promising validated isobutanol strategy under the simulated anaerobic conditions.

## FluxRETAP-Inspired Analysis

The `FluxRETAP_Isobutanol_Zmobilis/` directory contains a FluxRETAP-inspired analysis used to investigate metabolic flux changes associated with increasing isobutanol production.

The workflow compares low- and high-isobutanol production states using flux variability analysis (FVA).

The analysis includes:

- Increasing the required isobutanol production level
- Performing FVA at different production states
- Comparing reaction flux ranges between low- and high-isobutanol conditions
- Identifying reactions associated with changes in isobutanol production

The analysis showed decreased flux through ethanol-associated reactions as isobutanol production increased, including decreased flux associated with `ALCD2x`.

These results suggest competition between native ethanol production and the engineered isobutanol pathway and provide complementary information to the OptKnock analysis.

## Software Requirements

The computational analyses were performed using:

- MATLAB
- COBRA Toolbox
- Gurobi Optimizer

The SBML version of iZM_SDSU can also be imported into software that supports SBML-compatible constraint-based metabolic models.

## Model Formats

| File | Description |
|---|---|
| `iZM_SDSU.xml` | Final model in SBML format |
| `iZM_SDSU.mat` | Final model in MATLAB/COBRA Toolbox format |

The SBML version was checked by re-importing the exported model and confirming that the reaction, metabolite, and gene counts were preserved and that the model retained a feasible optimal FBA solution.

## Citation

If you use the iZM_SDSU model or computational workflows from this repository, please cite the associated master's thesis.

A permanent archived version and citation information will be added with the final repository release.

## Author

**Nhi Nguyen**  
M.S. Bioinformatics and Medical Informatics  
San Diego State University
