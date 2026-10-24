# FluxSprout models

All models are standard SBML Level 3 with the FBC package (gzipped), written and read back with COBRApy 0.32.1. COBRApy reads `.xml.gz` files directly; other tools (for example the COBRA Toolbox `readCbModel`) read them after `gunzip`.

```python
import cobra
ref = cobra.io.read_sbml_model('models/FluxSprout_reference.xml.gz')
tip = cobra.io.read_sbml_model('models/context/fastcore_default/ISV_tip.xml.gz')
tip.slim_optimize()
```

`python models/example_cobrapy.py` loads models and reproduces knockout results from the screen with `cobra.flux_analysis.single_gene_deletion`.

## Files
| File | Content |
| --- | --- |
| `FluxSprout_reference.xml.gz`, `FluxSprout_reference.json` | Reference network: Zebrafish-GEM restricted to the reactions that can carry flux under the embryonic medium (9,480 reactions, 5,578 metabolites, 2,170 genes). Same model in SBML and COBRApy JSON. |
| `context/<member>/<group>.xml.gz` | 108 context models: 4 ensemble members (`fastcore_loose`, `fastcore_default`, `fastcore_strict`, `gimme_default`) × 27 groups (14 endothelial and hemogenic states, 13 non-endothelial lineages). Each is the reference network with the reactions outside that model removed. |
| `model_index.csv` | One row per file: reactions, metabolites, genes, growth computed by COBRApy after reading the file, growth in the knockout screen, and the difference. |
| `example_cobrapy.py` | Loading, knockout and medium-change example. |

Group names follow the article; in file names `/` becomes `-`, spaces become `_` and dots are dropped (for example `LPM/cardiac` → `LPM-cardiac`, `Neural prog.` → `Neural_prog`).

## Settings stored in every model
* Objective: MAR13082 (generic human cell biomass, 45 mmol ATP gDW⁻¹ growth-associated maintenance).
* Non-growth ATP maintenance: MAR03964 lower bound 1 mmol gDW⁻¹ h⁻¹.
* Medium: exchange lower bounds set to the embryonic medium (Table S2; glucose 10, O2 20, glutamine 10, other amino acids 1, fatty acids 1, vitamins 0.1, ions 1000 mmol gDW⁻¹ h⁻¹). All other uptakes are closed, all secretions open. To change the medium, set the lower bound of the exchange reaction (negative = uptake), as in `example_cobrapy.py` and Table S2.
* The OXPHOS and PPase reaction sets used in the model audit are in `fluxsprout/fslib.py` (`OXPHOS6`, `PPI_LEGACY`).

The E-Flux metacell models are the reference network with expression-scaled bounds; they are built in section 8 of the notebook (`eflux_global` in `fluxsprout/pipeline.py`) and their flux read-outs are in `results/eflux_metacells.csv`.

## Rebuilding
`python fluxsprout/export_models.py` rewrites every file in this folder from `results/ctx_*.json` and the cached reference network, then checks that the growth COBRApy finds in each file equals the growth in the knockout screen.

The source model, Zebrafish-GEM, is from Wang et al. (SysBioChalmers); cite it with FluxSprout and keep its license terms with any redistribution.
