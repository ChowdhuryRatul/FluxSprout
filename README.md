# FluxSprout

**Single-cell metabolic modeling of vascular sprouting in zebrafish**

Mohammad Reza Zargar and Ratul Chowdhury, Department of Chemical and Biological Engineering and Nanovaccine Institute, Iowa State University

FluxSprout combines single-cell RNA sequencing of the 24 hpf zebrafish trunk with Zebrafish-GEM to build state-specific metabolic models of sprouting and established vessels. It applies four checks to every prediction:

1. A stoichiometric audit of the pathway labels.
2. An audit of thermodynamically infeasible loops.
3. Tests at both the marker and the flux level, controlled for cell-cycle phase.
4. A knockout screen across an ensemble of four extraction settings and all 13 non-endothelial lineages.

This repository contains the code, the models and every result behind the figures, tables and supplementary material of the article.

## Contents
| Path | Content |
| --- | --- |
| `FluxSprout_analysis.ipynb` | The full analysis, sections 1 to 14, from raw count matrices to every figure and table (executed, with outputs) |
| `fluxsprout/` | Python modules used by the notebook (listed below) |
| `models/` | All models in COBRApy-readable SBML: the reference network and the 108 context-specific models (see [`models/README.md`](models/README.md)) |
| `results/` | Cached results: cell annotations, metacells, context model reaction lists (`ctx_*.json`), knockout screens (`ko_*.csv`), E-Flux read-outs and result tables |
| `figures/` | Figs. 1 to 6, vector PDF and 600 dpi PNG |
| `supplementary/` | Figs. S1 to S4, Tables S1 to S5 (CSV) and Supplementary Data S1 (`FluxSprout_Supplementary_Data.xlsx`) |

Modules in `fluxsprout/`:
* `hlp.py`: HiGHS engine for FBA, FVA, pFBA, loop removal, gene knockouts, FASTCORE and GIMME.
* `pipeline.py`: context models, E-Flux and knockout screens.
* `fslib.py`: medium, objective and reaction sets.
* `ctxlib.py`: GPR evaluation and local-T2 activity calls.
* `figstyle.py`: figure style.
* `export_models.py`: SBML export.

## Using the models
No input data are needed to use the models; only COBRApy.
```python
import cobra
ref = cobra.io.read_sbml_model('models/FluxSprout_reference.xml.gz')
tip = cobra.io.read_sbml_model('models/context/fastcore_default/ISV_tip.xml.gz')
print(tip.slim_optimize())
```
Each model carries the embryonic medium, the maintenance requirement and the biomass objective used in the article. `models/example_cobrapy.py` reproduces knockout results from the screen with `cobra.flux_analysis.single_gene_deletion`.

## Reproducing the analysis
1. Install Python 3.11 or later and the pinned packages:
   ```
   pip install -r requirements.txt
   ```
2. Download the three input files into the repository root. They are not redistributed here.

   | File | Source | SHA-256 |
   | --- | --- | --- |
   | `GSM6138416_filtered_feature_bc_matrix_GFP+mCherry+.h5` | [GEO GSM6138416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6138416) | `2e60f2fd49e4ac52d33ae64e07a57033bf6c7f82aa152bb00a154cb1eaebdeb4` |
   | `GSM6138417_filtered_feature_bc_matrix_GFP+mCherry-.h5` | [GEO GSM6138417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6138417) | `3c2c4056606c70409910b544a862fbd72abec25b0807dd060f52579792aa6119` |
   | `Zebrafish-GEM.xml` | [SysBioChalmers/Zebrafish-GEM](https://github.com/SysBioChalmers/Zebrafish-GEM), `model/` folder | `3424f584f3f48711d79a0a5af475e0ad0413af10c93345a7bb76e53874a1e217` |

   Check a download with `shasum -a 256 <file>` (macOS) or `sha256sum <file>` (Linux). A different Zebrafish-GEM checksum means a different model release, and some numbers may change.
3. Run the notebook:
   ```
   jupyter lab FluxSprout_analysis.ipynb
   ```
   With `results/` present, the notebook runs in about 15 to 25 minutes. Quality control, scoring, the model audit and all plots are recomputed; the model ensemble and knockout screens are read from `results/`. Set `RECOMPUTE = True` in the first code cell to rebuild everything, which takes about 4 to 6 hours on 2 cores (the FASTCORE knockout screens dominate). Intermediate files are written to `cache/`.
4. Optionally, rewrite the SBML models: `python fluxsprout/export_models.py`. This needs `cache/` from step 3.

Cell labels are stored in `results/cell_annotations.csv`, because Leiden cluster numbering can change between library versions.

## Key settings
* **Medium:** yolk-like embryonic medium (Table S2), in mmol gDW⁻¹ h⁻¹:
  * glucose 10, O2 20, glutamine 10;
  * other amino acids 1, fatty acids 1, vitamins 0.1;
  * vitamin B12 and gamma-tocopherol, so that the biomass cofactor pool can be made.
* **Maintenance:** NGAM (MAR03964) at 1 mmol gDW⁻¹ h⁻¹.
* **Objective:** MAR13082, generic human biomass with 45 mmol ATP gDW⁻¹ growth-associated maintenance. MAR00021 and MAR04413 are used in the audit.
* **E-Flux:** bound = ŝ × U, with ŝ normalized to its 95th percentile and ε = 0.05. U = 50 is the primary setting and U = 100 the sensitivity setting.
* **Ensemble:** FASTCORE with local-T2 thresholds (15/65, 25/75 and 40/85 percentiles) and GIMME, over 27 groups (14 endothelial and hemogenic states, 13 non-endothelial lineages).
* **Calling rules:**
  * Two-lineage rule: sprouting inhibition ≥ 0.10, with erythroid and muscle spared.
  * All-lineage rule: all 13 non-endothelial lineages spared, and within-endothelium selectivity ≥ 0.10.

## Software
The reference run used Python 3.13 with these packages:
* COBRApy 0.32.1, highspy 1.15.1, python-libsbml 5.21.2;
* scanpy 1.12.4, anndata 0.13.4, leidenalg 0.12.0, igraph 1.0.0;
* statsmodels 0.15.0, scikit-learn 1.9.1, numpy 2.5.3, pandas 2.3.3, scipy 1.18.1;
* matplotlib 3.11.2.

The HiGHS engine in `hlp.py` is checked against COBRApy (FBA, FVA, pFBA, gene knockouts and blocked reactions) in section 1 of the notebook.

## Citation
If you use FluxSprout, please cite the article (citation will be added on publication) and the resources it builds on:
* Zebrafish-GEM: H. Wang *et al.*, *Proc. Natl. Acad. Sci. U.S.A.* 118, e2102344118 (2021).
* The single-cell data: GEO GSM6138416 and GSM6138417.

## Contact
Ratul Chowdhury (ratul@iastate.edu), Chowdhury Lab, Iowa State University.
