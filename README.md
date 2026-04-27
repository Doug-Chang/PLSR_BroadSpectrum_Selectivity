# Broad-Spectrum Antimicrobial Peptide Selectivity via PLSR

This repository contains the analysis pipeline for a **Partial Least Squares Regression (PLSR)** model that predicts and rationalizes the broad-spectrum selectivity of antimicrobial peptides (AMPs) against bacterial, fungal, and mammalian cell targets.

## Background

Broad-spectrum AMPs must kill pathogens across divergent phylogenetic clades while sparing host cells — a multidimensional optimization problem that resists intuitive design. This project uses PLSR to map peptide physicochemical descriptors (helicity, hydrophobicity, charge, etc.) to minimum inhibitory concentrations (MICs) and cytotoxicity values across nine species, enabling rational selectivity prediction.

## Repository structure

```
data/           — Raw experimental data (Excel)
scripts/        — Analysis scripts (run from the scripts/ directory)
output/         — Generated figures (gitignored, created on first run)
requirements.txt
```

## Analysis pipeline

Run scripts in order from the `scripts/` directory:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `cv_component_selection.py` | LOO cross-validation to choose number of PLS components; plots MSE and variance explained per component |
| 2 | `pca_analysis.py` | PCA of training vs. prediction sets to verify applicability domain |
| 3 | `pls_analysis.py` | Train 3-component PLSR model; generate loadings/scores biplots, prediction accuracy plots, VIP scores, and permutation importance |
| — | `species_correlation.py` | Pearson correlation heatmap between species MIC values (motivates multivariate modeling) |
| — | `plot_mic_curves.py` | Dose-response curves for fungal and bacterial species per peptide |
| — | `fit_hill_curves.py` | Hill curve fitting for mammalian IC50 values (3T3 and HUVEC) |

## Peptide descriptors

| Descriptor | Description |
|------------|-------------|
| `heli_100` / `heli_15` | Mean residue ellipticity at 100% / 15% TFE (helical propensity) |
| `ab ratio` | Amphipathicity ratio |
| `% Helicity` | Helical content by CD |
| `RT` | Reversed-phase HPLC retention time (hydrophobicity) |
| `ACPC` | α-Aminocyclopropa­ne carboxylic acid content |
| `charge` | Net charge at pH 7 |
| `MW` | Molecular weight |

## Output variables

MIC (log₂ μg/mL) against *E. coli*, *S. aureus*, *C. albicans*, *C. tropicalis*, *C. parapsilosis*, *C. glabrata*; IC50 (log₂ μg/mL) against 3T3 fibroblasts and HUVECs; HC10 (log₂ μg/mL) hemolysis.

## Citation

Chang, D. et al. Establishing Quantifiable Guidelines for Antimicrobial α/β-Peptide Design: A Partial Least-Squares Approach to Improve Antimicrobial Activity and Reduce Mammalian Cell Toxicity. *ACS Infectious Diseases* (2024). https://doi.org/10.1021/acsinfecdis.3c00468

## Requirements

```bash
pip install -r requirements.txt
```

Key dependencies: `scikit-learn`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `neutcurve`
