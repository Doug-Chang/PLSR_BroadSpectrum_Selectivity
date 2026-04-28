# Broad-Spectrum Antimicrobial α/β-Peptide Selectivity Prediction via PLSR and understanding design rules

This repository contains the analysis pipeline for the paper:

> Chang, D. et al. Establishing Quantifiable Guidelines for Antimicrobial α/β-Peptide Design: A Partial Least-Squares Approach to Improve Antimicrobial Activity and Reduce Mammalian Cell Toxicity. *ACS Infectious Diseases* (2024). https://doi.org/10.1021/acsinfecdis.3c00468

## Overview
- Predict antimicrobial peptide selectivity using PLSR
- Input: 8 peptide descriptors → Output: MIC + toxicity across species
- Reproduces figures and analysis from Chang et al. (2024)

## Background

Antimicrobial resistance is driving demand for peptide-based therapeutics with activity across divergent pathogens. α/β-Peptides — sequences that incorporate cyclic β-amino acids (here, ACPC: α-aminocyclopentane carboxylic acid) alongside standard α-residues — offer improved protease resistance and tunable secondary structure compared to all-α AMPs, but the high-dimensional relationship between their physicochemical properties and biological activity makes rational design difficult.

This work applies **Partial Least Squares Regression (PLSR)** to map eight peptide descriptors (helicity, hydrophobicity, charge, α/β ratio, etc.) to minimum inhibitory concentrations (MICs) and cytotoxicity values across nine microbial and mammalian targets simultaneously. The resulting 3-component model reveals which structural features drive broad-spectrum antimicrobial activity vs. mammalian cell toxicity, enabling quantitative selectivity prediction on a novel validation set.

## Quick start

```bash
git clone https://github.com/Doug-Chang/PLSR_BroadSpectrum_Selectivity.git
cd PLSR_BroadSpectrum_Selectivity
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd scripts
python pls_analysis.py
```

## Repository structure
data/           — Raw experimental data (Excel)
scripts/        — Analysis scripts (run from the scripts/ directory)
output/         — Generated figures (gitignored, created on first run)
requirements.txt
```

## Analysis pipeline

**All scripts must be run from inside the `scripts/` directory.** They use relative paths (`../data/`, `../output/`) that resolve correctly only from there. Running from the repo root will raise a `FileNotFoundError`.

```bash
cd scripts
python cv_component_selection.py
```

Run core model scripts in order:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `cv_component_selection.py` | LOO cross-validation over 1–8 components; plots MSE and X/Y variance explained to justify 3-component model |
| 2 | `pca_analysis.py` | PCA of training vs. prediction sets to confirm applicability domain overlap & diversity check |
| 3 | `pls_analysis.py` | Fit 3-component PLSR; generate loadings/scores biplots (2D & 3D), prediction accuracy plots, VIP scores, permutation importance, and selectivity heatmaps |

Supplementary / visualization scripts (order-independent):

| Script | Purpose |
|--------|---------|
| `species_correlation.py` | Pearson correlation heatmap across species MIC values — motivates multivariate over single-species modeling |
| `plot_mic_curves_BroadSpec.py` | Dose-response curves for all 6 pathogens (4 fungal, 2 bacterial) per training-set peptide |
| `plot_mic_curves_full_c_albi.py` | Full *C. albicans* MIC curves across the complete peptide library; groups of 5 peptides per figure, saved to `output/CA_mic_plots/` |
| `plot_ic50_curves.py` | Hill curve fitting for mammalian cytotoxicity (3T3 fibroblasts, HUVECs) using `neutcurve`; saved to `output/mammalIC50/` |
| `plot_HC10_curves.py` | Hemolysis dose-response curves (HC10) across the peptide library; groups of 5 peptides per figure, saved to `output/HC10_plots/` |

## Peptide descriptors

| Descriptor | Description |
|------------|-------------|
| `heli_100` / `heli_15` | Mean residue ellipticity at 100% / 15% TFE by CD (helical propensity in membrane-mimetic conditions) |
| `ab ratio` | α/β amino acid ratio |
| `% Helicity` | Helical content in aqueous buffer by CD |
| `RT` | Reversed-phase HPLC retention time (hydrophobicity proxy) |
| `ACPC` | Fraction of β-residues (α-aminocyclopentane carboxylic acid) |
| `charge` | Net charge at pH 7 |
| `MW` | Molecular weight (Da) |

## Output variables

| Variable | Description |
|----------|-------------|
| `MIC_CA`, `MIC_CT`, `MIC_CP`, `MIC_CG` | log₂ MIC (μg/mL) vs. *C. albicans*, *C. tropicalis*, *C. parapsilosis*, *C. glabrata* |
| `MIC_EC`, `MIC_SA` | log₂ MIC (μg/mL) vs. *E. coli*, *S. aureus* |
| `IC50_3T3`, `IC50_HUVEC` | log₂ IC50 (μg/mL) vs. 3T3 fibroblasts and HUVECs |
| `HC10` | log₂ concentration (μg/mL) at 10% hemolysis |

## Requirements

```bash
pip install -r requirements.txt
```

Key dependencies: `scikit-learn`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `neutcurve`
