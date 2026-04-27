from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

OUTPUT_DIR = Path('../output/MIC_plots')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df_CA = pd.read_excel('../data/viability_ca.xlsx', header=None, sheet_name=None)
df_CG = pd.read_excel('../data/viability_cg.xlsx', header=None, sheet_name=None)
df_CP = pd.read_excel('../data/viability_cp.xlsx', header=None, sheet_name=None)
df_CT = pd.read_excel('../data/viability_ct.xlsx', header=None, sheet_name=None)
df_EC = pd.read_excel('../data/viability_ec.xlsx', header=None, sheet_name=None)
df_SA = pd.read_excel('../data/viability_sa.xlsx', header=None, sheet_name=None)


def clean_data(df):
    return df.iloc[0:3, 15:27].dropna(axis=1)


for (kCA, vCA), (kCG, vCG), (kCP, vCP), (kCT, vCT), (kEC, vEC), (kSA, vSA) in zip(
        df_CA.items(), df_CG.items(), df_CP.items(), df_CT.items(), df_EC.items(), df_SA.items()):

    vCA, vCG, vCP, vCT, vEC, vSA = [clean_data(v) for v in (vCA, vCG, vCP, vCT, vEC, vSA)]

    fig, axs = plt.subplots(ncols=2, figsize=(11, 5))
    axs[0].set_title('Fungal Cell viability')
    axs[1].set_title('Bacterial Cell viability')
    fig.suptitle(f'Peptide {kCA}', size=16)

    axs[0].errorbar(vCA.iloc[0], vCA.iloc[1], vCA.iloc[2], marker='o', label='C. albicans',    capsize=2)
    axs[0].errorbar(vCG.iloc[0], vCG.iloc[1], vCG.iloc[2], marker='o', label='C. glabrata',    capsize=2)
    axs[0].errorbar(vCP.iloc[0], vCP.iloc[1], vCP.iloc[2], marker='o', label='C. parapsilosis', capsize=2)
    axs[0].errorbar(vCT.iloc[0], vCT.iloc[1], vCT.iloc[2], marker='o', label='C. tropicalis',  capsize=2)
    axs[1].errorbar(vEC.iloc[0], vEC.iloc[1], vEC.iloc[2], marker='o', label='E. coli',         capsize=2)
    axs[1].errorbar(vSA.iloc[0], vSA.iloc[1], vSA.iloc[2], marker='o', label='S. aureus',      capsize=2)

    for ax in axs.flat:
        ax.set_xscale('log')
        ax.legend()
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_box_aspect(1)
        ax.set_xlabel('Concentration (μg/mL)')
        ax.set_xlim(0.125, 1300)
        ax.set_ylim(-20, 200)
        ax.axhline(y=0, color='gray', linestyle='-')
    axs[0].set_ylabel('Cell viability (%)')

    fig.savefig(OUTPUT_DIR / f'peptide_{kCA}.png', bbox_inches='tight')
    plt.close(fig)
