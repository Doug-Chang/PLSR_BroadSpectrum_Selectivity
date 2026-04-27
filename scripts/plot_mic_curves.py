# Import the necessary libraries to read
# dataset and work on that
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

OUTPUT_DIR = Path('../output/MIC_plots')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#import MIC data sheets
df_3T3 = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum 3T3.xlsx', header=None, sheet_name = None)
df_HUVEC = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum HUVEC.xlsx', header=None, sheet_name = None)
df_CA = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum ca.xlsx', header=None, sheet_name = None)
df_CG = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum cg.xlsx', header=None, sheet_name = None)
df_CP = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum cp.xlsx', header=None, sheet_name = None)
df_CT = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum ct.xlsx', header=None, sheet_name = None)
df_EC = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum ec.xlsx', header=None, sheet_name = None)
df_SA = pd.read_excel('../data/Compiled list of percent cell viability broad spectrum sa.xlsx', header=None, sheet_name = None)

#Below function removes other info and saves only the MIC values for plotting
def cleanData (df):
    df = df.iloc[0:3,15:27]
    df = df.dropna(axis = 1) #remove 
    return df

#this is to plot all non-mammalian data
for (kCA,vCA),(kCG,vCG),(kCP,vCP),(kCT,vCT),(kEC,vEC),(kSA,vSA) in zip(df_CA.items(), df_CG.items(),df_CP.items(), df_CT.items(),df_EC.items(),df_SA.items()):
    
    vCA = cleanData(vCA)
    vCG = cleanData(vCG)
    vCP = cleanData(vCP)
    vCT = cleanData(vCT)
    vEC = cleanData(vEC)
    vSA = cleanData(vSA)
    
    #plotting formatting 
    fig, axs = plt.subplots(ncols = 2,figsize=(11, 5))
        
    axs[0].title.set_text('Fungal Cell viability')
    axs[1].title.set_text('Bacterial Cell viability')
    fig.suptitle('Peptide ' + kCA, size = 16) # peptide number, determined from key (peptide ID) of pandas dict
    
    axs[0].errorbar(vCA.iloc[0,:],vCA.iloc[1,:], vCA.iloc[2,:], marker = 'o', label = "C. albicans",capsize=2)
    axs[0].errorbar(vCG.iloc[0,:],vCG.iloc[1,:], vCG.iloc[2,:], marker = 'o',label = "C. glabrata",capsize=2)
    axs[0].errorbar(vCP.iloc[0,:],vCP.iloc[1,:], vCP.iloc[2,:], marker = 'o',label = "C. parapsilosis",capsize=2)
    axs[0].errorbar(vCT.iloc[0,:],vCT.iloc[1,:], vCT.iloc[2,:], marker = 'o',label = "C. tropicalis",capsize=2)
    
    axs[1].errorbar(vEC.iloc[0,:],vEC.iloc[1,:], vEC.iloc[2,:], marker = 'o',label = "E. coli",capsize=2)
    axs[1].errorbar(vSA.iloc[0,:],vSA.iloc[1,:], vSA.iloc[2,:], marker = 'o',label = "S. aureus",capsize=2)
    
    #formatting common things
    for ax in axs.flat:
        ax.set_xscale('log')
        ax.legend()
        ax.xaxis.set_major_formatter(ScalarFormatter())
        #ax.locator_params(axis='x', numticks=16)
        ax.set_box_aspect(1)
        ax.set_xlabel('Concentration (μg/mL)')
        ax.set_xlim(0.125,1300)
        ax.set_ylim(-20,200)
        ax.axhline(y=0, color='gray', linestyle='-')
        
    axs[0].set_ylabel('Cell viability (%)')
    fig.savefig(OUTPUT_DIR / f'peptide_{kCA}.png', bbox_inches='tight')
    plt.close(fig)

#to plot mammalian cell averaged viability curve data:
'''
#iterate through each spread sheet (peptide ID) and plot cell viability curves for all peptides
#k3T3..etc = keys of dataframe dictionary = peptide ID
#v3T3..etc = values of dataframe dictionary = matrix of peptide concentrations, cell viability values, and standard deviation
for (k3T3,v3T3), (kHUVEC,vHUVEC), (kCA,vCA),(kCG,vCG),(kCP,vCP),(kCT,vCT),(kEC,vEC),(kSA,vSA) in zip(df_3T3.items(), df_HUVEC.items(), df_CA.items(), df_CG.items(),df_CP.items(), df_CT.items(),df_EC.items(),df_SA.items()):
    
    v3T3 = cleanData(v3T3)
    vHUVEC = cleanData(vHUVEC)
    vCA = cleanData(vCA)
    vCG = cleanData(vCG)
    vCP = cleanData(vCP)
    vCT = cleanData(vCT)
    vEC = cleanData(vEC)
    vSA = cleanData(vSA)
    
    #plotting formatting 
    fig, axs = plt.subplots(ncols = 3,figsize=(15, 5))
        
    axs[0].title.set_text('Fungal Cell viability')
    axs[1].title.set_text('Bacterial Cell viability')
    axs[2].title.set_text('Mammalian Cell viability')
    fig.suptitle('Peptide ' + kCA, size = 16) # peptide number, determined from key (peptide ID) of pandas dict
    
    axs[0].errorbar(vCA.iloc[0,:],vCA.iloc[1,:], vCA.iloc[2,:], marker = 'o', label = "C. albicans",capsize=2)
    axs[0].errorbar(vCG.iloc[0,:],vCG.iloc[1,:], vCG.iloc[2,:], marker = 'o',label = "C. glabrata",capsize=2)
    axs[0].errorbar(vCP.iloc[0,:],vCP.iloc[1,:], vCP.iloc[2,:], marker = 'o',label = "C. parapsilosis",capsize=2)
    axs[0].errorbar(vCT.iloc[0,:],vCT.iloc[1,:], vCT.iloc[2,:], marker = 'o',label = "C. tropicalis",capsize=2)
    
    axs[1].errorbar(vEC.iloc[0,:],vEC.iloc[1,:], vEC.iloc[2,:], marker = 'o',label = "E. coli",capsize=2)
    axs[1].errorbar(vSA.iloc[0,:],vSA.iloc[1,:], vSA.iloc[2,:], marker = 'o',label = "S. aureus",capsize=2)
    
    axs[2].errorbar(v3T3.iloc[0,:],v3T3.iloc[1,:], v3T3.iloc[2,:], marker = 'o', label = "3T3",capsize=2)
    axs[2].errorbar(vHUVEC.iloc[0,:],vHUVEC.iloc[1,:], vHUVEC.iloc[2,:], marker = 'o',label = "HUVECs",capsize=2)
    
    #formatting common things
    for ax in axs.flat:
        ax.set_xscale('log')
        ax.legend()
        ax.xaxis.set_major_formatter(ScalarFormatter())
        #ax.locator_params(axis='x', numticks=16)
        ax.set_box_aspect(1)
        ax.set_xlabel('Concentration (μg/mL)')
        ax.set_xlim(0.5,1100)
        ax.set_ylim(-50,200)
        ax.axhline(y=0, color='gray', linestyle='-')
        
    axs[0].set_ylabel('Cell viability (%)')
'''

