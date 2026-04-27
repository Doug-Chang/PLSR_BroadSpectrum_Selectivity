# Import the necessary libraries to read
# dataset and work on that
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

#import MIC data sheets
df_3T3 = pd.read_excel('Compiled list of percent cell viability broad spectrum 3T3-IC50ver.xlsx', header=None, sheet_name = None)
df_HUVEC = pd.read_excel('Compiled list of percent cell viability broad spectrum HUVEC-IC50 ver.xlsx', header=None, sheet_name = None)

#Below function removes other info and saves only the MIC values for plotting

def cleanData (df):
    df = df.iloc[:,15:27]
    df = df.dropna(how='all')
#    df = df.dropna(axis = 1) #remove nan overlappers
    return df

#iterate through each spread sheet (peptide ID) and plot cell viability curves for all peptides
#k3T3..etc = keys of dataframe dictionary = peptide ID
#v3T3..etc = values of dataframe dictionary = matrix of peptide concentrations, cell viability values, and standard deviation
for (k3T3,v3T3), (kHUVEC,vHUVEC),  in zip(df_3T3.items(), df_HUVEC.items()):
    
    v3T3 = cleanData(v3T3)
    vHUVEC = cleanData(vHUVEC)
    
    #plotting formatting 
    fig, axs = plt.subplots(ncols = 2,figsize=(11, 5))
        
    axs[0].title.set_text('3T3 Cell viability')
    axs[1].title.set_text('HUVEC Cell viability')    
    fig.suptitle('Peptide ' + k3T3, size = 16) # peptide number, determined from key (peptide ID) of pandas dict
   
    count = 1
    for i in range(1, v3T3.shape[0]-1,2): #plot all experimental curves
        axs[0].errorbar(v3T3.iloc[0,:],v3T3.iloc[i,:], v3T3.iloc[i+1,:], marker = 'o', label = f"3T3_{count}",capsize=2)
        count+=1
    
    count = 1
    for i in range(1, vHUVEC.shape[0]-1,2): #plot all experimental curves
        axs[1].errorbar(vHUVEC.iloc[0,:],vHUVEC.iloc[i,:], vHUVEC.iloc[i+1,:], marker = 'o',label = f"HUVEC_{count}",capsize=2)
        count+=1
    
    #formatting common things
    for ax in axs.flat:
        ax.set_xscale('log')
        ax.legend()
        ax.xaxis.set_major_formatter(ScalarFormatter())
        #ax.locator_params(axis='x', numticks=16)
        ax.set_box_aspect(1)
        ax.set_xlabel('Concentration (μg/mL)')
        ax.set_xlim(0.1,1200)
        ax.set_ylim(-20,150)
        ax.axhline(y=0, color='gray', linestyle='-')
        
    axs[0].set_ylabel('Cell viability (%)')





