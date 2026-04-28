from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

OUTPUT_DIR = Path('../output/CA_mic_plots')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data from a CSV file
df = pd.read_excel('../data/c_albicans_data_full.xlsx')
df['Peptide'] = df['Peptide'].astype(str).str.zfill(3).apply(lambda x: '#' + x if x[0] != '#' else x)  # convert all peptide names to #01..etc
plt.rcParams["figure.dpi"] = 500

# Group data by peptide
groups = df.groupby('Peptide')

# Define a dictionary of marker shapes for each peptide
markers = ['o', 's', '^',  'd', '*']

# Set up counter variable
plot_counter = 0
marker_counter = 0 #set up marker counter
# Loop through groups
for peptide, data in groups:
    # Check if a new plot needs to be created
    if plot_counter % 5 == 0:
        fig, ax = plt.subplots()
        
    # Plot data on current axes
    ax.errorbar(data['con'], data['Viability'], yerr=data['sd'], label=peptide, marker = markers[marker_counter], capsize=2)
    
    # Increment counter
    plot_counter += 1
    marker_counter +=1
    
    # Set axis labels and title on current axes
    ax.set_xlabel('Concentration (μg/mL)',fontsize=15)
    ax.set_ylabel('Cell Viability (%)',fontsize=15)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlim(0.1,1024)
    # Show the legend on the last plot only
    if plot_counter % 5 == 0 or plot_counter == len(groups):
        ax.legend()
    
        plt.axhline(y=0, color='gray', linestyle='-')
        batch = plot_counter // 5
        fig.savefig(OUTPUT_DIR / f'batch_{batch:02d}.png', bbox_inches='tight')
        plt.close(fig)
        marker_counter = 0 #reset counter