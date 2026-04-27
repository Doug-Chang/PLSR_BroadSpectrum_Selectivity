# Import the necessary libraries to read
# dataset and work on that
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import neutcurve

OUTPUT_DIR = Path('../output/plotting_SI/mammalIC50')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#Set some pandas display options
pd.set_option('display.float_format', '{:.3g}'.format)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 500)

#import MIC data sheets
viabilities = pd.read_excel('../data/mammalian_hillcurve.xlsx', sheet_name=None)  # sera = cell type, virus = experimental replicates
ic50_storage = {}


for (pepName,data) in viabilities.items():
    
    data['fraction infectivity'] = data['fraction infectivity'].div(100) # normalize to 1
    
    fits = neutcurve.CurveFits(data, fixtop = False) #fit for 3T3 and HUVECs, here ***sera = cell type***
    
    
    ic50_storage[pepName] = fits.fitParams()
    
    fig, axes = fits.plotSera(xlabel='concentration (μg/mL)', ylabel='Normalized Cell viability',) #plots "sera" or cell type
    fig.suptitle('Peptide ' + pepName, size = 16) # peptide number, determined from key (peptide ID) of pandas dict
    fig.savefig(OUTPUT_DIR / f'peptide_{pepName}.png', bbox_inches='tight')
    plt.close(fig)


