import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns
from sklearn import linear_model
import math
import matplotlib.pyplot as plt
from statsmodels.stats.outliers_influence import variance_inflation_factor

labels = [ 'E. coli','S. aureus','C. albicans', 'C. tropicalis', 'C. parapsilosis', 'C. glabrata', '3T3 Fibroblasts', 'HUVECs','Red Blood Cells']
outputVariables = ['MIC_EC','MIC_SA','MIC_CA', 'MIC_CT','MIC_CP','MIC_CG','IC50_3T3','IC50_HUVEC','HC10']

df = pd.read_excel('../data/PearsonSpecies.xlsx')
#df = df.truncate(after = 17) 

#check indipendent variable colinearity

X = df[outputVariables]
X.columns = labels

# Correlation between different variables
corr = X.corr()

# Set up the matplotlib plot configuration
f, ax = plt.subplots(figsize=(12, 10))

# Generate a mask for upper traingle
mask = np.triu(np.ones_like(corr, dtype=bool))

# Configure a custom diverging colormap

cmap = sns.color_palette("coolwarm", as_cmap=True)
sns.set(font_scale=3, rc={'axes.facecolor':'white', 'figure.facecolor':'white'})

# Draw the heatmap
sns.heatmap(corr, annot=False, mask=(mask), cmap=cmap,vmin=-1, vmax=1,  cbar_kws={'label': 'Pearson correlation coefficient'})

#_ = sns.pairplot(X, kind="reg", diag_kind="kde", corner = True)