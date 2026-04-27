# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 13:53:23 2023

@author: Doug
PCA plot with labeled training and prediction

"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import scale 
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns

OUTPUT_DIR = Path('../output/PCA_plot')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_excel('../data/prevDataPepMR5.xlsx')
X = df[['ME100','ME15','ab ratio','Helicity','RT', 'ACPC','charge','MW']]
Y = df[['Set']]

target_names = {'Training','Prediction'}

sns.countplot(
    x='Set', 
    data=df)
plt.title('Training v.s. prediction size')
plt.savefig(OUTPUT_DIR / 'training_vs_prediction_count.png', bbox_inches='tight')
plt.close()

x_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=4) #determined from crossvalidation
 
pca_features = pca.fit_transform(x_scaled)
 
print('Shape before PCA: ', x_scaled.shape)
print('Shape after PCA: ', pca_features.shape)
 
pca_df = pd.DataFrame(
    data=pca_features, 
    columns=['PC1', 'PC2', 'PC3','PC4'])
pca_df['target'] = Y
pca_df = pca_df.reindex(index=pca_df.index[::-1])

import matplotlib.pyplot as plt 
import seaborn as sns
sns.set_style("ticks")
#setting colors
palette = {
    'Training ': 'tab:orange',
    'Prediction': 'tab:blue',
}

#plotting pca
g = sns.lmplot(
    x='PC1',
    y='PC2',
    data=pca_df,
    hue='target',
    palette =palette,
    fit_reg=False,
    legend=True
    )
g.axes[0, 0].invert_yaxis()
plt.title('2D PCA Graph')
plt.savefig(OUTPUT_DIR / 'pca_2d.png', bbox_inches='tight')
plt.close()