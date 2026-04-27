from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

OUTPUT_DIR = Path('../output/PCA_plot')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_excel('../data/prediction_set.xlsx')
X = df[['ME100', 'ME15', 'ab ratio', 'Helicity', 'RT', 'ACPC', 'charge', 'MW']]

sns.countplot(x='Set', data=df)
plt.title('Training vs. prediction set size')
plt.savefig(OUTPUT_DIR / 'training_vs_prediction_count.png', bbox_inches='tight')
plt.close()

x_scaled = StandardScaler().fit_transform(X)
pca_features = PCA(n_components=4).fit_transform(x_scaled)

pca_df = pd.DataFrame(pca_features, columns=['PC1', 'PC2', 'PC3', 'PC4'])
pca_df['target'] = df['Set'].values
pca_df = pca_df.iloc[::-1].reset_index(drop=True)

sns.set_style("ticks")
palette = {
    'Training ': 'tab:orange',
    'Prediction': 'tab:blue',
}

g = sns.lmplot(x='PC1', y='PC2', data=pca_df, hue='target', palette=palette,
               fit_reg=False, legend=True)
g.axes[0, 0].invert_yaxis()
plt.title('2D PCA')
plt.savefig(OUTPUT_DIR / 'pca_2d.png', bbox_inches='tight')
plt.close()
