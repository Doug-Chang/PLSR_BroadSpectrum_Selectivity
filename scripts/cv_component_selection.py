from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score
import seaborn as sns


OUTPUT_DIR = Path('../output/CV_results')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#FOR BROAD SPECTRUM
#import data and use only ones I need

df = pd.read_excel('../data/simpPepGeoMean.xlsx')

newData = df[['RT', 'MW','ACPC','charge','Helicity','ME100','ME15','ab ratio','MIC_CA', 'MIC_CT','MIC_CG','MIC_CP','IC50_3T3','IC50_HUVEC','HC10','MIC_SA','MIC_EC']]

#split to x and y
Y = newData[['MIC_CA', 'MIC_CT','MIC_CG','MIC_CP','IC50_3T3','IC50_HUVEC','HC10','MIC_SA','MIC_EC']]
X = newData[['ME100','ME15','ab ratio','Helicity','RT', 'ACPC','charge','MW']]

#Correlation matrix (PLS motivation)========================================================================
# # Correlation between different variables
# corr = X.corr()

# # Set up the matplotlib plot configuration
# f, ax = plt.subplots(figsize=(12, 10))

# # Generate a mask for upper traingle
# mask = np.triu(np.ones_like(corr, dtype=bool))

# # Configure a custom diverging colormap
# cmap = sns.color_palette("coolwarm", as_cmap=True)

# # Draw the heatmap
# sns.heatmap(corr, annot=True, mask = mask, cmap=cmap,vmin=-1, vmax=1,  cbar_kws={'label': 'Pearson correlation coefficient'})
# plt.savefig(OUTPUT_DIR / 'x_correlation_heatmap.png', bbox_inches='tight')

# _ = sns.pairplot(X, kind="reg", diag_kind="kde", corner = True)
# _.savefig(OUTPUT_DIR / 'x_pairplot.png', bbox_inches='tight')

#MSE CV ============================================================================================================
loo = LeaveOneOut()
X_arr = X.values
Y_arr = Y.values
n_components_range = np.arange(1, 9)
mse = []

for n_comp in n_components_range:
    fold_mse = []
    for train_idx, test_idx in loo.split(X_arr):
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        Y_train, Y_test = Y_arr[train_idx], Y_arr[test_idx]
        pls = PLSRegression(n_components=n_comp, scale=True, max_iter=1000)
        pls.fit(X_train, Y_train)
        Y_pred = pls.predict(X_test)
        fold_mse.append(np.mean((Y_test - Y_pred) ** 2))
    mse.append(np.mean(fold_mse))

#plot test MSE vs. number of components
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle('Component validation')

ax1.plot(n_components_range, mse)
ax1.set_xlabel('Number of PLS Components')
ax1.set_ylabel('MSE')
ax1.set_title('MSE vs Components')    

#PLS variance explained=========================================================================================

# X is a numpy ndarray with samples in rows and predictor variables in columns
# y is one-dimensional ndarray containing the response variable

n_comp = 8
my_plsr = PLSRegression(n_components=n_comp, scale=True, max_iter=1000)
my_plsr.fit(X, Y)
R2_arr = []
R2_sum_arr = []
r2_sum = 0
n_comp_plot = [*range(1,n_comp+1)]

Y_arr_full = Y.values
for i in range(0,n_comp):
        Y_pred_s = np.dot(my_plsr.x_scores_[:,i].reshape(-1,1), my_plsr.y_loadings_[:,i].reshape(-1,1).T)
        Y_pred = Y_pred_s * my_plsr._y_std + my_plsr._y_mean
        r2_comp = round(r2_score(Y_arr_full, Y_pred),3)
        r2_sum +=  r2_comp
        R2_sum_arr.append(r2_sum)
        print('R2 for %d component: %g' %(i+1,r2_comp))
        R2_arr.append(r2_comp)
print('R2 for all components (): %g' %r2_sum) #Sum of above
print('R2 for all components (): %g' %round(r2_score(Y_arr_full, my_plsr.predict(X)),3)) #Calcuted from PLSRegression's 'predict' function.


#for x variance
# Calculate variance in X explained by PLS model
variance_explained = np.zeros(my_plsr.n_components)
X_Scaled = (X.values - my_plsr._x_mean) / my_plsr._x_std
total_variance_in_x = np.var(X_Scaled, axis = 0)
R2_arr_x = []
R2_sum_arr_x = []
r2_sum_x = 0
for i in range(my_plsr.n_components):
    x_scores_i = my_plsr.x_scores_[:, i]
    variance_explained[i] = np.var(x_scores_i) / sum(total_variance_in_x)
    r2_sum_x +=  variance_explained[i]
    R2_sum_arr_x.append(r2_sum_x)
    R2_arr_x.append(variance_explained[i])
print("Variance in X explained by PLS model:")
for i, v in enumerate(variance_explained):
    print(f"Component {i+1}: {v*100:.2f}%")

#plot test MSE vs. number of components
ax2.plot(n_comp_plot,R2_arr_x, label ='X Variance per component', color = 'green', linestyle='dashed')
ax2.plot(n_comp_plot,R2_sum_arr_x, label = 'Cumulative X variance explained',  color = 'green')
ax2.plot(n_comp_plot,R2_arr, label ='Y Variance per component', color = 'red', linestyle='dashed')
ax2.plot(n_comp_plot,R2_sum_arr, label = 'Cumulative Y variance explained',color = 'red')
ax2.set_xlabel('Number of PLS Components')
ax2.set_ylabel('Variance Explained (R2)')
ax2.set_title('Variance vs Components')
ax2.legend()
plt.savefig(OUTPUT_DIR / 'component_validation.png', bbox_inches='tight')