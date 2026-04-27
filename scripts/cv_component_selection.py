from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score

OUTPUT_DIR = Path('../output/CV_results')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_excel('../data/training_set.xlsx')
data = df[['RT', 'MW', 'ACPC', 'charge', 'Helicity', 'ME100', 'ME15', 'ab ratio',
           'MIC_CA', 'MIC_CT', 'MIC_CG', 'MIC_CP', 'IC50_3T3', 'IC50_HUVEC', 'HC10', 'MIC_SA', 'MIC_EC']]

X = data[['ME100', 'ME15', 'ab ratio', 'Helicity', 'RT', 'ACPC', 'charge', 'MW']]
Y = data[['MIC_CA', 'MIC_CT', 'MIC_CG', 'MIC_CP', 'IC50_3T3', 'IC50_HUVEC', 'HC10', 'MIC_SA', 'MIC_EC']]

# ── Leave-one-out cross-validation to select number of PLS components ──────────

loo = LeaveOneOut()
X_arr, Y_arr = X.values, Y.values
n_components_range = np.arange(1, 9)
mse = []

for n_comp in n_components_range:
    fold_mse = []
    for train_idx, test_idx in loo.split(X_arr):
        pls = PLSRegression(n_components=n_comp, scale=True, max_iter=1000)
        pls.fit(X_arr[train_idx], Y_arr[train_idx])
        Y_pred = pls.predict(X_arr[test_idx])
        fold_mse.append(np.mean((Y_arr[test_idx] - Y_pred) ** 2))
    mse.append(np.mean(fold_mse))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle('Component validation')
ax1.plot(n_components_range, mse)
ax1.set_xlabel('Number of PLS Components')
ax1.set_ylabel('MSE')
ax1.set_title('LOO-CV MSE vs Components')

# ── Variance explained in X and Y per component ────────────────────────────────

n_comp = 8
pls_full = PLSRegression(n_components=n_comp, scale=True, max_iter=1000)
pls_full.fit(X, Y)

Y_arr_full = Y.values
R2_arr, R2_sum_arr = [], []
r2_sum = 0
for i in range(n_comp):
    Y_pred_s = np.dot(pls_full.x_scores_[:, i].reshape(-1, 1),
                      pls_full.y_loadings_[:, i].reshape(-1, 1).T)
    Y_pred = Y_pred_s * pls_full._y_std + pls_full._y_mean
    r2_comp = round(r2_score(Y_arr_full, Y_pred), 3)
    r2_sum += r2_comp
    R2_sum_arr.append(r2_sum)
    R2_arr.append(r2_comp)
    print(f'R² for component {i + 1}: {r2_comp:.3g}')
print(f'R² cumulative:     {r2_sum:.3g}')
print(f'R² from predict(): {round(r2_score(Y_arr_full, pls_full.predict(X)), 3):.3g}')

X_scaled = (X.values - pls_full._x_mean) / pls_full._x_std
total_var_x = np.var(X_scaled, axis=0).sum()
R2_arr_x, R2_sum_arr_x = [], []
r2_sum_x = 0
for i in range(n_comp):
    var_i = np.var(pls_full.x_scores_[:, i]) / total_var_x
    r2_sum_x += var_i
    R2_sum_arr_x.append(r2_sum_x)
    R2_arr_x.append(var_i)
    print(f'X variance by component {i + 1}: {var_i * 100:.2f}%')

n_comp_plot = list(range(1, n_comp + 1))
ax2.plot(n_comp_plot, R2_arr_x,     label='X variance per component', color='green', linestyle='dashed')
ax2.plot(n_comp_plot, R2_sum_arr_x, label='Cumulative X variance',    color='green')
ax2.plot(n_comp_plot, R2_arr,       label='Y variance per component', color='red',   linestyle='dashed')
ax2.plot(n_comp_plot, R2_sum_arr,   label='Cumulative Y variance',    color='red')
ax2.set_xlabel('Number of PLS Components')
ax2.set_ylabel('Variance Explained (R²)')
ax2.set_title('Variance vs Components')
ax2.legend()

plt.savefig(OUTPUT_DIR / 'component_validation.png', bbox_inches='tight')
plt.show()
