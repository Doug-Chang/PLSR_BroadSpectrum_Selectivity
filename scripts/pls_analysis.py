import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score
from sklearn.inspection import permutation_importance
from scipy.stats import pearsonr
import seaborn as sns

# ── Constants ─────────────────────────────────────────────────────────────────

VARIANCE_EXPLAINED = {
    "xR2": [15.0, 32.8, 22.5],  # % variance in X per component
    "yR2": [40.0, 21.5,  9.5],  # % variance in Y per component
}

DESCRIPTORS   = ['heli_100', 'heli_15', 'ab ratio', '% Helicity', 'RT', 'ACPC', 'charge', 'MW']
OUTPUT_VARS   = ['MIC_CA', 'MIC_CT', 'MIC_CG', 'MIC_CP', 'IC50_3T3', 'IC50_HUVEC', 'HC10', 'MIC_SA', 'MIC_EC']
SPECIES_LABELS = ['E. coli', 'S. aureus', 'C. albicans', 'C. tropicalis',
                  'C. parapsilosis', 'C. glabrata', '3T3 Fibroblasts', 'HUVECs', 'Red Blood Cells']

TOX_IDX      = [4, 5, 6]   # IC50_3T3, IC50_HUVEC, HC10
ACTIVITY_IDX = [0, 1, 2, 3, 7, 8]  # all MIC columns

MAKO_CMAP = ListedColormap(sns.color_palette("mako", 256))

plt.rcParams["figure.dpi"] = 144

TRAINING_DATA_PATH   = '../simpPepGeoMean.xlsx'
PREDICTION_DATA_PATH = '/home/dhc97/Peptide AMP Prediction Python/prevDataPepMR5.xlsx'

# ── Utility functions ──────────────────────────────────────────────────────────

def load_and_preprocess(path):
    df = pd.read_excel(path)
    df['ME100'] *= -1
    df['ME15']  *= -1
    return df.rename(columns={"Helicity": "% Helicity", "ME100": "heli_100", "ME15": "heli_15"})


def log2_selectivity(pred_MIC, pred_tox, actual_MIC, actual_tox):
    return actual_tox - actual_MIC, pred_tox - pred_MIC


def dilution_accuracy(predicted, actual):
    diff = np.abs(predicted - actual)
    n = actual.size
    return (sum(d < 1 for d in diff) / n * 100,
            sum(d < 2 for d in diff) / n * 100)


def calculate_pvalues(df):
    cols = df.columns
    pvalues = pd.DataFrame(np.nan, index=cols, columns=cols)
    for r in cols:
        for c in cols:
            if r != c:
                tmp = df[[r, c]].dropna()
                _, pval = pearsonr(tmp[r], tmp[c])
                pvalues.loc[r, c] = float(pval)
    return pvalues


def _calculate_vips(model):
    t, w, q = model.x_scores_, model.x_weights_, model.y_loadings_
    p, h = w.shape
    s = np.diag(t.T @ t @ q.T @ q)
    total_s = s.sum()
    vips = np.zeros(p)
    for i in range(p):
        weight = np.array([(w[i, j] / np.linalg.norm(w[:, j])) ** 2 for j in range(h)])
        vips[i] = np.sqrt(p * (s @ weight) / total_s)
    return vips


def _normalize(values, ref=None):
    """Scale values to [-1, 1] using ref range (defaults to values itself)."""
    ref = values if ref is None else ref
    lo, hi = ref.min(), ref.max()
    return 2 * (values - lo) / (hi - lo) - 1


def _comp_labels(dim):
    xR2, yR2 = VARIANCE_EXPLAINED["xR2"], VARIANCE_EXPLAINED["yR2"]
    return (f"Component {dim} ($R^2_X$ = {xR2[dim-1]}%, $R^2_Y$ = {yR2[dim-1]}%)")


# ── Plotting functions ─────────────────────────────────────────────────────────

def plot_prediction_accuracy(pred_MIC, pred_tox, actual_MIC, actual_tox, title=""):
    sel_actual, sel_pred = log2_selectivity(pred_MIC, pred_tox, actual_MIC, actual_tox)

    print(f"[{title}] MIC accuracy:        {dilution_accuracy(pred_MIC, actual_MIC)}")
    print(f"[{title}] Toxicity accuracy:   {dilution_accuracy(pred_tox, actual_tox)}")
    print(f"[{title}] Selectivity accuracy:{dilution_accuracy(sel_pred, sel_actual)}")

    ref = np.linspace(-3, 30, 150)
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 5))
    fig.suptitle(title)

    panels = [
        (ax1, actual_MIC,   pred_MIC,   "Actual log2 MIC",          "Predicted log2 MIC",          (0, 10)),
        (ax2, actual_tox,   pred_tox,   "Actual log2 Cytotoxicity",  "Predicted log2 Cytotoxicity", (0, 11)),
        (ax3, sel_actual,   sel_pred,   "Actual log2 Selectivity",   "Predicted log2 Selectivity",  (-2, 6)),
    ]
    for ax, x_data, y_data, xlabel, ylabel, lim in panels:
        x_data = np.array(x_data).flatten()
        y_data = np.array(y_data).flatten()
        ax.scatter(x_data, y_data, alpha=0.6)
        ax.plot(ref, ref,     'k',    alpha=0.5, label='Perfect')
        ax.plot(ref, ref + 1, '--k',  alpha=0.5, label='±2-fold')
        ax.plot(ref, ref - 1, '--k',  alpha=0.5)
        ax.plot(ref, ref + 2, '-.k',  alpha=0.5, label='±4-fold')
        ax.plot(ref, ref - 2, '-.k',  alpha=0.5)
        ax.annotate(f"R² = {r2_score(x_data, y_data):.3f}", xy=(0.05, 0.9), xycoords='axes fraction')
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(xlabel, size=13); ax.set_ylabel(ylabel, size=13)
        ax.set_box_aspect(1)
        ax.legend(loc='best', fontsize=9)


def plot_loadings_biplot(loadings_x, loadings_y, x_labels, y_labels):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot([], label="X loadings", color="darkslategrey")
    ax.plot([], label="Y loadings", color="crimson")

    for i, name in enumerate(x_labels):
        x, y = loadings_x.loc[i, 0], loadings_x.loc[i, 1]
        ax.arrow(0, 0, x, y, color='darkslategrey', width=0.007)
        ax.text(x, y, name, color='black', size=15,
                bbox=dict(boxstyle="round,pad=0.3", ec="k", fc="white", alpha=0.7))

    for i, name in enumerate(y_labels):
        x, y = loadings_y.loc[i, 0], loadings_y.loc[i, 1]
        ax.arrow(0, 0, x, y, color='crimson', width=0.007)
        ax.text(x, y, name, color='black', size=15,
                bbox=dict(boxstyle="round,pad=0.3", ec="k", fc="white", alpha=0.7))

    ax.set_title("PLSR X and Y Loadings", size=25, pad=15)
    ax.set_xlabel(_comp_labels(1), size=25)
    ax.set_ylabel(_comp_labels(2), size=25)
    ax.set_xlim([-1.2, 1.2]); ax.set_ylim([-1.2, 1.2])
    ax.set_aspect('equal')
    ax.legend(fontsize=20)


def plot_scores_biplot(scores_X, loadings_x, loadings_y, colorbar_values,
                       cbar_title, title, cbar_min=0, cbar_max=5, star_scores=None):
    sx = _normalize(scores_X[0])
    sy = _normalize(scores_X[1])

    fig, ax = plt.subplots(figsize=(10, 10))

    if star_scores is not None:
        star_x = _normalize(pd.Series([star_scores[0]]), ref=scores_X[0])[0]
        star_y = _normalize(pd.Series([star_scores[1]]), ref=scores_X[1])[0]
        ax.scatter(star_x, star_y, s=300, edgecolor='black', marker='*', zorder=5)
    else:
        ax.scatter(sx.iloc[0], sy.iloc[0], s=300, edgecolor='black', marker='*', zorder=5)

    sc = ax.scatter(sx, sy, s=100, c=colorbar_values, cmap=MAKO_CMAP, vmin=cbar_min, vmax=cbar_max)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_title, rotation=270, size=25, labelpad=20)

    ax.set_title(f"{title} vs Peptide Physicochemical Property", size=25, pad=15)
    ax.set_xlabel(_comp_labels(1), size=25)
    ax.set_ylabel(_comp_labels(2), size=25)
    ax.set_xlim([-1.2, 1.2]); ax.set_ylim([-1.2, 1.2])
    ax.set_aspect('equal')


def plot_3d_scores(scores_X, colorbar_values, cbar_title, title,
                   cbar_min=0, cbar_max=5, star_scores=None):
    sx = _normalize(scores_X[0])
    sy = _normalize(scores_X[1])
    sz = _normalize(scores_X[2])

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('white')

    if star_scores is not None:
        star_x = _normalize(pd.Series([star_scores[0]]), ref=scores_X[0])[0]
        star_y = _normalize(pd.Series([star_scores[1]]), ref=scores_X[1])[0]
        star_z = _normalize(pd.Series([star_scores[2]]), ref=scores_X[2])[0]
        ax.scatter(star_x, star_y, star_z, s=100, c='w', edgecolor='black', marker='*', zorder=5)
    else:
        ax.scatter(sx.iloc[0], sy.iloc[0], sz.iloc[0], s=100, c='w', edgecolor='black', marker='*', zorder=5)

    p = ax.scatter(sx, sy, sz, c=colorbar_values, cmap=MAKO_CMAP, vmin=cbar_min, vmax=cbar_max, edgecolor='k')
    fig.colorbar(p, ax=ax, fraction=0.03, pad=0.1).set_label(cbar_title, rotation=270, labelpad=20)

    ax.set_xlabel(_comp_labels(1), fontsize=7)
    ax.set_ylabel(_comp_labels(2), fontsize=7)
    ax.set_zlabel(_comp_labels(3), fontsize=7)
    ax.set_title(f"{title} vs Peptide Physicochemical Property", pad=15)


def plot_3d_loadings(scores_X, loadings_x, loadings_y, colorbar_values,
                     cbar_title, title, x_labels, y_labels, cbar_min=0, cbar_max=5):
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('white')

    ax.plot([], [], label="X loadings", color="darkslategrey")
    ax.plot([], [], label="Y loadings", color="crimson")

    origin8 = [0, 0, 0]
    X1, Y1, Z1 = zip(*[origin8] * len(x_labels))
    ax.quiver(X1, Y1, Z1,
              loadings_x.loc[:, 0], loadings_x.loc[:, 1], loadings_x.loc[:, 2],
              color='darkslategrey')
    for i, name in enumerate(x_labels):
        ax.text(*loadings_x.loc[i, :3], name, fontsize=7,
                bbox=dict(boxstyle="round,pad=0.3", ec="k", fc="white", alpha=0.7))

    origin9 = [0, 0, 0]
    X2, Y2, Z2 = zip(*[origin9] * len(y_labels))
    ax.quiver(X2, Y2, Z2,
              loadings_y.loc[:, 0], loadings_y.loc[:, 1], loadings_y.loc[:, 2],
              color='crimson')
    for i, name in enumerate(y_labels):
        ax.text(*loadings_y.loc[i, :3], name, fontsize=7,
                bbox=dict(boxstyle="round,pad=0.3", ec="k", fc="white", alpha=0.7))

    ax.set_xlabel(_comp_labels(1), fontsize=7)
    ax.set_ylabel(_comp_labels(2), fontsize=7)
    ax.set_zlabel(_comp_labels(3), fontsize=7)
    ax.set_title("PLSR X and Y Loadings", pad=15)
    ax.set_xlim3d(-0.75, 0.75); ax.set_ylim3d(-0.75, 0.75); ax.set_zlim3d(-0.75, 0.75)
    ax.legend()


# ── Load and train ─────────────────────────────────────────────────────────────

df = load_and_preprocess(TRAINING_DATA_PATH)
data = df[DESCRIPTORS + OUTPUT_VARS]
X_train = data[DESCRIPTORS]
Y = data[OUTPUT_VARS]

model = PLSRegression(n_components=3, scale=True)
model.fit(X_train, Y)
train_pred = model.predict(X_train)

# ── Training set selectivity ───────────────────────────────────────────────────

avg_tox_pred_tr  = train_pred[:, TOX_IDX].mean(axis=1)
avg_act_pred_tr  = train_pred[:, ACTIVITY_IDX].mean(axis=1)
pred_selectivity_tr = avg_tox_pred_tr - avg_act_pred_tr

avg_tox_true_tr  = Y.values[:, TOX_IDX].mean(axis=1)
avg_act_true_tr  = Y.values[:, ACTIVITY_IDX].mean(axis=1)
actual_selectivity_tr = avg_tox_true_tr - avg_act_true_tr

scores_train = pd.DataFrame(model.x_scores_)
loadings_x   = pd.DataFrame(model.x_loadings_)
loadings_y   = pd.DataFrame(model.y_loadings_)
training_star = scores_train.iloc[0]  # reference peptide carried into prediction plots

# ── Training set plots ─────────────────────────────────────────────────────────

plot_loadings_biplot(loadings_x, loadings_y, DESCRIPTORS, OUTPUT_VARS)

plot_scores_biplot(scores_train, loadings_x, loadings_y, actual_selectivity_tr,
                   cbar_title='Log2 Selectivity', title='Actual Broad Spectrum Selectivity',
                   cbar_min=-1, cbar_max=3.5)

plot_scores_biplot(scores_train, loadings_x, loadings_y, pred_selectivity_tr,
                   cbar_title='Log2 Selectivity', title='Predicted Broad Spectrum Selectivity',
                   cbar_min=-1, cbar_max=3.5)

plot_prediction_accuracy(avg_act_pred_tr, avg_tox_pred_tr,
                         avg_act_true_tr, avg_tox_true_tr, title='Broad Spectrum')
for i in range(9):
    plot_prediction_accuracy(train_pred[:, i], avg_tox_pred_tr,
                             Y.values[:, i], avg_tox_true_tr, title=OUTPUT_VARS[i])

# ── Load prediction set ────────────────────────────────────────────────────────

df2 = load_and_preprocess(PREDICTION_DATA_PATH)
X_pred  = df2[DESCRIPTORS]
Y_test  = df2[['MIC', 'Hemolysis']]
pep_ids = df2['# of Peptide'].astype(str).tolist()

preds = model.predict(X_pred)

HC10_sel_pred   = preds[:, 6] - preds[:, 0]   # HC10 - MIC_CA
HC10_sel_actual = Y_test.values[:, 1] - Y_test.values[:, 0]

avg_tox  = preds[:, TOX_IDX].mean(axis=1)
avg_act  = preds[:, ACTIVITY_IDX].mean(axis=1)
avg_sel  = avg_tox - avg_act

plot_prediction_accuracy(preds[:, 0], preds[:, 6],
                         Y_test.values[:, 0], Y_test.values[:, 1], title='C. albicans validation')

# ── Prediction set biplots ─────────────────────────────────────────────────────

scores_pred = pd.DataFrame(model.transform(X_pred))

plot_scores_biplot(scores_pred, loadings_x, loadings_y, HC10_sel_actual,
                   cbar_title='Log2 Selectivity', title='Actual C. albicans Selectivity',
                   cbar_min=-2, star_scores=training_star)

plot_scores_biplot(scores_pred, loadings_x, loadings_y, HC10_sel_pred,
                   cbar_title='Log2 Selectivity', title='Predicted C. albicans Selectivity',
                   cbar_min=-2, star_scores=training_star)

plot_scores_biplot(scores_pred, loadings_x, loadings_y, avg_sel,
                   cbar_title='$Log_2$ Broad Spectrum Selectivity',
                   title='Predicted Broad Spectrum Selectivity',
                   cbar_min=-1, cbar_max=3.5, star_scores=training_star)

plot_3d_scores(scores_pred, HC10_sel_pred,
               cbar_title='Log2 Selectivity', title='Predicted C. albicans Selectivity',
               cbar_min=-2, star_scores=training_star)

plot_3d_scores(scores_pred, HC10_sel_actual,
               cbar_title='Log2 Selectivity', title='Actual C. albicans Selectivity',
               cbar_min=-2, star_scores=training_star)

plot_3d_scores(scores_pred, avg_sel,
               cbar_title='Log2 Selectivity', title='Predicted Broad Spectrum Selectivity',
               cbar_min=-1, cbar_max=3.5, star_scores=training_star)

plot_3d_loadings(scores_pred, loadings_x, loadings_y, avg_sel,
                 cbar_title='Log2 Selectivity', title='Predicted Broad Spectrum Selectivity',
                 x_labels=DESCRIPTORS, y_labels=OUTPUT_VARS, cbar_min=-1, cbar_max=3.5)

# ── Prediction heatmap sorted by selectivity ───────────────────────────────────

REORDER = [8, 7, 0, 1, 3, 2, 4, 5, 6]  # rearrange columns to match species label order
preds_sorted = pd.DataFrame(preds[:, REORDER], index=pep_ids, columns=SPECIES_LABELS)
preds_sorted['avgSelectivity'] = avg_sel
preds_sorted = preds_sorted.sort_values('avgSelectivity', ascending=False)

plt.figure()
preds_labelled = preds_sorted[SPECIES_LABELS].reset_index(drop=True)
ax = sns.heatmap(preds_labelled, cmap='BuPu', yticklabels=False,
                 cbar_kws={'label': 'Peptide Toxicity Concentration (log2)'})
ax.set_ylabel('Peptides sorted by high to low BSI')

# ── Pearson correlation heatmap ────────────────────────────────────────────────

corr = preds_labelled.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
pval = calculate_pvalues(preds_sorted[SPECIES_LABELS])

plt.figure(figsize=(12, 10))
sns.heatmap(corr, mask=mask, cmap='coolwarm', vmin=-1, vmax=1, square=True,
            cbar_kws={'label': 'Pearson correlation coefficient'})

# ── Coefficient heatmap ────────────────────────────────────────────────────────

# x_rotations_ @ y_loadings_.T replicates the pre-sklearn-1.1 coef_ definition
coef = pd.DataFrame(
    (model.x_rotations_ @ model.y_loadings_.T)[:, REORDER],
    index=DESCRIPTORS, columns=SPECIES_LABELS
)
coef_sorted = coef.sort_values('E. coli')

plt.figure(figsize=(6, 4))
sns.heatmap(coef_sorted, cmap='coolwarm_r', vmin=-1.5, vmax=1.5,
            cbar_kws={'label': 'Coefficient'})

# ── Feature importance ─────────────────────────────────────────────────────────

vips = _calculate_vips(model)
vips_df = pd.DataFrame(vips, index=DESCRIPTORS, columns=['VIP values'])
vips_df = vips_df.sort_values('VIP values', ascending=False)

plt.figure()
vips_df.plot.bar(ax=plt.gca(), color='tab:purple', ylabel='VIP value', legend=False)

perm = permutation_importance(model, X_train, Y, n_repeats=1000)
perms_df = pd.DataFrame({
    'Permutation importance': perm.importances_mean,
    'std': perm.importances_std,
}, index=DESCRIPTORS).sort_values('Permutation importance', ascending=False)

plt.figure()
perms_df['Permutation importance'].plot(kind='bar', yerr=perms_df['std'])
plt.ylabel('Permutation importance')

plt.show()
