import numpy as np
import pandas as pd
from scipy.stats import t


def within_normalized_CI(df, x, y, group, conf=.95):
    """Compute within-subjects confidence intervals through
    normalization. Based on method described in Morey (2008)
    and Franz and Loftus (2012).

    df (pandas dataframe) : data
    x  (string)           : column name for within-subjects factor
    y  (string)           : column name for dependent variable
    group (string)        : column name for grouping variable (e.g., subject id)
    conf (float)          : confidence level
    """
    p = 1 - (1 - conf)/2.

    # number of within-subjects conditions
    M = df[x].unique().shape[0]

    # normalize observations
    normalize = lambda grp: grp[y] - grp[y].mean() + df[y].mean()
    df.loc[:,'y_norm'] = df.groupby(group).apply(normalize).reset_index()[y].values

    CI = []
    for i, rep in df.groupby(x):
        z = rep['y_norm'].values
        n = z.shape[0]
        cv = t.ppf(p, n)

        # standard error
        se = np.sqrt((M/float(M - 1)) * (1./(n*(n-1))) * np.sum((z - z.mean()) ** 2))

        # 95% CI
        CI.append([rep[x].values[0], cv * se])

    return pd.DataFrame(np.array(CI), columns=[x, '95-CI(within)'])


# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)



def add_letter_label(fig, x, y, label):
    fig.text(x, y, label, fontsize=12, fontname='Helvetica Neue', fontweight='bold')

