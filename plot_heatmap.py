
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def heatmap(prec, level=0, lim='a',cmap="coolwarm"):

    d = prec.shape[0]
    mask = np.zeros_like(prec)
    font_setting = {'fontsize': 5}

    for i in range(0, d):
        for j in range(0, d):
            if abs(prec[i][j]) <= level:
                mask[i][j] = True
    anc = prec.copy()
    for i in range(0, d):
        anc[i, i] = 0
    if lim == 'a':
        vmax = np.max([np.max(np.abs(anc)), 1e-1])
        vmin = -vmax
    elif lim == 'b':
        vmax = np.max([np.max(np.abs(anc)), 1e-1])
        vmin = 0
    else:
        vmax = lim[0]
        vmin = lim[1]
    ax = sns.heatmap(prec, annot=False, annot_kws=font_setting, cmap=cmap, vmin=vmin, vmax=vmax,
                square=True, mask=mask)
    ax.set_facecolor('.7')
    plt.savefig("./temp.png", dpi=600, bbox_inches='tight')
    plt.show()

