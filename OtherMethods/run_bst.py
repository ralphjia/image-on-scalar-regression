import sys
import multiprocessing

import numpy as np
import pandas as pd

from scipy.stats import norm as normal_rv
#from statsmodels.stats.multitest import multipletests

from utils import load_pickle, save_pickle
from fdaimage import fdaimage




def fit_bst(y, x, alpha=0.05):
    assert y.ndim == 3
    assert y.shape[0] == x.shape[0]
    img_shape = y.shape[1:]
    coords = np.stack(np.meshgrid(
        np.linspace(0, 1, img_shape[0]),
        np.linspace(0, 1, img_shape[1]),
        indexing='ij'), -1)
    effect, se, ci_lower, ci_upper, alpha_adj = fdaimage(
            y.reshape(y.shape[0], -1), x,
            coords.reshape(-1, coords.shape[-1]),
            alpha=alpha)
    effect = effect.reshape((-1,) + img_shape)
    ci_lower = ci_lower.reshape((-1,) + img_shape)
    ci_upper = ci_upper.reshape((-1,) + img_shape)
    # pval = normal_rv.sf(np.abs(effect / se)) * 2
    # sig = multipletests(pval.flatten(), alpha=alpha, method='fdr_bh')[0]
    # sig = sig.reshape(effect.shape)
    sig = np.logical_or(ci_lower > 0, ci_upper < 0)
    return effect, effect * sig


def fit_bst_3d(y, x, ncpu=None):
    assert y.ndim == 4
    assert y.shape[0] == x.shape[0]
    n_slices = y.shape[-1]
    y_slices = [y[..., i] for i in range(n_slices)]
    if ncpu is None:
        ncpu = min(n_slices, multiprocessing.cpu_count())
    print('ncpu:', ncpu)
    with multiprocessing.Pool(processes=ncpu) as pool:
        out = pool.starmap(fit_bst, zip(y_slices, [x]*n_slices))
    effect_dense = np.stack([e[0] for e in out], -1)
    effect = np.stack([e[1] for e in out], -1)
    return effect_dense, effect


if __name__ == '__main__':
    images = np.genfromtxt(
        "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/img.csv",
        delimiter=",", skip_header=1)
    covariates = np.genfromtxt(
        "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cov_mat.csv",
        delimiter=",", skip_header=1)

    images = images.transpose()

    coord = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/voxels.csv",
                          delimiter=",", skip_header=1)

    coord[:, 0] = coord[:, 0] - 1
    coord[:, 1] = coord[:, 1] - 1
    coord[:, 2] = coord[:, 2] - 1
    coord = coord.astype(np.int64)

    # true_beta4D = true_beta.reshape(3, 32, 32, 8, order="F")
    img4D = images.reshape(50, 20, 20, 8, order="F")

    # true_beta4D = np.zeros((np.shape(true_beta)[0],
    #                         len(np.unique(coord[:, 0])),
    #                         len(np.unique(coord[:, 1])),
    #                         len(np.unique(coord[:, 2]))
    #                         ))

    # for i in list(range(np.shape(true_beta)[0])):
    #     # img4D.flat[np.append(i * np.ones((np.shape(images)[1], 1)), coord, axis=1)] = images[i, :]
    #     np.put(true_beta4D,
    #            np.append(i * np.ones((np.shape(true_beta)[1], 1)), coord, axis=1).astype(int),
    #            true_beta[i, :],
    #            mode='wrap')

    # img4D = np.zeros((np.shape(images)[0],
    #                   len(np.unique(coord[:, 0])),
    #                   len(np.unique(coord[:, 1])),
    #                   len(np.unique(coord[:, 2]))
    #                   ))
    # img4D = img4D / 4

    # for i in list(range(np.shape(images)[0])):
    #     # img4D.flat[np.append(i * np.ones((np.shape(images)[1], 1)), coord, axis=1)] = images[i, :]
    #     np.put(img4D,
    #            np.append(i * np.ones((np.shape(images)[1], 1)), coord, axis=1).astype(int),
    #            images[i, :],
    #            mode='wrap')

    y = img4D
    x = covariates

    n_slices = 8
    stride = y.shape[-1] // n_slices
    maineff_pred_dense, maineff_pred = fit_bst_3d(y[:, ..., ::stride], x)
    maineff_pred_dense = np.repeat(maineff_pred_dense, stride, axis=-1)
    maineff_pred = np.repeat(maineff_pred, stride, axis=-1)
    maineffBST = maineff_pred.reshape(3, 3200, order="F")
    maineffBST = pd.DataFrame(maineffBST)
    maineffBST.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/BSTbeta.csv")
# mse = np.square(np.subtract(true_beta4D, maineff_pred)).mean()
# MSE_BST.append(mse)

# value = {'MSE' : MSE_BST}
# df = pd.DataFrame(value)

# saving the dataframe
# df.to_csv('MSE_BST1.csv')
