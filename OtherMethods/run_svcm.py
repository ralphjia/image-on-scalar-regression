import sys

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

from utils import load_pickle, save_pickle
from svcm import svcm


def fit_svcm(y, x, alpha=0.05):
    assert y.shape[0] == x.shape[0]
    effect, pval, __, __ = svcm(y, x, fromfile=False)
    sig = multipletests(pval.flatten(), alpha=alpha, method='fdr_bh')[0]
    sig = sig.reshape(effect.shape)
    return effect, effect * sig


MSE_SVCM = []
for num_obs in [50]:
    for sigma in [0.5, 1, 2, 4, 6]:
        for i in list(range(1,6)):

            images = np.genfromtxt(
                "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/img_num_obs_" + str(
                    num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
                delimiter=",", skip_header=1)
            covariates = np.genfromtxt(
                "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/cov_mat_num_obs_" + str(
                    num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
                delimiter=",", skip_header=1)

            images = images.transpose()

            coord = np.genfromtxt(
                "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/voxels_num_obs_" + str(
                    num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
                delimiter=",", skip_header=1)


            true_beta = np.genfromtxt(
                "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/true_beta_num_obs_" + str(
                    num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
                delimiter=",", skip_header=1).transpose()

            coord[:, 0] = coord[:, 0] - 1
            coord[:, 1] = coord[:, 1] - 1
            coord[:, 2] = coord[:, 2] - 1
            coord = coord.astype(np.int64)

            true_beta4D = true_beta.reshape(3, 32, 32, 8, order="F")
            img4D = images.reshape(num_obs,32,32,8, order="F")



            images = img4D

            n_slices = 5
            stride = images.shape[-1] // n_slices
            maineff_pred_dense, maineff_pred = fit_svcm(images[..., ::stride], covariates)
            maineff_pred_dense = np.repeat(maineff_pred_dense, stride, axis=-1)
            maineff_pred = np.repeat(maineff_pred, stride, axis=-1)

            maineffSVCM = maineff_pred.reshape(3, 8192, order="F")
            maineffSVCM = pd.DataFrame(maineffSVCM)
            maineffSVCM.to_csv(
                "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/SVCMbeta_num_obs_" + str(
                    num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv")
            mse = np.square(true_beta - maineffSVCM).mean()
            MSE_SVCM.append(mse)

value = {'MSE' : MSE_SVCM}
df = pd.DataFrame(value)

    # saving the dataframe
df.to_csv('MSE_SVCM.csv')


