
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from utils import load_pickle, save_pickle, fit_ls


def fit_mua(images, covariates, alpha=0.05):
    effect, __, pval = fit_ls(
            images.reshape(images.shape[0], -1), covariates)[:3]
    n_covariates = covariates.shape[1]
    img_shape = images.shape[1:]
    effect = effect.reshape((n_covariates,) + img_shape)
    sig = multipletests(pval.flatten(), alpha=alpha, method='fdr_bh')[0]
    sig = sig.reshape(effect.shape)
    return effect, effect * sig


images = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/img.csv",
                 delimiter=",", skip_header=1)
covariates = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cov_mat.csv",
                           delimiter=",", skip_header=1)

images = images.transpose()

effect = fit_mua(images, covariates)[0]

beta_MUA = pd.DataFrame(effect)
beta_MUA.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/MUAbeta.csv")

#MSE_MUA = []


# for num_obs in [50, 100]:
#     for sigma in [0.5, 1, 2, 4, 6]:
#         for i in list(range(1, 6)):
#             images = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/img_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                              delimiter=",", skip_header=1)
#             covariates = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/cov_mat_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                        delimiter=",", skip_header=1)

#             images = images.transpose()
#             effect = fit_mua(images, covariates)[0]
#             beta_MUA = pd.DataFrame(effect)
#             beta_MUA.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/MUAbeta_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv")
#             true_beta = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/true_beta_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                       delimiter=",", skip_header=1).transpose()
#             mse = np.square(np.subtract(true_beta, effect)).mean()
#             MSE_MUA.append(mse)
            
            
