import argparse
import random
from time import time

import numpy as np
import torch

from utils import write_text, load_pickle, save_pickle
from irrnn import fit_irrnn, get_coord

from metrics import (
        mean_squared_difference, true_positive, true_negative,
        false_positive, false_negative,
        true_discovery, true_omission,
        false_discovery, false_omission, rocauc)

import pandas as pd
import numpy as np

if __name__ == '__main__':
    MSE_IRRNN = []
    fp_IRRNN = []
    fn_IRRNN = []
    PPV_IRRNN = []
    NPV_IRRNN = []
    for sigma in [2, 4]:
        for outlier_ratio in [0.02, 0.05, 0.1, 0.2]:
            for outlier_mag in [10, 20, 50]:
                
                true_beta = np.genfromtxt(
                    "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OutlierData/true_beta_num_obs50sigma"+str(sigma)+"outlier_ratio"+str(outlier_ratio)+"outlier_mag"+str(outlier_mag)+".csv",
                    delimiter=",", skip_header=1)
                y = np.genfromtxt(
                    "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OutlierData/img_num_obs50sigma"+str(sigma)+"outlier_ratio"+str(outlier_ratio)+"outlier_mag"+str(outlier_mag)+".csv",
                    delimiter=",", skip_header=1)
                x = np.genfromtxt(
                    "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OutlierData/cov_mat_num_obs50sigma"+str(sigma)+"outlier_ratio"+str(outlier_ratio)+"outlier_mag"+str(outlier_mag)+".csv",
                    delimiter=",", skip_header=1)
                
                # true_beta = np.genfromtxt(
                #     "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/brain_test/true_beta_brain.csv",
                #     delimiter=",", skip_header=1)
                
                
                y = y.transpose()
                
                s = np.genfromtxt(
                    "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/voxels.csv",
                    delimiter=",", skip_header=1)
                
                
                
                img_shape = None
                n_voxels = y.shape[-1]
                batch_size = n_voxels // 16
                prefix = 'brain'
                print('Fitting irrnn...')
                t0 = time()
                pred = fit_irrnn(
                        x=x, y=y, s=s, img_shape=img_shape,
                        hidden_widths=(256,)*4,
                        activation='leaky',
                        alpha_threshold=0.05,
                        n_permute=100, lr=1e-3,
                        epochs=50, batch_size=batch_size,
                        max_iter=2, n_states=11,
                        alpha_states=0.5,
                        prefix=prefix, device='cpu', n_jobs=None)
                print(int(time() - t0), 'sec')
                maineff_pred = pred['maineff']
                maineff_pred = maineff_pred.transpose()
                
                mse = mean_squared_difference(true_beta, maineff_pred)
                falpos = false_positive(true_beta, maineff_pred)
                falneg = false_negative(true_beta, maineff_pred)
                NPV = true_omission(true_beta, maineff_pred)
                PPV = true_discovery(true_beta, maineff_pred)
                MSE_IRRNN = [MSE_IRRNN, mse]
                fp_IRRNN = [fp_IRRNN, falpos]
                fn_IRRNN = [fn_IRRNN, falneg]
                PPV_IRRNN = [PPV_IRRNN, PPV]
                NPV_IRRNN = [NPV_IRRNN, NPV]
                
                IRRNNbeta = pd.DataFrame(maineff_pred)
                IRRNNbeta.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OutlierData/IRRNNbeta_num_obs_50sigma"+str(sigma)+"outlier_ratio"+str(outlier_ratio)+"outlier_mag"+str(outlier_mag)+".csv")
            

# if __name__ == '__main__':
    
#     MSE_IRRNN = []
#     for i in list(range(1,6)):
#         for num_obs in [50]:
#             for sigma in [0.5, 1, 2, 4, 6]:

                
#                 y = np.genfromtxt(
#                     "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/img_num_obs_" + str(
#                         num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
#                     delimiter=",", skip_header=1)
#                 x = np.genfromtxt(
#                     "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/cov_mat_num_obs_" + str(
#                         num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
#                     delimiter=",", skip_header=1)
                
#                 y = y.transpose()
                
#                 s = np.genfromtxt(
#                     "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/voxels_num_obs_" + str(
#                         num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
#                     delimiter=",", skip_header=1)
                
                
#                 true_beta = np.genfromtxt(
#                     "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/true_beta_num_obs_" + str(
#                         num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv",
#                     delimiter=",", skip_header=1).transpose()
                
#                 img_shape = None
#                 n_voxels = y.shape[-1]
#                 batch_size = n_voxels // 16
#                 prefix = 'num_obs' + str(num_obs) + '_sigma' + str(sigma) + '_i' + str(i) 
#                 print('Fitting irrnn...')
#                 t0 = time()
#                 pred = fit_irrnn(
#                         x=x, y=y, s=s, img_shape=img_shape,
#                         hidden_widths=(256,)*4,
#                         activation='leaky',
#                         alpha_threshold=0.05,
#                         n_permute=100, lr=1e-3,
#                         epochs=50, batch_size=batch_size,
#                         max_iter=2, n_states=11,
#                         alpha_states=0.5,
#                         prefix=prefix, device='cpu', n_jobs=None)
#                 print(int(time() - t0), 'sec')
#                 maineff_pred = pred['maineff']
#                 IRRNNbeta = pd.DataFrame(maineff_pred)
#                 IRRNNbeta.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/IRRNNbeta_num_obs_" + str(
#                         num_obs) + "sigma_" + str(sigma) + "dist_noise_normali_" + str(i) + ".csv")
#                 mse = np.square(true_beta - maineff_pred).mean()
#                 MSE_IRRNN.append(mse)
                
#     value = {'MSE' : MSE_IRRNN}
#     df = pd.DataFrame(value)
   
#     # saving the dataframe
#     df.to_csv('MSE_IRRNN.csv')


