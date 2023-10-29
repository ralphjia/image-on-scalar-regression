#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 02:14:13 2023

@author: yutingduan
"""

import argparse

import numpy as np
import pandas as pd
from utils import load_pickle, save_pickle
from spm import spm, get_fwhm


if __name__ == '__main__':  
    images = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/img.csv",
                     delimiter=",", skip_header=1)
    covariates = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cov_mat.csv",
                               delimiter=",", skip_header=1)

    images = images.transpose()
    
    coord = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/voxels.csv",
                                delimiter=",", skip_header=1)
    
    fwhm = get_fwhm(
            images, covariates, coord, img_shape=None,
            alpha=0.05, n_cpus=None)
    effect = spm(images, covariates, coord, shape=None,fwhm=fwhm, alpha=0.05)[0]
    
    beta_SPM = pd.DataFrame(effect)
    beta_SPM.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/SPMbeta.csv")


# true_beta = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/data/true_beta_num_obs_50sigma_0.1dist_noise_normali_1.csv",
#                           delimiter=",", skip_header=1).transpose()
# mse = np.square(np.subtract(true_beta, effect)).mean()

# if __name__ == '__main__':  
#     MSE_SPM = []
#     for num_obs in [50]:
#         for sigma in [0.5, 1, 2, 4, 6]:
#             for i in list(range(1,6)):
    
#                 images = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/img_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                  delimiter=",", skip_header=1)
#                 covariates = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/cov_mat_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                            delimiter=",", skip_header=1)
        
#                 images = images.transpose()
                
#                 coord = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/voxels_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                            delimiter=",", skip_header=1)
                
#                 fwhm = get_fwhm(
#                         images, covariates, coord, img_shape=None,
#                         alpha=0.05, n_cpus=None)
#                 effect = spm(images, covariates, coord, shape=None,fwhm=fwhm, alpha=0.05)[0]
#                 beta_SPM = pd.DataFrame(effect)
#                 beta_SPM.to_csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/SPMbeta_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv")
#                 true_beta = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/true_beta_num_obs_"+str(num_obs)+"sigma_"+str(sigma)+"dist_noise_normali_"+str(i)+".csv",
#                                           delimiter=",", skip_header=1).transpose()
#                 mse = np.square(np.subtract(true_beta, effect)).mean()
#                 MSE_SPM.append(mse)
            
        

# SPM = np.asarray(MSE_SPM).T

# SPM.tofile('MSE_SPM.csv', sep = ',')


   