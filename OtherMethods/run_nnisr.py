import random
from time import time

import numpy as np
import torch

from utils import write_text, load_pickle, save_pickle
from nnisr import fit_nnisr, get_coord


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)



images = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/data/img_num_obs_50sigma_0.1dist_noise_normali_1.csv",
                  delimiter=",", skip_header=1)
covariates = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/data/cov_mat_num_obs_50sigma_0.1dist_noise_normali_1.csv",
                            delimiter=",", skip_header=1)

images = images.transpose()

coord = np.genfromtxt("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/data/voxels_num_obs_50sigma_0.1dist_noise_normali_1.csv",
                            delimiter=",", skip_header=1)
n_voxels = images.shape[-1]

set_seed(0)


pred = fit_nnisr(
        x=covariates, y=images, s=coord, img_shape=None,
        hidden_widths=(256,)*4,
        activation='leaky',
        alpha_threshold=0.05,
        n_permute=100, lr=1e-3,
        epochs=50, batch_size=4096,
        n_states=11, alpha_states=0.5,
        prefix='nnisr/', device='cpu', n_jobs=None)