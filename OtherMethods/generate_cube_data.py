import random

import numpy as np

from design import gen_data, get_coord
from utils import save_pickle
import pandas as pd

img_shape = (32, 32, 8)

# data = gen_data(
#             V_out=img_shape, N=50, Q=4,
#             beta_stn=1, omega_stn=0.5,
#             omega_itv=1.0,
#             noise_dist='gauss',
#             noise_var='wave', scale=1.0, cut=True)
s = get_coord(img_shape)
# coord = pd.DataFrame(s)
# coord.to_csv('/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/coord.csv')
#
# save_pickle(data, '/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/'
#                   +'V_out=(32,32,8), N=50, Q=4, beta_stn=1, omega_stn=0.5, omega_itv=1.0.pickle')

for N in [50, 100]:
    for beta_stn in [0.05, 0.2, 0.5, 1, 2]:
        data = gen_data(
            V_out=img_shape, N=N, Q=4,
            beta_stn=beta_stn, omega_stn=0.5,
            omega_itv=1.0,
            noise_dist='gauss',
            noise_var='wave', scale=1.0, cut=True)
        save_pickle(data, '/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/'
                    + 'N=' + str(N) + ', beta_stn=' + str(beta_stn) + '.pickle')
        y = data[0]
        y = y.reshape(y.shape[0], -1)
        img = pd.DataFrame(y)
        img.to_csv('/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/img_'
                    + 'N=' + str(N) + ', beta_stn=' + str(beta_stn) + '.csv')
        x = data[1]
        xx = pd.DataFrame(x)
        xx.to_csv('/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/xx_'
                   + 'N=' + str(N) + ', beta_stn=' + str(beta_stn) + '.csv')
        beta = data[2]
        betavv = beta.reshape(beta.shape[0], -1)
        betavv = pd.DataFrame(betavv)
        betavv.to_csv('/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/beta_'
                   + 'N=' + str(N) + ', beta_stn=' + str(beta_stn) + '.csv')


