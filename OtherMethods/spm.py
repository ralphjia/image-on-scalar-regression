import multiprocessing

import numpy as np
from scipy.stats import norm as norm_rv

from utils import fit_ls
from nnisr import fit_smooth


def prob_ec_nonzero(z, nsearch, fwhm):
    # worsley1992three
    return (
            nsearch / fwhm.prod() * (4 * np.log(2))**(3/2)
            * (2 * np.pi)**(-2) * (z**2 - 1) * np.exp(-0.5 * z**2))


def get_cutoff_rft(nsearch, fwhm, alpha, ngrid=1000):

    z_course = np.logspace(0, 2, ngrid+1)[1:]
    idx_course = np.where(
            prob_ec_nonzero(z_course, nsearch, fwhm) < alpha)[0][0]
    is_border = (idx_course > 0) and (idx_course < ngrid-1)
    if is_border:
        cutoff = z_course[idx_course]
    else:
        z_start, z_stop = z_course[idx_course - 1], z_course[idx_course+1]
        z = np.logspace(np.log10(z_start), np.log10(z_stop), ngrid+1)[1:]
        idx = np.where(prob_ec_nonzero(z, nsearch, fwhm) < alpha)[0][0]
        cutoff = z[idx]
    cutoff_min = norm_rv.ppf(alpha * 0.5) * (-1)
    cutoff = max(cutoff_min, cutoff)
    return cutoff


def spm(images, covariates, coord, shape, fwhm, alpha):
    if coord is None:
        size = np.prod(shape)
    else:
        size = coord.shape[0]
    stde = fwhm / np.sqrt(8*np.log(2))

    # smoothed = gaussian_filter(
    #         images,
    #         np.concatenate([[0], stde])).reshape(-1, size)
    smoothed = fit_smooth(
            value=images.T,
            coord=coord, img_shape=shape,
            filter_size=stde).T

    beta, se, = fit_ls(smoothed, covariates, alpha=alpha)[:2]
    score = np.abs(beta) / se

    if coord is None:
        fwhm = max(shape) * fwhm
    else:
        fwhm = (coord.max() - coord.min()) * fwhm

    cutoff_rft = get_cutoff_rft(size, fwhm, alpha)
    # cutoff_rft *= 1 - alpha
    sig_familywise = score > cutoff_rft
    return beta, beta * sig_familywise, cutoff_rft * se


def spm_cv(images, covariates, coord, shape, alpha, fwhm, n_folds=5):
    n_observations = images.shape[0]
    indices = list(range(n_observations))
    np.random.seed(0)
    np.random.shuffle(indices)
    ind_batches = np.array_split(indices, n_folds)
    loss_list = []
    for ind_bat in ind_batches:
        is_test = np.isin(np.arange(n_observations), ind_bat)
        effect = spm(
                images=images[~is_test],
                covariates=covariates[~is_test],
                coord=coord, shape=shape,
                fwhm=fwhm, alpha=alpha)[1]
        images_pred = covariates[is_test] @ effect
        loss = np.square(images_pred - images[is_test]).mean()
        loss_list.append(loss)
    loss = np.mean(loss_list)
    print('loss:', loss, 'fwhm:', fwhm)
    return loss


def get_fwhm(
        images, covariates, coord, img_shape, alpha,
        n_fwhm=20, fwhm_min=1e-3, fwhm_max=10, n_cpus=None, depth=1):
    fwhm_list = np.logspace(np.log10(fwhm_min), np.log10(fwhm_max), n_fwhm)
    if coord is None:
        img_ndim = len(img_shape)
    else:
        img_ndim = coord.shape[-1]
    fwhm_list = np.tile(fwhm_list, [img_ndim, 1]).T

    if n_cpus is None:
        n_cpus = min(n_fwhm, multiprocessing.cpu_count())
    print('n_cpus:', n_cpus)
    if n_cpus == 1:
        loss_list = [
                spm_cv(
                    images=images, covariates=covariates,
                    coord=coord, shape=img_shape, alpha=alpha, fwhm=fwhm)
                for fwhm in fwhm_list]
    else:
        with multiprocessing.Pool(processes=n_cpus) as pool:
            loss_list = pool.starmap(
                    spm_cv,
                    zip(
                        [images]*n_fwhm, [covariates]*n_fwhm,
                        [coord]*n_fwhm, [img_shape]*n_fwhm,
                        [alpha]*n_fwhm,
                        fwhm_list))

    i_best = np.argmin(loss_list)
    fwhm_best = fwhm_list[i_best]
    if i_best in (0, n_fwhm - 1):
        print('Warning: best parameter on boundary')
    else:
        if depth > 0:
            fwhm_best = get_fwhm(
                images=images, covariates=covariates, coord=coord,
                img_shape=img_shape, alpha=alpha, n_cpus=n_cpus,
                fwhm_min=fwhm_list[i_best-1][0],
                fwhm_max=fwhm_list[i_best+1][0],
                n_fwhm=n_fwhm, depth=depth-1)

    return fwhm_best
