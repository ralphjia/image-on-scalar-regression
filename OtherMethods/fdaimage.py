import numpy as np

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

from utils import erint


rpy2.robjects.numpy2ri.activate()
r = robjects.r
r['source']('fdaimage.R')


def get_coordinates(dim):
    d = np.prod(dim)
    coords = np.c_[np.unravel_index(np.arange(d).reshape(-1, 1), dim)]
    coords = coords / coords.max(0, keepdims=True)
    return coords


def fdaimage(
        Y, X, Z, alpha=0.05, nboot=10,
        boundary='rectangle', return_sigma=True, verbose=False):

    assert isinstance(Y, np.ndarray)
    assert isinstance(X, np.ndarray)
    assert Y.ndim == 2
    assert X.ndim == 2
    assert Z.ndim == 2
    assert X.shape[0] == Y.shape[0]
    assert Y.shape[1] == Z.shape[0]
    assert Z.shape[1] == 2

    out = r['fdaimage'](
            Y, X, Z, boundary=boundary, alpha=alpha, nboot=nboot,
            return_sigma=return_sigma, verbose=verbose, debug=False)
    beta, sigma, ci_lower, ci_upper, alpha_adj = out
    beta = np.array(beta).T
    sigma = np.array(sigma).T
    ci_lower = np.array(ci_lower).T
    ci_upper = np.array(ci_upper).T
    return beta, sigma, ci_lower, ci_upper, alpha_adj


def fdaimage3d(
        Y, X, Z, boundary='rectangle', return_sigma=True, verbose=False):

    assert isinstance(Y, np.ndarray)
    assert isinstance(X, np.ndarray)
    assert isinstance(Z, np.ndarray)
    assert Y.ndim == 2
    assert X.ndim == 2
    assert Z.ndim == 2
    assert X.shape[0] == Y.shape[0]
    assert Y.shape[1] == Z.shape[0]
    assert Z.shape[1] == 3

    beta = np.zeros((X.shape[1], Y.shape[1]))
    radi = np.zeros_like(beta)
    idxs, ys, xs, zs = [], [], [], []
    for w in np.unique(Z[:, 2]):
        idx = Z[:, 2] == w
        idxs.append(idx)
        ys.append(Y[:, idx])
        xs.append(X)
        zs.append(Z[idx, :2])

    for idx, y, x, z in zip(idxs, ys, xs, zs):
        if verbose:
            erint(Z[idx][0])
        b, r = fdaimage(
                y, x, z, boundary=boundary,
                return_sigma=return_sigma, verbose=verbose)
        beta[:, idx] = b
        if np.isfinite(r).all() and np.min(r) >= 0:
            radi[:, idx] = r
        else:
            radi[:, idx] = np.nan

    # from itertools import starmap
    # bs = list(starmap(fdaimage, zip(ys, xs, zs)))

    # for idx, b in zip(idxs, bs):
    #     beta[:, idx] = b

    return beta, radi


def test_fdaimage():
    V = (16, 16)
    N = 50
    Q = 3
    Y = np.random.randn(N, np.prod(V))
    X = np.random.randn(N, Q)
    Z = get_coordinates(V)
    beta, radi = fdaimage(Y, X, Z)
    assert beta.shape == (Q, np.prod(V))
    assert radi.shape == (Q, np.prod(V))


def test_fdaimage3d():
    V = (16, 16, 4)
    N = 50
    Q = 3
    Y = np.random.randn(N, np.prod(V))
    X = np.random.randn(N, Q)
    Z = get_coordinates(V)
    beta, radi = fdaimage3d(Y, X, Z)
    assert beta.shape == (Q, np.prod(V))
    assert radi.shape == (Q, np.prod(V))


def main():
    # test_fdaimage()
    test_fdaimage3d()


if __name__ == '__main__':
    main()
