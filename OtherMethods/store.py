import multiprocessing
from itertools import starmap

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
r = robjects.r
r['source']('store.R')


def stor_grid(Y, X, ss=None, ks=None, ncpu=None):
    if ss is None:
        ss = np.linspace(0.1, 0.9, 9)
    if ks is None:
        ks = list(range(1, 10))
    DD = np.array(Y.shape[1:])
    ntune = len(ks)
    assert len(ss) == ntune
    ss = [np.clip(s * DD, 1, DD-1) for s in ss]
    if ncpu is None:
        ncpu = min(ntune, multiprocessing.cpu_count())
    print('ncpu:', ncpu)
    if ncpu == 1:
        storout = list(starmap(stor, zip([Y]*ntune, [X]*ntune, ss, ks)))
    else:
        with multiprocessing.Pool(processes=ncpu) as pool:
            storout = pool.starmap(stor, zip([Y]*ntune, [X]*ntune, ss, ks))
    bhats = [e[0] for e in storout]
    bics = [e[1] for e in storout]
    idx_best = np.argmin(bics)
    if idx_best in [0, ntune-1]:
        print('Warning: best tuning parameter on boundaries')
    bhat_best = bhats[idx_best]
    return bhat_best, bhats, bics, ss, ks


def stor(Y, X, para, k):

    assert isinstance(Y, np.ndarray)
    assert isinstance(X, np.ndarray)
    assert (
            isinstance(para, np.ndarray)
            or isinstance(para, list)
            or isinstance(para, tuple))
    assert np.issubdtype(type(k), np.integer)
    assert X.shape[0] == Y.shape[0]
    assert X.ndim == 2
    assert np.size(para) == Y.ndim - 1
    assert np.ndim(para) == 1

    Y = [y for y in Y]
    out = r['mystor'](Y, X, para, k)
    bhat = out[list(out.names).index('B.hat')]
    bic = out[list(out.names).index('bic2')]
    bhat = np.array(bhat)
    bic = list(bic)[0]
    bhat = np.transpose(bhat, [-1] + list(range(bhat.ndim-1)))
    return bhat, bic


def test_stor():
    DD = (32, 32, 4)
    n = 50
    p = 3
    Y = np.random.randn(n, *DD)
    X = np.random.randn(n, p)
    para = [2, 2, 2]
    k = 10
    bhat, bic = stor(Y, X, para, k)
    assert bhat.shape == (p,) + DD
    assert bic >= 0


def test_stor_grid():
    DD = (32, 32, 4)
    n = 50
    p = 3
    Y = np.random.randn(n, *DD)
    X = np.random.randn(n, p)
    beta = np.random.randn(*DD[::-1], p)
    Y = (beta @ X.T).T
    Y += np.random.randn(*Y.shape)
    bhat, bhats, bics, ss, ks = stor_grid(Y, X)
    assert bhat.shape == (p,) + DD


def main():
    test_stor()
    test_stor_grid()


if __name__ == '__main__':
    main()
