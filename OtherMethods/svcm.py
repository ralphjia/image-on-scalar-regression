import numpy as np
from design import adj_variance

def svcm_dc(X, A, chunks):
    from scipy.stats import t as tdist
    N, Q = A.shape
    V = X.shape[1:]
    v = V[-1]
    kk = np.round(np.linspace(0, 1, (chunks-1)*2) * v).astype(int)
    kk = np.vstack([kk[:chunks], kk[(chunks-2):]]).T
    print(kk)
    # kk = [[0, v//2], [v//4, v//4*3], [v//2, v]]
    betas = np.full((len(kk), Q) + V, np.nan)
    tscores = np.full_like(betas, np.nan)
    # omegas = np.full((len(kk),) + X.shape, np.nan)
    for i, k in enumerate(kk):
        print(i)
        x = X[..., k[0]:k[1]]
        be, _, _, tsc = svcm(x, A, fromfile=False)
        betas[i, ..., k[0]:k[1]] = be
        tscores[i, ..., k[0]:k[1]] = tsc
        # omega[i, ..., k[0]:k[1]] = om
    beta = np.nanmean(betas, axis=0)
    tscore = np.nanmean(tscores, axis=0)
    pval = (1 - tdist.cdf(np.abs(tscore), df=N-Q)) * 2
    # omega = np.nanmean(omegas, axis=0)
    return beta, pval

def svcm(X, A, fromfile=True):
    if X.ndim != 4:
        raise ValueError('X must be a 4-D array')
    if A.ndim != 2:
        raise ValueError('A must be a 2-D array')
    if X.shape[0] != A.shape[0]:
        raise ValueError('X and A must have equal axis-0 length')

    X = np.transpose(X, [1, 2, 3, 0])

    if fromfile:
        from scipy.io import savemat, loadmat
        import subprocess
        randid = np.random.choice(int(1e6))
        matin = 'tmp/svcmin_{}.mat'.format(randid)
        matout = 'tmp/svcmout_{}.mat'.format(randid)
        savemat(matin, mdict={'imgData':X, 'xMatrix':A})
        subprocess.run(['./svcm.sh', matin, matout])
        mats = loadmat(matout)
        beta = mats['beta']
        pval = mats['pval']
        omega = mats['omega']
        tscore = mats['tscore']
    else:
        import matlab.engine
        eng = matlab.engine.start_matlab()
        nargout = 4
        X = matlab.double(X.tolist())
        A = matlab.double(A.tolist())
        beta, pval, omega, tscore = eng.svcm(X, A, nargout=nargout)
        beta = np.array(beta)
        pval = np.array(pval)
        omega = np.array(omega)
        tscore = np.array(tscore)

    beta = np.transpose(beta, [3, 0, 1, 2])
    pval = np.transpose(pval, [3, 0, 1, 2])
    omega = np.transpose(omega, [3, 0, 1, 2])
    tscore = np.transpose(tscore, [3, 0, 1, 2])

    return beta, pval, omega, tscore

def svcm_gendata(n, beta_stn=1.0, omega_stn=1.0, omega_itv=1.0, noise='gauss', scale=1.0, seed=None):
    assert isinstance(n, int)
    assert n > 0

    import matlab.engine
    eng = matlab.engine.start_matlab()
    nargout = 4
    if noise == 'gauss':
        noise = 1
    elif noise == 'chisq':
        noise = 2
    else:
        raise ValueError('Noise design not recognized')

    if seed is None:
        seed = 'shuffle'

    X, A, beta, omega = eng.svcm_gendata(n, noise, seed, nargout=nargout)

    X = np.array(X)
    A = np.array(A)
    beta = np.array(beta)
    omega = np.array(omega)
    X = np.transpose(X, [3, 0, 1, 2])
    beta = np.transpose(beta, [3, 0, 1, 2])
    omega = np.transpose(omega, [3, 0, 1, 2])
    V = X.shape[1:]
    VV = np.prod(V)
    X = X.reshape((-1, VV))
    beta = beta.reshape((-1, VV))
    omega = omega.reshape((-1, VV))

    noise = X - A @ beta - omega

    beta, omega, noise = adj_variance(A, beta, omega, noise, beta_stn, omega_stn, omega_itv, scale)

    X = A @ beta + omega + noise

    return X, A, beta, omega
