"""Numerical translation of paper_fig_v2.m; no MATLAB runtime required."""
import numpy as np
from scipy.optimize import least_squares

SPANS = np.array([1, 5, 11, 21, 31, 61, 121, 181, 241])
LABELS = np.array([1, 5, 10, 20, 30, 60, 120, 180, 240])
KEYS = ['cl', 'ov', 'bc']
NAMES = ['Clear sky', 'Overcast', 'Shallow cumulus']
COLORS = ['#0072bd', '#d95319', '#edb120']

def smooth_matlab(y, span):
    """Odd centered moving mean; endpoints use 1,3,5,... samples.

    Matches smooth(y,span) for the complete/all-missing series in this dataset.
    Reject partial missingness instead of silently guessing MATLAB's fallback.
    """
    y = np.asarray(y, dtype=float)
    if span == 1: return y.copy()
    n = y.shape[0]
    span = min(int(span), n - (n % 2 == 0))
    if span % 2 == 0: span -= 1
    count = np.isfinite(y).sum(axis=0)
    if np.any((count != 0) & (count != n)):
        raise ValueError('Partial missing series requires explicit smoothing policy.')
    i = np.arange(n)
    half = np.minimum(np.minimum(i, n - 1 - i), span // 2)
    lo, hi = i - half, i + half + 1
    total = np.concatenate([np.zeros_like(y[:1]), np.cumsum(np.nan_to_num(y), axis=0)])
    shape = (n,) + (1,) * (y.ndim - 1)
    result = (total[hi] - total[lo]) / (hi-lo).reshape(shape)
    return np.where(count == 0, np.nan, result)

def daily_std(y):
    """MATLAB nanstd(x,0,1), including singleton behavior."""
    count = np.isfinite(y).sum(axis=0)
    mean = np.nansum(y, axis=0) / np.maximum(count, 1)
    var = np.nansum((y-mean)**2, axis=0) / np.maximum(count-1, 1)
    return np.where(count > 0, np.sqrt(var), np.nan)

def mean_nan(y, axis=0):
    count = np.isfinite(y).sum(axis=axis)
    return np.divide(np.nansum(y, axis=axis), count,
                     out=np.full(np.shape(count), np.nan), where=count>0)

def std_groups(smoothed, groups, legacy=False):
    """Return class x span mean and between-day sample SD.

    Legacy keeps trailing day slots, as the MATLAB script actually does.
    """
    buffer = None
    means, errors, effective = [], [], []
    for ids in groups:
        current = smoothed[:, :, 0, ids]  # span, minute, day
        if legacy and buffer is not None and current.shape[-1] < buffer.shape[-1]:
            buffer = buffer.copy()
            buffer[..., :len(ids)] = current
        else:
            buffer = current.copy()
        s = daily_std(buffer.transpose(1, 0, 2))  # span, day
        means.append(mean_nan(s, axis=1))
        errors.append(daily_std(s.T))
        effective.append(buffer.shape[-1])
    return np.array(means), np.array(errors), effective, buffer

def fit_correlations(smoothed, site, groups):
    """36 pair correlations, negative values excluded before day averaging."""
    pairs = np.array([(i,j) for i in range(9) for j in range(i+1,9)])
    distance = (np.linalg.norm(site[:,pairs[:,0]]-site[:,pairs[:,1]],axis=0)*100.2-0.1627)*1000
    x = distance[None,:] / SPANS[:,None] / 60
    results = []
    for ids in groups:
        y = smoothed[:,:,:,ids]
        # Entire missing station-days remain NaN, matching corrcoef default.
        centered = y - np.mean(y,axis=1,keepdims=True)
        normalized = centered / np.sqrt(np.sum(centered**2,axis=1,keepdims=True))
        corr = np.sum(normalized[:,:,pairs[:,0],:] * normalized[:,:,pairs[:,1],:],axis=1)
        valid = np.isfinite(corr) & (corr >= 0)
        corr = np.where(valid,corr,np.nan)
        avg = mean_nan(corr,axis=2)
        keep = np.isfinite(avg)
        xx, yy = x[keep], avg[keep]
        # Fit in original rho space, NOT a log-linear regression.
        fits = [least_squares(lambda p:p[0]*xx**p[1]-yy,[0.8,b],
                              xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=10000)
                for b in [-0.2,0.0,0.2]]
        fit = min(fits,key=lambda f:np.sum(f.fun**2))
        if not fit.success: raise RuntimeError(fit.message)
        results.append(dict(coefficients=fit.x, x=x, rho=avg,
                            valid_days=valid.sum(axis=2), daily=corr,
                            sse=float(np.sum(fit.fun**2)), pairs=pairs, distance=distance))
    return results

def factor(site, coefficients, n, legacy=False):
    """Original analytical formula; legacy includes stale 9x9 cells for N=6.

    Self correlations are omitted to reproduce source methodology; this is
    explicitly documented, not presented as a corrected covariance model.
    """
    size = 9 if legacy else n
    d = np.linalg.norm(site[:,:size,None]-site[:,None,:size],axis=0)*100200
    x = d[:,:,None] / SPANS[None,None,:] / 60
    a,b = coefficients
    with np.errstate(divide='ignore',invalid='ignore'):
        p = np.where(x==0,np.nan,a*x**b)
    sij = np.nansum(p,axis=(0,1))
    sci = np.nansum(p[:,0,:],axis=0)
    return 1+sij/n**2-2*sci/n

def fig6_values(smoothed,site,groups,coefficients,legacy=False):
    values=[]
    # MATLAB x remains 18 days after z=1 and contaminates z=2 classes.
    buffer=None
    for n in [9,6]:
        means=[]
        for ids in groups:
            cur=smoothed[:,:,0,ids]
            if legacy and buffer is not None and len(ids)<buffer.shape[-1]:
                buffer=buffer.copy();buffer[...,:len(ids)]=cur
            else: buffer=cur.copy()
            means.append(mean_nan(daily_std(buffer.transpose(1,0,2)),axis=1))
        f=np.array([factor(site,c,n,legacy) for c in coefficients])
        if np.any(f<0): raise ValueError('Negative modeled variance; review model assumptions.')
        values.append(np.sqrt(f)*np.array(means))
    return np.array(values)
