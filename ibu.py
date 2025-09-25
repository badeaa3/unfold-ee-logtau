# Author: Patrick Komiske, 2020
# Edited by: Anthony Badea, 2025

import numpy as np

# Iterative Bayesian Unfolding, requires uniform binning (but det and mc can be different)
# data: measured histogram
# r: response matrix
# init: the prior
# it: number of iterations
def ibu(data, r, init, det_binwidth, mc_binwidth, it=10):
    
    # initialize the truth distribution to the prior
    phis = [init]
    
    # iterate the procedure
    for i in range(it):
        
        # update the estimate for the matrix m
        m = r * phis[-1]
        m /= (m.sum(axis=1)[:,np.newaxis] + 10**-50)

        # update the estimate for the truth distribution
        # the factors of binwidth show up here to change probabilities into probability densities
        phis.append(np.dot(m.T, data)*det_binwidth/mc_binwidth)
        
    return phis

# statistical uncertainty on the IBU distribution only from uncertainty on the prior
def ibu_unc(data_hist, response, mc_gen, binwidth_det, bins_mc, binwidth_mc, it=5, nresamples=20):
    
    rephis = []
    for resample in range(nresamples):
        
        # resample the weights
        reweights = np.random.poisson(1, size=len(mc_gen))

        # get the new generator-level histogram
        genobs_hist_rw = np.histogram(mc_gen, weights=reweights, bins=bins_mc, density=True)[0]

        # redo the IBU unfolding with this new prior (genobs_hist_rw)
        phi = ibu(data_hist, response, genobs_hist_rw, binwidth_det, binwidth_mc, it=it)[-1]

        # write down the phis
        rephis.append(phi)

    # return the standard deviation, bin-by-bin, as the uncertainty
    return np.std(np.asarray(rephis), axis=0)


def ibu_covariance(data, r, init, det_bw, mc_bw, cov_data, it=10):
    """
    Compute full covariance of the unfolded spectrum
    due to statistical (Poisson) uncertainties in 'data'.
    """
    data = np.asarray(data, dtype=float)
    n_det = len(data)
    n_mc  = len(init)

    # Unfold once to get the central value
    phi_final = ibu(data, r, init, det_bw, mc_bw, it)[-1]

    # Poisson covariance of the measured data
    # cov_data = np.diag(data)

    # Finite-difference Jacobian: J[i_truth, j_meas]
    J = np.zeros((n_mc, n_det))
    eps = 1e-6
    for j in range(n_det):
        # Small perturbation in bin j
        delta = np.zeros_like(data)
        delta[j] = max(np.sqrt(data[j]), 1.0) * 1e-4  # scale step to stats
        phi_up   = ibu(data + delta, r, init, det_bw, mc_bw, it)[-1]
        phi_down = ibu(np.clip(data - delta, 0, None), r, init, det_bw, mc_bw, it)[-1]
        J[:, j]  = (phi_up - phi_down) / (2 * delta[j])

    # Propagate: Cov_final = J Cov_data J^T
    cov_final = J @ cov_data @ J.T
    return phi_final, cov_final