# Author: Patrick Komiske, 2020
# Edited by: Anthony Badea, 2025

import numpy as np
import os

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

def ibu_wrapper(data, weights_data, mc_reco, mc_gen, weights_mc, binwidth, bins, it=5, output_directory=None):

    # data histogram
    data_hist = np.histogram(data, bins=bins, density=True)[0]
    
    # get the new generator-level histogram
    gen_hist = np.histogram(mc_gen, weights=weights_mc, bins=bins, density=True)[0]

    # make response
    response = np.histogram2d(mc_reco, mc_gen, bins=(bins, bins), weights=weights_mc)[0]
    response /= (response.sum(axis=0) + 10**-50)

    # redo the IBU unfolding with this new prior (genobs_hist_rw)
    phis = ibu(data_hist, response, gen_hist, binwidth, binwidth, it=it)

    # save histograms if outDir
    if output_directory != None:
        np.save(os.path.abspath(os.path.join(output_directory, "data_hist.npy")), data_hist)
        np.save(os.path.abspath(os.path.join(output_directory, "gen_hist.npy")), gen_hist)
        np.save(os.path.abspath(os.path.join(output_directory, "response.npy")), response)

    return phis
    
# statistical uncertainty on the IBU distribution only from uncertainty on the prior
def ibu_bootstrap(data, weights_data, mc_reco, mc_gen, weights_mc, binwidth, bins, it=5, boot="data", nresamples=20):
    
    rephis = []
    for resample in range(nresamples):

        if boot == "data":
            
            # resample the weights
            reweights = np.random.poisson(1, size=len(data)) * weights_data
            # print(weights_data)
            # print(reweights)
            # reweights *= weights_data
        
            # run ibu
            phi = ibu_wrapper(data, reweights, mc_reco, mc_gen, weights_mc, binwidth, bins, it=it)[-1]
    
            # write down the phis
            rephis.append(phi)
            
        elif boot == "mc":
            
            # resample the weights
            reweights = np.random.poisson(1, size=len(mc_gen)) * weights_mc
            # reweights *= weights_mc
        
            # run ibu
            phi = ibu_wrapper(data, weights_data, mc_reco, mc_gen, reweights, binwidth, bins, it=it)[-1]
    
            # write down the phis
            rephis.append(phi)
        
        else:
            print("Not a valid bootstrap.")
            return

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