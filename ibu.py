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
