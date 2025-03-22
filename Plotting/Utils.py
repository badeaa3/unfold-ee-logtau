import glob
import os
import numpy as np
import json
import matplotlib.pyplot as plt

def loadWeights(inPath):
    
    pathList = glob.glob(inPath)
    weights = {}
    for base_path in pathList:

        # pick up the weights
        omnifold_weights_path = glob.glob(os.path.join(base_path,"omnifold_weights*"))
        if len(omnifold_weights_path) > 0:
            omnifold_weights = np.load(omnifold_weights_path[0])
        else:
            print(f"No omnifold weights in {base_path}")
            continue

        # pick up the conf file
        conf = os.path.join(base_path,"conf.json")
        with open(conf) as f:
            conf = json.load(f)  
    
        # omnifold_weights_reco_path = glob.glob(os.path.join(base_path,"omnifold_weights_reco*"))[0]
        # omnifold_weights_reco = np.load(omnifold_weights_reco_path)
 
        starting_weights_path = glob.glob(os.path.join(base_path,"starting_weights*"))[0]
        with np.load(starting_weights_path) as f:
            starting_weights_mc = np.array(f['weights_mc'])

        omnifold_weights *= starting_weights_mc  

        # initialize the list if it doesn't exist, then append to it
        weights.setdefault(conf["SystematicVariation"], []).append((conf["i_ensemble_per_omnifold"], omnifold_weights))


    # sort according to the index of ensemble, i.e. submission order. Then stack in order of SystematicVariation 
    temp = []
    for key in sorted(weights.keys()):
        idx, w = zip(*sorted(weights[key], key=lambda pair: pair[0]))
        temp.append(np.array(w))
    weights = np.stack(temp, 0)

    return weights

def ensembleWeights(weights, N):
    temp = weights[:,:N*int(weights.shape[1]/N)] # extact enough to have full ensembles of N
    temp = temp.reshape(weights.shape[0], -1, N, weights.shape[-1]) # reshape to (nVariation, nEnsemble, nTrainingPerGroup, nEvents)
    # get number of groups and take median
    temp = np.median(temp, axis=2)
    return temp

def getStats(weights):
    nominal = weights[:,0] # nominal is the first
    std = np.std(weights, axis=1)
    std_div_nominal = std/nominal
    return std_div_nominal

# get 2004 aleph measurement
def loadALEPH2004Result(hepData = "/global/homes/b/badea/aleph/data/HEPData-ins636645-v1-Table_54.csv"):
    with open(hepData, 'r') as f:
        
        vals = []
        for row in f:
            if row.startswith('#'):
                continue
                
            if row.startswith('T'):
                print(row.strip())
            else:
                vals.append(row.strip().split(','))
                
    hepdata = np.asarray(vals, dtype=float)

    aleph_bins = 1 - np.append(1.0, hepdata[::-1,1])
    aleph_midbins = (aleph_bins[1:] + aleph_bins[:-1])/2
    aleph_binwidths = aleph_bins[1:] - aleph_bins[:-1]
    aleph_thrust = hepdata[::-1,3]
    aleph_thrust_errs = np.linalg.norm(hepdata[::-1,[-1,-3,-5]], axis=1)
    assert np.all(aleph_bins[1:] == 1 - hepdata[::-1,1]) and np.all(aleph_bins[:-1] == 1 - hepdata[::-1,2])

    log_bins_min = -8 # aleph reported linear binning down to 0 but can't do that for log, so must pick a lower bound. Found that beyond this no more stats
    # aleph_log_bins = np.log(aleph_bins + np.exp(log_bins_min)) # before we used this but this is confusing. For the lowest (1-T) bin with a lower bin edge of 0 we just want to modify the left bin edge but not the right. Instead use the below.
    aleph_log_bins = np.log(aleph_bins)
    aleph_log_bins[0] =  log_bins_min
    aleph_log_midbins = (aleph_log_bins[1:] + aleph_log_bins[:-1])/2
    aleph_log_binwidths = aleph_log_bins[1:] - aleph_log_bins[:-1]
    aleph_log_thrust = aleph_thrust * 0.01 / aleph_log_binwidths
    aleph_log_thrust_errs = aleph_thrust_errs * 0.01 / aleph_log_binwidths

    aleph = {
        "aleph_bins" : aleph_bins,
        "aleph_midbins" : aleph_midbins,
        "aleph_binwidths" : aleph_binwidths,
        "aleph_thrust" : aleph_thrust,
        "aleph_thrust_errs" : aleph_thrust_errs,
        "log_bins_min" : log_bins_min,
        "aleph_log_bins" : aleph_log_bins,
        "aleph_log_midbins" : aleph_log_midbins,
        "aleph_log_binwidths" : aleph_log_binwidths,
        "aleph_log_thrust" : aleph_log_thrust,
        "aleph_log_thrust_errs" : aleph_log_thrust_errs
    }
    return aleph

# function for getting histograms from observable values
def calc_hist(vals, bins=10, weights=None, density=True):
    
    if weights is None:
        weights = np.ones(vals.shape)
    
    # compute histogram
    hist, bins = np.histogram(vals, bins=bins, weights=weights)
    
    # compute which bins the values are in
    digits = np.digitize(vals, bins)

    # compute the errors per bin
    # note that lowest bin value that digitize returns is 1
    # hence the range in the following list comprehension should start at 1
    errs = np.asarray([np.linalg.norm(weights[digits==i]) for i in range(1, len(bins))])

    # handle normalization
    if density:
        binwidths = bins[1:] - bins[:-1]
        density_int = weights.sum() * binwidths # (bins[1] - bins[0])
        hist /= density_int
        errs /= density_int
        
    return hist, errs, bins

def plotThrust(style, inPlots, ratio_denom, epsilon = 1e-10 ):

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(6,6))
    plt.subplots_adjust(hspace=0)

    # plot error bar plots
    for key, plot in inPlots.items():
        if plot["plotType"] == "errorbar":
            # plot nominal
            ax1.errorbar(
                plot["x"], 
                plot["y"], 
                color = plot["color"], 
                label = plot["label"], 
                xerr = plot["xerr"], 
                yerr = plot["yerr"], 
                fmt='o', 
                lw=2, 
                capsize=1.5, 
                capthick=1, 
                markersize=1.5
            )
        elif plot["plotType"] == "stairs":
            ax1.stairs(
                plot["y"], 
                plot["x"], 
                label=plot["label"], 
                color=plot["color"],
                ls=plot["ls"]
            )

    # plot ratios
    for key, plot in inPlots.items():
        ratio_denom_idx = 0 if plot["y"].shape == ratio_denom[0].shape else 1
        plot["ratio_y"] = plot["y"] / (ratio_denom[ratio_denom_idx] + epsilon)
        if plot["plotType"] == "errorbar":
            plot["ratio_yerr"] = plot["yerr"] / (ratio_denom[ratio_denom_idx] + epsilon)
            ax2.errorbar(
                plot["x"], 
                plot["ratio_y"], 
                xerr = plot["xerr"], 
                yerr = plot["ratio_yerr"], 
                fmt = 'o', 
                color = plot["color"],  
                markersize = 4
            )
        elif plot["plotType"] == "stairs":
            ax2.plot(
                (plot["x"][:-1] + plot["x"][1:]) / 2, 
                plot["ratio_y"], 
                color = plot["color"],
                ls = plot["ls"]
            )

    # ratio horizontal line
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.5)  # Adding a horizontal line at y=1 for reference

    # legend
    ax1.legend(loc = style["legend_loc"], 
               bbox_to_anchor = style["legend_bbox"], 
               ncol = style["legend_ncol"],
               fontsize=11)
    
    # axis settings
    ax1.set_ylabel(style["ax1_ylabel"])
    ax1.set_yscale(style["ax1_yscale"])
    ax2.set_xlabel(style["ax2_xlabel"])
    ax2.set_xscale(style["ax2_xscale"])
    ax2.set_ylabel(style["ax2_ylabel"], fontsize=15)

    # set limits
    # ax1.set_ylim(0.2*10**-5, 10**0)
    ax1.set_ylim(style["ax1_ylim"][0], style["ax1_ylim"][1])
    ax2.set_xlim(style["bins"][0], style["bins"][-1])
    ax2.set_ylim(style["ax2_ylim"][0], style["ax2_ylim"][1])

    return fig, (ax1, ax2)

if __name__ == "__main__":

    ensemble = "/global/homes/b/badea/aleph/unfold-ee-logtau/results/training-745de56e/*/model_weights*"
    weights = loadWeights(ensemble)
    print(weights.shape)
    weights = ensembleWeights(weights, N=10)
    print(weights.shape)
    std_div_nominal = getStats(weights)
    print(std_div_nominal.shape, np.min(std_div_nominal), np.max(std_div_nominal))


    systematics = "/global/homes/b/badea/aleph/unfold-ee-logtau/results/training-439c81ff/*/model_weights*"
    weights = loadWeights(systematics)
    print(weights.shape)
    weights = ensembleWeights(weights, N=10)
    print(weights.shape)
    std_div_nominal = getStats(weights)
    print(std_div_nominal.shape, np.min(std_div_nominal), np.max(std_div_nominal))
    