import glob
import os
import numpy as np
import json
import matplotlib.pyplot as plt
import os
import PyPDF2
from PyPDF2 import PdfReader, PdfWriter, Transformation
import scipy.stats
import matplotlib.patches as patches

def plotBand(ax, x, y, syst_err, bin_widths, color, alpha=0.4):
    """
    Plot points with systematic uncertainty:
      * semi-transparent rectangles per bin
      * a continuous symmetric error band y ± syst_err

    Parameters
    ----------
    ax : matplotlib Axes
        Axis to draw on.
    x : array-like
        Bin centers.
    y : array-like
        Values at bin centers.
    syst_err : array-like
        Symmetric systematic uncertainties (same length as x).
    bin_widths : array-like
        Bin widths (same length as x).
    color : str
        Color for markers and boxes/band.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    syst_err = np.asarray(syst_err)
    bin_widths = np.asarray(bin_widths)

    for xi, yi, err, bw in zip(x, y, syst_err, bin_widths):
        rect = patches.Rectangle(
            (xi - bw/2, yi - err),   # bottom-left corner
            bw,                      # width
            2 * err,                 # height
            facecolor=color,
            alpha=alpha,
            edgecolor='none'
        )
        ax.add_patch(rect)
        
def chi2_ndf(y1, yerr1, y2, yerr2):
    """
    Compute χ²/ndf between two binned distributions with identical binning.

    Parameters
    ----------
    y1, yerr1 : array-like
        Values and uncertainties for dataset 1 (e.g. ALEPH).
    y2, yerr2 : array-like
        Values and uncertainties for dataset 2 (e.g. Unifold).

    Returns
    -------
    chi2, ndf, chi2_per_ndf : tuple of floats
        Total chi-square, number of degrees of freedom, and χ²/ndf.
    """
    y1, yerr1 = np.asarray(y1), np.asarray(yerr1)
    y2, yerr2 = np.asarray(y2), np.asarray(yerr2)

    # Combine uncertainties in quadrature
    sigma2 = yerr1**2 + yerr2**2

    # Avoid division by zero
    valid = sigma2 > 0
    chi2 = np.sum(((y1[valid] - y2[valid])**2) / sigma2[valid])
    ndf = np.count_nonzero(valid)
    p_value = scipy.stats.chi2.sf(chi2, ndf)  # survival function (1 - CDF)
    return chi2, ndf, chi2 / ndf, p_value

def pull_rms(y1, yerr1, y2, yerr2):
    dy = np.asarray(y1) - np.asarray(y2)
    sigma = np.sqrt(yerr1**2 + yerr2**2)
    pulls = dy / sigma
    return np.std(pulls), pulls

def global_shift(y1, yerr1, y2, yerr2):
    y1, y2 = np.asarray(y1), np.asarray(y2)
    yerr1, yerr2 = np.asarray(yerr1), np.asarray(yerr2)

    dy = y1 - y2
    sigma = np.sqrt(yerr1**2 + yerr2**2)
    w = 1 / sigma**2

    mean_shift = np.sum(w * dy) / np.sum(w)
    shift_sigma = np.sqrt(1 / np.sum(w))
    Z = mean_shift / shift_sigma

    return mean_shift, shift_sigma, Z
    
def mc_bin_quality(
    data_counts: np.ndarray,
    data_err: np.ndarray,
    mc_counts: np.ndarray,
    mc_err: np.ndarray,
    rel_mc_unc_max: float = 0.20,
    pull_abs_max: float = 3.0
) -> np.ndarray:
    """
    Return a boolean mask (same shape as input arrays) where each bin is True
    if the MC is statistically reliable AND compatible with data.

    Parameters
    ----------
    data_counts : array
        Data histogram bin contents. 
    data_err : array
        1σ statistical errors on the data bins.
    mc_counts : array
        MC histogram bin contents.
    mc_err : array
        1σ statistical errors on the MC bins (e.g. sqrt(sum w^2)).
    rel_mc_unc_max : float
        Maximum allowed relative MC uncertainty (default 0.20 = 20%).
    pull_abs_max : float
        Maximum allowed absolute pull (default 3σ).
    """
    data_counts = np.asarray(data_counts, dtype=float)
    data_err    = np.asarray(data_err, dtype=float)
    mc_counts   = np.asarray(mc_counts, dtype=float)
    mc_err      = np.asarray(mc_err, dtype=float)

    # Combined statistical uncertainty per bin
    sigma = np.sqrt(data_err**2 + mc_err**2)

    # Relative MC uncertainty
    rel_mc_unc = np.divide(mc_err, mc_counts,
                           out=np.zeros_like(mc_err),
                           where=mc_counts > 0)

    # Pull (z-score) for each bin
    pull = np.divide(mc_counts - data_counts, sigma,
                     out=np.zeros_like(mc_err),
                     where=sigma > 0)

    # Criteria
    good_rel_unc = rel_mc_unc <= rel_mc_unc_max
    good_pull    = np.abs(pull) <= pull_abs_max

    return good_rel_unc & good_pull
    
def watermark(
    in_file, # input file name
    out_file, # output file name
    scale=0.12, tx=44, ty=251,
    logo_fpath='./ee-logo.pdf',
    **kwargs
):

    # ensure out_plots_dir exists
    # os.makedirs(out_plots_dir, exist_ok=True)
    
    # open files for bare plot and the logo
    bare_plot = open(in_file, 'rb')
    logo = open(logo_fpath, 'rb')
    
    # extract pdf pages for bare plot and the logo
    plot_page = PyPDF2.PdfFileReader(bare_plot).getPage(0)
    logo_page = PyPDF2.PdfFileReader(logo).getPage(0)
    
    # add the watermark
    plot_page.mergeScaledTranslatedPage(logo_page, scale, tx, ty, expand=True)
    
    # create a pdf writer for the new plot
    out_plot_pdf = PyPDF2.PdfFileWriter()
    out_plot_pdf.addPage(plot_page)
    
    # write new plot to PDF
    out_plot = open(out_file, 'wb')
    out_plot_pdf.write(out_plot)
    
    # close all files
    bare_plot.close(); logo.close(); out_plot.close()
    
def loadWeightPaths(fileList): #file_pattern):
    # fileList = sorted(glob.glob(file_pattern))
    d = {}

    for file_path in fileList:
        with open(file_path) as f:
            conf = json.load(f)

        job_type = conf["job_type"]
        # data_key = conf["data"].split("_thrust_")[-1].split("_t.root")[0]
        data_key = conf["reco"].split("_thrust_")[-1].split("_t.root")[0]
        weight_path = os.path.dirname(file_path)

        if job_type not in d:
            d[job_type] = {}

        if data_key not in d[job_type]:
            d[job_type][data_key] = []

        d[job_type][data_key].append(weight_path)

    return d

def loadWeights(inPath):
    
    if type(inPath) == str:
        pathList = sorted(glob.glob(inPath))
    if type(inPath) == list:
        pathList = inPath

    weights = {}
    for base_path in pathList:

        # pick up the omnifold weights
        omnifold_weights_path = os.path.join(base_path, "omnifold_weights.npy")
        if os.path.isfile(omnifold_weights_path):
            omnifold_weights = np.load(omnifold_weights_path)
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
        if "SystematicVariation" in conf.keys():
            weights.setdefault(conf["SystematicVariation"], []).append((conf["i_ensemble_per_omnifold"], omnifold_weights))
        else:
            weights.setdefault(conf["data"], []).append((conf["i_ensemble_per_omnifold"], omnifold_weights))

    # sort according to the index of ensemble, i.e. submission order. Then stack in order of SystematicVariation 
    temp = []
    for key in sorted(weights.keys()):
        idx, w = zip(*sorted(weights[key], key=lambda pair: pair[0]))
        temp.append(np.array(w))
    weights = np.stack(temp, 0)

    return weights


def ensemblePredsThenGetWeight(w, N):
    f = w/(1+w) # go back to raw NN predictions from w = f(x)/(1-x) -> f(x) = w/(1+w)
    f = f[:,:N*int(f.shape[1]/N)] # extact enough to have full ensembles of N
    f = f.reshape(f.shape[0], -1, N, f.shape[-1]) # reshape to (nVariation, nEnsemble, nTrainingPerGroup, nEvents)
    # get number of groups and take median
    f = np.mean(f, axis=2)
    w = f/(1-f) # go back to weights
    return w

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
    aleph_thrust_errs_individual = hepdata[::-1,[4,6,8]]
    aleph_thrust_errs = np.linalg.norm(aleph_thrust_errs_individual, axis=1)
    assert np.all(aleph_bins[1:] == 1 - hepdata[::-1,1]) and np.all(aleph_bins[:-1] == 1 - hepdata[::-1,2])
    
    log_bins_min = -6 # aleph reported linear binning down to 0 but can't do that for log, so must pick a lower bound. Found that beyond this no more stats
    # aleph_log_bins = np.log(aleph_bins + np.exp(log_bins_min)) # before we used this but this is confusing. For the lowest (1-T) bin with a lower bin edge of 0 we just want to modify the left bin edge but not the right. Instead use the below.
    aleph_log_bins = np.log(aleph_bins)
    aleph_log_bins[0] =  log_bins_min
    aleph_log_midbins = (aleph_log_bins[1:] + aleph_log_bins[:-1])/2
    aleph_log_binwidths = aleph_log_bins[1:] - aleph_log_bins[:-1]
    aleph_log_thrust = aleph_thrust * aleph_binwidths[0] / aleph_log_binwidths # because the reported value were scaled by 1/bin width
    aleph_log_thrust_errs = aleph_thrust_errs * aleph_binwidths[0] / aleph_log_binwidths
    aleph_log_thrust_errs_individual = aleph_thrust_errs_individual * aleph_binwidths[0] / np.repeat(np.expand_dims(aleph_log_binwidths,1), aleph_thrust_errs_individual.shape[1], 1) 

    aleph = {
        "aleph_bins" : aleph_bins,
        "aleph_midbins" : aleph_midbins,
        "aleph_binwidths" : aleph_binwidths,
        "aleph_thrust" : aleph_thrust,
        "aleph_thrust_errs" : aleph_thrust_errs,
        "aleph_thrust_errs_individual" : aleph_thrust_errs_individual,
        "log_bins_min" : log_bins_min,
        "aleph_log_bins" : aleph_log_bins,
        "aleph_log_midbins" : aleph_log_midbins,
        "aleph_log_binwidths" : aleph_log_binwidths,
        "aleph_log_thrust" : aleph_log_thrust,
        "aleph_log_thrust_errs" : aleph_log_thrust_errs,
        "aleph_log_thrust_errs_individual" : aleph_log_thrust_errs_individual
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

def ratio_with_uncertainty(A, B, A_err=None, B_err=None):
    ratio = np.zeros_like(A, dtype=np.float64)
    ratio_err = np.zeros_like(A, dtype=np.float64)

    # Avoid division by zero
    nonzero = B != 0
    # ratio[nonzero] = A[nonzero] / B[nonzero]
    ratio = A / B
    
    if type(A_err) == np.ndarray and type(B_err) == np.ndarray:
        # ratio_err[nonzero] = ratio[nonzero] * np.sqrt(
        #     (A_err[nonzero] / A[nonzero])**2 + (B_err[nonzero] / B[nonzero])**2
        # )
        ratio_err = ratio * np.sqrt(
            (A_err / A)**2 + (B_err / B)**2
        )
    
    # if there's a divide by zero then set to positive infinity
    ratio = np.nan_to_num(ratio, posinf=np.inf)
    ratio_err = np.nan_to_num(ratio_err, posinf=np.inf)

    return ratio, ratio_err

def plotThrust(style, inPlots, ratio_denom, epsilon = 1e-10, header = r"ALEPH e$^{+}$e$^{-}$, $\sqrt{s}$ = 91.2 GeV"):

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3.5, 1]}, figsize=style["figsize"])
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
                fmt=plot.get("fmt", 'o'), 
                lw=plot.get("lw", 2), 
                capsize=plot.get("capsize", 3), 
                capthick=1, 
                markersize=plot.get("markersize", 1.5),
                alpha=plot.get("alpha", 1),
                markerfacecolor=plot.get("markerfacecolor", "auto")
            )
        elif plot["plotType"] == "stairs":
            ax1.stairs(
                plot["y"], 
                plot["x"], 
                label=plot["label"], 
                color=plot["color"],
                ls=plot["ls"],
                lw=plot.get("lw", 2),
                alpha=plot.get("alpha",1)
            )

    # plot ratios
    for key, plot in inPlots.items():
        if "noratio" in plot.keys() and plot["noratio"] == True:
            print(f"No ratio plot for {key}")
            continue
        ratio_denom_idx = 0 if plot["y"].shape == ratio_denom[0][0].shape else 1
        if "ratioidx" in plot.keys():
            ratio_denom_idx = plot["ratioidx"]
        # plot["ratio_y"] = plot["y"] / (ratio_denom[ratio_denom_idx] + epsilon)
        # get ratio
        if "yerr" in plot.keys():
            plot["ratio_y"], plot["ratio_yerr"] = ratio_with_uncertainty(plot["y"], ratio_denom[ratio_denom_idx][0], plot["yerr"], ratio_denom[ratio_denom_idx][1])
        else:
            plot["ratio_y"], _ = ratio_with_uncertainty(plot["y"], ratio_denom[ratio_denom_idx][0])
        # plot
        if plot["plotType"] == "errorbar":
            # plot["ratio_yerr"] = plot["yerr"] / (ratio_denom[ratio_denom_idx] + epsilon)
            ax2.errorbar(
                plot["x"], 
                plot["ratio_y"], 
                xerr = plot["xerr"], 
                yerr = plot["ratio_yerr"], 
                fmt = 'o', 
                color = plot["color"],  
                lw=plot.get("lw", 2), 
                capsize=plot.get("capsize", 3), 
                capthick=1, 
                markersize=plot.get("markersize", 1.5),
                alpha=plot.get("alpha", 1)
            )
        elif plot["plotType"] == "stairs":
            ax2.plot(
                (plot["x"][:-1] + plot["x"][1:]) / 2, 
                plot["ratio_y"], 
                color = plot["color"],
                ls = plot["ls"],
                lw = plot.get("lw", 2)
            )

    # ratio horizontal line
    ax2.axhline(y=1, color='black', linestyle='--', alpha=1, lw=1)  # Adding a horizontal line at y=1 for reference

    # legend
    ax1.legend(loc = style["legend_loc"], 
               bbox_to_anchor = style["legend_bbox"], 
               ncol = style["legend_ncol"],
               fontsize = style["legend_fontsize"],
               handletextpad=0.7,
               handlelength=0.8, 
               # handleheight=0.5, 
               # labelspacing=0.5, 
               # columnspacing=1.0
              )
       
    # axis settings
    ax1.set_ylabel(style["ax1_ylabel"], fontsize=style.get("ax1_ylabel_fs", 18), labelpad=8)
    ax1.set_yscale(style["ax1_yscale"])
    ax2.set_xlabel(style["ax2_xlabel"], fontsize=style.get("ax2_xlabel_fs", 18), labelpad=8)
    ax2.set_xscale(style["ax2_xscale"])
    ax2.set_ylabel(style["ax2_ylabel"], fontsize=style.get("ax2_ylabel_fs", 18))

    # set limits
    # ax1.set_ylim(0.2*10**-5, 10**0)
    if "ax1_ylim" in style.keys() and style["ax1_ylim"] is not None:
        ax1.set_ylim(style["ax1_ylim"][0], style["ax1_ylim"][1])
    if "ax2_xlim" in style.keys() and style["ax2_xlim"] is not None:
        ax2.set_xlim(style["ax2_xlim"][0], style["ax2_xlim"][1])
    else:
        ax2.set_xlim(style["bins"][0], style["bins"][-1])
    ax2.set_ylim(style["ax2_ylim"][0], style["ax2_ylim"][1])

    ax1.tick_params(axis='both', which='major', labelsize=style.get("ax1_tick_ls", 15))
    ax2.tick_params(axis='both', which='major', labelsize=style.get("ax2_tick_ls", 15))

    # top text
    ax1.text(0, 1.01, header, transform=ax1.transAxes, ha='left', va='bottom', fontsize=style["header_fontsize"])

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
    