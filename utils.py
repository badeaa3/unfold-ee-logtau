import numpy as np
import matplotlib.pyplot as plt

binning=np.linspace(0,0.5,20)

def HistRoutine(feed_dict,xlabel='',ylabel='',reference_name='Geant4',logy=False,binning=None,label_loc='best',plot_ratio=True,weights=None,uncertainty=None):
    assert reference_name in feed_dict.keys(), "ERROR: Don't know the reference distribution"
    
    fig,gs = SetGrid(ratio=plot_ratio) 
    ax0 = plt.subplot(gs[0])
    if plot_ratio:
        plt.xticks(fontsize=0)
        ax1 = plt.subplot(gs[1],sharex=ax0)

    
    if binning is None:
        binning = np.linspace(np.quantile(feed_dict[reference_name],0.0),np.quantile(feed_dict[reference_name],1),10)
        
    xaxis = [(binning[i] + binning[i+1])/2.0 for i in range(len(binning)-1)]
    reference_hist,_ = np.histogram(feed_dict[reference_name],bins=binning,density=True)
    for ip,plot in enumerate(feed_dict.keys()):
        if weights is not None:
            dist,_,_=ax0.hist(feed_dict[plot],bins=binning,label=plot,linestyle=line_style[plot],color=colors[plot],density=True,histtype="step",weights=weights[plot])
        else:
            dist,_,_=ax0.hist(feed_dict[plot],bins=binning,label=plot,linestyle=line_style[plot],color=colors[plot],density=True,histtype="step")
        
        if plot_ratio:
            if reference_name!=plot:
                ratio = np.ma.divide(dist,reference_hist).filled(0)
                ax1.plot(xaxis,ratio,color=colors[plot],marker='o',ms=10,lw=0,markerfacecolor='none',markeredgewidth=3)
                if uncertainty is not None:
                    for ibin in range(len(binning)-1):
                        xup = binning[ibin+1]
                        xlow = binning[ibin]
                        ax1.fill_between(np.array([xlow,xup]),
                                         uncertainty[ibin],-uncertainty[ibin], alpha=0.3,color='k')    
    if logy:
        ax0.set_yscale('log')

    ax0.legend(loc=label_loc,fontsize=16,ncol=1)        
    if plot_ratio:
        FormatFig(xlabel = "", ylabel = ylabel,ax0=ax0) 
        plt.ylabel('Ratio to Gen.')
        plt.axhline(y=1.0, color='r', linestyle='-',linewidth=1)
        # plt.axhline(y=10, color='r', linestyle='--',linewidth=1)
        # plt.axhline(y=-10, color='r', linestyle='--',linewidth=1)
        plt.ylim([0.8,1.2])
        plt.xlabel(xlabel)
    else:
        FormatFig(xlabel = xlabel, ylabel = ylabel,ax0=ax0) 
        
    return fig,ax0

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
        # print(binwidths)
        density_int = weights.sum() * binwidths # (bins[1] - bins[0])
        hist /= density_int
        errs /= density_int
        
    return hist, errs, bins

# function to propagate errors on division
def propagate_error_division(A, deltaA, B, deltaB):
    """
    Compute the propagated error for division of two quantities A and B with their uncertainties.

    Parameters:
    A (numpy array): Array of values for A.
    deltaA (numpy array): Array of uncertainties for A.
    B (numpy array): Array of values for B.
    deltaB (numpy array): Array of uncertainties for B.

    Returns:
    numpy array: Array of tuples containing (C, deltaC) for corresponding A, deltaA, B, deltaB.
    """
    C = A / B
    deltaC = C * np.sqrt((deltaA / A)**2 + (deltaB / B)**2)
    return deltaC