import ROOT
import os
import glob

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# custom
from PlotConfig import *
from style import *
SetALEPHStyle()

# pick up plot configurations
hist_name = "logtau"
pltConfig = plotConfig[hist_name]

def get_ratio(numeratorFilePath, denominatorFilePath, hist_name="logtau"):
    """
    Get the ratio of two histograms from different files.
    """
    numeratorFile = ROOT.TFile.Open(numeratorFilePath)
    denominatorFile = ROOT.TFile.Open(denominatorFilePath)
    if not numeratorFile or not denominatorFile:
        print(f"Error: Could not open files {numeratorFilePath} or {denominatorFilePath}.")
        return None
    
    num_hist = numeratorFile.Get(f"h_{hist_name}")
    denom_hist = denominatorFile.Get(f"h_{hist_name}")

    if not num_hist or not denom_hist:
        print(f"Warning: Histogram h_{hist_name} not found in one of the files.")
        return None

    # Ensure histograms have sum of weights squared enabled for error propagation
    num_hist.Sumw2()
    denom_hist.Sumw2()

    # normalize first
    num_hist.Scale(1.0 / (num_hist.GetBinWidth(1) * num_hist.Integral()))
    denom_hist.Scale(1.0 / (denom_hist.GetBinWidth(1) * denom_hist.Integral()))

    ratio_hist = num_hist.Clone("ratio")
    ratio_hist.Divide(denom_hist)

    ratio_hist.SetStats(0)  # Remove stats box
    ratio_hist.GetYaxis().SetTitleOffset(1.6)
    ratio_hist.GetXaxis().SetTitleOffset(1.2)
    ratio_hist.GetYaxis().SetMaxDigits(3)
    # rebin and then set divisions
    ratio_hist.Rebin(pltConfig["rebin"])
    ratio_hist.GetXaxis().SetNdivisions(pltConfig["Ndivisions"])
    # set style
    color = ROOT.kRed if "alephMCRecoAfterCutPaths" in numeratorFilePath else ROOT.kBlack
    ratio_hist.SetLineColor(color)
    ratio_hist.SetMarkerColor(color)
    ratio_hist.SetMarkerStyle(20)
    ratio_hist.SetMarkerSize(0.8)
    ratio_hist.SetLineWidth(2)
    ratio_hist.SetYTitle("Variation / Nominal") # (Nominal - Variation)/Nominal")
    ratio_hist.GetYaxis().SetRangeUser(0.02, 1.98) # -2.8,2.8)
    
    # Set directory to None to detach from file before closing
    ratio_hist.SetDirectory(0)
    
    # Close files
    numeratorFile.Close()
    denominatorFile.Close()

    return ratio_hist

if __name__ == "__main__":

    variation_labels = {
        "ech10" : "Total Charged Energy Variation", # E_{Ch} #geq 15 #rightarrow 10 GeV",
        "no_event_sel": "No Event Selection",
        "no_neutrals" : "Without Neutral Objects", # "Thrust w/o Neutral Objects",
        "nominal" : "Nominal",
        "ntpc7" : "Chg. Track NTPC Hit Variation", # "Chg. Track N_{TPC} #geq 4 #rightarrow 7",
        "pt04" : "Chg. Track pT Variation", # Chg. Track p_{T} #geq 0.2 #rightarrow 0.4 GeV",
        "with_met": "With MET Object", # Thrust w/ MET",
    }

    # output directory
    outDir = f"SystematicVariationsPlots"
    os.makedirs(outDir, exist_ok=True)

    # get list of files
    dataIdx = 1
    dataFileList = config[dataIdx]["file"].replace("nominal", "*")
    dataFileList = sorted(glob.glob(dataFileList))
    print(dataFileList)

    mcIdx = 3
    mcFileList = config[mcIdx]["file"].replace("nominal", "*")
    mcFileList = sorted(glob.glob(mcFileList))
    print(mcFileList)

    # loop over the others
    for dataFile, mcFile in zip(dataFileList, mcFileList):

        # print the files    
        print(dataFile)
        print(mcFile)

        # get variation name
        variation = dataFile.split("thrust_")[-1].split("_t.root")[0]
        print(variation)

        # get the ratios
        ratio_data = get_ratio(config[dataIdx]["file"], dataFile)
        ratio_mc = get_ratio(config[mcIdx]["file"], mcFile)

        # Create a new canvas for each histogram
        canvas = ALEPHCanvas(f"c_xyz")

        # legend
        legend = ALEPHLegend(loc=[0.5, 0.2, 0.7, 0.28], textsize=0.03) 

        # draw
        ratio_data.Draw("E1")
        ratio_mc.Draw("E1 same")

        # add to legend
        legend.AddEntry(ratio_data, "Uncorrected data", "pl")
        legend.AddEntry(ratio_mc, "Selected archived MC Reco.", "p")

        legend.Draw()

        # Draw a dashed light gray horizontal line at y=1
        line = ROOT.TLine(ratio_data.GetXaxis().GetXmin(), 1, ratio_data.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kGray + 1)
        line.SetLineStyle(2)  # Dashed line
        line.SetLineWidth(3)
        line.Draw()

        header = "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV" + f", {variation_labels[variation]}"
        myText(0.165, 0.97, header, size=0.03)
        
        # Save to individual PDF
        pdf_name = os.path.join(outDir, f"{hist_name}_check_{variation}.pdf")
        canvas.Print(pdf_name)
        print(f"Saved: {pdf_name}")

        # Cleanup
        canvas.Close()
