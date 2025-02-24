import ROOT
import os

from style import *
SetALEPHStyle()


# config
config = []
config.append({
    "file" : "../LEP1Data1994_recons_aftercut-MERGED_thrust.root",
    "legend" : "Data 1994"
})
config.append({
    "file" : "../alephMCRecoAfterCutPaths_1994_thrust.root",
    "legend" : "Pythia 6.1"
})

# individual plot settings
plotConfig = {}
plotConfig["cosTheta"] = {"SetLogy" : True, "rebin" : 5}
plotConfig["d0"] = {"SetLogy" : True, "rebin" : 10}
plotConfig["energy"] = {"SetLogy" : True, "rebin" : 5}
plotConfig["mass"] = {"SetLogy" : False, "rebin" : 5}
plotConfig["ntpc"] = {"SetLogy" : False, "rebin" : 1}
plotConfig["phi"] = {"SetLogy" : False, "rebin" : 5}
plotConfig["pmag"] = {"SetLogy" : True, "rebin" : 1}
plotConfig["pt"] = {"SetLogy" : True, "rebin" : 6}
plotConfig["z0"] = {"SetLogy" : True, "rebin" : 2}

# list of names to set log
SetLogy = ["cosTheta", "d0", "energy", "pmag", "pt", "z0"]

# pwflag
pwflags = ["Charged Tracks", "Charged Leptons 1", "Charged Leptons 2", "V0", "Photons", "Neutral Hadrons"]

# Open all files
files = [ROOT.TFile.Open(val["file"]) for val in config]

# Get the list of histograms in the first file (assuming all files have the same histograms)
hist_names = [key.GetName() for key in files[0].GetListOfKeys() if key.GetClassName().startswith("TH1")] # [:9]
# only look for pwflag hists for now and then do the event variables
hist_names = [hist for hist in hist_names if "pwflag" in hist]
print(hist_names)

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# # rebin settings
# rebinning = [5, 10, 5, 5, 1, 6, 2, 5, 2]

# Loop over histogram names and plot each set of histograms
for hist_name in hist_names:

    hist_list = [f.Get(hist_name) for f in files]
    print(hist_list)

    # get hist plot config
    pltType = [i for i in plotConfig.keys() if i in hist_name][0]

    # pick up pwflag
    pwflag = int(hist_name.split("_")[2].strip("pwflag"))
    outDir = f"pwflag{pwflag}"
    os.makedirs(outDir, exist_ok=True)

    if not all(hist_list):  # Skip if any file is missing the histogram
        print(f"Warning: Missing histogram {hist_name} in some files.")
        continue

    # Create a new canvas for each histogram
    canvas = ALEPHCanvas(hist_name)
    if any([i in hist_name for i in SetLogy]):
        canvas.SetLogy()

    # legend
    legend = ROOT.TLegend(0.68, 0.8, 0.88, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)

    # mc index
    mcI = 0

    # rebin factor
    rebin_factor = plotConfig[pltType]["rebin"]
    
    # Normalize histograms if needed (optional)
    for i, hist in enumerate(hist_list):
        
        # general settings
        max_bin_content = 0
        hist.SetStats(0)  # Remove stats box
        hist.GetYaxis().SetTitleOffset(1.6)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetMaxDigits(3)
        
        # rebin and then set divisions
        # rebin_factor = 1 if hist.GetBinWidth(0) == 1 else 2
        hist.Rebin(rebin_factor)
        # nBins = hist.GetNbinsX()
        # hist.GetXaxis().SetNdivisions(nBins, 1, 0)
        hist.GetXaxis().SetNdivisions(310)

        # create the y-axis title 
        XTitle = hist.GetXaxis().GetTitle()
        XBinWidth = hist.GetBinWidth(0)
        temp = XTitle.split("[")
        # get the name
        name = temp[0]
        if name == "Energy" : name = "E"
        if name == "Mass" : name = "m"
        # get the units
        units = f" [pb/{round(XBinWidth,2)}" + ("]" if len(temp) == 1 else " "+temp[1])
        # what we are plotting reduces to the ytitle in terms of sigma
        # 1/N dN/dX
        # = 1/(sigma * L) d(sigma*L)/dX
        # = 1/sigma d(sigma)/dX
        YTitle = "1/#sigma d#sigma/d" + temp[0] + units
        hist.SetYTitle(YTitle)

        # normalize
        if hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

        # save maximum
        max_bin_content = max(max_bin_content, hist.GetMaximum())
        
        # data set specific
        if "LEP1" in config[i]["file"]:
            # Data: Black dots with error bars
            hist.SetMarkerStyle(20)  # Solid circle
            hist.SetMarkerColor(ROOT.kBlack)
            hist.SetLineColor(ROOT.kBlack)
        elif "MC" in config[i]["file"]:
            # MC: Solid line histogram with color
            colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange]
            hist.SetLineColor(colors[mcI])
            hist.SetLineWidth(2)
            mcI += 1
        else:
            print(f"Warning: {fname} does not match 'data' or 'MC'. Defaulting to black.")
            hist.SetLineColor(ROOT.kBlack)

    # Set y-axis range to 1.5 times the maximum bin content
    if max_bin_content > 0:
        scale = 5 if any([i in hist_name for i in SetLogy]) else 1.5
        hist_list[0].SetMaximum(scale * max_bin_content)
    
    # draw hists
    for i, hist in enumerate(hist_list):
        draw_option = "PE" if "LEP1" in config[i]["file"] else "HIST"
        hist.Draw(draw_option if i==0 else draw_option + " SAME")
        # Add to legend
        legend.AddEntry(hist, config[i]["legend"], "p" if "LEP1" in config[i]["file"] else "l")

    legend.Draw()

    # By default tag in upper left. If in this list then put in lower left
    ALEPHTagLowerLeft = ["energy", "pmag", "pt"]
    top = 0.3 if any([i in hist_name for i in ALEPHTagLowerLeft]) else 0.87
    firstspace = 0.03
    space = 0.04
    ALEPHLabel(0.2, top)
    myText(0.2, top - firstspace - 0*space, "#sqrt{s} = 91.2 GeV, 44 pb ^{-1}", size=0.035)
    myText(0.2, top - firstspace - 1*space, pwflags[pwflag] + f" (pwflag = {pwflag})", size=0.035)

    # Save to individual PDF
    pdf_name = os.path.join(outDir, f"{hist_name}.pdf")
    canvas.Print(pdf_name)

    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()

# Close all files
for f in files:
    f.Close()
