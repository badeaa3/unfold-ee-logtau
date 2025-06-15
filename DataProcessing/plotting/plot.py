import ROOT
import os
import argparse


# Parse command line arguments
parser = argparse.ArgumentParser(description='Plot histograms from ROOT files')
parser.add_argument("-o", '--outDir', type=str, default='.', help='Output directory for plots (default: current directory)')
args = parser.parse_args()

# output directory
outDir = args.outDir
os.makedirs(outDir, exist_ok=True)

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# custom
from PlotConfig import *
from style import *
import argparse
SetALEPHStyle()

# Open all files
files = [ROOT.TFile.Open(val["file"]) for val in config]

# Get the list of histograms in the first file (assuming all files have the same histograms)
hist_names = [key.GetName() for key in files[0].GetListOfKeys() if key.GetClassName().startswith("TH1")]
# only look for pwflag hists for now and then do the event variables
# hist_names = [hist for hist in hist_names if "objSel1_" in hist]
# hist_names = [hist for hist in hist_names if "pwflag5" in hist]
print(hist_names)

# Loop over histogram names and plot each set of histograms
for hist_name in hist_names:

    # get hist plot config
    pltType = [i for i in plotConfig.keys() if i == hist_name.split("_")[-1]][0]
    pltConfig = plotConfig[pltType]

    # pick up plot type
    if "pwflag" in hist_name:
        pwflag = int(hist_name.split("_")[1].strip("pwflag"))


    # get list of hists
    hist_list = [f.Get(hist_name) for f in files]

    # Skip if any file is missing the histogram
    if not all(hist_list):  
        print(f"Warning: Missing histogram {hist_name} in some files.")

    # Create a new canvas for each histogram
    canvas = ALEPHCanvas(hist_name)
    if pltConfig["SetLogy"]:
        canvas.SetLogy()

    # legend
    legend = ALEPHLegend(loc=pltConfig['legloc'], textsize=0.025) 

    # Normalize histograms if needed (optional)
    goodHists = []
    for i, hist in enumerate(hist_list):
        
        # pwflag don't do no_event_sel?
        if "pwflag" in hist_name and "no_event_sel" in config[i]["file"]:
            goodHists.append(False)
            continue

        # general settings
        max_bin_content = 0
        try:
            hist.SetStats(0)  # Remove stats box
            goodHists.append(True)
        except:
            goodHists.append(False)
            continue
        hist.GetYaxis().SetTitleOffset(1.6)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetMaxDigits(3)
        # rebin and then set divisions
        hist.Rebin(pltConfig["rebin"])
        hist.GetXaxis().SetNdivisions(pltConfig["Ndivisions"])
        hist.GetYaxis().SetNdivisions(pltConfig["Ndivisions"])
        # set style
        hist.SetLineColor(config[i]["color"])
        hist.SetMarkerColor(config[i]["color"])
        hist.SetMarkerStyle(config[i]["MarkerStyle"])
        hist.SetLineWidth(config[i]["LineWidth"])
        hist.SetLineStyle(config[i]["LineStyle"])

        # create the y-axis title 
        XTitle = hist.GetXaxis().GetTitle()
        XBinWidth = hist.GetBinWidth(1)
        temp = XTitle.split("[")
        
        # get units
        if XBinWidth == int(XBinWidth):
            units = f"{int(XBinWidth)}"
        else:
            units = f"{round(XBinWidth, 2)}"
        hist.SetYTitle(pltConfig["YTitle"] + f" / {units}")

        # normalize
        if hist.Integral() > 0:
            hist.Sumw2()
            scale = XBinWidth * hist.Integral()
            hist.Scale(1.0 / scale)

        # save maximum
        max_bin_content = max(max_bin_content, hist.GetMaximum())

    # Set y-axis range to scale times the maximum bin content
    if max_bin_content > 0:
        first_idx = [i for i, val in enumerate(goodHists) if val][0]
        hist_list[first_idx].SetMaximum(pltConfig["scale_max_bin_content"] * max_bin_content)
    
    # draw hists and add legend
    for i, hist in enumerate(hist_list):
        if not goodHists[i]:
            print(f"Skipping: {files[i]}")
            continue
        hist.Draw(config[i]["DrawOption"] if i==0 else config[i]["DrawOption"] + " SAME")
        legend.AddEntry(hist, config[i]["legend"], config[i]["LegendDraw"])

    # draw legend
    legend.Draw()

    # set the ALEPH tag
    # firstspace = 0.05
    # space = 0.04
    # ALEPHLabel(0.2, pltConfig["ALEPHTagTop"])
    # myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 0*space, "#sqrt{s} = 91.2 GeV", size=0.035) # "#sqrt{s} = 91.2 GeV, 45 pb ^{-1}"
    # if "pwflag" in hist_name:
    #     myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 1*space, pwflags[pwflag], size=0.035)
    # elif "objSel" in hist_name:
    #     myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 1*space, f"Object Selection {objSel}", size=0.035)
    header = "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV"
    if "pwflag" in hist_name:
        header += f", {pwflags[pwflag]}"
    # elif "objSel" in hist_name:
    #     header += f", {objSel}"
    myText(0.165, 0.97, header, size=0.03)

    # Save to individual PDF
    pdf_name = os.path.join(outDir, f"{hist_name}.pdf")
    canvas.Print(pdf_name)
    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()

# Close all files
for f in files:
    f.Close()
