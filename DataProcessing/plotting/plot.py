import ROOT
import os

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# custom
from PlotConfig import *
from style import *
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
        pwflag = int(hist_name.split("_")[2].strip("pwflag"))
        outDir = f"pwflag{pwflag}"
        os.makedirs(outDir, exist_ok=True)
    elif "objSel" in hist_name:
        objSel = int(hist_name.split("_")[2].strip("objSel"))
        outDir = f"objSel{objSel}"
        os.makedirs(outDir, exist_ok=True)

    # get list of hists
    hist_list = [f.Get(hist_name) for f in files]
    # Skip if any file is missing the histogram
    if not all(hist_list):  
        print(f"Warning: Missing histogram {hist_name} in some files.")
        continue

    # Create a new canvas for each histogram
    canvas = ALEPHCanvas(hist_name)
    if pltConfig["SetLogy"]:
        canvas.SetLogy()

    # legend
    legend = ALEPHLegend() 

    # Normalize histograms if needed (optional)
    for i, hist in enumerate(hist_list):
        
        # general settings
        max_bin_content = 0
        hist.SetStats(0)  # Remove stats box
        hist.GetYaxis().SetTitleOffset(1.6)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetMaxDigits(3)
        # rebin and then set divisions
        hist.Rebin(pltConfig["rebin"])
        hist.GetXaxis().SetNdivisions(pltConfig["Ndivisions"])
        # set style
        hist.SetLineColor(config[i]["color"])
        hist.SetMarkerColor(config[i]["color"])
        hist.SetMarkerStyle(config[i]["MarkerStyle"])
        hist.SetLineWidth(config[i]["LineWidth"])

        # create the y-axis title 
        XTitle = hist.GetXaxis().GetTitle()
        XBinWidth = hist.GetBinWidth(0)
        temp = XTitle.split("[")
        # get the units
        units = "[pb/" 
        if XBinWidth == 1 and len(temp) != 1:
            units += temp[1]
        elif XBinWidth == 1 and len(temp) == 1:
            units += f"{int(XBinWidth)}]" # unit
        elif XBinWidth != 1 and len(temp) != 1:
            units += f"{round(XBinWidth,2)} {temp[1]}"
        else:
            units += f"{round(XBinWidth,2)}]"
        hist.SetYTitle(pltConfig["YTitle"] + " " + units)

        # normalize
        if hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

        # save maximum
        max_bin_content = max(max_bin_content, hist.GetMaximum())

    # Set y-axis range to scale times the maximum bin content
    if max_bin_content > 0:
        hist_list[0].SetMaximum(pltConfig["scale_max_bin_content"] * max_bin_content)
    
    # draw hists and add legend
    for i, hist in enumerate(hist_list):
        hist.Draw(config[i]["DrawOption"] if i==0 else config[i]["DrawOption"] + " SAME")
        legend.AddEntry(hist, config[i]["legend"], config[i]["LegendDraw"])

    # draw legend
    legend.Draw()

    # set the ALEPH tag
    firstspace = 0.03
    space = 0.04
    ALEPHLabel(0.2, pltConfig["ALEPHTagTop"])
    myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 0*space, "#sqrt{s} = 91.2 GeV, 45 pb ^{-1}", size=0.035)
    if "pwflag" in hist_name:
        myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 1*space, pwflags[pwflag], size=0.035)
    elif "objSel" in hist_name:
        myText(0.2, pltConfig["ALEPHTagTop"] - firstspace - 1*space, f"Object Selection {objSel}", size=0.035)

    # Save to individual PDF
    pdf_name = os.path.join(outDir, f"{hist_name}.pdf")
    canvas.Print(pdf_name)
    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()

# Close all files
for f in files:
    f.Close()
