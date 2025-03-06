import ROOT
import os

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# custom
from PlotConfig import *
from style import *
SetALEPHStyle()

def ratio_hist(A, B):
    if not A or not B:
        return None  # Ensure histograms are valid
    
    # Ensure histograms have sum of weights squared enabled for error propagation
    A.Sumw2()
    B.Sumw2()

    result = A.Clone("result") # Clone A to store the result, so we don't modify the original histogram
    # result.Add(B, -1) # Subtract B from A (A - B)
    # result.Divide(A) # Divide by A: (A - B) / A
    result.Divide(B)
    return result  # Return the new histogram

# pick up the plot config
pltConfig = plotConfig["thrust"]

# will compare the following nominal (1) vs no selection (0), nominal (1) for closure, and all variations (2-7)
comparisons = [
    [1, 0, "Variation = No Event Selection"],
    [1, 1, "Variation = Nominal"],
    [1, 2, "Variation = Chg. Track N_{TPC} #geq 4 #rightarrow 7"],
    [1, 3, "Variation = Chg. Track p_{T} #geq 0.2 #rightarrow 0.4 GeV"],
    [1, 4, "Variation = E_{Ch} #geq 15 #rightarrow 10 GeV"],
    [1, 5, "Variation = Thrust w/o Neutral Objects"],
    [1, 6, "Variation = Thrust w/ #vec{p}_{miss} as Object"],
    [1, 7, "Variation = E_{Vis} #geq 0 #rightarrow 0.5E_{cm}"],
    [1, 8, "Variation = MissP < 20 GeV"],
    [1, 9, "Variation = |#vec{p}_{miss}| < 20 GeV and Thrust w/#vec{p}_{miss}"],
    [1, 10, "Variation = Neutral Object E #geq 0.4 #rightarrow 0.8 GeV"],
]

for iC, name in enumerate(["Data"]): # , "MC"
    
    # Open all files
    f = ROOT.TFile.Open(config[iC]["file"]) # pick up the data

    # output directory
    outDir = f"ThrustVariations{name}"
    os.makedirs(outDir, exist_ok=True)

    # get hist plot config
    for i,j,description in comparisons:
        print(i,j, description)
        # Create a new canvas for each histogram
        canvas = ALEPHCanvas(f"c_hist_sel{i}_{j}_thrust")

        # legend
        # legend = ALEPHLegend()

        # pick up histograms
        hist_i = f.Get(f"t_hist_sel{i}_thrust")
        hist_j = f.Get(f"t_hist_sel{j}_thrust")
        # sumw2
        hist_i.Sumw2()
        hist_j.Sumw2()
        # normalize first
        hist_i.Scale(1.0 / (hist_i.GetBinWidth(1) * hist_i.Integral()))
        hist_j.Scale(1.0 / (hist_i.GetBinWidth(1) * hist_j.Integral()))
        # get ratio
        hist = ratio_hist(hist_j, hist_i) # returns hist_A/hist_B
        hist.SetStats(0)  # Remove stats box
        hist.GetYaxis().SetTitleOffset(1.6)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetMaxDigits(3)
        # rebin and then set divisions
        hist.Rebin(pltConfig["rebin"])
        hist.GetXaxis().SetNdivisions(pltConfig["Ndivisions"])
        # set style
        hist.SetLineColor(config[iC]["color"])
        hist.SetMarkerColor(config[iC]["color"])
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.8)
        hist.SetLineWidth(2)
        hist.SetYTitle("Variation/Nominal") # (Nominal - Variation)/Nominal")
        hist.GetYaxis().SetRangeUser(0.02, 1.98) # -2.8,2.8)
        # hist.GetXaxis().SetRangeUser(0.5, 1.3)
        hist.Draw()

        # legend.AddEntry(hist, "Data 1994", "p")
        # legend.Draw()

        # Draw a dashed light gray horizontal line at y=1
        line = ROOT.TLine(hist.GetXaxis().GetXmin(), 1, hist.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kGray + 1)
        line.SetLineStyle(2)  # Dashed line
        line.SetLineWidth(3)
        line.Draw()

        # set the ALEPH tag
        firstspace = 0.045
        space = 0.045
        left = 0.9
        ALEPHLabel(left, pltConfig["ALEPHTagTop"], align=31)
        myText(left, pltConfig["ALEPHTagTop"] - firstspace - 0*space, "#sqrt{s} = 91.2 GeV, 45 pb ^{-1}", size=0.035, align=31)
        myText(left, pltConfig["ALEPHTagTop"] - firstspace - 1*space, description, size=0.035, align=31)

        # Save to individual PDF
        pdf_name = os.path.join(outDir, f"sel{i}_{j}_thrust.pdf")
        canvas.Print(pdf_name)
        print(f"Saved: {pdf_name}")

        # Cleanup
        canvas.Close()
