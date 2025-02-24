import ROOT
from ROOT import gROOT
import AtlasStyle as ats
gROOT.LoadMacro("AtlasUtils.C")

ROOT.SetAtlasStyle()

def ALEPHLabel(x, y, text=None, color=ROOT.kBlack, size=0.055):
    """
    Draw an 'ATLAS' label on a ROOT canvas at position (x, y), with an optional text.

    Parameters:
    - x, y: Position in NDC coordinates (0-1)
    - text: Optional string to appear next to "ATLAS"
    - color: ROOT color (default is black)
    """
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)  # Bold ATLAS font
    l.SetTextColor(color)
    l.SetTextSize(size)
    
    # Compute delx using the gPad dimensions
    if ROOT.gPad:
        delx = 0.115 * 696 * ROOT.gPad.GetWh() / (472 * ROOT.gPad.GetWw())
    else:
        delx = 0.115  # Default value if no active pad

    # Draw "ATLAS"
    l.DrawLatex(x, y, "Archived ALEPH")

    # Draw additional text if provided
    if text:
        p = ROOT.TLatex()
        p.SetNDC()
        p.SetTextFont(42)  # Regular text font
        p.SetTextColor(color)
        p.DrawLatex(x + delx, y, text)

def myText(x, y, text, color=ROOT.kBlack, align=12, size=0.03):
    """
    Draws a text label on a ROOT canvas at position (x, y) in NDC coordinates.

    Parameters:
    - x, y: Position in NDC coordinates (0-1)
    - color: ROOT color (e.g., ROOT.kBlack, ROOT.kRed)
    - text: The text string to be displayed
    """
    l = ROOT.TLatex()
    l.SetTextAlign(align)
    l.SetTextSize(size) 
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x, y, text)

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

# list of names to set log
SetLogy = []

# Open all files
files = [ROOT.TFile.Open(val["file"]) for val in config]

# Get the list of histograms in the first file (assuming all files have the same histograms)
hist_names = [key.GetName() for key in files[0].GetListOfKeys() if key.GetClassName().startswith("TH1")]
# only look for pwflag hists for now and then do the event variables
hist_names = [hist for hist in hist_names if "pwflag" in hist]
print(hist_names)

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# Loop over histogram names and plot each set of histograms
for hist_name in hist_names:
    hist_list = [f.Get(hist_name) for f in files]

    if not all(hist_list):  # Skip if any file is missing the histogram
        print(f"Warning: Missing histogram {hist_name} in some files.")
        continue

    # Create a new canvas for each histogram
    canvas = ROOT.TCanvas(f"c_{hist_name}", hist_name, 50, 50, 600, 600)
    # if hist_name in SetLogy:
    # canvas.SetLogy()

    # Set line colors
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange]
    
    # legend
    legend = ROOT.TLegend(0.68, 0.8, 0.88, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    
    # Normalize histograms if needed (optional)
    for i, hist in enumerate(hist_list):
        
        # general settings
        max_bin_content = 0
        hist.SetStats(0)  # Remove stats box
        hist.GetYaxis().SetTitleOffset(1.4)
        hist.GetXaxis().SetTitleOffset(1.4)
        hist.GetYaxis().SetMaxDigits(3)
        
        # h1->SetYTitle("d#sigma_{jet}/dE_{T,jet} [fb/GeV]");

        # create the y-axis title 
        XTitle = hist.GetXaxis().GetTitle()
        XBinWidth = hist.GetBinWidth(0)
        temp = XTitle.split("[")
        units = f" [pb/{XBinWidth}" + ("]" if len(temp) == 1 else " "+temp[1])
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
            hist.SetLineColor(colors[i % len(colors)])
            hist.SetLineWidth(2)
        else:
            print(f"Warning: {fname} does not match 'data' or 'MC'. Defaulting to black.")
            hist.SetLineColor(ROOT.kBlack)

    # Set y-axis range to 1.5 times the maximum bin content
    if max_bin_content > 0:
        hist_list[0].SetMaximum(1.5 * max_bin_content)
    
    # draw hists
    for i, hist in enumerate(hist_list):
        draw_option = "PE" if "LEP1" in config[i]["file"] else "HIST"
        hist.Draw(draw_option if i==0 else draw_option + " SAME")
        # Add to legend
        legend.AddEntry(hist, config[i]["legend"], "p" if "LEP1" in config[i]["file"] else "l")

    legend.Draw()

    # Add "Archived ALEPH" text in upper-left corner
    ALEPHLabel(0.2, 0.87)
    myText(0.2, 0.83, "#sqrt{s} = 91.2 GeV, X pb^{-1}", size=0.04)

    # Save to individual PDF
    pdf_name = f"{hist_name}.pdf"
    canvas.Print(pdf_name)

    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()

# Close all files
for f in files:
    f.Close()
