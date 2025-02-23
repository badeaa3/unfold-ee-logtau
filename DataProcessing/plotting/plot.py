import ROOT
from ROOT import gROOT
import AtlasStyle as ats
gROOT.LoadMacro("AtlasUtils.C")

ROOT.SetAtlasStyle()

def ALEPHLabel(x, y, text=None, color=ROOT.kBlack):
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

# Open all files
files = [ROOT.TFile.Open(val["file"]) for val in config]

# Get the list of histograms in the first file (assuming all files have the same histograms)
hist_names = [key.GetName() for key in files[0].GetListOfKeys() if key.GetClassName().startswith("TH1")][:2]

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# Loop over histogram names and plot each set of histograms
for hist_name in hist_names:
    hist_list = [f.Get(hist_name) for f in files]

    if not all(hist_list):  # Skip if any file is missing the histogram
        print(f"Warning: Missing histogram {hist_name} in some files.")
        continue

    # Create a new canvas for each histogram
    canvas = ROOT.TCanvas(f"c_{hist_name}", hist_name, 800, 600)

    # Set line colors
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange]
    
    # legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.045)
    
    # Normalize histograms if needed (optional)
    for i, hist in enumerate(hist_list):
        
        # general settings
        max_bin_content = 0
        hist.SetStats(0)  # Remove stats box

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
    ALEPHLabel(0.2, 0.875); 
    # ROOT.myText(0.34,0.8,1,"Preliminary");

    # Save to individual PDF
    pdf_name = f"{hist_name}.pdf"
    canvas.Print(pdf_name)

    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()

# Close all files
for f in files:
    f.Close()
