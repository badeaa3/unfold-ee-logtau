import ROOT
import os

# Run in batch mode to suppress UI popups
ROOT.gROOT.SetBatch(True)

# custom
from PlotConfig import *
from style import *
SetALEPHStyle()
from ThrustVariationCheck import get_ratio

if __name__ == "__main__":

    # pick up plot configurations
    hist_name = "logtau"
    pltConfig = plotConfig[hist_name]

    variations = {
        "Archived MC" : "/data/abadea/e+e-/aleph/unfold-ee-logtau/DataProcessing/20250527/2/alephMCRecoAfterCutPaths_1994_thrust_no_event_sel_tgenBefore.root",
        "Pythia 8" : "/data/abadea/e+e-/aleph/unfold-ee-logtau/DataProcessing/20250607/0/LEP1_PYTHIA8_MC_TGenBefore_NoISR_thrust_no_event_sel_tgenBefore.root",
        "Herwig" : "/data/abadea/e+e-/aleph/unfold-ee-logtau/DataProcessing/20250607/0/Herwig_noISR_ALL_thrust_no_event_sel_tgenBefore.root",
        "Sherpa" : "/data/abadea/e+e-/aleph/unfold-ee-logtau/DataProcessing/20250607/0/Sherpa_noISR_ALL_thrust_no_event_sel_tgenBefore.root",
    }

    # Create a new canvas for each histogram
    canvas = ALEPHCanvas(f"c_xyz")

    # Create a legend
    legend = ALEPHLegend(loc=[0.58, 0.2, 0.78, 0.4], textsize=0.04) 

    ratios = []
    # Add variations to the legend
    for variation, file in variations.items():
        ratios.append(get_ratio(file, variations["Archived MC"]))
        # Set color cycling through predefined colors
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
        ratios[-1].SetMarkerColor(colors[len(ratios)-1 % len(colors)])
        ratios[-1].SetLineColor(colors[len(ratios)-1 % len(colors)])
        ratios[-1].SetYTitle("Variation / Archived MC")
        ratios[-1].Draw("E1" if len(ratios) == 1 else "E1 same")
        legend.AddEntry(ratios[-1], variation, "lep")

    legend.Draw()

    # Draw a dashed light gray horizontal line at y=1
    line = ROOT.TLine(ratios[0].GetXaxis().GetXmin(), 1, ratios[0].GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kGray + 1)
    line.SetLineStyle(2)  # Dashed line
    line.SetLineWidth(3)
    line.Draw()

    header = "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Generator Level Theory Variation"
    myText(0.165, 0.97, header, size=0.03)
        
    # Save to individual PDF
    pdf_name = os.path.join(f"logtau_check_gen_theory_variation.pdf")
    canvas.Print(pdf_name)
    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()
