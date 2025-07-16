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

    # zoom in and out versions
    zoomIn = False
    xrange = [-8, 0]
    if zoomIn:
        loc = [0.6, 0.2, 0.7, 0.3]
        yrange = [0.45, 1.9]
    else:
        loc = [0.6, 0.75, 0.75, 0.85]
        yrange = [0.1, 15]
    
    # Create a legend
    legend = ALEPHLegend(loc=loc, textsize=0.03)  # 

    ratios = []
    # Add variations to the legend
    for variation, file in variations.items():
        if variation == "Archived MC":
            continue
        ratios.append(get_ratio(file, variations["Archived MC"]))
        # Set color cycling through predefined colors
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
        ratios[-1].SetMarkerColor(colors[len(ratios)-1 % len(colors)])
        ratios[-1].SetMarkerSize(0.8)
        ratios[-1].SetMarkerStyle(20)
        ratios[-1].SetLineColor(colors[len(ratios)-1 % len(colors)])
        ratios[-1].SetYTitle("Variation / Archived MC")
        ratios[-1].GetYaxis().SetRangeUser(yrange[0], yrange[1])
        ratios[-1].GetXaxis().SetRangeUser(xrange[0], xrange[1])
        ratios[-1].Draw("E1" if len(ratios) == 1 else "E1 same")

        # Calculate and report agreement statistics
        y_vals = []
        for i in range(1, ratios[-1].GetNbinsX() + 1):
            y_vals.append(ratios[-1].GetBinContent(i))

        # Calculate mean and RMS deviation from 1
        mean_ratio = sum(y_vals) / len(y_vals)
        rms_deviation = (sum((y - 1)**2 for y in y_vals if y != 0) / len(y_vals))**0.5

        print(f"{variation} vs Archived MC:")
        print(f"  Mean ratio: {mean_ratio:.3f}")
        print(f"  RMS deviation from 1: {rms_deviation:.3f}")
        if rms_deviation < 0.05:
            print(f"  Agreement: Excellent (RMS < 5%)")
        elif rms_deviation < 0.1:
            print(f"  Agreement: Good (RMS < 10%)")
        else:
            print(f"  Agreement: Poor (RMS > 10%)")
        print()

        label = variation
        if variation != "Archived MC":
            label += f" (RMS {rms_deviation:.2f})"
        legend.AddEntry(ratios[-1], label, "lep")

    legend.Draw()

    # Draw a dashed light gray horizontal line at y=1
    line = ROOT.TLine(xrange[0], 1, xrange[1], 1)
    line.SetLineColor(ROOT.kGray + 1)
    line.SetLineStyle(2)  # Dashed line
    line.SetLineWidth(3)
    line.Draw()

    header = "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Generator Level Theory Variation"
    myText(0.165, 0.97, header, size=0.03)
        
    # Save to individual PDF
    extension = "ZoomIn" if zoomIn else "ZoomOut"
    pdf_name = os.path.join(f"logtau_check_gen_theory_variation_{extension}.pdf")
    canvas.Print(pdf_name)
    print(f"Saved: {pdf_name}")

    # Cleanup
    canvas.Close()
