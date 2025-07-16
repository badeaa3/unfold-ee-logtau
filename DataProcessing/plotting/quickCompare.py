from ThrustVariationCheck import get_ratio
from style import *
import os
import ROOT

outDir = "quickCompare"
os.makedirs(outDir, exist_ok=True)

# List of file pairs to compare (reference_file, comparison_file)
file_pairs = [
    ["../20250716/0/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root", "../20250625/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root"],
    ["../20250716/1/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root", "../20250625/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root"],
    ["../20250716/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root", "../20250625/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root"],
    ["../20250716/3/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root", "../20250625/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root"],
    ["../20250716/4/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root", "../20250625/2/LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root"]
]

# List to store canvases to prevent garbage collection
canvases = []
# List to store histograms to prevent garbage collection
histograms = []

hist_names = ["pwflag0_pt", "pwflag1_pt", "pwflag2_pt", "logtau"]
for hist_name in hist_names:
    for i, (numer_file, denom_file) in enumerate(file_pairs):
        ratio_hist = get_ratio(numer_file, denom_file, hist_name=hist_name)
        
        if ratio_hist:
            # Create a new canvas for each comparison
            canvas = ROOT.TCanvas(f"c{i+1}_{hist_name}", f"Ratio Comparison {i+1}", 800, 600)
            
            # Set histogram style
            ratio_hist.SetLineColor(ROOT.kBlue)
            ratio_hist.SetLineWidth(2)
            ratio_hist.SetTitle(f"Ratio Comparison {i+1}")
            
            # Draw histogram
            ratio_hist.Draw("HIST")
            
            # Check for note.txt in the numerator file directory
            numer_dir = os.path.dirname(numer_file)
            numer_note_file = os.path.join(numer_dir, "note.txt")
            
            text_y_position = 0.9  # Starting Y position for text
            
            if os.path.exists(numer_note_file):
                try:
                    with open(numer_note_file, 'r') as f:
                        numer_note_text = f.read()
                    
                    # Create text object for numerator and draw it on the canvas
                    numer_lines = numer_note_text.split('\n')
                    # Draw header
                    myText(0.2, text_y_position, "Numerator notes:", color=ROOT.kBlack, align=12, size=0.025)
                    current_y = text_y_position - 0.03
                    for j, line in enumerate(numer_lines):
                        if line.strip():  # Only draw non-empty lines
                            myText(0.22, current_y - (j * 0.03), f"- {line}", color=ROOT.kBlack, align=12, size=0.025)
                    text_y_position = current_y - (len([l for l in numer_lines if l.strip()]) * 0.03) - 0.03  # Move down for next text
                    
                except Exception as e:
                    print(f"Warning: Could not read note.txt from {numer_note_file}: {e}")
            
            # Check for note.txt in the denominator file directory
            denom_dir = os.path.dirname(denom_file)
            denom_note_file = os.path.join(denom_dir, "note.txt")
            
            if os.path.exists(denom_note_file):
                try:
                    with open(denom_note_file, 'r') as f:
                        denom_note_text = f.read()
                        print(denom_note_text)

                    # Create text object for denominator and draw it on the canvas
                    denom_lines = denom_note_text.split('\n')
                    # Draw header
                    myText(0.2, text_y_position, "Denominator notes:", color=ROOT.kBlack, align=12, size=0.025)
                    current_y = text_y_position - 0.03
                    for j, line in enumerate(denom_lines):
                        if line.strip():  # Only draw non-empty lines
                            myText(0.22, current_y - (j * 0.03), f"- {line}", color=ROOT.kBlack, align=12, size=0.025)

                except Exception as e:
                    print(f"Warning: Could not read note.txt from {denom_note_file}: {e}")
            
            # Save canvas
            # canvas.SaveAs(f"{outDir}/ratio_comparison_{i+1}.png")
            canvas.SaveAs(f"{outDir}/ratio_comparison_{i+1}_{hist_name}.pdf")
            
            # Store canvas to prevent garbage collection
            canvases.append(canvas)
            # Store histogram to prevent garbage collection
            histograms.append(ratio_hist)
            
            print(f"Plot {i+1} saved to {outDir}/")

print("All plots saved!")