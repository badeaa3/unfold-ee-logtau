import ROOT

def SetALEPHStyle():
    """Applies ALEPH style settings in PyROOT."""
    print("\nApplying ALEPH style settings...\n")
    aleph_style = aleph_style_settings()
    ROOT.gROOT.SetStyle("ALEPH")
    ROOT.gROOT.ForceStyle()

def aleph_style_settings():
    """Defines and returns the ALEPH style settings."""
    aleph_style = ROOT.TStyle("ALEPH", "ALEPH style")

    # Use plain black on white colors
    icol = 0  # WHITE
    aleph_style.SetFrameBorderMode(icol)
    aleph_style.SetFrameFillColor(icol)
    aleph_style.SetCanvasBorderMode(icol)
    aleph_style.SetCanvasColor(icol)
    aleph_style.SetPadBorderMode(icol)
    aleph_style.SetPadColor(icol)
    aleph_style.SetStatColor(icol)

    # Set the paper & margin sizes
    aleph_style.SetPaperSize(20, 26)

    # Set margin sizes
    aleph_style.SetPadTopMargin(0.05)
    aleph_style.SetPadRightMargin(0.05)
    aleph_style.SetPadBottomMargin(0.16)
    aleph_style.SetPadLeftMargin(0.16)

    # Set title offsets (for axis labels)
    aleph_style.SetTitleXOffset(1.4)
    aleph_style.SetTitleYOffset(1.4)

    # Use large fonts
    font = 42  # Helvetica
    tsize = 0.05
    aleph_style.SetTextFont(font)
    aleph_style.SetTextSize(tsize)
    aleph_style.SetLabelFont(font, "x")
    aleph_style.SetTitleFont(font, "x")
    aleph_style.SetLabelFont(font, "y")
    aleph_style.SetTitleFont(font, "y")
    aleph_style.SetLabelFont(font, "z")
    aleph_style.SetTitleFont(font, "z")

    aleph_style.SetLabelSize(tsize, "x")
    aleph_style.SetTitleSize(tsize, "x")
    aleph_style.SetLabelSize(tsize, "y")
    aleph_style.SetTitleSize(tsize, "y")
    aleph_style.SetLabelSize(tsize, "z")
    aleph_style.SetTitleSize(tsize, "z")

    # Use bold lines and markers
    aleph_style.SetMarkerStyle(20)
    aleph_style.SetMarkerSize(1.2)
    aleph_style.SetHistLineWidth(2)
    aleph_style.SetLineStyleString(2, "[12 12]")  # Postscript dashes

    # Get rid of X error bars and error bar caps
    aleph_style.SetEndErrorSize(0.)

    # Do not display standard histogram decorations
    aleph_style.SetOptTitle(0)
    aleph_style.SetOptStat(0)
    aleph_style.SetOptFit(0)

    # Put tick marks on top and RHS of plots
    aleph_style.SetPadTickX(1)
    aleph_style.SetPadTickY(1)

    return aleph_style

def ALEPHCanvas(name):
    return ROOT.TCanvas(f"c_{name}", name, 50, 50, 600, 600)

def ALEPHLegend():
    legend = ROOT.TLegend(0.68, 0.8, 0.88, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    return legend 


def ALEPHLabel(x, y, text=None, color=ROOT.kBlack, size=0.055):
    """
    Draw an 'ALEPH' label on a ROOT canvas at position (x, y), with an optional text.

    Parameters:
    - x, y: Position in NDC coordinates (0-1)
    - text: Optional string to appear next to "ALEPH"
    - color: ROOT color (default is black)
    """
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)  # Bold ALEPH font
    l.SetTextColor(color)
    l.SetTextSize(size)
    
    # Compute delx using the gPad dimensions
    if ROOT.gPad:
        delx = 0.115 * 696 * ROOT.gPad.GetWh() / (472 * ROOT.gPad.GetWw())
    else:
        delx = 0.115  # Default value if no active pad

    # Draw "ALEPH"
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