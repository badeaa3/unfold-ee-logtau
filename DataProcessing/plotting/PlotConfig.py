import ROOT

# config
config = []
config.append({
    "file" : "../LEP1Data1994_recons_aftercut-MERGED_thrust.root",
    "legend" : "Data 1994",
    "color" : ROOT.kBlack,
    "DrawOption" : "PE",
    "MarkerStyle" : 20,
    "LineWidth" : 2,
    "LegendDraw" : "p"
})
config.append({
    "file" : "../alephMCRecoAfterCutPaths_1994_thrust.root",
    "legend" : "Pythia 6.1",
    "DrawOption" : "HIST",
    "color" : ROOT.kRed,
    "MarkerStyle" : 20,
    "LineWidth" : 2,
    "LegendDraw" : "l"
})

# pwflag
pwflags = ["Charged Tracks", "Charged Leptons 1", "Charged Leptons 2", "V0", "Photons", "Neutral Hadrons"]

# ALEPH tag upper or lower
ALEPHTagUpperLeft = 0.87
ALEPHTagLowerLeft = 0.3

# individual plot settings
plotConfig = {}
# what we are plotting reduces to the ytitle in terms of sigma: 1/N dN/dX = 1/(sigma * L) d(sigma*L)/dX = 1/sigma d(sigma)/dX
plotConfig["cosTheta"] = {
    "SetLogy" : False, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dcos#theta", 
    "Ndivisions": 505, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 1.5
}
plotConfig["d0"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dd_{0}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 10
}
plotConfig["z0"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dz_{0}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 10
}
plotConfig["energy"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dE", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["mass"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dm", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 5
}
plotConfig["ntpc"] = {
    "SetLogy" : False, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dN_{TPC}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 1.5
}
plotConfig["phi"] = {
    "SetLogy" : False, "rebin" : 10, "YTitle": "1/#sigma d#sigma/d#phi", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 1.5
}
plotConfig["pmag"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/d|#vec{p}|", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["pt"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dp_{T}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["cosThetaSph"] = {
    "SetLogy" : False, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dcos#theta_{Sph}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 1.5
}
plotConfig["cosThetaThrust"] = {
    "SetLogy" : False, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dcos#theta_{Thr}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 1.5
}
plotConfig["eCh"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dE_{Ch}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["evis"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dE_{Vis}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["missP"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/d|#vec{p}_{miss}|", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["nneu"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dN_{Neu}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["ntrk"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dN_{Trk}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["ntrkPlusNeu"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dN_{Trk + Neu}", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagLowerLeft, "scale_max_bin_content" : 5
}
plotConfig["sphericity"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dSph", 
    "Ndivisions": 310, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 5
}
plotConfig["thrust"] = {
    "SetLogy" : True, "rebin" : 1, "YTitle": "1/#sigma d#sigma/dT", 
    "Ndivisions": 505, "ALEPHTagTop" : ALEPHTagUpperLeft, "scale_max_bin_content" : 5
}