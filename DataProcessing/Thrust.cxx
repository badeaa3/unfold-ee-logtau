/*
Author: Anthony Badea
Date: April 2, 2022
Analysis: MITHIG-MOD-20-001 Omnifold applied to ALEPH data
*/

// root dependencies
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TH1D.h"

// thrust code
#include "thrustTools.h"
#include "sphericityTools.h"

// c++ code
#include <vector>
#include <iostream>
#include <iomanip>


// track selection variation
std::map<std::string, float> getTrackVariation(
  /* charged track selections */
  int applySelection, // flag indicating if the selection should be applied to thrust and sphericity calculation
  int nTPCcut, // 2PC paper value: 4
  float chargedTracksAbsCosThCut, // 0.94
  float ptCut, // 0.2
  float d0Cut, // 2
  float z0Cut, // 10
  /* neutral track selections */
  float ECut, // 0.4
  float neutralTracksAbsCosThCut, // 0.98
  bool keepChargedTracks, // include charged tracks in calculation
  bool keepNeutralTracks,  // include neutral tracks in calculation
  bool doMET // include missing momentum vector in calculation
) {
  return std::map<std::string, float> {
    {"applyTrackSelection", applySelection},
    {"nTPCcut", nTPCcut},
    {"chargedTracksAbsCosThCut", chargedTracksAbsCosThCut},
    {"ptCut", ptCut},
    {"d0Cut", d0Cut},
    {"z0Cut", z0Cut},
    {"ECut", ECut},
    {"neutralTracksAbsCosThCut", neutralTracksAbsCosThCut},
    {"keepChargedTracks", keepChargedTracks},
    {"keepNeutralTracks", keepNeutralTracks},
    {"doMET", doMET},
  };
}

// event selection variation
std::map<std::string, float> getEventVariation(
  /* event selections */
  float TotalTrkEnergyCut, // 2PC paper value: 15
  float AbsCosSThetaCut, // 0.82
  float NTrkCut, // 5
  float NeuNchCut // 13
) {
  return std::map<std::string, float> {
    {"TotalTrkEnergyCut", TotalTrkEnergyCut},
    {"AbsCosSThetaCut", AbsCosSThetaCut},
    {"NTrkCut", NTrkCut},
    {"NeuNchCut", NeuNchCut}
  };
}

// progressbar used during the event loop
void pbftp(double time_diff, int nprocessed, int ntotal) {
  double rate = (double)(nprocessed + 1) / time_diff;
  std::cout << "\r > " << nprocessed << " / " << ntotal
            << " | "   << std::fixed << std::setprecision(1) << 100 * (double)(nprocessed) / (double)(ntotal) << "%"
            << " | "   << std::fixed << std::setprecision(1) << rate << "Hz"
            << " | "   << std::fixed << std::setprecision(1) << time_diff / 60 << "m elapsed"
            << " | "   << std::fixed << std::setprecision(1) << (double)(ntotal - nprocessed) / (rate * 60) << "m remaining";
  std::cout << std::flush;
}

/**
 * Main event loop to apply systematic variation of track and event selections to compute thrust.
 *
 * @param user supplies the input file.
 * Optionally the output file name (default to input file name with .root replaced with _thrust.root) and number of events to run over (default to all).
 * @return int and saves new root file with trees inside of it.
 */
int main(int argc, char* argv[]) {

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::string inFileName = "";
  std::string outFileName = "";
  int doNEvents = -1;
  bool debug = false;
  size_t divide    = 1; // how many divisions the file is being cut into, which division you are on, and how many events per division
  size_t thisdiv   = 0;
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "-i", 2) == 0) inFileName = argv[i + 1];
    if (strncmp(argv[i], "-o", 2) == 0) outFileName = argv[i + 1];
    // if (strncmp(argv[i], "-t", 2) == 0) treeName = argv[i + 1];
    if (strncmp(argv[i], "-n", 2) == 0) doNEvents = std::stoi(argv[i + 1]);
    if (strncmp(argv[i], "--debug", 7) == 0) debug = true;
    if (strncmp(argv[i], "--divide", 8) == 0) divide = (size_t)(std::atoi(argv[i + 1]));
    if (strncmp(argv[i], "--thisdiv", 9) == 0) thisdiv = (size_t)(std::atoi(argv[i + 1]));
  }
  if (inFileName == "") {
    std::cout << "No input file name provided. Exiting" << std::endl;
    return 0;
  }
  if (outFileName == "" && inFileName != "") {
    outFileName = inFileName;
    outFileName.erase(outFileName.length() - 5); // strip .root
    outFileName += "_thrust.root";
  }
  // set tree names
  std::vector<std::string> treeNames{"t"};
  if (inFileName.find("LEP1Data") == std::string::npos) {
    treeNames.push_back("tgen");
    treeNames.push_back("tgenBefore");
  }

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::unique_ptr<TFile> f (new TFile(inFileName.c_str(), "READ"));

  // declare variables
  int maxPart = 500;
  unsigned long long uniqueID;
  int nParticle;
  bool passesNTupleAfterCut;
  Short_t pwflag[maxPart];
  float theta[maxPart];
  float phi[maxPart];
  float pt[maxPart];
  float d0[maxPart];
  float z0[maxPart];
  Short_t ntpc[maxPart];
  float px[maxPart];
  float py[maxPart];
  float pz[maxPart];
  float pmag[maxPart];
  float mass[maxPart];

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Output File %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::unique_ptr<TFile> fout (new TFile(outFileName.c_str(), "RECREATE"));

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Track Variations %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::vector<std::map<std::string, float> > trackVariations; // vector of variations
  // trackVariations.push_back(getTrackVariation(0, 4, 0.94, 0.2, 2, 10, 0.4, 0.98, true, true, false)); // no track selection
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.98, true, true, false)); // 2PC track selection
  trackVariations.push_back(getTrackVariation(1, 7, 0.94, 0.2, 2, 10, 0.4, 0.98, true, true, false)); // ntpc varied from 4 -> 7
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.4, 2, 10, 0.4, 0.98, true, true, false)); // charged tracks pT cut from 0.2 to 0.4 GeV
  // trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.2, 0.98, true, true, false)); // neutral ECut variation scaling down by 1/2
  // trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.8, 0.98, true, true, false)); // neutral ECut variation scaling up by 2
  // trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.96, true, true, false)); // neutralTracksAbsCosThCut subtracting 0.02
  // trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 1.00, true, true, false)); // neutralTracksAbsCosThCut adding 0.02
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.98, true, false, false)); // charged tracks only thrust
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.98, false, true, false)); // neutral tracks only thrust
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.98, true, true, true)); // include met in thrust

  // vectors for selected objects
  std::vector<int> selectedParts;
  std::vector<std::vector<float> > selectedPx, selectedPy, selectedPz;
  std::vector<std::vector<Short_t> > selectedPwflag;
  // event level quantities
  TVector3 thrust;
  std::unique_ptr<Sphericity> spher;

  // save variation definitions to a tree
  std::unique_ptr<TTree> varDefs (new TTree("TrackVariationDefinitions", ""));
  int nTPCcut;
  float applyTrackSelection, chargedTracksAbsCosThCut, ptCut, d0Cut, z0Cut, ECut, neutralTracksAbsCosThCut;
  bool keepChargedTracks, keepNeutralTracks;
  varDefs->Branch("applyTrackSelection", &applyTrackSelection);
  varDefs->Branch("nTPCcut", &nTPCcut);
  varDefs->Branch("chargedTracksAbsCosThCut", &chargedTracksAbsCosThCut);
  varDefs->Branch("ptCut", &ptCut);
  varDefs->Branch("d0Cut", &d0Cut);
  varDefs->Branch("z0Cut", &z0Cut);
  varDefs->Branch("ECut", &ECut);
  varDefs->Branch("neutralTracksAbsCosThCut", &neutralTracksAbsCosThCut);
  varDefs->Branch("keepChargedTracks", &keepChargedTracks);
  varDefs->Branch("keepNeutralTracks", &keepNeutralTracks);

  // push back for each variation and save to tree
  for (unsigned int iV = 0; iV < trackVariations.size(); iV++) {
    selectedParts.push_back(0);
    selectedPx.push_back(std::vector<float>());
    selectedPy.push_back(std::vector<float>());
    selectedPz.push_back(std::vector<float>());
    selectedPwflag.push_back(std::vector<Short_t>());
    applyTrackSelection = trackVariations.at(iV)["applyTrackSelection"];
    nTPCcut = trackVariations.at(iV)["nTPCcut"];
    chargedTracksAbsCosThCut = trackVariations.at(iV)["chargedTracksAbsCosThCut"];
    ptCut = trackVariations.at(iV)["ptCut"];
    d0Cut = trackVariations.at(iV)["d0Cut"];
    z0Cut = trackVariations.at(iV)["z0Cut"];
    ECut = trackVariations.at(iV)["ECut"];
    neutralTracksAbsCosThCut = trackVariations.at(iV)["neutralTracksAbsCosThCut"];
    keepChargedTracks = trackVariations.at(iV)["keepChargedTracks"];
    keepNeutralTracks = trackVariations.at(iV)["keepNeutralTracks"];
    varDefs->Fill();
  }
  varDefs->Write();

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Event Variations %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::vector<std::map<std::string, float> > eventVariations; // vector of variations
  eventVariations.push_back(getEventVariation(15, 0.82, 5, 13)); // nominal selections
  eventVariations.push_back(getEventVariation(10, 0.82, 5, 13)); // total charged energy varied from 15 -> 10 GeV

  // save variation definitions to a tree
  std::unique_ptr<TTree> evtVarDefs (new TTree("EventVariationDefinitions", ""));
  float TotalTrkEnergyCut, AbsCosSThetaCut;
  int NTrkCut, NeuNchCut;
  evtVarDefs->Branch("TotalTrkEnergyCut", &TotalTrkEnergyCut);
  evtVarDefs->Branch("AbsCosSThetaCut", &AbsCosSThetaCut);
  evtVarDefs->Branch("NTrkCut", &NTrkCut);
  evtVarDefs->Branch("NeuNchCut", &NeuNchCut);
  for (unsigned int iEV = 0; iEV < eventVariations.size(); iEV++) {
    TotalTrkEnergyCut = eventVariations.at(iEV)["TotalTrkEnergyCut"];
    AbsCosSThetaCut = eventVariations.at(iEV)["AbsCosSThetaCut"];
    NTrkCut = eventVariations.at(iEV)["NTrkCut"];
    NeuNchCut = eventVariations.at(iEV)["NeuNchCut"];
    evtVarDefs->Fill();
  }
  evtVarDefs->Write();

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Data Tree %%%%%%%%%%%%%%%%%%%%%%%%%%#
  // std::vector<float> Thrust, TotalTrkEnergy, STheta;
  // std::vector<int> NTrk, Neu;
  // std::vector<std::vector<bool> > passEventSelection(eventVariations.size()); // allocate memory here without push_back to avoid copying which confuses tree->Branch  
    
  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Event Loop %%%%%%%%%%%%%%%%%%%%%%%%%%#

  // timekeeper
  std::chrono::time_point<std::chrono::system_clock> time_start;
  std::chrono::duration<double> elapsed_seconds;

  // loop over trees
  for (auto &tree : treeNames) {

    std::cout << TString::Format("Looping over tree: %s", tree.c_str()) << std::endl;

    // load input tree
    bool genTree = tree == "tgen" || tree == "tgenBefore";
    std::unique_ptr<TTree> t ((TTree*) f->Get(tree.c_str()));
    t->SetBranchAddress("uniqueID", &uniqueID);
    t->SetBranchAddress("nParticle", &nParticle);
    t->SetBranchAddress("passesNTupleAfterCut", &passesNTupleAfterCut);
    t->SetBranchAddress("pwflag", &pwflag);
    t->SetBranchAddress("theta", &theta);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("pt", &pt);
    t->SetBranchAddress("d0", &d0);
    t->SetBranchAddress("z0", &z0);
    t->SetBranchAddress("ntpc", &ntpc);
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("pmag", &pmag);
    t->SetBranchAddress("mass", &mass);

    // create output tree
    std::unique_ptr<TTree> tout (new TTree(tree.c_str(), ""));
    unsigned long long uniqueIDCopy; 
    std::vector<float> Thrust, TotalTrkEnergy, STheta, Sph, MissP, EVis, TTheta;
    std::vector<int> NTrk, Neu;
    std::vector<std::vector<bool> > passEventSelection(eventVariations.size()); // allocate memory here without push_back to avoid copying which confuses tree->Branch
    tout->Branch("uniqueID", &uniqueIDCopy);
    tout->Branch("Thrust", &Thrust);
    tout->Branch("TotalTrkEnergy", &TotalTrkEnergy);
    tout->Branch("NTrk", &NTrk);
    tout->Branch("Neu", &Neu);
    tout->Branch("STheta", &STheta);
    tout->Branch("Sphericity", &Sph);
    tout->Branch("MissP", &MissP);
    tout->Branch("EVis", &EVis);
    tout->Branch("TTheta", &TTheta);
    for (unsigned int iEV = 0; iEV < eventVariations.size(); iEV++) {
      tout->Branch(("passEventSelection_" + std::to_string(iEV)).c_str(), &passEventSelection.at(iEV));
    }

    // create object level histograms for each pwflag
    std::map<std::pair<int, std::string>, TH1D*> particleHists;
    for(int iP=0; iP <= 5; iP++){
      particleHists[{iP, "cosTheta"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "cosTheta").c_str(), ";cos#theta;Entries", 200, -1, 1);
      particleHists[{iP, "phi"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "phi").c_str(), ";#phi;Entries", 800, -4, 4);
      particleHists[{iP, "pt"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "pt").c_str(), ";p_{T} [GeV];Entries", 100, 0, 100);
      particleHists[{iP, "ntpc"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "ntpc").c_str(), ";N_{TPC};Entries", 31, -0.5, 30.5);
      particleHists[{iP, "d0"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "d0").c_str(), ";d_{0} [cm];Entries", 50, -2.5, 2.5);
      particleHists[{iP, "z0"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "z0").c_str(), ";z_{0} [cm];Entries", 300, -15, 15);
      particleHists[{iP, "pmag"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "pmag").c_str(), ";|#vec{p}| [GeV];Entries", 100, 0, 100);
      particleHists[{iP, "mass"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "mass").c_str(), ";Mass [GeV];Entries", 100, 0, 10);
      particleHists[{iP, "energy"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "energy").c_str(), ";Energy [GeV];Entries", 100, 0, 100);
    }

    // create event level histograms for each track selection variation
    std::map<std::pair<int, std::string>, TH1D*> eventHists;
    for (unsigned int iV = 0; iV < trackVariations.size(); iV++){
      if (genTree && iV > 0 ) break;
      eventHists[{iV, "ntrk"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "ntrk").c_str(), ";N_{Trk};Entries", 61, -0.5, 60.5);
      eventHists[{iV, "nneu"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "nneu").c_str(), ";N_{Neu};Entries", 51, -0.5, 50.5);
      eventHists[{iV, "ntrkPlusNeu"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "ntrkPlusNeu").c_str(), ";N_{Trk+Neu};Entries", 81, -0.5, 80.5);
      eventHists[{iV, "eCh"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "eCh").c_str(), ";E_{Ch} [GeV];Entries", 200, 0, 200);
      eventHists[{iV, "cosThetaSph"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "cosThetaSph").c_str(), ";cos#theta_{Sph};Entries", 100, -1, 1);
      eventHists[{iV, "sphericity"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "sphericity").c_str(), ";Sphericity;Entries", 100, 0, 1);
      eventHists[{iV, "thrust"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "thrust").c_str(), ";Thrust;Entries", 100, 0.5, 1);
      eventHists[{iV, "missP"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "missP").c_str(), ";|#vec{p}_{MET}| [GeV];Entries", 100, 0, 100);
      eventHists[{iV, "evis"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "evis").c_str(), ";E_{Vis} [GeV];Entries", 200, 0, 200);
      eventHists[{iV, "cosThetaThrust"}] = new TH1D( (tree + "_hist_objSel" + std::to_string(iV) + "_" + "cosThetaThrust").c_str(), ";cos#theta_{Thr};Entries", 100, -1, 1);
    }

    
    // interpret divide and thisdiv to event range
    int nEvents = t->GetEntries();
    int evtperdiv = nEvents / divide;
    int startevt  = evtperdiv * thisdiv;
    int endevt    = (divide == (thisdiv + 1)) ? nEvents : evtperdiv * (thisdiv + 1); // if the last division go till the end
    int ntotal    = endevt - startevt;
    // user says do nEvents
    if (doNEvents != -1) {
      startevt = 0;
      endevt = doNEvents;
      ntotal = doNEvents;
    }
    std::cout << TString::Format("Total events %d, Events per division %d, Start event %d, End event %d, Analysing %d events", nEvents, evtperdiv, startevt, endevt, ntotal) << std::endl;
    // start clock
    time_start = std::chrono::system_clock::now();
    for (int iE = startevt; iE < endevt; iE++) {

      // progressbar
      if (!debug) {
        elapsed_seconds = (std::chrono::system_clock::now() - time_start);
        // pbftp(elapsed_seconds.count(), iE + 1, nEvents);
        pbftp(elapsed_seconds.count(), iE + 1 - startevt, ntotal);
      }

      t->GetEntry(iE);

      // set uniqueID
      uniqueIDCopy = uniqueID;

      // reset variables
      TotalTrkEnergy.clear();
      NTrk.clear();
      Neu.clear();
      STheta.clear();
      Sph.clear();
      MissP.clear();
      EVis.clear();
      Thrust.clear();
      TTheta.clear();
      for (unsigned int iV = 0; iV < trackVariations.size(); iV++) {
        if (genTree && iV > 0 ) break;
        selectedParts.at(iV) = 0;
        selectedPx.at(iV).clear();
        selectedPy.at(iV).clear();
        selectedPz.at(iV).clear();
        selectedPwflag.at(iV).clear();
        TotalTrkEnergy.push_back(0);
	EVis.push_back(0);
        NTrk.push_back(0);
        Neu.push_back(0);
      }
      for (auto &ev : passEventSelection) ev.clear();

      // compute event selection variables
      for (int iP = 0; iP < nParticle; iP++) {

        if (debug) std::cout << TString::Format("iP %d, pwflag %d, theta %f, pt %f, d0 %f, z0 %f, ntpc %d", iP, pwflag[iP], theta[iP], pt[iP], d0[iP], z0[iP], ntpc[iP]) << std::endl;

	// compute the particle energy
	float energy = TMath::Sqrt(pmag[iP] * pmag[iP] + mass[iP] * mass[iP]);
	
	// fill histogram for all particles
	if(pwflag[iP] >= 0 && pwflag[iP] <= 5){
	  particleHists[{pwflag[iP], "cosTheta"}]->Fill(cos(theta[iP]));
	  particleHists[{pwflag[iP], "phi"}]->Fill(phi[iP]);
	  particleHists[{pwflag[iP], "pt"}]->Fill(pt[iP]);
	  particleHists[{pwflag[iP], "ntpc"}]->Fill(ntpc[iP]);
	  particleHists[{pwflag[iP], "d0"}]->Fill(d0[iP]);
	  particleHists[{pwflag[iP], "z0"}]->Fill(z0[iP]);
	  particleHists[{pwflag[iP], "pmag"}]->Fill(pmag[iP]);
	  particleHists[{pwflag[iP], "mass"}]->Fill(mass[iP]);
	  particleHists[{pwflag[iP], "energy"}]->Fill(energy);
	}
	  
        // loop over variations
        for (unsigned int iV = 0; iV < trackVariations.size(); iV++) {

          // if tgen or tgenbefore only do the first nominal variation
          if (genTree && iV > 0 ) break;
	  	  
          // count charged tracks
          bool chargedTrackSelections =
            (pwflag[iP] >= 0 && pwflag[iP] <= 2)
            && TMath::Abs(cos(theta[iP])) <= trackVariations.at(iV)["chargedTracksAbsCosThCut"]
            && pt[iP] >= trackVariations.at(iV)["ptCut"]
            && TMath::Abs(d0[iP]) <= trackVariations.at(iV)["d0Cut"]
            && TMath::Abs(z0[iP]) <= trackVariations.at(iV)["z0Cut"]
            && ntpc[iP] >= trackVariations.at(iV)["nTPCcut"];
          // count charged tracks
          if (chargedTrackSelections) {
            TotalTrkEnergy.at(iV) += energy;
	    EVis.at(iV) += energy;
            NTrk.at(iV) += 1;
            if (debug) std::cout << "Passed charged track selection" << std::endl;
          }

          // count neutral tracks
          bool neutralTrackSelections =
            (pwflag[iP] == 4 || pwflag[iP] == 5)
	    && energy >= trackVariations.at(iV)["ECut"]
            && TMath::Abs(cos(theta[iP])) <= trackVariations.at(iV)["neutralTracksAbsCosThCut"];
          if (neutralTrackSelections) {
	    EVis.at(iV) += energy;
            Neu.at(iV) += 1;
            if (debug) std::cout << "Passed neutral track selection" << std::endl;
          }

          // add to input list for thrust. check for -1 which indicates use all tracks for thrust
          bool keeptrack = tree == "tgen" || tree == "tgenBefore"; // generator level always keep
          keeptrack = keeptrack || trackVariations.at(iV)["applyTrackSelection"] == 0 ; // track selection should not be applied
	  // make choice if charged and neutral tracks are kept in thrust calculation and then verify their selection
	  keeptrack = keeptrack || (trackVariations.at(iV)["keepChargedTracks"] && chargedTrackSelections);
	  keeptrack = keeptrack	|| (trackVariations.at(iV)["keepNeutralTracks"] && neutralTrackSelections);
	  // save the track 
          if (keeptrack) {
            selectedParts.at(iV) += 1;
            selectedPx.at(iV).push_back(px[iP]);
            selectedPy.at(iV).push_back(py[iP]);
            selectedPz.at(iV).push_back(pz[iP]);
            selectedPwflag.at(iV).push_back(pwflag[iP]);
          }
          else {
            if (debug) std::cout << "Did not pass either charged or neutral track selection" << std::endl;
          }
        }
      }
      if (debug) std::cout << TString::Format("Nominal N Passing Charged Tracks %d, N Passing Neutral Tracks %d, and selected tracks %d", NTrk.at(0), Neu.at(0), selectedParts.at(0)) << std::endl;

      // compute event level variables
      for (unsigned int iV = 0; iV < trackVariations.size(); iV++) {

        // if tgen or tgenbefore only do the first nominal variation
        if (genTree && iV > 0 ) break;

	// calculate the missing momentum vector
	TVector3 met = TVector3(0, 0, 0);
	for (int t = 0; t < selectedParts.at(iV); t++) {
	  met += (TVector3(selectedPx.at(iV).at(t), selectedPy.at(iV).at(t), selectedPz.at(iV).at(t)));
	}
	met = -met;
	MissP.push_back(met.Mag());
  
        // sphericity
        spher = std::make_unique<Sphericity>(Sphericity(selectedParts.at(iV), selectedPx.at(iV).data(), selectedPy.at(iV).data(), selectedPz.at(iV).data(), selectedPwflag.at(iV).data(), false));
        STheta.push_back(spher->sphericityAxis().Theta());
	Sph.push_back(spher->sphericity());
      
        // thrust
	// TVector3 getThrust(int n, float *px, float *py, float *pz, THRUST::algorithm algo = THRUST::HERWIG, bool doWeight = false, bool doInvertWeight = false, float* weight = NULL, bool doMET = false, Short_t *pwflag = NULL)
        thrust = getThrust(
			   selectedParts.at(iV), // n
			   selectedPx.at(iV).data(), // px
			   selectedPy.at(iV).data(), // py
			   selectedPz.at(iV).data(), // pz 
			   THRUST::OPTIMAL, // THRUST::algorithm
			   false, // doWeight
			   false, // doInvertWeight
			   NULL, // weight
			   trackVariations.at(iV)["doMET"],  // doMET
			   NULL // pwflag
			   );
        Thrust.push_back(thrust.Mag());
	TTheta.push_back(thrust.Theta());
	
        // compute event selection passes
        for (unsigned int iEV = 0; iEV < eventVariations.size(); iEV++) {
          passEventSelection.at(iEV).push_back(
            passesNTupleAfterCut == 1
            && TotalTrkEnergy.at(iV) >= eventVariations.at(iEV)["TotalTrkEnergyCut"]
            && TMath::Abs(TMath::Cos(STheta.at(iV))) <= eventVariations.at(iEV)["AbsCosSThetaCut"]
            && NTrk.at(iV) >= eventVariations.at(iEV)["NTrkCut"]
            && (NTrk.at(iV) + Neu.at(iV)) >= eventVariations.at(iEV)["NeuNchCut"]
          );
        }

	// fill histograms
	eventHists[{iV, "ntrk"}]->Fill(NTrk.at(iV));
	eventHists[{iV, "nneu"}]->Fill(Neu.at(iV));
	eventHists[{iV, "ntrkPlusNeu"}]->Fill(NTrk.at(iV) + Neu.at(iV));
	eventHists[{iV, "eCh"}]->Fill(TotalTrkEnergy.at(iV));
	eventHists[{iV, "cosThetaSph"}]->Fill(TMath::Cos(STheta.at(iV)));
	eventHists[{iV, "sphericity"}]->Fill(Sph.at(iV));
	eventHists[{iV, "thrust"}]->Fill(Thrust.at(iV));
	eventHists[{iV, "missP"}]->Fill(MissP.at(iV));
	eventHists[{iV, "evis"}]->Fill(EVis.at(iV));
	eventHists[{iV, "cosThetaThrust"}]->Fill(TMath::Cos(TTheta.at(iV)));
	
      }
      if (debug) std::cout << TString::Format("Nominal STheta %f, Thrust %f", STheta.at(0), Thrust.at(0)) << std::endl;

      tout->Fill();
    }

    // write tree
    tout->Write();

    // write histograms
    for (auto &entry : particleHists) {
        entry.second->Write();
    }
    for (auto &entry : eventHists) {
        entry.second->Write();
    }
    
    std::cout << "\n" << std::endl;
  }

  return 1;
}
