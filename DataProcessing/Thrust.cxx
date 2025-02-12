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
  float neutralTracksAbsCosThCut // 0.98
) {
  return std::map<std::string, float> {
    {"applyTrackSelection", applySelection},
    {"nTPCcut", nTPCcut},
    {"chargedTracksAbsCosThCut", chargedTracksAbsCosThCut},
    {"ptCut", ptCut},
    {"d0Cut", d0Cut},
    {"z0Cut", z0Cut},
    {"ECut", ECut},
    {"neutralTracksAbsCosThCut", neutralTracksAbsCosThCut}
  };
}

// event selection variation
std::map<std::string, float> getEventVariation(
  /* event selections */
  float TotalChgEnergyCut, // 2PC paper value: 15
  float AbsCosSThetaCut, // 0.82
  float NTrkCut, // 5
  float NeuNchCut // 13
) {
  return std::map<std::string, float> {
    {"TotalChgEnergyCut", TotalChgEnergyCut},
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
  trackVariations.push_back(getTrackVariation(0, 4, 0.94, 0.2, 2, 10, 0.4, 0.98)); // no track selection
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.98)); // 2PC track selection
  for (int i = 5; i <= 7; i++){
    trackVariations.push_back(getTrackVariation(1, i, 0.94, 0.2, 2, 10, 0.4, 0.98)); // 2PC ntpc variations
  }
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.2, 0.98)); // neutral ECut variation scaling down by 1/2
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.8, 0.98)); // neutral ECut variation scaling up by 2
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 0.96)); // neutralTracksAbsCosThCut subtracting 0.02
  trackVariations.push_back(getTrackVariation(1, 4, 0.94, 0.2, 2, 10, 0.4, 1.00)); // neutralTracksAbsCosThCut adding 0.02

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
  varDefs->Branch("applyTrackSelection", &applyTrackSelection);
  varDefs->Branch("nTPCcut", &nTPCcut);
  varDefs->Branch("chargedTracksAbsCosThCut", &chargedTracksAbsCosThCut);
  varDefs->Branch("ptCut", &ptCut);
  varDefs->Branch("d0Cut", &d0Cut);
  varDefs->Branch("z0Cut", &z0Cut);
  varDefs->Branch("ECut", &ECut);
  varDefs->Branch("neutralTracksAbsCosThCut", &neutralTracksAbsCosThCut);

  // push back for each variation and save to tree
  for (int iV = 0; iV < trackVariations.size(); iV++) {
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
    varDefs->Fill();
  }
  varDefs->Write();

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Event Variations %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::vector<std::map<std::string, float> > eventVariations; // vector of variations
  // nominal values
  eventVariations.push_back(getEventVariation(15, 0.82, 5, 13));
  // total charged energy variation
  eventVariations.push_back(getEventVariation(10, 0.82, 5, 13));

  // save variation definitions to a tree
  std::unique_ptr<TTree> evtVarDefs (new TTree("EventVariationDefinitions", ""));
  float TotalChgEnergyCut, AbsCosSThetaCut;
  int NTrkCut, NeuNchCut;
  evtVarDefs->Branch("TotalChgEnergyCut", &TotalChgEnergyCut);
  evtVarDefs->Branch("AbsCosSThetaCut", &AbsCosSThetaCut);
  evtVarDefs->Branch("NTrkCut", &NTrkCut);
  evtVarDefs->Branch("NeuNchCut", &NeuNchCut);
  for (int iV = 0; iV < eventVariations.size(); iV++) {
    TotalChgEnergyCut = eventVariations.at(iV)["TotalChgEnergyCut"];
    AbsCosSThetaCut = eventVariations.at(iV)["AbsCosSThetaCut"];
    NTrkCut = eventVariations.at(iV)["NTrkCut"];
    NeuNchCut = eventVariations.at(iV)["NeuNchCut"];
    evtVarDefs->Fill();
  }
  evtVarDefs->Write();

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Data Tree %%%%%%%%%%%%%%%%%%%%%%%%%%#
  // std::vector<float> Thrust, TotalChgEnergy, STheta;
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
    std::vector<float> Thrust, TotalChgEnergy, STheta;
    std::vector<int> NTrk, Neu;
    std::vector<std::vector<bool> > passEventSelection(eventVariations.size()); // allocate memory here without push_back to avoid copying which confuses tree->Branch
    tout->Branch("uniqueID", &uniqueIDCopy);
    tout->Branch("Thrust", &Thrust);
    tout->Branch("TotalChgEnergy", &TotalChgEnergy);
    tout->Branch("NTrk", &NTrk);
    tout->Branch("Neu", &Neu);
    tout->Branch("STheta", &STheta);
    for (int iV = 0; iV < eventVariations.size(); iV++) {
      tout->Branch(("passEventSelection_" + std::to_string(iV)).c_str(), &passEventSelection.at(iV));
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
      TotalChgEnergy.clear();
      NTrk.clear();
      Neu.clear();
      STheta.clear();
      Thrust.clear();
      for (int iV = 0; iV < trackVariations.size(); iV++) {
        if (genTree && iV > 0 ) break;
        selectedParts.at(iV) = 0;
        selectedPx.at(iV).clear();
        selectedPy.at(iV).clear();
        selectedPz.at(iV).clear();
        selectedPwflag.at(iV).clear();
        TotalChgEnergy.push_back(0);
        NTrk.push_back(0);
        Neu.push_back(0);
      }
      for (auto &ev : passEventSelection) ev.clear();

      // compute event selection variables
      for (int iP = 0; iP < nParticle; iP++) {

        if (debug) std::cout << TString::Format("iP %d, pwflag %d, theta %f, pt %f, d0 %f, z0 %f, ntpc %d", iP, pwflag[iP], theta[iP], pt[iP], d0[iP], z0[iP], ntpc[iP]) << std::endl;

        // loop over variations
        for (int iV = 0; iV < trackVariations.size(); iV++) {

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
            TotalChgEnergy.at(iV) += TMath::Sqrt(pmag[iP] * pmag[iP] + mass[iP] * mass[iP]);
            NTrk.at(iV) += 1;
            if (debug) std::cout << "Passed charged track selection" << std::endl;
          }

          // count neutral tracks
          bool neutralTrackSelections =
            (pwflag[iP] == 4 || pwflag[iP] == 5)
            && TMath::Sqrt(pmag[iP] * pmag[iP] + mass[iP] * mass[iP]) >= trackVariations.at(iV)["ECut"]
            && TMath::Abs(cos(theta[iP])) <= trackVariations.at(iV)["neutralTracksAbsCosThCut"];
          if (neutralTrackSelections) {
            Neu.at(iV) += 1;
            if (debug) std::cout << "Passed neutral track selection" << std::endl;
          }

          // add to input list for thrust. check for -1 which indicates use all tracks for thrust
          bool keeptrack = tree == "tgen" || tree == "tgenBefore"; // generator level always keep
          keeptrack = keeptrack || trackVariations.at(iV)["applyTrackSelection"] == 0 ; // track selection should not be applied
          keeptrack = keeptrack || chargedTrackSelections || neutralTrackSelections; // passes selection
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
      for (int iV = 0; iV < trackVariations.size(); iV++) {

        // if tgen or tgenbefore only do the first nominal variation
        if (genTree && iV > 0 ) break;

        // sphericity
        spher = std::make_unique<Sphericity>(Sphericity(selectedParts.at(iV), selectedPx.at(iV).data(), selectedPy.at(iV).data(), selectedPz.at(iV).data(), selectedPwflag.at(iV).data(), false));
        STheta.push_back(spher->sphericityAxis().Theta());
        // thrust
        thrust = getThrust(selectedParts.at(iV), selectedPx.at(iV).data(), selectedPy.at(iV).data(), selectedPz.at(iV).data(), THRUST::OPTIMAL); //, false, false, pDataReader.weight);
        Thrust.push_back(thrust.Mag());

        // compute event selection passes
        for (int iEV = 0; iEV < eventVariations.size(); iEV++) {
          passEventSelection.at(iEV).push_back(
            passesNTupleAfterCut == 1
            && TotalChgEnergy.at(iV) >= eventVariations.at(iEV)["TotalChgEnergyCut"]
            && TMath::Abs(TMath::Cos(STheta.at(iV))) <= eventVariations.at(iEV)["AbsCosSThetaCut"]
            && NTrk.at(iV) >= eventVariations.at(iEV)["NTrkCut"]
            && (NTrk.at(iV) + Neu.at(iV)) >= eventVariations.at(iEV)["NeuNchCut"]
          );
        }

      }
      if (debug) std::cout << TString::Format("Nominal STheta %f, Thrust %f", STheta.at(0), Thrust.at(0)) << std::endl;

      tout->Fill();
    }

    tout->Write();
    std::cout << "\n" << std::endl;
  }

  return 1;
}