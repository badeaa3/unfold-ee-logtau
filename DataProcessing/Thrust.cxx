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

// default selection
std::map<std::string, float> getSelection(
 /* charged track selections */
 int nTPCcut, // number of hits in the TPC
 float chargedTracksAbsCosThCut, // cosine of polar angle for charged tracks
 float ptCut, // transverse momentum of charged tracks [GeV]
 float d0Cut, // transverse impact parameter of charged tracks [cm]
 float z0Cut, // lonitudinal impact parameter of charged tracks [cm]
 /* neutral particle selections */
 float ECut, // energy of neutral ECAL/HCAL objects [GeV]
 float neutralTracksAbsCosThCut, // cosine of polar angle for neutral objects
 /* event selections */
 float TotalTrkEnergyCut, // total charged energy [GeV]
 float AbsCosSThetaCut, // cosine of polar angle of sphericity axis 
 float NTrkCut, // number of charged tracks
 float NeuNchCut, // number of charged and neutral particles
 float EVisCut, // total visible energy [GeV]
 /* thrust variations */
 bool keepChargedTracks, // include charged tracks in calculation
 bool keepNeutralTracks,  // include neutral tracks in calculation
 bool doMET // include missing momentum vector in calculation
 ) {
  return std::map<std::string, float> {
    {"nTPCcut", nTPCcut},
    {"chargedTracksAbsCosThCut", chargedTracksAbsCosThCut},
    {"ptCut", ptCut},
    {"d0Cut", d0Cut},
    {"z0Cut", z0Cut},
    {"ECut", ECut},
    {"neutralTracksAbsCosThCut", neutralTracksAbsCosThCut},
    {"TotalTrkEnergyCut", TotalTrkEnergyCut},
    {"AbsCosSThetaCut", AbsCosSThetaCut},
    {"NTrkCut", NTrkCut},
    {"NeuNchCut", NeuNchCut},
    {"EVisCut", EVisCut},
    {"keepChargedTracks", keepChargedTracks},
    {"keepNeutralTracks", keepNeutralTracks},
    {"doMET", doMET},
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

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Selection Variations %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::vector<std::map<std::string, float> > selections; // vector of variations

  // default/nominal values
  int d_nTPC = 4;
  float d_AbsCosThetaChg = 0.94;
  float d_ptChg = 0.2; // GeV 
  float d_d0 = 2; // cm
  float d_z0 = 10; // cm
  /* neutral track selections */
  float d_ENeu = 0.4; // GeV
  float d_AbsCosThetaNeu = 0.98;
  /* event selections */
  float d_Ech = 15; // GeV 
  float d_AbsCosSTheta = 0.82;
  float d_NTrk = 5;
  float d_NTrkPlusNeu = 13;
  float d_EVis = 0; // GeV
  /* thrust variations */
  bool d_ThrCh = true; // include charged tracks in calculation
  bool d_ThrNeu = true;  // include neutral tracks in calculation
  bool d_ThrMissP = false; // include missing momentum vector in calculation

  // push back the selections
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, 0, 1, 0, 0, 0, d_ThrCh, d_ThrNeu, d_ThrMissP)); // no event selections
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, d_ThrNeu, d_ThrMissP)); // nominal
  selections.push_back(getSelection(7,      d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, d_ThrNeu, d_ThrMissP)); // ntpc 4 -> 7
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, 0.4,     d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, d_ThrNeu, d_ThrMissP)); // charged tracks pT 0.2 -> 0.4 GeV
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, 10,    d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, d_ThrNeu, d_ThrMissP)); // total charged energy 15 -> 10 GeV
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, false,    d_ThrMissP)); // thrust without neutral objects
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, d_EVis, d_ThrCh, d_ThrNeu, true));       // thrust with missing momentum vector as object
  selections.push_back(getSelection(d_nTPC, d_AbsCosThetaChg, d_ptChg, d_d0, d_z0, d_ENeu, d_AbsCosThetaNeu, d_Ech, d_AbsCosSTheta, d_NTrk, d_NTrkPlusNeu, 91.2/2, d_ThrCh, d_ThrNeu, d_ThrMissP)); // visible energy selection 0 -> 0.5*Ecm GeV
  
  // vectors for selected objects
  std::vector<int> selectedParts;
  std::vector<std::vector<float> > selectedPx, selectedPy, selectedPz;
  std::vector<std::vector<Short_t> > selectedPwflag;
  
  // event level quantities
  TVector3 thrust;
  std::unique_ptr<Sphericity> spher;

  // save variation definitions to a tree
  std::unique_ptr<TTree> varDefs (new TTree("Selections", ""));
  
  // selection variables
  int s_nTPC;
  float s_AbsCosThetaChg;
  float s_ptChg;
  float s_d0;
  float s_z0;
  /* neutral track selections */
  float s_ENeu;
  float s_AbsCosThetaNeu;
  /* event selections */
  float s_Ech;
  float s_AbsCosSTheta;
  float s_NTrk;
  float s_NTrkPlusNeu;
  float s_EVis;
  /* thrust variations */
  bool s_ThrCh;
  bool s_ThrNeu;
  bool s_ThrMissP;
  
  // set branches
  varDefs->Branch("nTPC", &s_nTPC);
  varDefs->Branch("AbsCosThetaChg", &s_AbsCosThetaChg);
  varDefs->Branch("ptChg", &s_ptChg);
  varDefs->Branch("d0", &s_d0);
  varDefs->Branch("z0", &s_z0);
  varDefs->Branch("ENeu", &s_ENeu);
  varDefs->Branch("AbsCosThetaNeu", &s_AbsCosThetaNeu);
  varDefs->Branch("Ech", &s_Ech);
  varDefs->Branch("AbsCosSTheta", &s_AbsCosSTheta);
  varDefs->Branch("NTrk", &s_NTrk);
  varDefs->Branch("NTrkPlusNeu", &s_NTrkPlusNeu);
  varDefs->Branch("EVis", &s_EVis);
  varDefs->Branch("ThrCh", &s_ThrCh);
  varDefs->Branch("ThrNeu", &s_ThrNeu);
  varDefs->Branch("ThrMissP", &s_ThrMissP);
  
  // push back for each variation and save to tree
  for (unsigned int iV = 0; iV < selections.size(); iV++) {

    // push back holders for selected particles to perform calculations
    selectedParts.push_back(0);
    selectedPx.push_back(std::vector<float>());
    selectedPy.push_back(std::vector<float>());
    selectedPz.push_back(std::vector<float>());
    selectedPwflag.push_back(std::vector<Short_t>());
      
    // write selections
    s_nTPC = selections.at(iV)["nTPCcut"];
    s_AbsCosThetaChg = selections.at(iV)["chargedTracksAbsCosThCut"];
    s_ptChg = selections.at(iV)["ptCut"];
    s_d0 = selections.at(iV)["d0Cut"];
    s_z0 = selections.at(iV)["z0Cut"];
    s_ENeu = selections.at(iV)["ECut"];
    s_AbsCosThetaNeu = selections.at(iV)["neutralTracksAbsCosThCut"];
    s_Ech = selections.at(iV)["TotalTrkEnergyCut"];
    s_AbsCosSTheta = selections.at(iV)["AbsCosSThetaCut"];
    s_NTrk = selections.at(iV)["NTrkCut"];
    s_NTrkPlusNeu = selections.at(iV)["NeuNchCut"];
    s_EVis = selections.at(iV)["EVisCut"];
    s_ThrCh = selections.at(iV)["keepChargedTracks"];
    s_ThrNeu = selections.at(iV)["keepNeutralTracks"];
    s_ThrMissP = selections.at(iV)["doMET"];

    // fill tree
    varDefs->Fill();
  }
  // write tree
  varDefs->Write();
    
  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Loop %%%%%%%%%%%%%%%%%%%%%%%%%%#

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
    std::vector<bool> passEventSelection;
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
    tout->Branch("passEventSelection", &passEventSelection);

    // create object level histograms for each pwflag
    std::map<std::pair<int, std::string>, TH1D*> hists;
    for(int iP=0; iP <= 5; iP++){
      hists[{iP, "cosTheta"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "cosTheta").c_str(), ";cos#theta;Entries", 200, -1, 1);
      hists[{iP, "phi"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "phi").c_str(), ";#phi;Entries", 800, -4, 4);
      hists[{iP, "pt"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "pt").c_str(), ";p_{T} [GeV];Entries", 100, 0, 100);
      hists[{iP, "ntpc"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "ntpc").c_str(), ";N_{TPC};Entries", 31, -0.5, 30.5);
      hists[{iP, "d0"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "d0").c_str(), ";d_{0} [cm];Entries", 50, -2.5, 2.5);
      hists[{iP, "z0"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "z0").c_str(), ";z_{0} [cm];Entries", 300, -15, 15);
      hists[{iP, "pmag"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "pmag").c_str(), ";|#vec{p}| [GeV];Entries", 100, 0, 100);
      hists[{iP, "mass"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "mass").c_str(), ";Mass [GeV];Entries", 100, 0, 10);
      hists[{iP, "energy"}] = new TH1D( (tree + "_hist_pwflag" + std::to_string(iP) + "_" + "energy").c_str(), ";Energy [GeV];Entries", 100, 0, 100);
    }

    // create event level histograms for each selection variation
    for (unsigned int iV = 0; iV < selections.size(); iV++){
      if (genTree && iV > 0 ) break;
      hists[{iV, "ntrk"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "ntrk").c_str(), ";N_{Trk};Entries", 61, -0.5, 60.5);
      hists[{iV, "nneu"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "nneu").c_str(), ";N_{Neu};Entries", 51, -0.5, 50.5);
      hists[{iV, "ntrkPlusNeu"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "ntrkPlusNeu").c_str(), ";N_{Trk+Neu};Entries", 81, -0.5, 80.5);
      hists[{iV, "eCh"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "eCh").c_str(), ";E_{Ch} [GeV];Entries", 200, 0, 200);
      hists[{iV, "cosThetaSph"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "cosThetaSph").c_str(), ";cos#theta_{Sph};Entries", 100, -1, 1);
      hists[{iV, "sphericity"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "sphericity").c_str(), ";Sphericity;Entries", 100, 0, 1);
      hists[{iV, "thrust"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "thrust").c_str(), ";Thrust;Entries", 100, 0.5, 1);
      hists[{iV, "missP"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "missP").c_str(), ";|#vec{p}_{MET}| [GeV];Entries", 100, 0, 100);
      hists[{iV, "evis"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "evis").c_str(), ";E_{Vis} [GeV];Entries", 200, 0, 200);
      hists[{iV, "cosThetaThrust"}] = new TH1D( (tree + "_hist_sel" + std::to_string(iV) + "_" + "cosThetaThrust").c_str(), ";cos#theta_{Thr};Entries", 100, -1, 1);
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
      passEventSelection.clear();
      for (unsigned int iV = 0; iV < selections.size(); iV++) {
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
      
      // loop over selections
      for (unsigned int iV = 0; iV < selections.size(); iV++) {

        // if tgen or tgenbefore only do the first nominal variation
        if (genTree && iV > 0 ) break;

        // loop over particles
        for (int iP = 0; iP < nParticle; iP++) {

          if (debug) std::cout << TString::Format("iP %d, pwflag %d, theta %f, pt %f, d0 %f, z0 %f, ntpc %d", iP, pwflag[iP], theta[iP], pt[iP], d0[iP], z0[iP], ntpc[iP]) << std::endl;

          // compute the particle energy
          float energy = TMath::Sqrt(pmag[iP] * pmag[iP] + mass[iP] * mass[iP]);
    
          // fill histogram for all particles on first pass
          if(iV == 0 && pwflag[iP] >= 0 && pwflag[iP] <= 5){
            hists[{pwflag[iP], "cosTheta"}]->Fill(cos(theta[iP]));
            hists[{pwflag[iP], "phi"}]->Fill(phi[iP]);
            hists[{pwflag[iP], "pt"}]->Fill(pt[iP]);
            hists[{pwflag[iP], "ntpc"}]->Fill(ntpc[iP]);
            hists[{pwflag[iP], "d0"}]->Fill(d0[iP]);
            hists[{pwflag[iP], "z0"}]->Fill(z0[iP]);
            hists[{pwflag[iP], "pmag"}]->Fill(pmag[iP]);
            hists[{pwflag[iP], "mass"}]->Fill(mass[iP]);
            hists[{pwflag[iP], "energy"}]->Fill(energy);
          }

          // charged track selection
          bool passChgTrkSel =
            (pwflag[iP] >= 0 && pwflag[iP] <= 2)
            && (TMath::Abs(cos(theta[iP])) <= selections.at(iV)["chargedTracksAbsCosThCut"])
            && (pt[iP] >= selections.at(iV)["ptCut"])
            && (TMath::Abs(d0[iP]) <= selections.at(iV)["d0Cut"])
            && (TMath::Abs(z0[iP]) <= selections.at(iV)["z0Cut"])
            && (ntpc[iP] >= selections.at(iV)["nTPCcut"]);
          
          // populate
          if (passChgTrkSel) {
            if (debug) std::cout << "Passed charged track selection" << std::endl;
            // increment values
            TotalTrkEnergy.at(iV) += energy;
	          EVis.at(iV) += energy;
            NTrk.at(iV) += 1;
            // add to input list for sphericity and thrust
            if(selections.at(iV)["keepChargedTracks"]){
              selectedParts.at(iV) += 1;
              selectedPx.at(iV).push_back(px[iP]);
              selectedPy.at(iV).push_back(py[iP]);
              selectedPz.at(iV).push_back(pz[iP]);
              selectedPwflag.at(iV).push_back(pwflag[iP]);
            }
          }

          // neutral particle selection
          bool passNeuPartSel = 
            (pwflag[iP] == 4 || pwflag[iP] == 5)
	          && (energy >= selections.at(iV)["ECut"])
            && (TMath::Abs(cos(theta[iP])) <= selections.at(iV)["neutralTracksAbsCosThCut"]);
          
          // populate
          if (passNeuPartSel) {
            if (debug) std::cout << "Passed neutral track selection" << std::endl;
            // increment values
            EVis.at(iV) += energy;
            Neu.at(iV) += 1;
            // add to input list for sphericity and thrust
            if(selections.at(iV)["keepNeutralTracks"]){
              selectedParts.at(iV) += 1;
              selectedPx.at(iV).push_back(px[iP]);
              selectedPy.at(iV).push_back(py[iP]);
              selectedPz.at(iV).push_back(pz[iP]);
              selectedPwflag.at(iV).push_back(pwflag[iP]);
            }
          }
        }

        // sphericity
        spher = std::make_unique<Sphericity>(Sphericity(selectedParts.at(iV), selectedPx.at(iV).data(), selectedPy.at(iV).data(), selectedPz.at(iV).data(), selectedPwflag.at(iV).data(), false));
        STheta.push_back(spher->sphericityAxis().Theta());
	      Sph.push_back(spher->sphericity());

        // calculate the missing momentum vector
	      TVector3 met = TVector3(0, 0, 0);
	      for (int t = 0; t < selectedParts.at(iV); t++) {
	        met += (TVector3(selectedPx.at(iV).at(t), selectedPy.at(iV).at(t), selectedPz.at(iV).at(t)));
	      }
	      met = -met;
	      MissP.push_back(met.Mag());

        // include missing momentum vector in thrust calculation
        if (selections.at(iV)["doMET"]) {
          selectedParts.at(iV) += 1;
          selectedPx.at(iV).push_back(met.X()); // X() same as Px() for TVector3
          selectedPy.at(iV).push_back(met.Y()); // Y() same as Py() for TVector3
          selectedPz.at(iV).push_back(met.Z()); // Z() same as Pz() for TVector3
          selectedPwflag.at(iV).push_back(-1); // not important for thrust but just to keep vectors the same length use -1 for missing momentum vector
        }

        // thrust
	      // TVector3 getThrust(int n, float *px, float *py, float *pz, THRUST::algorithm algo = THRUST::HERWIG, bool doWeight = false, bool doInvertWeight = false, float* weight = NULL, bool doMET = false, Short_t *pwflag = NULL)
        thrust = getThrust(selectedParts.at(iV), selectedPx.at(iV).data(), selectedPy.at(iV).data(), selectedPz.at(iV).data(), THRUST::OPTIMAL);
        Thrust.push_back(thrust.Mag());
	TTheta.push_back(thrust.Theta());

        // compute event selection passes
	bool eventSelection =
	  passesNTupleAfterCut == 1
	  && (TotalTrkEnergy.at(iV) >= selections.at(iV)["TotalTrkEnergyCut"])
          && (TMath::Abs(TMath::Cos(STheta.at(iV))) <= selections.at(iV)["AbsCosSThetaCut"])
          && (NTrk.at(iV) >= selections.at(iV)["NTrkCut"])
          && ((NTrk.at(iV) + Neu.at(iV)) >= selections.at(iV)["NeuNchCut"]);

	// append and fill histograms if selection passed
	passEventSelection.push_back(eventSelection);
	if(eventSelection){
	  hists[{iV, "ntrk"}]->Fill(NTrk.at(iV));
	  hists[{iV, "nneu"}]->Fill(Neu.at(iV));
	  hists[{iV, "ntrkPlusNeu"}]->Fill(NTrk.at(iV) + Neu.at(iV));
	  hists[{iV, "eCh"}]->Fill(TotalTrkEnergy.at(iV));
	  hists[{iV, "cosThetaSph"}]->Fill(TMath::Cos(STheta.at(iV)));
	  hists[{iV, "sphericity"}]->Fill(Sph.at(iV));
	  hists[{iV, "thrust"}]->Fill(Thrust.at(iV));
	  hists[{iV, "missP"}]->Fill(MissP.at(iV));
	  hists[{iV, "evis"}]->Fill(EVis.at(iV));
	  hists[{iV, "cosThetaThrust"}]->Fill(TMath::Cos(TTheta.at(iV)));
	}

      }

      // fill per event
      tout->Fill();
    }

    // write tree
    tout->Write();

    // write histograms
    for (auto &entry : hists) {
        entry.second->Write();
    }
    
    std::cout << "\n" << std::endl;

  }

  return 1;
}
