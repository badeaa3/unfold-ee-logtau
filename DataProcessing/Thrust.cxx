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
#include <sstream>
#include <chrono>

// Parse configuration string from command line (format: key1=value1;key2=value2;...)
// Returns both selection parameters and string parameters
struct ConfigParams {
  std::map<std::string, float> selections;
  std::map<std::string, std::string> strings;
};

ConfigParams parseConfigString(const std::string& configStr) {
  ConfigParams result;
  
  // Split by semicolon
  std::stringstream ss(configStr);
  std::string item;
  
  while (std::getline(ss, item, ';')) {
    // Find the equals sign
    size_t pos = item.find('=');
    if (pos == std::string::npos) continue;
    
    std::string key = item.substr(0, pos);
    std::string valueStr = item.substr(pos + 1);
    
    // Special handling for string parameters
    if (key == "inFileType" || key == "treeNames") {
      result.strings[key] = valueStr;
    } else {
      // Convert value to float for selection parameters
      float value;
      if (valueStr == "true") {
        value = true;
      } else if (valueStr == "false") {
        value = false;
      } else {
        value = std::stof(valueStr);
      }
      result.selections[key] = value;
    }
  }
  
  return result;
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
  std::string configStr = "";
  int doNEvents = -1;
  bool debug = false;
  size_t divide    = 1; // how many divisions the file is being cut into, which division you are on, and how many events per division
  size_t thisdiv   = 0;
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "-i", 2) == 0) inFileName = argv[i + 1];
    if (strncmp(argv[i], "-o", 2) == 0) outFileName = argv[i + 1];
    if (strncmp(argv[i], "-c", 2) == 0) configStr = argv[i + 1];
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
  
  // Parse configuration string from command line
  if (configStr.empty()) {
    std::cout << "No configuration string provided. Exiting" << std::endl;
    return 0;
  }
  ConfigParams config = parseConfigString(configStr);
  std::map<std::string, float> selMap = config.selections;
  
  // Extract file type and tree names from configuration
  std::string inFileType = config.strings["inFileType"];
  std::string treeNamesStr = config.strings["treeNames"];
  
  // Parse tree names (comma-separated)
  std::vector<std::string> treeNames;
  std::stringstream ss(treeNamesStr);
  std::string treeName;
  while (std::getline(ss, treeName, ',')) {
    treeNames.push_back(treeName);
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
  Short_t charge[maxPart];

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Output File %%%%%%%%%%%%%%%%%%%%%%%%%%#
  std::unique_ptr<TFile> fout (new TFile(outFileName.c_str(), "RECREATE"));

  // #%%%%%%%%%%%%%%%%%%%%%%%%%% Selection Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  // Save selection definition to a tree
  std::unique_ptr<TTree> varDefs (new TTree("Selection", ""));
  
  // Selection variables - using selMap directly
  int s_nTPC = selMap["nTPCcut"];
  float s_AbsCosThetaChg = selMap["chargedTracksAbsCosThCut"];
  float s_ptChg = selMap["ptCut"];
  float s_d0 = selMap["d0Cut"];
  float s_z0 = selMap["z0Cut"];
  float s_ENeu = selMap["ECut"];
  float s_AbsCosThetaNeu = selMap["neutralTracksAbsCosThCut"];
  float s_Ech = selMap["TotalTrkEnergyCut"];
  float s_AbsCosSTheta = selMap["AbsCosSThetaCut"];
  float s_NTrk = selMap["NTrkCut"];
  float s_NTrkPlusNeu = selMap["NeuNchCut"];
  float s_EVis = selMap["EVisCut"];
  float s_MissP = selMap["MissPCut"];
  bool s_ThrCh = selMap["keepChargedTracks"];
  bool s_ThrNeu = selMap["keepNeutralTracks"];
  bool s_ThrMissP = selMap["doMET"];
  
  // Set branches
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
  varDefs->Branch("MissP", &s_MissP);
  varDefs->Branch("ThrCh", &s_ThrCh);
  varDefs->Branch("ThrNeu", &s_ThrNeu);
  varDefs->Branch("ThrMissP", &s_ThrMissP);
  
  // Fill and write tree (single entry)
  varDefs->Fill();
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
    t->SetBranchAddress("charge", &charge);

    // event level quantities
    TVector3 thrust;
    std::unique_ptr<Sphericity> spher;
    
    // vectors for selected objects (single selection now)
    int selectedParts = 0;
    std::vector<float> selectedPx, selectedPy, selectedPz;
    std::vector<Short_t> selectedPwflag;
    
    // create output tree
    std::unique_ptr<TTree> tout (new TTree(tree.c_str(), ""));
    unsigned long long uniqueIDCopy; 
    float Thrust, TotalTrkEnergy, STheta, Sph, MissP, EVis, TTheta;
    int NTrk, Neu;
    bool passEventSelection;
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
    // save selected particles
    tout->Branch("NSelectedParticles", &selectedParts);
    tout->Branch("px", &selectedPx);
    tout->Branch("py", &selectedPy);
    tout->Branch("pz", &selectedPz);
    tout->Branch("pwflag", &selectedPwflag);

    std::vector<float> conversionElectronTheta, conversionElectronPhi, conversionElectronPt;
    tout->Branch("conversionElectronTheta", &conversionElectronTheta);
    tout->Branch("conversionElectronPhi", &conversionElectronPhi);
    tout->Branch("conversionElectronPt", &conversionElectronPt);
    
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

    // create event level histograms for single selection
    hists[{0, "ntrk"}] = new TH1D( (tree + "_hist_sel_ntrk").c_str(), ";N_{Trk};Entries", 61, -0.5, 60.5);
    hists[{0, "nneu"}] = new TH1D( (tree + "_hist_sel_nneu").c_str(), ";N_{Neu};Entries", 51, -0.5, 50.5);
    hists[{0, "ntrkPlusNeu"}] = new TH1D( (tree + "_hist_sel_ntrkPlusNeu").c_str(), ";N_{Trk+Neu};Entries", 81, -0.5, 80.5);
    hists[{0, "eCh"}] = new TH1D( (tree + "_hist_sel_eCh").c_str(), ";E_{Ch} [GeV];Entries", 200, 0, 200);
    hists[{0, "cosThetaSph"}] = new TH1D( (tree + "_hist_sel_cosThetaSph").c_str(), ";cos#theta_{Sph};Entries", 100, -1, 1);
    hists[{0, "sphericity"}] = new TH1D( (tree + "_hist_sel_sphericity").c_str(), ";Sphericity;Entries", 100, 0, 1);
    hists[{0, "thrust"}] = new TH1D( (tree + "_hist_sel_thrust").c_str(), ";Thrust;Entries", 100, 0.5, 1);
    hists[{0, "logtau"}] = new TH1D( (tree + "_hist_sel_logtau").c_str(), ";log(#tau);Entries", 100, -10, 0);
    hists[{0, "missP"}] = new TH1D( (tree + "_hist_sel_missP").c_str(), ";|#vec{p}_{MET}| [GeV];Entries", 100, 0, 100);
    hists[{0, "evis"}] = new TH1D( (tree + "_hist_sel_evis").c_str(), ";E_{Vis} [GeV];Entries", 200, 0, 200);
    hists[{0, "cosThetaThrust"}] = new TH1D( (tree + "_hist_sel_cosThetaThrust").c_str(), ";cos#theta_{Thr};Entries", 100, -1, 1);

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

      // reset variables for single selection
      TotalTrkEnergy = 0;
      NTrk = 0;
      Neu = 0;
      EVis = 0;
      conversionElectronTheta.clear();
      conversionElectronPhi.clear();
      conversionElectronPt.clear();
      selectedParts = 0;
      selectedPx.clear();
      selectedPy.clear();
      selectedPz.clear();
      selectedPwflag.clear();

      // loop over particles
      for (int iP = 0; iP < nParticle; iP++) {

        // debug printing
        if (debug) std::cout << TString::Format("iP %d, pwflag %d, theta %f, pt %f, d0 %f, z0 %f, ntpc %d, charge %i", iP, pwflag[iP], theta[iP], pt[iP], d0[iP], z0[iP], ntpc[iP], charge[iP]) << std::endl;

        // compute the particle energy
        float energy = TMath::Sqrt(pmag[iP] * pmag[iP] + mass[iP] * mass[iP]);

        // determine if good generator level
        bool goodGenPart = true;
        // neutral cleaning around phi = 0 for photon radiation along beam pipe
        if (!genTree) goodGenPart = false;
        if (genTree && inFileType == "ALEPHMC" && charge[iP] == 0 && std::abs(phi[iP]) <= 0.001 && pt[iP] > 0.00099 && pt[iP] < 0.001009){
          goodGenPart = false;
        }

        // determine if conversion electron
        bool isConversionElectron = true;
        //both electrons
        if( pwflag[iP] != 2 ) isConversionElectron = false;
        if( pwflag[iP-1] != 2 ) isConversionElectron = false;
        //opposite charge required
        if( charge[iP] != -(charge[iP-1])) isConversionElectron = false;
        //dtheta and dphi matching
        float conversionDPhi = 0.05;
        float conversionDTheta = 0.05;
        if( TMath::Abs(theta[iP] - theta[iP-1]) > conversionDTheta) isConversionElectron = false;
        if( TMath::ACos(TMath::Cos(phi[iP] - phi[iP-1])) > conversionDPhi) isConversionElectron = false;

        // apply reco level selections
        bool passChgTrkSel = false;
        bool passNeuPartSel = false;
        if (!genTree){

          // charged particle selections
          passChgTrkSel =
            (pwflag[iP] >= 0 && pwflag[iP] <= 2)
            && (TMath::Abs(cos(theta[iP])) <= selMap["chargedTracksAbsCosThCut"])
            && (pt[iP] >= selMap["ptCut"])
            && (TMath::Abs(d0[iP]) <= selMap["d0Cut"])
            && (TMath::Abs(z0[iP]) <= selMap["z0Cut"])
            && (ntpc[iP] >= selMap["nTPCcut"]);

          // neutral particle selections
          passNeuPartSel =
            (pwflag[iP] == 4 || pwflag[iP] == 5)
            && (energy >= selMap["ECut"])
            && (TMath::Abs(cos(theta[iP])) <= selMap["neutralTracksAbsCosThCut"]);
        }

        // fill and save
        bool saveParticle = false;
        if(goodGenPart) saveParticle = true;
        if(isConversionElectron){
          // if conversion electron then don't save particle and remove previous particle also
          saveParticle = false;
          // remove previous particle
          selectedParts -= 1;
          selectedPx.pop_back();
          selectedPy.pop_back();
          selectedPz.pop_back();
          selectedPwflag.pop_back();
          // conversion electrons
          conversionElectronTheta.push_back(theta[iP]);
          conversionElectronPhi.push_back(phi[iP]);
          conversionElectronPt.push_back(pt[iP]);
        }
        if(passChgTrkSel && selMap["keepChargedTracks"]){
          if (debug) std::cout << "Passed charged track selection" << std::endl;
          saveParticle = true;
          TotalTrkEnergy += energy;
          EVis += energy;
          NTrk += 1;
        }
        if(passNeuPartSel && selMap["keepNeutralTracks"]){
          if (debug) std::cout << "Passed neutral track selection" << std::endl;
          saveParticle = true;
          EVis += energy;
          Neu += 1;
        }

        // save the particle
        if (saveParticle){
          // save for event shape variables
          selectedParts += 1;
          selectedPx.push_back(px[iP]);
          selectedPy.push_back(py[iP]);
          selectedPz.push_back(pz[iP]);
          selectedPwflag.push_back(pwflag[iP]);
          // fill particle kinematic histograms
          if(!genTree) continue;
          if(inFileType == "PYTHIA8") continue; // no pwflag for pythia8, only pdgid
          if(!(pwflag[iP] >= 0 && pwflag[iP] <= 5)) continue; // only save pwflag 0-5
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

      }

      // sphericity
      spher = std::make_unique<Sphericity>(Sphericity(selectedParts, selectedPx.data(), selectedPy.data(), selectedPz.data(), selectedPwflag.data(), false));
      STheta = spher->sphericityAxis().Theta();
      Sph = spher->sphericity();

      // calculate the missing momentum vector
      TVector3 met = TVector3(0, 0, 0);
      for (int t = 0; t < selectedParts; t++) {
        met += (TVector3(selectedPx.at(t), selectedPy.at(t), selectedPz.at(t)));
      }
      met = -met;
      MissP = met.Mag();

      // include missing momentum vector in thrust calculation
      if (selMap["doMET"] > 0.5) {
        selectedParts += 1;
        selectedPx.push_back(met.X()); // X() same as Px() for TVector3
        selectedPy.push_back(met.Y()); // Y() same as Py() for TVector3
        selectedPz.push_back(met.Z()); // Z() same as Pz() for TVector3
        selectedPwflag.push_back(-1); // not important for thrust but just to keep vectors the same length use -1 for missing momentum vector
      }

      // thrust
      // TVector3 getThrust(int n, float *px, float *py, float *pz, THRUST::algorithm algo = THRUST::HERWIG, bool doWeight = false, bool doInvertWeight = false, float* weight = NULL, bool doMET = false, Short_t *pwflag = NULL)
      thrust = getThrust(selectedParts, selectedPx.data(), selectedPy.data(), selectedPz.data(), THRUST::OPTIMAL);
      Thrust = thrust.Mag();
      TTheta = thrust.Theta();

      // compute event selection passes
      bool eventSelection =
      passesNTupleAfterCut == 1
        && (TotalTrkEnergy >= selMap["TotalTrkEnergyCut"])
        && (TMath::Abs(TMath::Cos(STheta)) <= selMap["AbsCosSThetaCut"])
        && (NTrk >= selMap["NTrkCut"])
        && ((NTrk + Neu) >= selMap["NeuNchCut"])
        && (EVis >= selMap["EVisCut"])
        && (MissP < selMap["MissPCut"]);

      // gen always passes event selection
      eventSelection = eventSelection || genTree;
      
      // set passEventSelection and fill histograms if selection passed
      passEventSelection = eventSelection;
      if(eventSelection){
        hists[{0, "ntrk"}]->Fill(NTrk);
        hists[{0, "nneu"}]->Fill(Neu);
        hists[{0, "ntrkPlusNeu"}]->Fill(NTrk + Neu);
        hists[{0, "eCh"}]->Fill(TotalTrkEnergy);
        hists[{0, "cosThetaSph"}]->Fill(TMath::Cos(STheta));
        hists[{0, "sphericity"}]->Fill(Sph);
        hists[{0, "thrust"}]->Fill(Thrust);
        hists[{0, "logtau"}]->Fill(TMath::Log(1-Thrust));
        hists[{0, "missP"}]->Fill(MissP);
        hists[{0, "evis"}]->Fill(EVis);
        hists[{0, "cosThetaThrust"}]->Fill(TMath::Cos(TTheta));
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
