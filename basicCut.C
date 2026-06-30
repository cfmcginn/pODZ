//MISSING ELEMENTS: HF SELECTIONS FOR EVENTS

//c + cpp
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <unistd.h>

//ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/JSON_handler.h"

//Checking impact of basic cuts 
int basicCut()
{
  //Define some cut values here for ease of use
  const Float_t absVzMax = 15.0;
  const Float_t dMassFromPDG = 1.86484;//https://pdg.lbl.gov/2026/listings/rpp2026-list-D-zero.pdf
  const Float_t deltaDMass = 0.2;//200 mev around dmass for cut
  const Float_t maxDAbsY = 2.0;

  //D0 basic kinematics (secondary vertex / pointing related)
  const Float_t maxD0Alpha = 0.4;
  const Float_t minD0SecondaryDistOverErr = 3.5;
  const Float_t minD02ndVtxProb = 0.1;
  const Float_t d0MaxDeltaTheta = 0.3;
  
  const Float_t minDTrkPt = 1.0;  
  const Float_t maxDTrkAbsEta = 2.4;  
  
  //define some necessary filenames
  const std::string outFileName = "testBasicCut.root";
  std::string jsonName = "/home/cfm/Projects/pODZero/pO_golden.json";
  //  const std::string inFileName = "/home/cfm/MITHIG/tempPlots2026/Jun25/HiForestMINIAOD_numEvent50000.root";
  const std::string inFileName = "/home/cfm/MITHIG/tempPlots2026/Jun29/HiForestMINIAOD_numEvent50000.root";

  //Create my outfile and histogram for a test
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  const Int_t nD0Bins = 50;
  TH1D* outHist_p = new TH1D("d0Mass_h", ";D0 Mass [GeV];Counts in pO", nD0Bins, dMassFromPDG - deltaDMass, dMassFromPDG + deltaDMass); 
  
  //Variable declarations
  //hiTree
  UInt_t run, lumi;
  //  ULong64_t evt;
  Float_t vz;
  
  //hltTree
  Int_t HLT_OxyZeroBias_v1;
  Int_t HLT_MinimumBiasHF_OR_BptxAND_v1;

  //skimTree
  Int_t pprimaryVertexFilter;
  
  //dfinder
  const Int_t DsizeMax = 10000;  
  Int_t Dsize;
  Float_t Dmass[DsizeMax];
  Float_t Dpt[DsizeMax];
  Float_t Dy[DsizeMax];
  Float_t DsvpvDistance[DsizeMax];//svpv distance
  Float_t DsvpvDisErr[DsizeMax];//svpv distance uncertainty
  Float_t Dalpha[DsizeMax];//pointing of d0 relative secondary
  Float_t Ddtheta[DsizeMax];//opening angle of decay products
  Float_t Dchi2cl[DsizeMax];//probability secondary vtx
  //D0 track variables
  Float_t Dtrk1Pt[DsizeMax];
  Float_t Dtrk2Pt[DsizeMax];
  Float_t Dtrk1Eta[DsizeMax];
  Float_t Dtrk2Eta[DsizeMax];
  Bool_t Dtrk1highPurity[DsizeMax];
  Bool_t Dtrk2highPurity[DsizeMax];

  //Initialize the json for local use
  JSON_handler jsonIn(jsonName);  

  //Grab the root file and TTrees
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> ttreeNames = {"hltanalysis/HltTree",
					 "l1object/L1UpgradeFlatTree",
					 "ggHiNtuplizer/EventTree",
					 "muonAnalyzer/MuonTree",
					 "skimanalysis/HltTree",
					 "hiEvtAnalyzer/HiTree",
					 "zdcanalyzer/zdcdigi",
					 "ppTracks/trackTree",
					 "particleFlowAnalyser/pftree",
					 "ak0PFJetAnalyzer/t",
					 "Dfinder/ntDkpi"};

  //MAKE SURE YOU ADD TTREES TO THIS LIST ELSE YOUR CODE WILL NOT BE RELIABLE
  std::map<std::string,TTree*> ttreeList_p;

  //First pass of grabbing ttrees
  for(auto const & name : ttreeNames){
    ttreeList_p[name] = nullptr;//initialize
    ttreeList_p[name] = (TTree*)inFile_p->Get(name.c_str());

    ttreeList_p[name]->SetBranchStatus("*", 0);
  }

  //Assign a few additional pointers for convenience
  TTree* dTree_p = (TTree*)ttreeList_p["Dfinder/ntDkpi"];
  
  //Grab nEntries right away and do a quick nEntries check
  const Long64_t nEntries = (Long64_t)dTree_p->GetEntries();
  for(auto const & tree : ttreeList_p){
    Long64_t nEntriesTemp = (Long64_t)tree.second->GetEntries();

    if(nEntries != nEntriesTemp){
      std::cout << "TTree '" << tree.first << "' has nEntries mismatch (" << nEntriesTemp << ") with main TTree, '" << dTree_p->GetName() << "'. return 1" << std::endl;

      //Cleanup
      inFile_p->Close();
      delete inFile_p;
      
      delete outHist_p;
      outFile_p->Close();
      delete outFile_p;
      
      return 1;
    }
  }
  std::cout << "All nEntries passed" << std::endl;
  sleep(1);
  
  Int_t dTreeDSizeMax = dTree_p->GetMaximum("Dsize");
  if(dTreeDSizeMax >= DsizeMax){
    std::cerr << "ERROR: Dsize in file '" << inFileName << "', '" << dTreeDSizeMax << "', exceeds DsizeMax '" << DsizeMax << "'. Terminating" << std::endl;
    //clean inflie
    inFile_p->Close();
    delete inFile_p;

    //Clean outfile
    delete outHist_p;
    outFile_p->Close();
    delete outFile_p;

    return 1;
  }

  //hlt tree branches
  ttreeList_p["hltanalysis/HltTree"]->SetBranchStatus("HLT_OxyZeroBias_v1", 1);
  ttreeList_p["hltanalysis/HltTree"]->SetBranchStatus("HLT_MinimumBiasHF_OR_BptxAND_v1", 1);

  ttreeList_p["hltanalysis/HltTree"]->SetBranchAddress("HLT_OxyZeroBias_v1", &HLT_OxyZeroBias_v1);
  ttreeList_p["hltanalysis/HltTree"]->SetBranchAddress("HLT_MinimumBiasHF_OR_BptxAND_v1", &HLT_MinimumBiasHF_OR_BptxAND_v1);

  //skim tree branches
  ttreeList_p["skimanalysis/HltTree"]->SetBranchStatus("*", 0);
  ttreeList_p["skimanalysis/HltTree"]->SetBranchStatus("pprimaryVertexFilter", 1);

  ttreeList_p["skimanalysis/HltTree"]->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  
  //hi tree branches
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchStatus("*", 0);
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchStatus("run", 1);
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchStatus("lumi", 1);
  //  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchStatus("evt", 1);
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchStatus("vz", 1);
  
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchAddress("run", &run);
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchAddress("lumi", &lumi);
  //  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchAddress("evt", &evt);
  ttreeList_p["hiEvtAnalyzer/HiTree"]->SetBranchAddress("vz", &vz);

  //d tree branches
  dTree_p->SetBranchStatus("*", 0);
  dTree_p->SetBranchStatus("Dsize", 1);
  dTree_p->SetBranchStatus("Dmass", 1);
  dTree_p->SetBranchStatus("Dy", 1);
  dTree_p->SetBranchStatus("DsvpvDistance", 1);
  dTree_p->SetBranchStatus("DsvpvDisErr", 1);
  dTree_p->SetBranchStatus("Dalpha", 1);
  dTree_p->SetBranchStatus("Ddtheta", 1);
  dTree_p->SetBranchStatus("Dchi2cl", 1);
  dTree_p->SetBranchStatus("Dtrk1Pt", 1);
  dTree_p->SetBranchStatus("Dtrk2Pt", 1);
  dTree_p->SetBranchStatus("Dtrk1Eta", 1);
  dTree_p->SetBranchStatus("Dtrk2Eta", 1);
  dTree_p->SetBranchStatus("Dtrk1highPurity", 1);
  dTree_p->SetBranchStatus("Dtrk2highPurity", 1);

  dTree_p->SetBranchAddress("Dsize", &Dsize);
  dTree_p->SetBranchAddress("Dmass", Dmass);
  dTree_p->SetBranchAddress("Dy", Dy);
  dTree_p->SetBranchAddress("DsvpvDistance", DsvpvDistance);
  dTree_p->SetBranchAddress("DsvpvDisErr", DsvpvDisErr);
  dTree_p->SetBranchAddress("Dalpha", Dalpha);
  dTree_p->SetBranchAddress("Ddtheta", Ddtheta);
  dTree_p->SetBranchAddress("Dchi2cl", Dchi2cl);
  dTree_p->SetBranchAddress("Dtrk1Pt", Dtrk1Pt);
  dTree_p->SetBranchAddress("Dtrk2Pt", Dtrk2Pt);
  dTree_p->SetBranchAddress("Dtrk1Eta", Dtrk1Eta);
  dTree_p->SetBranchAddress("Dtrk2Eta", Dtrk2Eta);
  dTree_p->SetBranchAddress("Dtrk1highPurity", Dtrk1highPurity);
  dTree_p->SetBranchAddress("Dtrk2highPurity", Dtrk2highPurity);
  
  std::cout << "Running loop..." << std::endl;
  const Long64_t entryToPrint = 1000;

  //Declare some maps for stat and cut tracking
  std::map<UInt_t, std::map<UInt_t, Bool_t> > mapRunToLumi;

  //event selection map
  std::map<std::string, Long64_t> totalEventsPerEventCut;
  totalEventsPerEventCut["0: Total Events"] = nEntries;
  totalEventsPerEventCut["1: JSON Selection"] = 0;
  totalEventsPerEventCut["2: Trigger Selection"] = 0;
  totalEventsPerEventCut["3: pPrimaryVertex Selection"] = 0;
  totalEventsPerEventCut["4: fabs(vz) < 15.0 Selection"] = 0;
  totalEventsPerEventCut["5: Events containing D (all cut)"] = 0;

  //d selection map
  std::map<std::string, Long64_t> totalDPerD0Cut;
  totalDPerD0Cut["0: Total D0"] = 0;
  totalDPerD0Cut["1: D0 Mass Selection"] = 0;
  totalDPerD0Cut["2: D0 y Selection"] = 0;
  totalDPerD0Cut["3: D0 Trk Selection"] = 0;
  totalDPerD0Cut["4: D0 2nd Vtx Dist/Err Selection"] = 0;
  totalDPerD0Cut["5: D0 alpha Selection"] = 0;
  totalDPerD0Cut["6: D0 dTheta Selection"] = 0;
  totalDPerD0Cut["7: D0 2nd Vtx prob. Selection"] = 0;
  
  
  for(Long64_t entry = 0; entry < nEntries; ++entry){
    if(entry%entryToPrint == 0) std::cout << " Processing " << entry << "/" << nEntries << "..." << std::endl;
    for(auto const & tree : ttreeList_p){
      tree.second->GetEntry(entry);
    }

    //JSON sel., map fill, sel. fill
    if(!jsonIn.isGood(run, lumi)) continue;

    if(mapRunToLumi.count(run) == 0) mapRunToLumi[run] = {};
    (mapRunToLumi[run])[lumi] = true;

    ++(totalEventsPerEventCut["1: JSON Selection"]);
    
    //HLT Selections, just both triggers
    if(!HLT_MinimumBiasHF_OR_BptxAND_v1 && !HLT_OxyZeroBias_v1) continue;
    //    if(!HLT_MinimumBiasHF_OR_BptxAND_v1) continue;//just restrict to MB
    ++(totalEventsPerEventCut["2: Trigger Selection"]);

    //Vertex selections
    if(pprimaryVertexFilter == 0) continue;
    ++(totalEventsPerEventCut["3: pPrimaryVertex Selection"]);

    if(TMath::Abs(vz) >= absVzMax) continue;//selection is vz <, not < or equal to
    ++(totalEventsPerEventCut["4: fabs(vz) < 15.0 Selection"]);

    bool eventContainsDPassingCuts = false;
    totalDPerD0Cut["0: Total D0"] += Dsize;    
    for(Int_t dI = 0; dI < Dsize; ++dI){            
      //D mass cut
      if(TMath::Abs(dMassFromPDG - Dmass[dI]) > deltaDMass) continue;      
      ++(totalDPerD0Cut["1: D0 Mass Selection"]);      

      //D rapidity cut
      if(TMath::Abs(Dy[dI]) > maxDAbsY) continue;
      ++(totalDPerD0Cut["2: D0 y Selection"]);      
      
      //D trk cut block
      if(Dtrk1Pt[dI] < minDTrkPt) continue;
      if(Dtrk2Pt[dI] < minDTrkPt) continue;
      if(TMath::Abs(Dtrk1Eta[dI]) > maxDTrkAbsEta) continue;
      if(TMath::Abs(Dtrk2Eta[dI]) > maxDTrkAbsEta) continue;
      if(!Dtrk1highPurity[dI]) continue;
      if(!Dtrk2highPurity[dI]) continue;      
      ++(totalDPerD0Cut["3: D0 Trk Selection"]);      

      //D0 svpv/err sel
      if(DsvpvDistance[dI]/DsvpvDisErr[dI] < minD0SecondaryDistOverErr) continue;
      ++(totalDPerD0Cut["4: D0 2nd Vtx Dist/Err Selection"]);

      //D0 alpha cut
      if(Dalpha[dI] > maxD0Alpha) continue;
      ++(totalDPerD0Cut["5: D0 alpha Selection"]);

      //D0 dtheta cut
      if(Ddtheta[dI] > d0MaxDeltaTheta) continue;
      ++(totalDPerD0Cut["6: D0 dTheta Selection"]);
      
      //D0 2nd vtx prob
      if(Dchi2cl[dI] < minD02ndVtxProb) continue;
      ++(totalDPerD0Cut["7: D0 2nd Vtx prob. Selection"]);

      outHist_p->Fill(Dmass[dI]);
      
      eventContainsDPassingCuts = true;
    }

    //If d passes all cuts, add it to the counter
    if(eventContainsDPassingCuts) ++totalEventsPerEventCut["5: Events containing D (all cut)"];
  }

  inFile_p->Close();
  delete inFile_p;

  std::cout << "Loop complete." << std::endl;
  std::cout << "Run stats:" << std::endl;
  for(auto const & mapRun : mapRunToLumi){
    std::vector<UInt_t> lumiList;
    for(auto const & mapLumi : mapRun.second){
      lumiList.push_back(mapLumi.first);
    }
    std::sort(std::begin(lumiList), std::end(lumiList));

    if(lumiList.size() != 0) std::cout << " Run " << mapRun.first << ", lumi. section range " << lumiList[0] << "-" << lumiList[lumiList.size()-1] << std::endl;
    else std::cout << " Run " << mapRun.first << ", no lumi. sections!" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Event sel. stats: " << std::endl;
  for(auto const & perCut : totalEventsPerEventCut){
    std::cout << " " << perCut.first << ": " << perCut.second << std::endl;
  }
  
  std::cout << std::endl;
  std::cout << "D0 sel. stats: " << std::endl;
  for(auto const & perCut : totalDPerD0Cut){
    std::cout << " " << perCut.first << ": " << perCut.second << std::endl;
  }

  outFile_p->cd();

  outHist_p->Write("", TObject::kOverwrite);
  delete outHist_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}
