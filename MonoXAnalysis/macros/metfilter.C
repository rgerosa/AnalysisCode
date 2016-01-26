#include "eventlist.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

void metfilter(std::string inputFile, std::string outputFile) {

  // take the txt files with the event list
  string CMSSW_BASE = getenv("CMSSW_BASE");
  EventList listCSC  = readEventList(CMSSW_BASE+"/src/AnalysisCode/MonoXAnalysis/data/eventListMET/csc2015.txt");
  EventList listEESC = readEventList(CMSSW_BASE+"/src/AnalysisCode/MonoXAnalysis/data/eventListMET/ecalscn1043093.txt");
  EventList listMuonBadTrack = readEventList(CMSSW_BASE+"/src/AnalysisCode/MonoXAnalysis/data/eventListMET/muonBadTrack.txt");
  EventList listBadTracks    = readEventList(CMSSW_BASE+"/src/AnalysisCode/MonoXAnalysis/data/eventListMET/badResolutionTrack.txt");

  TFile* infile = new TFile(inputFile.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  
  unsigned char flagcscnew  = 1;    
  unsigned char flageescnew = 1;    
  unsigned char flagmuontrack = 1;    
  unsigned char flagbadtrack  = 1;    
  
  const char* cut = "";
  
  TFile* outfile = new TFile(outputFile.c_str(), "RECREATE");
  outfile->cd();
  TTree* outtree = frtree->CopyTree(cut);
  
  TBranch* bcscnew  = outtree->Branch("flagcscnew" , &flagcscnew , "flagcscnew/b");
  TBranch* beescnew = outtree->Branch("flageescnew", &flageescnew, "flageescnew/b");
  TBranch* bmuontrack = outtree->Branch("flagmuontrack", &flagmuontrack, "flagmuontrack/b");
  TBranch* bbadtrack = outtree->Branch("flagbadtrack", &flagbadtrack, "flagbadtrack/b");
  
  TBranch  *brun   = outtree->GetBranch("run");
  TBranch  *blumi  = outtree->GetBranch("lumi");
  TBranch  *bevent = outtree->GetBranch("event");
  
  unsigned int run   = 0;
  unsigned int lumi  = 0;
  unsigned int event = 0;
  
  brun  ->SetAddress(&run);
  blumi  ->SetAddress(&lumi);
  bevent->SetAddress(&event);
  
  for (Long64_t i = 0; i < outtree->GetEntries(); i++) {
    brun  ->GetEvent(i);
    bevent->GetEvent(i);
    blumi->GetEvent(i);

    flagcscnew  = 1;        
    flagmuontrack = 1;    
    flagbadtrack  = 1;    
    flageescnew = 1;    
    
    // list csc
    auto rItr1(listCSC.find(run));
    if (rItr1 != listCSC.end()) {
      auto eItr(rItr1->second.find(event));
      if (eItr != (rItr1->second).end()) {
	flagcscnew = 0;
      }
    }
    
    // list EESC
    auto rItr2(listEESC.find(run));
    if (rItr2 != listEESC.end()) {
      auto eItr(rItr2->second.find(event));
      if (eItr != (rItr2->second).end()) {
	flageescnew = 0;
      }
    }
    
    // list badMuTrack
    auto rItr3(listMuonBadTrack.find(run));
    if (rItr3 != listMuonBadTrack.end()) {
      auto eItr(rItr3->second.find(event));
      if (eItr != (rItr3->second).end()) {
	flagmuontrack = 0;
      }
    }

    // list bad tracks
    auto rItr4(listBadTracks.find(run));
    if (rItr4 != listBadTracks.end()) {
      auto eItr(rItr4->second.find(event));
      if (eItr != (rItr4->second).end()) {
	flagbadtrack = 0;
      }
    }
    
    bcscnew ->Fill();
    beescnew->Fill();
    bmuontrack->Fill();
    bbadtrack->Fill();

  }
  
  outfile->Write();
   
}

