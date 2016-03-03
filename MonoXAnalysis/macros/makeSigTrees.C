#include <vector>

#include "filters.C"

void makeSigTrees(string inputDIR, string outputDIR, bool isDIR = false, bool addBtag = true, bool storeGenTree = false) {

  vector<TFile*> fileList;
  vector<TTree*> treeList_gen;
  vector<TTree*> treeList_evt;
  vector<string> fileName; 

  if(isDIR){
    
    system(("ls "+inputDIR+" | grep root > file.tmp").c_str());
    
    ifstream infile("file.tmp");
    string line;
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);
	if(line == "") continue;
	fileName.push_back(line); 
	fileList.push_back(TFile::Open((inputDIR+"/"+line).c_str()));
	if(fileList.back()){
	  treeList_gen.push_back((TTree*)fileList.back()->Get("gentree/gentree"));
	  treeList_evt.push_back((TTree*)fileList.back()->Get("tree/tree"));
	}
	else{
	  fileList.pop_back();	  
	}
      }
    }
    
    infile.close();
    system("rm file.tmp");    
    system(("mkdir -p "+inputDIR+"/"+outputDIR).c_str());
  }
  else{ // dir is the file name
    fileName.push_back(inputDIR);
    fileList.push_back(TFile::Open(inputDIR.c_str()));
    if(fileList.back()){
      treeList_gen.push_back((TTree*)fileList.back()->Get("gentree/gentree"));
      treeList_evt.push_back((TTree*)fileList.back()->Get("tree/tree"));
    }
    else
      fileList.pop_back();	  
  }

  std::cout<<"####################################"<<std::endl;
  std::cout<<"make signal trees at posteriori    "<<std::endl;
  std::cout<<"####################################"<<std::endl;

  if(treeList_gen.size() != treeList_evt.size()) {
    cerr<<"exit --> different size between event tree list and gen tree list "<<endl;
    exit (EXIT_FAILURE);
  }


  for(size_t itree = 0; itree < treeList_gen.size(); itree++){
    
    TTreeReader treeReader(treeList_gen.at(itree));
    TTreeReaderValue<double> wgtsign (treeReader,"wgtsign");

    double wgtsum = 0.;
    float  nEvents   = 0.;
    double xsec      = 0.;
    while(treeReader.Next()){
      wgtsum += (*wgtsign);
      nEvents++;
    } 
    std::cout<<"Name "<<fileList.at(itree)->GetName()<<" XS "<<wgtsum/nEvents<<std::endl;
    xsec = (wgtsum/nEvents)*1000; // in femptobarn 

    TH1D* puRatio = pileupwgt(treeList_gen.at(itree));

    const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && (hltmet90 > 0 || hltmet120 > 0 || hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0)"; 

    TFile* outfile = NULL;
    if(isDIR)
      outfile = new TFile((inputDIR+"/"+outputDIR+"/sig_"+fileName.at(itree)).c_str(),"RECREATE");
    else
      outfile = new TFile(("sig_"+fileName.at(itree)).c_str(),"RECREATE");

    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree","tree");
    treedir->cd();
    // copy the tree applying a selection  
    treeList_evt.at(itree)->SetBranchStatus("xsec",0);
    TTree* outtree = treeList_evt.at(itree)->CopyTree(cut);
    // add a weight sum branch                                                                                                                                              
    TBranch* bwgtsum, *bwgtpileup, *bxsec;
    double wgtpileup = 1;

    bwgtsum    = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");
    bxsec      = outtree->Branch("xsec", &xsec, "xsec/D");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"sigfilter --> apply sumwgt and puweight"<<std::endl;
    while(myReader.Next()){
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();
      bxsec->Fill();
    }

    if(addBtag){
      TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
      TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
      eff_Num_b = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b");
      eff_Num_c = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c");
      eff_Num_ucsdg = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg");
      // take denominators                                                                                                                                                
      eff_Denom_b = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b");
      eff_Denom_c = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c");
      eff_Denom_ucsdg = (TH2F*) fileList.at(itree)->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg");

      // compute efficiency                                                                                                                                                   
      TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
      TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
      TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
      eff_b->Divide(eff_Denom_b);
      eff_c->Divide(eff_Denom_c);
      eff_ucsdg->Divide(eff_Denom_ucsdg);

      btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
    }

    // write tree in the file                                                                                                                                              
    outfile->cd();
    if(storeGenTree)
      treeList_gen.at(itree)->CloneTree()->Write();
    treedir->cd();
    outtree->Write();
    outfile->Close();
    
  }

  return;

}
