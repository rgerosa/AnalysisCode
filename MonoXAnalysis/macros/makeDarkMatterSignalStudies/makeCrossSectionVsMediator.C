#include "../CMS_lumi.h"

static float minBosonPt = 200;
static float brWtoQuarks = 0.66;
static float brZtoQuarks = 0.69;

class fileXSEC {

public:
  fileXSEC(){
    mediatorM_ = 0.;
    xsec_ = 0.;
    efficiency_ = 0.;
  };
  fileXSEC(double & mediatorM, double & xsec, double & efficiency):
    mediatorM_(mediatorM),
    xsec_(xsec),
    efficiency_(efficiency){}

  bool operator < (fileXSEC & a) const {
    if(mediatorM_ <= a.mediatorM_) return true;
    else return false;
  }

  double mediatorM_;
  double xsec_;
  double efficiency_;

};

void makeCrossSectionVsMediator(string inputDIR_monojet,
				string inputDIR_monoW,
				string inputDIR_monoZ,
				string outputDIR,
				string postfix,
				int    dmMass,
				bool   isDMF = true){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  vector<TTree*> chain_monojet;
  vector<TTree*> chain_monoW;
  vector<TTree*> chain_monoZ;
  vector<TFile*> file_monojet;
  vector<TFile*> file_monoW;
  vector<TFile*> file_monoZ;

  
  if(inputDIR_monojet != ""){
    system(("ls "+inputDIR_monojet+" | grep Mchi-"+to_string(dmMass)+"_gSM | grep root > list.temp ").c_str());
    ifstream file ("list.temp");
    if(file.is_open()){
      string line;
      while(!file.eof()){
	getline(file,line);
	if(line == "") continue;
	file_monojet.push_back(TFile::Open((inputDIR_monojet+"/"+line).c_str()));
	if(file_monojet.back())
	   chain_monojet.push_back((TTree*) file_monojet.back()->Get("gentree/tree"));
      }
    }
  }

  if(inputDIR_monoW != ""){
    system(("ls "+inputDIR_monoW+" | grep Mchi-"+to_string(dmMass)+"_gSM | grep root > list.temp ").c_str());
    ifstream file ("list.temp");
    if(file.is_open()){
      string line;
      while(!file.eof()){
	getline(file,line);
	if(line == "") continue;	
	file_monoW.push_back(TFile::Open((inputDIR_monoW+"/"+line).c_str()));
	if(file_monoW.back())
	  chain_monoW.push_back((TTree*) file_monoW.back()->Get("gentree/tree"));
      }
    }
  }

  if(inputDIR_monoZ != ""){
    system(("ls "+inputDIR_monoZ+" | grep Mchi-"+to_string(dmMass)+"_gSM | grep root > list.temp ").c_str());
    ifstream file ("list.temp");
    if(file.is_open()){
      string line;
      while(!file.eof()){
	getline(file,line);
	if(line == "") continue;
	file_monoZ.push_back(TFile::Open((inputDIR_monoZ+"/"+line).c_str()));
	if(file_monoZ.back())
	  chain_monoZ.push_back((TTree*) file_monoZ.back()->Get("gentree/tree"));
      }
    }
  }

  system("rm list.temp");

  cout<<"N monojet trees "<<file_monojet.size()<<endl;
  cout<<"N monoW   trees "<<file_monoW.size()<<endl;
  cout<<"N monoZ   trees "<<file_monoZ.size()<<endl;


  // calculate sum of weights
  vector<double> sumwgt_monojet;
  vector<double> sumwgt_monoW;
  vector<double> sumwgt_monoZ;

  cout<<"Sum of weights monojet "<<endl;
  for(auto tree : chain_monojet){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");  
    TTreeReaderValue<float> xsec (reader,"xsec");  
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");  
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");  
    TTreeReaderValue<float> dmpt (reader,"dmpt");  
    
    sumwgt_monojet.push_back(0);
    while(reader.Next())
      sumwgt_monojet.back() += *wgt;
  }

  // calculate sum of weights
  cout<<"Sum of weights monoW "<<endl;
  for(auto tree : chain_monoW){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");  
    TTreeReaderValue<float> xsec (reader,"xsec");  
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");  
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");  
    TTreeReaderValue<float> dmpt (reader,"dmpt");  
    
    sumwgt_monoW.push_back(0);
    while(reader.Next())
      sumwgt_monoW.back() += *wgt;
  }

  // calculate sum of weights
  cout<<"Sum of weights monoZ "<<endl;
  for(auto tree : chain_monoZ){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");  
    TTreeReaderValue<float> xsec (reader,"xsec");  
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");  
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");  
    TTreeReaderValue<float> dmpt (reader,"dmpt");  
    
    sumwgt_monoZ.push_back(0);
    while(reader.Next())
      sumwgt_monoZ.back() += *wgt;
  }
  
  
  /// xsec * efficiency
  vector<fileXSEC> xsec_monojet;
  vector<fileXSEC> xsec_monoW;
  vector<fileXSEC> xsec_monoZ;

  // new loop 
  int ifile = 0;
  cout<<"Loop on monojet chain "<<endl;
  for(auto tree : chain_monojet){

    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");

    fileXSEC temp;
    xsec_monojet.push_back(temp);    
    if(isDMF)
      xsec_monojet.back().xsec_ = sumwgt_monojet.at(ifile-1)/double(tree->GetEntries());
    
    int iEvent = 0;
    while(reader.Next()){
      iEvent++;
      if(not isDMF and iEvent == 1) xsec_monojet.back().xsec_ = *xsec/1000.;
      if(iEvent == 1){
	xsec_monojet.back().mediatorM_ = *samplemedM;      
	if(xsec_monojet.back().mediatorM_ < 1 or xsec_monojet.back().mediatorM_ > 14000){
	  TString name_tmp (file_monojet.at(ifile-1)->GetName());
	  name_tmp.ReplaceAll(inputDIR_monojet.c_str(),"");
	  vector<string> seglist;
	  stringstream name(name_tmp.Data());
	  string segment;	  
	  while(getline(name, segment, '-')){
	    seglist.push_back(segment);
	  }
	  TString name_tmp_2 (seglist.at(1).c_str());
	  name_tmp_2.ReplaceAll("_Mchi","");
	  xsec_monojet.back().mediatorM_ = stod(string(name_tmp_2));
	}
      }
      if(*dmpt < minBosonPt) continue;
      //event passing the selection
      xsec_monojet.back().efficiency_ += *wgt;
    }
  }
  
  ///// --- 
  for(size_t imonojet = 0; imonojet < xsec_monojet.size(); imonojet++)
    xsec_monojet.at(imonojet).efficiency_ = xsec_monojet.at(imonojet).efficiency_/sumwgt_monojet.at(imonojet);

  // new loop 
  ifile = 0;
  cout<<"Loop on monoW chain "<<endl;
  for(auto tree : chain_monoW){

    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");

    fileXSEC temp;
    xsec_monoW.push_back(temp);    
    int iEvent = 0;
    
    while(reader.Next()){
      iEvent++;
      if(iEvent == 1) xsec_monoW.back().xsec_ = *xsec*brWtoQuarks/1000.;
      if(iEvent == 1){
	xsec_monoW.back().mediatorM_ = *samplemedM;      	  
      }
      if(*dmpt < minBosonPt) continue;
      //event passing the selection
      xsec_monoW.back().efficiency_ += *wgt;
    }
  }
  
  ///// --- 
  for(size_t imonoW = 0; imonoW < xsec_monoW.size(); imonoW++)
    xsec_monoW.at(imonoW).efficiency_ = xsec_monoW.at(imonoW).efficiency_/sumwgt_monoW.at(imonoW);

  // new loop 
  ifile = 0;
  cout<<"Loop on monoZ chain "<<endl;
  for(auto tree : chain_monoZ){

    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");

    fileXSEC temp;
    xsec_monoZ.push_back(temp);        
    int iEvent = 0;

    while(reader.Next()){
      iEvent++;
      if(iEvent == 1) xsec_monoZ.back().xsec_ = *xsec*brZtoQuarks/1000.;
      if(iEvent == 1) xsec_monoZ.back().mediatorM_ = *samplemedM;      
      if(*dmpt < minBosonPt) continue;
      //event passing the selection
      xsec_monoZ.back().efficiency_ += *wgt;
    }
  }
  
  ///// --- 
  for(size_t imonoZ = 0; imonoZ < xsec_monoZ.size(); imonoZ++)
    xsec_monoZ.at(imonoZ).efficiency_ = xsec_monoZ.at(imonoZ).efficiency_/sumwgt_monoZ.at(imonoZ);


  std::sort(xsec_monojet.begin(),xsec_monojet.end());
  std::sort(xsec_monoW.begin(),xsec_monoW.end());
  std::sort(xsec_monoZ.begin(),xsec_monoZ.end());

  /// graphs
  TCanvas* canvas = new TCanvas("canvas","",625,600);
  canvas->cd();

  TGraph* monojet_xsec = new TGraph();
  TGraph* monoW_xsec = new TGraph();
  TGraph* monoZ_xsec = new TGraph();
  
  int ipoint = 0;
  for(auto xsec : xsec_monojet){
    monojet_xsec->SetPoint(ipoint,xsec.mediatorM_,xsec.xsec_*xsec.efficiency_);
    cout<<"ipoint "<<ipoint<<" med "<<xsec.mediatorM_<<" xsec "<<xsec.xsec_*xsec.efficiency_<<endl;
    ipoint++;
  }
  ipoint = 0;

  if(xsec_monoW.size() != 0){
    for(auto xsec : xsec_monoW){
      monoW_xsec->SetPoint(ipoint,xsec.mediatorM_,xsec.xsec_*xsec.efficiency_);
      cout<<"ipoint "<<ipoint<<" med "<<xsec.mediatorM_<<" xsec "<<xsec.xsec_*xsec.efficiency_<<endl;
      ipoint++;
    }
  }
  
  if(xsec_monoZ.size() != 0){ 
    ipoint = 0;
    for(auto xsec : xsec_monoZ){
      monoZ_xsec->SetPoint(ipoint,xsec.mediatorM_,xsec.xsec_*xsec.efficiency_);
      cout<<"ipoint "<<ipoint<<" med "<<xsec.mediatorM_<<" xsec "<<xsec.xsec_*xsec.efficiency_<<endl;
      ipoint++;
    }
  }

  double maximum = xsec_monojet.back().mediatorM_;
  if(xsec_monoW.size()!=0)
    maximum = max(maximum,xsec_monoW.back().mediatorM_);
  if(xsec_monoZ.size()!=0)
    maximum = max(maximum,xsec_monoZ.back().mediatorM_);

  TH1F* frame = new TH1F("frame","frame",1,0,maximum);
  frame->GetXaxis()->SetTitle("Mediator mass [GeV]");
  frame->GetYaxis()->SetTitle("#sigma (pb)");
  frame->SetFillColor(0);
  frame->SetFillStyle(0);

  if(xsec_monoW.size() != 0 and xsec_monoZ.size() != 0)
    frame->GetYaxis()->SetRangeUser(
				    min(TMath::MinElement(monojet_xsec->GetN(),monojet_xsec->GetY()),min(TMath::MinElement(monoW_xsec->GetN(),monoW_xsec->GetY()),TMath::MinElement(monoZ_xsec->GetN(),monoZ_xsec->GetY())))*0.05,
				    max(TMath::MaxElement(monojet_xsec->GetN(),monojet_xsec->GetY()),max(TMath::MaxElement(monoW_xsec->GetN(),monoW_xsec->GetY()),TMath::MaxElement(monoZ_xsec->GetN(),monoZ_xsec->GetY())))*100);
  else if(xsec_monoW.size() != 0 and xsec_monoZ.size() == 0)
    frame->GetYaxis()->SetRangeUser(
				    min(TMath::MinElement(monojet_xsec->GetN(),monojet_xsec->GetY()),TMath::MinElement(monoW_xsec->GetN(),monoW_xsec->GetY()))*0.05,
				    max(TMath::MaxElement(monojet_xsec->GetN(),monojet_xsec->GetY()),TMath::MaxElement(monoW_xsec->GetN(),monoW_xsec->GetY()))*100);
  else if(xsec_monoW.size() == 0 and xsec_monoZ.size() != 0)
    frame->GetYaxis()->SetRangeUser(
				    min(TMath::MinElement(monojet_xsec->GetN(),monojet_xsec->GetY()),TMath::MinElement(monoZ_xsec->GetN(),monoZ_xsec->GetY()))*0.05,
				    max(TMath::MaxElement(monojet_xsec->GetN(),monojet_xsec->GetY()),TMath::MaxElement(monoZ_xsec->GetN(),monoZ_xsec->GetY()))*100);
  else if(xsec_monoW.size() == 0 and xsec_monoZ.size() == 0)
    frame->GetYaxis()->SetRangeUser(TMath::MinElement(monojet_xsec->GetN(),monojet_xsec->GetY())*0.05,TMath::MaxElement(monojet_xsec->GetN(),monojet_xsec->GetY())*100);
  
  frame->Draw();
  CMS_lumi(canvas,"");

  monojet_xsec->SetLineColor(kBlack);
  monojet_xsec->SetLineWidth(2);
  monojet_xsec->SetMarkerColor(kBlack);
  monojet_xsec->SetMarkerStyle(20);
  monojet_xsec->SetMarkerSize(1);
  monojet_xsec->Draw("CPsame");

  monoW_xsec->SetLineColor(kRed);
  monoW_xsec->SetLineWidth(2);
  monoW_xsec->SetMarkerColor(kRed);
  monoW_xsec->SetMarkerStyle(20);
  monoW_xsec->SetMarkerSize(1);
  if(xsec_monoW.size() != 0)
    monoW_xsec->Draw("CPsame");
  
  monoZ_xsec->SetLineColor(kBlue);
  monoZ_xsec->SetLineWidth(2);
  monoZ_xsec->SetMarkerColor(kBlue);
  monoZ_xsec->SetMarkerStyle(20);
  monoZ_xsec->SetMarkerSize(1);
  if(xsec_monoZ.size() != 0)
    monoZ_xsec->Draw("CPsame");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(monojet_xsec,"mono-jet","LP");
  if(xsec_monoW.size() != 0)
    leg->AddEntry(monoW_xsec,"mono-W","LP");
  if(xsec_monoZ.size() != 0)
    leg->AddEntry(monoZ_xsec,"mono-Z","LP");
  leg->AddEntry((TObject*)(0),"p_{T}^{med} > 200 GeV","");
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/xsec_vs_mediator_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/xsec_vs_mediator_"+postfix+".png").c_str(),"png");
  
  TFile* outputFile = new TFile((outputDIR+"/xsec_vs_mediator_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();
  monojet_xsec->Write("monojet_xsec");
  monoW_xsec->Write("monoW_xsec");
  monoZ_xsec->Write("monoZ_xsec");
  outputFile->Close();
  
}
