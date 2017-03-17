#include "../CMS_lumi.h"

static float minBosonPt = 200;
vector<float> metBin    = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};

vector<int> colors = {1,2,4,7,6,5,10};

void makeDMSimpDMFMediatorPtComparison(string inputDIR_DMSimp,
				       string inputDIR_DMF,
				       string outputDIR,
				       string postfix,
				       vector<int>  medMass,
				       vector<int>  dmMass,
				       bool useXSECFromWgtForDMF = true,
				       float  lumi = 36
				       ){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  vector<TTree*> chain_dmsimp;
  vector<TTree*> chain_dmf;
  vector<TFile*> file_dmsimp;
  vector<TFile*> file_dmf;

  if(medMass.size() != dmMass.size()){
      cerr<<"Problem with mass point vectors --> med mass does not have same size of dm one --> exit "<<endl;
      return;
  }

  
  if(inputDIR_DMSimp != ""){
    for(size_t ipoint = 0; ipoint < medMass.size(); ipoint++){
      system(("ls "+inputDIR_DMSimp+" | grep Mphi-"+to_string(medMass.at(ipoint))+"_ | grep Mchi-"+to_string(dmMass.at(ipoint))+"_gSM | grep root > list.temp ").c_str());
      ifstream file ("list.temp");
      if(file.is_open()){
	string line;
	while(!file.eof()){
	  getline(file,line);
	  if(line == "") continue;
	  file_dmsimp.push_back(TFile::Open((inputDIR_DMSimp+"/"+line).c_str()));
	  if(file_dmsimp.back())
	    chain_dmsimp.push_back((TTree*) file_dmsimp.back()->Get("gentree/tree"));
	}
      }
    }
  }

  if(inputDIR_DMF != ""){
    for(size_t ipoint = 0; ipoint < medMass.size(); ipoint++){
      system(("ls "+inputDIR_DMF+" | grep Mphi-"+to_string(medMass.at(ipoint))+"_ | grep Mchi-"+to_string(dmMass.at(ipoint))+"_gSM | grep root > list.temp ").c_str());
      ifstream file ("list.temp");
      if(file.is_open()){
	string line;
	while(!file.eof()){
	  getline(file,line);
	  if(line == "") continue;
	  file_dmf.push_back(TFile::Open((inputDIR_DMF+"/"+line).c_str()));
	  if(file_dmf.back())
	    chain_dmf.push_back((TTree*) file_dmf.back()->Get("gentree/tree"));
	}
      }
    }
  }
  
  system("rm list.temp");

  cout<<"N DMSimp trees "<<file_dmsimp.size()<<endl;
  cout<<"N DMF    trees "<<file_dmf.size()<<endl;

  // calculate sum of weights
  vector<double> sumwgt_dmsimp;
  vector<double> sumwgt_dmf;

  cout<<"Sum of weights DMSimp "<<endl;
  for(auto tree : chain_dmsimp){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");  
    TTreeReaderValue<float> xsec (reader,"xsec");  
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");  
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");  
    TTreeReaderValue<float> dmpt (reader,"dmpt");  
    
    sumwgt_dmsimp.push_back(0);
    while(reader.Next())
      sumwgt_dmsimp.back() += *wgt;
  }

  cout<<"Sum of weights DMF "<<endl;
  for(auto tree : chain_dmf){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");  
    TTreeReaderValue<float> xsec (reader,"xsec");  
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");  
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");  
    TTreeReaderValue<float> dmpt (reader,"dmpt");  
    
    sumwgt_dmf.push_back(0);
    while(reader.Next())
      sumwgt_dmf.back() += *wgt;
  }

  
  /// xsec * efficiency
  vector<TH1F*> mediatorPt_dmsimp;
  vector<TH1F*> mediatorPt_dmf;

  // new loop 
  int ifile = 0;
  cout<<"Loop on DMSimp chain "<<endl;
  for(auto tree : chain_dmsimp){

    mediatorPt_dmsimp.push_back(new TH1F(Form("mediatorPt_mmed_%d_mdm_%d_dmsimp",medMass.at(ifile),dmMass.at(ifile)),"",metBin.size()-1,&metBin[0]));
    mediatorPt_dmsimp.back()->Sumw2();

    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");

    int iEvent = 0;
    while(reader.Next()){
      iEvent++;
      if(*dmpt < minBosonPt) continue;
      //event passing the selection
      mediatorPt_dmsimp.back()->Fill(*dmpt,*wgt*lumi*(*xsec)/(sumwgt_dmsimp.at(ifile-1)));
    }
  }

  
  // new loop 
  ifile = 0;
  cout<<"Loop on DMF chain "<<endl;
  for(auto tree : chain_dmf){

    mediatorPt_dmf.push_back(new TH1F(Form("mediatorPt_mmed_%d_mdm_%d_dmf",medMass.at(ifile),dmMass.at(ifile)),"",metBin.size()-1,&metBin[0]));
    mediatorPt_dmf.back()->Sumw2();

    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");

    double XSEC = 1000*sumwgt_dmf.at(ifile-1)/tree->GetEntries();

    while(reader.Next()){
      if(not useXSECFromWgtForDMF) XSEC = *xsec;
      if(*dmpt < minBosonPt) continue;
      //event passing the selection
      mediatorPt_dmf.back()->Fill(*dmpt,*wgt*lumi*(XSEC)/(sumwgt_dmf.at(ifile-1)));
    }
  }

  /// graphs
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  double max = -1;
  double min = 99999999;

  for(int idmsimp = 0; idmsimp < mediatorPt_dmsimp.size(); idmsimp++){
    if(mediatorPt_dmsimp.at(idmsimp)->GetMinimum() < min)
      min = mediatorPt_dmsimp.at(idmsimp)->GetMinimum();
    if(mediatorPt_dmsimp.at(idmsimp)->GetMaximum() > max)
      max = mediatorPt_dmsimp.at(idmsimp)->GetMaximum();
  }

  for(int idmf = 0; idmf < mediatorPt_dmf.size(); idmf++){
    if(mediatorPt_dmf.at(idmf)->GetMinimum() < min)
      min = mediatorPt_dmf.at(idmf)->GetMinimum();
    if(mediatorPt_dmf.at(idmf)->GetMinimum() > max)
      max = mediatorPt_dmf.at(idmf)->GetMaximum();
  }

  TH1F* frame = new TH1F("frame","frame",metBin.size()-1,&metBin[0]);
  frame->GetXaxis()->SetTitle("Mediator p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("Events");
  frame->SetFillColor(0);
  frame->SetFillStyle(0);
  frame->GetYaxis()->SetRangeUser(min*0.1,max*100);
  frame->Draw();
  CMS_lumi(canvas,Form("%1.f",lumi));

  TLegend leg (0.6,0.5,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);

  int icolor = 0;
  for(auto histo : mediatorPt_dmsimp){
    histo->SetLineColor(colors.at(icolor));
    histo->SetLineWidth(2);
    histo->SetMarkerColor(colors.at(icolor));
    histo->SetMarkerStyle(20);
    histo->Draw("hist same");
    leg.AddEntry(histo,Form("DMSimp m_{med} = %d, m_{dm} = %d ",medMass.at(icolor),dmMass.at(icolor)),"L");    
    icolor++;
  }

  icolor = 0;
  for(auto histo : mediatorPt_dmf){
    histo->SetLineColor(colors.at(icolor));
    histo->SetLineWidth(2);
    histo->SetMarkerColor(colors.at(icolor));
    histo->SetMarkerStyle(20);
    histo->Draw("EP same");
    leg.AddEntry(histo,Form("DMF m_{med} = %d, m_{dm} = %d ",medMass.at(icolor),dmMass.at(icolor)),"EP");    
    icolor++;
  }

  leg.Draw("same");

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/spectrumComparison_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/spectrumComparison_"+postfix+".pdf").c_str(),"pdf");

  // plot with ratio
  TCanvas* canvas2 = new TCanvas("canvas2","",600,700);
  canvas2->SetBottomMargin(0.3);
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  canvas2->cd();
  TH1* frame2 =  (TH1*) frame->Clone("frame2");
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw();
  CMS_lumi(canvas2,Form("%1.f",lumi));
  for(auto histo : mediatorPt_dmsimp){
    histo->Draw("hist same");
  }
  for(auto histo : mediatorPt_dmf){
    histo->Draw("EP same");
  }
  leg.Draw("same");

  pad2->Draw();
  pad2->cd();  

  frame2->GetYaxis()->SetTitle("DMSimp / DMF");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.2);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->Draw();

  vector<TH1*> ratios;
  float min_ratio = 1000000;
  float max_ratio = 0;

  for(size_t ihist = 0; ihist < mediatorPt_dmsimp.size(); ihist++){
    ratios.push_back((TH1*) mediatorPt_dmsimp.at(ihist)->Clone(Form("ratio_%d",int(ihist))));
    ratios.back()->Divide(mediatorPt_dmf.at(ihist));
    for(int iBin = 0; iBin < mediatorPt_dmsimp.at(ihist)->GetNbinsX()-4; iBin++){
      if(ratios.back()->GetBinContent(iBin+1) < min_ratio) min_ratio = ratios.back()->GetBinContent(iBin+1);
      if(ratios.back()->GetBinContent(iBin+1) > max_ratio) max_ratio = ratios.back()->GetBinContent(iBin+1);
    }
    ratios.back()->Draw("hist same");
  }
  frame2->GetYaxis()->SetRangeUser(min_ratio*0.8,max_ratio*1.2);
  
  canvas2->SetLogy();
  canvas2->SaveAs((outputDIR+"/spectrumComparison_"+postfix+"_withRatio.png").c_str(),"png");
  canvas2->SaveAs((outputDIR+"/spectrumComparison_"+postfix+"_withRatio.pdf").c_str(),"pdf");

  
  
}
 
