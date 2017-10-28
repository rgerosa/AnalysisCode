void makeYieldTableFromToysVBF(string referenceFileName, string directoryWithToys, string outputFileName){

  // take numbers for reference file
  TFile* referenceFile = TFile::Open(referenceFileName.c_str(),"READ");

  TH1F* zvvhist  = NULL;
  TH1F* wjetshist = NULL;
  TH1F* tophist  = NULL;
  TH1F* dibosonhist = NULL;
  TH1F* otherhist = NULL;
  TH1F* totalhist = NULL;

  zvvhist  = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/Znunu");
  zvvhist->Add((TH1F*)referenceFile->Get("shapes_fit_b/ch1/Znunu_EWK"));
  wjetshist = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/WJets");
  wjetshist->Add((TH1F*)referenceFile->Get("shapes_fit_b/ch1/WJets_EWK"));
  
  tophist = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/Top");

  dibosonhist = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/Dibosons");
  dibosonhist->Add((TH1F*)referenceFile->Get("shapes_fit_b/ch1/VGamma"));

  otherhist = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/QCD");
  otherhist->Add((TH1F*)referenceFile->Get("shapes_fit_b/ch1/ZJets"));

  totalhist = (TH1F*)referenceFile->Get("shapes_fit_b/ch1/total_background");

  TGraphAsymmErrors* data =  (TGraphAsymmErrors*)referenceFile->Get("shapes_fit_b/ch1/data");

  // Collect toys
  vector<TFile*> inputFile;
  vector<TH1F*> zvv_qcd;
  vector<TH1F*> zvv_ewk;
  vector<TH1F*> wjets_qcd;
  vector<TH1F*> wjets_ewk;

  system(("ls "+directoryWithToys+" | grep fitDiagnostics | grep root > file.list").c_str());
  ifstream infile;
  infile.open("file.list");
  string line;
  while(!infile.eof()){
    getline(infile,line);
    if(line == "" or line == "\n") continue;
    inputFile.push_back(TFile::Open((directoryWithToys+"/"+line).c_str(),"READ"));
    zvv_qcd.push_back((TH1F*) inputFile.back()->Get("shapes_fit_b/ch1/Znunu"));
    zvv_ewk.push_back((TH1F*) inputFile.back()->Get("shapes_fit_b/ch1/Znunu_EWK"));
    wjets_qcd.push_back((TH1F*) inputFile.back()->Get("shapes_fit_b/ch1/WJets"));
    wjets_ewk.push_back((TH1F*) inputFile.back()->Get("shapes_fit_b/ch1/WJets_EWK"));    
  }

  vector<double> zvv_variation;
  vector<double> wjets_variation;

  for(int iBin = 0; iBin < zvvhist->GetNbinsX(); iBin++){
    zvv_variation.push_back(0);
    wjets_variation.push_back(0);
  }
  
  for(size_t ihist = 0; ihist < zvv_qcd.size(); ihist++){ // loop on all histos
    for(int iBin = 0; iBin < zvvhist->GetNbinsX(); iBin++){ // loop on all bins
      // total Zvv in the toy - reference
      zvv_variation.at(iBin) += pow(((zvv_qcd.at(ihist)->GetBinContent(iBin+1)+zvv_ewk.at(ihist)->GetBinContent(iBin+1))-zvvhist->GetBinContent(iBin+1))*zvvhist->GetBinWidth(iBin+1),2);
      // total W+jets in the toy - reference
      wjets_variation.at(iBin) += pow(((wjets_qcd.at(ihist)->GetBinContent(iBin+1)+wjets_ewk.at(ihist)->GetBinContent(iBin+1))-wjetshist->GetBinContent(iBin+1))*wjetshist->GetBinWidth(iBin+1),2);
    }
  }			   
  infile.close();
  system("rm file.list");
  
  ofstream outputfile;
  outputfile.open(outputFileName.c_str());
  outputfile<<"$M_{jj}$ (GeV) & Observed & $Z \\rightarrow \\nu\\nu$+jets & $W \\rightarrow \\ell\\nu$+jets & Top & Dibosons & Other & Total Bkg. \\\\"<<endl;
  for(int ibin = 0; ibin < totalhist->GetNbinsX(); ibin++){
    double x,y;
    data->GetPoint(ibin,x,y);

    outputfile<<Form("%d-%d",int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)),int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2)))<<" & ";
    outputfile<<Form("%d",int(y*(int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2))-int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)))))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",zvvhist->GetBinContent(ibin+1)*zvvhist->GetBinWidth(ibin+1),sqrt(zvv_variation.at(ibin)/(zvv_qcd.size()+1)))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",wjetshist->GetBinContent(ibin+1)*wjetshist->GetBinWidth(ibin+1),sqrt(wjets_variation.at(ibin)/(wjets_qcd.size()+1)))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",tophist->GetBinContent(ibin+1)*tophist->GetBinWidth(ibin+1),tophist->GetBinError(ibin+1)*tophist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",dibosonhist->GetBinContent(ibin+1)*dibosonhist->GetBinWidth(ibin+1),dibosonhist->GetBinError(ibin+1)*dibosonhist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",otherhist->GetBinContent(ibin+1)*otherhist->GetBinWidth(ibin+1),otherhist->GetBinError(ibin+1)*otherhist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%.3f $\\pm$ %.3f",totalhist->GetBinContent(ibin+1)*totalhist->GetBinWidth(ibin+1),totalhist->GetBinError(ibin+1)*totalhist->GetBinWidth(ibin+1))<<" \\\\ ";
    outputfile<<"\n";
  }

  outputfile.close();
}
