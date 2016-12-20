{

  TFile* f = new TFile("PhotonSFandEffandPurity_Lumi7p6fb_16072016.root", "RECREATE");
  double binx[8] = {175, 195,215,245,285,335,385,1005};
  double biny[3] = {0, 1.4442, 2.5};
  TH2F* h2 = new TH2F("purity", "purity", 7, binx,2,biny);
  h2->SetBinContent(h2->GetXaxis()->FindBin(180.), h2->GetYaxis()->FindBin(1.),0.9488);
  h2->SetBinContent(h2->GetXaxis()->FindBin(200.), h2->GetYaxis()->FindBin(1.),0.9561);
  h2->SetBinContent(h2->GetXaxis()->FindBin(220.), h2->GetYaxis()->FindBin(1.),0.9583);
  h2->SetBinContent(h2->GetXaxis()->FindBin(260.), h2->GetYaxis()->FindBin(1.),0.9641);
  h2->SetBinContent(h2->GetXaxis()->FindBin(300.), h2->GetYaxis()->FindBin(1.),0.9659);
  h2->SetBinContent(h2->GetXaxis()->FindBin(350.), h2->GetYaxis()->FindBin(1.),0.9568);
  h2->SetBinContent(h2->GetXaxis()->FindBin(600.), h2->GetYaxis()->FindBin(1.),0.9449);

  h2->SetBinError(h2->GetXaxis()->FindBin(200.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.00082,2)+pow(0.04,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(270.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.00074,2)+pow(0.026,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(370.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.00065,2)+pow(0.041,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(370.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.00058,2)+pow(0.041,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(370.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.00063,2)+pow(0.041,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(370.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.0011,2)+pow(0.041,2)));
  h2->SetBinError(h2->GetXaxis()->FindBin(370.), h2->GetYaxis()->FindBin(1.),sqrt(pow(0.0015,2)+pow(0.041,2)));
  f->cd();
  h2->Write();
  f->Write();
  f->Close();
}
