#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include <iostream>

int mmed(double mh, int code){
    if (code == 800) return ((int)(mh-80000000000))/10000; 
    if (code == 801) return ((int)(mh-80100000000))/10000; 
    if (code == 805) return ((int)(mh-80500000000))/10000; 
    if (code == 806) return ((int)(mh-80600000000))/10000; 
    return -1;
}

int mdm(double mh, int code){
    if (code == 800) return (mh-80000000000)  - ( ((Int_t)(mh-80000000000))/10000 )*10000;
    if (code == 801) return (mh-80100000000)  - ( ((Int_t)(mh-80100000000))/10000 )*10000;
    if (code == 805) return (mh-80500000000)  - ( ((Int_t)(mh-80500000000))/10000 )*10000;
    if (code == 806) return (mh-80600000000)  - ( ((Int_t)(mh-80600000000))/10000 )*10000;
    return -1;
}

int code(double mh){
    return (int)(mh/100000000);
}

double const mtop = 173.34;
double const mb   = 4.18;

bool STOP=false;
double WidthDM(double mMED, double mF){
    double fac = (1- 4*mF*mF/(mMED*mMED));
    double vev = 246.;
    return TMath::Sqrt(fac)*mMED/(8*TMath::Pi());
}
double Width(double mMED, double mF){
    double fac = (1- 4*mF*mF/(mMED*mMED));
    double vev = 246.;
    return TMath::Sqrt(fac)*mF*mF*mMED/(8*TMath::Pi()*vev*vev);
}

double vecF(double mMED,double mDM){

    //if !(mMED>0) return 10;
    if (mMED<=0) return 10;
    if (2*mDM > 0.99*mMED) {

        mMED = 2*mDM; //return 10;
        STOP=true;
    }
    double iGamma = 1;
    if (mMED < 2.*mtop) iGamma = 3*Width(mMED,mb) + WidthDM(mMED,mDM) ;
    else iGamma = 3*Width(mMED,mb) + 3*Width(mMED,mtop) + WidthDM(mMED,mDM) ;

    double  lVal    = 0.5*1./TMath::Pi()*3;
    double  lDenom  = (mMED*mMED-4.*mDM*mDM)*(mMED*mMED-4.*mDM*mDM);
    lDenom += mMED*mMED*iGamma*iGamma;
    lVal /= lDenom;
    lVal *= mMED *mMED;
    //lVal *= iCoupl*iCoupl;
    lVal *= TMath::Sqrt(1.-(4.2*4.2)/mDM/mDM);
    //Adding Yukawa Copulings
    lVal *= (4.2/246.)*(4.2/246.);
    lVal *= (mDM/246.)*(mDM/246.);
    
    double c = 1.167E-17;
//std::cout <<  iGamma << " but if MMED = 125, mDM = 1, " << WidthDM(125,1)++3*Width(125,mb) << " " ;
    return c*lVal;
}

void majoranatodirac(TGraph* gr){
    double X;
    double Y;
    double fac = 2;
    for (int p=0;p<gr->GetN();p++){
        gr->GetPoint(p,X,Y);
        gr->SetPoint(p,X,fac*Y);
    }
}

TGraph * makeOBV(TGraph *Graph1){

    STOP=false;
    TGraph *gr = new TGraph();
    double X;
    double Y;
    int pp=0;
    Graph1->GetPoint(0,X,Y);
    for (double MDM=0.1;MDM<=Y;MDM+=0.1){

        gr->SetPoint(pp,MDM,vecF(X,MDM));
        pp++;
    }
    for (int p =0;p<Graph1->GetN();p++){
        Graph1->GetPoint(p,X,Y);
        if (!(X >0)) continue;
        
        //std::cout << " mMED = " << X << " mDM =  " << Y << std::endl;
        //std::cout << " width/x-sec = "<< vecF(X,Y) << std::endl;
        if (STOP) continue;
        gr->SetPoint(pp,Y,vecF(X,Y));
        pp++;
    }
    gr->GetXaxis()->SetTitle("m_{DM}");
    //gr->GetYaxis()->SetTitle("Annihilation Cross Section (cm^{3}/s)");
    gr->GetYaxis()->SetTitle("< #sigma v > (cm^{3}/s)");
    gr->SetName(Form("%s_DD",Graph1->GetName()));
    gr->SetLineStyle(Graph1->GetLineStyle());
    gr->SetLineColor(Graph1->GetLineColor());
    gr->SetLineWidth(Graph1->GetLineWidth());

    return gr;
}


void plotPseudoScalarDD() {

    TString coupling = "025";
    TString energy = "13";

    TFile *file = new TFile("/afs/cern.ch/user/r/rgerosa/public/xZeynep/Scans8April/higgsCombine_COMB.PseudoScalar_DD.root");
    TTree *tree = (TTree*)file->Get("limit");

    TFile* file2 = new TFile("FermiLat.root");
    TGraph* grdd = (TGraph*)file2->Get("FermiLAT8Y");   
 
    TGraph2D* grexp = new TGraph2D();
    TGraph2D* grobs = new TGraph2D();

    double mh;
    double limit;
    float quantile;
    
    tree->SetBranchAddress("mh",&mh);
    tree->SetBranchAddress("limit",&limit);
    tree->SetBranchAddress("quantileExpected",&quantile);
    
    int expcounter = 0;
    int obscounter = 0;
    for (int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
    
        if (quantile != 0.5 && quantile != -1) continue;
        int c = code(mh);
        int medmass = mmed(mh, c);
        int dmmass = mdm(mh, c);
        
        //hack the plot
        if (medmass < 100) continue;

        if (quantile == 0.5) {
            expcounter++;
            grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
        }
        if (quantile == -1) {
            obscounter++;
            grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
        }
    }
    tree->ResetBranchAddresses();

    TH2D* hexp = new TH2D("hexp", "", 245, 100, 5000, 75, 0, 1500);
    TH2D* hobs = new TH2D("hobs", "", 245, 100, 5000, 75, 0, 1500);

    for (int i = 1; i <= 245; i++) {
        for (int j = 1; j <= 75; j++) {
            hexp->SetBinContent(i, j, grexp->Interpolate(double(i+5)*20., double(j)*20.));
            hobs->SetBinContent(i, j, grobs->Interpolate(double(i+5)*20., double(j)*20.));
        }
    }
    TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
    TH2* hobs2 = (TH2*)hobs->Clone("hobs2");

    for (int i = 1; i <= hexp2->GetNbinsX(); i++) {
        for (int j = 1; j <= hexp2->GetNbinsY(); j++) {
            if (hexp2->GetBinContent(i, j) <= 0) hexp2->SetBinContent(i, j, 100.);
        }
    }

    for (int i = 1; i <= hobs2->GetNbinsX(); i++) {
        for (int j = 1; j <= hobs2->GetNbinsY(); j++) {
            if (hobs2->GetBinContent(i, j) <= 0) hobs2->SetBinContent(i, j, 100.);
        }
    }

    hexp2->SetContour(2);
    hexp2->SetContourLevel(1, 1);
    hobs2->SetContour(2);
    hobs2->SetContourLevel(1, 1);

    hexp2->Draw("contz list");
    gPad->Update();

    TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
    std::vector<double> lXE;
    std::vector<double> lYE;
    int lTotalContsE = lContoursE->GetSize();
    for(int i0 = 0; i0 < lTotalContsE; i0++){
        TList * pContLevel = (TList*)lContoursE->At(i0);
        TGraph *pCurv = (TGraph*)pContLevel->First();
        for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
            for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
                lXE.push_back(pCurv->GetX()[i2]); 
                lYE.push_back(pCurv->GetY()[i2]);
            }
            pCurv->SetLineColor(kBlack);                                                                                                            
            pCurv = (TGraph*)pContLevel->After(pCurv);                                                                                        
        }
    }
    if(lXE.size() == 0) {
        lXE.push_back(0); 
        lYE.push_back(0); 
    }
    TGraph *lTotalE = new TGraph(lXE.size(),&lXE[0],&lYE[0]);

    lTotalE->SetLineColor(1);
    lTotalE->SetLineWidth(3);
    lTotalE->SetLineStyle(2);

    hobs2->Draw("contz list");
    gPad->Update();

    TObjArray *lContours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
    std::vector<double> lX;
    std::vector<double> lY;
    int lTotalConts = lContours->GetSize();
    for(int i0 = 0; i0 < lTotalConts; i0++){
        TList * pContLevel = (TList*)lContours->At(i0);
        TGraph *pCurv = (TGraph*)pContLevel->First();
        for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
            for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
                lX.push_back(pCurv->GetX()[i2]);
                lY.push_back(pCurv->GetY()[i2]);
            }
            pCurv->SetLineColor(kBlack);
            pCurv = (TGraph*)pContLevel->After(pCurv);
        }
    }
    if(lX.size() == 0) {
        lX.push_back(0); 
        lY.push_back(0); 
    }
    TGraph *lTotal = new TGraph(lX.size(),&lX[0],&lY[0]);

    lTotal->SetLineColor(1);
    lTotal->SetLineWidth(3);

    TGraph *DDE_graph = makeOBV(lTotalE);
    TGraph *DD_graph  = makeOBV(lTotal);
    grdd->SetLineColor(kGreen+1);
    grdd->SetLineWidth(3);



    //TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 800);    
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->SetLogx();
    canvas->SetLogy();

    TH1* frame = canvas->DrawFrame(6., 5e-31, 300., 5e-20, "");
    //frame->GetYaxis()->SetTitle("Annihilation Cross Section (cm^{3}/s)");
    frame->GetYaxis()->SetTitle("<#sigma v> (cm^{3}/s)");
    frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.55);
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->CenterTitle();

    frame->Draw();
    DDE_graph->Draw("L SAME");
    DD_graph->Draw("L SAME");
    DDE_graph->SetLineColor(2);
    DD_graph->SetLineColor(2);

    grdd->Draw("L SAME");


    //gr8e->Draw("L SAME");
    //gr8o->Draw("L SAME");


    
    TLegend *leg = new TLegend(0.22,0.64,0.75,0.84,NULL,"brNDC");
    //TLegend *leg = new TLegend(0.22,0.68,0.80,0.88,NULL,"brNDC");
    leg->AddEntry(DDE_graph,"CMS median expected 90% CL","L");
    leg->AddEntry(DD_graph ,"CMS obs. 90% CL","L");
    leg->AddEntry(grdd ,"FermiLAT","L");
    
    leg->SetFillColor(0);
    leg->Draw("SAME");
    
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.035);
    tex->Draw();
    tex->DrawLatex(0.73,0.96,"35.9 fb^{-1} (13 TeV)");
    
    //TLatex* texCMS = new TLatex(0.22,0.96,"#bf{CMS} #it{Preliminary}");
    TLatex* texCMS = new TLatex(0.22,0.90,"#bf{CMS}");
    texCMS->SetNDC();
    texCMS->SetTextFont(45);
    texCMS->SetLineWidth(2);
    texCMS->SetTextSize(0.042); texCMS->Draw();
    //tex->DrawLatex(0.24,0.9,"All limits at 90% CL");
    tex->DrawLatex(0.22,0.86,"#bf{Pseudoscalar med, Dirac DM,  g_{q} = 1, g_{DM} = 1}");


    canvas->SaveAs("/afs/cern.ch/user/z/zdemirag/www/Monojet/moriond_80x/unblinding_march26/test_monojetv2/fits/scanDD_ps_g"+coupling+"_"+energy+"TeV_v1.pdf");
    canvas->SaveAs("/afs/cern.ch/user/z/zdemirag/www/Monojet/moriond_80x/unblinding_march26/test_monojetv2/fits/scanDD_ps_g"+coupling+"_"+energy+"TeV_v1.png");
    canvas->SaveAs("/afs/cern.ch/user/z/zdemirag/www/Monojet/moriond_80x/unblinding_march26/test_monojetv2/fits/scanDD_ps_g"+coupling+"_"+energy+"TeV_v1.C");


}

