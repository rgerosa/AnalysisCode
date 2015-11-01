#ifndef MAKEHIST_H
#define MAKEHIST_H

void makehist(TTree* tree, TH1* hist, bool isMC, int sample, double scale, TH2* sflhist=NULL, TH2* sfthist=NULL, TH1* rhist=NULL, TH1* khist=NULL) {
    double lumi = 1.28;

    TFile* pufile = new TFile("purwt.root");
    TH1*   puhist = (TH1*)pufile->Get("puhist");

    hist->Sumw2();

    TBranch  *bnvtx       = tree->GetBranch("nvtx");
    TBranch  *bxsec       = tree->GetBranch("xsec");
    TBranch  *bwgt        = tree->GetBranch("wgt");
    TBranch  *bwgtsum     = tree->GetBranch("wgtsum"); 

    TBranch  *bhltm       = tree->GetBranch("hltmet90");
    TBranch  *bhlte       = tree->GetBranch("hltsingleel");
    TBranch  *bhltp       = tree->GetBranch("hltphoton165");

    TBranch  *bfhbhe      = tree->GetBranch("flaghbheloose");
    TBranch  *bfhbiso     = tree->GetBranch("flaghbheiso");
    TBranch  *bfcsc       = tree->GetBranch("flagcsctight");
    TBranch  *bfeeb       = tree->GetBranch("flageebadsc");

    TBranch  *bnjets      = tree->GetBranch("njets");
    TBranch  *bnbjets     = tree->GetBranch("nbjetslowpt");
    TBranch  *bjetpt      = tree->GetBranch("signaljetpt");
    TBranch  *bchfrac     = tree->GetBranch("signaljetCHfrac");
    TBranch  *bnhfrac     = tree->GetBranch("signaljetNHfrac");
    TBranch  *bemfrac     = tree->GetBranch("signaljetEMfrac");

    TBranch  *bmet        = tree->GetBranch("t1pfmet");
    TBranch  *bmmet       = tree->GetBranch("t1mumet");
    TBranch  *bmmetphi    = tree->GetBranch("t1mumetphi");
    TBranch  *bemet       = tree->GetBranch("t1elmet");
    TBranch  *bemetphi    = tree->GetBranch("t1elmetphi");
    TBranch  *bpmet       = tree->GetBranch("t1phmet");
    TBranch  *bpmetphi    = tree->GetBranch("t1phmetphi");
    TBranch  *bjmmdphi    = tree->GetBranch("incjetmetdphimin");
    TBranch  *bjemdphi    = tree->GetBranch("incjetelmetdphimin");
    TBranch  *bjpmdphi    = tree->GetBranch("incjetphmetdphimin");

    TBranch  *bmu1id      = tree->GetBranch("mu1id");
    TBranch  *bmu2id      = tree->GetBranch("mu2id");
    TBranch  *bmu1pt      = tree->GetBranch("mu1pt");
    TBranch  *bmu2pt      = tree->GetBranch("mu2pt");
    TBranch  *bmu1eta     = tree->GetBranch("mu1eta");
    TBranch  *bmu2eta     = tree->GetBranch("mu2eta");
    TBranch  *bel1id      = tree->GetBranch("el1id");
    TBranch  *bel2id      = tree->GetBranch("el2id");
    TBranch  *bel1pt      = tree->GetBranch("el1pt");
    TBranch  *bel2pt      = tree->GetBranch("el2pt");
    TBranch  *bel1eta     = tree->GetBranch("el1eta");
    TBranch  *bel2eta     = tree->GetBranch("el2eta");
    TBranch  *bphid       = tree->GetBranch("phidm");
    TBranch  *bphpt       = tree->GetBranch("phpt");
    TBranch  *bpheta      = tree->GetBranch("pheta");

    TBranch  *bwzpt       = tree->GetBranch("wzpt");
    TBranch  *bzmass      = tree->GetBranch("zmass");
    TBranch  *bzpt        = tree->GetBranch("zpt");

    UInt_t   nvtx         = 0;
    Double_t xsec         = 0.0;
    Double_t wgt          = 0.0;
    Double_t wgtsum       = 0.0;

    UChar_t  hltm         = 0;
    UChar_t  hlte         = 0;
    UChar_t  hltp         = 0;

    UChar_t  fhbhe        = 0;
    UChar_t  fhbiso       = 0;
    UChar_t  fcsc         = 0;
    UChar_t  feeb         = 0;

    UInt_t   njets        = 0;
    UInt_t   nbjets       = 0;
    Double_t jetpt        = 0.0;
    Double_t chfrac       = 0.0;
    Double_t nhfrac       = 0.0;
    Double_t emfrac       = 0.0;

    Double_t pfmet        = 0.0;
    Double_t mmet         = 0.0;
    Double_t mmetphi      = 0.0;
    Double_t emet         = 0.0;
    Double_t emetphi      = 0.0;
    Double_t pmet         = 0.0;
    Double_t pmetphi      = 0.0;
    Double_t jmmdphi      = 0.0;
    Double_t jemdphi      = 0.0;
    Double_t jpmdphi      = 0.0;

    Int_t    mu1id        = 0;
    Int_t    mu2id        = 0;
    Double_t mu1pt        = 0.0;
    Double_t mu2pt        = 0.0;
    Double_t mu1eta       = 0.0;
    Double_t mu2eta       = 0.0;
    Int_t    el1id        = 0;
    Int_t    el2id        = 0;
    Double_t el1pt        = 0.0;
    Double_t el2pt        = 0.0;
    Double_t el1eta       = 0.0;
    Double_t el2eta       = 0.0;
    Int_t    phid         = 0;
    Double_t phpt         = 0.0;
    Double_t pheta        = 0.0;

    Double_t wzpt         = 0.0;
    Double_t zmass        = 0.0;
    Double_t zpt          = 0.0;

    bnvtx                 ->SetAddress(&nvtx);
    bxsec                 ->SetAddress(&xsec);
    bwgt                  ->SetAddress(&wgt);
    bwgtsum               ->SetAddress(&wgtsum);

    bhltm                 ->SetAddress(&hltm);
    bhlte                 ->SetAddress(&hlte);
    bhltp                 ->SetAddress(&hltp);

    bfhbhe                ->SetAddress(&fhbhe); 
    bfhbiso               ->SetAddress(&fhbiso);
    bfcsc                 ->SetAddress(&fcsc);
    bfeeb                 ->SetAddress(&feeb);

    bnjets                ->SetAddress(&njets);
    bnbjets               ->SetAddress(&nbjets);
    bjetpt                ->SetAddress(&jetpt);
    bchfrac               ->SetAddress(&chfrac);
    bnhfrac               ->SetAddress(&nhfrac);
    bemfrac               ->SetAddress(&emfrac);

    bmet                  ->SetAddress(&pfmet);
    bmmet                 ->SetAddress(&mmet);
    bmmetphi              ->SetAddress(&mmetphi);
    bemet                 ->SetAddress(&emet);
    bemetphi              ->SetAddress(&emetphi);
    bpmet                 ->SetAddress(&pmet);
    bpmetphi              ->SetAddress(&pmetphi);
    bjmmdphi              ->SetAddress(&jmmdphi);
    bjemdphi              ->SetAddress(&jemdphi);
    bjpmdphi              ->SetAddress(&jpmdphi);

    bmu1id                ->SetAddress(&mu1id);
    bmu2id                ->SetAddress(&mu2id);
    bmu1pt                ->SetAddress(&mu1pt);
    bmu2pt                ->SetAddress(&mu2pt);
    bmu1eta               ->SetAddress(&mu1eta);
    bmu2eta               ->SetAddress(&mu2eta);
    bel1id                ->SetAddress(&el1id);
    bel2id                ->SetAddress(&el2id);
    bel1pt                ->SetAddress(&el1pt);
    bel2pt                ->SetAddress(&el2pt);
    bel1eta               ->SetAddress(&el1eta);
    bel2eta               ->SetAddress(&el2eta);
    bphid                 ->SetAddress(&phid);
    bphpt                 ->SetAddress(&phpt);
    bpheta                ->SetAddress(&pheta);

    bwzpt                 ->SetAddress(&wzpt);
    bzmass                ->SetAddress(&zmass);
    bzpt                  ->SetAddress(&zpt);

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        bnvtx             ->GetEvent(i);
        bxsec             ->GetEvent(i);
        bwgt              ->GetEvent(i);
        bwgtsum           ->GetEvent(i);

        bhltm             ->GetEvent(i);
        bhlte             ->GetEvent(i);
        bhltp             ->GetEvent(i);

        bfhbhe            ->GetEvent(i);
        bfhbiso           ->GetEvent(i);
        bfcsc             ->GetEvent(i);
        bfeeb             ->GetEvent(i);

        bnjets            ->GetEvent(i);
        bnbjets           ->GetEvent(i);
        bjetpt            ->GetEvent(i);
        bchfrac           ->GetEvent(i);
        bnhfrac           ->GetEvent(i);
        bemfrac           ->GetEvent(i);

        bmet              ->GetEvent(i);
        bmmet             ->GetEvent(i);
        bmmetphi          ->GetEvent(i);
        bemet             ->GetEvent(i);
        bemetphi          ->GetEvent(i);
        bpmet             ->GetEvent(i);
        bpmetphi          ->GetEvent(i);
        bjmmdphi          ->GetEvent(i);
        bjemdphi          ->GetEvent(i);
        bjpmdphi          ->GetEvent(i);

        bmu1id            ->GetEvent(i);
        bmu2id            ->GetEvent(i);
        bmu1pt            ->GetEvent(i);
        bmu2pt            ->GetEvent(i);
        bmu1eta           ->GetEvent(i);
        bmu2eta           ->GetEvent(i);

        bel1id            ->GetEvent(i);
        bel2id            ->GetEvent(i);
        bel1pt            ->GetEvent(i);
        bel2pt            ->GetEvent(i);
        bel1eta           ->GetEvent(i);
        bel2eta           ->GetEvent(i);

        bphid             ->GetEvent(i);
        bphpt             ->GetEvent(i);
        bpheta            ->GetEvent(i);

        bwzpt             ->GetEvent(i);
        bzmass            ->GetEvent(i);
        bzpt              ->GetEvent(i);
        
        Double_t hlt = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) hlt = hltm;
        else if (sample == 3 || sample == 4)                hlt = hlte;
        else if (sample == 5)                               hlt = hltp;

        Double_t jmdphi = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) jmdphi = fabs(jmmdphi);
        else if (sample == 3 || sample == 4)                jmdphi = fabs(jemdphi);
        else if (sample == 5)                               jmdphi = fabs(jpmdphi);

        Double_t met = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) met = mmet;
        else if (sample == 3 || sample == 4)                met = emet;
        else if (sample == 5)                               met = pmet;

        Double_t puwgt = 0.;
        if (nvtx <= 35) puwgt = puhist->GetBinContent(nvtx);

        Int_t    id1   = 0;
        Int_t    id2   = 0;
        Double_t pt1   = 0.0;
        Double_t pt2   = 0.0;
        Double_t eta1  = 0.0;
        Double_t eta2  = 0.0;

        if (sample == 1 || sample == 2) {
            id1  = mu1id;
            id2  = mu2id;
            pt1  = mu1pt;
            pt2  = mu2pt;
            eta1 = fabs(mu1eta);
            eta2 = fabs(mu2eta);
        }
        else if (sample == 3 || sample == 4) {
            id1  = el1id;
            id2  = el2id;
            pt1  = el1pt;
            pt2  = el2pt;
            eta1 = fabs(el1eta);
            eta2 = fabs(el2eta);
        }
        else if (sample == 5) {
            id1  = 1.0;
            id2  = 1.0;
            pt1  = phpt;
            eta1 = fabs(pheta);
        }

        Double_t sfwgt = 1.0;
        if (isMC && sflhist && sfthist) {
            if (pt1 > 0.) {
                if (id1 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt1, eta1)); 
                else          sfwgt *= sfthist->GetBinContent(sflhist->FindBin(pt1, eta1)); 
            }
            if (pt2 > 0.) {
                if (id2 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt2, eta2)); 
                else          sfwgt *= sfthist->GetBinContent(sflhist->FindBin(pt2, eta2)); 
            }
        }

        Double_t rwgt = 1.0;
        if (rhist) rwgt = rhist->GetBinContent(rhist->FindBin(met));

        Double_t kwgt = 1.0;
        if (isMC && khist) kwgt = khist->GetBinContent(khist->FindBin(wzpt));

        if (hlt  == 0) continue;
        if (fhbhe == 0 || fhbiso == 0 || fcsc == 0 || feeb == 0) continue;
        if (chfrac < 0.1) continue;
        if (nhfrac > 0.8) continue;
        if (njets  < 1) continue;
        if (nbjets > 0) continue;
        if (jetpt  < 100.) continue;
        if (jmdphi < 0.5) continue;
        if (sample == 5 && fabs(pheta) > 1.4442) continue;
        if (sample == 4 && pfmet < 50.) continue;
        if (met < 200.) continue;

        double fillvar = met;
        //if (fillvar >= hist->GetBinContent(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());

        double evtwgt = 1.0;
        if (isMC) evtwgt = xsec*scale*lumi*wgt*puwgt*sfwgt*rwgt*kwgt/wgtsum;

        hist->Fill(fillvar, evtwgt);
    }

}

#endif
