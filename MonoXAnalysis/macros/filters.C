double sumwgt(TTree* tree) {
    TBranch *bweight = tree->GetBranch("wgtsign");

    double vweight  = 0.0;

    bweight->SetAddress(&vweight);

    double weightsum = 0.;
    for (int i = 0; i < tree->GetEntries(); i++) {
        bweight->GetEvent(i);
        weightsum += vweight;
    }

    return weightsum;
}

void sigfilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);

    const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && hltmet90 > 0";

    TFile* outfile = new TFile("sigtree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}


void zmmfilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);
    
    const char* cut = "nmuons == 2 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && zmass > 60 && zmass < 120 && mu1pt > 20 && (mu1id == 1 || mu2id == 1)";

    TFile* outfile = new TFile("zmmtree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}

void zeefilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);

    const char* cut = "nmuons == 0 && nelectrons == 2 && ntaus == 0 && nphotons == 0 && zeemass > 60 && zeemass < 120 && el1pt > 40 && (el1id == 1 || el2id == 1)";

    TFile* outfile = new TFile("zeetree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}

void wmnfilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);

    const char* cut = "nmuons == 1 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && mu1pt > 20 && mu1id == 1";

    TFile* outfile = new TFile("wmntree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}

void wenfilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);

    const char* cut = "nmuons == 0 && nelectrons == 1 && ntaus == 0 && nphotons == 0 && el1pt > 40 && el1id == 1";

    TFile* outfile = new TFile("wentree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}

void gamfilter() {
    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);

    const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 1 && phpt > 175 && phidm == 1";

    TFile* outfile = new TFile("gamtree.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    for (Long64_t i = 0; i < outtree->GetEntries(); i++)bwgtsum->Fill();

    outfile->Write();

}

