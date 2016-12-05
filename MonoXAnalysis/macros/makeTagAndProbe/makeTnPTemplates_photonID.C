{
    gSystem->Load("PDFs/RooErfExpPdf_cc.so");
    gROOT->ProcessLine(".L makeTnPTemplates.C");
    gROOT->ProcessLine("makeTnPTemplates(\"/home/rgerosa/MONOJET_ANALYSIS_2016_Data/TagAndProbe/TagAndProbe_MC/DYToEE_NNPDF30_13TeV-powheg-pythia8\",\"photon\",\"TemplatesNominal\",false)");    
}

