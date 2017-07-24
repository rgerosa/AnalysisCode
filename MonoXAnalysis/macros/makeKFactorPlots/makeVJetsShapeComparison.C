#include "../CMS_lumi.h"

enum class Sample {znn, zll, wjet};
enum class Category {monojet,VBF,VBFrelaxed};

// gen-level studies
void makeVJetsShapeComparison(string inputDIR_NLO, string inputDIR_LO, string inputKFactor, Sample sample, Category category, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  
  TChain* chain_nlo = new TChain("gentree/tree");
  TChain* chain_lo = new TChain("gentree/tree");
  chain_nlo->Add((inputDIR_NLO+"/*root").c_str());
  chain_lo->Add((inputDIR_LO+"/*root").c_str());


  
  makeSelection(chain_nlo,distribution_nlo,khists)

}
