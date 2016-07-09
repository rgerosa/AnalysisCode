#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "makehist.h"
#include "makeCorrHistograms.C"
#include "makeDataHistograms.C"
#include "makeSignalHistograms.C"

using namespace std;


// Run the final analysis on signal samples
void makeSignalTemplates(
			 const Category & category     = Category::monojet,  // 0 = inclusive mono-j, 1 = exclsuive mono-j, 2 V-tag HP ..
			 const double & lumi                   = 2.30, // 
			 const string & outDir                 = "", // output dir for template file
			 const string & templateSuffix         = "",  // suffix for the output file
			 vector<string> observables    = {"met"}, // 1D histo
			 vector<string> observables_2D = {},  // 2D histo
			 const string & interactionType        = "",  // can be Vector, Axial, Scalar, Pseudoscalar or All
			 const bool & doShapeSystematics       = false, // run all the met, b-tag shape variations
			 const int  & typeOfDMSignal           = 0,     // 0 means both mono-j and mono-V, 1 is mono-j, 2 is mono-V
			 const bool & runHiggsInvisible        = false // run Higgs invisible analysis
			 ) {

  system(("mkdir -p "+outDir).c_str());

  initializeBinning();  

  // find all possible mass pont to use in the analysis for each Model: Vector, Axial, Scalar and Pseudoscalar .. if onlyMonoJetSignal is true just use all the available mono-j signals in the directories
  vector<signalSample> signalMassPoint;
  if(not runHiggsInvisible){
    if(interactionType != "All")
      findAllPossibleMassPoints(signalMassPoint,interactionType,typeOfDMSignal);  
    else{
      findAllPossibleMassPoints(signalMassPoint,"Vector",typeOfDMSignal);  
      findAllPossibleMassPoints(signalMassPoint,"Axial", typeOfDMSignal);  
      findAllPossibleMassPoints(signalMassPoint,"Scalar",typeOfDMSignal);  
      findAllPossibleMassPoints(signalMassPoint,"PseudoScalar",typeOfDMSignal);  
    }
  }

  TFile outfile((outDir+"/templates_signal_"+interactionType+"_"+templateSuffix+".root").c_str(), "RECREATE");
  
  // signal region templates
  cout<<"start signal region shapes for signal"<<endl;
  if(not runHiggsInvisible){
    if(interactionType != "All")
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,interactionType,lumi,doShapeSystematics);
    else{
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Vector",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Axial",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Scalar",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Pseudoscalar",lumi,doShapeSystematics);
    }
  }
  else if(runHiggsInvisible){
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"110",{5.507E+04,4.434E+03,2.194E+03,1.309E+03},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"125",{4.414E+04,3.782E+03,1.373E+03,7.612E+02,1.227E+02},1);
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"150",{3.210E+04,3.239E+03,8.154E+02,5.279E+02},1);
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"200",{1.812E+04,2.282E+03,3.023E+02,2.054E+02},1);
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"300",{9.823E+03,1.256E+03,6.724E+01,4.132E+01},1);
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"400",{9.516E+04,7.580E+02,2.163E+01,1.273E+01},1);
    //    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"500",{4.538E+03,4.872E+02,8.621E+00,5.256E+00},1);
  }
  
  outfile.Close();
  
}

