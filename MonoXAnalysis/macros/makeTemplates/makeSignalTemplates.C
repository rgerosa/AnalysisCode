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
			 vector<string> observables            = {"met"}, // 1D histo
			 vector<string> observables_2D         = {},  // 2D histo
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
      findAllPossibleMassPoints(signalMassPoint,"Pseudoscalar",typeOfDMSignal);  
    }
  }

  string filename;
  if(interactionType == "")
    filename = outDir+"/templates_signal_"+templateSuffix+".root";
  else
    filename = outDir+"/templates_signal_"+interactionType+"_"+templateSuffix+".root";

  TFile outfile(filename.c_str(), "RECREATE");
  
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
    //signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"110",{5.579E+04,4.434E+03,1.335E+03,8.587E+02,1.309E+03,0.},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"125",{4.858E+04,3.782E+03,8.400E+02,5.328E+02,8.839E+02,1.227E+02},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"150",{3.129E+04,3.239E+03,5.037E+02,3.117E+02,5.279E+02},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"200",{1.694E+04,2.282E+03,1.899E+02,1.124E+02,2.054E+02},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"300",{6.590E+03,1.256E+03,4.348E+01,2.376E+01,4.132E+01},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"400",{3.160E+04,7.580E+02,1.432E+01,7.309E+00,1.273E+01},typeOfDMSignal);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"500",{1.709E+03,4.872E+02,5.825E+00,2.796E+00,5.256E+00},typeOfDMSignal);
  }
  
  outfile.Close();
  
}

