void checkFitStatus(string inputDirectory, string nameToGrep){

  system(("ls "+inputDirectory+" | grep root | grep "+nameToGrep+" > list.temp ").c_str());
  vector<TFile*> inputFiles;

  ifstream infile ("list.temp");
  if(infile.is_open()){
    string line;
    while(!infile.eof()){
      getline(infile,line);
      if(TString(line).Contains(".root"))
	inputFiles.push_back(TFile::Open((inputDirectory+"/"+line).c_str(),"READ"));
    }
  }
  infile.close();
  system("rm list.temp");

  int nbadfit_bonly = 0;
  int invalid_bonly = 0;
  int nbadfit_sb = 0;
  int invalid_sb = 0;

  for(auto file: inputFiles){
    RooFitResult* fit_b = (RooFitResult*) file->Get("fit_b");
    RooFitResult* fit_s = (RooFitResult*) file->Get("fit_s");
    if(fit_b == 0){ invalid_bonly++; continue;}
    if(fit_s == 0){ invalid_sb++; continue;}

    RooRealVar* mu = (RooRealVar*)fit_s->floatParsFinal().find("r");
    cout<<"Input file : "<<file->GetName()<<" --> b only fit status "<<fit_b->status()<<" s+b fit status "<<fit_s->status()<<" mu value "<<mu->getVal()<<" pm "<<mu->getError()<<endl;
    if(fit_b->status() != 0)
      nbadfit_bonly++;
    if(fit_s->status() != 0 and fit_s->status() != 1)
      nbadfit_sb++;
  }
  
  cout<<"#######: total fit "<<inputFiles.size()<<" bad b-only "<<nbadfit_bonly<<" bad s+b "<<nbadfit_sb<<" invalid b-only "<<invalid_bonly<<" invalid sb "<<invalid_sb<<endl;
        
}
