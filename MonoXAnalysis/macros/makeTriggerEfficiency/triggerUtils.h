// recoil binning                                                                                                                                                                                      
vector <float> bins_monojet_recoil   = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 225., 250., 275., 300., 325., 350., 400., 450., 500., 550., 650.,800., 1000., 1250};
vector <float> bins_vbf_recoil       = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};
vector <float> bins_monoV_recoil     = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};

// mjj                                                                                                                                                                                     
vector <float> bins_vbf_recoilvsmjj  = {0.,50.,75.,100.,125.,150.,175.,200., 225, 250., 300., 400., 550., 800., 1500};
vector <float> bins_vbf_mjj          = {0.,800.,1200.,1700.,3000.};

// detajj                                                                                                                                                                               
vector <float> bins_vbf_recoilvsdetajj = {0.,50.,75.,100.,125.,150.,175.,200., 225, 250., 300., 400., 550., 800., 1500};
vector <float> bins_vbf_detajj         = {0.,1.5,2.5,4.5,9};

// Eras for 2016 data                                                                                                                                                                    
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};


Double_t ApproxErf(Double_t arg) {

  static const double erflim = 5.0;
  if( arg >  erflim) return 1.0;
  if( arg < -erflim) return -1.0;
  
  return TMath::Erf(arg);
}

// Error function x CB
Double_t ErfCB(double *x, double *par) {

  double m = x[0];
  double m0 = par[0];
  double sigma = par[1];
  double alpha = par[2];
  double n = par[3];
  double norm = par[4];

  const double sqrtPiOver2 = 1.2533141373; // sqrt(pi/2)                                                                                                                       
  const double sqrt2 = 1.4142135624;

  Double_t sig = fabs((Double_t) sigma);
  Double_t t = (m - m0)/sig ;

  if (alpha < 0) t = -t;

  Double_t absAlpha = fabs(alpha / sig);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  Double_t b = absAlpha - n/absAlpha;

  Double_t leftArea = (1 + ApproxErf( absAlpha / sqrt2 )) * sqrtPiOver2 ;
  Double_t rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  Double_t area = leftArea + rightArea ;

  

  if ( t <= absAlpha ){
    double val = norm * (1 + ApproxErf( t / sqrt2 )) * sqrtPiOver2 / area ;
    if (val > 1) return 1;
    else return val;
  }
  else{
    double val = norm * (leftArea +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area ;
    if(val > 1 ) return 1;
    else return val;
  }
}

