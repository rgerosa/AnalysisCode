Double_t ApproxErf(Double_t arg) {
  static const double erflim = 5.0;
  if( arg >  erflim) return 1.0;
  if( arg < -erflim) return -1.0;

  return TMath::Erf(arg);
}

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

Double_t Logistic(double *x, double*par){

  double xx = x[0];
  double p0 = par[0]; 
  double p1 = par[1]; 
  double p2 = par[2]; 
  double p3 = par[3]; 
  double p4 = par[4]; 

  double denom = TMath::Power((1.+p1*exp(-p2*(xx-p3))),1./p4);
  double val = p0*1./denom;
  return val;		      
  
} 


