#include "makePhotonPurityPDFs.h"

#include "TMath.h"

ClassImp(RooPowPdf) 

RooPowPdf::RooPowPdf(const char *name, const char *title, 
		     RooAbsReal& _x,
		     RooAbsReal& _p0) :
RooAbsPdf(name,title), 
  x("x","x",this,_x),
  p0("p0","p0",this,_p0){} 


RooPowPdf::RooPowPdf(const RooPowPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  p0("p0",this,other.p0){} 



Double_t RooPowPdf::evaluate() const { 
  return TMath::Power(1+x, p0 )  ;   
} 
