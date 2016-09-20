#include "RooErfExpPdf.h"


ClassImp(RooErfExpPdf) 

RooErfExpPdf::RooErfExpPdf(const char *name, const char *title, 
			   RooAbsReal& _x,
			   RooAbsReal& _c,
			   RooAbsReal& _offset,
			   RooAbsReal& _width) :
RooAbsPdf(name,title), 
  x("x","x",this,_x),
  c("c","c",this,_c),
  offset("offset","offset",this,_offset),
  width("width","width",this,_width){ } 


RooErfExpPdf::RooErfExpPdf(const RooErfExpPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  c("c",this,other.c),
  offset("offset",this,other.offset),
  width("width",this,other.width){}



Double_t RooErfExpPdf::evaluate() const { 

  Double_t width_tmp = width; 
  if(width<1e-2){ width_tmp = 1e-2;}
  double c_tmp = c;
  if(c==0) c_tmp = -1e-7;
  return TMath::Exp(c_tmp*x)*(1.+TMath::Erf((x-offset)/width_tmp))/2. ;
} 

Int_t RooErfExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  { 

  if (matchArgs(allVars,analVars,x)) return 1 ; 
  return 0 ; 
} 

Double_t RooErfExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const  { 

  Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
  if (code==1) { 
    Double_t minTerm=0;
    Double_t maxTerm=0;
    if(c==0){ 
      Double_t delta=-1e-7;
      minTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
		 TMath::Erf((2*x.min(rangeName)-delta*width_tmp*width_tmp-
			     2*offset)/2/width_tmp) - 
		 TMath::Exp(delta*x.min(rangeName)) * 
		 TMath::Erf((x.min(rangeName)-offset)/width_tmp) - 
		 TMath::Exp(delta*x.min(rangeName)))/-2/delta;
      maxTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
		 TMath::Erf((2*x.max(rangeName)-delta*width_tmp*width_tmp-
			     2*offset)/2/width_tmp) - 
		 TMath::Exp(delta*x.max(rangeName)) * 
		 TMath::Erf((x.max(rangeName)-offset)/width_tmp) - 
		 TMath::Exp(delta*x.max(rangeName)))/-2/delta;
        
    }else{
      minTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
		 TMath::Erf((2*x.min(rangeName)-c*width_tmp*width_tmp-
			     2*offset)/2/width_tmp) - 
		 TMath::Exp(c*x.min(rangeName)) * 
		 TMath::Erf((x.min(rangeName)-offset)/width_tmp) - 
		 TMath::Exp(c*x.min(rangeName)))/-2/c;
      maxTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
		 TMath::Erf((2*x.max(rangeName)-c*width_tmp*width_tmp-
			     2*offset)/2/width_tmp) - 
		 TMath::Exp(c*x.max(rangeName)) * 
		 TMath::Erf((x.max(rangeName)-offset)/width_tmp) - 
		 TMath::Exp(c*x.max(rangeName)))/-2/c;
    }
    return (maxTerm-minTerm) ;
  } 
  return 0 ; 
} 
