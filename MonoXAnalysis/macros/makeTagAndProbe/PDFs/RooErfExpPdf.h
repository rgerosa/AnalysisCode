#ifndef ROOERFEXPPDF_H
#define ROOERFEXPPDF_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "RooMath.h"

class RooErfExpPdf : public RooAbsPdf {
 public:
  RooErfExpPdf() {} ;  // default constructor
  RooErfExpPdf(const char *name, const char *title,
	       RooAbsReal& _x,
	       RooAbsReal& _c, // slope of the exp
	       RooAbsReal& _offset, // offset of the erf
	       RooAbsReal& _width); // width of the erf

  RooErfExpPdf(const RooErfExpPdf& other, const char* name=0) ; // ctor

  virtual TObject* clone(const char* newname) const { return new RooErfExpPdf(*this,newname); } // clone 

  inline virtual ~RooErfExpPdf() { } // dtor

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ; // analytic integral
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

 protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy offset ;
  RooRealProxy width ;
  
  Double_t evaluate() const ; // evaluate method 

 private:

  ClassDef(RooErfExpPdf,1); // Your description goes here...
};

 
#endif
