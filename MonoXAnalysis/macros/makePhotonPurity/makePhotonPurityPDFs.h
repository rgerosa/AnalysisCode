#ifndef PhotonPurityPDFs
#define PhotonPurityPDFs

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

////// Pow Pdf 
class RooPowPdf : public RooAbsPdf {
 public:
  RooPowPdf() {} ; 
  RooPowPdf(const char *name, const char *title,
	    RooAbsReal& _x,
	    RooAbsReal& _p0);

  RooPowPdf(const RooPowPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const { return new RooPowPdf(*this,newname); }

  inline virtual ~RooPowPdf() { }
  
 protected:
  
  RooRealProxy x ;
  RooRealProxy p0 ;
  
  Double_t evaluate() const ;
  
 private:
  
  ClassDef(RooPowPdf,1) // Your description goes here...
    };

#endif
