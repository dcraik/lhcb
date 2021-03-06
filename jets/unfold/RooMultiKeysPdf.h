/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOMULTIKEYSPDF
#define ROOMULTIKEYSPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
 
class RooMultiKeysPdf : public RooAbsPdf {
public:
  RooMultiKeysPdf() {} ; 
  RooMultiKeysPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _cat,
	      RooArgList& _pdfs);
  RooMultiKeysPdf(const RooMultiKeysPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMultiKeysPdf(*this,newname); }
  inline virtual ~RooMultiKeysPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy cat ;
  RooListProxy pdfs ;

  mutable int idx;
  
  Double_t evaluate() const ;

private:

//  ClassDef(RooMultiKeysPdf,1) // Your description goes here...
};
 
#endif
