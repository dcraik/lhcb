/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOAPOLLONIOS
#define ROOAPOLLONIOS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooPromptShape : public RooAbsPdf {
public:
  RooPromptShape() {} ; 
  RooPromptShape(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mu,
	      RooAbsReal& _sigma,
	      RooAbsReal& _epsilon,
	      RooAbsReal& _rhoL,
	      RooAbsReal& _rhoR);
  RooPromptShape(const RooPromptShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPromptShape(*this,newname); }
  inline virtual ~RooPromptShape() { }

protected:

  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sigma ;
  RooRealProxy epsilon ;
  RooRealProxy rhoL ;
  RooRealProxy rhoR ;
  
  Double_t evaluate() const ;

private:

//  ClassDef(RooPromptShape,1) // Your description goes here...
};
 
#endif
