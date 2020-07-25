/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMultiKeysPdf.h" 
#include "RooAbsReal.h" 

#include "TMath.h"

//ClassImp(RooMultiKeysPdf) 

 RooMultiKeysPdf::RooMultiKeysPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _cat,
                        RooArgList& _pdfs) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   cat("cat","cat",this,_cat),
   pdfs("pdfs","pdfs",this),
   idx(-1)
 {
   pdfs.add(_pdfs);
 } 


 RooMultiKeysPdf::RooMultiKeysPdf(const RooMultiKeysPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   cat("cat",this,other.cat),
   pdfs("pdfs",this,other.pdfs),
   idx(other.idx)
 { 
 } 



 Double_t RooMultiKeysPdf::evaluate() const 
 { 
   if(idx<0 || idx>=pdfs.getSize()) {
     std::cout << "INFO in RooMultiKeysPdf::evaluate : Initialising index" << std::endl;
     idx = cat;

     if(cat!=idx || !cat.absArg()->getAttribute("Constant")) {
       std::cout << "WARNING in RooMultiKeysPdf::evaluate : Index must be a fixed integer" << std::endl;
       std::cout << "                                       Fixing to " << idx << std::endl;
     }
   }
   if(idx>=0 && idx<pdfs.getSize()) {
	   return static_cast<RooAbsPdf*>(pdfs.at(idx))->getValV();
   }

   return 0.0 ; 
 } 


