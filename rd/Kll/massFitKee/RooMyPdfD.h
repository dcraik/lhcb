/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooMyPdfD.rdl,v 1.9 2005/02/25 14:25:06 wverkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   UC San Diego,        raven@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MY_KEYS
#define ROO_MY_KEYS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooMyPdfD : public RooAbsPdf {
    public:
        RooMyPdfD() {}
        RooMyPdfD(const char *name, const char *title,
                RooAbsReal& mm, RooAbsReal& mc, RooAbsReal& a0, RooAbsReal& a1, RooAbsReal& a2, RooAbsReal& b0, RooAbsReal& b1, RooAbsReal& b2);//, RooAbsReal& mu, RooAbsReal& sigma);
        RooMyPdfD(const RooMyPdfD& other, const char* name=0);
        virtual TObject* clone(const char* newname) const {return new RooMyPdfD(*this,newname); }
        virtual ~RooMyPdfD();

    protected:

        RooRealProxy _mm ;
        RooRealProxy _mc ;
        RooRealProxy _a0 ;
        RooRealProxy _a1 ;
        RooRealProxy _a2 ;
        RooRealProxy _b0 ;
        RooRealProxy _b1 ;
        RooRealProxy _b2 ;
//        RooRealProxy _mu ;
//	RooRealProxy _sigma;
        Double_t evaluate() const;

    private:

	Double_t _lomm, _himm, _midmm;
	Double_t _lomc, _himc, _midmc;

        ClassDef(RooMyPdfD,2) // Non-Parametric KEYS PDF
};

#endif
