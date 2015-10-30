/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooMyPdfC.rdl,v 1.9 2005/02/25 14:25:06 wverkerke Exp $
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

class RooMyPdfC : public RooAbsPdf {
    public:
        RooMyPdfC() {}
        RooMyPdfC(const char *name, const char *title,
                RooAbsReal& mm, RooAbsReal& mc, RooAbsReal& alpha, 
		//RooAbsReal& beta);
		RooAbsReal& mu, RooAbsReal& sigma);
        RooMyPdfC(const RooMyPdfC& other, const char* name=0);
        virtual TObject* clone(const char* newname) const {return new RooMyPdfC(*this,newname); }
        virtual ~RooMyPdfC();

    protected:

        RooRealProxy _mm ;
        RooRealProxy _mc ;
        RooRealProxy _alpha ;
//        RooRealProxy _beta ;
        RooRealProxy _mu ;
	RooRealProxy _sigma;
        Double_t evaluate() const;

    private:

        ClassDef(RooMyPdfC,2) // Non-Parametric KEYS PDF
};

#endif
