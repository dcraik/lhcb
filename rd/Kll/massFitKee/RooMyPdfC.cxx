/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooMyPdfC.cc,v 1.19 2005/06/20 15:51:06 wverkerke Exp $
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
#include "RooFit.h"

#include <math.h>
#include <math.h>

#include "RooMyPdfC.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

ClassImp(RooMyPdfC)

    RooMyPdfC::RooMyPdfC(const char *name, const char *title,
            RooAbsReal& mm, RooAbsReal& mc, RooAbsReal& alpha, 
	    //RooAbsReal& beta) : 
	    RooAbsReal& mu, RooAbsReal& sigma) :
        RooAbsPdf(name,title),
        _mm(    "mm",    "", this, mm),
        _mc(    "mc",    "", this, mc),
        _alpha( "alpha", "", this, alpha),
//        _beta(  "beta",  "", this, beta)
        _mu(    "mu",    "", this, mu),
        _sigma( "sigma", "", this, sigma)
{
}


RooMyPdfC::RooMyPdfC(const RooMyPdfC& other, const char* name):
    RooAbsPdf(other,name), _mm("mm",this,other._mm), _mc("mc",this,other._mc), _alpha("alpha",this,other._alpha), 
	//_beta("beta",this,other._beta)
	_mu("mu",this,other._mu), _sigma("sigma",this,other._sigma) 
{
}

RooMyPdfC::~RooMyPdfC() {
}


Double_t RooMyPdfC::evaluate() const {
    if( _mc < _mm ) return 0.;

//    Double_t edgeTerm = TMath::Log(_mc/_mm);
    Double_t edgeTerm = TMath::ATan(_mc-_mm);

    Double_t mmPart = TMath::Exp(-_alpha*_mm);
//    Double_t mcPart = TMath::Exp(-_beta*_mc);
    Double_t mcPart = TMath::Gaus(_mc,_mu,_sigma);

    return edgeTerm * mmPart * mcPart;
}

