/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooMyPdfD.cc,v 1.19 2005/06/20 15:51:06 wverkerke Exp $
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

#include "RooMyPdfD.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

ClassImp(RooMyPdfD)

    RooMyPdfD::RooMyPdfD(const char *name, const char *title,
            RooAbsReal& mm, RooAbsReal& mc, RooAbsReal& a0, RooAbsReal& a1, RooAbsReal& a2, RooAbsReal& b0, RooAbsReal& b1, RooAbsReal& b2) : //, RooAbsReal& mu, RooAbsReal& sigma) :
        RooAbsPdf(name,title),
        _mm( "mm", "", this, mm),
        _mc( "mc", "", this, mc),
        _a0( "a0", "", this, a0),
        _a1( "a1", "", this, a1),
        _a2( "a2", "", this, a2),
        _b0( "b0", "", this, b0),
        _b1( "b1", "", this, b1),
        _b2( "b2", "", this, b2)//,
//        _mu(    "mu",    "", this, mu),
//        _sigma( "sigma", "", this, sigma)
{
    RooRealVar mm1= (RooRealVar&)(_mm.arg());
    _lomm = mm1.getMin();
    _himm = mm1.getMax();
    _midmm = (_himm + _lomm)/2.;
    RooRealVar mc1= (RooRealVar&)(_mc.arg());
    _lomc = mc1.getMin();
    _himc = mc1.getMax();
    _midmc = (_himc + _lomc)/2.;
}


RooMyPdfD::RooMyPdfD(const RooMyPdfD& other, const char* name):
    RooAbsPdf(other,name), _mm("mm",this,other._mm), _mc("mc",this,other._mc), _a0("a0",this,other._a0), _a1("a1",this,other._a1), _a2("a2",this,other._a2), _b0("b0",this,other._b0), _b1("b1",this,other._b1), _b2("b2",this,other._b2), _lomm(other._lomm), _himm(other._himm), _midmm(other._midmm), _lomc(other._lomc), _himc(other._himc), _midmc(other._midmc)//, _mu("mu",this,other._mu), _sigma("sigma",this,other._sigma) 
{
}

RooMyPdfD::~RooMyPdfD() {
}


Double_t RooMyPdfD::evaluate() const {
    if( _mc < _mm ) return 0.;

    Double_t mm = 2.*(_mm-_midmm)/(_himm-_lomm);
    Double_t mc = 2.*(_mc-_midmc)/(_himc-_lomc);

    Double_t edgeTerm = TMath::Log(_mc/_mm);

    Double_t mmPart = 1. + _a0*mm + _a1*mm*mm + _a2*mm*mm*mm; 
    Double_t mcPart = 1. + _b0*mc + _b1*mc*mc + _b2*mc*mc*mc; 
//    Double_t mcPart = TMath::Gaus(_mc,_mu,_sigma);

//    std::cout << mm << "\t" << mc << "\t" << mmPart << "\t" << mcPart << "\t" << edgeTerm << std::endl;
    return edgeTerm * mmPart * mcPart;
}

