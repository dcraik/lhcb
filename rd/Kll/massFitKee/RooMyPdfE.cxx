/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooMyPdfE.cc,v 1.19 2005/06/20 15:51:06 wverkerke Exp $
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

#include "RooMyPdfE.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

ClassImp(RooMyPdfE)

    RooMyPdfE::RooMyPdfE(const char *name, const char *title,
            RooAbsReal& mm, RooAbsReal& mc, RooAbsReal& muM, RooAbsReal& sigmaM, RooAbsReal& muC, RooAbsReal& sigmaC) :
        RooAbsPdf(name,title),
        _mm(     "mm",     "", this, mm),
        _mc(     "mc",     "", this, mc),
        _muM(    "muM",    "", this, muM),
        _sigmaM( "sigmaM", "", this, sigmaM),
        _muC(    "muC",    "", this, muC),
        _sigmaC( "sigmaC", "", this, sigmaC)
{
}


RooMyPdfE::RooMyPdfE(const RooMyPdfE& other, const char* name):
    RooAbsPdf(other,name), _mm("mm",this,other._mm), _mc("mc",this,other._mc), _muM("muM",this,other._muM), _sigmaM("sigmaM",this,other._sigmaM), _muC("muC",this,other._muC), _sigmaC("sigmaC",this,other._sigmaC) 
{
}

RooMyPdfE::~RooMyPdfE() {
}


Double_t RooMyPdfE::evaluate() const {
    if( _mc < _mm ) return 0.;

    Double_t mmPart = TMath::Gaus(_mm,    _muM,_sigmaM);
    Double_t mcPart = TMath::Gaus(_mc-_mm,_muC,_sigmaC);

    return mmPart * mcPart;
}

