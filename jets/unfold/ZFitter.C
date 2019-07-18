#include "ZFitter.h"

#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "RooAbsDataStore.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMsgService.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooPromptShape.h"

#include <boost/progress.hpp>

#include "cloneHists.h"
#include "outputFunctions.h"

void ZFitter::setInputs(TString data) {
	_dataFile  = data;

	_inputsSet=true;
}

bool ZFitter::fit(double& yieldZ, double& errorZ, uint binPT, uint binY) {
	if(binPT>0 && binPT<=static_cast<uint>(_ptBins->GetNbinsX())) {
		_jptmin = _ptBins->GetBinLowEdge(binPT);
		_jptmax = _ptBins->GetBinLowEdge(binPT+1);
	} else {
		_jptmin = _ptBins->GetBinLowEdge(1);
		_jptmax = _ptBins->GetBinLowEdge(_ptBins->GetNbinsX()+1);
	}
	if(binY>0 && binY<=static_cast<uint>(_yBins->GetNbinsX())) {
		_ymin = _yBins->GetBinLowEdge(binY);
		_ymax = _yBins->GetBinLowEdge(binY+1);
	} else {
		_ymin = _yBins->GetBinLowEdge(1);
		_ymax = _yBins->GetBinLowEdge(_yBins->GetNbinsX()+1);
	}

	std::cout << "INFO : fitting Z0 - mass fit to file " << _name << ", pT in (" << _jptmin << "," << _jptmax << "), y in (" << _ymin << "," << _ymax << ")" << std::endl;

	TString ptStr;
	ptStr+=_jptmin; ptStr+="-"; ptStr+=_jptmax;
	ptStr+=_ymin;  ptStr+="-"; ptStr+=_ymax;

	TFile* f0 = TFile::Open(_dataFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	RooRealVar ZM(  "ZM",  "ZM",  60e3,  120e3,  "MeV/#it{c}^{2}"); 
	ZM.setBins(100);
	RooRealVar JetPT( "JetPT",  "JetPT",  _jptmin,  _jptmax,  "MeV/#it{c}"); 
	RooRealVar    ZY( "ZY",     "ZY",     _ymin,    _ymax,  ""); 

	RooArgSet obs;
	obs.add(JetPT);
	obs.add(ZM);
	obs.add(ZY);

	RooDataSet* ds = new RooDataSet("ds","ds", obs, RooFit::Import(*t0));

	//fit model parameters
	RooRealVar zMass("zMass","m(Z)",90e3,80e3,100e3);
	RooRealVar zWidth1("zWidth1","width CB1",5e3,1e3,20e3);

	RooRealVar zRatio("zRatio","ratio", _ratio, 0.2, 5);
	RooFormulaVar zWidth2("zWidth2","width CB2","@0*@1",RooArgList(zWidth1,zRatio));
	RooRealVar zAlpha1("zAlpha1", "alpha1", _alpha1, 0.5, 5.0);
	RooRealVar zN1(    "zN1",     "N1", 2.);//2.0, 0.0, 5.0);
	RooRealVar zAlpha2("zAlpha2", "alpha2", _alpha2, -5.0, -0.5);
	RooRealVar zN2(    "zN2",     "N2", 2.);//2.0, 0.0, 5.0);
	RooRealVar zFrac("zFrac","frac", _frac, 0., 1.);

	if(_fixed) {
		zRatio.setConstant();
		zAlpha1.setConstant();
		zAlpha2.setConstant();
		zFrac.setConstant();
	}

	RooCBShape sig_cb1("sig_cb1","",ZM,zMass,zWidth1, zAlpha1, zN1);
	RooCBShape  sig_cb2("sig_cb2", "", ZM, zMass, zWidth2, zAlpha2, zN2); 
	RooAddPdf sig( "sig", "", RooArgList(sig_cb1,sig_cb2), RooArgList(zFrac) );
	//RooGaussian sig("sig","",ZM,zMass,zWidth);

	RooRealVar a("a","a",0., -0.1, 0.1);
	//RooExponential bkg("bkg","",ZM, a);
	RooPolynomial bkg("bkg","",ZM, RooArgList(a));

	double maxyield = ds->sumEntries();
	double fracsig = 0.9;
	std::cout << maxyield << "\t" << fracsig << std::endl;

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield Z",        fracsig *maxyield,    0.0,     1.1*maxyield);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", (1.-fracsig)*maxyield,    0.0,     1.1*maxyield);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sig,bkg), RooArgList(sigYield,bkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/Z_"+_name+"_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "sig" );
	sig_pdfs.push_back( "sig_cb1" );
	sig_pdfs.push_back( "sig_cb2" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "bkg" );

	plotFit(ZM, 60e3, 120e3, 100, ds, data_pdf, sig_pdfs, bkg_pdfs, "Z_"+_name+"_M_"+ptStr, "m");

	////print parameters
	RooArgList params;
	params.add(sigYield);
	params.add(bkgYield);
	params.add(zMass);
	params.add(zWidth1);
	params.add(zRatio);
	params.add(zAlpha1);
	params.add(zAlpha2);
	params.add(zFrac);
	//params.add(zN1);
	//params.add(zN2);
	params.add(a);
	printParams("Z_"+_name+"_params_"+ptStr+".dat",params);

	yieldZ= sigYield.getVal();
	errorZ= sigYield.getError();

	_ratio  = zRatio.getVal();
	_frac   = zFrac.getVal();
	_alpha1 = zAlpha1.getVal();
	_alpha2 = zAlpha2.getVal();

	return true;
}
