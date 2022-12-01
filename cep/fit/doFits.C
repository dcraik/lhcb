#include "doFits.h"

#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"

#include "RooExpWithSplineEff.h"
#include "RooExpWithEff.h"
#include "RooBkgShape.h"
#include "RooSpline.h"

#include "outputFunctions.h"

TString const Fitter::typeName(fitType type) {
	switch(type) {
		case fitSim:
			return "sim";
		case fitpp:
			return "pp";
		case fitppNoPrevEtCut:
			return "ppNoPrevEtCut";
		case fitpp5:
			return "pp5TeV";
		case fitpPb:
			return "pPb";
		case fitPbp:
			return "Pbp";
		case fitPbPb:
			return "PbPb";
		default:
			return "";
	}
}

TString const Fitter::optsName(uint options) {
	TString name="";
	if(options & optHiPTSq) {
		name+="-hiPTSq";
	}
	if(options & optHiHRC) {
		name+="-hiHRC";
	}
	if(options & optHiNTrk) {
		name+="-hiNTrk";
	}
	return name;
}

TString const Fitter::optsCut(uint options) {
	TString cutStr="";
	if(options & optHiPTSq) {
		cutStr +="&& phi_PTSq>3. && phi_PTSq<10.";
	} else {
		cutStr +="&& phi_PTSq>0.";// && phi_PTSq<5.";
	}
	//if(options & optHiHRC) {
	//	cutStr +="&& log_hrc_fom_v4>5. && log_hrc_fom_v4<10.";
	//} else {
	//	cutStr +="&& log_hrc_fom_v4>0. && log_hrc_fom_v4<5.";
	//}
	if(options & optHiNTrk) {
		cutStr +="&& (nVeloTracks>2. || nLongTracks>2.)";
	} else {
		cutStr +="&& nVeloTracks<=2. && nLongTracks<=2.";
	}
	return cutStr;
}

TTree* Fitter::getReducedTree(fitType which, uint options, bool massFit, /*bool hiPTSq, bool hiHRC, bool hiNTrk,*/ bool comb) {
	TTree *tin(0), *tout(0);
	TString dirName(""), cutStr("");

	switch(which) {
		case fitSim:
			dirName="sim";
			if(comb) return 0;
			tin = _simData;
			break;
		case fitpp:
		case fitppNoPrevEtCut:
			dirName="pp";
			if(comb) tin = _ppCmbData;
			else     tin = _ppData;
			break;
		case fitpp5:
			dirName="pp5";
			if(comb) tin = _pp5CmbData;
			else     tin = _pp5Data;
			break;
		case fitpPb:
			dirName="pPb";
			if(comb) tin = _ppbCmbData;
			else     tin = _ppbData;
			break;
		case fitPbp:
			dirName="Pbp";
			if(comb) tin = _pbpCmbData;
			else     tin = _pbpData;
			break;
		case fitPbPb:
			dirName="PbPb";
			if(comb) tin = _pbpbCmbData;
			else     tin = _pbpbData;
			break;
	}

	tin->SetBranchStatus("*", 0);
	tin->SetBranchStatus("phi_M", 1);
	tin->SetBranchStatus("phi_PTSq", 1);
	tin->SetBranchStatus("log_hrc_fom_v4", 1);
	tin->SetBranchStatus("L0DUTCK", 1);
	tin->SetBranchStatus("nVeloTracks", 1);
	tin->SetBranchStatus("nLongTracks", 1);
	tin->SetBranchStatus("Kplus_PT", 1);
	tin->SetBranchStatus("Kminus_PT", 1);
	tin->SetBranchStatus("Kplus_ProbNNk", 1);
	tin->SetBranchStatus("Kminus_ProbNNk", 1);
	tin->SetBranchStatus("Kplus_PIDK", 1);
	tin->SetBranchStatus("Kminus_PIDK", 1);

	if(massFit) {
		dirName+="-massFit";
		cutStr ="phi_M>1000. && phi_M<1050";
	} else {
		cutStr ="phi_M>";
		cutStr+=_valMu-3*_valSigma;
		cutStr+="&& phi_M<";
		cutStr+=_valMu+3*_valSigma;
	}
	dirName+=optsName(options);
	cutStr+=optsCut(options);
	//if(hiPTSq) {
	//	dirName+="-hiPTSq";
	//	cutStr +="&& phi_PTSq>3. && phi_PTSq<10.";
	//} else {
	//	cutStr +="&& phi_PTSq>0. && phi_PTSq<5.";
	//}
	//if(hiHRC) {
	//	dirName+="-hiHRC";
	//	cutStr +="&& log_hrc_fom_v3>5. && log_hrc_fom_v3<10.";
	//} else {
	//	cutStr +="&& log_hrc_fom_v3>0. && log_hrc_fom_v3<5.";
	//}
	//if(hiNTrk) {
	//	dirName+="-hiNTrk";
	//	cutStr +="&& (nVeloTracks>2. || nLongTracks>2.)";
	//} else {
	//	cutStr +="&& nVeloTracks<=2. && nLongTracks<=2.";
	//}
	if(which==fitpp) {
		cutStr +="&& L0DUTCK>0x1608";//0x1609,0x160e,0x160f,0x1611,0x1612
	} else if(which==fitppNoPrevEtCut) {
		dirName+="-noPrevEtCut";
		cutStr +="&& L0DUTCK<0x1606";//0x1603,0x1604,0x1605
	}
	if(comb) {
		dirName+="-comb";
	}

	if(!_fitTreeFile) return 0;
    	_fitTreeFile->cd();
	auto* dir = static_cast<TDirectory*>(_fitTreeFile->Get(dirName));
	if (!dir) {
		dir = _fitTreeFile->mkdir(dirName);
		dir = static_cast<TDirectory*>(_fitTreeFile->Get(dirName));
	}
	dir->cd();
	tout = tin->CopyTree(cutStr);
	tout->Write();

	tin->SetBranchStatus("*", 1);
	printf("Selected %d (%.1f%%) of %d entries\n",static_cast<int>(tout->GetEntries()),tout->GetEntries()*100./tin->GetEntries(),static_cast<int>(tin->GetEntries()));

	return tout;
}

bool Fitter::addGaussConstraint(RooArgSet& pdfs, RooLinkedList& fitCmds, RooAbsReal& var, double mean, double sigma) {
	TString name = var.GetName();

	RooRealVar* meanConstr  = new RooRealVar(name+"mean","",mean);
	RooRealVar* sigmaConstr = new RooRealVar(name+"sigma","",sigma);
	RooGaussian* constrain  = new RooGaussian(name+"constrain","",var,*meanConstr,*sigmaConstr);

	pdfs.add(*constrain);
	fitCmds.Add(new RooCmdArg{RooFit::Constrain(var)});

	return true;
}

void Fitter::fitEff() {
	TTree* t = _simData;
	if(!t) return;
	gSystem->RedirectOutput("log/eff_fit.log","w");

	TH1D den("den","",100,0.,10.);
	TH1D num("num","",100,0.,10.);
	TH1D eff("eff","",100,0.,10.);

	den.Sumw2();
	num.Sumw2();
	eff.Sumw2();

	t->Draw("phi_TRUEPT**2/1e6>>den");
	t->Draw("phi_PT**2/1e6>>num","phi_PT>0. && Kplus_P>5e3 && Kminus_P>5e3");

	eff.Divide(&num,&den);
	TH1* eff25 = eff.Rebin(4,"eff25");
	eff25->Scale(0.25);

	RooRealVar PTSq(  "phi_PTSq",  "phi_PTSq",  0., 10.,  "GeV^{2}/#it{c}^{2}"); 
	//RooSpline* spline = makeSpline("spline","spline",PTSq,std::vector<double>{0.,0.2,0.3,0.5,1.0,2.0,2.5,10.},std::vector<double>{0.,0.2,0.2,0.4,0.8,1.0,1.0,1.0},6u);
	RooSpline* spline = makeSpline("spline","spline",PTSq,5,eff25);

	RooRealVar norm("norm","norm",1.,0.,200.);

	RooAddPdf data_pdf("data_pdf","",RooArgList(*spline),RooArgList(norm));

	//TODO///
	//alt PDF
	RooRealVar a("a", "alpha", 0.); 
	RooRealVar x0("x0", "eff lower threshold",   _x0Val, 0.1, 0.5);
	RooRealVar x1("x1", "eff upper threshold",   _x1Val, 1.0, 5.0);
	RooRealVar gL("gL", "eff grad. left",        _gLVal, 2.0, 5.0);
	RooRealVar gR("gR", "eff grad. centre",      _gRVal, 0.2, 1.0);
	
	RooExpWithEff alt("alt", "alt", PTSq, a, x0, x1, gL, gR);
	RooAddPdf data_alt_pdf("data_alt_pdf","",RooArgList(alt),RooArgList(norm));
	/////////
	RooDataHist ds("ds",  "ds", RooArgList(PTSq), &eff);

	//fit y only
	data_alt_pdf.fitTo(ds);
	data_pdf.fitTo(ds);
	spline->fixXValues(false);
	spline->fixYValues();
	//fit x only
	data_pdf.fitTo(ds);
	spline->fixXValues();
	spline->fixYValues(false);
	//fit y only again
	data_pdf.fitTo(ds);
	RooPlot* frame = PTSq.frame();
	ds.plotOn(frame);
	data_pdf.plotOn(frame);
	data_alt_pdf.plotOn(frame);
	TCanvas c;
	frame->Draw();
	c.SaveAs("fig/fitSpline.pdf");

	spline->getXValues(_splineXVals);
	spline->getYValues(_splineYVals);
	printParams("dat/effFit.dat", &ds, &data_pdf);
	gSystem->RedirectOutput(0);
}

void Fitter::fitMass(fitType which, uint options, /*bool hiPTSq, bool hiHRC, bool hiNTrk,*/ bool comb) {
	TTree* t0(0);
	TString nameStr("_");
	nameStr+=typeName(which);
	nameStr+=optsName(options);
	//if(hiPTSq) nameStr+="_hiPTSq";
	//if(hiHRC)  nameStr+="_hiHRC";
	//if(hiNTrk)  nameStr+="_hiNTrk";
	if(comb)  nameStr+="_comb";

	t0=getReducedTree(which,options,true,/*hiPTSq,hiHRC,hiNTrk,*/comb);
	return;//TODO
	if(!t0) return;

	gSystem->RedirectOutput("log/mass_fit"+nameStr+".log","w");

	RooRealVar M(  "phi_M",  "phi_M",  1000.,  1050.,  "MeV/#it{c}^{2}"); 
	M.setBins(100);

	RooArgSet obs;
	obs.add(M);

	RooAbsData *ds(0);
	if(_doBinnedFit) {
		TH1D h0("h0","",200,1000.,1050.);

		t0->Draw("phi_M>>h0");
		ds = new RooDataHist("ds",  "ds",   obs, &h0);
	} else {
		ds = new RooDataSet("ds","ds", obs, RooFit::Import(*t0));
	}

	//fit model parameters
	RooRealVar sigMass ("sigMass", "sig mean",_valMu,1010.,1030.);
	//TODO//RooRealVar sigWidth1("sigWidth1","sig width",_valSigma,0.5,20.);

	//TODO//RooRealVar sigRatio("sigRatio","ratio of widths", 1.0);//_valR, 1.0, 1.5);
	//TODO//RooFormulaVar sigWidth2("sigWidth2","","@0*@1",RooArgList(sigWidth1,sigRatio));
	//TODO//RooRealVar sigAlpha1("sigAlpha1", "Alpha1", _valA1, 0.5, 5.0);
	//TODO//RooRealVar sigN1(    "sigN1",     "N1",     _valN1);
	//TODO//RooRealVar sigAlpha2("sigAlpha2", "Alpha2", _valA2, -5.0, -0.5);
	//TODO//RooRealVar sigN2(    "sigN2",     "N2",     _valN2);
	//TODO//RooRealVar sigFrac("sigFrac","fraction in CB1", _valF, 0., 1.);

	//TODO//RooCBShape  sig_cb1("sig_cb1", "", M, sigMass, sigWidth1, sigAlpha1, sigN1); 
	//TODO//RooCBShape  sig_cb2("sig_cb2", "", M, sigMass, sigWidth2, sigAlpha2, sigN2); 
	//TODO//RooAddPdf sig( "sig", "", RooArgList(sig_cb1,sig_cb2), RooArgList(sigFrac) );

	RooRealVar sigSigma("sigSigma","sig Gauss width",_valSigma,0.5,20.);
	RooRealVar sigGamma("sigGamma", "sig BW width",4.26);
	RooVoigtian  sig("sig", "", M, sigMass, sigGamma, sigSigma);

	//TODO//RooRealVar p0("p0","p0",0., -0.1, 0.1);
	//TODO//RooPolynomial bkg("bkg","",M, RooArgList(p0));
	RooRealVar c1("c1","c1",_valC1, -0.1, -0.0001);
	RooRealVar c2("c2","c2",_valC2, 900., 1000.);
	RooBkgShape bkg("bkg","",M, c1, c2);

	double maxyield = ds->sumEntries();
	double fracsig  = 0.9;

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield phi",      fracsig *maxyield,    0.0,     1.1*maxyield);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", (1.-fracsig)*maxyield,    0.0,     1.1*maxyield);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sig,bkg), RooArgList(sigYield,bkgYield) );

	if(comb) {
		sigYield.setVal(0);
		sigYield.setConstant();
		sigMass.setConstant();
		sigSigma.setConstant();
	}

	// -- fit model pdf to the dataset ----------------------------------------------
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));

	//update saved parameters
	_valMu = sigMass.getVal();
	_valSigma = sigSigma.getVal();
	//TODO//_valSigma = sigWidth1.getVal();
	//TODO//_valA1 = sigAlpha1.getVal();
	//TODO//_valA2 = sigAlpha2.getVal();
	//TODO//_valN1 = sigN1.getVal();
	//TODO//_valN2 = sigN2.getVal();
	//TODO////	_valR  = sigRatio.getVal();
	//TODO//_valF  = sigFrac.getVal();

	M.setRange("signal",_valMu-3*_valSigma,_valMu+3*_valSigma);
	M.setRange("sideLo",_valMu-6*_valSigma,_valMu-3*_valSigma);
	M.setRange("sideHi",_valMu+3*_valSigma,_valMu+6*_valSigma);
	M.setRange("full",  1000.,             1050.);
	printf("Using %.1f-%.1f for signal window\n", _valMu-3*_valSigma,_valMu+3*_valSigma);
	printf("Using %.1f-%.1f to estimate background\n", _valMu-6*_valSigma,_valMu+6*_valSigma);

	double fsig_1 = sig.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("signal"))->getVal();
	double fsig_2 = sig.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("sideLo"))->getVal();
	double fsig_3 = sig.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("sideHi"))->getVal();
	double fsig_0 = sig.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("full"))->getVal();

	double fbkg_1 = bkg.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("signal"))->getVal();
	double fbkg_2 = bkg.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("sideLo"))->getVal();
	double fbkg_3 = bkg.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("sideHi"))->getVal();
	double fbkg_0 = bkg.createIntegral(RooArgSet(M),RooFit::NormSet(M),RooFit::Range("full"))->getVal();

	std::cout << "yields in mass windows" << std::endl;
	std::cout << _valMu-6*_valSigma << "\t" << _valMu-3*_valSigma << "\t" << _valMu+3*_valSigma << "\t" << _valMu+6*_valSigma << std::endl;
	std::cout << sigYield.getVal()*fsig_1/fsig_0 << "\t" << bkgYield.getVal()*fbkg_1/fbkg_0 << std::endl;
	std::cout << sigYield.getVal()*(fsig_2+fsig_3)/fsig_0 << "\t" << bkgYield.getVal()*(fbkg_2+fbkg_3)/fbkg_0 << std::endl;
	std::cout << sigYield.getVal() << "\t" << bkgYield.getVal() << std::endl;

	if(comb) {
		_valC1 = c1.getVal();
		_valC2 = c2.getVal();
	} else {
		_bkgYield = bkgYield.getVal()*fbkg_1/fbkg_0;
		_bkgError = bkgYield.getError()*TMath::Sqrt(fbkg_1/fbkg_0);
	}

	//make plots
	std::vector<std::string> sig_pdfs;
	if(comb) sig_pdfs.push_back( "bkg" );
	else sig_pdfs.push_back( "sig" );
	std::vector<std::string> bkg_pdfs;
	if(!comb) bkg_pdfs.push_back( "bkg" );

	plotFit(M, 1000., 1050., 100, ds, data_pdf, sig_pdfs, bkg_pdfs, "fig/massFit"+nameStr, "#it{m}(#it{K}^{#plus}#it{K}^{#minus}) [MeV/#it{c}^{2}]");

	//////print parameters
	//RooArgList params;
	//params.add(sigYield);
	//params.add(bkgYield);
	//params.add(sigMass);
	//params.add(sigSigma);
	////TODO//params.add(sigWidth1);
	////TODO////	params.add(sigRatio);
	////TODO//params.add(sigAlpha1);
	////TODO//params.add(sigAlpha2);
	////TODO//params.add(sigFrac);
	////TODO//params.add(p0);
	//params.add(c1);
	//params.add(c2);
	//printParams("dat/massFit"+nameStr+".dat",params);
	printParams("dat/massFit"+nameStr+".dat", ds, &data_pdf);

	gSystem->RedirectOutput(0);
	return;
}

void Fitter::fitPtSq(fitType which, uint options) { //bool hiPTSq, bool hiHRC, bool hiNTrk) {
	TTree *t0(0), *t1(0);
	TString nameStr("_");
	nameStr+=typeName(which);
	nameStr+=optsName(options);
	//if(hiPTSq) nameStr+="_hiPTSq";
	//if(hiHRC)  nameStr+="_hiHRC";
	//if(hiNTrk)  nameStr+="_hiNTrk";

	double min(0.), max(5.);
	if(options & optHiPTSq) {
		min = 3.;
		max = 10.;
	}

	t0=getReducedTree(which,options,false,/*hiPTSq,hiHRC,hiNTrk,*/false);
	t1=getReducedTree(which,options,false,/*hiPTSq,hiHRC,hiNTrk,*/true);

	if(!t0) return;

	gSystem->RedirectOutput("log/ptsq_fit"+nameStr+".log","w");

	RooRealVar PTSq(  "phi_PTSq",  "phi_PTSq",  min,  max,  "GeV^{2}/#it{c}^{2}"); 
	PTSq.setBins(100);

	RooArgSet obs;
	obs.add(PTSq);

	RooAbsData *ds(0), *dsSS(0);
	if(_doBinnedFit) {
		TH1D h0("h0","",200,min,max);
		TH1D h1("h1","",200,min,max);

		t0->Draw("phi_PTSq>>h0");
		ds = new RooDataHist("ds",  "ds",   obs, &h0);
		if(t1) {
			t1->Draw("phi_PTSq>>h1");
			dsSS = new RooDataHist("dsSS","dsSS", obs, &h1);
		}
	} else {
		ds          = new RooDataSet("ds",  "ds",   obs, RooFit::Import(*t0));
		if(t1) dsSS = new RooDataSet("dsSS","dsSS", obs, RooFit::Import(*t1));
	}

	RooRealVar aSig("aSig", "alpha signal", -4.0, -20.0,  -1.); 
	RooRealVar aDis("aDis", "alpha proton dis.", _aDis, -5.,  -0.05); 

	RooAbsPdf *sig(0), *pDis(0);

	//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
	//   1  A            3.40597e-02   1.08719e-03   1.63542e-06  -1.47764e+00
	//   2  x0           3.31334e-01   1.81302e-02   1.35093e-05   1.44557e-01
	//   3  x1           3.50015e+00   1.42661e-01   4.84587e-04  -2.51537e-03
	//   4  g0           3.01892e+00   1.91345e-01   1.59954e-04   1.35935e-02
	//   5  g1           6.76119e-01   4.28981e-02   7.06607e-05  -2.39738e-02
	//
	//
	//   1  A            3.44250e-02   6.38658e-04   1.30861e-06  -2.07049e+00
	//   2  x0           3.13686e-01   9.89102e-03   1.15867e-05  -1.63405e-01
	//   3  x1           3.30997e+00   7.44192e-02   3.54055e-04  -7.81391e-03
	//   4  g0           3.08542e+00   1.12332e-01   1.52414e-04  -2.10145e-02
	//   5  g1           7.08739e-01   2.59935e-02   5.83149e-05  -8.85892e-02

	RooRealVar x0("x0", "eff lower threshold",   _x0Val, 0.1, 0.5);
	RooRealVar x1("x1", "eff upper threshold",   _x1Val, 1.0, 5.0);
	RooRealVar gL("gL", "eff grad. left",        _gLVal, 2.0, 5.0);
	RooRealVar gR("gR", "eff grad. centre",      _gRVal, 0.2, 1.0);


	double maxEff(0.);
	uint effMaxPoint(0);
	for(uint i=0; i<_splineYVals.size(); ++i) {
		if(_splineYVals.at(i)>maxEff) {
			maxEff = _splineYVals.at(i);
			effMaxPoint = i;
		}
	}
	RooSpline* effFunc = makeSpline("effFunc","effFunc",PTSq,_splineXVals,_splineYVals, effMaxPoint);
	effFunc->fixXValues();
	effFunc->fixYValues();
	//RooRealVar p0("p0","p0",0.);
	//RooRealVar p1("p1","p1",0.);//, -0.1, 0.1);
	//RooPolynomial effCor("effCor","",PTSq, RooArgList(p0));

	RooRealVar aDummy("aDummy", "alpha dummy", 0.); 
	RooExpWithEff eff2 ("eff2",  "eff2",             PTSq, aDummy, x0, x1, gL, gR);

	RooExponential sig1( "sig1",  "signal raw",             PTSq, aSig);
	RooExponential pDis1("pDis1", "proton disociation raw", PTSq, aDis);

	if(!_useSplineEff) {
		sig = new RooExpWithEff("sig",  "signal",             PTSq, aSig, x0, x1, gL, gR);
		pDis = new RooExpWithEff("pDis", "proton disociation", PTSq, aDis, x0, x1, gL, gR);
		//sig = new RooProdPdf("sig",   "signal", RooArgSet(sig1,eff2));//,effCor));
		//pDis = new RooProdPdf("pDis", "proton disociation", RooArgSet(pDis1,eff2));
	} else {
		sig = new RooExpWithSplineEff(*effFunc, aSig, "sig");
		pDis = new RooExpWithSplineEff(*effFunc, aDis, "pDis");
		//sig = new RooProdPdf("sig",   "signal", RooArgSet(sig1,*effFunc));//,effCor));
		//pDis = new RooProdPdf("pDis", "proton disociation", RooArgSet(pDis1,*effFunc));//,effCor));
	}

	RooAbsPdf* nonRes(0);
	if(_doBinnedFit) {
		nonRes = new RooHistPdf("nonRes","non-resonant", PTSq, *static_cast<RooDataHist*>(dsSS));
	} else {
		nonRes = new RooKeysPdf("nonRes","non-resonant", PTSq, *static_cast<RooDataSet*>(dsSS));
	}

	double maxyield = ds->sumEntries();

	RooRealVar sigYield(  "sigYield",  "yield CEP",                0.5*(maxyield-_bkgYield),    0.0,     1.1*maxyield);
	RooRealVar pDisYield(  "pDisYield",  "yield proton dis.",      0.5*(maxyield-_bkgYield),    0.0,     1.1*maxyield);
	RooRealVar nonResYield("nonResYield","yield non-res",          _bkgYield, 0.8*_bkgYield, 1.2*_bkgYield);

	if(options & optHiPTSq) {
		sigYield.setVal(0);
		sigYield.setConstant();
		aSig.setConstant();
	} else {
		aDis.setConstant();
	}

	// -- total PDF
	RooAddPdf data_pdf1( "data_pdf1",  "data_pdf1", RooArgList(*sig,*pDis,*nonRes), RooArgList(sigYield,pDisYield,nonResYield) );
	RooProdPdf* data_pdf(0);

	RooArgSet fitPdfs;
	RooLinkedList fitArgs;

	fitArgs.Add(new RooCmdArg{RooFit::Extended()});
	fitArgs.Add(new RooCmdArg{RooFit::Save()});
	fitArgs.Add(new RooCmdArg{RooFit::NumCPU(4)});
	fitArgs.Add(new RooCmdArg{RooFit::Range("FIT")});
	fitArgs.Add(new RooCmdArg{RooFit::SumW2Error(kTRUE)});

	fitPdfs.add(data_pdf1);

	addGaussConstraint(fitPdfs, fitArgs, nonResYield, _bkgYield, _bkgError);
	//RooRealVar meanConstrBkgYield("meanBkgYield","",_bkgYield);
	//RooRealVar sigmaConstrBkgYield("sigmaBkgYield","",_bkgError);
	//RooGaussian constrainBkgYield("constrainBkgYield","",nonResYield,meanConstrBkgYield,sigmaConstrBkgYield);
	//fitPdfs.add(constrainBkgYield);
	//fitArgs.Add(new RooCmdArg{RooFit::Constrain(nonResYield)});

	//RooRealVar meanConstrADis("meanADis","",_aDis);
	//RooRealVar sigmaConstrADis("sigmaADis","",_aDisError);
	//RooGaussian constrainADis("constrainADis","",aDis,meanConstrADis,sigmaConstrADis);

	if(!(options & optHiPTSq)) {
		addGaussConstraint(fitPdfs, fitArgs, aDis, _aDis, _aDisError);
		//fitPdfs.add(constrainADis);
		//fitArgs.Add(new RooCmdArg{RooFit::Constrain(aDis)});
	}

	//RooRealVar meanConstrx0("meanx0","",_x0Val);
	//RooRealVar sigmaConstrx0("sigmax0","",_x0Error);
	//RooGaussian constrainx0("constrainx0","",x0,meanConstrx0,sigmaConstrx0);
	//RooRealVar meanConstrx1("meanx1","",_x1Val);
	//RooRealVar sigmaConstrx1("sigmax1","",_x1Error);
	//RooGaussian constrainx1("constrainx1","",x1,meanConstrx1,sigmaConstrx1);
	//RooRealVar meanConstrgL("meangL","",_gLVal);
	//RooRealVar sigmaConstrgL("sigmagL","",_gLError);
	//RooGaussian constraingL("constraingL","",gL,meanConstrgL,sigmaConstrgL);
	//RooRealVar meanConstrgR("meangR","",_gRVal);
	//RooRealVar sigmaConstrgR("sigmagR","",_gRError);
	//RooGaussian constraingR("constraingR","",gR,meanConstrgR,sigmaConstrgR);
	if(!_useSplineEff) {
		if(_floatEffFunc) {
			addGaussConstraint(fitPdfs, fitArgs, x0, _x0Val, _x0Error);
			addGaussConstraint(fitPdfs, fitArgs, x1, _x1Val, _x1Error);
			addGaussConstraint(fitPdfs, fitArgs, gL, _gLVal, _gLError);
			addGaussConstraint(fitPdfs, fitArgs, gR, _gRVal, _gRError);
			//fitPdfs.add(constrainx0);
			//fitArgs.Add(new RooCmdArg{RooFit::Constrain(x0)});
			//fitPdfs.add(constrainx1);
			//fitArgs.Add(new RooCmdArg{RooFit::Constrain(x1)});
			//fitPdfs.add(constraingL);
			//fitArgs.Add(new RooCmdArg{RooFit::Constrain(gL)});
			//fitPdfs.add(constraingR);
			//fitArgs.Add(new RooCmdArg{RooFit::Constrain(gR)});
		} else {
			x0.setConstant();
			x1.setConstant();
			gL.setConstant();
			gR.setConstant();
		}
	}

	//data_pdf = new RooProdPdf("data_pdf", "data_pdf", RooArgSet(data_pdf1,constrainx0,constrainx1,constraingL,constraingR,constrainADis,constrainBkgYield));
	data_pdf = new RooProdPdf("data_pdf", "data_pdf", fitPdfs);
	/*RooFitResult * result =*/ data_pdf->fitTo(*ds, fitArgs);

	if(options & optHiPTSq) {
		_aDis  = aDis.getVal();
		_aDisError  = aDis.getError();
	}

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "sig" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "pDis" );
	bkg_pdfs.push_back( "nonRes" );

	plotFit(PTSq, min, max, 100, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/ptsqFit"+nameStr, "#it{p}_{T}^{2}(#it{K}^{#plus}#it{K}^{#minus}) [GeV^{2}/#it{c}^{2}]");

	//////print parameters
	//RooArgList params;
	//if(!hiPTSq && !hiHRC) params.add(sigYield);
	//params.add(pDisYield);
	//params.add(nonResYield);
	//if(hiPTSq || hiHRC) {
	//	params.add(aDis);
	//} else {
	//	params.add(aSig);
	//}
	//printParams("dat/ptsqFit"+nameStr+".dat",params);
	printParams("dat/ptsqFit"+nameStr+".dat", ds, data_pdf);

	gSystem->RedirectOutput(0);
	return;
}

void Fitter::setup() {
	_ppData = new TChain("DecayTree");
	_ppData->AddFile("/data/cep-phi/pp2phi.root");

	_pp5Data = new TChain("DecayTree");
	_pp5Data->AddFile("/data/cep-phi/pp52phi.root");

	_ppbData = new TChain("DecayTree");
	_ppbData->AddFile("/data/cep-phi/pPb2phi.root");

	_pbpData = new TChain("DecayTree");
	_pbpData->AddFile("/data/cep-phi/Pbp2phi.root");

	_pbpbData = new TChain("DecayTree");
	_pbpbData->AddFile("/data/cep-phi/PbPb2phi.root");

	_simData = new TChain("T");
	_simData->AddFile("/data/cep-phi/phi-sim.root");

	_ppCmbData = new TChain("DecayTree");
	_ppCmbData->AddFile("/data/cep-phi/pp2PiPi.root");

	_pp5CmbData = new TChain("DecayTree");
	_pp5CmbData->AddFile("/data/cep-phi/pp52PiPi.root");

	_ppbCmbData = new TChain("DecayTree");
	_ppbCmbData->AddFile("/data/cep-phi/pPb2PiPi.root");

	_pbpCmbData = new TChain("DecayTree");
	_pbpCmbData->AddFile("/data/cep-phi/Pbp2PiPi.root");

	_pbpbCmbData = new TChain("DecayTree");
	_pbpbCmbData->AddFile("/data/cep-phi/PbPb2PiPi.root");

	_fitTreeFile = TFile::Open("forFits.root", "RECREATE");
}

void Fitter::run() {
	if(_useSplineEff) fitEff();
	int stage(0);
	auto fits = std::vector<fitType>{/*fitppNoPrevEtCut,fitpp,fitpp5,*/fitpPb,fitPbp,fitPbPb};
	auto opts = std::vector<uint>{/*optHiPTSq|optHiHRC|optHiNTrk,optHiHRC|optHiNTrk,optHiPTSq|optHiNTrk,optHiNTrk,optHiPTSq|optHiHRC,optHiHRC,optHiPTSq,*/optNone};
	for( auto it=fits.begin(); it!=fits.end(); ++it) {
		TString name = typeName(*it);
		for( auto it2=opts.begin(); it2!=opts.end(); ++it2) {
			TString optsStr = optsName(*it2);
			printf("STAGE %02d : fit mass for %s, comb %s\n", ++stage, optsStr.Data(), name.Data());
			fitMass(*it,*it2,true);
			printf("STAGE %02d : fit mass for %s %s\n", ++stage, optsStr.Data(), name.Data());
			fitMass(*it,*it2);
			//TODO//printf("STAGE %02d : fit pT^2 for %s %s\n", ++stage, optsStr.Data(), name.Data());
			//TODO//fitPtSq(*it,*it2);
		}
		//printf("STAGE %02d : fit mass for high pT, high HRC, comb %s\n", ++stage, name.Data());
		//fitMass(*it,true,true,false,true);
		//printf("STAGE %02d : fit mass for high pT, high HRC %s\n", ++stage, name.Data());
		//fitMass(*it,true,true);
		//printf("STAGE %02d : fit pT^2 for high pT, high HRC %s\n", ++stage, name.Data());
		//fitPtSq(*it,true,true);
		//printf("STAGE %02d : fit mass for low pT, high HRC, comb %s\n", ++stage, name.Data());
		//fitMass(*it,false,true,false,true);
		//printf("STAGE %02d : fit mass for low pT, high HRC %s\n", ++stage, name.Data());
		//fitMass(*it,false,true);
		//printf("STAGE %02d : fit pT^2 for low pT, high HRC %s\n", ++stage, name.Data());
		//fitPtSq(*it,false,true);
		//printf("STAGE %02d : fit mass for high pT, low HRC, comb %s\n", ++stage, name.Data());
		//fitMass(*it,true,false,false,true);
		//printf("STAGE %02d : fit mass for high pT, low HRC %s\n", ++stage, name.Data());
		//fitMass(*it,true);
		//printf("STAGE %02d : fit pT^2 for high pT, low HRC %s\n", ++stage, name.Data());
		//fitPtSq(*it,true);
		//printf("STAGE %02d : fit mass of signal region, comb %s\n", ++stage, name.Data());
		//fitMass(*it,false,false,false,true);
		//printf("STAGE %02d : fit mass of signal region %s\n", ++stage, name.Data());
		//fitMass(*it);
		//printf("STAGE %02d : fit pT^2 of signal region %s\n", ++stage, name.Data());
		//fitPtSq(*it);
		printf("DONE: %s\n", name.Data());
	}
}

int main() {
	Fitter f;
	f.setup();
	f.run();
	return 0;
}
