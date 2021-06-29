#include "SimDFitter.h"

#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "RooAbsDataStore.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimWSTool.h"
#include "RooSimultaneous.h"
#include "RooThresholdCategory.h"
#include "RooWorkspace.h"

#include "RooMultiKeysPdf.h"
#include "RooPromptShape.h"

#include <boost/progress.hpp>

#include "AdaptBin.h"
#include "cloneHists.h"
#include "DatasetManager.h"
#include "outputFunctions.h"

void SimDFitter::setInputs(TString data, TString charm, TString beauty, TString eff, TString acc, bool isMC) {
	_dataFile  = data;
	_combFile  = data;
	_charmFile  = charm;
	_beautyFile  = beauty;
	_effFile   = eff;
	_accFile   = acc;
	_dataIsMC = isMC;

	_inputsSet=true;
	_effsSet = false;
	_initialised = false;
	_mcFitsDone = false;
}

void SimDFitter::setOptions(DFitterOptions& options) {
	_dptmin   = options.minpt;
	_dptmax   = options.maxpt;
	_skipSumW2Fits = options.skipSumw2Fits;
	_truthMatchData = static_cast<truthMatchType>(options.truthMatch);
	_nDPtBins = options.nDPtBins;
	_usePtFracBins = options.usePtFracBins;
	_combShapeSyst = options.combShapeSyst;
	_doEnhancedFits = options.doEnhancedFits;
	_allowMassWidthScale   = options.allowMassWidthScale;
	_allowPromptMeanShift  = options.allowPromptMeanShift;
	_allowPromptWidthScale = options.allowPromptWidthScale;
	_allowDisplcMeanShift  = options.allowDisplcMeanShift;
	_allowDisplcWidthScale = options.allowDisplcWidthScale;
	_enhanceMassWidthScale   = options.enhanceMassWidthScale;
	_enhancePromptMeanShift  = options.enhancePromptMeanShift;
	_enhancePromptWidthScale = options.enhancePromptWidthScale;
	_enhanceDisplcMeanShift  = options.enhanceDisplcMeanShift;
	_enhanceDisplcWidthScale = options.enhanceDisplcWidthScale;
	_splitMassWidthScale   = options.splitMassWidthScale;
	_splitPromptMeanShift  = options.splitPromptMeanShift;
	_splitPromptWidthScale = options.splitPromptWidthScale;
	_splitDisplcMeanShift  = options.splitDisplcMeanShift;
	_splitDisplcWidthScale = options.splitDisplcWidthScale;
	_splitBkgMassShape = options.splitBkgMassShape;
	_splitBkgIPShape = options.splitBkgIPShape;
	_allowLinearBkgMassShift=options.allowLinearBkgMassShift;
	_setMassWidthScale   = options.setMassWidthScale;
	_setPromptMeanShift  = options.setPromptMeanShift;
	_setPromptWidthScale = options.setPromptWidthScale;
	_setDisplcMeanShift  = options.setDisplcMeanShift;
	_setDisplcWidthScale = options.setDisplcWidthScale;
	_shiftFixedMassWidth   = options.shiftFixedMassWidth;
	_shiftFixedPromptMean  = options.shiftFixedPromptMean;
	_shiftFixedPromptWidth = options.shiftFixedPromptWidth;
	_shiftFixedDisplcMean  = options.shiftFixedDisplcMean;
	_shiftFixedDisplcWidth = options.shiftFixedDisplcWidth;
	_runToyFits = options.runToyFits;
	_useSimpleEffs = options.useSimpleEffs;
	_usePtBinCorrections = options.usePtBinCorrections;
	_useJetPtBinnedEffs = options.useJetPtBinnedEffs;
	std::cout << _nDPtBins << "!" << std::endl;//TODO
}

void SimDFitter::setDPtRange(double ptmin, double ptmax) {
	if(ptmin!=_dptmin || ptmax!=_dptmax) {
		_effsSet = false;
		_initialised = false;
		_mcFitsDone = false;
	}
	_dptmin   = ptmin;
	_dptmax   = ptmax;
}

void SimDFitter::setDPtBins(TH1D* ptBinsD) {
	_effsSet = false;
	_initialised = false;
	_mcFitsDone = false;
	_useFixedDBins = true;
	_nDPtBins = ptBinsD->GetNbinsX();

	if(ptBinsD) {
		_ptBinsD.clear();
		for( auto bin=1; bin<=ptBinsD->GetNbinsX()+1; ++bin) {
			_ptBinsD.push_back(ptBinsD->GetBinLowEdge(bin));
		}
	}
}

TString const SimDFitter::dFileName() {
	return gSaveDir+"/"+_whichD+"jets"+_name+".root";
}

void SimDFitter::initJetBinning(TH1D* ptBins, TH1D* yBins) {
	if(ptBins) {
		_ptBins.clear();
		for( auto bin=1; bin<=ptBins->GetNbinsX()+1; ++bin) {
			_ptBins.push_back(ptBins->GetBinLowEdge(bin));
		}
	}
	if(yBins) {
		_yBins.clear();
		for( auto bin=1; bin<=yBins->GetNbinsX()+1; ++bin) {
			_yBins.push_back(yBins->GetBinLowEdge(bin));
		}
	}
}

void SimDFitter::setupMasses() {
	if(_whichD=="D0" || _whichD=="K3PI") {
		_mCentre=1864.;
		_mLow=1784.;
		_mHigh=1944.;
		_mSBLow=1844.;
		_mSBHigh=1884.;
		_mSBLooseLow=1834.;
		_mSBLooseHigh=1894.;
		_decayLabel="#it{K}^{#minus}#pi^{#plus}";
	} else if(_whichD=="D") {
		_mCentre=1870.;
		_mLow=1790.;
		_mHigh=1950.;
		_mSBLow=1850.;
		_mSBHigh=1890.;
		_mSBLooseLow=1840.;
		_mSBLooseHigh=1900.;
		_decayLabel="#it{K}^{#minus}#pi^{#plus}#pi^{#plus}";
	}
}

bool SimDFitter::init() {
	if(ws) delete ws;
	ws = new RooWorkspace("ws");
	bool status(true);

	status&=loadDatasets();
	status&=setupModel();
	_initialised = true;
	_mcFitsDone = false;

	return status;
}

bool SimDFitter::addEffs() {
	if(!_inputsSet) return false;

	auto status(true);

	//if file already exists and we're not recreating then skip this step
	if(_recreateInputs || gSystem->AccessPathName(dFileName())) {
		status &= addEffs(0);
		status &= addEffs(4);
		status &= addEffs(5);
		status &= addEffs(7);
	}

	if(status) {
		_effsSet = true;
		_initialised = false;
		_mcFitsDone = false;
	}

	return status;
}

bool SimDFitter::addEffs(int flavour) {
	if(_useJetPtBinnedEffs) return addEffsPtBinned(flavour);

	std::cout << "INFO : adding " << _whichD << " efficiencies to file " << _name << ": " << flavour << std::endl;

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);
	TFile* f4 = TFile::Open("../efficiencies/"+_whichD+"4pTEffCorrections.root");
	TFile* f5 = TFile::Open("../efficiencies/"+_whichD+"5pTEffCorrections.root");

	bool isMC = _dataIsMC;

	DatasetManager* dm = DatasetManager::getInstance();
	switch(flavour) {
		case 0:
			dm->setDataset(_combFile);
			break;
		case 4:
			dm->setDataset(_charmFile);
			isMC=true;
			break;
		case 5:
			dm->setDataset(_beautyFile);
			isMC=true;
			break;
		case 7:
			dm->setDataset(_dataFile);
			break;
		default:
			return false;
	}

	//construct efficiency objects
	TH2D* hnumAR= dynamic_cast<TH2D*>(f3->Get("reco"));
	TH2D* hdenAR= dynamic_cast<TH2D*>(f3->Get("denom"));
	TH2D* hnumS4= dynamic_cast<TH2D*>(f2->Get("num"+_whichD+"4"));
	TH2D* hdenS4= dynamic_cast<TH2D*>(f2->Get("reco"+_whichD+"4"));
	TH2D* hnumS5= dynamic_cast<TH2D*>(f2->Get("num"+_whichD+"5"));
	TH2D* hdenS5= dynamic_cast<TH2D*>(f2->Get("reco"+_whichD+"5"));
	TH2D* hden4 = dynamic_cast<TH2D*>(f2->Get("denom"+_whichD+"4"));
	TH2D* hden5 = dynamic_cast<TH2D*>(f2->Get("denom"+_whichD+"5"));
//TODO//20201117 - old reco eff	TH2D* hnumAR= dynamic_cast<TH2D*>(f3->Get("num"));
//TODO//20201117 - old reco eff	TH2D* hdenAR= dynamic_cast<TH2D*>(f3->Get("denom"));
//TODO//20201117 - old reco eff	TH2D* hnumS4= dynamic_cast<TH2D*>(f2->Get("numD04"));
//TODO//20201117 - old reco eff	TH2D* hdenS4= dynamic_cast<TH2D*>(f2->Get("denomD04"));
//TODO//20201117 - old reco eff	TH2D* hnumS5= dynamic_cast<TH2D*>(f2->Get("numD05"));
//TODO//20201117 - old reco eff	TH2D* hdenS5= dynamic_cast<TH2D*>(f2->Get("denomD05"));

	if(_useSimpleEffs) {
		if(!hnumS4 || !hden4) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnumS4 <<" "<< hden4 << std::endl;
			std::cout << "4 histograms missing" << std::endl;
			return false;
		}
		if(!hnumS5 || !hden5) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnumS5 <<" "<< hden5 << std::endl;
			std::cout << "5 histograms missing" << std::endl;
			return false;
		}
	} else {
		if(!hnumAR || !hdenAR) {
			std::cout << "!" << _whichD << std::endl;
			f3->ls();
			std::cout << hnumAR <<" "<< hdenAR << std::endl;
			std::cout << "AccRec histograms missing" << std::endl;
			return false;
		}
		if(!hnumS4 || !hdenS4) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnumS4 <<" "<< hdenS4 << std::endl;
			std::cout << "SelPid4 histograms missing" << std::endl;
			return false;
		}
		if(!hnumS5 || !hdenS5) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnumS5 <<" "<< hdenS5 << std::endl;
			std::cout << "SelPid5 histograms missing" << std::endl;
			return false;
		}
	}

	std::vector<TH1D*> ptBinCorrections4;
	ptBinCorrections4.push_back(static_cast<TH1D*>(f4->Get("corrPtBin1")));
	ptBinCorrections4.push_back(static_cast<TH1D*>(f4->Get("corrPtBin2")));
	ptBinCorrections4.push_back(static_cast<TH1D*>(f4->Get("corrPtBin3")));
	ptBinCorrections4.push_back(static_cast<TH1D*>(f4->Get("corrPtBin4")));
	std::vector<TH1D*> ptBinCorrections5;
	ptBinCorrections5.push_back(static_cast<TH1D*>(f5->Get("corrPtBin1")));
	ptBinCorrections5.push_back(static_cast<TH1D*>(f5->Get("corrPtBin2")));
	ptBinCorrections5.push_back(static_cast<TH1D*>(f5->Get("corrPtBin3")));
	ptBinCorrections5.push_back(static_cast<TH1D*>(f5->Get("corrPtBin4")));

	TEfficiency heffAR(*hnumAR,*hdenAR);
	TEfficiency heffSP4(*hnumS4,*hdenS4);
	TEfficiency heffSP5(*hnumS5,*hdenS5);
	TEfficiency heffSimple4(*hnumS4,*hden4);
	TEfficiency heffSimple5(*hnumS5,*hden5);

	//compare against hand-rolled efficiencies
	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
//TODO	TH2D* hacc = dynamic_cast<TH2D*>(f2->Get("acceffD045"));
	TH2D* hacc4= dynamic_cast<TH2D*>(f2->Get("acceff"+_whichD+"4"));
	TH2D* hacc5= dynamic_cast<TH2D*>(f2->Get("acceff"+_whichD+"5"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
	TH2D* hrec4 = dynamic_cast<TH2D*>(f2->Get("receff"+_whichD+"4"));
	TH2D* hrec5 = dynamic_cast<TH2D*>(f2->Get("receff"+_whichD+"5"));
	//the following histogram give a vertex-position-dependent correction to the two-track efficiency
	//TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid(0);
	//if(dataIsMC) hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));
	//else hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	hpid = dynamic_cast<TH2D*>(f2->Get("pideff"+_whichD+"45"));
	TH2D* hsel4 = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"4"));
	TH2D* hsel5 = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"5"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045")); //superseded by RapidSim version
//TODO//20201117 - old reco eff	TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));//TODO//20201112 - switch back to old reco eff
//TODO//20201117 - old reco eff	TH2D* hrec4 = dynamic_cast<TH2D*>(f2->Get("receffD04"));//TODO//20201112 - switch back to old reco eff
//TODO//20201117 - old reco eff	TH2D* hrec5 = dynamic_cast<TH2D*>(f2->Get("receffD05"));//TODO//20201112 - switch back to old reco eff

	if(!hacc || !hrec /*|| !hcor*/ || !hpid || !hsel4 || !hsel5) {
		std::cout << "comparison histograms missing" << std::endl;
		return false;
	}

	//plot efficiency functions used
	TCanvas c1;
	hacc->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	hacc->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hrec->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	hrec->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	//hcor->GetXaxis()->SetTitle("#it{#rho}^{2}(#it{D}^{0}) [mm^{2}]");
	//hcor->GetYaxis()->SetTitle("#it{z}(#it{D}^{0}) [mm]");
	hsel4->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hsel4->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hsel5->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hsel5->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hpid->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hpid->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	if(_useSimpleEffs) {
		hacc4->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffAccPrompt.pdf");
		hacc5->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffAccDispl.pdf");
	} else {
		hacc->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffAcc.pdf");
	}
	if(_useSimpleEffs) {
		hrec4->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffRecPrompt.pdf");
		hrec5->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffRecDispl.pdf");
	} else {
		hrec->Draw("colz");
		c1.SaveAs(gSaveDir+"/"+_whichD+"EffRec.pdf");
	}
	//hcor->Draw("colz");
	//c1.SaveAs(gSaveDir+"/D0EffRecCor.pdf");
	hsel4->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"EffSelPrompt.pdf");
	hsel5->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"EffSelDispl.pdf");
	hpid->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"EffPID.pdf");

	std::vector<double> *vDM = new std::vector<double>();
	std::vector<double> *vDIPCHI2 = new std::vector<double>();
	std::vector<double> *vDPT = new std::vector<double>();
	std::vector<double> *vDPX = new std::vector<double>();
	std::vector<double> *vDPY = new std::vector<double>();
	std::vector<double> *vDPZ = new std::vector<double>();
	std::vector<double> *vDE = new std::vector<double>();
	std::vector<double> *vDX = new std::vector<double>();
	std::vector<double> *vDY = new std::vector<double>();
	std::vector<double> *vDZ = new std::vector<double>();
	std::vector<double> *vDh0P = new std::vector<double>();
	std::vector<double> *vDh0PT = new std::vector<double>();
	std::vector<double> *vDh0PX = new std::vector<double>();
	std::vector<double> *vDh0PY = new std::vector<double>();
	std::vector<double> *vDh0PZ = new std::vector<double>();
	std::vector<double> *vDh1P = new std::vector<double>();
	std::vector<double> *vDh1PT = new std::vector<double>();
	std::vector<double> *vDh1PX = new std::vector<double>();
	std::vector<double> *vDh1PY = new std::vector<double>();
	std::vector<double> *vDh1PZ = new std::vector<double>();
	std::vector<double> *vDh2P = new std::vector<double>();
	std::vector<double> *vDh2PT = new std::vector<double>();
	std::vector<double> *vDh2PX = new std::vector<double>();
	std::vector<double> *vDh2PY = new std::vector<double>();
	std::vector<double> *vDh2PZ = new std::vector<double>();
	std::vector<double> *vDh3P = new std::vector<double>();
	std::vector<double> *vDh3PT = new std::vector<double>();
	std::vector<double> *vDh3PX = new std::vector<double>();
	std::vector<double> *vDh3PY = new std::vector<double>();
	std::vector<double> *vDh3PZ = new std::vector<double>();
	std::vector<double> *vDh0PNN = new std::vector<double>();
	std::vector<double> *vDh1PNN = new std::vector<double>();
	std::vector<double> *vDh2PNN = new std::vector<double>();
	std::vector<double> *vDh3PNN = new std::vector<double>();
	std::vector<double> *vDh0WEIGHT = new std::vector<double>();
	std::vector<double> *vDh1WEIGHT = new std::vector<double>();
	std::vector<double> *vDh2WEIGHT = new std::vector<double>();
	std::vector<double> *vDh3WEIGHT = new std::vector<double>();
	std::vector<double> *vDTRUEIDX = new std::vector<double>();
	std::vector<double> *vDFROMB = new std::vector<double>();

	std::vector<double> *vTSVMCor = new std::vector<double>();
	std::vector<double> *vTSVN = new std::vector<double>();

	double JetPT;
	double JetEta;
	double JetTruePT;

	dm->setBranchAddress(""+_whichD+"M",           &vDM);
	dm->setBranchAddress(""+_whichD+"IPCHI2",      &vDIPCHI2);
	dm->setBranchAddress(""+_whichD+"PT",          &vDPT);
	dm->setBranchAddress(""+_whichD+"PX",          &vDPX);
	dm->setBranchAddress(""+_whichD+"PY",          &vDPY);
	dm->setBranchAddress(""+_whichD+"PZ",          &vDPZ);
	dm->setBranchAddress(""+_whichD+"E",           &vDE);
	dm->setBranchAddress(""+_whichD+"X",           &vDX);
	dm->setBranchAddress(""+_whichD+"Y",           &vDY);
	dm->setBranchAddress(""+_whichD+"Z",           &vDZ);
	if(_whichD=="D0") {
		dm->setBranchAddress("D0KP",          &vDh0P);
		dm->setBranchAddress("D0KPT",         &vDh0PT);
		dm->setBranchAddress("D0KPX",         &vDh0PX);
		dm->setBranchAddress("D0KPY",         &vDh0PY);
		dm->setBranchAddress("D0KPZ",         &vDh0PZ);
		dm->setBranchAddress("D0PIP",         &vDh1P);
		dm->setBranchAddress("D0PIPT",        &vDh1PT);
		dm->setBranchAddress("D0PIPX",        &vDh1PX);
		dm->setBranchAddress("D0PIPY",        &vDh1PY);
		dm->setBranchAddress("D0PIPZ",        &vDh1PZ);
		dm->setBranchAddress("D0KPNNK",       &vDh0PNN);
		dm->setBranchAddress("D0PIPNNPI",     &vDh1PNN);
		dm->setBranchAddress("D0KWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("D0PIWEIGHT",    &vDh1WEIGHT);
	} else if(_whichD=="D") {
		dm->setBranchAddress("DKP",          &vDh0P);
		dm->setBranchAddress("DKPT",         &vDh0PT);
		dm->setBranchAddress("DKPX",         &vDh0PX);
		dm->setBranchAddress("DKPY",         &vDh0PY);
		dm->setBranchAddress("DKPZ",         &vDh0PZ);
		dm->setBranchAddress("DPI1P",        &vDh1P);
		dm->setBranchAddress("DPI1PT",       &vDh1PT);
		dm->setBranchAddress("DPI1PX",       &vDh1PX);
		dm->setBranchAddress("DPI1PY",       &vDh1PY);
		dm->setBranchAddress("DPI1PZ",       &vDh1PZ);
		dm->setBranchAddress("DPI2P",        &vDh2P);
		dm->setBranchAddress("DPI2PT",       &vDh2PT);
		dm->setBranchAddress("DPI2PX",       &vDh2PX);
		dm->setBranchAddress("DPI2PY",       &vDh2PY);
		dm->setBranchAddress("DPI2PZ",       &vDh2PZ);
		dm->setBranchAddress("DKPNNK",       &vDh0PNN);
		dm->setBranchAddress("DPI1PNNPI",    &vDh1PNN);
		dm->setBranchAddress("DPI2PNNPI",    &vDh2PNN);
		dm->setBranchAddress("DKWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("DPI1WEIGHT",   &vDh1WEIGHT);
		dm->setBranchAddress("DPI2WEIGHT",   &vDh2WEIGHT);
	} else if(_whichD=="K3PI") {
		dm->setBranchAddress("K3PIKP",          &vDh0P);
		dm->setBranchAddress("K3PIKPT",         &vDh0PT);
		dm->setBranchAddress("K3PIKPX",         &vDh0PX);
		dm->setBranchAddress("K3PIKPY",         &vDh0PY);
		dm->setBranchAddress("K3PIKPZ",         &vDh0PZ);
		dm->setBranchAddress("K3PIPI1P",        &vDh1P);
		dm->setBranchAddress("K3PIPI1PT",       &vDh1PT);
		dm->setBranchAddress("K3PIPI1PX",       &vDh1PX);
		dm->setBranchAddress("K3PIPI1PY",       &vDh1PY);
		dm->setBranchAddress("K3PIPI1PZ",       &vDh1PZ);
		dm->setBranchAddress("K3PIPI2P",        &vDh2P);
		dm->setBranchAddress("K3PIPI2PT",       &vDh2PT);
		dm->setBranchAddress("K3PIPI2PX",       &vDh2PX);
		dm->setBranchAddress("K3PIPI2PY",       &vDh2PY);
		dm->setBranchAddress("K3PIPI2PZ",       &vDh2PZ);
		dm->setBranchAddress("K3PIPI3P",        &vDh3P);
		dm->setBranchAddress("K3PIPI3PT",       &vDh3PT);
		dm->setBranchAddress("K3PIPI3PX",       &vDh3PX);
		dm->setBranchAddress("K3PIPI3PY",       &vDh3PY);
		dm->setBranchAddress("K3PIPI3PZ",       &vDh3PZ);
		dm->setBranchAddress("K3PIKPNNK",       &vDh0PNN);
		dm->setBranchAddress("K3PIPI1PNNPI",    &vDh1PNN);
		dm->setBranchAddress("K3PIPI2PNNPI",    &vDh2PNN);
		dm->setBranchAddress("K3PIPI3PNNPI",    &vDh3PNN);
		dm->setBranchAddress("K3PIKWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("K3PIPI1WEIGHT",   &vDh1WEIGHT);
		dm->setBranchAddress("K3PIPI2WEIGHT",   &vDh2WEIGHT);
		dm->setBranchAddress("K3PIPI3WEIGHT",   &vDh3WEIGHT);
	}

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	dm->setBranchAddress("TSVMCor",       &vTSVMCor);
	dm->setBranchAddress("TSVN",          &vTSVN);
	if(isMC) {
		dm->setBranchAddress("JetTruePT",     &JetTruePT);
		dm->setBranchAddress(""+_whichD+"TRUEIDX",     &vDTRUEIDX);
		dm->setBranchAddress(""+_whichD+"FROMB",       &vDFROMB);
	}

	unsigned int nentries0 = dm->getEntries();

	TFile* fout = TFile::Open(dFileName(),"UPDATE");
	TString tname = "T";
	tname+=flavour;
	TTree* tout = new TTree(tname,"");

	double DM(0.), DPT(0.), DEta(0.), DLogIPChi2(0.), weight4(0.), weight5(0.), werror4(0.), werror5(0.);
	double weightAcc(0.), werrorAcc(0.), wbinAcc(0.);
	double weightSel4(0.), werrorSel4(0.), wbinSel4(0.);
	double weightSel5(0.), werrorSel5(0.), wbinSel5(0.);
	double enhanced(0.);
	//TODO//double DRhoSq(0.), DZ(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch(""+_whichD+"M",          &DM);
	tout->Branch(""+_whichD+"PT",         &DPT);
	tout->Branch(""+_whichD+"Eta",        &DEta);
	tout->Branch(""+_whichD+"LogIPChi2",  &DLogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);
	tout->Branch("werror4",      &werror4);
	tout->Branch("werror5",      &werror5);
	tout->Branch("weightAcc",    &weightAcc);
	tout->Branch("weightSel4",   &weightSel4);
	tout->Branch("weightSel5",   &weightSel5);
	tout->Branch("werrorAcc",    &werrorAcc);
	tout->Branch("werrorSel4",   &werrorSel4);
	tout->Branch("werrorSel5",   &werrorSel5);
	tout->Branch("wbinAcc",      &wbinAcc);
	tout->Branch("wbinSel4",     &wbinSel4);
	tout->Branch("wbinSel5",     &wbinSel5);
	tout->Branch("enhanced",     &enhanced);

	tout->Branch("JetTruePT", &JetTruePT);

	uint total(0), skipEff4(0), skipEff5(0), skipEff45(0);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	//for(unsigned int ientry=0; ientry<nentries0; ++ientry)
	while(dm->getNext()) {
		++progress;
		//t0->GetEntry(ientry);

		for(unsigned int s=0; s<vDM->size(); ++s) {
			//first check in our tight acceptance
			TVector3 DP (vDPX->at(s)  ,vDPY->at(s)  ,vDPZ->at(s));
			TVector3 DP0(vDh0PX->at(s) ,vDh0PY->at(s) ,vDh0PZ->at(s));
			TVector3 DP1(vDh1PX->at(s),vDh1PY->at(s),vDh1PZ->at(s));
			TVector3 DP2;//(vDh2PX->at(s),vDh2PY->at(s),vDh2PZ->at(s));
			TVector3 DP3;//(vDh3PX->at(s),vDh3PY->at(s),vDh3PZ->at(s));
			if(_whichD=="D" || _whichD=="K3PI") {
				DP2.SetXYZ(vDh2PX->at(s),vDh2PY->at(s),vDh2PZ->at(s));
				if(_whichD=="K3PI") {
					DP3.SetXYZ(vDh3PX->at(s),vDh3PY->at(s),vDh3PZ->at(s));
				}
			}

			if(flavour==7 && isMC && _truthMatchData!=truthMatchType::truthMatchOff) {
				if(_truthMatchData==truthMatchType::truthMatchPrompt) {
					if(vDTRUEIDX->at(s)==-1 || vDFROMB->at(s)==1) continue;
				}
				if(_truthMatchData==truthMatchType::truthMatchDispl) {
					if(vDTRUEIDX->at(s)==-1 || vDFROMB->at(s)==0) continue;
				}
			}

			if(!(DP0.Eta()>2. && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
			if(!(DP1.Eta()>2. && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
			if(_whichD=="D" || _whichD=="K3PI") {
				if(!(DP2.Eta()>2. && DP2.Eta()<4.5 && DP2.Pt()>500. && DP2.Mag()>5000.)) continue;
				if(_whichD=="K3PI") {
					if(!(DP3.Eta()>2. && DP3.Eta()<4.5 && DP3.Pt()>500. && DP3.Mag()>5000.)) continue;
				}
			}
			if(!(DP.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
			if((flavour==4 || flavour==5) && vDTRUEIDX->at(s)==-1) continue;
			if(flavour==5 && vDFROMB->at(s)==0) continue;
			if(flavour==4 && vDFROMB->at(s)==1) continue;
			if(flavour==0 && vDM->at(s)>_mSBLow && vDM->at(s)<_mSBHigh) continue;

			//check PID cuts
//TODO//20200420			if(!isMC && vDh0PNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			if(_whichD=="D0" || _whichD=="D" || _whichD=="K3PI") {
				if(vDh0PNN->at(s)<0.2 && DP0.Pt()<25000. && DP0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
				//if(!dataIsMC && vDh1PNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off
			}

			DM         = vDM->at(s);
			DLogIPChi2 = TMath::Log(vDIPCHI2->at(s));
			DPT        = DP.Pt();
			DEta       = DP.Eta();
			//TODO//DRhoSq     = vDX->at(s)*vDX->at(s) + vDY->at(s)*vDY->at(s);
			//TODO//DZ         = vDZ->at(s);

			//for data, check if we're in the enhanced c or b sample
			enhanced=0;
			if(vTSVMCor->size()>0) {
				if(vTSVMCor->at(0) < 2000. && vTSVN->at(0)<3.) enhanced=4;
				if(vTSVMCor->at(0) > 2000. && vTSVN->at(0)>2.) enhanced=5;
			}

			//get efficiency corrections
			double effacc4(0.), effacc5(0.), effrec4(0.), effrec5(0.), effcor(1.), effsel4(0.), effsel5(0.), effpid(0.);

			if(DPT>=100000.) DPT=99999.;
			effacc4 = effacc5 =      hacc ->GetBinContent(hacc ->FindBin(DPT/1000.,DEta));
			effrec4 = effrec5 =      hrec ->GetBinContent(hrec ->FindBin(DPT/1000.,DEta));
			if(_useSimpleEffs) {
				effacc4=      hacc4->GetBinContent(hacc4->FindBin(DPT      ,DEta));
				effacc5=      hacc5->GetBinContent(hacc5->FindBin(DPT      ,DEta));
				effrec4=      hrec4->GetBinContent(hrec4->FindBin(DPT      ,DEta));
				effrec5=      hrec5->GetBinContent(hrec5->FindBin(DPT      ,DEta));
			}
			//TODO//if(_useRhoZCorr) effcor = hcor ->GetBinContent(hcor ->FindBin(DRhoSq   ,DZ  ));
			effsel4=      hsel4->GetBinContent(hsel4->FindBin(DPT      ,DEta));
			effsel5=      hsel5->GetBinContent(hsel5->FindBin(DPT      ,DEta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(DPT      ,DEta));
//TODO//20200420			if(isMC || D0P0.Pt()>=25000. || D0P0.Mag()>=500000.) effpid=1.0; //no PID cut used in these cases
			if((_whichD=="D0" || _whichD=="D" || _whichD=="K3PI") && (DP0.Pt()>=25000. || DP0.Mag()>=500000.)) effpid=1.0; //no PID cut used in these cases

			double eff4=effacc4*effrec4*effcor*effsel4*effpid;
			double eff5=effacc5*effrec5*effcor*effsel5*effpid;

			++total;
			if(eff4<0.01) {
				++skipEff4;
			}
			if(eff5<0.01) {
				++skipEff5;
			}
			if(eff4<0.01 && eff5<0.01) {
				++skipEff45;
			}
			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << DPT << "\t" << DEta << "\t" << effacc << "\t" << effrec << "\t" << effcor << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
//TODO//				if(effcor==0) std::cout << DRhoSq << "\t" << DZ << std::endl;//TODO
				continue;
			}

			int binAR = heffAR.FindFixBin(DPT/1000.,DEta);
			int binSP4= heffSP4.FindFixBin(DPT,DEta);
			int binSP5= heffSP5.FindFixBin(DPT,DEta);
			int binSimple4= heffSimple4.FindFixBin(DPT,DEta);
			int binSimple5= heffSimple5.FindFixBin(DPT,DEta);
			double effAR = heffAR.GetEfficiency(binAR);
			double effSP4= heffSP4.GetEfficiency(binSP4);
			double effSP5= heffSP5.GetEfficiency(binSP5);
			double effSimple4= heffSimple4.GetEfficiency(binSP4);
			double effSimple5= heffSimple5.GetEfficiency(binSP5);
			//symmetrise errors
			double errAR = (heffAR.GetEfficiencyErrorUp(binAR)   + heffAR.GetEfficiencyErrorLow(binAR) )/2.;
			double errSP4= (heffSP4.GetEfficiencyErrorUp(binSP4) + heffSP4.GetEfficiencyErrorLow(binSP4))/2.;
			double errSP5= (heffSP5.GetEfficiencyErrorUp(binSP5) + heffSP5.GetEfficiencyErrorLow(binSP5))/2.;
			double errSimple4= (heffSimple4.GetEfficiencyErrorUp(binSimple4) + heffSimple4.GetEfficiencyErrorLow(binSimple4))/2.;
			double errSimple5= (heffSimple5.GetEfficiencyErrorUp(binSimple5) + heffSimple5.GetEfficiencyErrorLow(binSimple5))/2.;

			double fullEff4 = effAR*effSP4*effpid;
			double fullEff5 = effAR*effSP5*effpid;
			double fullErr4 = effAR*effSP4*effpid*TMath::Sqrt(errAR*errAR/effAR/effAR + errSP4*errSP4/effSP4/effSP4);
			double fullErr5 = effAR*effSP5*effpid*TMath::Sqrt(errAR*errAR/effAR/effAR + errSP5*errSP5/effSP5/effSP5);

			if(_useSimpleEffs) {
				fullEff4 = effSimple4*effpid;
				fullEff5 = effSimple5*effpid;
				fullErr4 = errSimple4*effpid;
				fullErr5 = errSimple5*effpid;
			}

			if(TMath::Abs(fullEff4-eff4)/fullEff4 > 0.01) {
				std::cout << "WARNING in SimDFitter::addEffs : prompt efficiency cross-check differs by more than 1%" << std::endl;
				std::cout << "                               : pT=" << DPT << ", eta=" << DEta << " effs:" << fullEff4 << " " << eff4 << std::endl;
				//std::cout << effacc*effrec << " " << effsel4 << std::endl;
				if(_useSimpleEffs) {
					std::cout << effacc4*effrec4*effsel4 << std::endl;
					std::cout << effSimple4 << std::endl;
				} else {
					std::cout << effacc4 << " " << effrec4*effsel4 << std::endl;
					std::cout << effAR << " " << effSP4 << std::endl;
				}
			}
			if(TMath::Abs(fullEff5-eff5)/fullEff5 > 0.01) {
				std::cout << "WARNING in SimDFitter::addEffs : disp efficiency cross-check differs by more than 1%" << std::endl;
				std::cout << "                               : pT=" << DPT << ", eta=" << DEta << " effs:" << fullEff5 << " " << eff5 << std::endl;
				//std::cout << effacc*effrec << " " << effsel5 << std::endl;
				if(_useSimpleEffs) {
					std::cout << effacc5*effrec5*effsel5 << std::endl;
					std::cout << effSimple5 << std::endl;
				} else {
					std::cout << effacc5 << " " << effrec5*effsel5 << std::endl;
					std::cout << effAR << " " << effSP5 << std::endl;
				}
			}

			if(_usePtBinCorrections) {
				int corrBin(0);
				if(JetPT>50e3) corrBin=3;
				else if(JetPT>30e3) corrBin=2;
				else if(JetPT>20e3) corrBin=1;
				else if(JetPT>15e3) corrBin=0;
				
				double corrDPtBin = ptBinCorrections4[corrBin]->FindBin(DPT);
				double effPtBinCorr4 = ptBinCorrections4[corrBin]->GetBinContent(corrDPtBin);
				double effPtBinCorr5 = ptBinCorrections5[corrBin]->GetBinContent(corrDPtBin);
				fullEff4 *= effPtBinCorr4;
				fullErr4 *= effPtBinCorr4;
				fullEff5 *= effPtBinCorr5;
				fullErr5 *= effPtBinCorr5;
			}

			weight4 = 1./fullEff4;
			weight5 = 1./fullEff5;
			werror4 = fullErr4/fullEff4/fullEff4;
			werror5 = fullErr5/fullEff5/fullEff5;

			//term-by-term weights and errors needed for overall uncertainty
			weightAcc = 1./effAR;
			werrorAcc = errAR/effAR/effAR;
			wbinAcc = binAR;

			weightSel4 = 1./effSP4/effpid;
			werrorSel4 = weightSel4*errSP4/effSP4;
			wbinSel4 = binSP4;

			weightSel5 = 1./effSP5/effpid;
			werrorSel5 = weightSel5*errSP5/effSP5;
			wbinSel5 = binSP5;

			if(isMC) {
				//put the PID efficiency into the weights for MC
				//TODO//testWithPIDWeightsOffForTrueD0Sample//if(D0P0.Pt()<25000. && D0P0.Mag()<500000.) {
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight4 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight5 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//}
				//currently no pion PID
				//if(D0P1.Pt()<25000. && D0P1.Mag()<500000.) {
				//	weight4 *= vD0PIWEIGHT->at(s);
				//	weight5 *= vD0PIWEIGHT->at(s);
				//}
				//scale weights for roughly continuous jet true pT
				//TODO//if(JetTruePT>50000.) {//TODO added weights back in
				//TODO//	weight4*=0.03;
				//TODO//	weight5*=0.02;
				//TODO//} else if(JetTruePT>20000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.4;
				//TODO//} else if(JetTruePT>15000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.5;
				//TODO//}
				//scale weights for roughly continuous jet true pT
				//TODO//breaks eff if comparing against unweighted MC//if(JetTruePT>50000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>20000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>15000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//}
			}       //TODO//breaks eff if comparing against unweighted MC//
                                //TODO//breaks eff if comparing against unweighted MC//
			tout->Fill();
			//TODO//20200420-keepAllD0s//break;//only keep one D0 candidate per entry
		}
	}

	std::cout << "SKIPPED " << skipEff4 << " " << skipEff5 << " " << skipEff45 << " " << total << std::endl;

	dm->reset();

	tout->AutoSave();
	fout->Close();

	return true;
}

bool SimDFitter::addEffsPtBinned(int flavour) {

	std::cout << "INFO : adding " << _whichD << " efficiencies to file " << _name << ": " << flavour << std::endl;

	TFile* f2 = TFile::Open("../efficiencies/Effs_ptbinned_16x8bins_up210319.root");

	bool isMC = _dataIsMC;

	DatasetManager* dm = DatasetManager::getInstance();
	switch(flavour) {
		case 0:
			dm->setDataset(_combFile);
			break;
		case 4:
			dm->setDataset(_charmFile);
			isMC=true;
			break;
		case 5:
			dm->setDataset(_beautyFile);
			isMC=true;
			break;
		case 7:
			dm->setDataset(_dataFile);
			break;
		default:
			return false;
	}

	//construct efficiency objects
	std::vector<TH2D*> hnum4, hnum5, hden4, hden5, heff4, heff5;
	for(uint i=1; i<=4; ++i) {
		TString formatStr = "%s"+_whichD+"%d_%d";
		hnum4.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"pid",4,i))));
		hnum5.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"pid",5,i))));
		hden4.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"denom",4,i))));
		hden5.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"denom",5,i))));
		if(!hnum4.back() || !hden4.back()) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnum4.back() <<" "<< hden4.back() << std::endl;
			std::cout << "4 histograms missing" << std::endl;
			return false;
		}
		if(!hnum5.back() || !hden5.back()) {
			std::cout << "!" << _whichD << std::endl;
			f2->ls();
			std::cout << hnum5.back() <<" "<< hden5.back() << std::endl;
			std::cout << "5 histograms missing" << std::endl;
			return false;
		}

		heff4.push_back(static_cast<TH2D*>(hnum4.back()->Clone(TString::Format(formatStr,"eff",4,i))));
		heff5.push_back(static_cast<TH2D*>(hnum5.back()->Clone(TString::Format(formatStr,"eff",5,i))));
		heff4.back()->Divide(hden4.back());
		heff5.back()->Divide(hden5.back());
	}

	std::vector<double> *vDM = new std::vector<double>();
	std::vector<double> *vDIPCHI2 = new std::vector<double>();
	std::vector<double> *vDPT = new std::vector<double>();
	std::vector<double> *vDPX = new std::vector<double>();
	std::vector<double> *vDPY = new std::vector<double>();
	std::vector<double> *vDPZ = new std::vector<double>();
	std::vector<double> *vDE = new std::vector<double>();
	std::vector<double> *vDX = new std::vector<double>();
	std::vector<double> *vDY = new std::vector<double>();
	std::vector<double> *vDZ = new std::vector<double>();
	std::vector<double> *vDh0P = new std::vector<double>();
	std::vector<double> *vDh0PT = new std::vector<double>();
	std::vector<double> *vDh0PX = new std::vector<double>();
	std::vector<double> *vDh0PY = new std::vector<double>();
	std::vector<double> *vDh0PZ = new std::vector<double>();
	std::vector<double> *vDh1P = new std::vector<double>();
	std::vector<double> *vDh1PT = new std::vector<double>();
	std::vector<double> *vDh1PX = new std::vector<double>();
	std::vector<double> *vDh1PY = new std::vector<double>();
	std::vector<double> *vDh1PZ = new std::vector<double>();
	std::vector<double> *vDh2P = new std::vector<double>();
	std::vector<double> *vDh2PT = new std::vector<double>();
	std::vector<double> *vDh2PX = new std::vector<double>();
	std::vector<double> *vDh2PY = new std::vector<double>();
	std::vector<double> *vDh2PZ = new std::vector<double>();
	std::vector<double> *vDh3P = new std::vector<double>();
	std::vector<double> *vDh3PT = new std::vector<double>();
	std::vector<double> *vDh3PX = new std::vector<double>();
	std::vector<double> *vDh3PY = new std::vector<double>();
	std::vector<double> *vDh3PZ = new std::vector<double>();
	std::vector<double> *vDh0PNN = new std::vector<double>();
	std::vector<double> *vDh1PNN = new std::vector<double>();
	std::vector<double> *vDh2PNN = new std::vector<double>();
	std::vector<double> *vDh3PNN = new std::vector<double>();
	std::vector<double> *vDh0WEIGHT = new std::vector<double>();
	std::vector<double> *vDh1WEIGHT = new std::vector<double>();
	std::vector<double> *vDh2WEIGHT = new std::vector<double>();
	std::vector<double> *vDh3WEIGHT = new std::vector<double>();
	std::vector<double> *vDTRUEIDX = new std::vector<double>();
	std::vector<double> *vDFROMB = new std::vector<double>();

	std::vector<double> *vTSVMCor = new std::vector<double>();
	std::vector<double> *vTSVN = new std::vector<double>();

	double JetPT;
	double JetEta;
	double JetTruePT;

	dm->setBranchAddress(""+_whichD+"M",           &vDM);
	dm->setBranchAddress(""+_whichD+"IPCHI2",      &vDIPCHI2);
	dm->setBranchAddress(""+_whichD+"PT",          &vDPT);
	dm->setBranchAddress(""+_whichD+"PX",          &vDPX);
	dm->setBranchAddress(""+_whichD+"PY",          &vDPY);
	dm->setBranchAddress(""+_whichD+"PZ",          &vDPZ);
	dm->setBranchAddress(""+_whichD+"E",           &vDE);
	dm->setBranchAddress(""+_whichD+"X",           &vDX);
	dm->setBranchAddress(""+_whichD+"Y",           &vDY);
	dm->setBranchAddress(""+_whichD+"Z",           &vDZ);
	if(_whichD=="D0") {
		dm->setBranchAddress("D0KP",          &vDh0P);
		dm->setBranchAddress("D0KPT",         &vDh0PT);
		dm->setBranchAddress("D0KPX",         &vDh0PX);
		dm->setBranchAddress("D0KPY",         &vDh0PY);
		dm->setBranchAddress("D0KPZ",         &vDh0PZ);
		dm->setBranchAddress("D0PIP",         &vDh1P);
		dm->setBranchAddress("D0PIPT",        &vDh1PT);
		dm->setBranchAddress("D0PIPX",        &vDh1PX);
		dm->setBranchAddress("D0PIPY",        &vDh1PY);
		dm->setBranchAddress("D0PIPZ",        &vDh1PZ);
		dm->setBranchAddress("D0KPNNK",       &vDh0PNN);
		dm->setBranchAddress("D0PIPNNPI",     &vDh1PNN);
		dm->setBranchAddress("D0KWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("D0PIWEIGHT",    &vDh1WEIGHT);
	} else if(_whichD=="D") {
		dm->setBranchAddress("DKP",          &vDh0P);
		dm->setBranchAddress("DKPT",         &vDh0PT);
		dm->setBranchAddress("DKPX",         &vDh0PX);
		dm->setBranchAddress("DKPY",         &vDh0PY);
		dm->setBranchAddress("DKPZ",         &vDh0PZ);
		dm->setBranchAddress("DPI1P",        &vDh1P);
		dm->setBranchAddress("DPI1PT",       &vDh1PT);
		dm->setBranchAddress("DPI1PX",       &vDh1PX);
		dm->setBranchAddress("DPI1PY",       &vDh1PY);
		dm->setBranchAddress("DPI1PZ",       &vDh1PZ);
		dm->setBranchAddress("DPI2P",        &vDh2P);
		dm->setBranchAddress("DPI2PT",       &vDh2PT);
		dm->setBranchAddress("DPI2PX",       &vDh2PX);
		dm->setBranchAddress("DPI2PY",       &vDh2PY);
		dm->setBranchAddress("DPI2PZ",       &vDh2PZ);
		dm->setBranchAddress("DKPNNK",       &vDh0PNN);
		dm->setBranchAddress("DPI1PNNPI",    &vDh1PNN);
		dm->setBranchAddress("DPI2PNNPI",    &vDh2PNN);
		dm->setBranchAddress("DKWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("DPI1WEIGHT",   &vDh1WEIGHT);
		dm->setBranchAddress("DPI2WEIGHT",   &vDh2WEIGHT);
	} else if(_whichD=="K3PI") {
		dm->setBranchAddress("K3PIKP",          &vDh0P);
		dm->setBranchAddress("K3PIKPT",         &vDh0PT);
		dm->setBranchAddress("K3PIKPX",         &vDh0PX);
		dm->setBranchAddress("K3PIKPY",         &vDh0PY);
		dm->setBranchAddress("K3PIKPZ",         &vDh0PZ);
		dm->setBranchAddress("K3PIPI1P",        &vDh1P);
		dm->setBranchAddress("K3PIPI1PT",       &vDh1PT);
		dm->setBranchAddress("K3PIPI1PX",       &vDh1PX);
		dm->setBranchAddress("K3PIPI1PY",       &vDh1PY);
		dm->setBranchAddress("K3PIPI1PZ",       &vDh1PZ);
		dm->setBranchAddress("K3PIPI2P",        &vDh2P);
		dm->setBranchAddress("K3PIPI2PT",       &vDh2PT);
		dm->setBranchAddress("K3PIPI2PX",       &vDh2PX);
		dm->setBranchAddress("K3PIPI2PY",       &vDh2PY);
		dm->setBranchAddress("K3PIPI2PZ",       &vDh2PZ);
		dm->setBranchAddress("K3PIPI3P",        &vDh3P);
		dm->setBranchAddress("K3PIPI3PT",       &vDh3PT);
		dm->setBranchAddress("K3PIPI3PX",       &vDh3PX);
		dm->setBranchAddress("K3PIPI3PY",       &vDh3PY);
		dm->setBranchAddress("K3PIPI3PZ",       &vDh3PZ);
		dm->setBranchAddress("K3PIKPNNK",       &vDh0PNN);
		dm->setBranchAddress("K3PIPI1PNNPI",    &vDh1PNN);
		dm->setBranchAddress("K3PIPI2PNNPI",    &vDh2PNN);
		dm->setBranchAddress("K3PIPI3PNNPI",    &vDh3PNN);
		dm->setBranchAddress("K3PIKWEIGHT",     &vDh0WEIGHT);
		dm->setBranchAddress("K3PIPI1WEIGHT",   &vDh1WEIGHT);
		dm->setBranchAddress("K3PIPI2WEIGHT",   &vDh2WEIGHT);
		dm->setBranchAddress("K3PIPI3WEIGHT",   &vDh3WEIGHT);
	}

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	dm->setBranchAddress("TSVMCor",       &vTSVMCor);
	dm->setBranchAddress("TSVN",          &vTSVN);
	if(isMC) {
		dm->setBranchAddress("JetTruePT",     &JetTruePT);
		dm->setBranchAddress(""+_whichD+"TRUEIDX",     &vDTRUEIDX);
		dm->setBranchAddress(""+_whichD+"FROMB",       &vDFROMB);
	}

	unsigned int nentries0 = dm->getEntries();

	TFile* fout = TFile::Open(dFileName(),"UPDATE");
	TString tname = "T";
	tname+=flavour;
	TTree* tout = new TTree(tname,"");

	double DM(0.), DPT(0.), DEta(0.), DLogIPChi2(0.), weight4(0.), weight5(0.), werror4(0.), werror5(0.);
	double weightAcc(0.), werrorAcc(0.), wbinAcc(0.);
	double weightSel4(0.), werrorSel4(0.), wbinSel4(0.);
	double weightSel5(0.), werrorSel5(0.), wbinSel5(0.);
	double enhanced(0.);
	//TODO//double DRhoSq(0.), DZ(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch(""+_whichD+"M",          &DM);
	tout->Branch(""+_whichD+"PT",         &DPT);
	tout->Branch(""+_whichD+"Eta",        &DEta);
	tout->Branch(""+_whichD+"LogIPChi2",  &DLogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);
	tout->Branch("werror4",      &werror4);
	tout->Branch("werror5",      &werror5);
	tout->Branch("weightAcc",    &weightAcc);
	tout->Branch("weightSel4",   &weightSel4);
	tout->Branch("weightSel5",   &weightSel5);
	tout->Branch("werrorAcc",    &werrorAcc);
	tout->Branch("werrorSel4",   &werrorSel4);
	tout->Branch("werrorSel5",   &werrorSel5);
	tout->Branch("wbinAcc",      &wbinAcc);
	tout->Branch("wbinSel4",     &wbinSel4);
	tout->Branch("wbinSel5",     &wbinSel5);
	tout->Branch("enhanced",     &enhanced);

	tout->Branch("JetTruePT", &JetTruePT);

	uint total(0), skipEff4(0), skipEff5(0), skipEff45(0);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	//for(unsigned int ientry=0; ientry<nentries0; ++ientry)
	while(dm->getNext()) {
		++progress;
		//t0->GetEntry(ientry);

		for(unsigned int s=0; s<vDM->size(); ++s) {
			//first check in our tight acceptance
			TVector3 DP (vDPX->at(s)  ,vDPY->at(s)  ,vDPZ->at(s));
			TVector3 DP0(vDh0PX->at(s) ,vDh0PY->at(s) ,vDh0PZ->at(s));
			TVector3 DP1(vDh1PX->at(s),vDh1PY->at(s),vDh1PZ->at(s));

			if(flavour==7 && isMC && _truthMatchData!=truthMatchType::truthMatchOff) {
				if(_truthMatchData==truthMatchType::truthMatchPrompt) {
					if(vDTRUEIDX->at(s)==-1 || vDFROMB->at(s)==1) continue;
				}
				if(_truthMatchData==truthMatchType::truthMatchDispl) {
					if(vDTRUEIDX->at(s)==-1 || vDFROMB->at(s)==0) continue;
				}
			}

			if(!(DP0.Eta()>2. && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
			if(!(DP1.Eta()>2. && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
			if(!(DP.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
			if((flavour==4 || flavour==5) && vDTRUEIDX->at(s)==-1) continue;
			if(flavour==5 && vDFROMB->at(s)==0) continue;
			if(flavour==4 && vDFROMB->at(s)==1) continue;
			if(flavour==0 && vDM->at(s)>_mSBLow && vDM->at(s)<_mSBHigh) continue;

			//check PID cuts
//TODO//20200420			if(!isMC && vDh0PNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			if(_whichD=="D0" || _whichD=="D") {
				if(vDh0PNN->at(s)<0.2 && DP0.Pt()<25000. && DP0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
				//if(!dataIsMC && vDh1PNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off
			}

			DM         = vDM->at(s);
			DLogIPChi2 = TMath::Log(vDIPCHI2->at(s));
			DPT        = DP.Pt();
			DEta       = DP.Eta();
			//TODO//DRhoSq     = vDX->at(s)*vDX->at(s) + vDY->at(s)*vDY->at(s);
			//TODO//DZ         = vDZ->at(s);

			//for data, check if we're in the enhanced c or b sample
			enhanced=0;
			if(vTSVMCor->size()>0) {
				if(vTSVMCor->at(0) < 2000. && vTSVN->at(0)<3.) enhanced=4;
				if(vTSVMCor->at(0) > 2000. && vTSVN->at(0)>2.) enhanced=5;
			}

			//get efficiencies
			uint jptbin(0);
			if(JetPT>50e3) jptbin=3;
			else if(JetPT>30e3) jptbin=2;
			else if(JetPT>20e3) jptbin=1;
			else if(JetPT>15e3) jptbin=0;
			//else continue;
			if(DPT>=100000.) DPT=99999.;
			double eff4 = heff4[jptbin]->GetBinContent(heff4[jptbin]->FindBin(DPT,DEta));
			double eff5 = heff5[jptbin]->GetBinContent(heff5[jptbin]->FindBin(DPT,DEta));
			double err4 = heff4[jptbin]->GetBinError(heff4[jptbin]->FindBin(DPT,DEta));
			double err5 = heff5[jptbin]->GetBinError(heff5[jptbin]->FindBin(DPT,DEta));

			++total;
			if(eff4<0.01) {
				++skipEff4;
			}
			if(eff5<0.01) {
				++skipEff5;
			}
			if(eff4<0.01 && eff5<0.01) {
				++skipEff45;
			}
			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << DPT << "\t" << DEta << "\t" << effacc << "\t" << effrec << "\t" << effcor << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
//TODO//				if(effcor==0) std::cout << DRhoSq << "\t" << DZ << std::endl;//TODO
				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;
			werror4 = err4/eff4/eff4;
			werror5 = err5/eff5/eff5;

			if(isMC) {
				//put the PID efficiency into the weights for MC
				//TODO//testWithPIDWeightsOffForTrueD0Sample//if(D0P0.Pt()<25000. && D0P0.Mag()<500000.) {
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight4 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight5 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//}
				//currently no pion PID
				//if(D0P1.Pt()<25000. && D0P1.Mag()<500000.) {
				//	weight4 *= vD0PIWEIGHT->at(s);
				//	weight5 *= vD0PIWEIGHT->at(s);
				//}
				//scale weights for roughly continuous jet true pT
				//TODO//if(JetTruePT>50000.) {//TODO added weights back in
				//TODO//	weight4*=0.03;
				//TODO//	weight5*=0.02;
				//TODO//} else if(JetTruePT>20000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.4;
				//TODO//} else if(JetTruePT>15000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.5;
				//TODO//}
				//scale weights for roughly continuous jet true pT
				//TODO//breaks eff if comparing against unweighted MC//if(JetTruePT>50000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>20000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>15000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//}
			}       //TODO//breaks eff if comparing against unweighted MC//
                                //TODO//breaks eff if comparing against unweighted MC//
			tout->Fill();
			//TODO//20200420-keepAllD0s//break;//only keep one D0 candidate per entry
		}
	}

	std::cout << "SKIPPED " << skipEff4 << " " << skipEff5 << " " << skipEff45 << " " << total << std::endl;

	dm->reset();

	tout->AutoSave();
	fout->Close();

	return true;
}

bool SimDFitter::testEffs(int flavour) {
	if(_useJetPtBinnedEffs) return testEffsPtBinned(flavour);
	if(!_dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << _name << std::endl;
	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset(_dataFile);

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);
	TFile* f4(0);
	if(flavour==5) f4 = TFile::Open("../efficiencies/"+_whichD+"5pTEffCorrections.root");
	else f4 = TFile::Open("../efficiencies/"+_whichD+"4pTEffCorrections.root");

	//TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pideff"+_whichD+"45"));
	TH2D *hacc(0), *hrec(0), *hsel(0);
	if(_flavour==5) hsel = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"5"));
	else hsel = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"4"));

	double unitsScale(1.);
	if(!_useSimpleEffs) {
		hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
		hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
		unitsScale=1e-3;
	} else {
		hacc  = dynamic_cast<TH2D*>(f2->Get("acceff"+_whichD+"45"));
		hrec  = dynamic_cast<TH2D*>(f2->Get("receff"+_whichD+"45"));
	}

	std::vector<TH1D*> ptBinCorrections;
	ptBinCorrections.push_back(static_cast<TH1D*>(f4->Get("corrPtBin1")));
	ptBinCorrections.push_back(static_cast<TH1D*>(f4->Get("corrPtBin2")));
	ptBinCorrections.push_back(static_cast<TH1D*>(f4->Get("corrPtBin3")));
	ptBinCorrections.push_back(static_cast<TH1D*>(f4->Get("corrPtBin4")));

	if(!hacc || !hrec /*|| !hcor*/ || !hpid || !hsel) return false;

	TH2D* effModel = cloneTH2D("effModel", hpid);
	TH2D* denData  = cloneTH2D("denData", hpid);
	TH2D* numData  = cloneTH2D("numData", hpid);
	TH2D* effData  = cloneTH2D("effData", hpid);
	std::vector<TH2D*> denDataJPTBins;
	std::vector<TH2D*> numDataJPTBins;
	std::vector<TH2D*> effDataJPTBins;
	std::vector<TH2D*> denDataJTRUEPTBins;
	std::vector<TH2D*> numDataJTRUEPTBins;
	std::vector<TH2D*> effDataJTRUEPTBins;
	for(uint i=1; i<_ptBins.size(); ++i) {
		denDataJPTBins.push_back(cloneTH2D(TString::Format("denDataJPTBin%d",i), hpid));
		numDataJPTBins.push_back(cloneTH2D(TString::Format("numDataJPTBin%d",i), hpid));
		effDataJPTBins.push_back(cloneTH2D(TString::Format("effDataJPTBin%d",i), hpid));
		denDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("denDataJTRUEPTBin%d",i), hpid));
		numDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("numDataJTRUEPTBin%d",i), hpid));
		effDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("effDataJTRUEPTBin%d",i), hpid));
	}

	//add over/underflow to the pT binning
	std::vector<double> ptBinsWithOF = {0.};
	for(uint i=0; i<_ptBins.size(); ++i) ptBinsWithOF.push_back(_ptBins.at(i));
	ptBinsWithOF.push_back(110e3);
	for(uint i=0; i<ptBinsWithOF.size(); ++i) std::cout << ptBinsWithOF.at(i) << std::endl;//TODO

	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
	std::vector<TString> stages;
	stages.push_back("all");
	stages.push_back("geometric");
	stages.push_back("reconstruction");
	stages.push_back("selection");
	stages.push_back("PID");
	std::vector<TH1D*> truDtruePT;//the true distribution at this stage
	truDtruePT.push_back(cloneTH1D("tru0truept", ptBinsWithOF));
	truDtruePT.push_back(cloneTH1D("tru1truept", ptBinsWithOF));
	truDtruePT.push_back(cloneTH1D("tru2truept", ptBinsWithOF));
	truDtruePT.push_back(cloneTH1D("tru3truept", ptBinsWithOF));
	truDtruePT.push_back(cloneTH1D("tru4truept", ptBinsWithOF));
	std::vector<TH1D*> effDtruePT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effDtruePT.push_back(cloneTH1D("eff0truept", ptBinsWithOF));
	effDtruePT.push_back(cloneTH1D("eff1truept", ptBinsWithOF));
	effDtruePT.push_back(cloneTH1D("eff2truept", ptBinsWithOF));
	effDtruePT.push_back(cloneTH1D("eff3truept", ptBinsWithOF));
	effDtruePT.push_back(cloneTH1D("eff4truept", ptBinsWithOF));
	std::vector<TH1D*> fitDtruePT;//sideband-subtracted distribution at this stage
	fitDtruePT.push_back(cloneTH1D("fit0truept", ptBinsWithOF));
	fitDtruePT.push_back(cloneTH1D("fit1truept", ptBinsWithOF));
	fitDtruePT.push_back(cloneTH1D("fit2truept", ptBinsWithOF));
	fitDtruePT.push_back(cloneTH1D("fit3truept", ptBinsWithOF));
	fitDtruePT.push_back(cloneTH1D("fit4truept", ptBinsWithOF));
	std::vector<TH1D*> fndDtruePT;//truth-matched distribution at this stage
	fndDtruePT.push_back(cloneTH1D("fnd0truept", ptBinsWithOF));
	fndDtruePT.push_back(cloneTH1D("fnd1truept", ptBinsWithOF));
	fndDtruePT.push_back(cloneTH1D("fnd2truept", ptBinsWithOF));
	fndDtruePT.push_back(cloneTH1D("fnd3truept", ptBinsWithOF));
	fndDtruePT.push_back(cloneTH1D("fnd4truept", ptBinsWithOF));
	//binned in true jet pT
	std::vector<TH1D*> truDrecoPT;
	truDrecoPT.push_back(cloneTH1D("tru0recopt", ptBinsWithOF));
	truDrecoPT.push_back(cloneTH1D("tru1recopt", ptBinsWithOF));
	truDrecoPT.push_back(cloneTH1D("tru2recopt", ptBinsWithOF));
	truDrecoPT.push_back(cloneTH1D("tru3recopt", ptBinsWithOF));
	truDrecoPT.push_back(cloneTH1D("tru4recopt", ptBinsWithOF));
	std::vector<TH1D*> effDrecoPT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effDrecoPT.push_back(cloneTH1D("eff0recopt", ptBinsWithOF));
	effDrecoPT.push_back(cloneTH1D("eff1recopt", ptBinsWithOF));
	effDrecoPT.push_back(cloneTH1D("eff2recopt", ptBinsWithOF));
	effDrecoPT.push_back(cloneTH1D("eff3recopt", ptBinsWithOF));
	effDrecoPT.push_back(cloneTH1D("eff4recopt", ptBinsWithOF));
	std::vector<TH1D*> fitDrecoPT;
	fitDrecoPT.push_back(cloneTH1D("fit0recopt", ptBinsWithOF));
	fitDrecoPT.push_back(cloneTH1D("fit1recopt", ptBinsWithOF));
	fitDrecoPT.push_back(cloneTH1D("fit2recopt", ptBinsWithOF));
	fitDrecoPT.push_back(cloneTH1D("fit3recopt", ptBinsWithOF));
	fitDrecoPT.push_back(cloneTH1D("fit4recopt", ptBinsWithOF));
	std::vector<TH1D*> fndDrecoPT;
	fndDrecoPT.push_back(cloneTH1D("fnd0recopt", ptBinsWithOF));
	fndDrecoPT.push_back(cloneTH1D("fnd1recopt", ptBinsWithOF));
	fndDrecoPT.push_back(cloneTH1D("fnd2recopt", ptBinsWithOF));
	fndDrecoPT.push_back(cloneTH1D("fnd3recopt", ptBinsWithOF));
	fndDrecoPT.push_back(cloneTH1D("fnd4recopt", ptBinsWithOF));
	//binned in both reco and true jet pT
	std::vector<TH2D*> truD2DPT;
	truD2DPT.push_back(cloneTH2D("tru02Dpt", ptBinsWithOF, ptBinsWithOF));
	truD2DPT.push_back(cloneTH2D("tru12Dpt", ptBinsWithOF, ptBinsWithOF));
	truD2DPT.push_back(cloneTH2D("tru22Dpt", ptBinsWithOF, ptBinsWithOF));
	truD2DPT.push_back(cloneTH2D("tru32Dpt", ptBinsWithOF, ptBinsWithOF));
	truD2DPT.push_back(cloneTH2D("tru42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> effD2DPT;
	effD2DPT.push_back(cloneTH2D("eff02Dpt", ptBinsWithOF, ptBinsWithOF));
	effD2DPT.push_back(cloneTH2D("eff12Dpt", ptBinsWithOF, ptBinsWithOF));
	effD2DPT.push_back(cloneTH2D("eff22Dpt", ptBinsWithOF, ptBinsWithOF));
	effD2DPT.push_back(cloneTH2D("eff32Dpt", ptBinsWithOF, ptBinsWithOF));
	effD2DPT.push_back(cloneTH2D("eff42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> fitD2DPT;
	fitD2DPT.push_back(cloneTH2D("fit02Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD2DPT.push_back(cloneTH2D("fit12Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD2DPT.push_back(cloneTH2D("fit22Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD2DPT.push_back(cloneTH2D("fit32Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD2DPT.push_back(cloneTH2D("fit42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> fndD2DPT;
	fndD2DPT.push_back(cloneTH2D("fnd02Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD2DPT.push_back(cloneTH2D("fnd12Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD2DPT.push_back(cloneTH2D("fnd22Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD2DPT.push_back(cloneTH2D("fnd32Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD2DPT.push_back(cloneTH2D("fnd42Dpt", ptBinsWithOF, ptBinsWithOF));

	//entries in this vector are true and reco for each stage
	std::vector<TH2D*> dKinematics;
	dKinematics.push_back(cloneTH2D("dKineTrue0", hsel));
	dKinematics.push_back(cloneTH2D("dKineReco0", hsel));
	dKinematics.push_back(cloneTH2D("dKineTrue1", hsel));
	dKinematics.push_back(cloneTH2D("dKineReco1", hsel));
	dKinematics.push_back(cloneTH2D("dKineTrue2", hsel));
	dKinematics.push_back(cloneTH2D("dKineReco2", hsel));
	dKinematics.push_back(cloneTH2D("dKineTrue3", hsel));
	dKinematics.push_back(cloneTH2D("dKineReco3", hsel));
	dKinematics.push_back(cloneTH2D("dKineTrue4", hsel));
	dKinematics.push_back(cloneTH2D("dKineReco4", hsel));
	for(unsigned int i=0; i<dKinematics.size(); ++i) {
		dKinematics[i]->Sumw2();
	}
	std::vector<TH1D*> dKinematicsPulls;
	dKinematicsPulls.push_back(new TH1D("pulls1", "", 30, -1., 1.));
	dKinematicsPulls.push_back(new TH1D("pulls2", "", 30, -1., 1.));
	dKinematicsPulls.push_back(new TH1D("pulls3", "", 30, -1., 1.));
	dKinematicsPulls.push_back(new TH1D("pulls4", "", 30, -1., 1.));

	TH1D* closure_true = cloneTH1D("closure_true", _ptBins);
	TH1D* closure_reco = cloneTH1D("closure_reco", _ptBins);

	//also divide distribution into bins of true and reco pT
	std::vector<TH2D*> dKinematicsTRBinned;
	std::vector<TH2D*> dKinematicsTRBinned_recD_tru;
	std::vector<TH2D*> dKinematicsTRBinned_recD_rec;
	for(uint i=0; i<ptBinsWithOF.size(); ++i) {
		for(uint j=0; j<ptBinsWithOF.size()-1; ++j) {
			TString hName = TString::Format("dKineJPTBin%d%d",i,j);
			TString hTitle = TString::Format("true%d reco%d",j,i);
			if(j==ptBinsWithOF.size()-2) {
				hTitle.Form("reco%d",i);
			}
			if(i==ptBinsWithOF.size()-1) {
				hTitle.Form("true%d",j);
			}
			if(j==ptBinsWithOF.size()-2 && i==ptBinsWithOF.size()-1) {
				hTitle = "all";
			}
			dKinematicsTRBinned.push_back(cloneTH2D(hName, hsel));
			dKinematicsTRBinned.back()->Sumw2();
			dKinematicsTRBinned.back()->SetTitle(hTitle);
			dKinematicsTRBinned.back()->GetXaxis()->SetTitle("pT(D)");
			dKinematicsTRBinned.back()->GetYaxis()->SetTitle("#eta(D)");
			
			dKinematicsTRBinned_recD_tru.push_back(cloneTH2D(hName+"_recD_tru", hsel));
			dKinematicsTRBinned_recD_tru.back()->Sumw2();
			dKinematicsTRBinned_recD_tru.back()->SetTitle(hTitle);
			dKinematicsTRBinned_recD_tru.back()->GetXaxis()->SetTitle("pT(D)");
			dKinematicsTRBinned_recD_tru.back()->GetYaxis()->SetTitle("#eta(D)");
			
			dKinematicsTRBinned_recD_rec.push_back(cloneTH2D(hName+"_recD_rec", hsel));
			dKinematicsTRBinned_recD_rec.back()->Sumw2();
			dKinematicsTRBinned_recD_rec.back()->SetTitle(hTitle);
			dKinematicsTRBinned_recD_rec.back()->GetXaxis()->SetTitle("pT(D)");
			dKinematicsTRBinned_recD_rec.back()->GetYaxis()->SetTitle("#eta(D)");
		}
	}

	for(int i=0; i<5; ++i) {
		truDtruePT.at(i)->Sumw2();
		effDtruePT.at(i)->Sumw2();
		fitDtruePT.at(i)->Sumw2();
		fndDtruePT.at(i)->Sumw2();
		truDrecoPT.at(i)->Sumw2();
		effDrecoPT.at(i)->Sumw2();
		fitDrecoPT.at(i)->Sumw2();
		fndDrecoPT.at(i)->Sumw2();
		truD2DPT.at(i)->Sumw2();
		effD2DPT.at(i)->Sumw2();
		fitD2DPT.at(i)->Sumw2();
		fndD2DPT.at(i)->Sumw2();
		truD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		effD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		fitD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		fndD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		truD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		effD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		fitD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		fndD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
	}

	std::vector<double> *DM = new std::vector<double>();
	std::vector<double> *DIPCHI2 = new std::vector<double>();
	std::vector<double> *DPT = new std::vector<double>();
	std::vector<double> *DPX = new std::vector<double>();
	std::vector<double> *DPY = new std::vector<double>();
	std::vector<double> *DPZ = new std::vector<double>();
	std::vector<double> *DE = new std::vector<double>();
	std::vector<double> *DX = new std::vector<double>();
	std::vector<double> *DY = new std::vector<double>();
	std::vector<double> *DZ = new std::vector<double>();
	std::vector<double> *Dh0P = new std::vector<double>();
	std::vector<double> *Dh0PT = new std::vector<double>();
	std::vector<double> *Dh0PX = new std::vector<double>();
	std::vector<double> *Dh0PY = new std::vector<double>();
	std::vector<double> *Dh0PZ = new std::vector<double>();
	std::vector<double> *Dh1P = new std::vector<double>();
	std::vector<double> *Dh1PT = new std::vector<double>();
	std::vector<double> *Dh1PX = new std::vector<double>();
	std::vector<double> *Dh1PY = new std::vector<double>();
	std::vector<double> *Dh1PZ = new std::vector<double>();
	std::vector<double> *Dh2P = new std::vector<double>();
	std::vector<double> *Dh2PT = new std::vector<double>();
	std::vector<double> *Dh2PX = new std::vector<double>();
	std::vector<double> *Dh2PY = new std::vector<double>();
	std::vector<double> *Dh2PZ = new std::vector<double>();
	std::vector<double> *Dh0WEIGHT = new std::vector<double>();
	std::vector<double> *Dh1WEIGHT = new std::vector<double>();
	std::vector<double> *Dh2WEIGHT = new std::vector<double>();
	std::vector<double> *DTRUEIDX = new std::vector<double>();
	std::vector<double> *DTRUETRK0 = new std::vector<double>();
	std::vector<double> *DTRUETRK1 = new std::vector<double>();
	std::vector<double> *DTRUETRK2 = new std::vector<double>();

	std::vector<double> *TRUEDID = new std::vector<double>();
	std::vector<double> *TRUEDPX = new std::vector<double>();
	std::vector<double> *TRUEDPY = new std::vector<double>();
	std::vector<double> *TRUEDPZ = new std::vector<double>();
	std::vector<double> *TRUEDE = new std::vector<double>();
	std::vector<double> *TRUEDX = new std::vector<double>();
	std::vector<double> *TRUEDY = new std::vector<double>();
	std::vector<double> *TRUEDZ = new std::vector<double>();
	std::vector<double> *TRUEDFROMB = new std::vector<double>();
	std::vector<double> *TRUEDTRK0P = new std::vector<double>();
	std::vector<double> *TRUEDTRK0PT = new std::vector<double>();
	std::vector<double> *TRUEDTRK0INACC = new std::vector<double>();
	std::vector<double> *TRUEDTRK0RECO = new std::vector<double>();
	std::vector<double> *TRUEDTRK1P = new std::vector<double>();
	std::vector<double> *TRUEDTRK1PT = new std::vector<double>();
	std::vector<double> *TRUEDTRK1INACC = new std::vector<double>();
	std::vector<double> *TRUEDTRK1RECO = new std::vector<double>();
	std::vector<double> *TRUEDTRK0IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK1IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK2IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK3IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK0ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK1ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK2ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK3ID = new std::vector<double>();

	double JetPT;
	double JetEta;
	double JetTruePT;

	dm->setBranchAddress(""+_whichD+"M",           &DM);
	dm->setBranchAddress(""+_whichD+"IPCHI2",      &DIPCHI2);
	dm->setBranchAddress(""+_whichD+"PT",          &DPT);
	dm->setBranchAddress(""+_whichD+"PX",          &DPX);
	dm->setBranchAddress(""+_whichD+"PY",          &DPY);
	dm->setBranchAddress(""+_whichD+"PZ",          &DPZ);
	dm->setBranchAddress(""+_whichD+"E",           &DE);
	dm->setBranchAddress(""+_whichD+"X",           &DX);
	dm->setBranchAddress(""+_whichD+"Y",           &DY);
	dm->setBranchAddress(""+_whichD+"Z",           &DZ);
	if(_whichD=="D0") {
		dm->setBranchAddress("D0KP",          &Dh0P);
		dm->setBranchAddress("D0KPT",         &Dh0PT);
		dm->setBranchAddress("D0KPX",         &Dh0PX);
		dm->setBranchAddress("D0KPY",         &Dh0PY);
		dm->setBranchAddress("D0KPZ",         &Dh0PZ);
		dm->setBranchAddress("D0PIP",         &Dh1P);
		dm->setBranchAddress("D0PIPT",        &Dh1PT);
		dm->setBranchAddress("D0PIPX",        &Dh1PX);
		dm->setBranchAddress("D0PIPY",        &Dh1PY);
		dm->setBranchAddress("D0PIPZ",        &Dh1PZ);
		dm->setBranchAddress("D0KWEIGHT",     &Dh0WEIGHT);
		dm->setBranchAddress("D0PIWEIGHT",    &Dh1WEIGHT);
		dm->setBranchAddress("D0TRUEIDX",     &DTRUEIDX);
		dm->setBranchAddress("D0TRUETRK0",    &DTRUETRK0);
		dm->setBranchAddress("D0TRUETRK1",    &DTRUETRK1);
	} if(_whichD=="D") {
		dm->setBranchAddress("DKP",          &Dh0P);
		dm->setBranchAddress("DKPT",         &Dh0PT);
		dm->setBranchAddress("DKPX",         &Dh0PX);
		dm->setBranchAddress("DKPY",         &Dh0PY);
		dm->setBranchAddress("DKPZ",         &Dh0PZ);
		dm->setBranchAddress("DPI1P",        &Dh1P);
		dm->setBranchAddress("DPI1PT",       &Dh1PT);
		dm->setBranchAddress("DPI1PX",       &Dh1PX);
		dm->setBranchAddress("DPI1PY",       &Dh1PY);
		dm->setBranchAddress("DPI1PZ",       &Dh1PZ);
		dm->setBranchAddress("DPI2P",        &Dh2P);
		dm->setBranchAddress("DPI2PT",       &Dh2PT);
		dm->setBranchAddress("DPI2PX",       &Dh2PX);
		dm->setBranchAddress("DPI2PY",       &Dh2PY);
		dm->setBranchAddress("DPI2PZ",       &Dh2PZ);
		dm->setBranchAddress("DKWEIGHT",     &Dh0WEIGHT);
		dm->setBranchAddress("DPI1WEIGHT",   &Dh1WEIGHT);
		dm->setBranchAddress("DPI2WEIGHT",   &Dh2WEIGHT);
		dm->setBranchAddress("DTRUEIDX",     &DTRUEIDX);
		dm->setBranchAddress("DTRUETRK0",    &DTRUETRK0);
		dm->setBranchAddress("DTRUETRK1",    &DTRUETRK1);
		dm->setBranchAddress("DTRUETRK2",    &DTRUETRK2);
	}

	dm->setBranchAddress("TRUEDID",       &TRUEDID);
	dm->setBranchAddress("TRUEDPX",       &TRUEDPX);
	dm->setBranchAddress("TRUEDPY",       &TRUEDPY);
	dm->setBranchAddress("TRUEDPZ",       &TRUEDPZ);
	dm->setBranchAddress("TRUEDE",        &TRUEDE);
	dm->setBranchAddress("TRUEDX",        &TRUEDX);
	dm->setBranchAddress("TRUEDY",        &TRUEDY);
	dm->setBranchAddress("TRUEDZ",        &TRUEDZ);
	dm->setBranchAddress("TRUEDFROMB",    &TRUEDFROMB);
	dm->setBranchAddress("TRUEDTRK0P",    &TRUEDTRK0P);
	dm->setBranchAddress("TRUEDTRK0PT",   &TRUEDTRK0PT);
	dm->setBranchAddress("TRUEDTRK0INACC",&TRUEDTRK0INACC);
	dm->setBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO);
	dm->setBranchAddress("TRUEDTRK1P",    &TRUEDTRK1P);
	dm->setBranchAddress("TRUEDTRK1PT",   &TRUEDTRK1PT);
	dm->setBranchAddress("TRUEDTRK1INACC",&TRUEDTRK1INACC);
	dm->setBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO);
	dm->setBranchAddress("TRUEDTRK0IDX",  &TRUEDTRK0IDX);
	dm->setBranchAddress("TRUEDTRK1IDX",  &TRUEDTRK1IDX);
	dm->setBranchAddress("TRUEDTRK2IDX",  &TRUEDTRK2IDX);
	dm->setBranchAddress("TRUEDTRK3IDX",  &TRUEDTRK3IDX);
	dm->setBranchAddress("TRUEDTRK0ID",  &TRUEDTRK0ID);
	dm->setBranchAddress("TRUEDTRK1ID",  &TRUEDTRK1ID);
	dm->setBranchAddress("TRUEDTRK2ID",  &TRUEDTRK2ID);
	dm->setBranchAddress("TRUEDTRK3ID",  &TRUEDTRK3ID);

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	dm->setBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = dm->getEntries();

	double dpt(0.), deta(0.);

	for(int i=1; i<=effModel->GetNbinsX(); ++i) {
		for(int j=1; j<=effModel->GetNbinsY(); ++j) {
			dpt = effModel->GetXaxis()->GetBinCenter(i);
			deta = effModel->GetYaxis()->GetBinCenter(j);
			double effacc(0.), effrec(0.), effsel(0.), effpid(0.), eff(0.);
			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));
			eff = effacc*effrec*effsel*effpid;
			//std::cout << i << " " << j << "\t" << dpt << "\t" << deta << "\t" << eff << std::endl;
			effModel->SetBinContent(i,j,eff);
		}
	}

	boost::progress_display progress( nentries0 );
	//for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
	while(dm->getNext()) {
		++progress;
		//t0->GetEntry(ientry);

		std::vector<int> trueDs, foundDs;

		int jrptbin(0), jtptbin(0);
		for(uint i=0; i<_ptBins.size(); ++i) {
			if(JetPT>_ptBins[i]) ++jrptbin;
			if(JetTruePT>_ptBins[i]) ++jtptbin;
		}

		//put anything out of range into the overflow bin
		if(JetPT>_ptBins.back()) JetPT = _ptBins.back();
		if(JetTruePT>_ptBins.back()) JetTruePT = _ptBins.back();

		for(unsigned int d=0; d<TRUEDID->size(); ++d) {
			if(_whichD=="D0") {
				//only use D0->Kpi with Pt>5GeV
				if(TMath::Abs(TRUEDID->at(d))!=421) continue;
				if(TRUEDTRK0IDX->at(d)==-1) continue;
				if(TRUEDTRK2IDX->at(d)!=-1) continue;
			} else if(_whichD=="D") {
				//only use D->Kpipi with Pt>5GeV
				if(TMath::Abs(TRUEDID->at(d))!=411) continue;
				if(TRUEDTRK0IDX->at(d)==-1) continue;
				if(TRUEDTRK2IDX->at(d)==-1) continue;
				if(TRUEDTRK3IDX->at(d)!=-1) continue;
			}
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;

			TVector3 TRUEDP (TRUEDPX->at(d)  ,TRUEDPY->at(d)  ,TRUEDPZ->at(d));

			dpt        = TRUEDP.Pt();
			deta       = TRUEDP.Eta();
			//double rhoSq = TRUEDX->at(d)*TRUEDX->at(d) + TRUEDY->at(d)*TRUEDY->at(d);
			//double z = TRUEDZ->at(d);
			//if(rhoSq>=100) rhoSq=99.9;//TODO

			//get efficiency corrections
			double effacc(0.), effrec(0.), effsel(0.), effpid(0.);

			if(dpt>=100000.) dpt=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
			//apply correction to the reco efficiency
			//if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			denData->Fill(dpt,deta);
			if(jrptbin>0 && jrptbin<=denDataJPTBins.size()) denDataJPTBins[jrptbin-1]->Fill(dpt,deta);
			if(jtptbin>0 && jtptbin<=denDataJPTBins.size()) denDataJTRUEPTBins[jtptbin-1]->Fill(dpt,deta);

			truDrecoPT.at(0)->Fill(JetPT);
			truDtruePT.at(0)->Fill(JetTruePT);
			truD2DPT.at(0)->Fill(JetTruePT,JetPT);
			dKinematics[0]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			//fill kinematics separately for each true/reco jet pT bin
			int ibin = truDrecoPT.at(0)->FindBin(JetTruePT) - 1;
			int jbin = truDrecoPT.at(0)->FindBin(JetPT) - 1;
			dKinematicsTRBinned.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
			dKinematicsTRBinned.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
			dKinematicsTRBinned.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
			dKinematicsTRBinned.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size()-1)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all

			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;

			truDrecoPT.at(1)->Fill(JetPT);
			truDtruePT.at(1)->Fill(JetTruePT);
			truD2DPT.at(1)->Fill(JetTruePT,JetPT);
			if(effacc>0.) {
				effDrecoPT.at(0)->Fill(JetPT,1./effacc);
				effDtruePT.at(0)->Fill(JetTruePT,1./effacc);
				effD2DPT.at(0)->Fill(JetTruePT,JetPT,1./effacc);
			} else {
				std::cout << "ACC=" << effacc << " at " << dpt << "," << deta << std::endl;
			}
			dKinematics[2]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			truDrecoPT.at(2)->Fill(JetPT);
			truDtruePT.at(2)->Fill(JetTruePT);
			truD2DPT.at(2)->Fill(JetTruePT,JetPT);
			if(effrec>0.) {
				effDrecoPT.at(1)->Fill(JetPT,1./effrec);
				effDtruePT.at(1)->Fill(JetTruePT,1./effrec);
				effD2DPT.at(1)->Fill(JetTruePT,JetPT,1./effrec);
			} else {
				std::cout << "REC=" << effrec << " at " << dpt << "," << deta << std::endl;
			}
			dKinematics[4]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			for(unsigned int s=0; s<DTRUEIDX->size(); ++s) {
				if(DTRUEIDX->at(s)==d) {
					TVector3 DP (DPX->at(s)  ,DPY->at(s)  ,DPZ->at(s));
					TVector3 DP0(Dh0PX->at(s) ,Dh0PY->at(s) ,Dh0PZ->at(s));
					TVector3 DP1(Dh1PX->at(s),Dh1PY->at(s),Dh1PZ->at(s));

					if(!(DP0.Eta()>2.0 && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
					if(!(DP1.Eta()>2.0 && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
					if(!(DP.Pt()>_dptmin)) continue;
					if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
					//if(!(D0P.Eta()>2.5&&D0P.Eta()<4.0)) continue;//TODO test restricting the eta(D0) range (turn on in two places 1/2)

					truDrecoPT.at(3)->Fill(JetPT);
					truDtruePT.at(3)->Fill(JetTruePT);
					truD2DPT.at(3)->Fill(JetTruePT,JetPT);
					if(effsel>0.) {
						effDrecoPT.at(2)->Fill(JetPT,1./effsel);
						effDtruePT.at(2)->Fill(JetTruePT,1./effsel);
						effD2DPT.at(2)->Fill(JetTruePT,JetPT,1./effsel);
					} else {
						std::cout << "SEL=" << effsel << " at " << dpt << "," << deta << std::endl;
					}
					dKinematics[6]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

					//use weights for PID
					double weight = Dh0WEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
					if(DP0.Pt()>25000. || DP0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
					truDrecoPT.at(4)->Fill(JetPT,weight);
					truDtruePT.at(4)->Fill(JetTruePT,weight);
					truD2DPT.at(4)->Fill(JetTruePT,JetPT);
					numData->Fill(dpt,deta,weight);
					if(jrptbin>0 && jrptbin<=numDataJPTBins.size()) numDataJPTBins[jrptbin-1]->Fill(dpt,deta,weight);
					if(jtptbin>0 && jtptbin<=numDataJPTBins.size()) numDataJTRUEPTBins[jtptbin-1]->Fill(dpt,deta,weight);
					if(effpid>0.) {
						effDrecoPT.at(3)->Fill(JetPT,weight/effpid);
						effDtruePT.at(3)->Fill(JetTruePT,weight/effpid);
						effD2DPT.at(3)->Fill(JetTruePT,JetPT,1./effpid);
					} else {
						std::cout << "PID=" << effpid << " at " << dpt << "," << deta << std::endl;
					}
					dKinematics[8]->Fill(TRUEDP.Pt(),TRUEDP.Eta(),weight);

					trueDs.push_back(d*100+s);

					//fill kinematics separately for each true/reco jet pT bin
					dKinematicsTRBinned_recD_tru.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
					dKinematicsTRBinned_recD_tru.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
					dKinematicsTRBinned_recD_tru.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
					dKinematicsTRBinned_recD_tru.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size() -1)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all

					//fill kinematics separately for each true/reco jet pT bin
					dKinematicsTRBinned_recD_rec.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());
					dKinematicsTRBinned_recD_rec.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());//row sums
					dKinematicsTRBinned_recD_rec.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());//col sums
					dKinematicsTRBinned_recD_rec.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size() -1)->Fill(DP.Pt(),DP.Eta());//all

					break;
				}
			}
		}
		for(unsigned int s=0; s<DM->size(); ++s) {
			//first check in our tight acceptance
			TVector3 DP (DPX->at(s)  ,DPY->at(s)  ,DPZ->at(s));
			TVector3 DP0(Dh0PX->at(s),Dh0PY->at(s),Dh0PZ->at(s));
			TVector3 DP1(Dh1PX->at(s),Dh1PY->at(s),Dh1PZ->at(s));
			//double rhoSq = D0X->at(s)*D0X->at(s) + D0Y->at(s)*D0Y->at(s);
			//double z = D0Z->at(s);
			//if(rhoSq>=100) rhoSq=99.9;//TODO

			if(!(DP0.Eta()>2.0 && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
			if(!(DP1.Eta()>2.0 && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
			if(!(DP.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
			//if(!(D0P.Eta()>2.2&&D0P.Eta()<4.0)) continue;//TODO test restrcting the eta(D0) range (turn on in two places 2/2)

			//check PID cuts
			double weight = Dh0WEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
			if(DP0.Pt()>25000. || DP0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT

			dpt        = DP.Pt();
			deta       = DP.Eta();

			//get efficiency corrections
			double effacc(0.), effrec(0.), effsel(0.), effpid(0.), corrPt(0.);

			if(dpt>=100000.) dpt=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
			//apply correction to the reco efficiency
			//if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			uint corrBin(0);
			if(JetPT>50e3) corrBin=3;
			else if(JetPT>30e3) corrBin=2;
			else if(JetPT>20e3) corrBin=1;
			else if(JetPT>15e3) corrBin=0;
			corrPt = ptBinCorrections[corrBin]->GetBinContent(ptBinCorrections[corrBin]->FindBin(DP.Pt()));

			double eff=effacc*effrec*corrPt*effsel*effpid;

			if(eff<0.01) {
				//TODO//std::cout << dpt << "\t" << deta << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
				//if(D0TRUEIDX->at(s)>-1) {
				//	int d = D0TRUEIDX->at(s);
				//	if(TRUEDTRK0IDX->at(d)!=-1 && TRUEDTRK2IDX->at(d)==-1 && TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)>=5000.*5000.) {
				//		if(TRUEDTRK0INACC->at(d)==1 && TRUEDTRK1INACC->at(d)==1) {
				//			if(TRUEDTRK0RECO->at(d)==1 && TRUEDTRK1RECO->at(d)==1) {
				//				std::cout << "LOST TRUE D0" << std::endl;
				//std::cout << dpt << "\t" << deta << "\t" << eff << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
				//std::cout << rhoSq << "\t" << z << "\t" << hrec ->GetBinContent(hrec ->FindBin(dpt/1000.,deta)) << "\t" << hcor ->GetBinContent(hcor ->FindBin(rhoSq,z)) << std::endl;
				//			}
				//		}
				//	}
				//}
				continue;
			}

			double sbSubSwitch(0.);
			if(TMath::Abs(DM->at(s)-_mCentre)<40.) sbSubSwitch = 1.;
			else if(TMath::Abs(DM->at(s)-_mCentre)<80.) sbSubSwitch = -1.;

			fitDrecoPT.at(4)->Fill(JetPT,weight*sbSubSwitch);
			fitDrecoPT.at(3)->Fill(JetPT,weight*sbSubSwitch/effpid);
			fitDrecoPT.at(2)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel));
			fitDrecoPT.at(1)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitDrecoPT.at(0)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			fitDtruePT.at(4)->Fill(JetTruePT,weight*sbSubSwitch);
			fitDtruePT.at(3)->Fill(JetTruePT,weight*sbSubSwitch/effpid);
			fitDtruePT.at(2)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel));
			fitDtruePT.at(1)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitDtruePT.at(0)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			fitD2DPT.at(4)->Fill(JetTruePT,JetPT,weight*sbSubSwitch);
			fitD2DPT.at(3)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/effpid);
			fitD2DPT.at(2)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel));
			fitD2DPT.at(1)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitD2DPT.at(0)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			dKinematics[9]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch);
			dKinematics[7]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid));
			dKinematics[5]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel));
			dKinematics[3]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec));
			dKinematics[1]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			//truth match to real, accepted and reconstructed D0
			if(DTRUEIDX->at(s)<0) continue;
			int d = DTRUEIDX->at(s);
			if(TRUEDTRK0IDX->at(d)==-1) continue;
			if(TRUEDTRK2IDX->at(d)!=-1) continue;
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;
			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			fndDrecoPT.at(4)->Fill(JetPT,weight);
			fndDrecoPT.at(3)->Fill(JetPT,weight/effpid);
			fndDrecoPT.at(2)->Fill(JetPT,weight/(effpid*effsel));
			fndDrecoPT.at(1)->Fill(JetPT,weight/(effpid*effsel*effrec*corrPt));
			fndDrecoPT.at(0)->Fill(JetPT,weight/(effpid*effsel*effrec*corrPt*effacc));

			fndDtruePT.at(4)->Fill(JetTruePT,weight);
			fndDtruePT.at(3)->Fill(JetTruePT,weight/effpid);
			fndDtruePT.at(2)->Fill(JetTruePT,weight/(effpid*effsel));
			fndDtruePT.at(1)->Fill(JetTruePT,weight/(effpid*effsel*effrec*corrPt));
			fndDtruePT.at(0)->Fill(JetTruePT,weight/(effpid*effsel*effrec*corrPt*effacc));

			fndD2DPT.at(4)->Fill(JetTruePT,JetPT,weight);
			fndD2DPT.at(3)->Fill(JetTruePT,JetPT,weight/effpid);
			fndD2DPT.at(2)->Fill(JetTruePT,JetPT,weight/(effpid*effsel));
			fndD2DPT.at(1)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec*corrPt));
			fndD2DPT.at(0)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec*corrPt*effacc));

			foundDs.push_back(d*100+s);

		//TODO	break;//only keep one D0 candidate per entry
		}
		for(unsigned int i=0; i<trueDs.size(); ++i) {
			bool matched(false);
			for(unsigned int j=0; j<foundDs.size(); ++j) {
				if(foundDs[j]==trueDs[i]) {
					matched=true;
					break;
				}
			}
			if(!matched) std::cout << progress.count()-1 << ": truD " << trueDs[i] << " has no matching fndD" << std::endl;
		}
		for(unsigned int i=0; i<foundDs.size(); ++i) {
			bool matched(false);
			for(unsigned int j=0; j<trueDs.size(); ++j) {
				if(foundDs[i]==trueDs[j]) {
					matched=true;
					break;
				}
			}
			if(!matched) std::cout << progress.count()-1 << ": fndD " << foundDs[i] << " has no matching truD" << std::endl;
		}
	}

	effData->Divide(numData,denData);
	for(uint i=1; i<_ptBins.size(); ++i) {
		effDataJPTBins[i-1]->Divide(numDataJPTBins[i-1],denDataJPTBins[i-1]);
		effDataJTRUEPTBins[i-1]->Divide(numDataJTRUEPTBins[i-1],denDataJTRUEPTBins[i-1]);
	}

	printf("bins of true jet pT\n");
	for(int i=0; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int j=1; j<=truDtruePT.at(i)->GetNbinsX(); ++j) {
			printf("%6.0f-%6.0f\t",truDtruePT.at(i)->GetBinLowEdge(j),truDtruePT.at(i)->GetBinLowEdge(j+1));
			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f",truDtruePT.at(i)->GetBinContent(j),effDtruePT.at(i)->GetBinContent(j),fndDtruePT.at(i)->GetBinContent(j),fitDtruePT.at(i)->GetBinContent(j));
			if(i>0) {
				printf("\t%5.3f", effDtruePT.at(i-1)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f\t%5.3f", truDtruePT.at(i)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j),
						         fndDtruePT.at(i)->GetBinContent(j)/fndDtruePT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f", (fndDtruePT.at(i)->GetBinContent(j)/fndDtruePT.at(i-1)->GetBinContent(j))/
						  (truDtruePT.at(i)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j)));
			} else {
				printf("\t%5.3f", (fndDtruePT.at(4)->GetBinContent(j)/fndDtruePT.at(0)->GetBinContent(j))/
						  (truDtruePT.at(4)->GetBinContent(j)/truDtruePT.at(0)->GetBinContent(j)));
				if(j>1 && j<truDtruePT.at(i)->GetNbinsX()) {
				closure_true->SetBinContent(j-1,(fndDtruePT.at(4)->GetBinContent(j)/fndDtruePT.at(0)->GetBinContent(j))/
						                (truDtruePT.at(4)->GetBinContent(j)/truDtruePT.at(0)->GetBinContent(j)));
				}
			}
			printf("\n");
		}
	}

	printf("bins of reco jet pT\n");
	for(int i=0; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int j=1; j<=truDrecoPT.at(i)->GetNbinsX(); ++j) {
			printf("%6.0f-%6.0f\t",truDrecoPT.at(i)->GetBinLowEdge(j),truDrecoPT.at(i)->GetBinLowEdge(j+1));
			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f\t",truDrecoPT.at(i)->GetBinContent(j),effDrecoPT.at(i)->GetBinContent(j),fndDrecoPT.at(i)->GetBinContent(j),fitDrecoPT.at(i)->GetBinContent(j));
			if(i>0) {
				printf("\t%5.3f", effDrecoPT.at(i-1)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f\t%5.3f", truDrecoPT.at(i)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j),
						         fndDrecoPT.at(i)->GetBinContent(j)/fndDrecoPT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f", (fndDrecoPT.at(i)->GetBinContent(j)/fndDrecoPT.at(i-1)->GetBinContent(j))/
						  (truDrecoPT.at(i)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j)));
			} else {
				printf("\t%5.3f", (fndDrecoPT.at(4)->GetBinContent(j)/fndDrecoPT.at(0)->GetBinContent(j))/
						  (truDrecoPT.at(4)->GetBinContent(j)/truDrecoPT.at(0)->GetBinContent(j)));
				if(j>1 && j<truDrecoPT.at(i)->GetNbinsX()) {
				closure_reco->SetBinContent(j-1,(fndDrecoPT.at(4)->GetBinContent(j)/fndDrecoPT.at(0)->GetBinContent(j))/
						                (truDrecoPT.at(4)->GetBinContent(j)/truDrecoPT.at(0)->GetBinContent(j)));
				}
			}
			printf("\n");
		}
	}

	printf("bins of (X) true and (Y) reco jet pT\neff(found)/eff(truth)\n");
	for(int i=1; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
			printf("\t%3.0f-%3.0f\t",truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
		}
		printf("\n");
		for(int j=1; j<=truD2DPT.at(i)->GetNbinsY(); ++j) {
			printf("%3.0f-%3.0f",truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
			for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
				double ratio = (fndD2DPT.at(i)->GetBinContent(k,j)/fndD2DPT.at(i-1)->GetBinContent(k,j))/
				               (truD2DPT.at(i)->GetBinContent(k,j)/truD2DPT.at(i-1)->GetBinContent(k,j));
				double error = ratio*TMath::Sqrt(TMath::Power(fndD2DPT.at(i  )->GetBinError(k,j)/fndD2DPT.at(i  )->GetBinContent(k,j),2.)
						             //+ TMath::Power(fndD02DPT.at(i-1)->GetBinError(k,j)/fndD02DPT.at(i-1)->GetBinContent(k,j),2.)
						               + TMath::Power(truD2DPT.at(i  )->GetBinError(k,j)/truD2DPT.at(i  )->GetBinContent(k,j),2.));
						             //+ TMath::Power(truD02DPT.at(i-1)->GetBinError(k,j)/truD02DPT.at(i-1)->GetBinContent(k,j),2.));
				printf("\t%5.3f+/-%5.3f", ratio, error);
			}
			printf("\n");
		}
	}

	printf("bins of (X) true and (Y) reco jet pT\nN(found)\n");
	for(int i=1; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
			printf("\t%3.0f-%3.0f\t",truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
		}
		printf("\n");
		for(int j=1; j<=truD2DPT.at(i)->GetNbinsY(); ++j) {
			printf("%3.0f-%3.0f",truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
			for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
				printf("\t%6.0f\t", fndD2DPT.at(i)->GetBinContent(k,j));
			}
			printf("\n");
		}
	}

	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 0.00, 0.00};
	Double_t greens[NRGBs] = { 0.00, 0.50, 1.00, 0.50, 0.00};
	Double_t blues[NRGBs]  = { 0.00, 0.00, 1.00, 1.00, 1.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);
	TCanvas c1;
	c1.SetLogx();
	dKinematics[9]->Divide(dKinematics[7]);
	dKinematics[8]->Divide(dKinematics[6]);
	dKinematics[9]->Divide(dKinematics[8]);
	dKinematics[9]->SetMinimum(0.);
	dKinematics[9]->SetMaximum(2.);
	dKinematics[9]->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"KinePIDRatio.pdf");
	dKinematics[7]->Divide(dKinematics[5]);
	dKinematics[6]->Divide(dKinematics[4]);
	dKinematics[7]->Divide(dKinematics[6]);
	dKinematics[7]->SetMinimum(0.);
	dKinematics[7]->SetMaximum(2.);
	dKinematics[7]->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineSELRatio.pdf");
	dKinematics[5]->Divide(dKinematics[3]);
	dKinematics[4]->Divide(dKinematics[2]);
	dKinematics[5]->Divide(dKinematics[4]);
	dKinematics[5]->SetMinimum(0.);
	dKinematics[5]->SetMaximum(2.);
	dKinematics[5]->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineRECRatio.pdf");
	dKinematics[3]->Divide(dKinematics[1]);
	dKinematics[2]->Divide(dKinematics[0]);
	dKinematics[3]->Divide(dKinematics[2]);
	dKinematics[3]->SetMinimum(0.);
	dKinematics[3]->SetMaximum(2.);
	dKinematics[3]->Draw("colz");
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineACCRatio.pdf");
	//for(unsigned int i=0; i<dKinematics.size(); ++i) {
	//	dKinematics[i]->Draw("colz");
	//	TString plotname("D0Kinematics");
	//	plotname+=i;
	//	plotname+=".pdf";
	//	c1.SaveAs(gSaveDir+"/"+plotname);
	//}
	
	for(int i=1; i<=dKinematics[3]->GetNbinsX(); ++i) {
		for(int j=1; j<=dKinematics[3]->GetNbinsY(); ++j) {
			if(dKinematics[3]->GetBinContent(i,j)>0) {
				dKinematicsPulls[0]->Fill((dKinematics[3]->GetBinContent(i,j)-1.)/dKinematics[3]->GetBinError(i,j));
			}
			if(dKinematics[5]->GetBinContent(i,j)>0) {
				dKinematicsPulls[1]->Fill((dKinematics[5]->GetBinContent(i,j)-1.)/dKinematics[5]->GetBinError(i,j));
			}
			if(dKinematics[7]->GetBinContent(i,j)>0) {
				dKinematicsPulls[2]->Fill((dKinematics[7]->GetBinContent(i,j)-1.)/dKinematics[7]->GetBinError(i,j));
			}
			if(dKinematics[9]->GetBinContent(i,j)>0) {
				dKinematicsPulls[3]->Fill((dKinematics[9]->GetBinContent(i,j)-1.)/dKinematics[9]->GetBinError(i,j));
			}
		}
	}
	c1.SetLogx(0);
	dKinematicsPulls[0]->Draw();
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineACCPulls.pdf");
	dKinematicsPulls[1]->Draw();
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineRECPulls.pdf");
	dKinematicsPulls[2]->Draw();
	c1.SaveAs(gSaveDir+"/"+_whichD+"KineSELPulls.pdf");
	dKinematicsPulls[3]->Draw();
	c1.SaveAs(gSaveDir+"/"+_whichD+"KinePIDPulls.pdf");

	const Int_t NRGBs2 = 5;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds2[NRGBs2]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
	Double_t greens2[NRGBs2] = { 1.00, 0.95, 0.50, 0.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);

	c1.Clear();
	c1.Divide(4,5);
	for(int i=0; i<5; ++i) {
		double max = 1.1*TMath::Max(TMath::Max(truD2DPT.at(i)->GetMaximum(),effD2DPT.at(i)->GetMaximum()),TMath::Max(fndD2DPT.at(i)->GetMaximum(),fitD2DPT.at(i)->GetMaximum()));
		truD2DPT.at(i)->SetMinimum(0.);
		truD2DPT.at(i)->SetMaximum(max);
		fndD2DPT.at(i)->SetMinimum(0.);
		fndD2DPT.at(i)->SetMaximum(max);
		fitD2DPT.at(i)->SetMinimum(0.);
		fitD2DPT.at(i)->SetMaximum(max);
		effD2DPT.at(i)->SetMinimum(0.);
		effD2DPT.at(i)->SetMaximum(max);
		c1.cd(4*i+1);
		truD2DPT.at(i)->Draw("colz");
		c1.cd(4*i+2);
		fndD2DPT.at(i)->Draw("colz");
		c1.cd(4*i+3);
		fitD2DPT.at(i)->Draw("colz");
		c1.cd(4*i+4);
		effD2DPT.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/"+_whichD+"2DEffMaps.pdf");

	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);
	c1.Clear();
	c1.Divide(5,3);
	for(int i=0; i<5; ++i) {
		fndD2DPT.at(i)->Divide(truD2DPT.at(i));
		fitD2DPT.at(i)->Divide(truD2DPT.at(i));
		effD2DPT.at(i)->Divide(truD2DPT.at(i));
		double max = 2.0;
		fndD2DPT.at(i)->SetMinimum(0.);
		fndD2DPT.at(i)->SetMaximum(max);
		fitD2DPT.at(i)->SetMinimum(0.);
		fitD2DPT.at(i)->SetMaximum(max);
		effD2DPT.at(i)->SetMinimum(0.);
		effD2DPT.at(i)->SetMaximum(max);
		c1.cd(i+1);
		fndD2DPT.at(i)->Draw("colz");
		c1.cd(i+6);
		fitD2DPT.at(i)->Draw("colz");
		c1.cd(i+11);
		effD2DPT.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/2DEffMapsRatios.pdf");

	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);
	c1.Clear();
	c1.Divide(6,7);
	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
		c1.cd(i+1);
		dKinematicsTRBinned.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin.pdf");
	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
		c1.cd(i+1);
		dKinematicsTRBinned_recD_tru.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin_recD_truKine.pdf");
	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
		c1.cd(i+1);
		dKinematicsTRBinned_recD_rec.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin_recD_recKine.pdf");

	TString closureOutName = gSaveDir+"/effClosure";
	closureOutName+=flavour;
	if(_useSimpleEffs) closureOutName+="_simpleEffs";
	closureOutName+=".root";
	TFile* fout = TFile::Open(closureOutName,"RECREATE");
	closure_true->Write();
	closure_reco->Write();
	effModel->Write();
	effData->Write();
	numData->Write();
	denData->Write();
	for(uint i=0; i<_ptBins.size()-1; ++i) {
		effDataJPTBins[i]->Write();
		numDataJPTBins[i]->Write();
		denDataJPTBins[i]->Write();
		effDataJTRUEPTBins[i]->Write();
		numDataJTRUEPTBins[i]->Write();
		denDataJTRUEPTBins[i]->Write();
	}
	fout->Close();

	dm->reset();

	if(!_useSimpleEffs) {
		//rerun with simple effs
		_useSimpleEffs=true;
		testEffs(flavour);
		_useSimpleEffs=false;
	}

	return true;
}

bool SimDFitter::testEffsPtBinned(int flavour) {
	if(!_dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << _name << std::endl;
	std::cout << "TODO : not written yet" << std::endl;
//	DatasetManager* dm = DatasetManager::getInstance();
//	dm->setDataset(_dataFile);
//
//	TFile* f2 = TFile::Open("../efficiencies/Effs_ptbinned_16x8bins_up210319.root");
//
//	//construct efficiency objects
//	std::vector<TH2D*> hpid4, hpid5, hnum4, hnum5, hden4, hden5, heffS4, heffS5, heffP4, heffP5, heff4, heff5;
//	for(uint i=1; i<=4; ++i) {
//		TString formatStr = "%s"+_whichD+"%d_%d";
//		hpid4.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"pid",4,i))));
//		hpid5.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"pid",5,i))));
//		hnum4.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"num",4,i))));
//		hnum5.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"num",5,i))));
//		hden4.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"denom",4,i))));
//		hden5.push_back(static_cast<TH2D*>(f2->Get(TString::Format(formatStr,"denom",5,i))));
//		if(!hpid4.back() || !hnum4.back() || !hden4.back()) {
//			std::cout << "!" << _whichD << std::endl;
//			f2->ls();
//			std::cout << hpid4.back() <<" "<< << hnum4.back() <<" "<< hden4.back() << std::endl;
//			std::cout << "4 histograms missing" << std::endl;
//			return false;
//		}
//		if(!hpid5.back() || !hnum5.back() || !hden5.back()) {
//			std::cout << "!" << _whichD << std::endl;
//			f2->ls();
//			std::cout << hpid5.back() <<" "<< << hnum5.back() <<" "<< hden5.back() << std::endl;
//			std::cout << "5 histograms missing" << std::endl;
//			return false;
//		}
//
//		heffS4.push_back(static_cast<TH2D*>(hnum4.back()->Clone(TString::Format(formatStr,"effS",4,i))));
//		heffS5.push_back(static_cast<TH2D*>(hnum5.back()->Clone(TString::Format(formatStr,"effS",5,i))));
//		heffS4.back()->Divide(hden4.back());
//		heffS5.back()->Divide(hden5.back());
//
//		heffP4.push_back(static_cast<TH2D*>(hpid4.back()->Clone(TString::Format(formatStr,"effP",4,i))));
//		heffP5.push_back(static_cast<TH2D*>(hpid5.back()->Clone(TString::Format(formatStr,"effP",5,i))));
//		heffP4.back()->Divide(hden4.back());
//		heffP5.back()->Divide(hden5.back());
//
//		heff4.push_back(static_cast<TH2D*>(hpid4.back()->Clone(TString::Format(formatStr,"eff",4,i))));
//		heff5.push_back(static_cast<TH2D*>(hpid5.back()->Clone(TString::Format(formatStr,"eff",5,i))));
//		heff4.back()->Divide(hden4.back());
//		heff5.back()->Divide(hden5.back());
//	}
//
//	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pideff"+_whichD+"45"));
//	TH2D *hacc(0), *hrec(0), *hsel(0);
//	if(_flavour==5) hsel = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"5"));
//	else hsel = dynamic_cast<TH2D*>(f2->Get("seleff"+_whichD+"4"));
//
//	double unitsScale(1.);
//	if(!_useSimpleEffs) {
//		hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
//		hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
//		unitsScale=1e-3;
//	} else {
//		hacc  = dynamic_cast<TH2D*>(f2->Get("acceff"+_whichD+"45"));
//		hrec  = dynamic_cast<TH2D*>(f2->Get("receff"+_whichD+"45"));
//	}
//
//	if(!hacc || !hrec /*|| !hcor*/ || !hpid || !hsel) return false;
//
//	TH2D* effModel = cloneTH2D("effModel", hpid);
//	TH2D* denData  = cloneTH2D("denData", hpid);
//	TH2D* numData  = cloneTH2D("numData", hpid);
//	TH2D* effData  = cloneTH2D("effData", hpid);
//	std::vector<TH2D*> denDataJPTBins;
//	std::vector<TH2D*> numDataJPTBins;
//	std::vector<TH2D*> effDataJPTBins;
//	std::vector<TH2D*> denDataJTRUEPTBins;
//	std::vector<TH2D*> numDataJTRUEPTBins;
//	std::vector<TH2D*> effDataJTRUEPTBins;
//	for(uint i=1; i<_ptBins.size(); ++i) {
//		denDataJPTBins.push_back(cloneTH2D(TString::Format("denDataJPTBin%d",i), hpid));
//		numDataJPTBins.push_back(cloneTH2D(TString::Format("numDataJPTBin%d",i), hpid));
//		effDataJPTBins.push_back(cloneTH2D(TString::Format("effDataJPTBin%d",i), hpid));
//		denDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("denDataJTRUEPTBin%d",i), hpid));
//		numDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("numDataJTRUEPTBin%d",i), hpid));
//		effDataJTRUEPTBins.push_back(cloneTH2D(TString::Format("effDataJTRUEPTBin%d",i), hpid));
//	}
//
//	//add over/underflow to the pT binning
//	std::vector<double> ptBinsWithOF = {0.};
//	for(uint i=0; i<_ptBins.size(); ++i) ptBinsWithOF.push_back(_ptBins.at(i));
//	ptBinsWithOF.push_back(110e3);
//	for(uint i=0; i<ptBinsWithOF.size(); ++i) std::cout << ptBinsWithOF.at(i) << std::endl;//TODO
//
//	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
//	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
//	std::vector<TString> stages;
//	stages.push_back("all");
//	stages.push_back("geometric");
//	stages.push_back("reconstruction");
//	stages.push_back("selection");
//	stages.push_back("PID");
//	std::vector<TH1D*> truDtruePT;//the true distribution at this stage
//	truDtruePT.push_back(cloneTH1D("tru0truept", ptBinsWithOF));
//	truDtruePT.push_back(cloneTH1D("tru4truept", ptBinsWithOF));
//	std::vector<TH1D*> effDtruePT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
//	effDtruePT.push_back(cloneTH1D("eff0truept", ptBinsWithOF));
//	effDtruePT.push_back(cloneTH1D("eff4truept", ptBinsWithOF));
//	std::vector<TH1D*> fitDtruePT;//sideband-subtracted distribution at this stage
//	fitDtruePT.push_back(cloneTH1D("fit0truept", ptBinsWithOF));
//	fitDtruePT.push_back(cloneTH1D("fit4truept", ptBinsWithOF));
//	std::vector<TH1D*> fndDtruePT;//truth-matched distribution at this stage
//	fndDtruePT.push_back(cloneTH1D("fnd0truept", ptBinsWithOF));
//	fndDtruePT.push_back(cloneTH1D("fnd4truept", ptBinsWithOF));
//	//binned in true jet pT
//	std::vector<TH1D*> truDrecoPT;
//	truDrecoPT.push_back(cloneTH1D("tru0recopt", ptBinsWithOF));
//	truDrecoPT.push_back(cloneTH1D("tru4recopt", ptBinsWithOF));
//	std::vector<TH1D*> effDrecoPT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
//	effDrecoPT.push_back(cloneTH1D("eff0recopt", ptBinsWithOF));
//	effDrecoPT.push_back(cloneTH1D("eff4recopt", ptBinsWithOF));
//	std::vector<TH1D*> fitDrecoPT;
//	fitDrecoPT.push_back(cloneTH1D("fit0recopt", ptBinsWithOF));
//	fitDrecoPT.push_back(cloneTH1D("fit4recopt", ptBinsWithOF));
//	std::vector<TH1D*> fndDrecoPT;
//	fndDrecoPT.push_back(cloneTH1D("fnd0recopt", ptBinsWithOF));
//	fndDrecoPT.push_back(cloneTH1D("fnd4recopt", ptBinsWithOF));
//	//binned in both reco and true jet pT
//	std::vector<TH2D*> truD2DPT;
//	truD2DPT.push_back(cloneTH2D("tru02Dpt", ptBinsWithOF, ptBinsWithOF));
//	truD2DPT.push_back(cloneTH2D("tru42Dpt", ptBinsWithOF, ptBinsWithOF));
//	std::vector<TH2D*> effD2DPT;
//	effD2DPT.push_back(cloneTH2D("eff02Dpt", ptBinsWithOF, ptBinsWithOF));
//	effD2DPT.push_back(cloneTH2D("eff42Dpt", ptBinsWithOF, ptBinsWithOF));
//	std::vector<TH2D*> fitD2DPT;
//	fitD2DPT.push_back(cloneTH2D("fit02Dpt", ptBinsWithOF, ptBinsWithOF));
//	fitD2DPT.push_back(cloneTH2D("fit42Dpt", ptBinsWithOF, ptBinsWithOF));
//	std::vector<TH2D*> fndD2DPT;
//	fndD2DPT.push_back(cloneTH2D("fnd02Dpt", ptBinsWithOF, ptBinsWithOF));
//	fndD2DPT.push_back(cloneTH2D("fnd42Dpt", ptBinsWithOF, ptBinsWithOF));
//
//	//entries in this vector are true and reco for each stage
//	std::vector<TH2D*> dKinematics;
//	dKinematics.push_back(cloneTH2D("dKineTrue0", hsel));
//	dKinematics.push_back(cloneTH2D("dKineReco0", hsel));
//	dKinematics.push_back(cloneTH2D("dKineTrue4", hsel));
//	dKinematics.push_back(cloneTH2D("dKineReco4", hsel));
//	for(unsigned int i=0; i<dKinematics.size(); ++i) {
//		dKinematics[i]->Sumw2();
//	}
//	std::vector<TH1D*> dKinematicsPulls;
//	dKinematicsPulls.push_back(new TH1D("pulls1", "", 30, -1., 1.));
//	dKinematicsPulls.push_back(new TH1D("pulls4", "", 30, -1., 1.));
//
//	TH1D* closure_true = cloneTH1D("closure_true", _ptBins);
//	TH1D* closure_reco = cloneTH1D("closure_reco", _ptBins);
//
//	//also divide distribution into bins of true and reco pT
//	std::vector<TH2D*> dKinematicsTRBinned;
//	std::vector<TH2D*> dKinematicsTRBinned_recD_tru;
//	std::vector<TH2D*> dKinematicsTRBinned_recD_rec;
//	for(uint i=0; i<ptBinsWithOF.size(); ++i) {
//		for(uint j=0; j<ptBinsWithOF.size()-1; ++j) {
//			TString hName = TString::Format("dKineJPTBin%d%d",i,j);
//			TString hTitle = TString::Format("true%d reco%d",j,i);
//			if(j==ptBinsWithOF.size()-2) {
//				hTitle.Form("reco%d",i);
//			}
//			if(i==ptBinsWithOF.size()-1) {
//				hTitle.Form("true%d",j);
//			}
//			if(j==ptBinsWithOF.size()-2 && i==ptBinsWithOF.size()-1) {
//				hTitle = "all";
//			}
//			dKinematicsTRBinned.push_back(cloneTH2D(hName, hsel));
//			dKinematicsTRBinned.back()->Sumw2();
//			dKinematicsTRBinned.back()->SetTitle(hTitle);
//			dKinematicsTRBinned.back()->GetXaxis()->SetTitle("pT(D)");
//			dKinematicsTRBinned.back()->GetYaxis()->SetTitle("#eta(D)");
//			
//			dKinematicsTRBinned_recD_tru.push_back(cloneTH2D(hName+"_recD_tru", hsel));
//			dKinematicsTRBinned_recD_tru.back()->Sumw2();
//			dKinematicsTRBinned_recD_tru.back()->SetTitle(hTitle);
//			dKinematicsTRBinned_recD_tru.back()->GetXaxis()->SetTitle("pT(D)");
//			dKinematicsTRBinned_recD_tru.back()->GetYaxis()->SetTitle("#eta(D)");
//			
//			dKinematicsTRBinned_recD_rec.push_back(cloneTH2D(hName+"_recD_rec", hsel));
//			dKinematicsTRBinned_recD_rec.back()->Sumw2();
//			dKinematicsTRBinned_recD_rec.back()->SetTitle(hTitle);
//			dKinematicsTRBinned_recD_rec.back()->GetXaxis()->SetTitle("pT(D)");
//			dKinematicsTRBinned_recD_rec.back()->GetYaxis()->SetTitle("#eta(D)");
//		}
//	}
//
//	for(int i=0; i<2; ++i) {
//		truDtruePT.at(i)->Sumw2();
//		effDtruePT.at(i)->Sumw2();
//		fitDtruePT.at(i)->Sumw2();
//		fndDtruePT.at(i)->Sumw2();
//		truDrecoPT.at(i)->Sumw2();
//		effDrecoPT.at(i)->Sumw2();
//		fitDrecoPT.at(i)->Sumw2();
//		fndDrecoPT.at(i)->Sumw2();
//		truD2DPT.at(i)->Sumw2();
//		effD2DPT.at(i)->Sumw2();
//		fitD2DPT.at(i)->Sumw2();
//		fndD2DPT.at(i)->Sumw2();
//		truD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
//		effD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
//		fitD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
//		fndD2DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
//		truD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
//		effD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
//		fitD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
//		fndD2DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
//	}
//
//	std::vector<double> *DM = new std::vector<double>();
//	std::vector<double> *DIPCHI2 = new std::vector<double>();
//	std::vector<double> *DPT = new std::vector<double>();
//	std::vector<double> *DPX = new std::vector<double>();
//	std::vector<double> *DPY = new std::vector<double>();
//	std::vector<double> *DPZ = new std::vector<double>();
//	std::vector<double> *DE = new std::vector<double>();
//	std::vector<double> *DX = new std::vector<double>();
//	std::vector<double> *DY = new std::vector<double>();
//	std::vector<double> *DZ = new std::vector<double>();
//	std::vector<double> *Dh0P = new std::vector<double>();
//	std::vector<double> *Dh0PT = new std::vector<double>();
//	std::vector<double> *Dh0PX = new std::vector<double>();
//	std::vector<double> *Dh0PY = new std::vector<double>();
//	std::vector<double> *Dh0PZ = new std::vector<double>();
//	std::vector<double> *Dh1P = new std::vector<double>();
//	std::vector<double> *Dh1PT = new std::vector<double>();
//	std::vector<double> *Dh1PX = new std::vector<double>();
//	std::vector<double> *Dh1PY = new std::vector<double>();
//	std::vector<double> *Dh1PZ = new std::vector<double>();
//	std::vector<double> *Dh2P = new std::vector<double>();
//	std::vector<double> *Dh2PT = new std::vector<double>();
//	std::vector<double> *Dh2PX = new std::vector<double>();
//	std::vector<double> *Dh2PY = new std::vector<double>();
//	std::vector<double> *Dh2PZ = new std::vector<double>();
//	std::vector<double> *Dh0WEIGHT = new std::vector<double>();
//	std::vector<double> *Dh1WEIGHT = new std::vector<double>();
//	std::vector<double> *Dh2WEIGHT = new std::vector<double>();
//	std::vector<double> *DTRUEIDX = new std::vector<double>();
//	std::vector<double> *DTRUETRK0 = new std::vector<double>();
//	std::vector<double> *DTRUETRK1 = new std::vector<double>();
//	std::vector<double> *DTRUETRK2 = new std::vector<double>();
//
//	std::vector<double> *TRUEDID = new std::vector<double>();
//	std::vector<double> *TRUEDPX = new std::vector<double>();
//	std::vector<double> *TRUEDPY = new std::vector<double>();
//	std::vector<double> *TRUEDPZ = new std::vector<double>();
//	std::vector<double> *TRUEDE = new std::vector<double>();
//	std::vector<double> *TRUEDX = new std::vector<double>();
//	std::vector<double> *TRUEDY = new std::vector<double>();
//	std::vector<double> *TRUEDZ = new std::vector<double>();
//	std::vector<double> *TRUEDFROMB = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0P = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0PT = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0INACC = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0RECO = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1P = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1PT = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1INACC = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1RECO = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0IDX = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1IDX = new std::vector<double>();
//	std::vector<double> *TRUEDTRK2IDX = new std::vector<double>();
//	std::vector<double> *TRUEDTRK3IDX = new std::vector<double>();
//	std::vector<double> *TRUEDTRK0ID = new std::vector<double>();
//	std::vector<double> *TRUEDTRK1ID = new std::vector<double>();
//	std::vector<double> *TRUEDTRK2ID = new std::vector<double>();
//	std::vector<double> *TRUEDTRK3ID = new std::vector<double>();
//
//	double JetPT;
//	double JetEta;
//	double JetTruePT;
//
//	dm->setBranchAddress(""+_whichD+"M",           &DM);
//	dm->setBranchAddress(""+_whichD+"IPCHI2",      &DIPCHI2);
//	dm->setBranchAddress(""+_whichD+"PT",          &DPT);
//	dm->setBranchAddress(""+_whichD+"PX",          &DPX);
//	dm->setBranchAddress(""+_whichD+"PY",          &DPY);
//	dm->setBranchAddress(""+_whichD+"PZ",          &DPZ);
//	dm->setBranchAddress(""+_whichD+"E",           &DE);
//	dm->setBranchAddress(""+_whichD+"X",           &DX);
//	dm->setBranchAddress(""+_whichD+"Y",           &DY);
//	dm->setBranchAddress(""+_whichD+"Z",           &DZ);
//	if(_whichD=="D0") {
//		dm->setBranchAddress("D0KP",          &Dh0P);
//		dm->setBranchAddress("D0KPT",         &Dh0PT);
//		dm->setBranchAddress("D0KPX",         &Dh0PX);
//		dm->setBranchAddress("D0KPY",         &Dh0PY);
//		dm->setBranchAddress("D0KPZ",         &Dh0PZ);
//		dm->setBranchAddress("D0PIP",         &Dh1P);
//		dm->setBranchAddress("D0PIPT",        &Dh1PT);
//		dm->setBranchAddress("D0PIPX",        &Dh1PX);
//		dm->setBranchAddress("D0PIPY",        &Dh1PY);
//		dm->setBranchAddress("D0PIPZ",        &Dh1PZ);
//		dm->setBranchAddress("D0KWEIGHT",     &Dh0WEIGHT);
//		dm->setBranchAddress("D0PIWEIGHT",    &Dh1WEIGHT);
//		dm->setBranchAddress("D0TRUEIDX",     &DTRUEIDX);
//		dm->setBranchAddress("D0TRUETRK0",    &DTRUETRK0);
//		dm->setBranchAddress("D0TRUETRK1",    &DTRUETRK1);
//	} if(_whichD=="D") {
//		dm->setBranchAddress("DKP",          &Dh0P);
//		dm->setBranchAddress("DKPT",         &Dh0PT);
//		dm->setBranchAddress("DKPX",         &Dh0PX);
//		dm->setBranchAddress("DKPY",         &Dh0PY);
//		dm->setBranchAddress("DKPZ",         &Dh0PZ);
//		dm->setBranchAddress("DPI1P",        &Dh1P);
//		dm->setBranchAddress("DPI1PT",       &Dh1PT);
//		dm->setBranchAddress("DPI1PX",       &Dh1PX);
//		dm->setBranchAddress("DPI1PY",       &Dh1PY);
//		dm->setBranchAddress("DPI1PZ",       &Dh1PZ);
//		dm->setBranchAddress("DPI2P",        &Dh2P);
//		dm->setBranchAddress("DPI2PT",       &Dh2PT);
//		dm->setBranchAddress("DPI2PX",       &Dh2PX);
//		dm->setBranchAddress("DPI2PY",       &Dh2PY);
//		dm->setBranchAddress("DPI2PZ",       &Dh2PZ);
//		dm->setBranchAddress("DKWEIGHT",     &Dh0WEIGHT);
//		dm->setBranchAddress("DPI1WEIGHT",   &Dh1WEIGHT);
//		dm->setBranchAddress("DPI2WEIGHT",   &Dh2WEIGHT);
//		dm->setBranchAddress("DTRUEIDX",     &DTRUEIDX);
//		dm->setBranchAddress("DTRUETRK0",    &DTRUETRK0);
//		dm->setBranchAddress("DTRUETRK1",    &DTRUETRK1);
//		dm->setBranchAddress("DTRUETRK2",    &DTRUETRK2);
//	}
//
//	dm->setBranchAddress("TRUEDID",       &TRUEDID);
//	dm->setBranchAddress("TRUEDPX",       &TRUEDPX);
//	dm->setBranchAddress("TRUEDPY",       &TRUEDPY);
//	dm->setBranchAddress("TRUEDPZ",       &TRUEDPZ);
//	dm->setBranchAddress("TRUEDE",        &TRUEDE);
//	dm->setBranchAddress("TRUEDX",        &TRUEDX);
//	dm->setBranchAddress("TRUEDY",        &TRUEDY);
//	dm->setBranchAddress("TRUEDZ",        &TRUEDZ);
//	dm->setBranchAddress("TRUEDFROMB",    &TRUEDFROMB);
//	dm->setBranchAddress("TRUEDTRK0P",    &TRUEDTRK0P);
//	dm->setBranchAddress("TRUEDTRK0PT",   &TRUEDTRK0PT);
//	dm->setBranchAddress("TRUEDTRK0INACC",&TRUEDTRK0INACC);
//	dm->setBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO);
//	dm->setBranchAddress("TRUEDTRK1P",    &TRUEDTRK1P);
//	dm->setBranchAddress("TRUEDTRK1PT",   &TRUEDTRK1PT);
//	dm->setBranchAddress("TRUEDTRK1INACC",&TRUEDTRK1INACC);
//	dm->setBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO);
//	dm->setBranchAddress("TRUEDTRK0IDX",  &TRUEDTRK0IDX);
//	dm->setBranchAddress("TRUEDTRK1IDX",  &TRUEDTRK1IDX);
//	dm->setBranchAddress("TRUEDTRK2IDX",  &TRUEDTRK2IDX);
//	dm->setBranchAddress("TRUEDTRK3IDX",  &TRUEDTRK3IDX);
//	dm->setBranchAddress("TRUEDTRK0ID",  &TRUEDTRK0ID);
//	dm->setBranchAddress("TRUEDTRK1ID",  &TRUEDTRK1ID);
//	dm->setBranchAddress("TRUEDTRK2ID",  &TRUEDTRK2ID);
//	dm->setBranchAddress("TRUEDTRK3ID",  &TRUEDTRK3ID);
//
//	dm->setBranchAddress("JetPT",         &JetPT);
//	dm->setBranchAddress("JetEta",        &JetEta);
//	dm->setBranchAddress("JetTruePT",     &JetTruePT);
//
//	unsigned int nentries0 = dm->getEntries();
//
//	double dpt(0.), deta(0.);
//
//	for(int i=1; i<=effModel->GetNbinsX(); ++i) {
//		for(int j=1; j<=effModel->GetNbinsY(); ++j) {
//			dpt = effModel->GetXaxis()->GetBinCenter(i);
//			deta = effModel->GetYaxis()->GetBinCenter(j);
//			double effacc(0.), effrec(0.), effsel(0.), effpid(0.), eff(0.);
//			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
//			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
//			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
//			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));
//			eff = effacc*effrec*effsel*effpid;
//			//std::cout << i << " " << j << "\t" << dpt << "\t" << deta << "\t" << eff << std::endl;
//			effModel->SetBinContent(i,j,eff);
//		}
//	}
//
//	boost::progress_display progress( nentries0 );
//	//for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
//	while(dm->getNext()) {
//		++progress;
//		//t0->GetEntry(ientry);
//
//		std::vector<int> trueDs, foundDs;
//
//		int jrptbin(0), jtptbin(0);
//		for(uint i=0; i<_ptBins.size(); ++i) {
//			if(JetPT>_ptBins[i]) ++jrptbin;
//			if(JetTruePT>_ptBins[i]) ++jtptbin;
//		}
//
//		//put anything out of range into the overflow bin
//		if(JetPT>_ptBins.back()) JetPT = _ptBins.back();
//		if(JetTruePT>_ptBins.back()) JetTruePT = _ptBins.back();
//
//		for(unsigned int d=0; d<TRUEDID->size(); ++d) {
//			if(_whichD=="D0") {
//				//only use D0->Kpi with Pt>5GeV
//				if(TMath::Abs(TRUEDID->at(d))!=421) continue;
//				if(TRUEDTRK0IDX->at(d)==-1) continue;
//				if(TRUEDTRK2IDX->at(d)!=-1) continue;
//			} else if(_whichD=="D") {
//				//only use D->Kpipi with Pt>5GeV
//				if(TMath::Abs(TRUEDID->at(d))!=411) continue;
//				if(TRUEDTRK0IDX->at(d)==-1) continue;
//				if(TRUEDTRK2IDX->at(d)==-1) continue;
//				if(TRUEDTRK3IDX->at(d)!=-1) continue;
//			}
//			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;
//
//			TVector3 TRUEDP (TRUEDPX->at(d)  ,TRUEDPY->at(d)  ,TRUEDPZ->at(d));
//
//			dpt        = TRUEDP.Pt();
//			deta       = TRUEDP.Eta();
//			//double rhoSq = TRUEDX->at(d)*TRUEDX->at(d) + TRUEDY->at(d)*TRUEDY->at(d);
//			//double z = TRUEDZ->at(d);
//			//if(rhoSq>=100) rhoSq=99.9;//TODO
//
//			//get efficiency corrections
//			double effacc(0.), effrec(0.), effsel(0.), effpid(0.);
//
//			if(dpt>=100000.) dpt=99999.;
//			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
//			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
//			//apply correction to the reco efficiency
//			//if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
//			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
//			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));
//
//			denData->Fill(dpt,deta);
//			if(jrptbin>0 && jrptbin<=denDataJPTBins.size()) denDataJPTBins[jrptbin-1]->Fill(dpt,deta);
//			if(jtptbin>0 && jtptbin<=denDataJPTBins.size()) denDataJTRUEPTBins[jtptbin-1]->Fill(dpt,deta);
//
//			truDrecoPT.at(0)->Fill(JetPT);
//			truDtruePT.at(0)->Fill(JetTruePT);
//			truD2DPT.at(0)->Fill(JetTruePT,JetPT);
//			dKinematics[0]->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//
//			//fill kinematics separately for each true/reco jet pT bin
//			int ibin = truDrecoPT.at(0)->FindBin(JetTruePT) - 1;
//			int jbin = truDrecoPT.at(0)->FindBin(JetPT) - 1;
//			dKinematicsTRBinned.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//			dKinematicsTRBinned.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
//			dKinematicsTRBinned.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
//			dKinematicsTRBinned.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size()-1)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all
//
//			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
//
//			truDrecoPT.at(1)->Fill(JetPT);
//			truDtruePT.at(1)->Fill(JetTruePT);
//			truD2DPT.at(1)->Fill(JetTruePT,JetPT);
//			if(effacc>0.) {
//				effDrecoPT.at(0)->Fill(JetPT,1./effacc);
//				effDtruePT.at(0)->Fill(JetTruePT,1./effacc);
//				effD2DPT.at(0)->Fill(JetTruePT,JetPT,1./effacc);
//			} else {
//				std::cout << "ACC=" << effacc << " at " << dpt << "," << deta << std::endl;
//			}
//			dKinematics[2]->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//
//			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;
//
//			truDrecoPT.at(2)->Fill(JetPT);
//			truDtruePT.at(2)->Fill(JetTruePT);
//			truD2DPT.at(2)->Fill(JetTruePT,JetPT);
//			if(effrec>0.) {
//				effDrecoPT.at(1)->Fill(JetPT,1./effrec);
//				effDtruePT.at(1)->Fill(JetTruePT,1./effrec);
//				effD2DPT.at(1)->Fill(JetTruePT,JetPT,1./effrec);
//			} else {
//				std::cout << "REC=" << effrec << " at " << dpt << "," << deta << std::endl;
//			}
//			dKinematics[4]->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//
//			for(unsigned int s=0; s<DTRUEIDX->size(); ++s) {
//				if(DTRUEIDX->at(s)==d) {
//					TVector3 DP (DPX->at(s)  ,DPY->at(s)  ,DPZ->at(s));
//					TVector3 DP0(Dh0PX->at(s) ,Dh0PY->at(s) ,Dh0PZ->at(s));
//					TVector3 DP1(Dh1PX->at(s),Dh1PY->at(s),Dh1PZ->at(s));
//
//					if(!(DP0.Eta()>2.0 && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
//					if(!(DP1.Eta()>2.0 && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
//					if(!(DP.Pt()>_dptmin)) continue;
//					if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
//					//if(!(D0P.Eta()>2.5&&D0P.Eta()<4.0)) continue;//TODO test restricting the eta(D0) range (turn on in two places 1/2)
//
//					truDrecoPT.at(3)->Fill(JetPT);
//					truDtruePT.at(3)->Fill(JetTruePT);
//					truD2DPT.at(3)->Fill(JetTruePT,JetPT);
//					if(effsel>0.) {
//						effDrecoPT.at(2)->Fill(JetPT,1./effsel);
//						effDtruePT.at(2)->Fill(JetTruePT,1./effsel);
//						effD2DPT.at(2)->Fill(JetTruePT,JetPT,1./effsel);
//					} else {
//						std::cout << "SEL=" << effsel << " at " << dpt << "," << deta << std::endl;
//					}
//					dKinematics[6]->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//
//					//use weights for PID
//					double weight = Dh0WEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
//					if(DP0.Pt()>25000. || DP0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
//					truDrecoPT.at(4)->Fill(JetPT,weight);
//					truDtruePT.at(4)->Fill(JetTruePT,weight);
//					truD2DPT.at(4)->Fill(JetTruePT,JetPT);
//					numData->Fill(dpt,deta,weight);
//					if(jrptbin>0 && jrptbin<=numDataJPTBins.size()) numDataJPTBins[jrptbin-1]->Fill(dpt,deta,weight);
//					if(jtptbin>0 && jtptbin<=numDataJPTBins.size()) numDataJTRUEPTBins[jtptbin-1]->Fill(dpt,deta,weight);
//					if(effpid>0.) {
//						effDrecoPT.at(3)->Fill(JetPT,weight/effpid);
//						effDtruePT.at(3)->Fill(JetTruePT,weight/effpid);
//						effD2DPT.at(3)->Fill(JetTruePT,JetPT,1./effpid);
//					} else {
//						std::cout << "PID=" << effpid << " at " << dpt << "," << deta << std::endl;
//					}
//					dKinematics[8]->Fill(TRUEDP.Pt(),TRUEDP.Eta(),weight);
//
//					trueDs.push_back(d*100+s);
//
//					//fill kinematics separately for each true/reco jet pT bin
//					dKinematicsTRBinned_recD_tru.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
//					dKinematicsTRBinned_recD_tru.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
//					dKinematicsTRBinned_recD_tru.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
//					dKinematicsTRBinned_recD_tru.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size() -1)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all
//
//					//fill kinematics separately for each true/reco jet pT bin
//					dKinematicsTRBinned_recD_rec.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());
//					dKinematicsTRBinned_recD_rec.at((ptBinsWithOF.size()-2) + jbin*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());//row sums
//					dKinematicsTRBinned_recD_rec.at(ibin + (ptBinsWithOF.size()-1)*(ptBinsWithOF.size()-1))->Fill(DP.Pt(),DP.Eta());//col sums
//					dKinematicsTRBinned_recD_rec.at((ptBinsWithOF.size()-1)*ptBinsWithOF.size() -1)->Fill(DP.Pt(),DP.Eta());//all
//
//					break;
//				}
//			}
//		}
//		for(unsigned int s=0; s<DM->size(); ++s) {
//			//first check in our tight acceptance
//			TVector3 DP (DPX->at(s)  ,DPY->at(s)  ,DPZ->at(s));
//			TVector3 DP0(Dh0PX->at(s),Dh0PY->at(s),Dh0PZ->at(s));
//			TVector3 DP1(Dh1PX->at(s),Dh1PY->at(s),Dh1PZ->at(s));
//			//double rhoSq = D0X->at(s)*D0X->at(s) + D0Y->at(s)*D0Y->at(s);
//			//double z = D0Z->at(s);
//			//if(rhoSq>=100) rhoSq=99.9;//TODO
//
//			if(!(DP0.Eta()>2.0 && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
//			if(!(DP1.Eta()>2.0 && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
//			if(!(DP.Pt()>_dptmin)) continue;
//			if(!(_dptmax==-1 || DP.Pt() < _dptmax) ) continue;
//			//if(!(D0P.Eta()>2.2&&D0P.Eta()<4.0)) continue;//TODO test restrcting the eta(D0) range (turn on in two places 2/2)
//
//			//check PID cuts
//			double weight = Dh0WEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
//			if(DP0.Pt()>25000. || DP0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
//
//			dpt        = DP.Pt();
//			deta       = DP.Eta();
//
//			//get efficiency corrections
//			double effacc(0.), effrec(0.), effsel(0.), effpid(0.), corrPt(0.);
//
//			if(dpt>=100000.) dpt=99999.;
//			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt*unitsScale,deta));
//			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt*unitsScale,deta));
//			//apply correction to the reco efficiency
//			//if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
//			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
//			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));
//
//			uint corrBin(0);
//			if(JetPT>50e3) corrBin=3;
//			else if(JetPT>30e3) corrBin=2;
//			else if(JetPT>20e3) corrBin=1;
//			else if(JetPT>15e3) corrBin=0;
//			corrPt = ptBinCorrections[corrBin]->GetBinContent(ptBinCorrections[corrBin]->FindBin(DP.Pt()));
//
//			double eff=effacc*effrec*corrPt*effsel*effpid;
//
//			if(eff<0.01) {
//				//TODO//std::cout << dpt << "\t" << deta << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
//				//if(D0TRUEIDX->at(s)>-1) {
//				//	int d = D0TRUEIDX->at(s);
//				//	if(TRUEDTRK0IDX->at(d)!=-1 && TRUEDTRK2IDX->at(d)==-1 && TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)>=5000.*5000.) {
//				//		if(TRUEDTRK0INACC->at(d)==1 && TRUEDTRK1INACC->at(d)==1) {
//				//			if(TRUEDTRK0RECO->at(d)==1 && TRUEDTRK1RECO->at(d)==1) {
//				//				std::cout << "LOST TRUE D0" << std::endl;
//				//std::cout << dpt << "\t" << deta << "\t" << eff << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
//				//std::cout << rhoSq << "\t" << z << "\t" << hrec ->GetBinContent(hrec ->FindBin(dpt/1000.,deta)) << "\t" << hcor ->GetBinContent(hcor ->FindBin(rhoSq,z)) << std::endl;
//				//			}
//				//		}
//				//	}
//				//}
//				continue;
//			}
//
//			double sbSubSwitch(0.);
//			if(TMath::Abs(DM->at(s)-_mCentre)<40.) sbSubSwitch = 1.;
//			else if(TMath::Abs(DM->at(s)-_mCentre)<80.) sbSubSwitch = -1.;
//
//			fitDrecoPT.at(4)->Fill(JetPT,weight*sbSubSwitch);
//			fitDrecoPT.at(3)->Fill(JetPT,weight*sbSubSwitch/effpid);
//			fitDrecoPT.at(2)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel));
//			fitDrecoPT.at(1)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
//			fitDrecoPT.at(0)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));
//
//			fitDtruePT.at(4)->Fill(JetTruePT,weight*sbSubSwitch);
//			fitDtruePT.at(3)->Fill(JetTruePT,weight*sbSubSwitch/effpid);
//			fitDtruePT.at(2)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel));
//			fitDtruePT.at(1)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec));
//			fitDtruePT.at(0)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));
//
//			fitD2DPT.at(4)->Fill(JetTruePT,JetPT,weight*sbSubSwitch);
//			fitD2DPT.at(3)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/effpid);
//			fitD2DPT.at(2)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel));
//			fitD2DPT.at(1)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
//			fitD2DPT.at(0)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));
//
//			dKinematics[9]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch);
//			dKinematics[7]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid));
//			dKinematics[5]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel));
//			dKinematics[3]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec));
//			dKinematics[1]->Fill(DP.Pt(),DP.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec*effacc));
//
//			//truth match to real, accepted and reconstructed D0
//			if(DTRUEIDX->at(s)<0) continue;
//			int d = DTRUEIDX->at(s);
//			if(TRUEDTRK0IDX->at(d)==-1) continue;
//			if(TRUEDTRK2IDX->at(d)!=-1) continue;
//			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;
//			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
//			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;
//
//			fndDrecoPT.at(4)->Fill(JetPT,weight);
//			fndDrecoPT.at(3)->Fill(JetPT,weight/effpid);
//			fndDrecoPT.at(2)->Fill(JetPT,weight/(effpid*effsel));
//			fndDrecoPT.at(1)->Fill(JetPT,weight/(effpid*effsel*effrec*corrPt));
//			fndDrecoPT.at(0)->Fill(JetPT,weight/(effpid*effsel*effrec*corrPt*effacc));
//
//			fndDtruePT.at(4)->Fill(JetTruePT,weight);
//			fndDtruePT.at(3)->Fill(JetTruePT,weight/effpid);
//			fndDtruePT.at(2)->Fill(JetTruePT,weight/(effpid*effsel));
//			fndDtruePT.at(1)->Fill(JetTruePT,weight/(effpid*effsel*effrec*corrPt));
//			fndDtruePT.at(0)->Fill(JetTruePT,weight/(effpid*effsel*effrec*corrPt*effacc));
//
//			fndD2DPT.at(4)->Fill(JetTruePT,JetPT,weight);
//			fndD2DPT.at(3)->Fill(JetTruePT,JetPT,weight/effpid);
//			fndD2DPT.at(2)->Fill(JetTruePT,JetPT,weight/(effpid*effsel));
//			fndD2DPT.at(1)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec*corrPt));
//			fndD2DPT.at(0)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec*corrPt*effacc));
//
//			foundDs.push_back(d*100+s);
//
//		//TODO	break;//only keep one D0 candidate per entry
//		}
//		for(unsigned int i=0; i<trueDs.size(); ++i) {
//			bool matched(false);
//			for(unsigned int j=0; j<foundDs.size(); ++j) {
//				if(foundDs[j]==trueDs[i]) {
//					matched=true;
//					break;
//				}
//			}
//			if(!matched) std::cout << progress.count()-1 << ": truD " << trueDs[i] << " has no matching fndD" << std::endl;
//		}
//		for(unsigned int i=0; i<foundDs.size(); ++i) {
//			bool matched(false);
//			for(unsigned int j=0; j<trueDs.size(); ++j) {
//				if(foundDs[i]==trueDs[j]) {
//					matched=true;
//					break;
//				}
//			}
//			if(!matched) std::cout << progress.count()-1 << ": fndD " << foundDs[i] << " has no matching truD" << std::endl;
//		}
//	}
//
//	effData->Divide(numData,denData);
//	for(uint i=1; i<_ptBins.size(); ++i) {
//		effDataJPTBins[i-1]->Divide(numDataJPTBins[i-1],denDataJPTBins[i-1]);
//		effDataJTRUEPTBins[i-1]->Divide(numDataJTRUEPTBins[i-1],denDataJTRUEPTBins[i-1]);
//	}
//
//	printf("bins of true jet pT\n");
//	for(int i=0; i<5; ++i) {
//		printf("stage %d (%s)\n",i,stages[i].Data());
//		for(int j=1; j<=truDtruePT.at(i)->GetNbinsX(); ++j) {
//			printf("%6.0f-%6.0f\t",truDtruePT.at(i)->GetBinLowEdge(j),truDtruePT.at(i)->GetBinLowEdge(j+1));
//			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f",truDtruePT.at(i)->GetBinContent(j),effDtruePT.at(i)->GetBinContent(j),fndDtruePT.at(i)->GetBinContent(j),fitDtruePT.at(i)->GetBinContent(j));
//			if(i>0) {
//				printf("\t%5.3f", effDtruePT.at(i-1)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j));
//				printf("\t%5.3f\t%5.3f", truDtruePT.at(i)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j),
//						         fndDtruePT.at(i)->GetBinContent(j)/fndDtruePT.at(i-1)->GetBinContent(j));
//				printf("\t%5.3f", (fndDtruePT.at(i)->GetBinContent(j)/fndDtruePT.at(i-1)->GetBinContent(j))/
//						  (truDtruePT.at(i)->GetBinContent(j)/truDtruePT.at(i-1)->GetBinContent(j)));
//			} else {
//				printf("\t%5.3f", (fndDtruePT.at(4)->GetBinContent(j)/fndDtruePT.at(0)->GetBinContent(j))/
//						  (truDtruePT.at(4)->GetBinContent(j)/truDtruePT.at(0)->GetBinContent(j)));
//				if(j>1 && j<truDtruePT.at(i)->GetNbinsX()) {
//				closure_true->SetBinContent(j-1,(fndDtruePT.at(4)->GetBinContent(j)/fndDtruePT.at(0)->GetBinContent(j))/
//						                (truDtruePT.at(4)->GetBinContent(j)/truDtruePT.at(0)->GetBinContent(j)));
//				}
//			}
//			printf("\n");
//		}
//	}
//
//	printf("bins of reco jet pT\n");
//	for(int i=0; i<5; ++i) {
//		printf("stage %d (%s)\n",i,stages[i].Data());
//		for(int j=1; j<=truDrecoPT.at(i)->GetNbinsX(); ++j) {
//			printf("%6.0f-%6.0f\t",truDrecoPT.at(i)->GetBinLowEdge(j),truDrecoPT.at(i)->GetBinLowEdge(j+1));
//			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f\t",truDrecoPT.at(i)->GetBinContent(j),effDrecoPT.at(i)->GetBinContent(j),fndDrecoPT.at(i)->GetBinContent(j),fitDrecoPT.at(i)->GetBinContent(j));
//			if(i>0) {
//				printf("\t%5.3f", effDrecoPT.at(i-1)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j));
//				printf("\t%5.3f\t%5.3f", truDrecoPT.at(i)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j),
//						         fndDrecoPT.at(i)->GetBinContent(j)/fndDrecoPT.at(i-1)->GetBinContent(j));
//				printf("\t%5.3f", (fndDrecoPT.at(i)->GetBinContent(j)/fndDrecoPT.at(i-1)->GetBinContent(j))/
//						  (truDrecoPT.at(i)->GetBinContent(j)/truDrecoPT.at(i-1)->GetBinContent(j)));
//			} else {
//				printf("\t%5.3f", (fndDrecoPT.at(4)->GetBinContent(j)/fndDrecoPT.at(0)->GetBinContent(j))/
//						  (truDrecoPT.at(4)->GetBinContent(j)/truDrecoPT.at(0)->GetBinContent(j)));
//				if(j>1 && j<truDrecoPT.at(i)->GetNbinsX()) {
//				closure_reco->SetBinContent(j-1,(fndDrecoPT.at(4)->GetBinContent(j)/fndDrecoPT.at(0)->GetBinContent(j))/
//						                (truDrecoPT.at(4)->GetBinContent(j)/truDrecoPT.at(0)->GetBinContent(j)));
//				}
//			}
//			printf("\n");
//		}
//	}
//
//	printf("bins of (X) true and (Y) reco jet pT\neff(found)/eff(truth)\n");
//	for(int i=1; i<5; ++i) {
//		printf("stage %d (%s)\n",i,stages[i].Data());
//		for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
//			printf("\t%3.0f-%3.0f\t",truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
//		}
//		printf("\n");
//		for(int j=1; j<=truD2DPT.at(i)->GetNbinsY(); ++j) {
//			printf("%3.0f-%3.0f",truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
//			for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
//				double ratio = (fndD2DPT.at(i)->GetBinContent(k,j)/fndD2DPT.at(i-1)->GetBinContent(k,j))/
//				               (truD2DPT.at(i)->GetBinContent(k,j)/truD2DPT.at(i-1)->GetBinContent(k,j));
//				double error = ratio*TMath::Sqrt(TMath::Power(fndD2DPT.at(i  )->GetBinError(k,j)/fndD2DPT.at(i  )->GetBinContent(k,j),2.)
//						             //+ TMath::Power(fndD02DPT.at(i-1)->GetBinError(k,j)/fndD02DPT.at(i-1)->GetBinContent(k,j),2.)
//						               + TMath::Power(truD2DPT.at(i  )->GetBinError(k,j)/truD2DPT.at(i  )->GetBinContent(k,j),2.));
//						             //+ TMath::Power(truD02DPT.at(i-1)->GetBinError(k,j)/truD02DPT.at(i-1)->GetBinContent(k,j),2.));
//				printf("\t%5.3f+/-%5.3f", ratio, error);
//			}
//			printf("\n");
//		}
//	}
//
//	printf("bins of (X) true and (Y) reco jet pT\nN(found)\n");
//	for(int i=1; i<5; ++i) {
//		printf("stage %d (%s)\n",i,stages[i].Data());
//		for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
//			printf("\t%3.0f-%3.0f\t",truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD2DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
//		}
//		printf("\n");
//		for(int j=1; j<=truD2DPT.at(i)->GetNbinsY(); ++j) {
//			printf("%3.0f-%3.0f",truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD2DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
//			for(int k=1; k<=truD2DPT.at(i)->GetNbinsX(); ++k) {
//				printf("\t%6.0f\t", fndD2DPT.at(i)->GetBinContent(k,j));
//			}
//			printf("\n");
//		}
//	}
//
//	gStyle->SetPalette(1,0);
//	const Int_t NRGBs = 5;
//	const Int_t NCont = 255;
//	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
//	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 0.00, 0.00};
//	Double_t greens[NRGBs] = { 0.00, 0.50, 1.00, 0.50, 0.00};
//	Double_t blues[NRGBs]  = { 0.00, 0.00, 1.00, 1.00, 1.00};
//	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
//	gStyle->SetNumberContours(NCont);
//	TCanvas c1;
//	c1.SetLogx();
//	dKinematics[9]->Divide(dKinematics[7]);
//	dKinematics[8]->Divide(dKinematics[6]);
//	dKinematics[9]->Divide(dKinematics[8]);
//	dKinematics[9]->SetMinimum(0.);
//	dKinematics[9]->SetMaximum(2.);
//	dKinematics[9]->Draw("colz");
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KinePIDRatio.pdf");
//	dKinematics[7]->Divide(dKinematics[5]);
//	dKinematics[6]->Divide(dKinematics[4]);
//	dKinematics[7]->Divide(dKinematics[6]);
//	dKinematics[7]->SetMinimum(0.);
//	dKinematics[7]->SetMaximum(2.);
//	dKinematics[7]->Draw("colz");
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineSELRatio.pdf");
//	dKinematics[5]->Divide(dKinematics[3]);
//	dKinematics[4]->Divide(dKinematics[2]);
//	dKinematics[5]->Divide(dKinematics[4]);
//	dKinematics[5]->SetMinimum(0.);
//	dKinematics[5]->SetMaximum(2.);
//	dKinematics[5]->Draw("colz");
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineRECRatio.pdf");
//	dKinematics[3]->Divide(dKinematics[1]);
//	dKinematics[2]->Divide(dKinematics[0]);
//	dKinematics[3]->Divide(dKinematics[2]);
//	dKinematics[3]->SetMinimum(0.);
//	dKinematics[3]->SetMaximum(2.);
//	dKinematics[3]->Draw("colz");
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineACCRatio.pdf");
//	//for(unsigned int i=0; i<dKinematics.size(); ++i) {
//	//	dKinematics[i]->Draw("colz");
//	//	TString plotname("D0Kinematics");
//	//	plotname+=i;
//	//	plotname+=".pdf";
//	//	c1.SaveAs(gSaveDir+"/"+plotname);
//	//}
//	
//	for(int i=1; i<=dKinematics[3]->GetNbinsX(); ++i) {
//		for(int j=1; j<=dKinematics[3]->GetNbinsY(); ++j) {
//			if(dKinematics[3]->GetBinContent(i,j)>0) {
//				dKinematicsPulls[0]->Fill((dKinematics[3]->GetBinContent(i,j)-1.)/dKinematics[3]->GetBinError(i,j));
//			}
//			if(dKinematics[5]->GetBinContent(i,j)>0) {
//				dKinematicsPulls[1]->Fill((dKinematics[5]->GetBinContent(i,j)-1.)/dKinematics[5]->GetBinError(i,j));
//			}
//			if(dKinematics[7]->GetBinContent(i,j)>0) {
//				dKinematicsPulls[2]->Fill((dKinematics[7]->GetBinContent(i,j)-1.)/dKinematics[7]->GetBinError(i,j));
//			}
//			if(dKinematics[9]->GetBinContent(i,j)>0) {
//				dKinematicsPulls[3]->Fill((dKinematics[9]->GetBinContent(i,j)-1.)/dKinematics[9]->GetBinError(i,j));
//			}
//		}
//	}
//	c1.SetLogx(0);
//	dKinematicsPulls[0]->Draw();
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineACCPulls.pdf");
//	dKinematicsPulls[1]->Draw();
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineRECPulls.pdf");
//	dKinematicsPulls[2]->Draw();
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KineSELPulls.pdf");
//	dKinematicsPulls[3]->Draw();
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KinePIDPulls.pdf");
//
//	const Int_t NRGBs2 = 5;
//	const Int_t NCont2 = 255;
//	Double_t stops2[NRGBs2]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
//	Double_t reds2[NRGBs2]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
//	Double_t greens2[NRGBs2] = { 1.00, 0.95, 0.50, 0.00, 0.00};
//	Double_t blues2[NRGBs2]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
//	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
//	gStyle->SetNumberContours(NCont2);
//
//	c1.Clear();
//	c1.Divide(4,5);
//	for(int i=0; i<5; ++i) {
//		double max = 1.1*TMath::Max(TMath::Max(truD2DPT.at(i)->GetMaximum(),effD2DPT.at(i)->GetMaximum()),TMath::Max(fndD2DPT.at(i)->GetMaximum(),fitD2DPT.at(i)->GetMaximum()));
//		truD2DPT.at(i)->SetMinimum(0.);
//		truD2DPT.at(i)->SetMaximum(max);
//		fndD2DPT.at(i)->SetMinimum(0.);
//		fndD2DPT.at(i)->SetMaximum(max);
//		fitD2DPT.at(i)->SetMinimum(0.);
//		fitD2DPT.at(i)->SetMaximum(max);
//		effD2DPT.at(i)->SetMinimum(0.);
//		effD2DPT.at(i)->SetMaximum(max);
//		c1.cd(4*i+1);
//		truD2DPT.at(i)->Draw("colz");
//		c1.cd(4*i+2);
//		fndD2DPT.at(i)->Draw("colz");
//		c1.cd(4*i+3);
//		fitD2DPT.at(i)->Draw("colz");
//		c1.cd(4*i+4);
//		effD2DPT.at(i)->Draw("colz");
//	}
//	c1.SaveAs(gSaveDir+"/"+_whichD+"2DEffMaps.pdf");
//
//	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
//	gStyle->SetNumberContours(NCont);
//	c1.Clear();
//	c1.Divide(5,3);
//	for(int i=0; i<5; ++i) {
//		fndD2DPT.at(i)->Divide(truD2DPT.at(i));
//		fitD2DPT.at(i)->Divide(truD2DPT.at(i));
//		effD2DPT.at(i)->Divide(truD2DPT.at(i));
//		double max = 2.0;
//		fndD2DPT.at(i)->SetMinimum(0.);
//		fndD2DPT.at(i)->SetMaximum(max);
//		fitD2DPT.at(i)->SetMinimum(0.);
//		fitD2DPT.at(i)->SetMaximum(max);
//		effD2DPT.at(i)->SetMinimum(0.);
//		effD2DPT.at(i)->SetMaximum(max);
//		c1.cd(i+1);
//		fndD2DPT.at(i)->Draw("colz");
//		c1.cd(i+6);
//		fitD2DPT.at(i)->Draw("colz");
//		c1.cd(i+11);
//		effD2DPT.at(i)->Draw("colz");
//	}
//	c1.SaveAs(gSaveDir+"/2DEffMapsRatios.pdf");
//
//	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
//	gStyle->SetNumberContours(NCont2);
//	c1.Clear();
//	c1.Divide(6,7);
//	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
//		c1.cd(i+1);
//		dKinematicsTRBinned.at(i)->Draw("colz");
//	}
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin.pdf");
//	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
//		c1.cd(i+1);
//		dKinematicsTRBinned_recD_tru.at(i)->Draw("colz");
//	}
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin_recD_truKine.pdf");
//	for(uint i=0; i<(ptBinsWithOF.size()-1)*ptBinsWithOF.size(); ++i) {
//		c1.cd(i+1);
//		dKinematicsTRBinned_recD_rec.at(i)->Draw("colz");
//	}
//	c1.SaveAs(gSaveDir+"/"+_whichD+"KinematicsByJPTBin_recD_recKine.pdf");
//
//	TString closureOutName = gSaveDir+"/effClosure";
//	closureOutName+=flavour;
//	if(_useSimpleEffs) closureOutName+="_simpleEffs";
//	closureOutName+=".root";
//	TFile* fout = TFile::Open(closureOutName,"RECREATE");
//	closure_true->Write();
//	closure_reco->Write();
//	effModel->Write();
//	effData->Write();
//	numData->Write();
//	denData->Write();
//	for(uint i=0; i<_ptBins.size()-1; ++i) {
//		effDataJPTBins[i]->Write();
//		numDataJPTBins[i]->Write();
//		denDataJPTBins[i]->Write();
//		effDataJTRUEPTBins[i]->Write();
//		numDataJTRUEPTBins[i]->Write();
//		denDataJTRUEPTBins[i]->Write();
//	}
//	fout->Close();
//
//	dm->reset();
//
//	if(!_useSimpleEffs) {
//		//rerun with simple effs
//		_useSimpleEffs=true;
//		testEffs(flavour);
//		_useSimpleEffs=false;
//	}

	return true;
}

bool SimDFitter::fitD(double& yield4, double& error4, double& yield5, double& error5, int binPT, int binY, uint effType) {
	if(!_inputsSet) return false;
	if(!_effsSet) {
		std::cout << "INFO in SimDFitter::fitD: D efficiencies not calculated yet" << std::endl;
		std::cout << "                          running now" << std::endl;
		if(!addEffs()) return false;
	}

	//find the pT and y limits of the requested bin
	double ptmin(0.), ptmax(100e3), ymin(0.), ymax(10.);
	if(binPT<0 || binPT>=_ptBins.size()) {
		ptmin = _ptBins[0];
		ptmax = _ptBins[_ptBins.size()-1];
	} else {
		ptmin = _ptBins[binPT];
		ptmax = _ptBins[binPT+1];
	}
	if(binY<0 || binY>=_yBins.size()) {
		ymin = _yBins[0];
		ymax = _yBins[_yBins.size()-1];
	} else {
		ymin = _yBins[binY];
		ymax = _yBins[binY+1];
	}

	//check if we've changed bin (need to recreate D pT bins)
	if(_jptmin != ptmin || _jptmax != ptmax ||
	   _ymin != ymin || _ymax != ymax) {
		_jptmin = ptmin;
		_jptmax = ptmax;
		_ymin = ymin;
		_ymax = ymax;
		_initialised = false;
	}

	if(effType==4) {
		_flavour = 4;
		_useEffs = true;
	} else if(effType==5) {
		_flavour = 5;
		_useEffs = true;
	} else {
		_useEffs = false;
	}
	
	if(!doFits()) return false;
	std::cout << _yield4 << "\t" << _error4 << "\t" << _yield5 << "\t" << _error5 << std::endl;//TODO

	yield5 = _yield5;
	error5 = _error5;
	yield4 = _yield4;
	error4 = _error4;

	return true;
}

//bool SimDFitter::makeInput(int flav) {
//	TString fname = "/data/dijets/dijets_";
//	if(flav==0) {
//		fname+="2016.root";
//	} else {
//		fname+="sim";
//		fname+=flav;
//		fname+=".root";
//	}
//
//	TFile* f = TFile::Open(fname);
//	if(!f) return false;
//
//	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
//	if(!t) return false;
//
//	std::vector<double> *vD0M = new std::vector<double>();
//	std::vector<double> *vD0IPCHI2 = new std::vector<double>();
//	std::vector<double> *vD0PT = new std::vector<double>();
//	std::vector<double> *vD0PX = new std::vector<double>();
//	std::vector<double> *vD0PY = new std::vector<double>();
//	std::vector<double> *vD0PZ = new std::vector<double>();
//	std::vector<double> *vD0E = new std::vector<double>();
//	std::vector<double> *vD0X = new std::vector<double>();
//	std::vector<double> *vD0Y = new std::vector<double>();
//	std::vector<double> *vD0Z = new std::vector<double>();
//	std::vector<double> *vD0KP = new std::vector<double>();
//	std::vector<double> *vD0KPT = new std::vector<double>();
//	std::vector<double> *vD0KPX = new std::vector<double>();
//	std::vector<double> *vD0KPY = new std::vector<double>();
//	std::vector<double> *vD0KPZ = new std::vector<double>();
//	std::vector<double> *vD0PIP = new std::vector<double>();
//	std::vector<double> *vD0PIPT = new std::vector<double>();
//	std::vector<double> *vD0PIPX = new std::vector<double>();
//	std::vector<double> *vD0PIPY = new std::vector<double>();
//	std::vector<double> *vD0PIPZ = new std::vector<double>();
//	std::vector<double> *vD0KPNNK = new std::vector<double>();
//	std::vector<double> *vD0PIPNNPI = new std::vector<double>();
//	std::vector<double> *vD0KWEIGHT = new std::vector<double>();
//	std::vector<double> *vD0PIWEIGHT = new std::vector<double>();
//	std::vector<double> *vD0TRUEIDX = new std::vector<double>();
//	std::vector<double> *vD0FROMB = new std::vector<double>();
//
//	double JetPT;
//	double JetEta;
//	double JetTruePT;
//
//	t->SetBranchAddress("D0M",           &vD0M);
//	t->SetBranchAddress("D0IPCHI2",      &vD0IPCHI2);
//	t->SetBranchAddress("D0PT",          &vD0PT);
//	t->SetBranchAddress("D0PX",          &vD0PX);
//	t->SetBranchAddress("D0PY",          &vD0PY);
//	t->SetBranchAddress("D0PZ",          &vD0PZ);
//	t->SetBranchAddress("D0E",           &vD0E);
//	t->SetBranchAddress("D0X",           &vD0X);
//	t->SetBranchAddress("D0Y",           &vD0Y);
//	t->SetBranchAddress("D0Z",           &vD0Z);
//	t->SetBranchAddress("D0KP",          &vD0KP);
//	t->SetBranchAddress("D0KPT",         &vD0KPT);
//	t->SetBranchAddress("D0KPX",         &vD0KPX);
//	t->SetBranchAddress("D0KPY",         &vD0KPY);
//	t->SetBranchAddress("D0KPZ",         &vD0KPZ);
//	t->SetBranchAddress("D0PIP",         &vD0PIP);
//	t->SetBranchAddress("D0PIPT",        &vD0PIPT);
//	t->SetBranchAddress("D0PIPX",        &vD0PIPX);
//	t->SetBranchAddress("D0PIPY",        &vD0PIPY);
//	t->SetBranchAddress("D0PIPZ",        &vD0PIPZ);
//	t->SetBranchAddress("D0KPNNK",       &vD0KPNNK);
//	t->SetBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);
//	t->SetBranchAddress("D0KWEIGHT",     &vD0KWEIGHT);
//	t->SetBranchAddress("D0PIWEIGHT",    &vD0PIWEIGHT);
//	if(flav!=0) {
//		t->SetBranchAddress("D0TRUEIDX",     &vD0TRUEIDX);
//		t->SetBranchAddress("D0FROMB",       &vD0FROMB);
//	}
//
//	t->SetBranchAddress("JetPT",         &JetPT);
//	t->SetBranchAddress("JetEta",        &JetEta);
//	if(flav!=0) {
//		t->SetBranchAddress("JetTruePT",     &JetTruePT);
//	}
//
//	unsigned int nentries0 = t->GetEntries();
//
//	TFile* fout = TFile::Open("d0InputDatasets.root","UPDATE");
//	TString tname = "T";
//	tname+=flav;
//	TTree* tout = new TTree(tname,"");
//
//	double D0M(0.), D0PT(0.), D0Eta(0.), D0LogIPChi2(0.);
//
//	tout->Branch("JetPT",        &JetPT);
//	tout->Branch("JetEta",       &JetEta);
//	tout->Branch("JetTruePT",    &JetTruePT);
//	tout->Branch("D0M",          &D0M);
//	tout->Branch("D0PT",         &D0PT);
//	tout->Branch("D0Eta",        &D0Eta);
//	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
//
//	//first make non-vector tree for fits
//	boost::progress_display progress( nentries0 );
//	//for(unsigned int ientry=0; ientry<nentries0; ++ientry)
//	for(uint i=0; i<nentries0; ++i) {
//		t->GetEntry(i);
//		++progress;
//		//t0->GetEntry(ientry);
//
//		for(unsigned int s=0; s<vD0M->size(); ++s) {
//			//first check in our tight acceptance
//			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
//			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
//			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));
//
//			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
//			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
//			//if(!(D0P.Pt()>_dptmin)) continue;
//			//if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;
//			if(flav!=0 && vD0TRUEIDX->at(s)==-1) continue;
//			if(flav==5 && vD0FROMB->at(s)==0) continue;
//			if(flav==4 && vD0FROMB->at(s)==1) continue;
//			if(flav==0 && vD0M->at(s)>1844. && vD0M->at(s)<1884.) continue;
//
//			//check PID cuts
//			//if(!_dataIsMC && vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
//			//if(!dataIsMC && vD0PIPNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off
//
//			D0M         = vD0M->at(s);
//			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
//			D0PT        = D0P.Pt();
//			D0Eta       = D0P.Eta();
//
//			//for the combinatorial, we need to resample the mass to remove signal window
//			if(flav==0) D0M = gRandom->Uniform(1784.,1944.);
//
//			tout->Fill();
//			break;//only keep one D0 candidate per entry
//		}
//	}
//
//	tout->AutoSave();
//	fout->Close();
//
//	return true;
//}

bool SimDFitter::loadDatasets() {
	TFile* f = TFile::Open(dFileName());
	//if(!f) {
	//	if(!makeInput(0)) return false;
	//	if(!makeInput(4)) return false;
	//	if(!makeInput(5)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//}
	if(!f) return false;

	TTree* t0 = dynamic_cast<TTree*>(f->Get("T0"));
	//if(!t0) {
	//	if(!makeInput(0)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t0 = dynamic_cast<TTree*>(f->Get("T0"));
	//}
	TTree* t4 = dynamic_cast<TTree*>(f->Get("T4"));
	//if(!t4) {
	//	if(!makeInput(4)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t4 = dynamic_cast<TTree*>(f->Get("T4"));
	//}
	TTree* t5 = dynamic_cast<TTree*>(f->Get("T5"));
	//if(!t5) {
	//	if(!makeInput(5)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t5 = dynamic_cast<TTree*>(f->Get("T5"));
	//}
	TTree* t7 = dynamic_cast<TTree*>(f->Get("T7"));
	if(!t0 || !t4 || !t5 || !t7) return false;

	std::cout << "?" << t7->GetEntries() << std::endl;//TODO
	//setup equal bins in pT(D)
	if(!_useFixedDBins) setupDPtBins(t7);

	RooRealVar DM        (""+_whichD+"M",         ""+_whichD+"M",       _mLow, _mHigh, "MeV/#it{c}^{2}"); 
	RooRealVar DLOGIPCHI2(""+_whichD+"LogIPChi2", ""+_whichD+"LogIPChi2", -5.,   15., ""); 
	RooRealVar DPT       (""+_whichD+"PT",        ""+_whichD+"PT",         0.,   1e5, "MeV/#it{c}"); 
	RooRealVar JetPT      ("JetPT",       "JetPT",        0.,   1e5, "MeV/#it{c}"); 
	RooRealVar JetEta     ("JetEta",      "JetEta",       0.,   10., ""); 
	RooRealVar weight4    ("weight4",     "weight4",      0.,   30., ""); 
	RooRealVar weight5    ("weight5",     "weight5",      0.,   30., ""); 
	RooRealVar enhanced   ("enhanced",    "enhanced",    -1.,    6., "");

	RooThresholdCategory* ptCat(0);

	TString binName="PTbin";
	binName+=_ptBinsD.size()-2;
	if(!_usePtFracBins) {
		ptCat = new RooThresholdCategory("ptCat","ptCat", DPT, binName, _ptBinsD.size()-2);
		for(uint i=1; i<_ptBinsD.size()-1; ++i) {
			binName="PTbin";
			binName+=i-1;
			ptCat->addThreshold(_ptBinsD[i],binName.Data(),i-1);
		}
	} else {
		RooFormulaVar* DPTFrac = new RooFormulaVar(""+_whichD+"PTFrac","","@0/@1",RooArgList(DPT,JetPT));
		ptCat = new RooThresholdCategory("ptCat","ptCat", *DPTFrac, binName, _nDPtBins-1);
		for(uint i=1; i<_nDPtBins; ++i) {
			binName="PTbin";
			binName+=i-1;
			ptCat->addThreshold(static_cast<double>(i)/_nDPtBins,binName.Data(),i-1);
		}
	}

	binName="JetPTbin";
	binName+=_ptBins.size()-2;
	RooThresholdCategory jetptCat("jetptCat","jetptCat", JetPT, binName, _ptBins.size()-2);
	for(uint i=1; i<_ptBins.size()-1; ++i) {
		binName="JetPTbin";
		binName+=i-1;
		jetptCat.addThreshold(_ptBins[i],binName.Data(),i-1);
	}

	binName="JetYbin";
	binName+=_yBins.size()-2;
	RooThresholdCategory jetyCat("jetyCat","jetyCat", JetEta, binName, _yBins.size()-2);
	for(uint i=1; i<_yBins.size()-1; ++i) {
		binName="JetYbin";
		binName+=i-1;
		jetyCat.addThreshold(_yBins[i],binName.Data(),i-1);
	}

	RooThresholdCategory enhancedCat("enhancedCat","enhancedCat", enhanced, "enhancedB", 2);
	enhancedCat.addThreshold(3.5, "unenhanced", 0);
	enhancedCat.addThreshold(4.5, "enhancedC", 1);

	RooCategory typeCat("typeCat","typeCat");
	typeCat.defineType("data",3);
	typeCat.defineType("charm",1);
	typeCat.defineType("beauty",2);
	typeCat.defineType("sideband",0);

	RooCmdArg sidebandCutCmd   = RooCmdArg::none();
	if(_combShapeSyst) {
		TString sbCutStr = _whichD+"M<";
		sbCutStr+=_mSBLooseLow;
		sbCutStr+=" || "+_whichD+"M>";
		sbCutStr+=_mSBLooseHigh;
		sidebandCutCmd = RooFit::Cut(sbCutStr);
	}

	RooDataSet cData("cData", "cData", RooArgList(DLOGIPCHI2,DM,DPT,JetPT,JetEta,weight4,weight5,enhanced), RooFit::Import(*t4));
	RooDataSet bData("bData", "bData", RooArgList(DLOGIPCHI2,DM,DPT,JetPT,JetEta,weight4,weight5,enhanced), RooFit::Import(*t5));
	RooDataSet sData("sData", "sData", RooArgList(DLOGIPCHI2,DM,DPT,JetPT,JetEta,weight4,weight5,enhanced), RooFit::Import(*t0), sidebandCutCmd);
	RooDataSet dData("dData", "dData", RooArgList(DLOGIPCHI2,DM,DPT,JetPT,JetEta,weight4,weight5,enhanced), RooFit::Import(*t7));
	RooDataSet data("data", "data", RooArgList(DLOGIPCHI2,DM,DPT,JetPT,JetEta,weight4,weight5,enhanced), RooFit::Index(typeCat), RooFit::Import("charm",cData), RooFit::Import("beauty",bData), RooFit::Import("sideband",sData), RooFit::Import("data",dData));
	data.addColumn(*ptCat);
	data.addColumn(jetptCat);
	data.addColumn(jetyCat);
	data.addColumn(enhancedCat);
	sData.addColumn(*ptCat);
	dData.addColumn(*ptCat);

	ws->import(cData);
	ws->import(bData);
	ws->import(sData);
	ws->import(dData);
	ws->import(data);
	nEntries = data.sumEntries();

	return true;
}

bool SimDFitter::setupDPtBins(TTree* t) {
	//TODO
	//_ptBinsD.clear();
	//_ptBinsD.push_back(0);
	//_ptBinsD.push_back(7e3);
	//_ptBinsD.push_back(10e3);
	//_ptBinsD.push_back(15e3);
	//_ptBinsD.push_back(20e3);
	//_ptBinsD.push_back(100e3);
	//return true;
	//TODO
	uint nbins(500);
	TH1D pt_pdf("pt_pdf","",nbins,0.,100e3);
	//TH1D pt_cdf("pt_cdf","",nbins,0.,100e3);

	TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f",_jptmin,_jptmax,_ymin,_ymax);
	t->Draw(""+_whichD+"PT>>pt_pdf",cutStr);
	double total = pt_pdf.GetEntries();
	double sum(0.);
	double step = total/_nDPtBins;
	std::cout << _nDPtBins << "?" << std::endl;//TODO

	_ptBinsD.clear();
	_ptBinsD.push_back(0.);
	for(uint iBin=1; iBin<=nbins; ++iBin) {
		sum += pt_pdf.GetBinContent(iBin);
		if(sum >= _ptBinsD.size()*step) {
			std::cout << iBin << " " << pt_pdf.GetBinLowEdge(iBin) << " " << sum << " " << total << " " << step << " " << _ptBinsD.size()*step << std::endl;//TODO
			_ptBinsD.push_back(pt_pdf.GetBinLowEdge(iBin+1));
		}
		//pt_cdf.SetBinContent(i,sum);
	}
	//pt_cdf.Scale(1./sum);
	_ptBinsD.back() = 1e5;//TODO

	//TCanvas c;
	//pt_cdf.Draw();
	//c.SaveAs("test.pdf");//TODO
	for(uint iBin=0; iBin<_ptBinsD.size(); ++iBin) {
		std::cout << _ptBinsD[iBin] << std::endl;
	}

	return true;
}

bool SimDFitter::setupModel() {
	RooRealVar* DM = ws->var(""+_whichD+"M");
	RooRealVar* DLOGIPCHI2 = ws->var(""+_whichD+"LogIPChi2");

	//mass part
	RooRealVar dMean("dMean","$\\mu_m$",valMassMean,_mSBLow,_mSBHigh);
	RooRealVar dWidthBase("dWidthBase","$\\sigma_C$",valMassWidth,4.,40.);
	RooRealVar dWidthScale("dWidthScale","$\\sigma_C^{\\rm data}/\\sigma_C^{\\rm sim}$",1.,0.5,2.);
	RooFormulaVar dWidth("dWidth","","@0*@1",RooArgList(dWidthBase,dWidthScale));

	RooRealVar dRatio("dRatio","$r$", valMassRatio, 1.0, 5);
	RooFormulaVar dWidthG("dWidthG","","@0*@1",RooArgList(dWidth,dRatio));
	RooRealVar dAlpha("dAlpha", "$\\alpha$", valMassAlpha, 0.5, 5.0);
	RooRealVar dN(    "dN",     "$N$",     valMassN);
	RooRealVar dFrac("dFrac","$f$", valMassFrac, 0., 1.);

	RooGaussian signalMass_gauss("signalMass_gauss","",*DM,dMean,dWidthG);
	RooCBShape  signalMass_cb("signalMass_cb", "", *DM, dMean, dWidth, dAlpha, dN); 
	RooAddPdf signalMass( "signalMass", "", RooArgList(signalMass_cb,signalMass_gauss), RooArgList(dFrac) );

	RooAbsReal* p0;
	if(_allowLinearBkgMassShift) {
		RooRealVar* p0offset = new RooRealVar("p0offset","offsetp0",-4., -10., 10.);
		RooRealVar* p0grad   = new RooRealVar("p0grad","gradp0",0., -10., 10.);
		RooRealVar* avePt     = new RooRealVar("avePt","avePt",0.);
		RooRealVar* scalep0   = new RooRealVar("scalep0","scalep0",1e-4);
		p0 = new RooFormulaVar("p0","$p_0$","(@0*@1 + @2)*@3", RooArgList(*avePt,*p0grad,*p0offset,*scalep0));
	} else {
		//RooRealVar p0("p0","$p_0$",-4e-4, -1, 1);
		p0 = new RooRealVar("p0","$p_0$",0., -1, 10);//works
	}
	RooPolynomial bkgrndMass("bkgrndMass","",*DM, RooArgList(*p0));

	//IP part
	RooRealVar promptMeanBase("promptMeanBase",  "$\\mu^{\\rm prompt}$",     valPromptMean,0.5,1.);
	RooRealVar promptMeanShift("promptMeanShift","$\\mu^{\\rm prompt data}-\\mu^{\\rm prompt sim}$",0.,-1.,1.);
	RooFormulaVar promptMean("promptMean","","@0+@1",RooArgList(promptMeanBase,promptMeanShift));
	RooRealVar promptWidthBase("promptWidthBase",  "$\\sigma^{\\rm prompt}$",  valPromptWidth,0.5,1.5);
	RooRealVar promptWidthScale("promptWidthScale","$\\sigma^{\\rm prompt data}/\\sigma^{\\rm prompt sim}$",1.,0.2,2.);
	RooFormulaVar promptWidth("promptWidth","","@0*@1",RooArgList(promptWidthBase,promptWidthScale));
	RooRealVar promptAsym("promptAsym",          "$\\epsilon^{\\rm prompt}$",valPromptAsym,-0.5,0.);
	RooRealVar promptRhoL("promptRhoL",          "$\\rho_L^{\\rm prompt}$",  valPromptRhoL,0.5,2.);
	RooRealVar promptRhoR("promptRhoR",          "$\\rho_R^{\\rm prompt}$",  valPromptRhoR);//,0.5,2.);
	RooPromptShape promptIP("promptIP","",*DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displcMeanBase("displcMeanBase",  "$\\mu^{\\rm displc}$",     valDisplMean,2.,12.);
	RooRealVar displcMeanShift("displcMeanShift","$\\mu^{\\rm displc data}-\\mu^{\\rm displc sim}$",0.,-1.,1.);
	RooFormulaVar displcMean("displcMean","","@0+@1",RooArgList(displcMeanBase,displcMeanShift));
	RooRealVar displcWidthBase("displcWidthBase","$\\sigma^{\\rm displc}$",  valDisplWidth,0.5,6.);
	RooRealVar displcWidthScale("displcWidthScale","$\\sigma^{\\rm displc data}/\\sigma^{\\rm displc sim}$",1.,0.5,2.);
	RooFormulaVar displcWidth("displcWidth","","@0*@1",RooArgList(displcWidthBase,displcWidthScale));
	RooRealVar displcAsym("displcAsym",  "$\\epsilon^{\\rm displc}$",valDisplAsym,-1.,1.);
	RooRealVar displcRhoL("displcRhoL",  "$\\rho_L^{\\rm displc}$",  valDisplRhoL);
	RooRealVar displcRhoR("displcRhoR",  "$\\rho_R^{\\rm displc}$",  valDisplRhoR);
	RooPromptShape displcIP("displcIP","",*DLOGIPCHI2,displcMean,displcWidth,displcAsym,displcRhoL,displcRhoR);

	RooAbsPdf* bkgrndIP(0);
	RooRealVar bkgrndBinIdx("bkgrndBinIdx","pt bin index",0.,-1.,_ptBinsD.size());
	RooArgList pdfs("pdfs");

	if(!_splitBkgIPShape) {
		bkgrndIP = new RooKeysPdf("bkgrndIP","",*DLOGIPCHI2,*dynamic_cast<RooDataSet*>(ws->data("sData")));
	} else {
		RooDataSet* sData = dynamic_cast<RooDataSet*>(ws->data("sData"));
		const RooArgSet* vars = sData->get();
		for(uint iBin=0; iBin<_ptBinsD.size()-1; ++iBin) {
			TString nameStr = TString::Format("_%d",iBin);
			TString cutStr = TString::Format("ptCat==%d",iBin);
			RooDataSet* ds = new RooDataSet("sData"+nameStr, "", sData, *vars, cutStr);
			RooKeysPdf* pdf = new RooKeysPdf("bkgrndIP"+nameStr, "", *DLOGIPCHI2, *ds);
			pdfs.add(*pdf);
		}
		bkgrndIP = new RooMultiKeysPdf("bkgrndIP", "", *DLOGIPCHI2, bkgrndBinIdx, pdfs);
	}

	RooProdPdf prompt("prompt","",RooArgSet(signalMass,promptIP));
	RooProdPdf displc ("displc", "",RooArgSet(signalMass,displcIP));
	RooProdPdf bkgrnd("bkgrnd","",RooArgSet(bkgrndMass,*bkgrndIP));

	// -- yields
	RooRealVar signalYield( "signalYield", "$N_{\\rm sig}$",   -2.,  nEntries);
	RooRealVar promptYield( "promptYield", "$N_{\\rm prmpt}$", -2.,  nEntries);
	RooRealVar displcYield( "displcYield", "$N_{\\rm displ}$", -2.,  nEntries);
	RooRealVar bkgrndYield( "bkgrndYield", "$N_{\\rm bkgnd}$", -2.,  nEntries);

	RooAddPdf mass_pdf( "mass_pdf",  "mass_pdf", RooArgList(signalMass,bkgrndMass), RooArgList(signalYield,bkgrndYield) );
	RooAddPdf ip_pdf(   "ip_pdf",    "ip_pdf",   RooArgList(promptIP,displcIP,*bkgrndIP), RooArgList(promptYield,displcYield,bkgrndYield) );
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displc,bkgrnd), RooArgList(promptYield,displcYield,bkgrndYield) );

	ws->importClassCode(RooPromptShape::Class());
	ws->importClassCode(RooKeysPdf::Class());
	ws->importClassCode(RooMultiKeysPdf::Class());
	ws->import(mass_pdf);
	ws->import(ip_pdf);
	ws->import(data_pdf,RooFit::RenameConflictNodes("_2d"));//RooFit::RecycleConflictNodes());

	TString massFitSplitPars("dWidthBase,dRatio,dAlpha,dFrac");
	TString ipFitSplitPars("promptMeanBase,promptWidthBase,promptAsym,promptRhoL,displcMeanBase,displcWidthBase,displcAsym");

	if(_splitMassWidthScale) massFitSplitPars += ",dWidthScale";
	if(_splitPromptMeanShift) ipFitSplitPars += ",promptMeanShift";
	if(_splitPromptWidthScale) ipFitSplitPars += ",promptWidthScale";
	if(_splitDisplcMeanShift) ipFitSplitPars += ",displcMeanShift";
	if(_splitDisplcWidthScale) ipFitSplitPars += ",displcWidihScale";
	if(_splitBkgMassShape && !_allowLinearBkgMassShift) massFitSplitPars += ",p0";
	if(_allowLinearBkgMassShift) massFitSplitPars += ",avePt";
	if(_splitBkgIPShape) ipFitSplitPars += ",bkgrndBinIdx";

	RooSimWSTool sct(*ws);
	sct.build("sim_mass_pdf","mass_pdf",
	  	  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	  	  RooFit::SplitParam(massFitSplitPars+",signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	sct.build("sim_ip_pdf","ip_pdf",
//TODO//20200420	  	  RooFit::SplitParam("promptMean,promptWidthBase,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,bkgrndBinIdx,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	  	  RooFit::SplitParam(ipFitSplitPars+",promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	sct.build("sim_data_pdf","data_pdf",
		  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,promptMean,promptWidthBase,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//TODO//20200420		  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,promptMean,promptWidthBase,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,bkgrndBinIdx,promptYield,displcYield,bkgrndYield,p0,dMean","ptCat"));//,jetptCat,jetyCat"));
		  RooFit::SplitParam(massFitSplitPars+","+ipFitSplitPars+",promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));

	if(_splitBkgIPShape) {
		for(uint iBin=0; iBin<_ptBinsD.size()-1; ++iBin) {
			RooRealVar* binIdx = ws->var(TString::Format("bkgrndBinIdx_PTbin%d",iBin));
			binIdx->setConstant(true);
			binIdx->setVal(iBin);
		}
	}
	if(_allowLinearBkgMassShift) {
		RooDataSet* dData = dynamic_cast<RooDataSet*>(ws->data("dData"));
		RooRealVar* DPT  = dynamic_cast<RooRealVar*>(ws->var("DPT"));
		for(uint iBin=0; iBin<_ptBinsD.size()-1; ++iBin) {
			RooRealVar* avePt = ws->var(TString::Format("avePt_PTbin%d",iBin));
			//TString cutStr = TString::Format("ptCat==%d",iBin);
			TString cutStr = TString::Format("ptCat==%d && JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f",iBin,_jptmin,_jptmax,_ymin,_ymax);
			avePt->setConstant(true);
			avePt->setVal(dData->mean(*DPT,cutStr)/1e5);
			std::cout << cutStr << " " << avePt->getVal() << std::endl;//TODO
		}
	}

	//ws->Print();//TODO
	//TODO//initPars();
	//return false;

	return true;
}

//bool SimDFitter::setupDPtBins() {
//	
//	TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f",_jptmin,_jptmax,_ymin,_ymax);
//	TString nameStr = TString::Format("_%.1f-%.1f_%.3f-%.3f",_jptmin,_jptmax,_ymin,_ymax);
//
//	RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
//	RooRealVar* D0PT = ws->var("D0PT");
//	const RooArgSet* vars = data->get();
//	RooDataSet* ds = new RooDataSet("data"+nameStr, "", data, *vars, cutStr);
//
//	TString binName="PTbin";
//	binName+=_ptBinsD.size()-2;
//	RooThresholdCategory ptCat("ptCat"+nameStr,"ptCat", *D0PT, binName, _ptBinsD.size()-2);
//	for(uint i=1; i<_ptBinsD.size()-1; ++i) {
//		binName="PTbin";
//		binName+=i-1;
//		ptCat.addThreshold(_ptBinsD[i],binName.Data(),i-1);
//	}
//
//	ds->addColumn(ptCat);
//	RooSimWSTool sct(*ws);
//	sct.build("sim_mass_pdf"+nameStr,"mass_pdf",
//	  	  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//	  	  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,signalYield,bkgrndYield","ptCat"+nameStr));
//	sct.build("sim_ip_pdf"+nameStr,"ip_pdf",
//	  	  RooFit::SplitParam("promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"+nameStr));
//	sct.build("sim_data_pdf"+nameStr,"data_pdf",
//		  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//		  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"+nameStr));
//
//	ws->Print();
//	//TODO//initPars();
//
//	return true;
//}

void SimDFitter::initPars() {
//signal mass
	ws->var(         "dAlpha_PTbin0"  )->setVal(  1.9180e+00);// +/-  6.72e-02
	ws->var(         "dAlpha_PTbin1"  )->setVal(  1.8576e+00);// +/-  4.03e-02
	ws->var(         "dAlpha_PTbin2"  )->setVal(  2.0967e+00);// +/-  6.33e-02
	ws->var(         "dAlpha_PTbin3"  )->setVal(  2.0526e+00);// +/-  8.37e-02
	ws->var(         "dAlpha_PTbin4"  )->setVal(  2.1923e+00);// +/-  1.05e-01
	ws->var(          "dFrac_PTbin0"  )->setVal(  9.3494e-01);// +/-  1.24e-02
	ws->var(          "dFrac_PTbin1"  )->setVal(  9.4841e-01);// +/-  8.24e-03
	ws->var(          "dFrac_PTbin2"  )->setVal(  9.2646e-01);// +/-  8.37e-03
	ws->var(          "dFrac_PTbin3"  )->setVal(  8.8369e-01);// +/-  1.61e-02
	ws->var(          "dFrac_PTbin4"  )->setVal(  7.9370e-01);// +/-  1.95e-02
	ws->var(          "dMean_PTbin0"  )->setVal(  1.8656e+03);// +/-  1.44e-01
	ws->var(          "dMean_PTbin1"  )->setVal(  1.8655e+03);// +/-  1.01e-01
	ws->var(          "dMean_PTbin2"  )->setVal(  1.8655e+03);// +/-  1.03e-01
	ws->var(          "dMean_PTbin3"  )->setVal(  1.8654e+03);// +/-  1.47e-01
	ws->var(          "dMean_PTbin4"  )->setVal(  1.8655e+03);// +/-  1.43e-01
	ws->var(         "dRatio_PTbin0"  )->setVal(  3.5897e+00);// +/-  3.58e-01
	ws->var(         "dRatio_PTbin1"  )->setVal(  3.4258e+00);// +/-  3.14e-01
	ws->var(         "dRatio_PTbin2"  )->setVal(  3.5701e+00);// +/-  2.29e-01
	ws->var(         "dRatio_PTbin3"  )->setVal(  2.8388e+00);// +/-  1.73e-01
	ws->var(         "dRatio_PTbin4"  )->setVal(  2.4915e+00);// +/-  8.27e-02
	ws->var(         "dWidth_PTbin0"  )->setVal(  7.5904e+00);// +/-  1.47e-01
	ws->var(         "dWidth_PTbin1"  )->setVal(  8.0534e+00);// +/-  1.09e-01
	ws->var(         "dWidth_PTbin2"  )->setVal(  8.9742e+00);// +/-  1.11e-01
	ws->var(         "dWidth_PTbin3"  )->setVal(  9.5189e+00);// +/-  1.69e-01
	ws->var(         "dWidth_PTbin4"  )->setVal(  1.1444e+01);// +/-  2.16e-01
	//ws->var(    "signalYield_PTbin0"  )->setVal(  1.4054e+04);// +/-  2.28e+02
	//ws->var(    "signalYield_PTbin1"  )->setVal(  2.7708e+04);// +/-  3.00e+02
	//ws->var(    "signalYield_PTbin2"  )->setVal(  3.2729e+04);// +/-  3.19e+02
	//ws->var(    "signalYield_PTbin3"  )->setVal(  1.8350e+04);// +/-  2.31e+02
	//ws->var(    "signalYield_PTbin4"  )->setVal(  3.6071e+04);// +/-  3.39e+02
//bkgrnd mass
	//ws->var(    "bkgrndYield_PTbin0"  )->setVal(  4.8410e+03);// +/-  1.31e+02
	//ws->var(    "bkgrndYield_PTbin1"  )->setVal(  7.6104e+03);// +/-  1.55e+02
	//ws->var(    "bkgrndYield_PTbin2"  )->setVal(  7.6209e+03);// +/-  1.50e+02
	//ws->var(    "bkgrndYield_PTbin3"  )->setVal(  4.0486e+03);// +/-  1.09e+02
	//ws->var(    "bkgrndYield_PTbin4"  )->setVal(  5.3398e+03);// +/-  1.29e+02
	ws->var(             "p0_PTbin0"  )->setVal( -4.8470e-04);// +/-  4.24e-06
	ws->var(             "p0_PTbin1"  )->setVal( -4.7094e-04);// +/-  5.41e-06
	ws->var(             "p0_PTbin2"  )->setVal( -4.6175e-04);// +/-  6.86e-06
	ws->var(             "p0_PTbin3"  )->setVal( -4.5688e-04);// +/-  1.07e-05
	ws->var(             "p0_PTbin4"  )->setVal( -4.1021e-04);// +/-  2.87e-05
//prompt IP
	ws->var(     "promptAsym_PTbin0"  )->setVal( -2.44938e-01);// +/-  1.94264e-02
	ws->var(     "promptAsym_PTbin1"  )->setVal( -3.85341e-01);// +/-  1.12410e-02
	ws->var(     "promptAsym_PTbin2"  )->setVal( -2.32890e-01);// +/-  1.20989e-02
	ws->var(     "promptAsym_PTbin3"  )->setVal( -1.58127e-01);// +/-  1.57671e-02
	ws->var(     "promptAsym_PTbin4"  )->setVal( -2.20983e-01);// +/-  1.23541e-02
	ws->var(     "promptMean_PTbin0"  )->setVal(  7.41476e-01);// +/-  3.22027e-02
	ws->var(     "promptMean_PTbin1"  )->setVal(  8.69684e-01);// +/-  2.10475e-02
	ws->var(     "promptMean_PTbin2"  )->setVal(  6.84181e-01);// +/-  1.58891e-02
	ws->var(     "promptMean_PTbin3"  )->setVal(  6.49780e-01);// +/-  2.59778e-02
	ws->var(     "promptMean_PTbin4"  )->setVal(  7.91027e-01);// +/-  1.74429e-02
	ws->var(     "promptRhoL_PTbin0"  )->setVal(  1.26607e+00);// +/-  4.62087e-02
	ws->var(     "promptRhoL_PTbin1"  )->setVal(  1.63344e+00);// +/-  5.11881e-02
	ws->var(     "promptRhoL_PTbin2"  )->setVal(  1.23796e+00);// +/-  3.79378e-02
	ws->var(     "promptRhoL_PTbin3"  )->setVal(  1.12382e+00);// +/-  2.89408e-02
	ws->var(     "promptRhoL_PTbin4"  )->setVal(  1.32783e+00);// +/-  3.59866e-02
	ws->var(     "promptRhoR_PTbin0"  )->setVal(  1.90330e+00);// +/-  5.90527e-02
	ws->var(     "promptRhoR_PTbin1"  )->setVal(  1.56755e+00);// +/-  3.65293e-02
	ws->var(     "promptRhoR_PTbin2"  )->setVal(  1.65797e+00);// +/-  3.21594e-02
	ws->var(     "promptRhoR_PTbin3"  )->setVal(  1.87261e+00);// +/-  5.23660e-02
	ws->var(     "promptRhoR_PTbin4"  )->setVal(  1.73098e+00);// +/-  4.19735e-02
	ws->var(    "promptWidth_PTbin0"  )->setVal(  1.12137e+00);// +/-  1.50698e-02
	ws->var(    "promptWidth_PTbin1"  )->setVal(  1.17254e+00);// +/-  8.66453e-03
	ws->var(    "promptWidth_PTbin2"  )->setVal(  1.09664e+00);// +/-  1.11622e-02
	ws->var(    "promptWidth_PTbin3"  )->setVal(  1.07137e+00);// +/-  1.07803e-02
	ws->var(    "promptWidth_PTbin4"  )->setVal(  1.15567e+00);// +/-  9.26588e-03
	//ws->var(    "promptYield_PTbin0"  )->setVal(  2.31630e+03);// +/-  2.69955e+01
	//ws->var(    "promptYield_PTbin1"  )->setVal(  8.68934e+03);// +/-  5.83562e+01
	//ws->var(    "promptYield_PTbin2"  )->setVal(  1.25728e+04);// +/-  7.68802e+01
	//ws->var(    "promptYield_PTbin3"  )->setVal(  7.29075e+03);// +/-  5.99772e+01
	//ws->var(    "promptYield_PTbin4"  )->setVal(  1.87655e+04);// +/-  8.82830e+01
//disp IP
	ws->var(     "displcAsym_PTbin0"  )->setVal(  1.1835e-02);// +/-  3.07e-02
	ws->var(     "displcAsym_PTbin1"  )->setVal( -7.3169e-02);// +/-  2.52e-02
	ws->var(     "displcAsym_PTbin2"  )->setVal( -2.1166e-01);// +/-  2.45e-02
	ws->var(     "displcAsym_PTbin3"  )->setVal( -1.6608e-01);// +/-  3.28e-02
	ws->var(     "displcAsym_PTbin4"  )->setVal( -2.5108e-01);// +/-  2.66e-02
	ws->var(     "displcMean_PTbin0"  )->setVal(  3.6616e+00);// +/-  1.23e-01
	ws->var(     "displcMean_PTbin1"  )->setVal(  4.1797e+00);// +/-  1.02e-01
	ws->var(     "displcMean_PTbin2"  )->setVal(  4.8858e+00);// +/-  9.70e-02
	ws->var(     "displcMean_PTbin3"  )->setVal(  4.7789e+00);// +/-  1.29e-01
	ws->var(     "displcMean_PTbin4"  )->setVal(  4.8852e+00);// +/-  1.01e-01
	ws->var(    "displcWidth_PTbin0"  )->setVal(  2.2864e+00);// +/-  2.94e-02
	ws->var(    "displcWidth_PTbin1"  )->setVal(  2.3217e+00);// +/-  2.31e-02
	ws->var(    "displcWidth_PTbin2"  )->setVal(  2.3108e+00);// +/-  2.44e-02
	ws->var(    "displcWidth_PTbin3"  )->setVal(  2.2657e+00);// +/-  2.97e-02
	ws->var(    "displcWidth_PTbin4"  )->setVal(  2.2072e+00);// +/-  2.47e-02
	//ws->var(    "displcYield_PTbin0"  )->setVal(  1.1203e+04);// +/-  2.04e+02
	//ws->var(    "displcYield_PTbin1"  )->setVal(  1.6159e+04);// +/-  2.28e+02
	//ws->var(    "displcYield_PTbin2"  )->setVal(  1.5315e+04);// +/-  2.20e+02
	//ws->var(    "displcYield_PTbin3"  )->setVal(  8.4982e+03);// +/-  1.57e+02
	//ws->var(    "displcYield_PTbin4"  )->setVal(  1.2530e+04);// +/-  1.97e+02
//bkgrnd IP
	//ws->var(    "bkgrndYield_PTbin0"  )->setVal(  4.8410e+03);// +/-  1.31e+02
	//ws->var(    "bkgrndYield_PTbin1"  )->setVal(  7.6104e+03);// +/-  1.55e+02
	//ws->var(    "bkgrndYield_PTbin2"  )->setVal(  7.6209e+03);// +/-  1.50e+02
	//ws->var(    "bkgrndYield_PTbin3"  )->setVal(  4.0486e+03);// +/-  1.09e+02
	//ws->var(    "bkgrndYield_PTbin4"  )->setVal(  5.3398e+03);// +/-  1.29e+02
}

bool SimDFitter::fit(fitType whichFit, fitSet whichSet, uint ptBinMin, uint ptBinMax, unsigned int nToys) {
	bool status=true;

	RooRealVar* DM = ws->var(""+_whichD+"M");
	RooRealVar* DIP = ws->var(""+_whichD+"LogIPChi2");

	TString name="_"+_whichD;
	TString typeName="";
	TString cutStr="ptCat>=";
	cutStr+=ptBinMin;
	cutStr+=" && ptCat<";
	cutStr+=ptBinMax;
	cutStr+=" && JetPT>=";
	cutStr+=_jptmin;
	cutStr+=" && JetPT<";
	cutStr+=_jptmax;

	TString weightName="";
	if(_useEffs) {
		if(_flavour==4) weightName="weight4";
		if(_flavour==5) weightName="weight5";
	}

	switch(whichSet) {
		case fitSet::fitJustCharm:
			cutStr+=" && typeCat==1";
			name+="_prompt";
			break;
		case fitSet::fitJustBeauty:
			cutStr+=" && typeCat==2";
			name+="_displ";
			break;
		case fitSet::fitJustSignal:
			cutStr+=" && typeCat>0 && typeCat<3";
			name+="_signal";
			break;
		case fitSet::fitJustSideband:
			cutStr+=" && typeCat==0";
			name+="_bkgrnd";
			break;
		case fitSet::fitData:
			cutStr+=" && JetEta>=";
			cutStr+=_ymin;
			cutStr+=" && JetEta<";
			cutStr+=_ymax;
			cutStr+=" && typeCat==3";
			name+="_data";
			break;
		case fitSet::fitDataEnhancedC:
			cutStr+=" && JetEta>=";
			cutStr+=_ymin;
			cutStr+=" && JetEta<";
			cutStr+=_ymax;
			cutStr+=" && typeCat==3";
			cutStr+=" && enhancedCat==1";
			name+="_dataEnhancedC";
			break;
		case fitSet::fitDataEnhancedB:
			cutStr+=" && JetEta>=";
			cutStr+=_ymin;
			cutStr+=" && JetEta<";
			cutStr+=_ymax;
			cutStr+=" && typeCat==3";
			cutStr+=" && enhancedCat==2";
			name+="_dataEnhancedB";
			break;
		case fitSet::fitToys:
			cutStr+=" && JetEta>=";
			cutStr+=_ymin;
			cutStr+=" && JetEta<";
			cutStr+=_ymax;
			cutStr+=" && typeCat==3";
			name+="_toys";
			break;
	}
	name+="_";
	name+=_jptmin;
	name+="-";
	name+=_jptmax;
	name+="_";
	name+=_ymin;
	name+="-";
	name+=_ymax;
	if(_useEffs) {
		name+="_effWeight";
		if(_flavour==4) name+="_fit4";
		if(_flavour==5) name+="_fit5";
	}

	RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
	const RooArgSet* vars = data->get();
	RooDataSet* ds(0);
	if(_useEffs) ds = new RooDataSet("ds", "ds", data, *vars, cutStr, weightName);
	else ds = new RooDataSet("ds", "ds", data, *vars, cutStr);
	//RooDataSet* ds = dynamic_cast<RooDataSet*>(ws->data("data")->reduce(cutStr));
	name+="_";
	name+=_ptBinsD[ptBinMin];
	name+="-";
	name+=_ptBinsD[ptBinMax];
	std::cout << ds->sumEntries() << std::endl;//TODO
	vars->Print();//TODO

	RooAbsPdf* data_pdf(0);

	TString toFloatStr, toFixStr, toZeroStr, toOneStr;
	switch(whichFit) {
		case fitType::fitMass:
			data_pdf = ws->pdf("sim_mass_pdf");
			typeName = "Mass";
			switch(whichSet) {
				case fitSet::fitJustCharm:
				case fitSet::fitJustBeauty:
				case fitSet::fitJustSignal:
					toFloatStr="signalYield*,dMean*,dWidthBase*,dRatio*,dAlpha*,dFrac*";
					toFixStr="p0*,dN*";
					toZeroStr="bkgrndYield*";
					toOneStr="dWidthScale*";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*,p0*";
					toFixStr="dMean*,dWidthBase*,dWidthScale*,dRatio*,dAlpha*,dFrac*,dN*";
					toZeroStr="signalYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitDataEnhancedC:
				case fitSet::fitDataEnhancedB:
				case fitSet::fitToys:
					toFloatStr="signalYield*,bkgrndYield*,dMean*,dWidthScale*,p0*";
					toFixStr="dWidthBase*,dRatio*,dAlpha*,dFrac*,dN*";
					toZeroStr="";
					toOneStr="";
					break;
			}
			break;
		case fitType::fitIP:
			data_pdf = ws->pdf("sim_ip_pdf");
			typeName = "IP";
			switch(whichSet) {
				case fitSet::fitJustCharm:
					toFloatStr="promptYield*,promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*";
					toFixStr="promptRhoR*,displcMeanBase*,displcMeanShift*,displcWidthBase*,displcWidthScale*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="displcYield*,bkgrndYield*,promptMeanShift*";
					toOneStr="promptWidthScale*";
					break;
				case fitSet::fitJustBeauty:
					toFloatStr="displcYield*,displcMeanBase*,displcWidthBase*,displcAsym*";
					toFixStr="promptMeanBase*,promptMeanShift*,promptWidthBase*,promptWidthScale*,promptAsym*,promptRhoL*,promptRhoR*,displcRho*,bkgrndBinIdx*";
					toZeroStr="promptYield*,bkgrndYield*,displcMeanShift*";
					toOneStr="displcWidthScale*";
					break;
				case fitSet::fitJustSignal:
					toFloatStr="promptYield*,displcYield*";
					toFixStr="promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="bkgrndYield*,promptMeanShift*,displcMeanShift*";
					toOneStr="promptWidthScale*,displcWidthScale*";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*";
					toFixStr="promptMeanBase*,promptMeanShift*,promptWidthBase*,promptWidthScale*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcMeanShift*,,displcWidthBase*,displcWidthScale*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="promptYield*,displcYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitToys:
				case fitSet::fitDataEnhancedC:
				case fitSet::fitDataEnhancedB:
					toFloatStr="promptYield*,displcYield*,bkgrndYield*,promptMeanShift*,promptWidthScale*,displcMeanShift*,displcWidthScale*";
					toFixStr="promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="";
					toOneStr="";
					break;
				//case fitSet::fitDataEnhancedC:
				//	toFloatStr="promptYield*,displcYield*,bkgrndYield*,promptMeanBase*,promptWidthBase*";
				//	toFixStr="promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,displcMeanShift*,displcWidthScale*,bkgrndBinIdx*";
				//	toZeroStr="promptMeanShift*";
				//	toOneStr="promptWidthScale*";
				//	break;
				//case fitSet::fitDataEnhancedB:
				//	toFloatStr="promptYield*,displcYield*,bkgrndYield*,displcMeanBase*,displcWidthBase*";
				//	toFixStr="promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcAsym*,displcRho*,promptMeanShift*,promptWidthScale*,bkgrndBinIdx*";
				//	toZeroStr="displcMeanShift*";
				//	toOneStr="displcWidthScale*";
				//	break;
			}
			break;
		case fitType::fit2D:
			data_pdf = ws->pdf("sim_data_pdf");
			typeName = "2D";
			switch(whichSet) {
				case fitSet::fitJustCharm:
					toFloatStr="promptYield*,dMean*,dWidthBase*,promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,p0*,promptRhoR*,displcMeanBase*,displcMeanShift*,displcWidthBase*,displcWidthSclae*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="displcYield*,bkgrndYield*,promptMeanShift*";
					toOneStr="dWidthScale*,promptWidthScale*";
					break;
				case fitSet::fitJustBeauty:
					toFloatStr="displcYield*,dMean*,dWidthBase*,displcMeanBase*,displcWidthBase*,displcAsym*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,p0*,promptMeanBase*,promptMeanShift*,promptWidthBase*,promptWidthScale*,promptAsym*,promptRhoL*,promptRhoR*,displcRho*,bkgrndBinIdx*";
					toZeroStr="promptYield*,bkgrndYield*,displcMeanShift*";
					toOneStr="dWidthScale*,displcWidthScale*";
					break;
				case fitSet::fitJustSignal:
					toFloatStr="promptYield*,displcYield*,dMean*,dWidthBase*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,p0*,promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="bkgrndYield*,promptMeanShift*,displcMeanShift*";
					toOneStr="dWidthScale*,promptWidthScale*,displcWidthScale*";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*,p0*";
					toFixStr="dMean*,dWidthBase*,dWidthScale*,dRatio*,dAlpha*,dFrac*,dN*,promptMeanBase*,promptMeanShift*,promptWidthBase*,promptWidthScale*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcMeanShift*,displcWidthBase*,displcWidthScale*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="promptYield*,displcYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitToys:
				case fitSet::fitDataEnhancedC:
				case fitSet::fitDataEnhancedB:
					toFloatStr="promptYield*,displcYield*,bkgrndYield*,dMean*,dWidthScale*,p0*,promptMeanShift*,promptWidthScale*,displcMeanShift*,displcWidthScale*";
					toFixStr="dWidthBase*,dRatio*,dAlpha*,dFrac*,dN*,promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,bkgrndBinIdx*";
					toZeroStr="";
					toOneStr="";
					break;
				//case fitSet::fitDataEnhancedC:
				//	toFloatStr="promptYield*,displcYield*,bkgrndYield*,dMean*,dWidthScale*,p0*,promptMeanBase*,promptWidthBase*";
				//	toFixStr="dWidthBase*,dRatio*,dAlpha*,dFrac*,dN*,promptAsym*,promptRhoL*,promptRhoR*,displcMeanBase*,displcWidthBase*,displcAsym*,displcRho*,displcMeanShift*,displcWidthScale*,bkgrndBinIdx*";
				//	toZeroStr="promptMeanShift*";
				//	toOneStr="promptWidthScale*";
				//	break;
				//case fitSet::fitDataEnhancedB:
				//	toFloatStr="promptYield*,displcYield*,bkgrndYield*,dMean*,dWidthScale*,p0*,displcMeanBase*,displcWidthBase*";
				//	toFixStr="dWidthBase*,dRatio*,dAlpha*,dFrac*,dN*,promptMeanBase*,promptWidthBase*,promptAsym*,promptRhoL*,promptRhoR*,displcAsym*,displcRho*,promptMeanShift*,promptWidthScale*,bkgrndBinIdx*";
				//	toZeroStr="displcMeanShift*";
				//	toOneStr="displcWidthScale*";
				//	break;
			}
			break;
	}

	//unused PT bins
	for(uint ibin=0; ibin<_ptBinsD.size()-1; ++ibin) {
		if(ibin<ptBinMin || ibin>=ptBinMax) {
			toFixStr+=",*PTbin";
			toFixStr+=ibin;
			toZeroStr+=",*Yield*PTbin";
			toZeroStr+=ibin;
		}
	}

	//std::cout << toFloatStr << std::endl;//TODO
	//std::cout << toFixStr << std::endl;//TODO
	//std::cout << toZeroStr << std::endl;//TODO

	RooArgSet*  allPars   = data_pdf->getParameters(*ds);
	RooArgList* floatPars = static_cast<RooArgList*>(allPars->selectByName(toFloatStr));
	RooArgList* fixPars   = static_cast<RooArgList*>(allPars->selectByName(toFixStr));
	RooArgList* zeroPars  = static_cast<RooArgList*>(allPars->selectByName(toZeroStr));
	RooArgList* onePars   = static_cast<RooArgList*>(allPars->selectByName(toOneStr));
	RooArgList* yield4Pars= static_cast<RooArgList*>(allPars->selectByName("*promptYield*"));
	RooArgList* yield5Pars= static_cast<RooArgList*>(allPars->selectByName("*displcYield*"));

	for ( int i=0; i<floatPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(floatPars->at(i))->setConstant(false);
	}
	for ( int i=0; i<fixPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(fixPars->at(i))->setConstant(true);
	}
	for ( int i=0; i<onePars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(onePars->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(onePars->at(i))->setVal(1.0);
	}
	for ( int i=0; i<zeroPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(zeroPars->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(zeroPars->at(i))->setVal(0.0);
	}

	//check for unused or overlapping parameters
	RooArgList missingPars(*allPars);
	missingPars.remove(*floatPars);
	missingPars.remove(*fixPars);
	missingPars.remove(*zeroPars);
	missingPars.remove(*onePars);
	if(!missingPars.empty()) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were not specified" << std::endl;
		missingPars.printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*fixPars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*fixPars))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*zeroPars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed to zero" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*zeroPars))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*onePars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed to one" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*onePars))->printMultiline(std::cout,0xFF);
	}
	if(zeroPars->overlaps(*onePars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be fixed to both zero and one" << std::endl;
		static_cast<RooArgList*>(zeroPars->selectCommon(*onePars))->printMultiline(std::cout,0xFF);
	}

	//set sensible starting values for yields
	for(uint ibin=ptBinMin; ibin<ptBinMax; ++ibin) {
		//extract the yield in each bin from the dataset
		TString ptCutStr="ptCat==";
		ptCutStr+=ibin;
		double binYield = ds->sumEntries(ptCutStr);
		std::cout << "entries" << ibin << " " << binYield << std::endl;//TODO

		if(binYield>=1) {
			TString yieldNameStr="*Yield*PTbin";
			yieldNameStr+=ibin;
			RooArgList* yields= static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(yieldNameStr)->selectByAttrib("Constant",kFALSE));
			for ( int i=0; i<yields->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(yields->at(i))->setVal(binYield/yields->getSize());
				//std::cout << dynamic_cast<RooRealVar*>(yields->at(i))->getTitle() << " set to " << dynamic_cast<RooRealVar*>(yields->at(i))->getValV() << std::endl;
			}
		} else {
			std::cout << "Bin " << ibin << " has yield of " << binYield << ". Fixing to 0..." << std::endl;
			toFixStr+=",*PTbin";
			toFixStr+=ibin;
			toZeroStr+=",*Yield*PTbin";
			toZeroStr+=ibin;
		}
		if(binYield<200) {//TODO
		//if(binYield<20) {//TODO
			std::cout << "Bin " << ibin << " has yield of " << binYield << ". Fixing combinatorial gradient..." << std::endl;
			toZeroStr+=",p0_PTbin";
			toZeroStr+=ibin;
		}
		if(binYield<20) {
			std::cout << "Bin " << ibin << " has yield of " << binYield << ". Fixing data/MC corrections..." << std::endl;
			toOneStr+=",*Scale_PTbin";
			toOneStr+=ibin;
			toZeroStr+=",*Shift_PTbin";
			toZeroStr+=ibin;
			//in case of enhanced data fits, also fix the base means and widths for low stats
			if(whichSet==fitSet::fitDataEnhancedC || whichSet==fitSet::fitDataEnhancedB) {
				toFixStr+=",*Base_PTbin";
			}
		}
	}

	RooArgList* fixPars2   = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(toFixStr));
	RooArgList* zeroPars2  = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(toZeroStr));
	RooArgList* onePars2   = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(toOneStr));

	if(floatPars->overlaps(*fixPars2)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to fixed due to low bin populations" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*fixPars2))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*zeroPars2)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to zeroed due to low bin populations" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*zeroPars2))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*onePars2)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to fixed to one due to low bin populations" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*onePars2))->printMultiline(std::cout,0xFF);
	}

	//repeat in case more have been added
	for ( int i=0; i<fixPars2->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(fixPars2->at(i))->setConstant(true);
	}
	for ( int i=0; i<zeroPars2->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(zeroPars2->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(zeroPars2->at(i))->setVal(0.0);
	}
	for ( int i=0; i<onePars2->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(onePars2->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(onePars2->at(i))->setVal(1.0);
	}

	//set values of any shift or scale parameters that were given in the fit options
	if(whichSet==fitSet::fitData || whichSet==fitSet::fitDataEnhancedC || whichSet==fitSet::fitDataEnhancedB || whichSet==fitSet::fitToys) {
		if(_setMassWidthScale!=-999.) {
			RooArgList* setPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("dWidthScale*"));
			std::cout << "INFO in SimDFitter::fit: setting "+_whichD+" mass width scale factor to " << _setMassWidthScale << std::endl;
			for ( int i=0; i<setPars->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(setPars->at(i))->setVal(_setMassWidthScale);
			}
		}
		if(_setPromptMeanShift!=-999.) {
			RooArgList* setPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("promptMeanShift*"));
			std::cout << "INFO in SimDFitter::fit: setting prompt IP mean shift factor to " << _setPromptMeanShift << std::endl;
			for ( int i=0; i<setPars->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(setPars->at(i))->setVal(_setPromptMeanShift);
			}
		}
		if(_setPromptWidthScale!=-999.) {
			RooArgList* setPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("promptWidthScale*"));
			std::cout << "INFO in SimDFitter::fit: setting prompt IP width scale factor to " << _setPromptWidthScale << std::endl;
			for ( int i=0; i<setPars->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(setPars->at(i))->setVal(_setPromptWidthScale);
			}
		}
		if(_setDisplcMeanShift!=-999.) {
			RooArgList* setPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("displcMeanShift*"));
			std::cout << "INFO in SimDFitter::fit: setting displaced IP mean shift factor to " << _setDisplcMeanShift << std::endl;
			for ( int i=0; i<setPars->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(setPars->at(i))->setVal(_setDisplcMeanShift);
			}
		}
		if(_setDisplcWidthScale!=-999.) {
			RooArgList* setPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("displcWidthScale*"));
			std::cout << "INFO in SimDFitter::fit: setting displaced IP width scale factor to " << _setDisplcWidthScale << std::endl;
			for ( int i=0; i<setPars->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(setPars->at(i))->setVal(_setDisplcWidthScale);
			}
		}
	}

	//finally, remove any shift or scale parameters that were turned off in the fit options
	TString extraFixedZeroStr(""), extraFixedOneStr(""), extraFixedStr("");
	if(whichSet==fitSet::fitDataEnhancedC) {
		if(!_enhanceMassWidthScale)   extraFixedOneStr+=",dWidthScale*";
		if(!_enhancePromptMeanShift)  extraFixedZeroStr+=",promptMeanShift*";
		if(!_enhancePromptWidthScale) extraFixedOneStr+=",promptWidthScale*";
		extraFixedStr+=",displcMeanShift*";
		extraFixedStr+=",displcWidthScale*";
	} else if(whichSet==fitSet::fitDataEnhancedB) {
		if(!_enhanceMassWidthScale)   extraFixedOneStr+=",dWidthScale*";
		if(!_enhanceDisplcMeanShift)  extraFixedZeroStr+=",displcMeanShift*";
		if(!_enhanceDisplcWidthScale) extraFixedOneStr+=",displcWidthScale*";
		extraFixedStr+=",promptMeanShift*";
		extraFixedStr+=",promptWidthScale*";
	} else {
		if(!_allowMassWidthScale)   extraFixedOneStr+=",dWidthScale*";
		if(!_allowPromptMeanShift)  extraFixedZeroStr+=",promptMeanShift*";
		if(!_allowPromptWidthScale) extraFixedOneStr+=",promptWidthScale*";
		if(!_allowDisplcMeanShift)  extraFixedZeroStr+=",displcMeanShift*";
		if(!_allowDisplcWidthScale) extraFixedOneStr+=",displcWidthScale*";
	}
	RooArgList* fixedParsExtra  = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(extraFixedStr));
	if(floatPars->overlaps(*fixedParsExtra)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to fixed due to configured fit options" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*fixedParsExtra))->printMultiline(std::cout,0xFF);
	}
	for ( int i=0; i<fixedParsExtra->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(fixedParsExtra->at(i))->setConstant(true);
	}
	RooArgList* zeroParsExtra  = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(extraFixedZeroStr));
	if(floatPars->overlaps(*zeroParsExtra)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to zeroed due to configured fit options" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*zeroParsExtra))->printMultiline(std::cout,0xFF);
	}
	for ( int i=0; i<zeroParsExtra->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(zeroParsExtra->at(i))->setConstant(true);
		//dynamic_cast<RooRealVar*>(zeroParsExtra->at(i))->setVal(0.0);
	}
	RooArgList* oneParsExtra  = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(extraFixedOneStr));
	if(floatPars->overlaps(*oneParsExtra)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to fixed to one due to configured fit options" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*oneParsExtra))->printMultiline(std::cout,0xFF);
	}
	for ( int i=0; i<oneParsExtra->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(oneParsExtra->at(i))->setConstant(true);
		//dynamic_cast<RooRealVar*>(oneParsExtra->at(i))->setVal(1.0);
	}

	//Print final list of floated parameters before fitting
	//std::cout << "Fit parameters:" << std::endl;
	//static_cast<RooArgList*>(allPars->selectByAttrib("Constant",kFALSE))->printMultiline(std::cout, 0xFF);
	
	double minCovQual(3);
	minCovQual=2;//allow not pos def//TODO
	gSystem->RedirectOutput(gSaveDir+"/log/fitFull"+typeName+name+".log","w");
	RooAbsReal* nll = data_pdf->createNLL(*ds, RooFit::Extended(), RooFit::NumCPU(4), RooFit::Range("FIT"));
	RooMinuit m(*nll);
	m.setPrintEvalErrors(-1);
	std::cout << "FIT STAGE 1" << std::endl;//TODO
	m.migrad();
	RooFitResult* r1 = m.save();
	std::cout << "STATUS=" << r1->status() << ", COVQUAL=" << r1->covQual() << std::endl;//TODO

	std::cout << "FIT STAGE 2" << std::endl;//TODO
	m.hesse();
	RooFitResult* r2 = m.save();
	std::cout << "STATUS=" << r2->status() << ", COVQUAL=" << r2->covQual() << std::endl;//TODO
	RooFitResult *r3(0), *r4(0), *r5(0), *r6(0);
	if(r2->status()!=0 || r2->covQual()<minCovQual) {
		std::cout << "Hesse failed STATUS=" << r2->status() << ", COVQUAL=" << r2->covQual() << std::endl;
		std::cout << "Retrying fit" << std::endl;
		std::cout << "FIT STAGE 3" << std::endl;//TODO
		m.migrad();
		r3 = m.save();
		std::cout << "STATUS=" << r3->status() << ", COVQUAL=" << r3->covQual() << std::endl;//TODO
		std::cout << "FIT STAGE 4" << std::endl;//TODO
		m.hesse();
		r4 = m.save();
		std::cout << "STATUS=" << r4->status() << ", COVQUAL=" << r4->covQual() << std::endl;//TODO
		if(r4->status()!=0 || r4->covQual()<minCovQual) {
			std::cout << "Hesse failed STATUS=" << r4->status() << ", COVQUAL=" << r4->covQual() << std::endl;
			//if(r4->status()!=0 || r4->covQual()<2) {//TODO allow covQual==2 for now
			status=false;
			//}
		}
	}
	//if we're using efficiency weights then re-run the fit using weights to calculate uncertainties
	if(_useEffs && !_skipSumW2Fits) {
		std::cout << "FIT STAGE 5" << std::endl;//TODO
		RooFitResult* r5 = data_pdf->fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(true), RooFit::PrintLevel(3));//TODO
		std::cout << "STATUS=" << r5->status() << ", COVQUAL=" << r5->covQual() << std::endl;//TODO
		status = (r5->status()==0);
		if(r5->status()!=0) {
			std::cout << "FIT STAGE 6" << std::endl;//TODO
			r6 = data_pdf->fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(true));
			std::cout << "STATUS=" << r6->status() << ", COVQUAL=" << r6->covQual() << std::endl;//TODO
			status = (r6->status()==0);
		}
	}
	gSystem->RedirectOutput(0);
	//static_cast<RooArgList*>(allPars->selectByAttrib("Constant",kFALSE))->printLatex(std::cout, 1, "NEYU", 2);
	delete nll;
	std::ofstream logfile(gSaveDir+"/log/fit"+typeName+name+".log",std::ofstream::out);
	r1->printMultiline(logfile,0);
	r2->printMultiline(logfile,0);
	if(r3) r3->printMultiline(logfile,0);
	if(r4) r4->printMultiline(logfile,0);
	if(r5) r5->printMultiline(logfile,0);
	if(r6) r6->printMultiline(logfile,0);
	logfile.close();

	if(true) { //status) {
		if(whichSet==fitSet::fitToys) {
			//do toy fits
			const RooCategory* cat = dynamic_cast<const RooCategory*>(&dynamic_cast<RooSimultaneous*>(data_pdf)->indexCat());
			RooMCStudy mcs(*data_pdf, RooArgSet(*DM,*DIP,*cat), RooFit::Extended(), RooFit::FitOptions(RooFit::PrintEvalErrors(-1), RooFit::PrintLevel(2), RooFit::Warnings(kFALSE)),RooFit::Silence());
			mcs.generateAndFit(nToys, ds->sumEntries(),true,gSaveDir+"/dtoys/toy_"+name+"_%04d.dat");
			//TODO//mcs.fit(nToys, gSaveDir+"/dtoys/toy_"+name+"_%04d.dat");
			mcs.fitParDataSet().write(gSaveDir+"/dtoys/allToys_"+name+".dat");
			std::cout << "Toy parameters:" << std::endl;
			mcs.fitParDataSet().printMultiline(std::cout,0);
			TCanvas c1;
			TString varName;
			RooRealVar* var(0);
			for(uint ibin=ptBinMin; ibin<ptBinMax; ++ibin) {
				varName = "promptYield_PTbin";
				varName+=ibin;
				var = ws->var(varName);
				RooPlot* promptYieldPlot = mcs.plotPull(*var,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-3,3));
				promptYieldPlot->SetTitle("");
				promptYieldPlot->SetXTitle("Pull(N_{prompt})");
				promptYieldPlot->Draw();
				c1.SaveAs(gSaveDir+"/fig/pull_"+varName+"_"+name+".pdf");

				varName = "displcYield_PTbin";
				varName+=ibin;
				var = ws->var(varName);
				RooPlot* displYieldPlot  = mcs.plotPull(*var,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-3,3));
				displYieldPlot->SetTitle("");
				displYieldPlot->SetXTitle("Pull(N_{displ.})");
				displYieldPlot->Draw();
				c1.SaveAs(gSaveDir+"/fig/pull_"+varName+"_"+name+".pdf");
			}
		}
		//gSystem->RedirectOutput(0);
		if(whichSet!=fitSet::fitToys) {
//			double yieldS = signalYield.getVal();
//			double yieldB = bkgrndYield.getVal();
//			std::cout << yieldS << "\t" << yieldB << std::endl;
			
			std::vector<std::string> sig_pdfs;
			if(whichFit==fitType::fitMass && whichSet!=fitSet::fitJustSideband) {
				sig_pdfs.push_back( "signal*" );
			}
			if(whichFit!=fitType::fitMass && whichSet!=fitSet::fitJustSideband && whichSet!=fitSet::fitJustBeauty) {
				sig_pdfs.push_back( "prompt*" );
			}
			if(whichFit!=fitType::fitMass && whichSet!=fitSet::fitJustSideband && whichSet!=fitSet::fitJustCharm) {
				sig_pdfs.push_back( "displc*" );
			}
			std::vector<std::string> bkg_pdfs;
			if(whichSet!=fitSet::fitJustCharm && whichSet!=fitSet::fitJustBeauty && whichSet!=fitSet::fitJustSignal) {
				bkg_pdfs.push_back( "bkgrnd*" );
			}

			//for single-bin fits just plot that bin
			//int plotBins(-1);
			//if(ptBinMax==ptBinMin+1) {
			//	plotBins=ptBinMin;
			//}

			switch(whichFit) {
				case fitType::fitMass:
					//plotFit(*DM,     1784., 1944., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fitMass"+name, "#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]", plotBins);
					makeMassPlots(ds,name,ptBinMin,ptBinMax);
					break;
				case fitType::fitIP:
					//plotFit(*DIP,      -5.,   15., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fitIP"+name, "log #chi_{IP}^{2}", plotBins);
					makeIPPlots(ds,name,ptBinMin,ptBinMax);
					break;
				case fitType::fit2D:
					//plotFit(*DM,     1784., 1944., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fit2D_Mass"+name, "#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]", plotBins);
					//plotFit(*DIP,      -5.,   15., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fit2D_IP"+name, "log #chi_{IP}^{2}", plotBins);
					make2DPlots(ds,name,ptBinMin,ptBinMax);
					break;
			}
		}

		if(whichFit==fitType::fit2D && whichSet!=fitSet::fitToys) {
			_yield4=0.;
			_error4=0.;
			_yield5=0.;
			_error5=0.;
			_names4.clear();
			_names5.clear();
			_yields4.clear();
			_yields5.clear();
			_errors4.clear();
			_errors5.clear();
			for ( int i=0; i<yield4Pars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(yield4Pars->at(i));
				//std::cout << v->getTitle() << " " << v->getVal() << " " << v->getError() << std::endl;//TODO
				_yield4+= v->getVal();
				_error4+= v->getError()*v->getError();
				_names4.push_back(v->getTitle());
				_yields4.push_back(v->getVal());
				_errors4.push_back(v->getError());
			}
			_error4 = TMath::Sqrt(_error4);
			//std::cout << "..." << std::endl;//TODO
			for ( int i=0; i<yield5Pars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(yield5Pars->at(i));
				//std::cout << v->getTitle() << " " << v->getVal() << " " << v->getError() << std::endl;//TODO
				_yield5+= v->getVal();
				_error5+= v->getError()*v->getError();
				_names5.push_back(v->getTitle());
				_yields5.push_back(v->getVal());
				_errors5.push_back(v->getError());
			}
			_error5 = TMath::Sqrt(_error5);
		}

		printParams("dat/params"+typeName+name+".dat",ds,data_pdf);
		printAllParams("dat/paramsAll"+typeName+name+".dat",ds,data_pdf);

		//if we're applying a pull to any of our fixed parameters then do it here after printing
		if(whichFit == fitType::fitMass && whichSet==fitSet::fitJustSignal
		   && TMath::Abs(_shiftFixedMassWidth)>1e-2) {
			std::cout << "Shifting mass width by parameters by " << _shiftFixedMassWidth << " sigma:" << std::endl;
			RooArgList* shiftPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("dWidthBase*"));
			for ( int i=0; i<shiftPars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(shiftPars->at(i));
				std::cout << v->getVal() << " +/- " << v->getError() << " ->";
				v->setVal(v->getVal() + _shiftFixedMassWidth*v->getError());
				std::cout << v->getVal() << std::endl;
			}
		}
		if(whichFit == fitType::fitIP && whichSet==fitSet::fitJustCharm
		   && TMath::Abs(_shiftFixedPromptMean)>1e-2) {
			std::cout << "Shifting prompt mean by parameters by " << _shiftFixedPromptMean << " sigma:" << std::endl;
			RooArgList* shiftPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("promptMeanBase*"));
			for ( int i=0; i<shiftPars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(shiftPars->at(i));
				std::cout << v->getVal() << " +/- " << v->getError() << " ->";
				v->setVal(v->getVal() + _shiftFixedPromptMean*v->getError());
				std::cout << v->getVal() << std::endl;
			}
		}
		if(whichFit == fitType::fitIP && whichSet==fitSet::fitJustCharm
		   && TMath::Abs(_shiftFixedPromptWidth)>1e-2) {
			std::cout << "Shifting prompt width by parameters by " << _shiftFixedPromptWidth << " sigma:" << std::endl;
			RooArgList* shiftPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("promptWidthBase*"));
			for ( int i=0; i<shiftPars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(shiftPars->at(i));
				std::cout << v->getVal() << " +/- " << v->getError() << " ->";
				v->setVal(v->getVal() + _shiftFixedPromptWidth*v->getError());
				std::cout << v->getVal() << std::endl;
			}
		}
		if(whichFit == fitType::fitIP && whichSet==fitSet::fitJustBeauty
		   && TMath::Abs(_shiftFixedDisplcMean)>1e-2) {
			std::cout << "Shifting displaced mean by parameters by " << _shiftFixedDisplcMean << " sigma:" << std::endl;
			RooArgList* shiftPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("displcMeanBase*"));
			for ( int i=0; i<shiftPars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(shiftPars->at(i));
				std::cout << v->getVal() << " +/- " << v->getError() << " ->";
				v->setVal(v->getVal() + _shiftFixedDisplcMean*v->getError());
				std::cout << v->getVal() << std::endl;
			}
		}
		if(whichFit == fitType::fitIP && whichSet==fitSet::fitJustBeauty
		   && TMath::Abs(_shiftFixedDisplcWidth)>1e-2) {
			std::cout << "Shifting displaced width by parameters by " << _shiftFixedDisplcWidth << " sigma:" << std::endl;
			RooArgList* shiftPars = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName("displcWidthBase*"));
			for ( int i=0; i<shiftPars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(shiftPars->at(i));
				std::cout << v->getVal() << " +/- " << v->getError() << " ->";
				v->setVal(v->getVal() + _shiftFixedDisplcWidth*v->getError());
				std::cout << v->getVal() << std::endl;
			}
		}
	}

	delete ds;

	return status;
}

void SimDFitter::makeMassPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	name.ReplaceAll(".","");
	RooRealVar* DM = ws->var(""+_whichD+"M");

	uint nBins1D(80);
	uint nBins(160);

	double mmin(_mLow), mmax(_mHigh);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);

	//pT-integrated histograms
	TH1* dataHist   = static_cast<TH1*>(ds->createHistogram("dataHist",*DM,RooFit::Binning(nBins,mmin,mmax)));
	TH1* modelHist  = new TH1D("modelHist", "",nBins,mmin,mmax);
	TH1* signalHist = new TH1D("signalHist","",nBins,mmin,mmax);
	TH1* bkgrndHist = new TH1D("bkgrndHist","",nBins,mmin,mmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}("+_decayLabel+") [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	bkgrndHist->SetLineColor(kBlue);

	//store integrated and pT slice histograms
	std::vector<TH1*> dataHists {dataHist};
	std::vector<TH1*> modelHists {modelHist};
	std::vector<TH1*> signalHists {signalHist};
	std::vector<TH1*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* signal = ws->pdf(TString("signalMass_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrndMass_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!signal) signal = ws->pdf(TString("signalMass"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrndMass"));

		TH1* signalSlice = static_cast<TH1*>(signal->createHistogram(TString::Format("signalSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax)));
		TH1* bkgrndSlice = static_cast<TH1*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax)));
		TH1* dataSlice   = static_cast<TH1*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH1* modelSlice  = new TH1D(TString::Format("modelSlice%d",ibin), "",nBins,mmin,mmax);

		signalSlice->Scale(ws->var(TString::Format("signalYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(signalSlice);
		modelSlice->Add(bkgrndSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("#it{m}("+_decayLabel+") [MeV/#it{c}^{2}]");
			hist->GetYaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH1* dataHist50 = dataHists[i]->Rebin( nBins/nBins1D, "dataHist50");

		//Mass plots
		TH1* dataHistMass   = dataHist50;
		TH1* signalHistMass = signalHists[i];
		TH1* bkgrndHistMass = bkgrndHists[i];
		TH1* modelHistMass  = modelHists[i];

		signalHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());

		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (mmax-mmin)/nBins1D));
		bkgrndHistMass->SetLineStyle(kDashed);
		signalHistMass->SetFillColor(kRed);
		signalHistMass->SetFillStyle(3350);

		dataHistMass->Draw("E1 X0 P");
		if(signalHistMass->Integral()>0) signalHistMass->Draw("HIST C SAME");
		if(bkgrndHistMass->Integral()>0) bkgrndHistMass->Draw("HIST C SAME");
		modelHistMass->Draw( "HIST C SAME");
		dataHistMass->Draw("E1 X0 P SAME");

		TPaveText label(0.65,0.8,0.9,0.9,"BRNDC");
		if(i==0) {
			label.AddText(TString::Format(
				"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
				_jptmin/1e3,
				_jptmax/1e3));
		} else {
			if(_usePtFracBins) {
				label.AddText(TString::Format(
					"#it{f}_{#it{p}_{T}}^{#it{D}} #in [%.2f,%.2f]",
					static_cast<double>(i-1)/_nDPtBins,
					static_cast<double>(i)/_nDPtBins));
			} else {
				label.AddText(TString::Format(
					"#it{p}_{T}(#it{D}) #in [%.0f,%.0f] GeV/#it{c}",
					_ptBinsD[i-1],
					_ptBinsD[i]));
			}
		}
		label.SetFillColor(0);
		label.SetTextAlign(12);
		label.SetBorderSize(0);
		label.Draw();
		if(i==0) can.SaveAs(gSaveDir+"/fig/fitMass"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fitMass"+name+"PTbin"+(i-1)+".pdf");
	}
}

void SimDFitter::makeIPPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	name.ReplaceAll(".","");
	RooRealVar* DIP = ws->var(""+_whichD+"LogIPChi2");

	uint nBins1D(80);
	uint nBins(160);

	double ipmin(-5.), ipmax(15.);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);

	//pT-integrated histograms
	TH1* dataHist   = static_cast<TH1*>(ds->createHistogram("dataHist",*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
	TH1* modelHist  = new TH1D("modelHist", "",nBins,ipmin,ipmax);
	TH1* signalHist = new TH1D("signalHist","",nBins,ipmin,ipmax);
	TH1* promptHist = new TH1D("promptHist","",nBins,ipmin,ipmax);
	TH1* displcHist = new TH1D("displcHist","",nBins,ipmin,ipmax);
	TH1* bkgrndHist = new TH1D("bkgrndHist","",nBins,ipmin,ipmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,promptHist,displcHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("log #chi_{IP}^{2}");
		hist->GetYaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	promptHist->SetLineColor(kGreen+2);
	displcHist->SetLineColor(kMagenta);
	bkgrndHist->SetLineColor(kBlue);

	//store integrated and pT slice histograms
	std::vector<TH1*> dataHists {dataHist};
	std::vector<TH1*> modelHists {modelHist};
	std::vector<TH1*> signalHists {signalHist};
	std::vector<TH1*> promptHists {promptHist};
	std::vector<TH1*> displcHists {displcHist};
	std::vector<TH1*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* prompt = ws->pdf(TString("promptIP_PTbin")+ibin);
		RooAbsPdf* displc = ws->pdf(TString("displcIP_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrndIP_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!prompt) prompt = ws->pdf(TString("promptIP"));
		if(!displc) displc = ws->pdf(TString("displcIP"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrndIP"));

		TH1* promptSlice = static_cast<TH2*>(prompt->createHistogram(TString::Format("promptSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* displcSlice = static_cast<TH2*>(displc->createHistogram(TString::Format("displcSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* bkgrndSlice = static_cast<TH2*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* dataSlice   = static_cast<TH2*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DIP,RooFit::Binning(nBins,ipmin,ipmax),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH1* modelSlice  = new TH1D(TString::Format("modelSlice%d",ibin), "",nBins,ipmin,ipmax);
		TH1* signalSlice = new TH1D(TString::Format("signalSlice%d",ibin),"",nBins,ipmin,ipmax);

		promptSlice->Scale(ws->var(TString::Format("promptYield_PTbin%d",ibin))->getVal());
		displcSlice->Scale(ws->var(TString::Format("displcYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(promptSlice);
		modelSlice->Add(displcSlice);
		modelSlice->Add(bkgrndSlice);
		signalSlice->Add(promptSlice);
		signalSlice->Add(displcSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		promptHist->Add(promptSlice);
		displcHist->Add(displcSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,promptSlice,displcSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("log #chi_{IP}^{2}");
			hist->GetYaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		promptSlice->SetLineColor(kGreen+2);
		displcSlice->SetLineColor(kMagenta);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		promptHists.push_back(promptSlice);
		displcHists.push_back(displcSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH1* dataHist50 = dataHists[i]->Rebin( nBins/nBins1D,"dataHist50");

		//IP plots
		TH1* dataHistIP   = dataHist50;
		TH1* signalHistIP = signalHists[i];
		TH1* promptHistIP = promptHists[i];
		TH1* displcHistIP = displcHists[i];
		TH1* bkgrndHistIP = bkgrndHists[i];
		TH1* modelHistIP  = modelHists[i];

		signalHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		promptHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		displcHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		bkgrndHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		modelHistIP ->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());

		dataHistIP->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f)", (ipmax-ipmin)/nBins1D));
		bkgrndHistIP->SetLineStyle(kDashed);
		signalHistIP->SetLineStyle(kDotted);
		promptHistIP->SetFillColor(kGreen+2);
		displcHistIP->SetFillColor(kMagenta);
		promptHistIP->SetFillStyle(3345);
		displcHistIP->SetFillStyle(3354);

		dataHistIP->Draw("E1 X0 P");
		if(promptHistIP->Integral()>0 && displcHistIP->Integral()>0) signalHistIP->Draw("HIST C SAME");
		if(promptHistIP->Integral()>0) promptHistIP->Draw("HIST C SAME");
		if(displcHistIP->Integral()>0) displcHistIP->Draw("HIST C SAME");
		if(bkgrndHistIP->Integral()>0) bkgrndHistIP->Draw("HIST C SAME");
		modelHistIP->Draw( "HIST C SAME");
		dataHistIP->Draw("E1 X0 P SAME");

		TPaveText label(0.65,0.8,0.9,0.9,"BRNDC");
		if(i==0) {
			label.AddText(TString::Format(
				"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
				_jptmin/1e3,
				_jptmax/1e3));
		} else {
			if(_usePtFracBins) {
				label.AddText(TString::Format(
					"#it{f}_{#it{p}_{T}}^{#it{D}} #in [%.2f,%.2f]",
					static_cast<double>(i-1)/_nDPtBins,
					static_cast<double>(i)/_nDPtBins));
			} else {
				label.AddText(TString::Format(
					"#it{p}_{T}(#it{D}) #in [%.0f,%.0f] GeV/#it{c}",
					_ptBinsD[i-1],
					_ptBinsD[i]));
			}
		}
		label.SetFillColor(0);
		label.SetTextAlign(12);
		label.SetBorderSize(0);
		label.Draw();
		if(i==0) can.SaveAs(gSaveDir+"/fig/fitIP"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fitIP"+name+"PTbin"+(i-1)+".pdf");
	}
}

void SimDFitter::make2DPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	name.ReplaceAll(".","");
	RooRealVar* DM = ws->var(""+_whichD+"M");
	RooRealVar* DIP = ws->var(""+_whichD+"LogIPChi2");

	uint nBins1D(50);
	uint nBins2D(30);
	uint nBins(600);

	double mmin(_mLow), mmax(_mHigh), ipmin(-5.), ipmax(15.);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
	Double_t greens[NRGBs] = { 1.00, 0.95, 0.50, 0.00, 0.00};
	Double_t blues[NRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);

	//pT-integrated histograms
	TH2* dataHist   = static_cast<TH2*>(ds->createHistogram("dataHist",*DM,RooFit::Binning(nBins,mmin,mmax),RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
	TH2* modelHist  = new TH2D("modelHist", "",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* signalHist = new TH2D("signalHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* promptHist = new TH2D("promptHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* displcHist = new TH2D("displcHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* bkgrndHist = new TH2D("bkgrndHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,promptHist,displcHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}("+_decayLabel+") [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
		hist->GetZaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	promptHist->SetLineColor(kGreen+2);
	displcHist->SetLineColor(kMagenta);
	bkgrndHist->SetLineColor(kBlue);
	promptHist->SetContour(5);
	displcHist->SetContour(5);
	bkgrndHist->SetContour(5);

	//store integrated and pT slice histograms
	std::vector<TH2*> dataHists {dataHist};
	std::vector<TH2*> modelHists {modelHist};
	std::vector<TH2*> signalHists {signalHist};
	std::vector<TH2*> promptHists {promptHist};
	std::vector<TH2*> displcHists {displcHist};
	std::vector<TH2*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* prompt = ws->pdf(TString("prompt_PTbin")+ibin);
		RooAbsPdf* displc = ws->pdf(TString("displc_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrnd_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!prompt) prompt = ws->pdf(TString("prompt"));
		if(!displc) displc = ws->pdf(TString("displc"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrnd"));

		TH2* promptSlice = static_cast<TH2*>(prompt->createHistogram(TString::Format("promptSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* displcSlice = static_cast<TH2*>(displc->createHistogram(TString::Format("displcSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* bkgrndSlice = static_cast<TH2*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* dataSlice   = static_cast<TH2*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax)),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH2* modelSlice  = new TH2D(TString::Format("modelSlice%d",ibin), "",nBins,mmin,mmax,nBins,ipmin,ipmax);
		TH2* signalSlice = new TH2D(TString::Format("signalSlice%d",ibin),"",nBins,mmin,mmax,nBins,ipmin,ipmax);

		promptSlice->Scale(ws->var(TString::Format("promptYield_PTbin%d",ibin))->getVal());
		displcSlice->Scale(ws->var(TString::Format("displcYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(promptSlice);
		modelSlice->Add(displcSlice);
		modelSlice->Add(bkgrndSlice);
		signalSlice->Add(promptSlice);
		signalSlice->Add(displcSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		promptHist->Add(promptSlice);
		displcHist->Add(displcSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,promptSlice,displcSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("#it{m}("+_decayLabel+") [MeV/#it{c}^{2}]");
			hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
			hist->GetZaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		promptSlice->SetLineColor(kGreen+2);
		displcSlice->SetLineColor(kMagenta);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		promptHists.push_back(promptSlice);
		displcHists.push_back(displcSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	//make pull histograms
	TH2* dataHist10 = dataHist->Rebin2D(nBins/10,nBins/10,"dataHist10");
	TH2* modelHist10 = modelHist->Rebin2D(nBins/10,nBins/10,"modelHist10");

	//make pull histograms with adaptive binning
	int nBinAdapt = TMath::Nint(dataHist10->GetSumOfWeights()/30.); //target 30 entries per bins
	if(nBinAdapt<16) nBinAdapt=16;//set a reasonable minimum
	AdaptBin binning("pullHist",nBinAdapt,mmin,mmax,ipmin,ipmax);
	//binning.setVerbose();//TODO
	binning.loadDataFromHist(dataHist);
	TH2Poly* dataPoly = binning.getHisto("dataPoly");
	TH2Poly* modelPoly = binning.getHisto("modelPoly");
	modelPoly->SetTitle("");
	modelPoly->GetXaxis()->SetTitle("#it{m}("+_decayLabel+") [MeV/#it{c}^{2}]");
	modelPoly->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
	modelPoly->GetZaxis()->SetTitle("");

	//create pull histogram
	modelHist10->Scale(dataHist10->Integral()/modelHist10->Integral());
	for(int ibin=1; ibin<=10; ++ibin) {
		for(int jbin=1; jbin<=10; ++jbin) {
			double D = dataHist10->GetBinContent(ibin,jbin);
			double F = modelHist10->GetBinContent(ibin,jbin);
			double E = dataHist10->GetBinError(ibin,jbin);
			modelHist10->SetBinContent(ibin,jbin,(D-F)/E);
		}
	}

	double x, y, n;
	//fill adaptive binning histograms
	for(uint ibin=1; ibin<=nBins; ++ibin) {
		for(uint jbin=1; jbin<=nBins; ++jbin) {
			x = dataHist->GetXaxis()->GetBinCenter(ibin);
			y = dataHist->GetYaxis()->GetBinCenter(jbin);
			n = dataHist->GetBinContent( ibin, jbin );
			dataPoly->Fill(x, y, n);
			x = modelHist->GetXaxis()->GetBinCenter(ibin);
			y = modelHist->GetYaxis()->GetBinCenter(jbin);
			n = modelHist->GetBinContent( ibin, jbin );
			modelPoly->Fill(x, y, n);
		}
	}
	for(int ibin=1; ibin<=dataPoly->GetNumberOfBins(); ++ibin) {
		double D = dataPoly->GetBinContent(ibin);
		double F = modelPoly->GetBinContent(ibin);
		double E = TMath::Sqrt(D);
		modelPoly->SetBinContent(ibin,(D-F)/E);
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH2* dataHist50 = dataHists[i]->Rebin2D( nBins/nBins1D, nBins/nBins1D,"dataHist50");

		//Mass plots
		TH1* dataHistMass   = dataHist50->ProjectionX("dataHistMass");
		TH1* signalHistMass = signalHists[i]->ProjectionX("signalHistMass");
		TH1* promptHistMass = promptHists[i]->ProjectionX("promptHistMass");
		TH1* displcHistMass = displcHists[i]->ProjectionX("displcHistMass");
		TH1* bkgrndHistMass = bkgrndHists[i]->ProjectionX("bkgrndHistMass");
		TH1* modelHistMass  = modelHists[i]->ProjectionX("modelHistMass");

		signalHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		promptHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		displcHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());

		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (mmax-mmin)/nBins1D));
		bkgrndHistMass->SetLineStyle(kDashed);
		signalHistMass->SetLineStyle(kDotted);
		promptHistMass->SetFillColor(kGreen+2);
		displcHistMass->SetFillColor(kMagenta);
		promptHistMass->SetFillStyle(3345);
		displcHistMass->SetFillStyle(3354);

		dataHistMass->Draw("E1 X0 P");
		signalHistMass->Draw("HIST C SAME");
		promptHistMass->Draw("HIST C SAME");
		displcHistMass->Draw("HIST C SAME");
		bkgrndHistMass->Draw("HIST C SAME");
		modelHistMass->Draw( "HIST C SAME");
		dataHistMass->Draw("E1 X0 P SAME");

		TPaveText label(0.65,0.8,0.9,0.9,"BRNDC");
		if(i==0) {
			label.AddText(TString::Format(
				"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
				_jptmin/1e3,
				_jptmax/1e3));
		} else {
			if(_usePtFracBins) {
				label.AddText(TString::Format(
					"#it{f}_{#it{p}_{T}}^{#it{D}} #in [%.2f,%.2f]",
					static_cast<double>(i-1)/_nDPtBins,
					static_cast<double>(i)/_nDPtBins));
			} else {
				label.AddText(TString::Format(
					"#it{p}_{T}(#it{D}) #in [%.0f,%.0f] GeV/#it{c}",
					_ptBinsD[i-1],
					_ptBinsD[i]));
			}
		}
		label.SetFillColor(0);
		label.SetTextAlign(12);
		label.SetBorderSize(0);
		label.Draw();
		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+"PTbin"+(i-1)+".pdf");

		//IP plots
		TH1* dataHistIP   = dataHist50->ProjectionY("dataHistIP");
		TH1* signalHistIP = signalHists[i]->ProjectionY("signalHistIP");
		TH1* promptHistIP = promptHists[i]->ProjectionY("promptHistIP");
		TH1* displcHistIP = displcHists[i]->ProjectionY("displcHistIP");
		TH1* bkgrndHistIP = bkgrndHists[i]->ProjectionY("bkgrndHistIP");
		TH1* modelHistIP  = modelHists[i]->ProjectionY("modelHistIP");

		signalHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		promptHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		displcHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		bkgrndHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		modelHistIP ->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());

		dataHistIP->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f)", (ipmax-ipmin)/nBins1D));
		bkgrndHistIP->SetLineStyle(kDashed);
		signalHistIP->SetLineStyle(kDotted);
		promptHistIP->SetFillColor(kGreen+2);
		displcHistIP->SetFillColor(kMagenta);
		promptHistIP->SetFillStyle(3345);
		displcHistIP->SetFillStyle(3354);

		dataHistIP->Draw("E1 X0 P");
		signalHistIP->Draw("HIST C SAME");
		promptHistIP->Draw("HIST C SAME");
		displcHistIP->Draw("HIST C SAME");
		bkgrndHistIP->Draw("HIST C SAME");
		modelHistIP->Draw( "HIST C SAME");
		dataHistIP->Draw("E1 X0 P SAME");
		label.Draw();
		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+"PTbin"+(i-1)+".pdf");
	}
	
	//2D plots
	TH2* dataHist30 = dataHist->Rebin2D( nBins/nBins2D, nBins/nBins2D,"dataHist30");
	can.SetTickx();
	can.SetTicky();
	//TODO
	//std::cout << dataHist30->GetMaximum() << " " << dataHist30->GetMinimum() << " " << dataHist30->GetEntries() << std::endl;
	dataHist30->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	promptHist->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	displcHist->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	bkgrndHist->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	//END TODO

	promptHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/promptPDF"+name+".pdf");
	promptHist->Draw("A col");
	can.SaveAs(gSaveDir+"/fig/promptPDF"+name+"_noaxis.png");
	can.SetFillColorAlpha(kWhite,0.);
	promptHist->Draw("AXIS");
	can.SaveAs(gSaveDir+"/fig/promptPDF"+name+"_axis.pdf");
	can.SetFillColor(kWhite);

	displcHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/displcPDF"+name+".pdf");
	displcHist->Draw("A col");
	can.SaveAs(gSaveDir+"/fig/displcPDF"+name+"_noaxis.png");
	can.SetFillColorAlpha(kWhite,0.);
	displcHist->Draw("AXIS");
	can.SaveAs(gSaveDir+"/fig/displcPDF"+name+"_axis.pdf");
	can.SetFillColor(kWhite);

	bkgrndHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/bkgrndPDF"+name+".pdf");
	bkgrndHist->Draw("A col");
	can.SaveAs(gSaveDir+"/fig/bkgrndPDF"+name+"_noaxis.png");
	can.SetFillColorAlpha(kWhite,0.);
	bkgrndHist->Draw("AXIS");
	can.SaveAs(gSaveDir+"/fig/bkgrndPDF"+name+"_axis.pdf");
	can.SetFillColor(kWhite);

	can.SetRightMargin(0.14);
	dataHist30->Draw("colz");
	//signalHist->Draw("cont2 same");
	promptHist->Draw("cont2 same");
	displcHist->Draw("cont2 same");
	bkgrndHist->Draw("cont2 same");

	TPaveText label(0.53,0.8,0.82,0.9,"BRNDC");
	label.AddText(TString::Format(
		"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
		_jptmin/1e3,
		_jptmax/1e3));
	label.SetFillColor(0);
	label.SetTextAlign(12);
	label.SetBorderSize(0);
	label.Draw();
	can.SaveAs(gSaveDir+"/fig/data"+name+".pdf");

	can.SetRightMargin(0.10);
	const Int_t NRGBs2 = 4;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.45, 0.55, 1.00};
	Double_t reds2[NRGBs2]   = { 0.00, 1.00, 1.00, 1.00};
	Double_t greens2[NRGBs2] = { 0.00, 1.00, 1.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 1.00, 1.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);
	modelHist10->SetMinimum(-4.);
	modelHist10->SetMaximum( 4.);
	modelHist10->Draw("colz");
	can.SaveAs(gSaveDir+"/fig/pull"+name+".pdf");
	modelPoly->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	modelPoly->SetMinimum(-4.);
	modelPoly->SetMaximum( 4.);
	modelPoly->Draw("colz");
	can.SaveAs(gSaveDir+"/fig/pull_adapbin"+name+".pdf");
}

bool SimDFitter::doFits() {
	if(!_initialised) init();
	ws->writeToFile("testWS.root");

	bool status(true);
	int nRetries(0);
	int nAttemptsMC(1);
	int nAttemptsData(2);
	if(!_mcFitsDone) {
		do {
			std::cout << "Fitting mass: simulated signal" << std::endl;
			status=fit(fitType::fitMass, fitSet::fitJustSignal,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		//do {
		//	std::cout << "Fitting mass: simulated background" << std::endl;
		//	status=fit(fitType::fitMass, fitSet::fitJustSideband,0,_ptBinsD.size()-1);
		//	if(status) break;
		//	std::cout << "Failed: retrying..." << std::endl;
		//	++nRetries;
		//} while(!status && nRetries<nAttemptsMC);

		if(!status) {
			std::cout << "mass fit failed, skipping IP fit" << std::endl;
			return status;
		}

		do {
			std::cout << "Fitting IP: simulated prompt" << std::endl;
			//for(uint ibin=0; ibin<_ptBinsD.size()-1; ++ibin) {
			//	status=fit(fitType::fitIP, fitSet::fitJustCharm,ibin,ibin+1);
			//	std::cout << ibin << " " << status << std::endl;//TODO
			//}
			status=fit(fitType::fitIP, fitSet::fitJustCharm,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping further fits" << std::endl;
			return status;
		}
		do {
			std::cout << "Fitting IP: simulated displaced" << std::endl;
			status=fit(fitType::fitIP, fitSet::fitJustBeauty,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping further fits" << std::endl;
			return status;
		}
		//do {
		//	std::cout << "Fitting IP: simulated background" << std::endl;
		//	status=fit(fitType::fitIP, fitSet::fitJustSideband,0,_ptBinsD.size()-1);
		//	if(status) break;
		//	std::cout << "Failed: retrying..." << std::endl;
		//	++nRetries;
		//} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping 2D fit" << std::endl;
			return status;
		}
		_mcFitsDone=true;
	}
	nRetries=0;//reset allowed retries for fit to data
	if(_doEnhancedFits) {
		do {
			std::cout << "Fitting 2D: b-enhanced data" << std::endl;
			status=fit(fitType::fit2D, fitSet::fitDataEnhancedB,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << nRetries << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsData);
		do {
			std::cout << "Fitting 2D: c-enhance data" << std::endl;
			status=fit(fitType::fit2D, fitSet::fitDataEnhancedC,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << nRetries << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsData);
	}
	do {
		std::cout << "Fitting 2D: data" << std::endl;
		status=fit(fitType::fit2D, fitSet::fitData,0,_ptBinsD.size()-1);
		if(status) break;
		std::cout << "Failed: retrying..." << nRetries << std::endl;
		++nRetries;
	} while(!status && nRetries<nAttemptsData);
	if(_runToyFits) status&=fit(fitType::fit2D, fitSet::fitToys,0,_ptBinsD.size()-1,100);
	return status;
}

double SimDFitter::getAveWeight(double flavour) {
	TFile* f = TFile::Open(dFileName());
	if(!f) return 1.;

	TTree* t7 = static_cast<TTree*>(f->Get("T7"));
	if(!t7) return 1.;

	double mmin(_mSBLooseLow), mmax(_mSBLooseHigh);
	double ipmin(-5.), ipmax(10.);

	TString cutStr = TString::Format("(JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f)",_jptmin,_jptmax,_ymin,_ymax);
	TString weightStr;

	if(flavour==4) {
		weightStr="(weight4)*";
		ipmax= 3;
	} else if(flavour==5) {
		weightStr="(weight5)*";
		ipmin= 3;
	}
	TH2D sum("sum","",1,mmin,mmax,1,ipmin,ipmax);
	TH2D count("count","",1,mmin,mmax,1,ipmin,ipmax);
	t7->Draw(""+_whichD+"LogIPChi2:"+_whichD+"M>>sum",weightStr+cutStr);
	t7->Draw(""+_whichD+"LogIPChi2:"+_whichD+"M>>count",cutStr);

	std::cout << dFileName() << " " << t7->GetEntries() << std::endl;//TODO
	std::cout << "HERE " << ""+_whichD+"LogIPChi2:"+_whichD+"M>>sum" << " " << cutStr << " " << weightStr << std::endl;//TODO
	std::cout << sum.GetEntries() << " " << sum.GetSumOfWeights() << " " << sum.GetBinContent(1,1) << std::endl;//TODO
	std::cout << count.GetEntries() << " " << count.GetSumOfWeights() << " " << count.GetBinContent(1,1) << std::endl;//TODO

	return sum.GetBinContent(1,1)/count.GetBinContent(1,1);
}

//calculate uncertainty on the average weight
//when the acceptance and selection binning schemes are different need a toy approach
double SimDFitter::getAveWeightError(double flavour) {
	TFile* f = TFile::Open(dFileName());
	if(!f) return 1.;

	TTree* t7 = static_cast<TTree*>(f->Get("T7"));
	if(!t7) return 1.;

	if(flavour!=4 && flavour!=5) return 1.;

	double JetPT(0.), JetEta(0.), DM(0.), DLogIPChi2(0.);
	double weightAcc(0.), werrorAcc(0.), wbinAcc(0.);
	double weightSel(0.), werrorSel(0.), wbinSel(0.);
	int maxAcc(0);
	double mmin(_mSBLooseLow), mmax(_mSBLooseHigh);
	double ipmin(-5.), ipmax(10.);

	t7->SetBranchAddress("JetPT",        &JetPT);
	t7->SetBranchAddress("JetEta",       &JetEta);
	t7->SetBranchAddress(""+_whichD+"M",          &DM);
	t7->SetBranchAddress(""+_whichD+"LogIPChi2",  &DLogIPChi2);

	t7->SetBranchAddress("weightAcc",    &weightAcc);
	t7->SetBranchAddress("werrorAcc",    &werrorAcc);
	t7->SetBranchAddress("wbinAcc",      &wbinAcc);

	maxAcc = t7->GetMaximum("wbinAcc");
	if(flavour==4) {
		t7->SetBranchAddress("weightSel4",   &weightSel);
		t7->SetBranchAddress("werrorSel4",   &werrorSel);
		t7->SetBranchAddress("wbinSel4",     &wbinSel);

		ipmax= 3;
	} else {
		t7->SetBranchAddress("weightSel5",   &weightSel);
		t7->SetBranchAddress("werrorSel5",   &werrorSel);
		t7->SetBranchAddress("wbinSel5",     &wbinSel);

		ipmin= 3;
	}

	//if simple efficiencies then we can skip the toy method an sum total uncertainties in quadrature
	std::cout << "A0" << std::endl;//TODO
	if(_useSimpleEffs) {
	std::cout << "A1" << std::endl;//TODO
		TString cutStr = TString::Format("(JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f)",_jptmin,_jptmax,_ymin,_ymax);
		TString weightStr;

		if(flavour==4) {
			weightStr="(werror4)*(werror4)*";
		} else if(flavour==5) {
			weightStr="(werror5)*(werror5)*";
		}

		TH2D sum("sum","",1,mmin,mmax,1,ipmin,ipmax);
		TH2D count("count","",1,mmin,mmax,1,ipmin,ipmax);
		t7->Draw(""+_whichD+"LogIPChi2:"+_whichD+"M>>sum",weightStr+cutStr);
		t7->Draw(""+_whichD+"LogIPChi2:"+_whichD+"M>>count",cutStr);

	std::cout << "A2 " << sum.GetBinContent(1,1) << " " << count.GetBinContent(1,1) << std::endl;//TODO
		return TMath::Sqrt(sum.GetBinContent(1,1))/count.GetBinContent(1,1);
	}

	//first load information from TTree into maps
	std::map<int, int> counts;
	std::map<int, int> accBins;
	std::map<int, int> selBins;

	std::map<int, double> accWeights;
	std::map<int, double> accErrors;
	std::map<int, double> accPerturbed;

	std::map<int, double> selWeights;
	std::map<int, double> selErrors;
	std::map<int, double> selPerturbed;

	int bin(0);
	boost::progress_display progress_setup( t7->GetEntries() );
	for(int ievt=0; ievt<t7->GetEntries(); ++ievt) {
		++progress_setup;
		t7->GetEntry(ievt);
		if(JetPT<_jptmin || JetPT>=_jptmax) continue;
		if(JetEta<_ymin  || JetEta>=_ymax) continue;
		if(DM<mmin || DM>mmax) continue;
		if(DLogIPChi2<ipmin || DLogIPChi2>=ipmax) continue;

		bin = wbinAcc + wbinSel*maxAcc;
		if(counts.count(bin)==0) {
			counts[bin] = 1;
			accBins[bin] = wbinAcc;
			selBins[bin] = wbinSel;

			if(accWeights.count(wbinAcc)==0) {
				accWeights[wbinAcc] = weightAcc;
				accErrors[wbinAcc]  = werrorAcc;
				accPerturbed[wbinAcc] = 0;
			}

			if(selWeights.count(wbinSel)==0) {
				selWeights[wbinSel] = weightSel;
				selErrors[wbinSel]  = werrorSel;
				selPerturbed[wbinSel] = 0;
			}
		} else {
			counts[bin] += 1;
		}
	}

	double aveWeight = getAveWeight(flavour);
	TH1D toyWeights("toyWeights","",100,0.95*aveWeight,1.05*aveWeight);

	//now throw toys and calculate the weight in each case
	int ntoys=100;
	boost::progress_display progress_run( ntoys );
	for(int itoy=0; itoy<ntoys; ++itoy) {
		++progress_run;
		//generate new pulls for all bins
		for( auto & [bin, val] : accPerturbed ) val = gRandom->Gaus(accWeights[bin],accErrors[bin]);
		for( auto & [bin, val] : selPerturbed ) val = gRandom->Gaus(selWeights[bin],selErrors[bin]);

		double num(0.), denom(0.);
		int binA(0), binS(0);

		for( auto const& [bin, count] : counts ) {
			binA = accBins[bin];
			binS = selBins[bin];
			denom += count;
			num += count*accPerturbed[binA]*selPerturbed[binS];
		}
		toyWeights.Fill(num/denom);
	}
	TCanvas c;
	toyWeights.Draw();
	c.SaveAs(TString::Format(gSaveDir+"/aveDEffWeights_%.0f-%.0f_%.0f-%.0f.pdf",_jptmin,_jptmax,_ymin,_ymax));
	//std::cout << aveWeight << " " << toyWeights.GetMean() << " " << toyWeights.GetStdDev() << std::endl;

	return toyWeights.GetStdDev();
}

void SimDFitter::printYields() {
	for(uint i=0; i<_yields4.size(); ++i) {
		printf("%30s %6.0f +/- %3.0f\n", _names4[i].Data(), _yields4[i], _errors4[i]);
	}
	std::cout << std::endl;
	for(uint i=0; i<_yields5.size(); ++i) {
		printf("%30s %6.0f +/- %3.0f\n", _names5[i].Data(), _yields5[i], _errors5[i]);
	}
}

//int main() {
//	RooMsgService::instance().setSilentMode(true);
//	RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
//	gErrorIgnoreLevel = kError;
//	gSaveDir="simFit/";
//	SimDFitter s;
//	s.doFits();
//}
