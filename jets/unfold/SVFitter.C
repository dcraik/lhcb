#include "SVFitter.h"

#include <vector>
#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMsgService.h"
#include "RooNDKeysPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooThresholdCategory.h"

#include "outputFunctions.h"
#include "DatasetManager.h"

void SVFitter::setInputs(TString light, TString charm, TString beauty, TString data, bool lightIsMC, bool charmIsMC, bool beautyIsMC) {
	_lightInputFile  = light;
	_charmInputFile  = charm;
	_beautyInputFile = beauty;
	_dataInputFile   = data;
	_lightIsMC  = lightIsMC;
	_charmIsMC  = charmIsMC;
	_beautyIsMC = beautyIsMC;

	_inputsSet=true;
}

void SVFitter::setInputWeightings(bool light, bool charm, bool beauty) {
	_lightIsWeighted = light;
	_charmIsWeighted = charm;
	_beautyIsWeighted = beauty;
}

void SVFitter::setOptions(SVFitterOptions& options) {
	_nmbins               = options.nMCorBins;
	_mmin                 = options.minMCor;
	_mmax                 = options.maxMCor;
	_ntbins               = options.nNTrkBins;
	_recreateInputs       = options.rerunTemplates;
	_usePtBinnedTemplates = options.usePtBinnedTemplates;
	_useBinnedTemplates = options.useBinnedTemplates;
}

void SVFitter::setSVBinning(int nmbins, double mmin, double mmax, int ntbins) {
	_nmbins = nmbins;
	_mmin   = mmin;
	_mmax   = mmax;
	_ntbins = ntbins;
}

void SVFitter::makeSVFitHists(int which) {
	//printf("...\t%d\t%f\t%f\t%d\t%d\t%d\n",_nmbins,_mmin,_mmax,_ntbins,2,2+_ntbins);
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return;
	}
	//TFile* f(0);
	TFile* fout(0);
	TString fname;
	DatasetManager* dm = DatasetManager::getInstance();

	bool isMC(false), isWeighted(false);

	if(which==0) {
		std::cout << "INFO : making 2D SV histograms: light" << std::endl;
		dm->setDataset(_lightInputFile);
		fname = templateFileName("0");
		if(_lightIsMC) isMC=true;
		if(_lightIsWeighted) isWeighted=true;
		//fout = TFile::Open(templateFileName("0"),"RECREATE");
	} else if(which==4) {
		std::cout << "INFO : making 2D SV histograms: charm" << std::endl;
		dm->setDataset(_charmInputFile);
		fname = templateFileName("4");
		if(_charmIsMC) isMC=true;
		if(_charmIsWeighted) isWeighted=true;
		//fout = TFile::Open(templateFileName("4"),"RECREATE");
	} else if(which==5) {
		std::cout << "INFO : making 2D SV histograms: beauty" << std::endl;
		dm->setDataset(_beautyInputFile);
		fname = templateFileName("5");
		if(_beautyIsMC) isMC=true;
		if(_beautyIsWeighted) isWeighted=true;
		//fout = TFile::Open(templateFileName("5"),"RECREATE");
	} else if(which==7) {
		std::cout << "INFO : making 2D SV histograms: data" << std::endl;
		dm->setDataset(_dataInputFile);
		fname = templateFileName("D");
		//fout = TFile::Open(templateFileName("D"),"RECREATE");
	}

	if(!_recreateInputs && !gSystem->AccessPathName(fname)) {
		std::cout << "INFO in SVFitter::makeSVFitHists: histograms already made" << std::endl;
		return;
	}

	fout = TFile::Open(fname,"RECREATE");
	TTree* tout = new TTree("T","");
	//if(!f) return;
	//TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	//if(!t) return;

	int npt=1;
	if(_ptBins) {
		npt = _ptBins->GetNbinsX();
	}
	int nzy=1;
	if(_yBins) {
		nzy = _yBins->GetNbinsX();
	}

	std::vector<TH2D*> template2D;
	for(int i=0; i<npt; ++i) {
		for(int j=0; j<nzy; ++j) {
			TString name;
			if(_ptBins) {
				name+="_";
				name+=_ptBins->GetXaxis()->GetBinLowEdge(i+1);
				name+="-";
				name+=_ptBins->GetXaxis()->GetBinUpEdge(i+1);
			} else {
				name+="_0-200000";
			}
			if(_yBins) {
				name+="_";
				name+=_yBins->GetXaxis()->GetBinLowEdge(j+1);
				name+="-";
				name+=_yBins->GetXaxis()->GetBinUpEdge(j+1);
			} else {
				name+="_0-10";
			}
			template2D.push_back(new TH2D("template2D"+name,"",_nmbins,_mmin,_mmax,_ntbins,2,2+_ntbins));
			template2D[i*nzy+j]->Sumw2();
		}
	}

	std::vector<double>* SVN = new std::vector<double>();
	std::vector<double>* SVMCor = new std::vector<double>();
	int NPV;//TODO
	double JetPT;
	//TODO//double JetTRUEc;//TODO debug
	//TODO//double JetTRUEb;//TODO debug
	double ZE, ZPZ;
	double ZY(3.);
	double mcor, ntrk;
	int ntrkI;
	double weight(1.);
	bool vectorInput(true);

	dm->setBranchAddress("NPV", &NPV);//TODO
	dm->setBranchAddress("JetPT", &JetPT);
	//TODO//dm->setBranchAddress("JetTRUEc", &JetTRUEc);//TODO debug
	//TODO//dm->setBranchAddress("JetTRUEb", &JetTRUEb);//TODO debug
	
	//check whether we have SVs stored in vectors or scalars
	TString btype = dm->getBranchType("SVMCor");
	if(btype!="vector<double>") {
		vectorInput=false;
	}

	if(!vectorInput) {
		dm->setBranchAddress("SVMCor", &mcor);
		dm->setBranchAddress("SVN", &ntrk);
	} else {
		dm->setBranchAddress("SVMCor", &SVMCor);
		dm->setBranchAddress("SVN", &SVN);
	}

	if(isMC) {
		if(_yBins) {
			dm->setBranchAddress("ZTRUEE", &ZE);
			dm->setBranchAddress("ZTRUEPZ", &ZPZ);
		}
	} else {
		if(_yBins) {
			dm->setBranchAddress("ZE", &ZE);
			dm->setBranchAddress("ZPZ", &ZPZ);
		}
	}
	if(isWeighted) {
		dm->setBranchAddress("weight", &weight);
	}

	//TODO debug
	//TODO//int countc(0), countb(0);
	//TODO debug

	tout->Branch("JetPT",  &JetPT);
	tout->Branch("ZY",     &ZY);
	tout->Branch("SVMCor", &mcor);
	tout->Branch("SVN",    &ntrkI);
	tout->Branch("weight", &weight);

	int jetPtBin(1), zyBin(1); //defautls to 1 to work when no binning used
	int bin;

	std::cout << dm->getEntries() << std::endl;//TODO
	boost::progress_display progress(dm->getEntries());
	//for(int i=0; i<t->GetEntries(); ++i) 
	while(dm->getNext()) {
		++progress;
		//t->GetEntry(i);
		//if(which==7 && NPV==1) continue;//TODO
		if(_ptBins) {
			jetPtBin=_ptBins->FindBin(JetPT);
			if(jetPtBin<1 || jetPtBin>_ptBins->GetNbinsX()) continue;
		}
		if(_yBins) {
			ZY = 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ));
			zyBin = _yBins->FindBin(ZY);
			if(zyBin<1 || zyBin>_yBins->GetNbinsX()) continue;
		}
		bin=(jetPtBin-1)*nzy+(zyBin-1);
		if(vectorInput) {//extract first SV if our input is in vector form
			if(SVN->size()<1) continue;
			ntrk = SVN->at(0);
			mcor = SVMCor->at(0);
		}
//		if(weightPt) weight = jetTruePtWeights->GetBinContent(jetTruePtWeights->FindBin(JetTruePT));
		//if(mcor>_mmax-0.1) mcor=_mmax-0.1;
		if(mcor>=_mmax) continue;//TODO remove overflow bin
		if(ntrk>1+_ntbins) ntrk=_ntbins+1;
		if(mcor<_mmin || ntrk<2) continue;//second case should be impossible
		//TODO//if(JetTRUEb) ++countb; //TODO debug
		//TODO//if(JetTRUEc) ++countc; //TODO debug
		template2D[bin]->Fill(mcor,ntrk,weight);
		ntrkI = ntrk;
		tout->Fill();
	}

	//TODO//std::cout << countc << " " << countb << std::endl;//TODO debug

	dm->reset();

	for(uint i=0; i<template2D.size(); ++i) {
		template2D[i]->Write();
	}
	tout->Write();
	fout->Close();
}

//function to fit features for a single sample
bool SVFitter::fitSV(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY){
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return false;
	}

	double minPT(0.);
	double maxPT(200000.);
	if(_ptBins!=0) {
		minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
		maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
	}
	double minY(0.), maxY(10.);
	if(_yBins!=0) {
		minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
		maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
	}

	std::cout << "INFO : fitting SV - 1D corrected mass and Ntrk fit for " << _name << ", pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;
	ptStr+="_"; ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* fh0 = TFile::Open(templateFileName("0"));
	TFile* fh4 = TFile::Open(templateFileName("4"));
	TFile* fh5 = TFile::Open(templateFileName("5"));
	TFile* fhd = TFile::Open(templateFileName("D"));
	TH2D* template2D0 = dynamic_cast<TH2D*>(fh0->Get("template2D_"+ptStr));
	TH2D* template2D4 = dynamic_cast<TH2D*>(fh4->Get("template2D_"+ptStr));
	TH2D* template2D5 = dynamic_cast<TH2D*>(fh5->Get("template2D_"+ptStr));
	TH2D* template2Dd = dynamic_cast<TH2D*>(fhd->Get("template2D_"+ptStr));

	if(!template2D0 || !template2D4 || !template2D5 || !template2Dd) return false;
	if(template2Dd->GetEntries()==0) {
		return false;
	}

	int nbinM = template2D0->GetNbinsX();
	int nbinN = template2D0->GetNbinsY();
	int nbin  = nbinM+nbinN;
	double ntrkScale = nbinN/static_cast<double>(nbinM);
	double scale = 1./(1.+ntrkScale);

	TH1D light("light","",nbin,0,nbin);
	TH1D charm("charm","",nbin,0,nbin);
	TH1D beaut("beaut","",nbin,0,nbin);
	TH1D data ("data", "",nbin,0,nbin);

	double nQ, nC, nB, nD;
	for(int i=0; i<nbinM; ++i) {
		for(int j=0; j<nbinN; ++j) {
			nQ = template2D0->GetBinContent(i+1,j+1);
			nC = template2D4->GetBinContent(i+1,j+1);
			nB = template2D5->GetBinContent(i+1,j+1);
			nD = template2Dd->GetBinContent(i+1,j+1);
			light.SetBinContent(i+1,light.GetBinContent(i+1)+nQ);
			charm.SetBinContent(i+1,charm.GetBinContent(i+1)+nC);
			beaut.SetBinContent(i+1,beaut.GetBinContent(i+1)+nB);
			data .SetBinContent(i+1,data .GetBinContent(i+1)+nD);
			light.SetBinContent(nbinM+j+1,light.GetBinContent(nbinM+j+1)+nQ);
			charm.SetBinContent(nbinM+j+1,charm.GetBinContent(nbinM+j+1)+nC);
			beaut.SetBinContent(nbinM+j+1,beaut.GetBinContent(nbinM+j+1)+nB);
			data .SetBinContent(nbinM+j+1,data .GetBinContent(nbinM+j+1)+nD);
		}
	}

	for(int j=0; j<nbinN; ++j) {
		light.SetBinContent(nbinM+j+1,light.GetBinContent(nbinM+j+1)*ntrkScale);
		charm.SetBinContent(nbinM+j+1,charm.GetBinContent(nbinM+j+1)*ntrkScale);
		beaut.SetBinContent(nbinM+j+1,beaut.GetBinContent(nbinM+j+1)*ntrkScale);
		data .SetBinContent(nbinM+j+1,data .GetBinContent(nbinM+j+1)*ntrkScale);
	}

	// -- variables from datasets
	RooRealVar SVComb(  "SVComb",  "SVComb",  0, nbin,  ""); 
	SVComb.setBins(nbin);

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.4*data.Integral(), 0., data.Integral());
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*data.Integral(), 0., data.Integral());
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.3*data.Integral(), 0., data.Integral());

	RooDataHist histSVMSVNB("histCombB", "histCombB", RooArgList(SVComb), &beaut);
	RooDataHist histSVMSVNC("histCombC", "histCombC", RooArgList(SVComb), &charm);
	RooDataHist histSVMSVNQ("histCombQ", "histCombQ", RooArgList(SVComb), &light);

	// -- simulation PDFs for each category
	RooHistPdf pdfB( "pdfB", "pdfB", RooArgSet(SVComb), histSVMSVNB ); 
	RooHistPdf pdfC( "pdfC", "pdfC", RooArgSet(SVComb), histSVMSVNC ); 
	RooHistPdf pdfQ( "pdfQ", "pdfQ", RooArgSet(SVComb), histSVMSVNQ ); 

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(pdfB, pdfC, pdfQ), RooArgList(yieldB, yieldC, yieldQ) );

	// -- add all feature observables to dataset
	RooArgSet obs;
	obs.add(SVComb);

	RooDataHist dh("dh", "dh", RooArgList(SVComb), &data);

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/SVComb_"+_name+"_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
	gSystem->RedirectOutput(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	TString plotName = "SVComb_"+_name+"_"; plotName+=ptStr;
	plotFit(SVComb, 0, nbin, nbin, &dh, data_pdf, sig_pdfs, bkg_pdfs, plotName, "M_{cor}, N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_SVcomb_"+_name+"_"; paramsName+=ptStr; paramsName+=".dat";
	printParams(paramsName,params);

	double Ntot;

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV()*scale;
	NC = yieldC.getValV()*scale;
	NQ = yieldQ.getValV()*scale;
	eB = yieldB.getError()*scale;
	eC = yieldC.getError()*scale;
	eQ = yieldQ.getError()*scale;
	Ntot*=scale;
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ);

	TString saveName=gSaveDir+"/SVFitHists";
	saveName+=_name+"_";
	saveName+=ptStr;
	saveName+=".root";
	TFile* fsave = TFile::Open(saveName,"RECREATE");
	light.SetDirectory(fsave);
	light.Write();
	charm.SetDirectory(fsave);
	charm.Write();
	beaut.SetDirectory(fsave);
	beaut.Write();
	data .SetDirectory(fsave);
	data .Write();
	fsave->Close();

	return true;
}

//function to fit features for a single sample
bool SVFitter::fitSVSim(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY){
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return false;
	}
	
	double minPT(0.);
	double maxPT(200000.);
	if(_ptBins!=0) {
		minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
		maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
	}
	double minY(0.), maxY(10.);
	if(_yBins!=0) {
		minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
		maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
	}

	std::cout << "INFO : fitting SV - binned 2D corrected mass and Ntrk fit for " << _name << ", pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;
	ptStr+="_"; ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* fh0 = TFile::Open(templateFileName("0"));
	TFile* fh4 = TFile::Open(templateFileName("4"));
	TFile* fh5 = TFile::Open(templateFileName("5"));
	TFile* fhd = TFile::Open(templateFileName("D"));
	TTree* t0 = dynamic_cast<TTree*>(fh0->Get("T"));
	TTree* t4 = dynamic_cast<TTree*>(fh4->Get("T"));
	TTree* t5 = dynamic_cast<TTree*>(fh5->Get("T"));
	TTree* td = dynamic_cast<TTree*>(fhd->Get("T"));

	if(!t0 || !t4 || !t5) return false;
	if(!td) return false;
	if(td->GetEntries()==0) {
		return false;
	}
	TH1D hSVN0("hSVN0","",_ntbins,2.,_ntbins+2.);
	TH1D hSVN4("hSVN4","",_ntbins,2.,_ntbins+2.);
	TH1D hSVN5("hSVN5","",_ntbins,2.,_ntbins+2.);
	TH1D hSVND("hSVND","",_ntbins,2.,_ntbins+2.);
	TH1D hSVNF("hSVNF","",_ntbins,2.,_ntbins+2.);

	TString cutStrData = TString::Format("JetPT>=%f && JetPT<%f && ZY>=%f && ZY<%f && SVMCor>=%f && SVMCor<%f", minPT, maxPT, minY, maxY, _mmin, _mmax);
	TString cutStrTemplates = TString::Format("JetPT>=%f && JetPT<%f && ZY>=%f && ZY<%f && SVMCor>=%f && SVMCor<%f", 10e3, 100e3, 0., 10., _mmin, _mmax);
	if(_usePtBinnedTemplates) cutStrTemplates = cutStrData;
	//TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && SVMCor>=%f && SVMCor<%f", minPT, maxPT, _mmin, _mmax);
	//std::cout << cutStr << std::endl;
	t0->Draw("SVN>>hSVN0",cutStrTemplates);
	t4->Draw("SVN>>hSVN4",cutStrTemplates);
	t5->Draw("SVN>>hSVN5",cutStrTemplates);
	td->Draw("SVN>>hSVND",cutStrData);

	//std::cout << hSVN0.GetEntries() << " " << hSVN4.GetEntries() << " " << hSVN5.GetEntries() << " " << hSVND.GetEntries() << std::endl;

	hSVN0.Scale(1./hSVN0.Integral());
	hSVN4.Scale(1./hSVN4.Integral());
	hSVN5.Scale(1./hSVN5.Integral());

	// -- variables from datasets
	RooRealVar JetPT("JetPT",   "JetPT",   1e4, 1e5,  "MeV/#it{c}"); 
	RooRealVar ZY   ("ZY",      "ZY",      0., 10.,   ""); 
	RooRealVar  SVM ("SVMCor",  "SVMCor",  _mmin, _mmax,  "MeV/#it{c}^{2}"); 
	SVM.setBins(_nmbins);
	//RooRealVar  SVN ("SVN",     "SVN",     2, _ntbins+2,  ""); 
	//SVN.setBins(_ntbins);
	RooRealVar weight("weight",  "",  0., 2.); 

	RooCategory SVN(  "SVN",  "SVN");
	for(int i=2; i<_ntbins+2; ++i) {
		//std::cout << TString::Format("SVN%d",i) << " " << i << std::endl;//TODO
		SVN.defineType(TString::Format("SVN%d",i),i);
	}

	RooDataSet ds("ds", "ds", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*td), RooFit::WeightVar(weight), RooFit::Cut(cutStrData));
	RooDataSet d0("d0", "d0", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t0), RooFit::WeightVar(weight), RooFit::Cut(cutStrTemplates));
	RooDataSet d4("d4", "d4", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t4), RooFit::WeightVar(weight), RooFit::Cut(cutStrTemplates));
	RooDataSet d5("d5", "d5", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t5), RooFit::WeightVar(weight), RooFit::Cut(cutStrTemplates));

	//std::cout << d0.sumEntries() << " " << d4.sumEntries() << " " << d5.sumEntries() << " " << ds.sumEntries() << std::endl;
	
	double Ntot = ds.sumEntries();

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.5*Ntot, -0.02*Ntot, Ntot);
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*Ntot, -0.02*Ntot, Ntot);
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.2*Ntot, -0.02*Ntot, Ntot);

	RooSimultaneous data_pdf("data_pdf","data_pdf",SVN);
	for(int i=0; i<_ntbins; ++i) {
		//slices of the datasets
		RooDataSet* d0slice = static_cast<RooDataSet*>(d0.reduce(RooFit::Name(TString::Format("d0_%d",i+2)),RooFit::WeightVar(weight),RooFit::Cut(TString::Format("SVN==%d",i+2))));
		RooDataSet* d4slice = static_cast<RooDataSet*>(d4.reduce(RooFit::Name(TString::Format("d4_%d",i+2)),RooFit::WeightVar(weight),RooFit::Cut(TString::Format("SVN==%d",i+2))));
		RooDataSet* d5slice = static_cast<RooDataSet*>(d5.reduce(RooFit::Name(TString::Format("d5_%d",i+2)),RooFit::WeightVar(weight),RooFit::Cut(TString::Format("SVN==%d",i+2))));

		RooAbsPdf *pdfQ(0), *pdfC(0), *pdfB(0);
		if(_useBinnedTemplates) {
			//create binned datasets
			RooDataHist* d0sliceHist = new RooDataHist(TString::Format("d0Hist_%d",i+2),"",RooArgList(SVM),*d0slice);
			RooDataHist* d4sliceHist = new RooDataHist(TString::Format("d4Hist_%d",i+2),"",RooArgList(SVM),*d4slice);
			RooDataHist* d5sliceHist = new RooDataHist(TString::Format("d5Hist_%d",i+2),"",RooArgList(SVM),*d5slice);
			//create histogram PDFs
			pdfQ = new RooHistPdf(TString::Format("pdfQ_%d",i+2),"",SVM,*d0sliceHist);
			pdfC = new RooHistPdf(TString::Format("pdfC_%d",i+2),"",SVM,*d4sliceHist);
			pdfB = new RooHistPdf(TString::Format("pdfB_%d",i+2),"",SVM,*d5sliceHist);
		} else {
			//create KDE PDFs
			pdfQ = new RooKeysPdf(TString::Format("pdfQ_%d",i+2),"",SVM,*d0slice);
			pdfC = new RooKeysPdf(TString::Format("pdfC_%d",i+2),"",SVM,*d4slice);
			pdfB = new RooKeysPdf(TString::Format("pdfB_%d",i+2),"",SVM,*d5slice);
		}

		//determine yield fractions
		RooRealVar* fracQ = new RooRealVar(TString::Format("fracQ%d",i+2),"",d0slice->sumEntries()/d0.sumEntries());
		RooRealVar* fracC = new RooRealVar(TString::Format("fracC%d",i+2),"",d4slice->sumEntries()/d4.sumEntries());
		RooRealVar* fracB = new RooRealVar(TString::Format("fracB%d",i+2),"",d5slice->sumEntries()/d5.sumEntries());

		RooFormulaVar* partYieldQ = new RooFormulaVar(TString::Format("partYieldQ%d",i+2),"","@0*@1",RooArgList(yieldQ,*fracQ));
		RooFormulaVar* partYieldC = new RooFormulaVar(TString::Format("partYieldC%d",i+2),"","@0*@1",RooArgList(yieldC,*fracC));
		RooFormulaVar* partYieldB = new RooFormulaVar(TString::Format("partYieldB%d",i+2),"","@0*@1",RooArgList(yieldB,*fracB));

		//create total PDF
		RooAddPdf* pdf = new RooAddPdf( TString::Format("data_pdf_%d",i+2),  "", RooArgList(*pdfB, *pdfC, *pdfQ), RooArgList(*partYieldB, *partYieldC, *partYieldQ) );
		data_pdf.addPdf(*pdf, TString::Format("SVN%d",i+2));
		//std::cout << data_pdf.getPdf(TString::Format("SVN%d",i+2))->GetName() << std::endl;//TODO
	}

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/SV2D_"+_name+"_"+ptStr+"_fits.log","w");
	RooFitResult* r(0);
	int attempt(0), nattempts(1);
	while(true) {
		std::cout << "ATTEMPT " << attempt << std::endl;
		std::cout << "STARTING NQ=" << yieldQ.getVal() << ", NC=" << yieldC.getVal() << ", NB=" << yieldB.getVal() << std::endl;
		if(r) delete r;
		r = data_pdf.fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
		std::cout << "STATUS=" << r->status() << ", COVQUAL=" << r->covQual() << std::endl;
		if(r->status()==0 && r->covQual()==3) break;
		++attempt;
		if(attempt>=nattempts) break;

		//otherwise rethrow starting values
		double a = gRandom->Rndm();
		double b = gRandom->Rndm();
		double c = gRandom->Rndm();

		yieldQ.setVal(Ntot*a/(a+b+c));
		yieldB.setVal(Ntot*b/(a+b+c));
		yieldC.setVal(Ntot*c/(a+b+c));
	}
	gSystem->RedirectOutput(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ*" );
	sig_pdfs.push_back( "pdfB*" );
	sig_pdfs.push_back( "pdfC*" );
	std::vector<std::string> bkg_pdfs;

	plotFit(SVM, _mmin, _mmax, _nmbins, &ds, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+_name+"_"+ptStr, "M_{cor}");
	//plotFit(SVN, 2, 2+_ntbins, _ntbins, &ds, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_SV2D_"+_name+"_"; paramsName+=ptStr; paramsName+=".dat";
	printParams(paramsName,params);

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV();
	NC = yieldC.getValV();
	NQ = yieldQ.getValV();
	eB = yieldB.getError();
	eC = yieldC.getError();
	eQ = yieldQ.getError();
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ);

	hSVN0.Scale(NQ);
	hSVN4.Scale(NC);
	hSVN5.Scale(NB);
	hSVNF.Add(&hSVN0);
	hSVNF.Add(&hSVN4);
	hSVNF.Add(&hSVN5);

	hSVN0.SetLineColor(kBlue);
	hSVN4.SetLineColor(kGreen+2);
	hSVN5.SetLineColor(kRed);
	hSVNF.SetLineColor(kBlack);
	hSVND.SetMarkerStyle(kFullCircle);
	hSVND.SetLineColor(kBlack);

	hSVNF.SetMinimum(0.);
	hSVNF.SetMaximum(1.1*TMath::Max(hSVNF.GetMaximum(),hSVND.GetMaximum()));
	hSVNF.GetXaxis()->SetTitle("N_{trk}");
	hSVNF.GetYaxis()->SetTitle("Events");
	hSVNF.SetTitle("");
	hSVNF.GetXaxis()->SetLabelOffset(0.02);
	hSVNF.GetXaxis()->SetTitleOffset(1.18);

	hSVNF.GetXaxis()->SetBinLabel(1,"2");
	hSVNF.GetXaxis()->SetBinLabel(2,"3");
	hSVNF.GetXaxis()->SetBinLabel(3,"4+");

	TCanvas c;
	hSVNF.Draw("HIST");
	hSVN5.Draw("HIST SAME");
	hSVN4.Draw("HIST SAME");
	hSVN0.Draw("HIST SAME");
	hSVND.Draw("SAME E1 X0 P");
	c.SaveAs(gSaveDir+"/SVN_"+_name+"_"+ptStr+".pdf");

	return true;
}

//function to fit features for a single sample
bool SVFitter::fitSV2D(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY){
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return false;
	}
	
	double minPT(0.);
	double maxPT(200000.);
	if(_ptBins!=0) {
		minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
		maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
	}
	double minY(0.), maxY(10.);
	if(_yBins!=0) {
		minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
		maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
	}

	std::cout << "INFO : fitting SV - binned 2D corrected mass and Ntrk fit for " << _name << ", pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;
	ptStr+="_"; ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* fh0 = TFile::Open(templateFileName("0"));
	TFile* fh4 = TFile::Open(templateFileName("4"));
	TFile* fh5 = TFile::Open(templateFileName("5"));
	TFile* fhd = TFile::Open(templateFileName("D"));
	TH2D* template2D0 = dynamic_cast<TH2D*>(fh0->Get("template2D_"+ptStr));
	TH2D* template2D4 = dynamic_cast<TH2D*>(fh4->Get("template2D_"+ptStr));
	TH2D* template2D5 = dynamic_cast<TH2D*>(fh5->Get("template2D_"+ptStr));
	TH2D* template2Dd = dynamic_cast<TH2D*>(fhd->Get("template2D_"+ptStr));

	if(!template2D0 || !template2D4 || !template2D5 || !template2Dd) return false;
	if(template2Dd->GetEntries()==0) {
		return false;
	}

	// -- variables from datasets
	RooRealVar SVM(  "SVM",  "SVM",  _mmin, _mmax,  ""); 
	RooRealVar SVN(  "SVN",  "SVN",  2, 2+_ntbins,  ""); 
	SVM.setBins(_nmbins);
	SVN.setBins(_ntbins);
	
	double Ntot = template2Dd->Integral();

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.5*Ntot, 0., Ntot);
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*Ntot, 0., Ntot);
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.2*Ntot, 0., Ntot);

	RooDataHist histSVMSVNB("histB", "histB", RooArgList(SVM,SVN), template2D5);
	RooDataHist histSVMSVNC("histC", "histC", RooArgList(SVM,SVN), template2D4);
	RooDataHist histSVMSVNQ("histQ", "histQ", RooArgList(SVM,SVN), template2D0);

	// -- simulation PDFs for each category
	RooHistPdf pdfB( "pdfB", "pdfB", RooArgSet(SVM,SVN), histSVMSVNB ); 
	RooHistPdf pdfC( "pdfC", "pdfC", RooArgSet(SVM,SVN), histSVMSVNC ); 
	RooHistPdf pdfQ( "pdfQ", "pdfQ", RooArgSet(SVM,SVN), histSVMSVNQ ); 

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(pdfB, pdfC, pdfQ), RooArgList(yieldB, yieldC, yieldQ) );

	// -- add all feature observables to dataset
	RooArgSet obs;
	obs.add(SVM);
	obs.add(SVN);

	RooDataHist dh("dh", "dh", RooArgList(SVM,SVN), template2Dd);

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/SV2D_"+_name+"_"+ptStr+"_fits.log","w");
	RooFitResult* r(0);
	int attempt(0), nattempts(1);
	while(true) {
		std::cout << "ATTEMPT " << attempt << std::endl;
		std::cout << "STARTING NQ=" << yieldQ.getVal() << ", NC=" << yieldC.getVal() << ", NB=" << yieldB.getVal() << std::endl;
		if(r) delete r;
		r = data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
		std::cout << "STATUS=" << r->status() << ", COVQUAL=" << r->covQual() << std::endl;
		if(r->status()==0 && r->covQual()==3) break;
		++attempt;
		if(attempt>=nattempts) break;

		//otherwise rethrow starting values
		double a = gRandom->Rndm();
		double b = gRandom->Rndm();
		double c = gRandom->Rndm();

		yieldQ.setVal(Ntot*a/(a+b+c));
		yieldB.setVal(Ntot*b/(a+b+c));
		yieldC.setVal(Ntot*c/(a+b+c));
	}
	gSystem->RedirectOutput(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	plotFit(SVM, _mmin, _mmax, _nmbins, &dh, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+_name+"_"+ptStr, "M_{cor}");
	plotFit(SVN, 2, 2+_ntbins, _ntbins, &dh, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_SV2D_"+_name+"_"; paramsName+=ptStr; paramsName+=".dat";
	printParams(paramsName,params);

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV();
	NC = yieldC.getValV();
	NQ = yieldQ.getValV();
	eB = yieldB.getError();
	eC = yieldC.getError();
	eQ = yieldQ.getError();
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ);

	return true;
}

//function to fit features for a single sample
bool SVFitter::fitSV2DUB(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY){
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return false;
	}

	double minPT(0.);
	double maxPT(200000.);
	if(_ptBins!=0) {
		minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
		maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
	}
	double minY(0.), maxY(10.);
	if(_yBins!=0) {
		minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
		maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
	}

	std::cout << "INFO : fitting SV - unbinned 2D corrected mass and Ntrk fit for " << _name << ", pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;
	ptStr+="_"; ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* fh0 = TFile::Open(templateFileName("0"));
	TFile* fh4 = TFile::Open(templateFileName("4"));
	TFile* fh5 = TFile::Open(templateFileName("5"));
	TFile* fhd = TFile::Open(templateFileName("D"));
	TTree* t0 = dynamic_cast<TTree*>(fh0->Get("T"));
	TTree* t4 = dynamic_cast<TTree*>(fh4->Get("T"));
	TTree* t5 = dynamic_cast<TTree*>(fh5->Get("T"));
	TTree* td = dynamic_cast<TTree*>(fhd->Get("T"));

	if(!t0 || !t4 || !t5 || !td) return false;
	if(td->GetEntries()==0) {
		return false;
	}

	// -- variables from datasets
	RooRealVar JetPT("JetPT",   "JetPT",   minPT, maxPT,  "");
	RooRealVar ZY   ("ZY",      "ZY",      minY,  maxY,   "");
	RooRealVar SVM  ("SVMCor",  "SVMCor",   _mmin, _mmax,  ""); 
	RooRealVar SVN  ("SVN",     "SVN",     2, 2+_ntbins,  ""); 
	SVM.setBins(_nmbins);
	SVN.setBins(_ntbins);

	RooDataSet light ("light",  "", RooArgList(JetPT,ZY,SVM,SVN), RooFit::Import(*t0));
	RooDataSet charm ("charm",  "", RooArgList(JetPT,ZY,SVM,SVN), RooFit::Import(*t4));
	RooDataSet beauty("beauty", "", RooArgList(JetPT,ZY,SVM,SVN), RooFit::Import(*t5));
	RooDataSet data  ("data",   "", RooArgList(JetPT,ZY,SVM,SVN), RooFit::Import(*td));
	
	double Ntot = data.sumEntries();

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.5*Ntot, 0., Ntot);
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*Ntot, 0., Ntot);
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.2*Ntot, 0., Ntot);

	RooNDKeysPdf pdfB("pdfB", "", RooArgList(SVM,SVN), beauty);
	RooNDKeysPdf pdfC("pdfC", "", RooArgList(SVM,SVN), charm);
	RooNDKeysPdf pdfQ("pdfQ", "", RooArgList(SVM,SVN), light);

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(pdfB, pdfC, pdfQ), RooArgList(yieldB, yieldC, yieldQ) );

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/SV2DUB_"+_name+"_"+ptStr+"_fits.log","w");
	RooFitResult* r(0);
	int attempt(0), nattempts(1);
	while(true) {
		std::cout << "ATTEMPT " << attempt << std::endl;
		std::cout << "STARTING NQ=" << yieldQ.getVal() << ", NC=" << yieldC.getVal() << ", NB=" << yieldB.getVal() << std::endl;
		if(r) delete r;
		r = data_pdf.fitTo( data, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
		std::cout << "STATUS=" << r->status() << ", COVQUAL=" << r->covQual() << std::endl;
		if(r->status()==0 && r->covQual()==3) break;
		++attempt;
		if(attempt>=nattempts) break;

		//otherwise rethrow starting values
		double a = gRandom->Rndm();
		double b = gRandom->Rndm();
		double c = gRandom->Rndm();

		yieldQ.setVal(Ntot*a/(a+b+c));
		yieldB.setVal(Ntot*b/(a+b+c));
		yieldC.setVal(Ntot*c/(a+b+c));
	}
	gSystem->RedirectOutput(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	plotFit(SVM, _mmin, _mmax, _nmbins, &data, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+_name+"_"+ptStr, "M_{cor}");
	plotFit(SVN, 2, 2+_ntbins, _ntbins, &data, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_SV2DUB_"+_name+"_"; paramsName+=ptStr; paramsName+=".dat";
	printParams(paramsName,params);

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV();
	NC = yieldC.getValV();
	NQ = yieldQ.getValV();
	eB = yieldB.getError();
	eC = yieldC.getError();
	eQ = yieldQ.getError();
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ);

	return true;
}

TString SVFitter::templateFileName(TString which) {
	return gSaveDir+ "/svFitTemplates_"+_name+which+".root";
}
