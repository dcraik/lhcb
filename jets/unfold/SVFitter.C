#include "SVFitter.h"

#include <vector>
#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooRealVar.h"

#include "outputFunctions.h"

void SVFitter::setInputs(TString light, TString charm, TString beauty, TString data) {
	_lightInputFile  = light;
	_charmInputFile  = charm;
	_beautyInputFile = beauty;
	_dataInputFile   = data;

	_inputsSet=true;
}

void SVFitter::setInputWeightings(bool light, bool charm, bool beauty) {
	_lightIsWeighted = light;
	_charmIsWeighted = charm;
	_beautyIsWeighted = beauty;
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
	TFile* f(0);
	TFile* fout(0);

	if(which==0) {
		std::cout << "INFO : making 2D SV histograms: light" << std::endl;
		f = TFile::Open(_lightInputFile);
		fout = TFile::Open(templateFileName("0"),"RECREATE");
	} else if(which==4) {
		std::cout << "INFO : making 2D SV histograms: charm" << std::endl;
		f = TFile::Open(_charmInputFile);
		fout = TFile::Open(templateFileName("4"),"RECREATE");
	} else if(which==5) {
		std::cout << "INFO : making 2D SV histograms: beauty" << std::endl;
		f = TFile::Open(_beautyInputFile);
		fout = TFile::Open(templateFileName("5"),"RECREATE");
	} else if(which==7) {
		std::cout << "INFO : making 2D SV histograms: data" << std::endl;
		f = TFile::Open(_dataInputFile);
		fout = TFile::Open(templateFileName("D"),"RECREATE");
	}
	if(!f) return;
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

	if(!t) return;

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
	double ZE, ZPZ;
	double mcor, ntrk;
	double weight(1.);
	bool vectorInput(true);

	t->SetBranchAddress("NPV", &NPV);//TODO
	t->SetBranchAddress("JetPT", &JetPT);
	if(which==4 || which==5) {
		t->SetBranchAddress("SVMCor", &mcor);
		t->SetBranchAddress("SVN", &ntrk);
		if(_yBins) {
			t->SetBranchAddress("ZTRUEE", &ZE);
			t->SetBranchAddress("ZTRUEPZ", &ZPZ);
		}
		t->SetBranchAddress("weight", &weight);
		vectorInput=false;
	} else {
		t->SetBranchAddress("SVMCor", &SVMCor);
		t->SetBranchAddress("SVN", &SVN);
		if(_yBins) {
			t->SetBranchAddress("ZE", &ZE);
			t->SetBranchAddress("ZPZ", &ZPZ);
		}
	}

	double ZY;
	int jetPtBin(1), zyBin(1); //defautls to 1 to work when no binning used
	int bin;

	boost::progress_display progress(t->GetEntries());
	for(int i=0; i<t->GetEntries(); ++i) {
		++progress;
		t->GetEntry(i);
		//if(which==7 && NPV==1) continue;//TODO
		if(_ptBins) {
			jetPtBin=_ptBins->FindBin(JetPT);
			if(jetPtBin<1 || jetPtBin>_ptBins->GetNbinsX()) continue;
		}
		ZY = 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ));
		if(_yBins) {
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
		if(mcor>_mmax-0.1) mcor=_mmax-0.1;
		if(ntrk>1+_ntbins) ntrk=_ntbins+1;
		if(mcor<_mmin || ntrk<2) continue;//second case should be impossible
		template2D[bin]->Fill(mcor,ntrk,weight);
	}

	for(uint i=0; i<template2D.size(); ++i) {
		template2D[i]->Write();
	}
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

TString SVFitter::templateFileName(TString which) {
	return gSaveDir+ "/svFitTemplates_"+_name+which+".root";
}
