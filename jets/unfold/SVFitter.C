#include "SVFitter.h"

#include <vector>
#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMCStudy.h"
#include "RooModulatedHistPdf.h"
#include "RooMsgService.h"
#include "RooNDKeysPdf.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooThresholdCategory.h"

#include "outputFunctions.h"
#include "DatasetManager.h"

void SVFitter::setInputs(TString light, TString charm, TString beauty, TString data, bool lightIsMC, bool charmIsMC, bool beautyIsMC, bool dataIsMC, TString training, TString back) {
	_lightInputFile  = light;
	_charmInputFile  = charm;
	_beautyInputFile = beauty;
	_dataInputFile   = data;
	_trainingInputFile = training;
	_backwardsInputFile = back;
	_lightIsMC  = lightIsMC;
	_charmIsMC  = charmIsMC;
	_beautyIsMC = beautyIsMC;
	_dataIsMC  = dataIsMC;

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
	_useBinnedTemplates   = options.useBinnedTemplates;
	_smearMCorShapes      = options.smearMCorShapes;
	_additiveMCorCorrectionFactor = options.additiveMCorCorrectionFactor;
	_lightYieldFloat      = options.lightYieldFloat;
	_lightYieldScale      = options.lightYieldScale;
	_runToyFits           = options.runToyFits; 
}

void SVFitter::setSVBinning(int nmbins, double mmin, double mmax, int ntbins) {
	_nmbins = nmbins;
	_mmin   = mmin;
	_mmax   = mmax;
	_ntbins = ntbins;
}

void SVFitter::makeSVFitHists() {
	makeSVFitHists(0);
	makeSVFitHists(4);
	makeSVFitHists(5);
	makeSVFitHists(7);
	makeSVFitHists(8);
	makeSVFitHists(9);
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
		if(_dataIsMC) isMC=true;
		//fout = TFile::Open(templateFileName("D"),"RECREATE");
	} else if(which==8) {
		if(_trainingInputFile=="") {
			std::cout << "INFO : skipping 2D SV histogram: training data - will use data file instead" << std::endl;
			return;
		}
		std::cout << "INFO : making 2D SV histograms: training data" << std::endl;
		dm->setDataset(_trainingInputFile);
		fname = templateFileName("DE");
	} else if(which==9) {
		if(_backwardsInputFile=="") {
			std::cout << "INFO : skipping 2D SV histogram: backwards-tagged data" << std::endl;
			return;
		}
		std::cout << "INFO : making 2D SV histograms: backwards-tagged data" << std::endl;
		dm->setDataset(_backwardsInputFile);
		fname = templateFileName("DB");
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
	//TODO don't bin in rapidity for mistag
	if(_yBins && which!=0) {
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
			if(_yBins && which!=0) {
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
	std::vector<double>* TSVN = new std::vector<double>();
	std::vector<double>* TSVMCor = new std::vector<double>();
	int NPV;//TODO
	double JetPT;
	//TODO//double JetTRUEc;//TODO debug
	//TODO//double JetTRUEb;//TODO debug
	double ZE, ZPZ;
	double ZY(3.);
	double mcor, ntrk;
	int ntrkI;
	double weight(1.);
	double enhanced(0.);
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
		dm->setBranchAddress("TSVMCor", &TSVMCor);
		dm->setBranchAddress("TSVN", &TSVN);
	}

	if(isMC) {
		if(_yBins && which!=0) {
			dm->setBranchAddress("ZTRUEE", &ZE);
			dm->setBranchAddress("ZTRUEPZ", &ZPZ);
		}
	} else {
		if(_yBins && which!=0) {
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
	tout->Branch("enhanced", &enhanced);

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
		if(_yBins && which!=0) {
			ZY = 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ));
			zyBin = _yBins->FindBin(ZY);
			if(zyBin<1 || zyBin>_yBins->GetNbinsX()) continue;
		}
		bin=(jetPtBin-1)*nzy+(zyBin-1);
		//if(which==0) std::cout << "!" << jetPtBin << " " << zyBin << " " << npt << " " << nzy << " " << bin << std::endl;//TODO
		enhanced = 0;
		if(vectorInput) {//extract first SV if our input is in vector form
			if(SVN->size()<1) continue;
			ntrk = SVN->at(0);
			mcor = SVMCor->at(0);
			if(TSVN->size()>0) {
				if(TSVN->at(0)>2 && TSVMCor->at(0)>2000) enhanced = 5;
				if(TSVN->at(0)==2 && TSVMCor->at(0)<2000) enhanced = 4;
			}
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

bool SVFitter::fit(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY, TString enhancedTemplatesFile) {
	_finalRun = false;
	if(enhancedTemplatesFile!="" && loadedTemplateFile!=enhancedTemplatesFile) {
		_pdfCache.clear();
		_fracYieldsCache.clear();

		std::map<std::string,TH1*>* pdfMap = new std::map<std::string,TH1*>();
		std::map<std::string,double>* yieldMap = new std::map<std::string,double>();
		_pdfCache.push_back(pdfMap);
		_fracYieldsCache.push_back(yieldMap);

		std::cout << "INFO in SVFitter::fit : Using enhanced templates from " << enhancedTemplatesFile << std::endl;
		//fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY, 0, "noEnhancement");
		TFile* fin = TFile::Open(enhancedTemplatesFile);
		if(fin) {
			TH2D* hC = static_cast<TH2D*>(fin->Get("hC"));
			TH2D* hB = static_cast<TH2D*>(fin->Get("hB"));
			if(hC && hB) {
				for(uint nbin=0; nbin<_ntbins; ++nbin) {
					TH1 *hC1(0), *hB1(0);
					double fyC(-1.), fyB(-1.);
					TString nameC = TString::Format("pdfC_%d",nbin+2);
					TString nameB = TString::Format("pdfB_%d",nbin+2);

					//if no PDF map or any PDFs are already in the back map then make a new map
					if(pdfMap->count(nameC.Data()) || pdfMap->count(nameB.Data())) {
						pdfMap = new std::map<std::string,TH1*>();
						_pdfCache.push_back(pdfMap);
						yieldMap = new std::map<std::string,double>();
						_fracYieldsCache.push_back(yieldMap);
					}

					hC1 = hC->ProjectionX(nameC,nbin+1,nbin+1);
					hB1 = hB->ProjectionX(nameB,nbin+1,nbin+1);

					fyC = hC1->Integral() / hC->Integral();
					fyB = hB1->Integral() / hB->Integral();
					
					pdfMap->insert(std::pair<TString,TH1*>(nameC.Data(),hC1));
					yieldMap->insert(std::pair<TString,double>(nameC.Data(),fyC));
					pdfMap->insert(std::pair<TString,TH1*>(nameB.Data(),hB1));
					yieldMap->insert(std::pair<TString,double>(nameB.Data(),fyB));
				}
			} else {
				std::cout << "WARNING in SVFitter::fit : could not load enhanced templates" << std::endl
					  << "                             performing fit without enhancement" << std::endl;
			}
	//		fin->Close();
		} else {
			std::cout << "WARNING in SVFitter::fit : could not load file for enhanced templates" << std::endl
				  << "                             performing fit without enhancement" << std::endl;
		}
		for(uint nbin=0; nbin<_ntbins; ++nbin) {
			TString nameC = TString::Format("pdfC_%d",nbin+2);
			TString nameB = TString::Format("pdfB_%d",nbin+2);

			std::cout << _pdfCache.back()->at(nameC.Data())->Integral() << " " << _pdfCache.back()->at(nameB.Data())->Integral() << std::endl;//TODO
			std::cout << _fracYieldsCache.back()->at(nameC.Data()) << " " << _fracYieldsCache.back()->at(nameB.Data()) << std::endl;//TODO
		}
		loadedTemplateFile = enhancedTemplatesFile;
	}
	_finalRun = true;
	return fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY);
}

bool SVFitter::enhanceTemplatesAndFit(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY, uint nPass) {
	_finalRun = false;
	std::cout << "INFO in SVFitter::enhanceTemplates : Starting iterative correction for c/b templates" << std::endl;
	double NCtest, NBtest;
	double NCenhanceTest, NBenhanceTest;
	//_cShiftVals.clear();
	//_bShiftVals.clear();
	_fracYieldsCache.clear();
	_pdfCache.clear();

	fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY, 0, "noEnhancement");
	NCtest = NC;
	NBtest = NB;

	for(int iPass=0; iPass<nPass; ++iPass) {
		TString passName = TString::Format("pass%d",iPass);

		if(!fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY, 5, passName)) return false;
		//for(auto shift : _bShiftVals.back()) std::cout << shift << " ";
		//std::cout << std::endl;
		if(TMath::Abs(NBenhanceTest-NB)/NB < 0.01) {
			std::cout << "b results stablised " << NB << " " << NBenhanceTest << std::endl;
			break;
		}
		NBenhanceTest = NB;
		if(!fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY, 4, passName)) return false;
		//for(auto shift : _cShiftVals.back()) std::cout << shift << " ";
		//std::cout << std::endl;
		if(TMath::Abs(NCenhanceTest-NC)/NC < 0.01) {
			std::cout << "c results stablised " << NC << " " << NCenhanceTest << std::endl;
			break;
		}
		NCenhanceTest = NC;

		if(!fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY, 0, passName)) return false;
		if(TMath::Abs(NCtest-NC)/NC < 0.01 && TMath::Abs(NBtest-NB)/NB < 0.01) {
			std::cout << "results stablised " << NC << " " << NCtest << " " <<  NB << " " << NBtest << std::endl;
			break;
		}
		NCtest = NC;
		NBtest = NB;
	}

	TFile* outFile = TFile::Open(gSaveDir+TString::Format("/enhancedTemplates_%d_%d.root",binPT,binY),"RECREATE");

	TH2D hC("hC","",_nmbins,_mmin,_mmax,_ntbins,2.,2.+_ntbins);
	TH2D hB("hB","",_nmbins,_mmin,_mmax,_ntbins,2.,2.+_ntbins);
	for(uint nbin=0; nbin<_ntbins; ++nbin) {
		TH1 *hC1(0), *hB1(0);
		double fyC(-1.), fyB(-1.);
		TString nameC = TString::Format("pdfC_%d",nbin+2);
		TString nameB = TString::Format("pdfB_%d",nbin+2);

		for(auto it = _pdfCache.rbegin(); it!=_pdfCache.rend(); ++it) {
			auto pdfMap = *it;
			if(pdfMap) {
				if(!hC1 && pdfMap->count(nameC.Data())) {
					hC1 = pdfMap->at(nameC.Data());
				}
				if(!hB1 && pdfMap->count(nameB.Data())) {
					hB1 = pdfMap->at(nameB.Data());
				}
				if(hC1 && hB1) break;
			}
		}
		for(auto it = _fracYieldsCache.rbegin(); it!=_fracYieldsCache.rend(); ++it) {
			auto yieldMap = *it;
			if(yieldMap) {
				if(fyC<0 && yieldMap->count(nameC.Data())) {
					fyC = yieldMap->at(nameC.Data());
				}
				if(fyB<0 && yieldMap->count(nameB.Data())) {
					fyB = yieldMap->at(nameB.Data());
				}
				if(fyC>=0 && fyB>=0) break;
			}
		}
		double scaleC = fyC/hC1->Integral();
		double scaleB = fyB/hB1->Integral();

		for(uint mbin=0; mbin<_nmbins; ++mbin) {
			hC.SetBinContent(mbin+1,nbin+1,scaleC*hC1->GetBinContent(mbin+1));
			hB.SetBinContent(mbin+1,nbin+1,scaleB*hB1->GetBinContent(mbin+1));
		}
	}
	hC.Write();
	hB.Write();
	outFile->Close();

	_finalRun=true;
	return fitSVSim(NB, eB, NC, eC, NQ, eQ, binPT, binY);
}

//function to fit features for a single sample
bool SVFitter::fitSVComb(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY){
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

	TString plotLabel = TString::Format(
				"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
				minPT,
				maxPT);
	if(minY>2. || maxY<6.) plotLabel+=TString::Format(
						"\n#it{y}(#it{Z}) #in [%.2f,%.2f]",
						minY,
						maxY);
	TString plotName = "SVComb_"+_name+"_"; plotName+=ptStr;
	plotFit(SVComb, 0, nbin, nbin, &dh, data_pdf, sig_pdfs, bkg_pdfs, plotName, "M_{cor}, N_{trk}",plotLabel);

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
bool SVFitter::fitSVSim(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY, uint enhance, TString nameMod){
	if(!_inputsSet) {
		std::cout << "ERROR in SVFitter::makeSVFitHists: input files not set yet" << std::endl;
		return false;
	}
	
	double minPT(0.);
	double maxPT(200000.);
	if(_ptBins!=0) {
		if(binPT==0) {
			minPT = _ptBins->GetXaxis()->GetXmin();
			maxPT = _ptBins->GetXaxis()->GetXmax();
		} else {
			minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
			maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
		}
	}
	double minY(0.), maxY(10.);
	if(_yBins!=0) {
		if(binY==0) {
			minY = _yBins->GetXaxis()->GetXmin();
			maxY = _yBins->GetXaxis()->GetXmax();
		} else {
			minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
			maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
		}
	}
	bool allowBShifts(false), allowCShifts(false);

	std::cout << "INFO : fitting SV - binned 2D corrected mass and Ntrk fit for " << _name << ", pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << "), (enhance=" << enhance << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;
	ptStr+="_"; ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* fh0 = TFile::Open(templateFileName("0"));
	TFile* fh4 = TFile::Open(templateFileName("4"));
	TFile* fh5 = TFile::Open(templateFileName("5"));
	TFile* fhd(0);

	if(enhance==4 || enhance==5) fhd = TFile::Open(templateFileName("DE"));
	else fhd = TFile::Open(templateFileName("D"));

	TFile* fhb = TFile::Open(templateFileName("DB"));

	TTree* t0 = dynamic_cast<TTree*>(fh0->Get("T"));
	TTree* t4 = dynamic_cast<TTree*>(fh4->Get("T"));
	TTree* t5 = dynamic_cast<TTree*>(fh5->Get("T"));
	TTree* td = dynamic_cast<TTree*>(fhd->Get("T"));
	TTree* tb(0);
	if(fhb) tb = dynamic_cast<TTree*>(fhb->Get("T"));

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
	TString cutStrTemplates = TString::Format("JetPT>=%f && JetPT<%f && ZY>=%f && ZY<%f && SVMCor>=%f && SVMCor<%f", 20e3, 30e3, 0., 10., _mmin, _mmax);//TODO
	TString cutStrMisID     = TString::Format("JetPT>=%f && JetPT<%f && ZY>=%f && ZY<%f && SVMCor>=%f && SVMCor<%f", 20e3, 30e3, 0., 10., _mmin, _mmax);//TODO
	if(_usePtBinnedTemplates) {
		cutStrTemplates = cutStrData;
		cutStrMisID     = TString::Format("JetPT>=%f && JetPT<%f && ZY>=%f && ZY<%f && SVMCor>=%f && SVMCor<%f", minPT, maxPT, 0., 10., _mmin, _mmax);//TODO
	}

	if(enhance==4) {
		cutStrData += " && enhanced==4";
		ptStr+="enhanceC";
		allowCShifts=true;
		//_cShiftVals.push_back(std::vector<double>());
	} else if(enhance==5) {
		cutStrData += " && enhanced==5";
		ptStr+="enhanceB";
		allowBShifts=true;
		//_bShiftVals.push_back(std::vector<double>());
	}
	ptStr+=nameMod;
	//TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && SVMCor>=%f && SVMCor<%f", minPT, maxPT, _mmin, _mmax);
	//std::cout << cutStr << std::endl;
	td->Draw("SVN>>hSVND",cutStrData);

	//if we have backwards-tagged data and the yield is not being floated then fix the light background yield
	double nLightBack(-1);
	if(tb && !_lightYieldFloat) {
		nLightBack = tb->GetEntries(cutStrData) * _lightYieldScale;
		std::cout << "Fixing light yield to " << nLightBack << std::endl;
	}

	//std::cout << hSVN0.GetEntries() << " " << hSVN4.GetEntries() << " " << hSVN5.GetEntries() << " " << hSVND.GetEntries() << std::endl;

	// -- variables from datasets
	RooRealVar JetPT("JetPT",   "JetPT",   1e4, 1e5,  "MeV/#it{c}"); 
	RooRealVar ZY   ("ZY",      "ZY",      0., 10.,   ""); 
	RooRealVar  SVM ("SVMCor",  "SVMCor",  _mmin, _mmax,  "MeV/#it{c}^{2}"); 
	SVM.setBins(_nmbins);
	//RooRealVar  SVN ("SVN",     "SVN",     2, _ntbins+2,  ""); 
	//SVN.setBins(_ntbins);
	RooRealVar weight("weight",  "",  0., 2.); 
	RooRealVar enhanced("enhanced",  "",  0., 6.); 

	RooCategory SVN(  "SVN",  "SVN");
	for(int i=2; i<_ntbins+2; ++i) {
		//std::cout << TString::Format("SVN%d",i) << " " << i << std::endl;//TODO
		SVN.defineType(TString::Format("SVN%d",i),i);
	}

	RooDataSet ds("ds", "ds", RooArgList(JetPT,ZY,SVM,SVN,weight,enhanced), RooFit::Import(*td), RooFit::WeightVar(weight), RooFit::Cut(cutStrData));
	RooDataSet d0("d0", "d0", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t0), RooFit::WeightVar(weight), RooFit::Cut(cutStrMisID));
	RooDataSet d4("d4", "d4", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t4), RooFit::WeightVar(weight), RooFit::Cut(cutStrTemplates));
	RooDataSet d5("d5", "d5", RooArgList(JetPT,ZY,SVM,SVN,weight), RooFit::Import(*t5), RooFit::WeightVar(weight), RooFit::Cut(cutStrTemplates));

	//TODO probably won't work with RooCategory
	//save raw PDFs
	//TH2* rawHistData   = static_cast<TH2*>(ds->createHistogram("rawHistData",  *SVM,RooFit::Binning(_nmbins,_mmin,_mmax),RooFit::YVar(*SVN,RooFit::Binning(_ntbins,2.,2.+_ntbins))));
	//TH2* rawHistLight  = static_cast<TH2*>(ds->createHistogram("rawHistLight", *SVM,RooFit::Binning(_nmbins,_mmin,_mmax),RooFit::YVar(*SVN,RooFit::Binning(_ntbins,2.,2.+_ntbins))));
	//TH2* rawHistCharm  = static_cast<TH2*>(ds->createHistogram("rawHistCharm", *SVM,RooFit::Binning(_nmbins,_mmin,_mmax),RooFit::YVar(*SVN,RooFit::Binning(_ntbins,2.,2.+_ntbins))));
	//TH2* rawHistBeauty = static_cast<TH2*>(ds->createHistogram("rawHistBeauty",*SVM,RooFit::Binning(_nmbins,_mmin,_mmax),RooFit::YVar(*SVN,RooFit::Binning(_ntbins,2.,2.+_ntbins))));

	//std::cout << d0.sumEntries() << " " << d4.sumEntries() << " " << d5.sumEntries() << " " << ds.sumEntries() << std::endl;
	
	double Ntot = ds.sumEntries();

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.5*Ntot, -0.02*Ntot, Ntot);
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*Ntot, -0.02*Ntot, Ntot);
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.2*Ntot, -0.02*Ntot, Ntot);

	if(nLightBack>=0) {
		//TODO//yieldQ.setConstant(true);
		yieldQ.setVal(nLightBack);
	}

	RooRealVar SVMmin("SVMmin","",_mmin);
	RooRealVar SVMmax("SVMmax","",_mmax);
	RooFormulaVar* SVMNorm = new RooFormulaVar("SVMNorm","","2*(@0-@1)/(@2-@1) - 1",RooArgList(SVM,SVMmin,SVMmax));

	RooSimultaneous data_pdf("data_pdf","data_pdf",SVN);
	RooRealVar *shiftBPrev(0), *shiftCPrev(0);
	std::vector<RooFormulaVar*> partYieldsQ, partYieldsC, partYieldsB;
	std::vector<RooRealVar*> shiftsC, shiftsB;
	for(int i=0; i<_ntbins; ++i) {
		TString nameC = TString::Format("pdfC_%d",i+2);
		TString nameB = TString::Format("pdfB_%d",i+2);

		//slices of the datasets
		RooDataSet* d0slice = static_cast<RooDataSet*>(d0.reduce(RooFit::Name(TString::Format("d0_%d",i+2)),RooFit::Cut(TString::Format("SVN==%d",i+2))));
		RooDataSet* d4slice = static_cast<RooDataSet*>(d4.reduce(RooFit::Name(TString::Format("d4_%d",i+2)),RooFit::Cut(TString::Format("SVN==%d",i+2))));
		RooDataSet* d5slice = static_cast<RooDataSet*>(d5.reduce(RooFit::Name(TString::Format("d5_%d",i+2)),RooFit::Cut(TString::Format("SVN==%d",i+2))));

		//determine yield fractions
		RooRealVar* fracQ = new RooRealVar(TString::Format("fracQ%d",i+2),TString::Format("fracQ%d",i+2),d0slice->sumEntries()/d0.sumEntries());
		RooRealVar* fracC = new RooRealVar(TString::Format("fracC%d",i+2),TString::Format("fracC%d",i+2),d4slice->sumEntries()/d4.sumEntries());
		RooRealVar* fracB = new RooRealVar(TString::Format("fracB%d",i+2),TString::Format("fracB%d",i+2),d5slice->sumEntries()/d5.sumEntries());

		//grab and update the fractional yields
		bool foundFYC(false), foundFYB(false);
		for(auto it = _fracYieldsCache.rbegin(); it!=_fracYieldsCache.rend(); ++it) {
			auto fracYieldsMap = *it;
			if(fracYieldsMap) {
				if(!foundFYC && fracYieldsMap->count(nameC.Data())) {
					fracC->setVal(fracYieldsMap->at(nameC.Data()));
					foundFYC=true;
				}
				if(!foundFYB && fracYieldsMap->count(nameB.Data())) {
					fracB->setVal(fracYieldsMap->at(nameB.Data()));
					foundFYB=true;
				}
				if(foundFYC && foundFYB) break;
			}
		}

		RooAbsPdf *pdfQ(0), *pdfC(0), *pdfB(0), *rawpdfC(0), *rawpdfB(0), *smearC(0), *smearB(0);
		if(_useBinnedTemplates) {
			RooDataHist *d4sliceHist(0), *d5sliceHist(0);
			double smearCVal(10.), smearBVal(10.);

			TString pdfBaseName="pdf";
			if(_smearMCorShapes) {
				pdfBaseName="rawpdf";

				//std::cout << "START" << std::endl;//TODO
				//check whether we've cached smearing factors
				bool foundSC(false), foundSB(false);
				for(auto it = _smearCache.rbegin(); it!=_smearCache.rend(); ++it) {
					auto smearMap = *it;
					//std::cout << smearMap->size() << "!" << std::endl;//TODO
					if(smearMap) {
						if(!foundSC && smearMap->count("smearC")) {
							smearCVal = smearMap->at("smearC");
							//std::cout << "C" << smearCVal << std::endl;//TODO
							foundSC=true;
						}
						if(!foundSB && smearMap->count("smearB")) {
							smearBVal = smearMap->at("smearB");
							//std::cout << "B" << smearBVal << std::endl;//TODO
							foundSB=true;
						}
						if(foundSC && foundSB) break;
					}
				}
				//std::cout << "END" << std::endl;//TODO
			} else{
				//check whether we've cached any corrected PDFs
				for(auto it = _pdfCache.rbegin(); it!=_pdfCache.rend(); ++it) {
					auto pdfMap = *it;
					if(pdfMap) {
						if(!d4sliceHist && pdfMap->count(nameC.Data())) {
							TH1* h = pdfMap->at(nameC.Data());
							d4sliceHist = new RooDataHist(TString::Format("d4Hist_%d",i+2),"",RooArgList(SVM),h);
							//std::cout << "Loaded saved c template" << std::endl;
						}
						if(!d5sliceHist && pdfMap->count(nameB.Data())) {
							TH1* h = pdfMap->at(nameB.Data());
							d5sliceHist = new RooDataHist(TString::Format("d5Hist_%d",i+2),"",RooArgList(SVM),h);
							//std::cout << "Loaded saved b template" << std::endl;
						}
						if(d4sliceHist && d5sliceHist) break;
					}
				}
			}

			//create binned datasets
			if(!d4sliceHist) d4sliceHist = new RooDataHist(TString::Format("d4Hist_%d",i+2),"",RooArgList(SVM),*d4slice);
			if(!d5sliceHist) d5sliceHist = new RooDataHist(TString::Format("d5Hist_%d",i+2),"",RooArgList(SVM),*d5slice);

			//create histogram PDFs
			rawpdfC = new RooHistPdf(pdfBaseName+TString::Format("C_%d",i+2),"",SVM,*d4sliceHist,2);
			rawpdfB = new RooHistPdf(pdfBaseName+TString::Format("B_%d",i+2),"",SVM,*d5sliceHist,2);
			//pdfC = new RooHistPdf(TString::Format("pdfC_%d",i+2),"",SVM,*d4sliceHist,2);
			//pdfB = new RooHistPdf(TString::Format("pdfB_%d",i+2),"",SVM,*d5sliceHist,2);

			if(_smearMCorShapes) {
				RooRealVar* meanSmear = new RooRealVar("meanSmear","meanSmear",0.);
				RooRealVar* widthSmearC = new RooRealVar("widthSmearC","smearC",smearCVal, 1., 5000.);
				if(!allowCShifts) widthSmearC->setConstant();
				RooRealVar* widthSmearB = new RooRealVar("widthSmearB","smearB",smearBVal, 1., 5000.);
				if(!allowBShifts) widthSmearB->setConstant();
				smearC = new RooGaussian("smearC","", SVM, *meanSmear, *widthSmearC);
				smearB = new RooGaussian("smearC","", SVM, *meanSmear, *widthSmearB);
				SVM.setBins(100000, "cache");
				pdfC = new RooFFTConvPdf(TString::Format("pdfC_%d",i+2),"", SVM, *rawpdfC, *smearC);
				pdfB = new RooFFTConvPdf(TString::Format("pdfB_%d",i+2),"", SVM, *rawpdfB, *smearB);
			} else {
				pdfC = rawpdfC;
				pdfB = rawpdfB;
			}
		
		} else {
			//create KDE PDFs
			pdfC = new RooKeysPdf(TString::Format("pdfC_%d",i+2),"",SVM,*d4slice);
			pdfB = new RooKeysPdf(TString::Format("pdfB_%d",i+2),"",SVM,*d5slice);
		}
		//always use KDE for light shape
		pdfQ = new RooKeysPdf(TString::Format("pdfQ_%d",i+2),"",SVM,*d0slice);

		RooRealVar* shiftC(0);
		if(allowCShifts && i<_ntbins-1) {
			shiftC = new RooRealVar(TString::Format("shiftC%d_%d",i+2,i+3),"",0.,-1.,1.);
			shiftsC.push_back(shiftC);
		//} else if(_cShiftVals.size()>0 && _cShiftVals.back().size() > i) {
		//	shiftC = new RooRealVar(TString::Format("shiftC%d_%d",i+2,i+3),"",_cShiftVals.back()[i]);
		//	shiftsC.push_back(shiftC);
		}
		RooRealVar* shiftB(0);
		if(allowBShifts && i<_ntbins-1) {
			shiftB = new RooRealVar(TString::Format("shiftB%d_%d",i+2,i+3),"",0.,-1.,1.);
			shiftsB.push_back(shiftB);
		//} else if(_bShiftVals.size()>0 && _bShiftVals.back().size() > i) {
		//	shiftB = new RooRealVar(TString::Format("shiftB%d_%d",i+2,i+3),"",_bShiftVals.back()[i]);
		//	shiftsB.push_back(shiftB);
		}

		RooFormulaVar* partYieldQ = new RooFormulaVar(TString::Format("partYieldQ_%d",i+2),"","@0*@1",RooArgList(yieldQ,*fracQ));
		RooFormulaVar* partYieldC(0);
		if(shiftCPrev && shiftC) {
			partYieldC = new RooFormulaVar(TString::Format("partYieldC_%d",i+2),"","@0*(@1+@2-@3)",RooArgList(yieldC,*fracC,*shiftCPrev,*shiftC));
		} else if(shiftC) {
			partYieldC = new RooFormulaVar(TString::Format("partYieldC_%d",i+2),"","@0*(@1-@2)",RooArgList(yieldC,*fracC,*shiftC));
		} else if(shiftCPrev) {
			partYieldC = new RooFormulaVar(TString::Format("partYieldC_%d",i+2),"","@0*(@1+@2)",RooArgList(yieldC,*fracC,*shiftCPrev));
		} else {
			partYieldC = new RooFormulaVar(TString::Format("partYieldC_%d",i+2),"","@0*@1",RooArgList(yieldC,*fracC));
		}
		RooFormulaVar* partYieldB(0);
		if(shiftBPrev && shiftB) {
			partYieldB = new RooFormulaVar(TString::Format("partYieldB_%d",i+2),"","@0*(@1+@2-@3)",RooArgList(yieldB,*fracB,*shiftBPrev,*shiftB));
		} else if(shiftB) {
			partYieldB = new RooFormulaVar(TString::Format("partYieldB_%d",i+2),"","@0*(@1-@2)",RooArgList(yieldB,*fracB,*shiftB));
		} else if(shiftBPrev) {
			partYieldB = new RooFormulaVar(TString::Format("partYieldB_%d",i+2),"","@0*(@1+@2)",RooArgList(yieldB,*fracB,*shiftBPrev));
		} else {
			partYieldB = new RooFormulaVar(TString::Format("partYieldB_%d",i+2),"","@0*@1",RooArgList(yieldB,*fracB));
		}
		partYieldsQ.push_back(partYieldQ);
		partYieldsC.push_back(partYieldC);
		partYieldsB.push_back(partYieldB);
		shiftCPrev = shiftC;
		shiftBPrev = shiftB;

		//create total PDF
		RooAddPdf* pdf = new RooAddPdf( TString::Format("data_pdf_%d",i+2),  "", RooArgList(*pdfB, *pdfC, *pdfQ), RooArgList(*partYieldB, *partYieldC, *partYieldQ) );
		data_pdf.addPdf(*pdf, TString::Format("SVN%d",i+2));
		//std::cout << data_pdf.getPdf(TString::Format("SVN%d",i+2))->GetName() << std::endl;//TODO
	}
	//if we don't have any versions of the PDFs cached then cache uncorrected versions
	if(_pdfCache.empty()) cachePDFs(&SVM, &ds, &data_pdf, "pdfB*,pdfC*");

	RooProdPdf* data_pdf_constr(0);
	if(nLightBack>=0) {
		RooRealVar* yieldQmean = new RooRealVar("yieldQmean","",nLightBack);
		RooRealVar* yieldQsigma = new RooRealVar("yieldQsigma","",TMath::Sqrt(nLightBack));
		RooGaussian* light_constr = new RooGaussian("light_constr","",yieldQ,*yieldQmean,*yieldQsigma);
		data_pdf_constr = new RooProdPdf("data_pdf_constr","",RooArgList(data_pdf,*light_constr));
	}

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/log/SV2D_"+_name+"_"+ptStr+"_fits.log","w");
	RooFitResult* r(0);
	int attempt(0), nattempts(1);
	while(true) {
		std::cout << "ATTEMPT " << attempt << std::endl;
		std::cout << "STARTING NQ=" << yieldQ.getVal() << ", NC=" << yieldC.getVal() << ", NB=" << yieldB.getVal() << std::endl;
		if(r) delete r;
		if(nLightBack>=0) {
			r = data_pdf_constr->fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::Constrain(RooArgSet(yieldQ)));//TODO
		} else {
			r = data_pdf.fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
		}
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

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ*" );
	sig_pdfs.push_back( "pdfB*" );
	sig_pdfs.push_back( "pdfC*" );
	std::vector<std::string> bkg_pdfs;

	plotFit(SVM, _mmin, _mmax, _nmbins, &ds, data_pdf, sig_pdfs, bkg_pdfs, "svfig/SVM_"+_name+"_"+ptStr, "M_{cor} [MeV/#it{c}^{2}]","");
	//plotFit(SVN, 2, 2+_ntbins, _ntbins, &ds, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}","");
	makeMassPlots(&ds, &data_pdf, _name+"_"+ptStr, binPT, binY, 2, 2+_ntbins);

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "dat/params_SV2D_"+_name+"_"; paramsName+=ptStr; paramsName+=".dat";
	printParams(paramsName,params);
	TString allParamsName = "dat/paramsAll_SV2D_"+_name+"_"; allParamsName+=ptStr; allParamsName+=".dat";
	printAllParams(allParamsName,&ds,&data_pdf);

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV();
	NC = yieldC.getValV();
	NQ = yieldQ.getValV();
	eB = yieldB.getError();
	eC = yieldC.getError();
	eQ = yieldQ.getError();
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ);
	
	double maxCorrC(0.), maxCorrB(0.);
	if(allowCShifts) {
		//for(auto shiftC : shiftsC) {
		//	_cShiftVals.back().push_back(shiftC->getVal());
		//}
		maxCorrC = calcCorrection(&SVM, &ds, &data_pdf, "pdfC*", "corrC"+_name+"_"+ptStr);
	}
	if(allowBShifts) {
		//for(auto shiftB : shiftsB) {
		//	_bShiftVals.back().push_back(shiftB->getVal());
		//}
		maxCorrB = calcCorrection(&SVM, &ds, &data_pdf, "pdfB*", "corrB"+_name+"_"+ptStr);
	}
	//std::cout << maxCorrC << " " << maxCorrB << std::endl;

	for(int ibin=0; ibin<_ntbins; ++ibin) {
		hSVN0.SetBinContent(ibin+1, partYieldsQ[ibin]->getValV());
		hSVN4.SetBinContent(ibin+1, partYieldsC[ibin]->getValV());
		hSVN5.SetBinContent(ibin+1, partYieldsB[ibin]->getValV());
		hSVN0.SetBinError(ibin+1, partYieldsQ[ibin]->getPropagatedError(*r));
		hSVN4.SetBinError(ibin+1, partYieldsC[ibin]->getPropagatedError(*r));
		hSVN5.SetBinError(ibin+1, partYieldsB[ibin]->getPropagatedError(*r));
	}

	hSVNF.Add(&hSVN0);
	hSVNF.Add(&hSVN4);
	hSVNF.Add(&hSVN5);

	hSVN0.SetLineColor(kBlue);
	hSVN4.SetLineColor(kGreen+2);
	hSVN5.SetLineColor(kRed);
	hSVNF.SetLineColor(kBlack);
	hSVND.SetMarkerStyle(kFullCircle);
	hSVND.SetLineColor(kBlack);
	hSVN0.SetFillStyle(3245);
	hSVN0.SetFillColor(kBlue);
	hSVN4.SetFillStyle(3644);
	hSVN4.SetFillColor(kGreen+2);
	hSVN5.SetFillStyle(3454);
	hSVN5.SetFillColor(kRed);

	hSVNF.SetMinimum(0.);
	hSVNF.SetMaximum(1.1*TMath::Max(hSVNF.GetMaximum(),hSVND.GetMaximum()));
	hSVNF.GetXaxis()->SetTitle("N_{trk}");
	hSVNF.GetYaxis()->SetTitle("Candidates");
	hSVNF.SetTitle("");
	hSVNF.GetXaxis()->SetLabelOffset(0.02);
	hSVNF.GetXaxis()->SetTitleOffset(1.18);
	hSVNF.GetYaxis()->SetTitleOffset(1.1);

	hSVNF.GetXaxis()->SetBinLabel(1,"2");
	hSVNF.GetXaxis()->SetBinLabel(2,"3");
	hSVNF.GetXaxis()->SetBinLabel(3,"4");//4+

	TCanvas c;
	c.SetLeftMargin(0.16);
	c.SetBottomMargin(0.18);

	hSVNF.Draw("HIST");
	hSVN5.Draw("HIST SAME");
	hSVN4.Draw("HIST SAME");
	hSVN0.Draw("HIST SAME");
	hSVND.Draw("SAME E1 X0 P");

	TPaveText label(0.65,0.7,0.9,0.9,"BRNDC");
	if(binPT>0) {
		label.AddText(TString::Format(
			"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
			minPT/1e3,
			maxPT/1e3));
	}
	if(binY>0) {
		label.AddText(TString::Format(
			"#it{y}(#it{Z}) #in [%.2f,%.2f]",
			minY,
			maxY));
	}
	label.SetFillColor(0);
	label.SetTextAlign(12);
	label.SetBorderSize(0);
	label.Draw();

	TString name = _name+"_"+ptStr;
	name.ReplaceAll(".","");
	c.SaveAs(gSaveDir+"/svfig/SVN_"+name+".pdf");

	if(_runToyFits && _finalRun) {
		//do toy fits
		RooMCStudy mcs(data_pdf, RooArgSet(SVM,SVN), RooFit::Extended(), RooFit::FitOptions(RooFit::PrintEvalErrors(-1), RooFit::PrintLevel(2), RooFit::Warnings(kFALSE)),RooFit::Silence());
		mcs.generateAndFit(_nToys, ds.sumEntries(),true,gSaveDir+"/svtoys/toy_"+_name+ptStr+"_%04d.dat");
		//TODO//mcs.fit(nToys, gSaveDir+"/svtoys/toy_"+_name+ptStr+"_%04d.dat");
		mcs.fitParDataSet().write(gSaveDir+"/svtoys/allToys_"+_name+ptStr+".dat");
		std::cout << "Toy parameters:" << std::endl;
		mcs.fitParDataSet().printMultiline(std::cout,0);

		TCanvas c1;
		RooPlot* cYieldPlot = mcs.plotPull(yieldC,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-3,3));
		cYieldPlot->SetTitle("");
		cYieldPlot->SetXTitle("Pull(#it{N}_{#it{c}})");
		cYieldPlot->Draw();
		c1.SaveAs(gSaveDir+"/svfig/pull_charmYield"+_name+ptStr+".pdf");
		RooPlot* bYieldPlot = mcs.plotPull(yieldB,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-3,3));
		bYieldPlot->SetTitle("");
		bYieldPlot->SetXTitle("Pull(#it{N}_{#it{b}})");
		bYieldPlot->Draw();
		c1.SaveAs(gSaveDir+"/svfig/pull_beautyYield_"+_name+ptStr+".pdf");

	}

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

	plotFit(SVM, _mmin, _mmax, _nmbins, &dh, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+_name+"_"+ptStr, "M_{cor}","");
	plotFit(SVN, 2, 2+_ntbins, _ntbins, &dh, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}","");

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

	plotFit(SVM, _mmin, _mmax, _nmbins, &data, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+_name+"_"+ptStr, "M_{cor}","");
	plotFit(SVN, 2, 2+_ntbins, _ntbins, &data, data_pdf, sig_pdfs, bkg_pdfs, "SVN_"+_name+"_"+ptStr, "N_{trk}","");

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
	//check whether we have a separate enhanced data sample for training
	//data-driven correction. If not, use the main data sample
	if(which=="DE" && _trainingInputFile=="") which="D";
	return gSaveDir+ "/svFitTemplates_"+_name+which+".root";
}

double SVFitter::calcCorrection(RooRealVar* var, RooAbsData* ds, RooAbsPdf* pdf, TString toCorrect, TString plotName) {
	std::cout << "INFO in SVFitter::calcCorrection: calculating correction for " << toCorrect << std::endl;
	TH1 *hData(0), *hFit(0), *hComp(0), *hErr(0), *hCorr(0), *hNew(0);

	//plotName+="pass";
	//plotName+=_pdfCache.size();

	if(!ds || !pdf) return 0;
	RooArgList fullCompList;

	pdf->branchNodeServerList(&fullCompList);
	//for(int i=0; i<fullCompList.getSize(); ++i) std::cout << fullCompList.at(i)->GetName() << std::endl;
	RooArgList* fitList = static_cast<RooArgList*>(fullCompList.selectByName("pdf*"));
	RooArgList* compList = static_cast<RooArgList*>(fullCompList.selectByName(toCorrect));
	if(fitList->getSize()<1) return 0;
	if(compList->getSize()<1) return 0;

	hData = ds ->createHistogram("hData",*var);

	for(int i=0; i<fitList->getSize(); ++i) {
		RooAbsPdf* comp = static_cast<RooAbsPdf*>(fitList->at(i));
		TString name = comp->GetName();
		name.ReplaceAll("pdf","partYield");
		RooAbsReal* yieldPar = static_cast<RooRealVar*>(fullCompList.find(name));
		if(!yieldPar) {
			std::cout << "Yield " << name << " not found!" << std::endl;
			continue;
		}

		if(hFit) {
			TH1* hTemp = comp->createHistogram("hTemp",*var);
			hTemp->Scale(yieldPar->getValV()/hTemp->Integral());
			hFit->Add(hTemp);
			delete hTemp;
		} else {
			hFit = comp->createHistogram("hFit",*var);
			hFit->Scale(yieldPar->getValV()/hFit->Integral());
		}
	}

	for(int i=0; i<compList->getSize(); ++i) {
		RooAbsPdf* comp = static_cast<RooAbsPdf*>(compList->at(i));
		TString name = comp->GetName();
		name.ReplaceAll("pdf","partYield");
		RooAbsReal* yieldPar = static_cast<RooRealVar*>(fullCompList.find(name));
		if(!yieldPar) {
			std::cout << "Yield " << name << " not found!" << std::endl;
			continue;
		}

		if(hComp) {
			TH1* hTemp = comp->createHistogram("hTemp",*var);
			hTemp->Scale(yieldPar->getValV()/hTemp->Integral());
			hComp->Add(hTemp);
			delete hTemp;
		} else {
			hComp = comp->createHistogram("hComp",*var);
			hComp->Scale(yieldPar->getValV()/hComp->Integral());
		}
	}

	hErr = static_cast<TH1*>(hData->Clone("hErr"));
	hErr->Add(hFit,-1.);

	hCorr = static_cast<TH1*>(hErr->Clone("hCorr"));

	if(_additiveMCorCorrectionFactor) {
		hCorr->Scale(1./hComp->Integral());
	} else {
		hCorr->Add(hComp);
		hCorr->Divide(hComp);
	}

	double min(-1.), max(-1.);
	bool foundStart(false);
	for(int ibin=1; ibin<=hCorr->GetNbinsX(); ++ibin) {
		//hCorr->SetBinError(ibin,1./TMath::Sqrt(hComp->GetBinContent(ibin)));
		if(!_additiveMCorCorrectionFactor && (hComp->GetBinContent(ibin)<0.05*hComp->GetMaximum() || hComp->GetBinContent(ibin)<30.)) {
			hCorr->SetBinContent(ibin,1.);
		} else {
			if(!foundStart) {
				min = hCorr->GetBinLowEdge(ibin);
				foundStart = true;
			} else {
				max = hCorr->GetBinLowEdge(ibin+1);
			}
		}
	}

	double retVal(0.);

	TF1* func(0);
	if(foundStart) {
		if(_additiveMCorCorrectionFactor) {
			func = new TF1("func","pol7",min,max);
		} else {
			func = new TF1("func","pol5",min,max);
		}
		hCorr->Fit(func,"Q");

		retVal = TMath::Max(func->GetMaximum()-1.,1.-func->GetMinimum());
	}

	//cache the corrected PDFs
	cachePDFs(var, ds, pdf, toCorrect, func);

	if(foundStart) {
		//correct total histogram for plotting
		hNew = static_cast<TH1*>(hComp->Clone("hNew"));
		for(int ibin=1; ibin<=hNew->GetNbinsX(); ++ibin) {
			double x = hNew->GetBinCenter(ibin);
			if(x<min) x=min;
			if(x>max) x=max;
			double corr = func->Eval(x);
			if(_additiveMCorCorrectionFactor) {
				hNew->SetBinContent(ibin, hNew->GetBinContent(ibin) + corr*hComp->Integral());
			} else {
				hNew->SetBinContent(ibin, hNew->GetBinContent(ibin) * corr);
			}
		}
	}

	//make plots
	TCanvas c;
	hData->SetMaximum(1.1*TMath::Max(hData->GetMaximum(),hFit->GetMaximum()));
	hData->Draw("EP");
	hFit->SetLineColor(kRed);
	hFit->Draw("hist same");
	hComp->SetLineColor(kGreen+2);
	hComp->Draw("hist same");
	if(foundStart) {
		hNew->SetLineColor(kMagenta);
		hNew->Draw("hist same");
	}
	c.SaveAs(gSaveDir+"/svfig/"+plotName+"Components.pdf");

	hErr->GetYaxis()->SetTitle("Residual / (100 MeV/#it{c}^{2})");
	hErr->Draw("hist");
	c.SaveAs(gSaveDir+"/svfig/"+plotName+"Difference.pdf");

	if(foundStart) {
		if(!_additiveMCorCorrectionFactor) {
			hCorr->SetMinimum(0.);
			hCorr->SetMaximum(2.);
			hCorr->GetYaxis()->SetTitle("Correction factor");
		} else {
			hCorr->GetYaxis()->SetTitle("Residual");
		}
		func->SetLineColor(kRed);
		hCorr->Draw("HIST");//"EP");
		func->Draw("same");
		c.SaveAs(gSaveDir+"/svfig/"+plotName+"CorrectionFit.pdf");
	}

	delete hData;
	delete hFit;
	delete hComp;
	delete hErr;
	delete hCorr;
	if(foundStart) {
		delete hNew;
		delete func;
	}

	return retVal;
}

void SVFitter::cachePDFs(RooRealVar* var, RooAbsData* ds, RooAbsPdf* pdf, TString toCorrect, TF1* corrFunc) {
	//std::cout << "caching" << std::endl;//TODO
	std::map<std::string,TH1*>* pdfMap(0);
	std::map<std::string,double>* yieldMap(0);
	std::map<std::string,double>* smearMap(0);
	if(_pdfCache.size()>0) {
		pdfMap = _pdfCache.back();
		yieldMap = _fracYieldsCache.back();
	} else {
		pdfMap = new std::map<std::string,TH1*>();
		_pdfCache.push_back(pdfMap);
		yieldMap = new std::map<std::string,double>();
		_fracYieldsCache.push_back(yieldMap);
	}

	//now make new histograms for all requested PDFs
	RooArgList fullCompList;
	pdf->branchNodeServerList(&fullCompList);
	RooArgList* compList = static_cast<RooArgList*>(fullCompList.selectByName(toCorrect));
	RooArgList* parList  = static_cast<RooArgList*>(pdf->getParameters(*ds)->selectByName("*"));

	for(int icomp=0; icomp<compList->getSize(); ++icomp) {
		RooAbsPdf* comp = static_cast<RooAbsPdf*>(compList->at(icomp));
		TString name = comp->GetName();

		//if no PDF map or any PDFs are already in the back map then make a new map
		if(pdfMap->count(name.Data())) {
			pdfMap = new std::map<std::string,TH1*>();
			_pdfCache.push_back(pdfMap);
			yieldMap = new std::map<std::string,double>();
			_fracYieldsCache.push_back(yieldMap);
		}

		//get yield fraction
		TString yname = name;
		yname.ReplaceAll("pdf","partYield");
		RooAbsReal* yieldPar = static_cast<RooRealVar*>(fullCompList.find(yname));
		RooAbsReal* denomPar(0);
		if(yname.BeginsWith("partYieldB")) denomPar = static_cast<RooRealVar*>(parList->find("yieldB"));
		if(yname.BeginsWith("partYieldC")) denomPar = static_cast<RooRealVar*>(parList->find("yieldC"));
		if(!yieldPar || !denomPar) {
			std::cout << "Yield " << yname << " not found!" << std::endl;
			return;
		}

		TString hName = name;
		hName += _pdfCache.size();

		TH1* h = comp->createHistogram(hName,*var);
		double intFact = h->Integral();
		for(int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
			double x = h->GetBinCenter(ibin);
			if(corrFunc && x<corrFunc->GetXmin()) x=corrFunc->GetXmin();
			if(corrFunc && x>corrFunc->GetXmax()) x=corrFunc->GetXmax();
			double corr(1.);
			if(_additiveMCorCorrectionFactor) {
				corr = 0.;
			}
			if(corrFunc) corr = corrFunc->Eval(x);
			if(_additiveMCorCorrectionFactor) {
				h->SetBinContent(ibin, h->GetBinContent(ibin) + corr*intFact);
			}
			else {
				h->SetBinContent(ibin, h->GetBinContent(ibin) * corr);
			}
		}
		pdfMap->insert(std::pair<TString,TH1*>(name.Data(),h));
		yieldMap->insert(std::pair<TString,double>(name.Data(),yieldPar->getValV()/denomPar->getValV()));
	}

	//get smear factors
	smearMap = new std::map<std::string,double>();
	_smearCache.push_back(smearMap);
	RooRealVar* smear = static_cast<RooRealVar*>(parList->find("widthSmearC"));
	if(smear) {
		smearMap->insert(std::pair<TString,double>("smearC",smear->getValV()));
		//std::cout << "caching smearC=" << smear->getValV() << std::endl;//TODO
	}
	smear = static_cast<RooRealVar*>(parList->find("widthSmearB"));
	if(smear) {
		smearMap->insert(std::pair<TString,double>("smearB",smear->getValV()));
		//std::cout << "caching smearB=" << smear->getValV() << std::endl;//TODO
	}
}

void SVFitter::makeMassPlots(RooAbsData* ds, RooAbsPdf* pdf, TString name, int binPT, int binY, int minNBin, int maxNBin) {
	RooArgList fullCompList;
	pdf->branchNodeServerList(&fullCompList);
	RooArgList* compList = static_cast<RooArgList*>(fullCompList.selectByName("*"));
	RooArgList* parList  = static_cast<RooArgList*>(pdf->getParameters(*ds)->selectByName("*"));
	RooArgList* obsList  = static_cast<RooArgList*>(pdf->getObservables(*ds)->selectByName("*"));

	//for(int icomp=0; icomp<compList->getSize(); ++icomp) {
	//	RooAbsPdf* comp = static_cast<RooAbsPdf*>(compList->at(icomp));
	//	TString name = comp->GetName();
	//	std::cout << "HERE" << name << std::endl;//TODO
	//}

	//for(int ipar=0; ipar<parList->getSize(); ++ipar) {
	//	RooAbsPdf* par = static_cast<RooAbsPdf*>(parList->at(ipar));
	//	TString name = par->GetName();
	//	std::cout << "PARHERE" << name << std::endl;//TODO
	//}
	//for(int iobs=0; iobs<obsList->getSize(); ++iobs) {
	//	RooAbsPdf* obs = static_cast<RooAbsPdf*>(obsList->at(iobs));
	//	TString name = obs->GetName();
	//	std::cout << "OBSHERE" << name << std::endl;//TODO
	//}
	name.ReplaceAll(".","");
	RooRealVar* SVMCor = static_cast<RooRealVar*>(obsList->find("SVMCor"));

	uint nBins1D = _nmbins;
	uint nBins2D = 30;
	uint nBins = _nmbins*2;

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

	//Ntrk-integrated histograms
	TH1* dataHist   = static_cast<TH1*>(ds->createHistogram("dataHist",*SVMCor,RooFit::Binning(nBins, _mmin, _mmax)));
	TH1* modelHist  = new TH1D("modelHist", "",nBins, _mmin, _mmax);
	TH1* charmHist = new TH1D("charmHist","",nBins, _mmin, _mmax);
	TH1* beautyHist = new TH1D("beautyHist","",nBins, _mmin, _mmax);
	TH1* bkgrndHist = new TH1D("bkgrndHist","",nBins, _mmin, _mmax);

	double* mbins = new double[nBins+2];
	double* tbins = new double[_ntbins+3];
	for(int i=0; i<=nBins; ++i) {
		mbins[i+1] = _mmin + i*_mmax/nBins;
	}
	mbins[0] = _mmin-0.01;
	for(int i=0; i<=_ntbins; ++i) {
		tbins[i+1] = i+2;
	}
	tbins[0] = 1.99;
	tbins[_ntbins+2] = 2.+_ntbins+0.01;

	//extra 2 bins in Ntrks for plotting purposes (narrow border bins to make contours plot nicely)
	TH2* dataHist2D   = new TH2D("dataHist2D",  "",nBins, _mmin, _mmax, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* modelHist2D  = new TH2D("modelHist2D", "",nBins+1, mbins, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* charmHist2D  = new TH2D("charmHist2D", "",nBins+1, mbins, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* beautyHist2D = new TH2D("beautyHist2D","",nBins+1, mbins, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* bkgrndHist2D = new TH2D("bkgrndHist2D","",nBins+1, mbins, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* pullHist2D   = new TH2D("pullHist2D","",  10, _mmin, _mmax, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* pullHist2D_D = new TH2D("pullHist2D_D","",  10, _mmin, _mmax, _ntbins+2, tbins);//2., 2.+_ntbins);
	TH2* pullHist2D_M = new TH2D("pullHist2D_M","",  10, _mmin, _mmax, _ntbins+2, tbins);//2., 2.+_ntbins);

	//axis labels
	for(auto hist : {dataHist,modelHist,charmHist,beautyHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}_{cor} [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("");
	}

	//axis labels
	for(auto hist : {dataHist2D,modelHist2D,charmHist2D,beautyHist2D,bkgrndHist2D}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}_{cor} [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("#it{N}_{trk}");
		hist->GetZaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	charmHist->SetLineColor(kGreen+2);
	beautyHist->SetLineColor(kRed);
	bkgrndHist->SetLineColor(kBlue);

	dataHist2D->SetLineColor(kBlack);
	dataHist2D->SetMarkerStyle(20);
	modelHist2D->SetLineColor(kBlack);
	charmHist2D->SetLineColor(kGreen+2);
	beautyHist2D->SetLineColor(kRed);
	bkgrndHist2D->SetLineColor(kBlue);
	charmHist2D->SetContour(5);
	beautyHist2D->SetContour(5);
	bkgrndHist2D->SetContour(5);

	//store integrated and Ntrk slice histograms
	std::vector<TH1*> dataHists {dataHist};
	std::vector<TH1*> modelHists {modelHist};
	std::vector<TH1*> charmHists {charmHist};
	std::vector<TH1*> beautyHists {beautyHist};
	std::vector<TH1*> bkgrndHists {bkgrndHist};

	//load all Ntrk slices
	for(int ntrk=minNBin; ntrk<maxNBin; ++ntrk) {
		RooAbsPdf* charm  = static_cast<RooAbsPdf*>(compList->find(TString::Format("pdfC_%d",ntrk)));
		RooAbsPdf* beauty = static_cast<RooAbsPdf*>(compList->find(TString::Format("pdfB_%d",ntrk)));
		RooAbsPdf* bkgrnd = static_cast<RooAbsPdf*>(compList->find(TString::Format("pdfQ_%d",ntrk)));

		//if no parameters are split then only a single PDF will exist
		if(!charm)  charm  = static_cast<RooAbsPdf*>(compList->find("pdfC"));
		if(!beauty) beauty = static_cast<RooAbsPdf*>(compList->find("pdfB"));
		if(!bkgrnd) bkgrnd = static_cast<RooAbsPdf*>(compList->find("pdfQ"));

		TH1* charmSlice = static_cast<TH1*>(charm->createHistogram(TString::Format("charmSlice%d",ntrk),*SVMCor,RooFit::Binning(nBins,_mmin,_mmax)));
		TH1* beautySlice = static_cast<TH1*>(beauty->createHistogram(TString::Format("beautySlice%d",ntrk),*SVMCor,RooFit::Binning(nBins,_mmin,_mmax)));
		TH1* bkgrndSlice = static_cast<TH1*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ntrk),*SVMCor,RooFit::Binning(nBins,_mmin,_mmax)));
		TH1* dataSlice   = static_cast<TH1*>(ds->createHistogram(    TString::Format("dataSlice%d",ntrk),  *SVMCor,RooFit::Binning(nBins,_mmin,_mmax),
					                                                             RooFit::Cut(TString::Format("SVN==%d",ntrk))));

		TH1* modelSlice  = new TH1D(TString::Format("modelSlice%d",ntrk), "",nBins,_mmin,_mmax);

		charmSlice->Scale( static_cast<RooRealVar*>(parList->find(TString::Format("fracC%d",ntrk)))->getVal() * static_cast<RooRealVar*>(parList->find("yieldC"))->getVal());
		beautySlice->Scale(static_cast<RooRealVar*>(parList->find(TString::Format("fracB%d",ntrk)))->getVal() * static_cast<RooRealVar*>(parList->find("yieldB"))->getVal());
		bkgrndSlice->Scale(static_cast<RooRealVar*>(parList->find(TString::Format("fracQ%d",ntrk)))->getVal() * static_cast<RooRealVar*>(parList->find("yieldQ"))->getVal());

		//make composites
		modelSlice->Add(charmSlice);
		modelSlice->Add(beautySlice);
		modelSlice->Add(bkgrndSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		charmHist->Add(charmSlice);
		beautyHist->Add(beautySlice);
		bkgrndHist->Add(bkgrndSlice);

		//add to 2D histograms
		int jbin = ntrk;
		for(int ibin=1; ibin<nBins; ++ibin) {
			dataHist2D  ->SetBinContent(ibin,jbin, dataSlice  ->GetBinContent(ibin));
			modelHist2D ->SetBinContent(ibin+1,jbin, modelSlice ->GetBinContent(ibin));
			charmHist2D ->SetBinContent(ibin+1,jbin, charmSlice ->GetBinContent(ibin));
			beautyHist2D->SetBinContent(ibin+1,jbin, beautySlice->GetBinContent(ibin));
			bkgrndHist2D->SetBinContent(ibin+1,jbin, bkgrndSlice->GetBinContent(ibin));
		}

		//setup style
		for(auto hist : {dataSlice,modelSlice,charmSlice,beautySlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("#it{m}_{cor} [MeV/#it{c}^{2}]");
			hist->GetYaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		charmSlice->SetLineColor(kGreen+2);
		beautySlice->SetLineColor(kRed);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		charmHists.push_back(charmSlice);
		beautyHists.push_back(beautySlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	double x, y, n;
	//fill adaptive binning histograms
	for(uint ibin=1; ibin<=dataHist2D->GetNbinsX(); ++ibin) {
		for(uint jbin=1; jbin<=dataHist2D->GetNbinsY(); ++jbin) {
			x = dataHist2D->GetXaxis()->GetBinCenter(ibin);
			y = dataHist2D->GetYaxis()->GetBinCenter(jbin);
			n = dataHist2D->GetBinContent( ibin, jbin );
			//std::cout << "D " << ibin << " " << jbin << " " << x << " " << y << " " << n << std::endl;//TODO
			pullHist2D_D->Fill(x, y, n);
			//model has extra underflow bin to make contours plot
			x = modelHist2D->GetXaxis()->GetBinCenter(ibin+1);
			y = modelHist2D->GetYaxis()->GetBinCenter(jbin);
			n = modelHist2D->GetBinContent( ibin+1, jbin );
			//std::cout << "M " << ibin << " " << jbin << " " << x << " " << y << " " << n << std::endl;//TODO
			pullHist2D_M->Fill(x, y, n);
		}
	}
	for(uint ibin=1; ibin<=pullHist2D_D->GetNbinsX(); ++ibin) {
		for(uint jbin=1; jbin<=pullHist2D_D->GetNbinsY(); ++jbin) {
			double D = pullHist2D_D->GetBinContent(ibin,jbin);
			double F = pullHist2D_M->GetBinContent(ibin,jbin);
			double E = TMath::Sqrt(D);
			//std::cout << ibin << " " << jbin << " " << D << " " << F << " " << E << std::endl;//TODO
			pullHist2D->SetBinContent(ibin,jbin,(D-F)/E);
		}
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH1* dataHist50 = dataHists[i]->Rebin( nBins/nBins1D, "dataHist50");

		//Mass plots
		TH1* dataHistMass   = dataHist50;
		TH1* charmHistMass = charmHists[i];
		TH1* beautyHistMass = beautyHists[i];
		TH1* bkgrndHistMass = bkgrndHists[i];
		TH1* modelHistMass  = modelHists[i];

		charmHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		beautyHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());

		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (_mmax-_mmin)/nBins1D));
		bkgrndHistMass->SetFillStyle(3245);
		charmHistMass->SetFillColor(kGreen+2);
		charmHistMass->SetFillStyle(3644);
		beautyHistMass->SetFillColor(kRed);
		beautyHistMass->SetFillStyle(3454);

		dataHistMass->Draw("E1 X0 P");
		if(charmHistMass->Integral()>0) charmHistMass->Draw("HIST C SAME");
		if(beautyHistMass->Integral()>0) beautyHistMass->Draw("HIST C SAME");
		if(bkgrndHistMass->Integral()>0) bkgrndHistMass->Draw("HIST C SAME");
		modelHistMass->Draw( "HIST C SAME");
		dataHistMass->Draw("E1 X0 P SAME");

		TPaveText label(0.65,0.7,0.9,0.9,"BRNDC");
		//if(i==0) {
			if(binPT>0) {
				double minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
				double maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
				label.AddText(TString::Format(
					"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
					minPT/1e3,
					maxPT/1e3));
			}
			if(binY>0) {
				double minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
				double maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
				label.AddText(TString::Format(
					"#it{y}(#it{Z}) #in [%.2f,%.2f]",
					minY,
					maxY));
			}
		//} else {
		if(i!=0) {
			label.AddText(TString::Format(
				"#it{N}_{trk} = %.0f",
				static_cast<double>(1+i)));
		}
		label.SetFillColor(0);
		label.SetTextAlign(12);
		label.SetBorderSize(0);
		label.Draw();
		if(i==0) can.SaveAs(gSaveDir+"/svfig/fitSVM"+name+".pdf");
		else can.SaveAs(gSaveDir+"/svfig/fitSVM"+name+"NTrkBin"+(i-1)+".pdf");
	}
	
	//2D plots
	TH2* dataHist2D30 = dataHist2D->Rebin2D( nBins/nBins2D, 1,"dataHist2D30");
	dataHist2D30->GetYaxis()->SetBinLabel(2,"2");
	dataHist2D30->GetYaxis()->SetBinLabel(3,"3");
	dataHist2D30->GetYaxis()->SetBinLabel(4,"4");//4+
	can.SetTickx();
	can.SetTicky();
	dataHist2D30->Draw("col");
	charmHist2D->Draw("cont2 same");
	beautyHist2D->Draw("cont2 same");
	bkgrndHist2D->Draw("cont2 same");

	TPaveText label(0.65,0.7,0.9,0.9,"BRNDC");
	if(binPT>0) {
		double minPT = _ptBins->GetXaxis()->GetBinLowEdge(binPT);
		double maxPT = _ptBins->GetXaxis()->GetBinUpEdge(binPT);
		label.AddText(TString::Format(
			"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
			minPT/1e3,
			maxPT/1e3));
	}
	if(binY>0) {
		double minY = _yBins->GetXaxis()->GetBinLowEdge(binY);
		double maxY = _yBins->GetXaxis()->GetBinUpEdge(binY);
		label.AddText(TString::Format(
			"#it{y}(#it{Z}) #in [%.2f,%.2f]",
			minY,
			maxY));
	}
	label.SetFillColor(0);
	label.SetTextAlign(12);
	label.SetBorderSize(0);
	label.Draw();
	can.SaveAs(gSaveDir+"/svfig/fitSV2D"+name+".pdf");

	can.SetRightMargin(0.10);
	const Int_t NRGBs2 = 4;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.45, 0.55, 1.00};
	Double_t reds2[NRGBs2]   = { 0.00, 1.00, 1.00, 1.00};
	Double_t greens2[NRGBs2] = { 0.00, 1.00, 1.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 1.00, 1.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);
	pullHist2D->SetMinimum(-4.);
	pullHist2D->SetMaximum( 4.);
	pullHist2D->Draw("colz");
	can.SaveAs(gSaveDir+"/svfig/pull"+name+".pdf");

	delete dataHist;
	delete modelHist;
	delete charmHist;
	delete beautyHist;
	delete bkgrndHist;
	delete dataHist2D;
	delete modelHist2D;
	delete charmHist2D;
	delete beautyHist2D;
	delete bkgrndHist2D;
	delete pullHist2D;
	delete pullHist2D_D;
	delete pullHist2D_M;
}

//void SVFitter::make2DPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
//	name.ReplaceAll(".","");
//	RooRealVar* DM = ws->var("D0M");
//	RooRealVar* DIP = ws->var("D0LogIPChi2");
//
//	uint nBins1D(50);
//	uint nBins2D(30);
//	uint nBins(150);
//
//	double mmin(1784.), mmax(1944.), ipmin(-5.), ipmax(15.);
//
//	//plotting style
//	setLHCbStyle();
//	gStyle->SetOptStat(0);
//	gStyle->SetPalette(1,0);
//	const Int_t NRGBs = 5;
//	const Int_t NCont = 255;
//	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
//	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
//	Double_t greens[NRGBs] = { 1.00, 0.95, 0.50, 0.00, 0.00};
//	Double_t blues[NRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
//	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
//	gStyle->SetNumberContours(NCont);
//
//	//pT-integrated histograms
//	TH2* dataHist   = static_cast<TH2*>(ds->createHistogram("dataHist",*DM,RooFit::Binning(nBins,mmin,mmax),RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
//	TH2* modelHist  = new TH2D("modelHist", "",nBins,mmin,mmax,nBins,ipmin,ipmax);
//	TH2* signalHist = new TH2D("signalHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
//	TH2* promptHist = new TH2D("promptHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
//	TH2* displcHist = new TH2D("displcHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
//	TH2* bkgrndHist = new TH2D("bkgrndHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
//
//	//axis labels
//	for(auto hist : {dataHist,modelHist,signalHist,promptHist,displcHist,bkgrndHist}) {
//		hist->SetTitle("");
//		hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
//		hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
//		hist->GetZaxis()->SetTitle("");
//	}
//
//	dataHist->SetLineColor(kBlack);
//	dataHist->SetMarkerStyle(20);
//	modelHist->SetLineColor(kBlack);
//	signalHist->SetLineColor(kRed);
//	promptHist->SetLineColor(kGreen+2);
//	displcHist->SetLineColor(kMagenta);
//	bkgrndHist->SetLineColor(kBlue);
//	promptHist->SetContour(5);
//	displcHist->SetContour(5);
//	bkgrndHist->SetContour(5);
//
//	//store integrated and pT slice histograms
//	std::vector<TH2*> dataHists {dataHist};
//	std::vector<TH2*> modelHists {modelHist};
//	std::vector<TH2*> signalHists {signalHist};
//	std::vector<TH2*> promptHists {promptHist};
//	std::vector<TH2*> displcHists {displcHist};
//	std::vector<TH2*> bkgrndHists {bkgrndHist};
//
//	//load all pT slices
//	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
//		RooAbsPdf* prompt = ws->pdf(TString("prompt_PTbin")+ibin);
//		RooAbsPdf* displc = ws->pdf(TString("displc_PTbin")+ibin);
//		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrnd_PTbin")+ibin);
//
//		//if no parameters are split then only a single PDF will exist
//		if(!prompt) prompt = ws->pdf(TString("prompt"));
//		if(!displc) displc = ws->pdf(TString("displc"));
//		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrnd"));
//
//		TH2* promptSlice = static_cast<TH2*>(prompt->createHistogram(TString::Format("promptSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
//					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
//		TH2* displcSlice = static_cast<TH2*>(displc->createHistogram(TString::Format("displcSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
//					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
//		TH2* bkgrndSlice = static_cast<TH2*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
//					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
//		TH2* dataSlice   = static_cast<TH2*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DM,RooFit::Binning(nBins,mmin,mmax),
//					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax)),
//					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));
//
//		TH2* modelSlice  = new TH2D(TString::Format("modelSlice%d",ibin), "",nBins,mmin,mmax,nBins,ipmin,ipmax);
//		TH2* signalSlice = new TH2D(TString::Format("signalSlice%d",ibin),"",nBins,mmin,mmax,nBins,ipmin,ipmax);
//
//		promptSlice->Scale(ws->var(TString::Format("promptYield_PTbin%d",ibin))->getVal());
//		displcSlice->Scale(ws->var(TString::Format("displcYield_PTbin%d",ibin))->getVal());
//		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());
//
//		//make composites
//		modelSlice->Add(promptSlice);
//		modelSlice->Add(displcSlice);
//		modelSlice->Add(bkgrndSlice);
//		signalSlice->Add(promptSlice);
//		signalSlice->Add(displcSlice);
//
//		//add to integrated histograms
//		modelHist->Add(modelSlice);
//		signalHist->Add(signalSlice);
//		promptHist->Add(promptSlice);
//		displcHist->Add(displcSlice);
//		bkgrndHist->Add(bkgrndSlice);
//
//		//setup style
//		for(auto hist : {dataSlice,modelSlice,signalSlice,promptSlice,displcSlice,bkgrndSlice}) {
//			hist->SetTitle("");
//			hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
//			hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
//			hist->GetZaxis()->SetTitle("");
//		}
//
//		dataSlice->SetLineColor(kBlack);
//		dataSlice->SetMarkerStyle(20);
//		modelSlice->SetLineColor(kBlack);
//		signalSlice->SetLineColor(kRed);
//		promptSlice->SetLineColor(kGreen+2);
//		displcSlice->SetLineColor(kMagenta);
//		bkgrndSlice->SetLineColor(kBlue);
//
//		dataHists.push_back(dataSlice);
//		modelHists.push_back(modelSlice);
//		signalHists.push_back(signalSlice);
//		promptHists.push_back(promptSlice);
//		displcHists.push_back(displcSlice);
//		bkgrndHists.push_back(bkgrndSlice);
//	}
//
//	//make pull histograms
//	TH2* dataHist10 = dataHist->Rebin2D(nBins/10,nBins/10,"dataHist10");
//	TH2* modelHist10 = modelHist->Rebin2D(nBins/10,nBins/10,"modelHist10");
//
//	//make pull histograms with adaptive binning
//	int nBinAdapt = TMath::Nint(dataHist10->GetSumOfWeights()/30.); //target 30 entries per bins
//	if(nBinAdapt<16) nBinAdapt=16;//set a reasonable minimum
//	AdaptBin binning("pullHist",nBinAdapt,mmin,mmax,ipmin,ipmax);
//	binning.setVerbose();//TODO
//	binning.loadDataFromHist(dataHist);
//	TH2Poly* dataPoly = binning.getHisto("dataPoly");
//	TH2Poly* modelPoly = binning.getHisto("modelPoly");
//	modelPoly->SetTitle("");
//	modelPoly->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
//	modelPoly->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
//	modelPoly->GetZaxis()->SetTitle("");
//
//	//create pull histogram
//	modelHist10->Scale(dataHist10->Integral()/modelHist10->Integral());
//	for(int ibin=1; ibin<=10; ++ibin) {
//		for(int jbin=1; jbin<=10; ++jbin) {
//			double D = dataHist10->GetBinContent(ibin,jbin);
//			double F = modelHist10->GetBinContent(ibin,jbin);
//			double E = dataHist10->GetBinError(ibin,jbin);
//			modelHist10->SetBinContent(ibin,jbin,(D-F)/E);
//		}
//	}
//
//	double x, y, n;
//	//fill adaptive binning histograms
//	for(uint ibin=1; ibin<=nBins; ++ibin) {
//		for(uint jbin=1; jbin<=nBins; ++jbin) {
//			x = dataHist->GetXaxis()->GetBinCenter(ibin);
//			y = dataHist->GetYaxis()->GetBinCenter(jbin);
//			n = dataHist->GetBinContent( ibin, jbin );
//			dataPoly->Fill(x, y, n);
//			x = modelHist->GetXaxis()->GetBinCenter(ibin);
//			y = modelHist->GetYaxis()->GetBinCenter(jbin);
//			n = modelHist->GetBinContent( ibin, jbin );
//			modelPoly->Fill(x, y, n);
//		}
//	}
//	for(int ibin=1; ibin<=dataPoly->GetNumberOfBins(); ++ibin) {
//		double D = dataPoly->GetBinContent(ibin);
//		double F = modelPoly->GetBinContent(ibin);
//		double E = TMath::Sqrt(D);
//		modelPoly->SetBinContent(ibin,(D-F)/E);
//	}
//
//	TCanvas can;
//
//	for(uint i=0; i<dataHists.size(); ++i) {
//		TH2* dataHist50 = dataHists[i]->Rebin2D( nBins/nBins1D, nBins/nBins1D,"dataHist50");
//
//		//Mass plots
//		TH1* dataHistMass   = dataHist50->ProjectionX("dataHistMass");
//		TH1* signalHistMass = signalHists[i]->ProjectionX("signalHistMass");
//		TH1* promptHistMass = promptHists[i]->ProjectionX("promptHistMass");
//		TH1* displcHistMass = displcHists[i]->ProjectionX("displcHistMass");
//		TH1* bkgrndHistMass = bkgrndHists[i]->ProjectionX("bkgrndHistMass");
//		TH1* modelHistMass  = modelHists[i]->ProjectionX("modelHistMass");
//
//		signalHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
//		promptHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
//		displcHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
//		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
//		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
//
//		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (mmax-mmin)/nBins1D));
//		bkgrndHistMass->SetLineStyle(kDashed);
//		signalHistMass->SetLineStyle(kDotted);
//		promptHistMass->SetFillColor(kGreen+2);
//		displcHistMass->SetFillColor(kMagenta);
//		promptHistMass->SetFillStyle(3345);
//		displcHistMass->SetFillStyle(3354);
//
//		dataHistMass->Draw("E1 X0 P");
//		signalHistMass->Draw("HIST C SAME");
//		promptHistMass->Draw("HIST C SAME");
//		displcHistMass->Draw("HIST C SAME");
//		bkgrndHistMass->Draw("HIST C SAME");
//		modelHistMass->Draw( "HIST C SAME");
//		dataHistMass->Draw("E1 X0 P SAME");
//
//		TPaveText label(0.65,0.8,0.9,0.9,"BRNDC");
//		if(i==0) {
//			label.AddText(TString::Format(
//				"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
//				_jptmin/1e3,
//				_jptmax/1e3));
//		} else {
//			if(_usePtFracBins) {
//				label.AddText(TString::Format(
//					"#it{f}_{#it{p}_{T}}^{#it{D}} #in [%.2f,%.2f]",
//					static_cast<double>(i-1)/_nDPtBins,
//					static_cast<double>(i)/_nDPtBins));
//			} else {
//				label.AddText(TString::Format(
//					"#it{p}_{T}(#it{D}) #in [%.0f,%.0f] GeV/#it{c}",
//					_ptBinsD[i-1],
//					_ptBinsD[i]));
//			}
//		}
//		label.SetFillColor(0);
//		label.SetTextAlign(12);
//		label.SetBorderSize(0);
//		label.Draw();
//		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+".pdf");
//		else can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+"PTbin"+(i-1)+".pdf");
//
//		//IP plots
//		TH1* dataHistIP   = dataHist50->ProjectionY("dataHistIP");
//		TH1* signalHistIP = signalHists[i]->ProjectionY("signalHistIP");
//		TH1* promptHistIP = promptHists[i]->ProjectionY("promptHistIP");
//		TH1* displcHistIP = displcHists[i]->ProjectionY("displcHistIP");
//		TH1* bkgrndHistIP = bkgrndHists[i]->ProjectionY("bkgrndHistIP");
//		TH1* modelHistIP  = modelHists[i]->ProjectionY("modelHistIP");
//
//		signalHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
//		promptHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
//		displcHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
//		bkgrndHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
//		modelHistIP ->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
//
//		dataHistIP->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f)", (ipmax-ipmin)/nBins1D));
//		bkgrndHistIP->SetLineStyle(kDashed);
//		signalHistIP->SetLineStyle(kDotted);
//		promptHistIP->SetFillColor(kGreen+2);
//		displcHistIP->SetFillColor(kMagenta);
//		promptHistIP->SetFillStyle(3345);
//		displcHistIP->SetFillStyle(3354);
//
//		dataHistIP->Draw("E1 X0 P");
//		signalHistIP->Draw("HIST C SAME");
//		promptHistIP->Draw("HIST C SAME");
//		displcHistIP->Draw("HIST C SAME");
//		bkgrndHistIP->Draw("HIST C SAME");
//		modelHistIP->Draw( "HIST C SAME");
//		dataHistIP->Draw("E1 X0 P SAME");
//		label.Draw();
//		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+".pdf");
//		else can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+"PTbin"+(i-1)+".pdf");
//	}
//	
//	//2D plots
//	TH2* dataHist30 = dataHist->Rebin2D( nBins/nBins2D, nBins/nBins2D,"dataHist30");
//	can.SetTickx();
//	can.SetTicky();
//	dataHist30->Draw("col");
//	//signalHist->Draw("cont2 same");
//	promptHist->Draw("cont2 same");
//	displcHist->Draw("cont2 same");
//	bkgrndHist->Draw("cont2 same");
//
//	TPaveText label(0.65,0.8,0.9,0.9,"BRNDC");
//	label.AddText(TString::Format(
//		"#it{p}_{T}(#it{j}) #in [%.0f,%.0f] GeV/#it{c}",
//		_jptmin/1e3,
//		_jptmax/1e3));
//	label.SetFillColor(0);
//	label.SetTextAlign(12);
//	label.SetBorderSize(0);
//	label.Draw();
//	can.SaveAs(gSaveDir+"/fig/data"+name+".pdf");
//	promptHist->Draw("col");
//	can.SaveAs(gSaveDir+"/fig/promptPDF"+name+".pdf");
//	displcHist->Draw("col");
//	can.SaveAs(gSaveDir+"/fig/displcPDF"+name+".pdf");
//	bkgrndHist->Draw("col");
//	can.SaveAs(gSaveDir+"/fig/bkgrndPDF"+name+".pdf");
//
//	can.SetRightMargin(0.10);
//	const Int_t NRGBs2 = 4;
//	const Int_t NCont2 = 255;
//	Double_t stops2[NRGBs2]  = { 0.00, 0.45, 0.55, 1.00};
//	Double_t reds2[NRGBs2]   = { 0.00, 1.00, 1.00, 1.00};
//	Double_t greens2[NRGBs2] = { 0.00, 1.00, 1.00, 0.00};
//	Double_t blues2[NRGBs2]  = { 1.00, 1.00, 1.00, 0.00};
//	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
//	gStyle->SetNumberContours(NCont2);
//	modelHist10->SetMinimum(-4.);
//	modelHist10->SetMaximum( 4.);
//	modelHist10->Draw("colz");
//	can.SaveAs(gSaveDir+"/fig/pull"+name+".pdf");
//	modelPoly->SetMinimum(-4.);
//	modelPoly->SetMaximum( 4.);
//	modelPoly->Draw("colz");
//	can.SaveAs(gSaveDir+"/fig/pull_adapbin"+name+".pdf");
//}
