#include <vector>
#include <fstream>

#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooMsgService.h"

#include "outputFunctions.h"
#include "DFitter.h"
#include "SVFitter.h"
#include "MCJets.h"

//globals to save passing these around
double d0minpt(5000);
double d0maxpt(-1);
TString savedir("output");
bool useSimpleEff(false);
bool useRhoZEffCor(true);
bool weightRecPt(false);

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
bool oldSim = false;
TString charmSimFile  = "/data/dijets/for_yandex_data_new_14X.root";
TString beautySimFile = "/data/dijets/for_yandex_data_new_15X.root";
//TString lightSimFile  = "/data/dijets/for_yandex_data_new_101.root";
//TString dataFile      = "/data/dijets/for_yandex_data_new_100.root";
//bool oldSim = false;
//TString charmSimFile  = "/data/dijets/dijets_sim4.root";
//TString beautySimFile = "/data/dijets/dijets_sim5.root";
TString lightSimFile  = "/data/dijets/dijets_back_2016.root";
TString dataFile      = "/data/dijets/dijets_2016.root";

//efficiency inputs - update these with latest version numbers
TString simpleEffFile = "../efficiencies/SimpleEffs_2XX_16x8bins_up190216.root";
TString effFile = "../efficiencies/D0Effs_2XX_16x8bins_up190213.root";
TString accFile = "../efficiencies/D0AccEffNewUp190205.root";

TString lightHistFile = "svFitHists0.root";
TString charmHistFile = "svFitHists4.root";
TString beautyHistFile = "svFitHists5.root";
TString dataHistFile = "svFitHistsD.root";

bool dataIsMC(false);
bool dataIsResampledMC(false);

double getPtCorrFactor(MCJets::jetType type, double ptMin, double ptMax) {
	TFile* f(0);
	TTree* t(0);
	TString cut1;
	TString cut2;
	double corr(1.);

	cut1 = "JetPT>"; cut1+=ptMin; cut1+=" && JetPT<"; cut1+=ptMax;
	cut2 = "JetPT>"; cut2+=ptMin; cut2+=" && JetPT<"; cut2+=ptMax;

	switch(type) {
		case MCJets::jetRecoD04:
			std::cout << "getting pT correction factor for c->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return corr;//no correction needed
		case MCJets::jetRecoSV4:
			std::cout << "getting pT correction factor for c->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(charmSimFile);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000. && SVM[0]";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && SVM[0]";
			break;
		case MCJets::jetRecoD05:
			std::cout << "getting pT correction factor for b->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(beautySimFile);
			//cut1+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEBPT>5000.";
			//cut2+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEDPT>5000.";
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case MCJets::jetRecoSV5:
			std::cout << "getting pT correction factor for b->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(beautySimFile);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000. && SVM[0]";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && SVM[0]";
			break;
		default:
			//TODO not covered
			return 1.;
	}

	if(!f) return 1.;
	t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return 1.;

	corr = t->GetEntries(cut1)/static_cast<double>(t->GetEntries(cut2));
	std::cout << "Correction factor: " << corr << std::endl;
	f->Close();
	return corr;
}

void makeTrainTestSamples(int sample) {
	std::cout << "INFO : making samples for MC study " << sample << std::endl;
	TRandom3 rand(1000+sample);
	TString name = "_"; name+=sample;
	TFile* fin4 = TFile::Open(charmSimFile);
	TTree* tin4 = dynamic_cast<TTree*>(fin4->Get("T"));

	TFile* fout40 = TFile::Open("/tmp/dcraik/for_yandex_data_new_14X_train"+name+".root","RECREATE");
	TTree* tout40 = tin4->CloneTree(0);

	TFile* fout41 = TFile::Open("/tmp/dcraik/for_yandex_data_new_14X_test"+name+".root","RECREATE");
	TTree* tout41 = tin4->CloneTree(0);

	TFile* fin5 = TFile::Open(beautySimFile);
	TTree* tin5 = dynamic_cast<TTree*>(fin5->Get("T"));

	TFile* fout50 = TFile::Open("/tmp/dcraik/for_yandex_data_new_15X_train"+name+".root","RECREATE");
	TTree* tout50 = tin5->CloneTree(0);

	TFile* fout51 = TFile::Open("/tmp/dcraik/for_yandex_data_new_15X_test"+name+".root","RECREATE");
	TTree* tout51 = tin5->CloneTree(0);

	boost::progress_display progress( tin4->GetEntries()+tin5->GetEntries() );
	for(int ientry=0; ientry<tin4->GetEntries(); ++ientry) {
		++progress;
		tin4->GetEntry(ientry);
		if(rand.Rndm()<0.5) tout40->Fill();
		else tout41->Fill();
	}

	for(int ientry=0; ientry<tin5->GetEntries(); ++ientry) {
		++progress;
		tin5->GetEntry(ientry);
		if(rand.Rndm()<0.5) tout50->Fill();
		else tout51->Fill();
	}

	tout40->AutoSave();
	fout40->Close();
	tout41->AutoSave();
	fout41->Close();
	tout50->AutoSave();
	fout50->Close();
	tout51->AutoSave();
	fout51->Close();

	//now merge the test samples
	TChain* tout45 = new TChain("T");
	tout45->Add("/tmp/dcraik/for_yandex_data_new_14X_test"+name+".root");
	tout45->Add("/tmp/dcraik/for_yandex_data_new_15X_test"+name+".root");
	tout45->Merge("/tmp/dcraik/for_yandex_data_new"+name+".root");

	//update the standard locations and set the dataIsMC flag
	charmSimFile = "/tmp/dcraik/for_yandex_data_new_14X_train"+name+".root";
	beautySimFile = "/tmp/dcraik/for_yandex_data_new_15X_train"+name+".root";
	dataFile = "/tmp/dcraik/for_yandex_data_new"+name+".root";
	dataIsMC=true;
	dataIsResampledMC=true;
}

bool getTruth(TH1D* trueD04, TH1D* trueD05, TH1D* trueSV4, TH1D* trueSV5, bool useTruePT=false) {
	std::cout << "INFO : getting truth information" << std::endl;
	if(!trueD04 || !trueD05 || !trueSV4 || !trueSV5) return false;

	TFile* f = TFile::Open(dataFile);
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

	double JetPT;
	double JetTruePT;
	double JetTrueD0;
	double JetTrueDSV;
	double JetTrueBSV;

	std::vector<double>* TRUEDID = new std::vector<double>();
	std::vector<double>* TRUEDFROMB = new std::vector<double>();

	std::vector<double>* SVN = new std::vector<double>();

	t->SetBranchAddress("JetPT",      &JetPT);
	t->SetBranchAddress("JetTruePT",  &JetTruePT);
	t->SetBranchAddress("JetTRUED0",  &JetTrueD0);
	t->SetBranchAddress("JetTRUEDSV", &JetTrueDSV);
	t->SetBranchAddress("JetTRUEBSV", &JetTrueBSV);

	t->SetBranchAddress("TRUEDID",    &TRUEDID);
	t->SetBranchAddress("TRUEDFROMB", &TRUEDFROMB);

	t->SetBranchAddress("SVN",        &SVN);

	boost::progress_display progress( t->GetEntries() );
	for(int i=0; i<t->GetEntries(); ++i) {
		++progress;
		t->GetEntry(i);
		if(useTruePT) JetPT = JetTruePT;

		double weight(1.);
		if(JetTruePT>50000.) {
			weight=0.007;
		} else if(JetTruePT>20000.) {
			weight=0.10;
		} else if(JetTruePT>15000.) {
			weight=0.25;
		}

		if(JetTrueD0) {
			for(unsigned int id=0; id<TRUEDID->size(); ++id) {
				if(TMath::Abs(TRUEDID->at(id))!=421.) continue;
				if(TRUEDFROMB->at(id)) {
					trueD05->Fill(JetPT,weight);
				} else {
					trueD04->Fill(JetPT,weight);
				}
				break;
			}
		}
		if(JetTrueDSV) {
			if(SVN->size()>0) trueSV4->Fill(JetPT,weight);
		}
		if(JetTrueBSV) {
			if(SVN->size()>0) trueSV5->Fill(JetPT,weight);
		}
	}

	return true;
}

int main(int argc, char** argv) {
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	//0-99 are MC closure tests, 100 is dijet tagged jets, 101 is jets with backwards SVs, 102 is J/psi tagged jets
	TString file="100";//"201X";
	if(argc>1) file = argv[1];
	if(argc>2) savedir = argv[2];
	if(argc>3) d0minpt = atoi(argv[3]);//5000
	if(argc>4) d0maxpt = atoi(argv[4]);//-1

	gSystem->Exec("mkdir -p "+savedir);
	gSaveDir = savedir;

	if(file.IsDec() && atoi(file)<100) {
		std::cout << "Running MC closure test " << file << " (" << atoi(file) << ")" << std::endl;
		makeTrainTestSamples(atoi(file));
	} else {
		if(file.BeginsWith("sim")) {// || file.BeginsWith("14") || file.BeginsWith("15") || file.BeginsWith("16") || file.BeginsWith("45")) {
			dataIsMC=true;
		}
//		dataFile = "/data/dijets/dijets_"+file+".root";
//		dataHistFile = "svFitHistsD"+file+".root";
	}

	unsigned int npt=3;
	//double* binsPt  = new double[npt +1]{17500.,20000.,30000.,50000.};
	double* binsPt  = new double[npt +1]{15000.,20000.,30000.,50000.};
	//double* binsPt  = new double[npt +1]{15000.,20000.,30000.,100000.};
	//unsigned int npt=4;
	//double* binsPt  = new double[npt +1]{10000.,15000.,20000.,30000.,100000.};

	//truth histograms for MC studies
	TH1D trueD04("trueD04","",npt,binsPt);
	TH1D trueD05("trueD05","",npt,binsPt);
	TH1D trueSV4("trueSV4","",npt,binsPt);
	TH1D trueSV5("trueSV5","",npt,binsPt);
	TH1D trueUnfldD04("trueUnfldD04","",npt,binsPt);
	TH1D trueUnfldD05("trueUnfldD05","",npt,binsPt);
	TH1D trueUnfldSV4("trueUnfldSV4","",npt,binsPt);
	TH1D trueUnfldSV5("trueUnfldSV5","",npt,binsPt);

	DFitter d0fit("jjDFit"+file, &trueD04);
	if(useSimpleEff) {
		d0fit.setInputs(dataFile,simpleEffFile,"",true,dataIsMC,useRhoZEffCor);
	} else {
		d0fit.setInputs(dataFile,effFile,accFile,false,dataIsMC,useRhoZEffCor);
	}
	d0fit.addEffs();

	//weight MC for continuous true jet pT
	unsigned int nmcpt=4;
	double* ptInputWeightBins  = new double[nmcpt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* jetTruePtWeights4 = new TH1D("jetTruePtWeights4","",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights5 = new TH1D("jetTruePtWeights5","",nmcpt,ptInputWeightBins);
	if(oldSim) {
		jetTruePtWeights4->SetBinContent(1,1.);
		jetTruePtWeights4->SetBinContent(2,0.3);
		jetTruePtWeights4->SetBinContent(3,0.12);
		jetTruePtWeights4->SetBinContent(4,0.008);
		jetTruePtWeights5->SetBinContent(1,1.);
		jetTruePtWeights5->SetBinContent(2,0.2);
		jetTruePtWeights5->SetBinContent(3,0.08);
		jetTruePtWeights5->SetBinContent(4,0.006);
	} else {
		jetTruePtWeights4->SetBinContent(1,1.);
		jetTruePtWeights4->SetBinContent(2,0.6);
		jetTruePtWeights4->SetBinContent(3,0.6);
		jetTruePtWeights4->SetBinContent(4,0.07);
		jetTruePtWeights5->SetBinContent(1,1.);
		jetTruePtWeights5->SetBinContent(2,0.6);
		jetTruePtWeights5->SetBinContent(3,0.5);
		jetTruePtWeights5->SetBinContent(4,0.02);
	}
	//TH2D test("test","",npt,binsPt,1,0.,10.);
	MCJets wmc("jj");
	wmc.setInputs("",charmSimFile,beautySimFile,"",d0fit.dFileName());
	wmc.setInputTruePtWeights(MCJets::jetRecoD04,jetTruePtWeights4);
	wmc.setInputTruePtWeights(MCJets::jetRecoD05,jetTruePtWeights5);
	wmc.setInputTruePtWeights(MCJets::jetRecoSV4,jetTruePtWeights4);
	wmc.setInputTruePtWeights(MCJets::jetRecoSV5,jetTruePtWeights5);
	if(!wmc.weightMC(MCJets::jetRecoD04)) return 1;
	if(!wmc.weightMC(MCJets::jetRecoD05)) return 1;
	if(!wmc.weightMC(MCJets::jetRecoSV4/*,&test*/)) return 1;
	if(!wmc.weightMC(MCJets::jetRecoSV5)) return 1;

	//for(int i=1; i<=npt; ++i) {
	//	std::cout << test.GetBinContent(i,1) << "+/-" << test.GetBinError(i,1) << std::endl;
	//}

	if(dataIsMC) {
		getTruth(&trueD04,&trueD05,&trueSV4,&trueSV5);
		getTruth(&trueUnfldD04,&trueUnfldD05,&trueUnfldSV4,&trueUnfldSV5,true);
		//for the real MC samples, test D0 efficiencies and return
		//for permutations, compare truth results to extracted yields
		if(!dataIsResampledMC) {
			d0fit.testEffs(4);
			d0fit.testEffs(5);
			return 0;//TODO
		}
	}

	SVFitter svfit("jjFit"+file, &trueD04);
	svfit.setInputs(lightSimFile,wmc.outputName(MCJets::jetRecoSV4),wmc.outputName(MCJets::jetRecoSV5),dataFile);
	svfit.setSVBinning(96,500.,10000.,3);
	svfit.makeSVFitHists(0);
	svfit.makeSVFitHists(4);
	svfit.makeSVFitHists(5);
	svfit.makeSVFitHists(7);

	//first do D0 fits for denominators
	TH1D recoD04("recoD04","",npt,binsPt);
	TH1D recoD05("recoD05","",npt,binsPt);

	double yield(0.), error(0.), corr(1.);
	for(int i=1; i<=recoD04.GetNbinsX(); ++i) {
		if(d0fit.fitD(4,yield,error,i)) {
			corr = getPtCorrFactor(MCJets::jetRecoD04,recoD04.GetBinLowEdge(i),recoD04.GetBinLowEdge(i+1));
			recoD04.SetBinContent(i,yield*corr);
			recoD04.SetBinError(  i,error*corr);
		}
		if(d0fit.fitD(5,yield,error,i)) {
			corr = getPtCorrFactor(MCJets::jetRecoD05,recoD05.GetBinLowEdge(i),recoD05.GetBinLowEdge(i+1));
			recoD05.SetBinContent(i,yield*corr);
			recoD05.SetBinError(  i,error*corr);
		}
	}

	//do SV fits for numerator
	TH1D recoSV4("recoSV4","",npt,binsPt);
	TH1D recoSV5("recoSV5","",npt,binsPt);
	for(int i=1; i<=recoSV4.GetNbinsX(); ++i) {
		double nB(0.), eB(0.), nC(0.), eC(0.), nQ(0.), eQ(0.);
		if(!svfit.fitSV(nB,eB,nC,eC,nQ,eQ,i)) continue;
		double corrC = getPtCorrFactor(MCJets::jetRecoSV4,recoSV4.GetBinLowEdge(i),recoSV4.GetBinLowEdge(i+1));
		double corrB = getPtCorrFactor(MCJets::jetRecoSV5,recoSV5.GetBinLowEdge(i),recoSV5.GetBinLowEdge(i+1));
		recoSV4.SetBinContent(i,nC*corrC);
		recoSV4.SetBinError(  i,eC*corrC);
		recoSV5.SetBinContent(i,nB*corrB);
		recoSV5.SetBinError(  i,eB*corrB);
	}

	//do unfolding
	TH1D* unfoldedD04 = wmc.unfold(&recoD04, MCJets::jetRecoD04);
	TH1D* unfoldedD05 = wmc.unfold(&recoD05, MCJets::jetRecoD05);
	TH1D* unfoldedSV4 = wmc.unfold(&recoSV4, MCJets::jetRecoSV4);
	TH1D* unfoldedSV5 = wmc.unfold(&recoSV5, MCJets::jetRecoSV5);

	if(!unfoldedD04 || !unfoldedD05) return 1;
	if(!unfoldedSV4 || !unfoldedSV5) return 1;

	//print results
	double bfD0 = 0.0389;
	double errBFD0 = 0.0004;
	double ffc2D0 = 0.542;
	double errFFc2D0 = TMath::Sqrt(.024*.024 + .007*.007);
	double ffb2D0 = 0.587;
	double errFFb2D0 = TMath::Sqrt(.021*.021 + .008*.008);

	double bfffErr4 = TMath::Sqrt( (errFFc2D0/ffc2D0)*(errFFc2D0/ffc2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );
	double bfffErr5 = TMath::Sqrt( (errFFb2D0/ffb2D0)*(errFFb2D0/ffb2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );

	//D0 results
	std::cout << "D0 reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinError(i)/recoD04.GetBinContent(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinError(i)/recoD05.GetBinContent(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinContent(i) << " +/- " << recoD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinContent(i) << " +/- " << recoD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04->GetBinContent(i) << " +/- " << unfoldedD04->GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05->GetBinContent(i) << " +/- " << unfoldedD05->GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) << " +/- " << trueD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) << " +/- " << trueD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) << " +/- " << trueUnfldD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) << " +/- " << trueUnfldD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	//D0 results scaled
	std::cout << "jet reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << recoD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << recoD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << unfoldedD04->GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << unfoldedD05->GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueUnfldD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueUnfldD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	//SV results
	std::cout << "SV reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV4.GetBinContent(i) << " +/- " << recoSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV5.GetBinContent(i) << " +/- " << recoSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV4->GetBinContent(i) << " +/- " << unfoldedSV4->GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV5->GetBinContent(i) << " +/- " << unfoldedSV5->GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV4.GetBinContent(i) << " +/- " << trueSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV5.GetBinContent(i) << " +/- " << trueSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV4.GetBinContent(i) << " +/- " << trueUnfldSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV5.GetBinContent(i) << " +/- " << trueUnfldSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	std::ofstream fout("toysSVResults.log",std::ofstream::app);
	fout << file << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV4.GetBinContent(i) << "\t" << recoSV4.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV5.GetBinContent(i) << "\t" << recoSV5.GetBinError(i) << "\t";
	fout << std::endl;
	fout.close();

	//ratios
	std::vector<double> ratioRec4;
	std::vector<double> ratioRec5;
	std::vector<double> ratioUnfld4;
	std::vector<double> ratioUnfld5;
	std::vector<double> ratioTrue4;
	std::vector<double> ratioTrue5;
	std::vector<double> ratioTrueUnfld4;
	std::vector<double> ratioTrueUnfld5;
	std::vector<double> errorRec4;
	std::vector<double> errorRec5;
	std::vector<double> errorUnfld4;
	std::vector<double> errorUnfld5;
	for (unsigned int i=1; i<=npt; ++i) {
		ratioRec4.push_back(recoSV4.GetBinContent(i) / (recoD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioRec5.push_back(recoSV5.GetBinContent(i) / (recoD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioUnfld4.push_back(unfoldedSV4->GetBinContent(i) / (unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioUnfld5.push_back(unfoldedSV5->GetBinContent(i) / (unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioTrue4.push_back(trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrue5.push_back(trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioTrueUnfld4.push_back(trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrueUnfld5.push_back(trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		errorRec4.push_back(ratioRec4[i-1]*TMath::Sqrt(TMath::Power(recoSV4.GetBinError(i)/recoSV4.GetBinContent(i),2)+TMath::Power(recoD04.GetBinError(i)/recoD04.GetBinContent(i),2)));
		errorRec5.push_back(ratioRec5[i-1]*TMath::Sqrt(TMath::Power(recoSV5.GetBinError(i)/recoSV5.GetBinContent(i),2)+TMath::Power(recoD05.GetBinError(i)/recoD05.GetBinContent(i),2)));
		errorUnfld4.push_back(ratioUnfld4[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV4->GetBinError(i)/unfoldedSV4->GetBinContent(i),2)+TMath::Power(unfoldedD04->GetBinError(i)/unfoldedD04->GetBinContent(i),2)));
		errorUnfld5.push_back(ratioUnfld5[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV5->GetBinError(i)/unfoldedSV5->GetBinContent(i),2)+TMath::Power(unfoldedD05->GetBinError(i)/unfoldedD05->GetBinContent(i),2)));
	}
	std::cout << "ratios reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioRec4[i] << "+/-" << errorRec4[i] << "+/-" << ratioRec4[i]*bfffErr4 << "\t";
	std::cout << std::endl;                                                              
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioRec5[i] << "+/-" << errorRec5[i] << "+/-" << ratioRec5[i]*bfffErr5 << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioUnfld4[i] << "+/-" << errorUnfld4[i] << "+/-" << ratioUnfld4[i]*bfffErr4 << "\t";
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioUnfld5[i] << "+/-" << errorUnfld5[i] << "+/-" << ratioUnfld5[i]*bfffErr5 << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	//plots
	recoD04.SetMinimum(0.);
	recoD04.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04.GetMaximum(),recoD05.GetMaximum()),TMath::Max(unfoldedD04->GetMaximum(),unfoldedD05->GetMaximum())));
	recoD04.SetLineColor(kBlue);
	recoD05.SetLineColor(kRed);
	unfoldedD04->SetLineColor(kBlue);
	unfoldedD05->SetLineColor(kRed);
	unfoldedD04->SetLineStyle(kDashed);
	unfoldedD05->SetLineStyle(kDashed);
	trueD04.SetLineColor(kBlue);
	trueD05.SetLineColor(kRed);
	trueD04.SetLineStyle(kDotted);
	trueD05.SetLineStyle(kDotted);

	recoSV4.SetMinimum(0.);
	recoSV4.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04.GetMaximum(),recoD05.GetMaximum()),TMath::Max(unfoldedD04->GetMaximum(),unfoldedD05->GetMaximum())));
	recoSV4.SetLineColor(kBlue);
	recoSV5.SetLineColor(kRed);
	unfoldedSV4->SetLineColor(kBlue);
	unfoldedSV5->SetLineColor(kRed);
	unfoldedSV4->SetLineStyle(kDashed);
	unfoldedSV5->SetLineStyle(kDashed);
	trueSV4.SetLineColor(kBlue);
	trueSV5.SetLineColor(kRed);
	trueSV4.SetLineStyle(kDotted);
	trueSV5.SetLineStyle(kDotted);

	TCanvas c;
	recoD04.Draw();
	unfoldedD04->Draw("same");
	recoD05.Draw("same");
	unfoldedD05->Draw("same");
	if(dataIsMC) trueD04.Draw("same");
	if(dataIsMC) trueD05.Draw("same");
	c.SaveAs(savedir+"/d0Unfolding"+file+".pdf");

	recoSV4.Draw();
	unfoldedSV4->Draw("same");
	recoSV5.Draw("same");
	unfoldedSV5->Draw("same");
	if(dataIsMC) trueSV4.Draw("same");
	if(dataIsMC) trueSV5.Draw("same");
	c.SaveAs(savedir+"/svUnfolding"+file+".pdf");

	return 0;
}
