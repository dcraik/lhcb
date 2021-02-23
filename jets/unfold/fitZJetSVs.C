#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <boost/progress.hpp>
#include <boost/program_options.hpp>

#include "outputFunctions.h"
#include "DatasetManager.h"
#include "SVFitter.h"
#include "ZFitter.h"
#include "MCJets.h"

namespace po = boost::program_options;

//globals to save passing these around
//TString savedir("output-fitZj-statOnly-24bins-newFit-new-new-fixSVSel-fixJetPVMatch2-ptCorr-rerun");
TString savedir("output-fitZj-statOnly-24bins-newFit-new-new-fixSVSel-fixJetPVMatch2-ptCorr-rerun-new-rerun-overflowBin");

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
//TString charmSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_14X.root";
//TString beautySimFile = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_15X.root";
//TString lightSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_101.root";
//TString dataFile      = "/tmp/dcraik/skimmed.root";
TString charmSimFile  = "/data/zjet/zjet_sim4.root";
TString beautySimFile = "/data/zjet/zjet_sim5.root";
TString lightSimFile  = "/data/zjet/zjet_sim0.root";
TString lightDataFile  = "/data/zjet/zjet_back_201X.root";
TString lightSVFile  = "/data/zjet/zjet_sim0_sv.root";
TString dataFile      = "/data/zjet/zjet_201X.root";
TString ssDataFile    = "/data/zjet/zjet_ss_201X.root";
TString trainFile      = "/data/dijets/dijets_2016.root";
//TString tagEffFile = "dijets_20200519_4ptbins/taggingEffs.root";
//TString templateFilesPattern = "enhancedTemplates_%d_0.root";
TString tagEffFile = "taggingEffs.root";
TString templateFilesPattern = "enhancedTemplates_%d_0.root";

//TString lightHistFile = "svFitHists0_zj.root";
//TString charmHistFile = "svFitHists4_zj.root";
//TString beautyHistFile = "svFitHists5_zj.root";
//TString dataHistFile = "svFitHistsD_zj.root";

//unfolding
TString unfoldingMode("bayes");
double unfoldingReg(1.);
double jetEnergyScale(0.95);
double jetEnergySmear(0.14);

void fillJetsHist(TH2D* h, bool isMC) {
	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset("data");
	//TFile* f = TFile::Open(dataFile);
	//if(!f) return;
	//TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	//if(!t) return;

	double JetPT, ZPZ, ZE;

	int NPV;//TODO
	dm->setBranchAddress("NPV", &NPV);//TODO
	dm->setBranchAddress("JetPT", &JetPT);
	if(isMC) {
		dm->setBranchAddress("ZTRUEPZ",   &ZPZ);
		dm->setBranchAddress("ZTRUEE",    &ZE);
	} else {
		dm->setBranchAddress("ZPZ",   &ZPZ);
		dm->setBranchAddress("ZE",    &ZE);
	}

	std::cout << dm->getEntries() << std::endl;
	//for(int i=0; i<t->GetEntries(); ++i) {
	while(dm->getNext()) {
		//t->GetEntry(i);
		//if(NPV==1) continue;//TODO
		//std::cout << JetPT << "\t" << 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)) << "\t" << h->FindBin(JetPT,0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ))) << std::endl;
		h->Fill(JetPT,0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)));
	}
	//TString cutStr="0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)):JetPT>>";
	//cutStr+=h->GetName();
	//t->Draw(cutStr,"","goff");

	dm->reset();
}

int main(int argc, char** argv) {
	//print arguments used so we have them in the log file
	for(int iarg=0; iarg<argc; ++iarg) {
		printf("%s ",argv[iarg]);
	}
	setLHCbStyle();
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kError;
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

	std::vector<double> ptBinBoundaries;
	std::vector<double> yBinBoundaries;
	SVFitterOptions svOpts;
	bool useBackDataForLightShape;
	bool splitUnfoldingInY;
	TString inputDir;
	// Declare the supported options.
	po::options_description general_options("General");
	general_options.add_options()
	    ("help", "produce help message")
	    ("data", po::value<std::string>(), "input dataset [201X]")
	    ("dir", po::value<std::string>(), "save directory [output]")
	    ("pt-bins", po::value<std::vector<double>>(&ptBinBoundaries)->multitoken()->zero_tokens()->composing()->default_value(std::vector<double>{15e3,20e3,30e3,50e3,100e3}, "15e3, 20e3, 30e3, 50e3, 100e3"), "boundaries for pT bins [MeV]")
	    ("y-bins", po::value<std::vector<double>>(&yBinBoundaries)->multitoken()->zero_tokens()->composing()->default_value(std::vector<double>{2.,2.5,3.,3.5,4.5}, "2., 2.5, 3., 3.5, 4.5"), "boundaries for rapidity bins")
	    ("input-dir", po::value<std::string>(), "directory containing tagging efficiencies file and enhanced SV templates (defaults to dir)")
	    ("tag-eff-file", po::value<std::string>(), "input tagging efficiencies")
	;
	po::options_description svfit_options("SV fit options");
	svfit_options.add_options()
	    ("sv-nmcorbins", po::value<int>(&svOpts.nMCorBins)->default_value(94), "number of bins to use for MCor in SV fit (for binned fits and plots)")
	    ("sv-minmcor", po::value<double>(&svOpts.minMCor)->default_value(600), "minimum value of corrected mass to use in SV fits (in MeV)")
	    ("sv-maxmcor", po::value<double>(&svOpts.maxMCor)->default_value(10000), "maximum value of corrected mass to use in SV fits (in MeV)")
	    ("sv-nntrkbins", po::value<int>(&svOpts.nNTrkBins)->default_value(3), "number of bins to use for NTrk in SV fit")
	    ("sv-binnedtemplates", po::value<bool>(&svOpts.useBinnedTemplates)->default_value(false)->implicit_value(true), "use SV templates binned in MCor")
	    ("sv-ptbintemplates", po::value<bool>(&svOpts.usePtBinnedTemplates)->default_value(false)->implicit_value(true), "bin simulated samples in jet pT to obtain SV fit templates")
	    ("sv-templates-files", po::value<std::string>(), "input enhanced SV templates")
	    ("sv-mistag-shape-from-back", po::value<bool>(&useBackDataForLightShape)->default_value(false)->implicit_value(true), "use the back-tagged data file for the mis-tag SV shape (instead of simulation)")
	    ("sv-light-yield-float", po::value<bool>(&svOpts.lightYieldFloat)->default_value(false)->implicit_value(true), "float the mistag yields")
	    ("sv-light-yield-scale", po::value<double>(&svOpts.lightYieldScale)->default_value(1.), "scale the fixed mistag yields by a constant factor")
	;
	po::options_description unfold_options("Unfolding options");
	unfold_options.add_options()
	    ("unfold-mode", po::value<std::string>(), "unfolding method (invert, [bayes], ids, svd, tunfold}")
	    ("unfold-reg", po::value<double>(), "unfolding regularisation parameter (-1 for method default) [-1]")
	    ("jet-energy-scale", po::value<double>(&jetEnergyScale), "scale factor to correct jet energy in simulation [0.95]")
	    ("jet-energy-smear", po::value<double>(&jetEnergySmear), "factor to smear simulated jet energy by [0.14]")
	    ("unfold-rapiditybins", po::value<bool>(&splitUnfoldingInY)->default_value(false)->implicit_value(true), "whether to train unfolding separately for each bin of Z rapidity")
	;

	po::options_description desc("Options");
	desc.add(general_options).add(svfit_options).add(unfold_options);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    
	
	if (vm.count("help")) {
		std::cout << desc << std::endl;
	    return 1;
	}

	printf("Full configuration:\n");
	for ( auto& [k,v] : vm) {
		if(v.value().type() == typeid(std::string)) {
			printf("%35s:\t\"%s\"\n", k.data(), v.as<std::string>().data());
		} else if(v.value().type() == typeid(bool)) {
			printf("%35s:\t% 7d\n", k.data(), v.as<bool>());
		} else if(v.value().type() == typeid(int)) {
			printf("%35s:\t% 7d\n", k.data(), v.as<int>());
		} else if(v.value().type() == typeid(uint)) {
			printf("%35s:\t% 7d\n", k.data(), v.as<uint>());
		} else if(v.value().type() == typeid(double)) {
			printf("%35s:\t% 7g\n", k.data(), v.as<double>());
		} else if(v.value().type() == typeid(std::vector<std::string>)) {
			printf("%35s:\t", k.data());
			auto vv = v.as<std::vector<std::string>>();
			for( auto i : vv) {
				printf(" %s", i.data());
			}
			printf("\n");
		} else if(v.value().type() == typeid(std::vector<int>)) {
			printf("%35s:\t", k.data());
			auto vv = v.as<std::vector<int>>();
			for( auto i : vv) {
				printf(" %d", i);
			}
			printf("\n");
		} else if(v.value().type() == typeid(std::vector<double>)) {
			printf("%35s:\t", k.data());
			auto vv = v.as<std::vector<double>>();
			for( auto i : vv) {
				printf(" %g", i);
			}
			printf("\n");
		}
	}

	TString file="201X";
	if (vm.count("data")) {
		file = vm["data"].as<std::string>();
		std::cout << "Running over dataset:" << file << std::endl;
	}
	dataFile = "/data/zjet/zjet_"+file+".root";

	if (vm.count("dir")) {
		savedir = vm["dir"].as<std::string>();
		std::cout << "Saving to directory:" << savedir << std::endl;
	}
	if (vm.count("input-dir")) {
		inputDir = vm["input-dir"].as<std::string>();
		std::cout << "Reading inputs from directory:" << inputDir << std::endl;
	} else {
		inputDir = savedir;
	}
	if (vm.count("tag-eff-file")) {
		tagEffFile = vm["tag-eff-file"].as<std::string>();
		std::cout << "Using tagging efficiencies:" << tagEffFile << std::endl;
	}
	if (vm.count("sv-templates-files")) {
		templateFilesPattern = vm["sv-templates-files"].as<std::string>();
		std::cout << "Using enhanced templates:" << templateFilesPattern << std::endl;
	}
	if (vm.count("unfold-mode")) {
		unfoldingMode = vm["unfold-mode"].as<std::string>();
		std::cout << "Using unfolding mode: " << unfoldingMode << std::endl;
	}
	if (vm.count("unfold-reg")) {
		unfoldingReg = vm["unfold-reg"].as<double>();
		std::cout << "Unfolding regularisation parameter set to " << unfoldingReg << std::endl;
	}

	gSaveDir = savedir;
	gSystem->Exec("mkdir -p "+savedir+"/fig");
	gSystem->Exec("mkdir -p "+savedir+"/svfig");
	gSystem->Exec("mkdir -p "+savedir+"/dat");
	gSystem->Exec("mkdir -p "+savedir+"/log");

	double nB(0.), eB(0.), nC(0.), eC(0.), nQ(0.), eQ(0.);//, nZ(0.), eZ(0.);

	if(ptBinBoundaries.size()<2) {
		std::cout << "Must define at least 2 pT bin boundaries" << std::endl;
		return 1;
	}
	if(yBinBoundaries.size()<2) {
		std::cout << "Must define at least 2 rapidity bin boundaries" << std::endl;
		return 1;
	}
	unsigned int npt = ptBinBoundaries.size()-1;
	//const int npt(4);
	//double ptBounds[npt+1] = {15000.,20000.,30000.,50000.,100000.};
	//const int nzy(2);
	//double zyBounds[nzy+1] = {2.0,3.0,4.5};
	//const int nzy(4);
	unsigned int nzy = yBinBoundaries.size()-1;
	//double zyBounds[nzy+1] = {2.0,2.5,3.0,3.5,4.5};//4.0,4.5};
	//double zyBoundsOneBin[2] = {2.0,4.5};
	std::vector<double> yBinBoundsOneBin;
	yBinBoundsOneBin.push_back(yBinBoundaries[0]);
	yBinBoundsOneBin.push_back(yBinBoundaries[nzy]);

	TH1D ptScheme("ptScheme","",npt,ptBinBoundaries.data());
	TH1D zyScheme("zyScheme","",nzy,yBinBoundaries.data());

	TH2D hC("hC","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hB("hB","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hJ("hJ","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hEc("hEc","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hEb("hEb","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hE4("hE4","",npt,ptBinBoundaries.data(),1,yBinBoundsOneBin.data());
	TH2D hE5("hE5","",npt,ptBinBoundaries.data(),1,yBinBoundsOneBin.data());
	TH2D hRc("hRc","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURc("hURc","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRb("hRb","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURb("hURb","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hCE("hCE","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUCE("hUCE","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBE("hBE","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUBE("hUBE","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	//TH2D hZ("hZ","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	//TH2D hCSS("hCSS","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	//TH2D hBSS("hBSS","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	//histograms for propagating uncertainties separately
	TH2D hEcFull("hEcFull","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBFcErr("hBFcErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hEcNoErr("hEcNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hCNoErr("hCNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hJNoErr("hJNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUCNoErr("hUCNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUJNoErr("hUJNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRc_statErr("hRc_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURc_statErr("hURc_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRc_effErr("hRc_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURc_effErr("hURc_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRc_bfErr("hRc_bfErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURc_bfErr("hURc_bfErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	TH2D hEbFull("hEbFull","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBFbErr("hBFbErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hEbNoErr("hEbNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBNoErr("hBNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUBNoErr("hUBNoErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRb_statErr("hRb_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURb_statErr("hURb_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRb_effErr("hRb_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURb_effErr("hURb_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hRb_bfErr("hRb_bfErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hURb_bfErr("hURb_bfErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	TH2D hCE_statErr("hCE_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUCE_statErr("hUCE_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBE_statErr("hBE_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUBE_statErr("hUBE_statErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	TH2D hCE_effErr("hCE_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUCE_effErr("hUCE_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hBE_effErr("hBE_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hUBE_effErr("hUBE_effErr","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	//truth histograms
	TH2D hT("hT","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hT4("hT4","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hT5("hT5","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTSV4("hTSV4","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTSV5("hTSV5","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTU("hTU","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTU4("hTU4","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTU5("hTU5","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTUSV4("hTUSV4","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());
	TH2D hTUSV5("hTUSV5","",npt,ptBinBoundaries.data(),nzy,yBinBoundaries.data());

	TH1D hRc20up("hRc20up","",nzy,yBinBoundaries.data());
	TH1D hURc20up("hURc20up","",nzy,yBinBoundaries.data());
	TH1D hRb20up("hRb20up","",nzy,yBinBoundaries.data());
	TH1D hURb20up("hURb20up","",nzy,yBinBoundaries.data());

	hC.Sumw2();
	hJ.Sumw2();
	hEc.Sumw2();
	hEb.Sumw2();
	hE4.Sumw2();
	hE5.Sumw2();
	hRc.Sumw2();
	hURc.Sumw2();
	hRc_statErr.Sumw2();
	hURc_statErr.Sumw2();
	hRc_effErr.Sumw2();
	hURc_effErr.Sumw2();
	hRc_bfErr.Sumw2();
	hURc_bfErr.Sumw2();
	hRb_statErr.Sumw2();
	hURb_statErr.Sumw2();
	hRb_effErr.Sumw2();
	hURb_effErr.Sumw2();
	hRb_bfErr.Sumw2();
	hURb_bfErr.Sumw2();
	hCE.Sumw2();
	hUCE.Sumw2();
	hBE.Sumw2();
	hUBE.Sumw2();
	//hZ.Sumw2();
	//hCSS.Sumw2();
	//hBSS.Sumw2();

	hRc20up.Sumw2();
	hURc20up.Sumw2();
	hRb20up.Sumw2();
	hURb20up.Sumw2();

	//unsigned int nmcpt=4;
	//double* mcPtInputWeightBins  = new double[nmcpt +1]{10000.,15000.,20000.,50000.,100000.};

	//TH1D* jetTruePtWeights = new TH1D("jetTruePtWeights","",nmcpt,mcPtInputWeightBins);
	//jetTruePtWeights->SetBinContent(1,1.);
	//jetTruePtWeights->SetBinContent(2,0.6);
	//jetTruePtWeights->SetBinContent(3,0.5);
	//jetTruePtWeights->SetBinContent(4,0.04);

	bool dataIsMC(false);
	DatasetManager* dm = DatasetManager::getInstance();
	if(file.IsDec() && atoi(file)<100) {
		std::cout << "Running MC closure test " << file << " (" << atoi(file) << ")" << std::endl;
		gRandom->SetSeed(atoi(file));
		dataIsMC=true;

		dm->loadDataset("!light", "T",std::vector<TString>{lightSimFile},110000);
		dm->loadDataset("!charm", "T",std::vector<TString>{charmSimFile},5000);
		dm->loadDataset("!beauty","T",std::vector<TString>{beautySimFile},2000);
		dm->buildDataset("data",   "T",std::vector<TString>{"!charm","!beauty","!light"});
	} else {
		dm->loadDataset("light", "T",std::vector<TString>{lightSimFile});
		dm->loadDataset("charm", "T",std::vector<TString>{charmSimFile});
		dm->loadDataset("beauty","T",std::vector<TString>{beautySimFile});
		dm->loadDataset("data",  "T",std::vector<TString>{dataFile});
		dm->loadDataset("ssdata",  "T",std::vector<TString>{ssDataFile});
		dm->loadDataset("back",  "T",std::vector<TString>{lightDataFile});
		//dm->loadDataset("trainSVs","T",std::vector<TString>{trainFile});
	}

	MCJets mcj("Zj");
	mcj.setInputs("light","charm","beauty","","");
	//mcj.setInputTruePtWeights(MCJets::jetRecoSV4,jetTruePtWeights);
	//mcj.setInputTruePtWeights(MCJets::jetRecoSV5,jetTruePtWeights);
	if(!mcj.weightMC(MCJets::jetRecoSV4,&hE4)) return 1;
	if(!mcj.weightMC(MCJets::jetRecoSV5,&hE5)) return 1;

	if(dataIsMC) {
		mcj.getTruth("data",&hT,&hT4,&hT5,&hTSV4,&hTSV5);
		mcj.getTruth("data",&hTU,&hTU4,&hTU5,&hTUSV4,&hTUSV5,true);
	}

	bool usingDataEffsC(false);
	bool usingDataEffsB(false);
	double bfffErr(0.);
	double bfffErrB(0.);
	TFile* feff = TFile::Open(inputDir+"/"+tagEffFile);
	if(feff) {
		TH1D* h = static_cast<TH1D*>(feff->Get("tagEff4"));
		if(h) {
			if(h->GetNbinsX() == static_cast<int>(npt)) {
				for(uint i=0; i<npt; ++i) {
					for(uint j=0; j<nzy; ++j) {
						hEc.SetBinContent(i+1,j+1,h->GetBinContent(i+1));
						hEc.SetBinError(i+1,j+1,h->GetBinError(i+1));
						hEcNoErr.SetBinContent(i+1,j+1,h->GetBinContent(i+1));
						hEcNoErr.SetBinError(i+1,j+1,0.);
					}
				}
				usingDataEffsC=true;
				std::cout << "Will use tagging efficiencies from input data file" << std::endl;
			}
		}
		h = static_cast<TH1D*>(feff->Get("tagEff5"));
		if(h) {
			if(h->GetNbinsX() == static_cast<int>(npt)) {
				for(uint i=0; i<npt; ++i) {
					for(uint j=0; j<nzy; ++j) {
						hEb.SetBinContent(i+1,j+1,h->GetBinContent(i+1));
						hEb.SetBinError(i+1,j+1,h->GetBinError(i+1));
						hEbNoErr.SetBinContent(i+1,j+1,h->GetBinContent(i+1));
						hEbNoErr.SetBinError(i+1,j+1,0.);
					}
				}
				usingDataEffsB=true;
				std::cout << "Will use tagging efficiencies from input data file" << std::endl;
			}
		}
		TH1D* hbfff = static_cast<TH1D*>(feff->Get("bfffErr4"));
		if(hbfff) {
			bfffErr = hbfff->GetBinError(1);
			for(uint i=0; i<npt; ++i) {
				for(uint j=0; j<nzy; ++j) {
					hBFcErr.SetBinContent(i+1,j+1,1.);
					hBFcErr.SetBinError(i+1,j+1,bfffErr);
				}
			}
		}
		hbfff = static_cast<TH1D*>(feff->Get("bfffErr5"));
		if(hbfff) {
			bfffErrB = hbfff->GetBinError(1);
			for(uint i=0; i<npt; ++i) {
				for(uint j=0; j<nzy; ++j) {
					hBFbErr.SetBinContent(i+1,j+1,1.);
					hBFbErr.SetBinError(i+1,j+1,bfffErrB);
				}
			}
		}
		feff->Close();
	}

	//for(int j=1; j<=nzy; ++j) {
	//	hEc.SetBinContent(1,j,0.191);
	//	hEc.SetBinError(1,j,0.);//0.021);
	//	hEc.SetBinContent(2,j,0.237);
	//	hEc.SetBinError(2,j,0.);//0.016);
	//	hEc.SetBinContent(3,j,0.241);
	//	hEc.SetBinError(3,j,0.);//0.022);
	//}

	fillJetsHist(&hJ,dataIsMC);
	for(uint i=0; i<npt; ++i) {
		for(uint j=0; j<nzy; ++j) {
			hJNoErr.SetBinContent(i+1,j+1,hJ.GetBinContent(i+1,j+1));
			hJNoErr.SetBinError(i+1,j+1,0.);
		}
	}

	//ZFitter zfit("ZFit", &ptScheme, &zyScheme);
	//zfit.setInputs(dataFile);
	//zfit.fit(nZ,eZ);
	//zfit.fixShape();

	SVFitter svfit("ZjFit", &ptScheme, &zyScheme);
	dm->loadDataset("charmSV", "T",std::vector<TString>{mcj.outputName(MCJets::jetRecoSV4)});
	dm->loadDataset("beautySV","T",std::vector<TString>{mcj.outputName(MCJets::jetRecoSV5)});
	if(useBackDataForLightShape) {
		svfit.setInputs("back","charmSV","beautySV","data",dataIsMC,true,true,dataIsMC,"","back");
	} else {
		svfit.setInputs("light","charmSV","beautySV","data",true,true,true,dataIsMC,"","back");
	}
	svfit.setInputWeightings(false,true,true);
	svfit.setOptions(svOpts);
	//svfit.setSVBinning(18,500.,5000.,3);
	//svfit.setSVBinning(94,600.,10000.,3);
	svfit.makeSVFitHists();

	//SVFitter sssvfit("SSFit", &ptScheme, &zyScheme);
	//sssvfit.setInputs(lightDataFile,mcj.outputName(MCJets::jetRecoSV4),mcj.outputName(MCJets::jetRecoSV5),ssDataFile);
	//if(useBackDataForLightShape) {
	//	sssvfit.setInputs("back","charmSV","beautySV","ssdata",dataIsMC,true,true,dataIsMC,"","back");
	//} else {
	//	sssvfit.setInputs("light","charmSV","beautySV","ssdata",true,true,true,dataIsMC,"","back");
	//}
	//sssvfit.setInputWeightings(false,true,true);
	//sssvfit.setOptions(svOpts);
	////svfit.setSVBinning(18,500.,5000.,3);
	//sssvfit.makeSVFitHists();

	for(uint i=1; i<=npt; ++i) {
		for(uint j=1; j<=nzy; ++j) {
			if(!usingDataEffsC) {
				hEc.SetBinContent(i,j,hE4.GetBinContent(i,1));
				hEc.SetBinError  (i,j,hE4.GetBinError  (i,1));
			}
			if(!usingDataEffsB) {
				hEb.SetBinContent(i,j,hE5.GetBinContent(i,1));
				hEb.SetBinError  (i,j,hE5.GetBinError  (i,1));
			}
			TString templateFile = "";
			if(templateFilesPattern!="") {
				templateFile = TString::Format(inputDir+"/"+templateFilesPattern.Data(),i);
			}
			if(svfit.fit(nB,eB,nC,eC,nQ,eQ,i,j,templateFile)) {
				double corrC = mcj.getPtCorrFactor(MCJets::jetRecoSV4,ptBinBoundaries[i-1],ptBinBoundaries[i]);
				double corrB = mcj.getPtCorrFactor(MCJets::jetRecoSV5,ptBinBoundaries[i-1],ptBinBoundaries[i]);
				hC.SetBinContent(i,j,corrC*nC);
				hC.SetBinError(  i,j,corrC*eC);
				hCNoErr.SetBinContent(i,j,corrC*nC);
				hCNoErr.SetBinError(  i,j,0.);
				hB.SetBinContent(i,j,corrB*nB);
				hB.SetBinError(  i,j,corrB*eB);
				hBNoErr.SetBinContent(i,j,corrB*nB);
				hBNoErr.SetBinError(  i,j,0.);
			}
			//if(sssvfit.fit(nB,eB,nC,eC,nQ,eQ,i,j,"")) {
			//	hCSS.SetBinContent(i,j,nC);
			//	hCSS.SetBinError(  i,j,eC);
			//	hBSS.SetBinContent(i,j,nB);
			//	hBSS.SetBinError(  i,j,eB);
			//}
			//if(zfit.fit(nZ,eZ,i,j)) {
			//	hZ.SetBinContent(i,j,nZ);
			//	hZ.SetBinError(  i,j,eZ);
			//}
		}
	}

	mcj.setUnfoldingMethod(unfoldingMode,unfoldingReg);
	mcj.setUnfoldingScaleSmear(jetEnergyScale,jetEnergySmear);
	TH2D* hUC = mcj.unfold(&hC, MCJets::jetRecoSV4, splitUnfoldingInY);
	TH2D* hUB = mcj.unfold(&hB, MCJets::jetRecoSV5, splitUnfoldingInY);
	TH2D* hUJ = mcj.unfold(&hJ, MCJets::jetAll0, splitUnfoldingInY);

	for(uint i=0; i<npt; ++i) {
		for(uint j=0; j<nzy; ++j) {
			hUBNoErr.SetBinContent(i+1,j+1,hUB->GetBinContent(i+1,j+1));
			hUBNoErr.SetBinError(i+1,j+1,0.);
			hUCNoErr.SetBinContent(i+1,j+1,hUC->GetBinContent(i+1,j+1));
			hUCNoErr.SetBinError(i+1,j+1,0.);
			hUJNoErr.SetBinContent(i+1,j+1,hUJ->GetBinContent(i+1,j+1));
			hUJNoErr.SetBinError(i+1,j+1,0.);
		}
	}

	//hZ.Divide(&hJ);
	//hCSS.Divide(&hC);
	//hBSS.Divide(&hB);

	hRc.Divide(&hC,&hJ);
	hRc.Divide(&hEc);
	hURc.Divide(hUC,hUJ);
	hURc.Divide(&hEc);

	hRb.Divide(&hB,&hJ);
	hRb.Divide(&hEb);
	hURb.Divide(hUB,hUJ);
	hURb.Divide(&hEb);

	//versions with only one source of uncertainty
	hRc_statErr.Divide(&hC,&hJ);
	hRc_statErr.Divide(&hEcNoErr);
	hURc_statErr.Divide(hUC,hUJ);
	hURc_statErr.Divide(&hEcNoErr);

	hRc_effErr.Divide(&hCNoErr,&hJNoErr);
	hRc_effErr.Divide(&hEc);
	hURc_effErr.Divide(&hUCNoErr,&hUJNoErr);
	hURc_effErr.Divide(&hEc);

	hRc_bfErr.Divide(&hCNoErr,&hJNoErr);
	hRc_bfErr.Divide(&hEcNoErr);
	hRc_bfErr.Multiply(&hBFcErr);
	hURc_bfErr.Divide(&hUCNoErr,&hUJNoErr);
	hURc_bfErr.Divide(&hEcNoErr);
	hURc_bfErr.Multiply(&hBFcErr);

	hEcFull.Multiply(&hEc,&hBFcErr);

	hRb_statErr.Divide(&hB,&hJ);
	hRb_statErr.Divide(&hEbNoErr);
	hURb_statErr.Divide(hUB,hUJ);
	hURb_statErr.Divide(&hEbNoErr);

	hRb_effErr.Divide(&hBNoErr,&hJNoErr);
	hRb_effErr.Divide(&hEb);
	hURb_effErr.Divide(&hUBNoErr,&hUJNoErr);
	hURb_effErr.Divide(&hEb);

	hRb_bfErr.Divide(&hBNoErr,&hJNoErr);
	hRb_bfErr.Divide(&hEbNoErr);
	hRb_bfErr.Multiply(&hBFbErr);
	hURb_bfErr.Divide(&hUBNoErr,&hUJNoErr);
	hURb_bfErr.Divide(&hEbNoErr);
	hURb_bfErr.Multiply(&hBFbErr);

	hEbFull.Multiply(&hEb,&hBFcErr);

	//versions for pT>20 integrated numbers
	hCE.Divide(&hC,&hEc);
	hBE.Divide(&hB,&hEb);
	hUCE.Divide(hUC,&hEc);
	hUBE.Divide(hUB,&hEb);

	hCE_statErr.Divide(&hC,&hEcNoErr);
	hBE_statErr.Divide(&hB,&hEbNoErr);
	hUCE_statErr.Divide(hUC,&hEcNoErr);
	hUBE_statErr.Divide(hUB,&hEbNoErr);

	hCE_effErr.Divide(&hCNoErr,&hEc);
	hBE_effErr.Divide(&hBNoErr,&hEb);
	hUCE_effErr.Divide(&hUCNoErr,&hEc);
	hUBE_effErr.Divide(&hUBNoErr,&hEb);

	TH1D* hCE20up  = hCE .ProjectionY("hCE20up", 2,npt);
	TH1D* hBE20up  = hBE .ProjectionY("hBE20up", 2,npt);
	TH1D* hJ20up   = hJ  .ProjectionY("hJ20up",  2,npt);
	TH1D* hUCE20up = hUCE.ProjectionY("hUCE20up",2,npt);
	TH1D* hUBE20up = hUBE.ProjectionY("hUBE20up",2,npt);
	TH1D* hUJ20up  = hUJ->ProjectionY("hUJ20up", 2,npt);

	TH1D* hCE20up_statErr  = hCE_statErr .ProjectionY("hCE20up_statErr", 2,npt);
	TH1D* hBE20up_statErr  = hBE_statErr .ProjectionY("hBE20up_statErr", 2,npt);
	TH1D* hUCE20up_statErr = hUCE_statErr.ProjectionY("hUCE20up_statErr",2,npt);
	TH1D* hUBE20up_statErr = hUBE_statErr.ProjectionY("hUBE20up_statErr",2,npt);

	TH1D* hCE20up_effErr  = hCE_effErr .ProjectionY("hCE20up_effErr", 2,npt);
	TH1D* hBE20up_effErr  = hBE_effErr .ProjectionY("hBE20up_effErr", 2,npt);
	TH1D* hUCE20up_effErr = hUCE_effErr.ProjectionY("hUCE20up_effErr",2,npt);
	TH1D* hUBE20up_effErr = hUBE_effErr.ProjectionY("hUBE20up_effErr",2,npt);

	//TH1D* hPtCE  = hCE .ProjectionX("hPtCE");
	//TH1D* hPtJ   = hJ  .ProjectionX("hPtJ");
	//TH1D* hPtUCE = hUCE.ProjectionX("hPtUCE");
	//TH1D* hPtUJ  = hUJ->ProjectionX("hPtUJ");

	hRc20up.Divide(hCE20up,hJ20up);
	hRb20up.Divide(hBE20up,hJ20up);
	hURc20up.Divide(hUCE20up,hUJ20up);
	hURb20up.Divide(hUBE20up,hUJ20up);

	double Nc20up_statErr(0.);
	double Nb20up_statErr(0.);
	double Nj20up_statErr(0.);
	double UNc20up_statErr(0.);
	double UNb20up_statErr(0.);
	double UNj20up_statErr(0.);

	double Nc20up_effErr(0.);
	double Nb20up_effErr(0.);
	double UNc20up_effErr(0.);
	double UNb20up_effErr(0.);

	double Nc20up_bfErr(0.);
	double Nb20up_bfErr(0.);
	double UNc20up_bfErr(0.);
	double UNb20up_bfErr(0.);

	double Nc20up = hCE20up_statErr->IntegralAndError(1,nzy, Nc20up_statErr);
	double Nb20up = hBE20up_statErr->IntegralAndError(1,nzy, Nb20up_statErr);
	double Nj20up = hJ20up         ->IntegralAndError(1,nzy, Nj20up_statErr);
	double UNc20up = hUCE20up_statErr->IntegralAndError(1,nzy, UNc20up_statErr);
	double UNb20up = hUBE20up_statErr->IntegralAndError(1,nzy, UNb20up_statErr);
	double UNj20up = hUJ20up         ->IntegralAndError(1,nzy, UNj20up_statErr);

	Nc20up_effErr = Nc20up * hCE20up_effErr->GetBinError(1) / hCE20up_effErr->GetBinContent(1);
	Nb20up_effErr = Nb20up * hBE20up_effErr->GetBinError(1) / hBE20up_effErr->GetBinContent(1);
	UNc20up_effErr = UNc20up * hUCE20up_effErr->GetBinError(1) / hUCE20up_effErr->GetBinContent(1);
	UNb20up_effErr = UNb20up * hUBE20up_effErr->GetBinError(1) / hUBE20up_effErr->GetBinContent(1);

	Nc20up_bfErr = bfffErr*Nc20up;
	Nb20up_bfErr = bfffErrB*Nb20up;
	UNc20up_bfErr = bfffErr*UNc20up;
	UNb20up_bfErr = bfffErrB*UNb20up;

	double Rc20up =  Nc20up  / Nj20up;
	double Rb20up =  Nb20up  / Nj20up;
	double URc20up = UNc20up / UNj20up;
	double URb20up = UNb20up / UNj20up;

	double Rc20up_statErr  = Rc20up  * TMath::Sqrt(TMath::Power( Nc20up_statErr/ Nc20up,2.) + TMath::Power( Nj20up_statErr/ Nj20up,2.));
	double Rb20up_statErr  = Rb20up  * TMath::Sqrt(TMath::Power( Nb20up_statErr/ Nb20up,2.) + TMath::Power( Nj20up_statErr/ Nj20up,2.));
	double URc20up_statErr = URc20up * TMath::Sqrt(TMath::Power(UNc20up_statErr/UNc20up,2.) + TMath::Power(UNj20up_statErr/UNj20up,2.));
	double URb20up_statErr = URb20up * TMath::Sqrt(TMath::Power(UNb20up_statErr/UNb20up,2.) + TMath::Power(UNj20up_statErr/UNj20up,2.));

	double Rc20up_effErr  = Rc20up  *  Nc20up_effErr/ Nc20up;
	double Rb20up_effErr  = Rb20up  *  Nb20up_effErr/ Nb20up;
	double URc20up_effErr = URc20up * UNc20up_effErr/UNc20up;
	double URb20up_effErr = URb20up * UNb20up_effErr/UNb20up;

	double Rc20up_bfErr  = Rc20up  *  Nc20up_bfErr/ Nc20up;
	double Rb20up_bfErr  = Rb20up  *  Nb20up_bfErr/ Nb20up;
	double URc20up_bfErr = URc20up * UNc20up_bfErr/UNc20up;
	double URb20up_bfErr = URb20up * UNb20up_bfErr/UNb20up;

	hUJ20up->SetLineColor(kRed);
	hUJ20up->SetMarkerColor(kRed);
	hUCE20up->SetLineColor(kRed);
	hUCE20up->SetMarkerColor(kRed);
	hUCE20up->SetLineStyle(kDashed);
	hCE20up->SetLineStyle(kDashed);
	hURc20up.SetLineColor(kRed);
	hURc20up.SetMarkerColor(kRed);

	//hPtUJ->SetLineColor(kRed);
	//hPtUJ->SetMarkerColor(kRed);
	//hPtUCE->SetLineColor(kRed);
	//hPtUCE->SetMarkerColor(kRed);
	//hPtUCE->SetLineStyle(kDashed);
	//hPtCE->SetLineStyle(kDashed);

	TCanvas c;
	//hRc20up.SetMaximum(1.1*TMath::Max(hRc20up.GetMaximum(),hURc20up.GetMaximum()));
	//hRc20up.SetMinimum(0.);
	//hRc20up.GetXaxis()->SetTitle("#it{y}(#it{Z})");
	//hRc20up.GetYaxis()->SetTitle("#it{f}_{#it{c}}");
	//hRc20up.Draw();
	//hURc20up.Draw("same");
	//c.SaveAs(savedir+"/fracC_20up.pdf");
	//hJ20up->SetMaximum(1.1*TMath::Max(hJ20up->GetMaximum(),hUJ20up->GetMaximum()));
	//hJ20up->SetMinimum(0.);
	//hJ20up->GetXaxis()->SetTitle("#it{y}(#it{Z})");
	//hJ20up->GetYaxis()->SetTitle("#it{N}_{jet}");
	//hJ20up->Draw();
	//hUJ20up->Draw("same");
	//hCE20up->Draw("same");
	//hUCE20up->Draw("same");
	//c.SaveAs(savedir+"/NJNC_20up.pdf");
	//hPtJ->SetMaximum(1.1*TMath::Max(hPtJ->GetMaximum(),hPtUJ->GetMaximum()));
	//hPtJ->SetMinimum(0.);
	//hPtJ->GetXaxis()->SetTitle("#it{p}_{T}(#it{j})");
	//hPtJ->GetYaxis()->SetTitle("#it{N}_{jet}");
	//hPtJ->Draw();
	//hPtUJ->Draw("same");
	//hPtCE->Draw("same");
	//hPtUCE->Draw("same");
	//c.SaveAs(savedir+"/NJNC_pT.pdf");

	TH1D* hURc1520 = hURc.ProjectionY("hURc1520",1,1);
	TH1D* hURc2030 = hURc.ProjectionY("hURc2030",2,2);
	TH1D* hURc3050 = hURc.ProjectionY("hURc3050",3,3);
	TH1D* hURc50up = hURc.ProjectionY("hURc50up",4,4);

	//add correlated uncertainties after projections
	TH1D* hURc1520Full = static_cast<TH1D*>(hURc1520->Clone("hURc1520Full"));
	TH1D* hURc2030Full = static_cast<TH1D*>(hURc2030->Clone("hURc2030Full"));
	TH1D* hURc3050Full = static_cast<TH1D*>(hURc3050->Clone("hURc3050Full"));
	TH1D* hURc50upFull = static_cast<TH1D*>(hURc50up->Clone("hURc50upFull"));
	TH1D* hURc20upFull = static_cast<TH1D*>(hURc20up .Clone("hURc20upFull"));
	for(uint ybin=1; ybin<=nzy; ++ybin) {
		hURc1520Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURc1520->GetBinError(ybin),2.) + TMath::Power(bfffErr*hURc1520->GetBinContent(ybin),2.)));
		hURc2030Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURc2030->GetBinError(ybin),2.) + TMath::Power(bfffErr*hURc2030->GetBinContent(ybin),2.)));
		hURc3050Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURc3050->GetBinError(ybin),2.) + TMath::Power(bfffErr*hURc3050->GetBinContent(ybin),2.)));
		hURc50upFull->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURc50up->GetBinError(ybin),2.) + TMath::Power(bfffErr*hURc50up->GetBinContent(ybin),2.)));
		hURc20upFull->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURc20up .GetBinError(ybin),2.) + TMath::Power(bfffErr*hURc20up .GetBinContent(ybin),2.)));
		//std::cout << hURc1520Full->GetBinError(ybin) << " " << hURc1520->GetBinError(ybin) << " " << bfffErr*hURc1520->GetBinContent(ybin) << std::endl;
	}

	hURc20up.SetLineColor(kBlack);
	hURc20up.SetMarkerColor(kBlack);
	hURc50up->SetLineColor(kMagenta);
	hURc50up->SetMarkerColor(kMagenta);
	hURc3050->SetLineColor(kRed);
	hURc3050->SetMarkerColor(kRed);
	hURc2030->SetLineColor(kGreen+2);
	hURc2030->SetMarkerColor(kGreen+2);
	hURc1520->SetLineColor(kBlue);
	hURc1520->SetMarkerColor(kBlue);
	hURc20upFull->SetLineColor(kBlack);
	hURc20upFull->SetMarkerColor(kBlack);
	hURc50upFull->SetLineColor(kMagenta);
	hURc50upFull->SetMarkerColor(kMagenta);
	hURc3050Full->SetLineColor(kRed);
	hURc3050Full->SetMarkerColor(kRed);
	hURc2030Full->SetLineColor(kGreen+2);
	hURc2030Full->SetMarkerColor(kGreen+2);
	hURc1520Full->SetLineColor(kBlue);
	hURc1520Full->SetMarkerColor(kBlue);
	hURc1520->SetMinimum(0.);
	hURc1520->SetMaximum(0.15);
	hURc1520->GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	hURc1520->GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{c}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	hURc1520->Draw("E1 P");
	hURc1520Full->Draw("E1 P SAME");
	if(npt>2) {
		hURc50up->Draw("E1 P same");
		hURc50upFull->Draw("E1 P same");
		hURc3050->Draw("E1 P same");
		hURc3050Full->Draw("E1 P same");
		hURc2030->Draw("E1 P same");
		hURc2030Full->Draw("E1 P same");
	}
	hURc20up.Draw("E1 P same");
	hURc20upFull->Draw("E1 P same");
	c.SaveAs(savedir+"/fracCharmZj.pdf");
	hURc1520->Draw("E1 P");
	hURc1520Full->Draw("E1 P SAME");
	if(npt>2) {
		hURc50up->Draw("E1 P same");
		hURc50upFull->Draw("E1 P same");
		hURc3050->Draw("E1 P same");
		hURc3050Full->Draw("E1 P same");
	}
	hURc2030->Draw("E1 P same");
	hURc2030Full->Draw("E1 P same");
	c.SaveAs(savedir+"/fracCharmZj_sep.pdf");
	//hURc1520->Draw("E1 P");
	//hURc1520Full->Draw("E1 P SAME");
	hURc20up.SetMinimum(0.);
	hURc20up.SetMaximum(0.15);
	hURc20up.GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	hURc20up.GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{c}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	hURc20up.Draw("E1 P");
	hURc20upFull->Draw("E1 P same");
	c.SaveAs(savedir+"/fracCharmZj_combined.pdf");

	TH1D* hURb1520 = hURb.ProjectionY("hURb1520",1,1);
	TH1D* hURb2030 = hURb.ProjectionY("hURb2030",2,2);
	TH1D* hURb3050 = hURb.ProjectionY("hURb3050",3,3);
	TH1D* hURb50up = hURb.ProjectionY("hURb50up",4,4);

	//add correlated uncertainties after projections
	TH1D* hURb1520Full = static_cast<TH1D*>(hURb1520->Clone("hURb1520Full"));
	TH1D* hURb2030Full = static_cast<TH1D*>(hURb2030->Clone("hURb2030Full"));
	TH1D* hURb3050Full = static_cast<TH1D*>(hURb3050->Clone("hURb3050Full"));
	TH1D* hURb50upFull = static_cast<TH1D*>(hURb50up->Clone("hURb50upFull"));
	TH1D* hURb20upFull = static_cast<TH1D*>(hURb20up .Clone("hURb20upFull"));
	for(uint ybin=1; ybin<=nzy; ++ybin) {
		hURb1520Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURb1520->GetBinError(ybin),2.) + TMath::Power(bfffErrB*hURb1520->GetBinContent(ybin),2.)));
		hURb2030Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURb2030->GetBinError(ybin),2.) + TMath::Power(bfffErrB*hURb2030->GetBinContent(ybin),2.)));
		hURb3050Full->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURb3050->GetBinError(ybin),2.) + TMath::Power(bfffErrB*hURb3050->GetBinContent(ybin),2.)));
		hURb50upFull->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURb50up->GetBinError(ybin),2.) + TMath::Power(bfffErrB*hURb50up->GetBinContent(ybin),2.)));
		hURb20upFull->SetBinError(ybin, TMath::Sqrt(TMath::Power(hURb20up .GetBinError(ybin),2.) + TMath::Power(bfffErrB*hURb20up .GetBinContent(ybin),2.)));
	}

	hURb20up.SetLineColor(kBlack);
	hURb20up.SetMarkerColor(kBlack);
	hURb50up->SetLineColor(kMagenta);
	hURb50up->SetMarkerColor(kMagenta);
	hURb3050->SetLineColor(kRed);
	hURb3050->SetMarkerColor(kRed);
	hURb2030->SetLineColor(kGreen+2);
	hURb2030->SetMarkerColor(kGreen+2);
	hURb1520->SetLineColor(kBlue);
	hURb1520->SetMarkerColor(kBlue);
	hURb20upFull->SetLineColor(kBlack);
	hURb20upFull->SetMarkerColor(kBlack);
	hURb50upFull->SetLineColor(kMagenta);
	hURb50upFull->SetMarkerColor(kMagenta);
	hURb3050Full->SetLineColor(kRed);
	hURb3050Full->SetMarkerColor(kRed);
	hURb2030Full->SetLineColor(kGreen+2);
	hURb2030Full->SetMarkerColor(kGreen+2);
	hURb1520Full->SetLineColor(kBlue);
	hURb1520Full->SetMarkerColor(kBlue);
	hURb1520->SetMinimum(0.);
	hURb1520->SetMaximum(0.15);
	hURb1520->GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	hURb1520->GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{b}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	hURb1520->Draw("E1 P");
	hURb1520Full->Draw("E1 P SAME");
	if(npt>2) {
		hURb50up->Draw("E1 P same");
		hURb50upFull->Draw("E1 P same");
		hURb3050->Draw("E1 P same");
		hURb3050Full->Draw("E1 P same");
		hURb2030->Draw("E1 P same");
		hURb2030Full->Draw("E1 P same");
	}
	hURb20up.Draw("E1 P same");
	hURb20upFull->Draw("E1 P same");
	c.SaveAs(savedir+"/fracBeautyZj.pdf");
	hURb1520->Draw("E1 P");
	hURb1520Full->Draw("E1 P SAME");
	if(npt>2) {
		hURb50up->Draw("E1 P same");
		hURb50upFull->Draw("E1 P same");
		hURb3050->Draw("E1 P same");
		hURb3050Full->Draw("E1 P same");
	}
	hURb2030->Draw("E1 P same");
	hURb2030Full->Draw("E1 P same");
	c.SaveAs(savedir+"/fracBeautyZj_sep.pdf");
	//hURb1520->Draw("E1 P");
	//hURb1520Full->Draw("E1 P SAME");
	hURb20up.SetMinimum(0.);
	hURb20up.SetMaximum(0.15);
	hURb20up.GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	hURb20up.GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{b}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	hURb20up.Draw("E1 P");
	hURb20upFull->Draw("E1 P same");
	c.SaveAs(savedir+"/fracBeautyZj_combined.pdf");

	//hEc.Draw("colz text45");
	//c.SaveAs(savedir+"/charmEffs.pdf");
	//hC.Draw("colz text45");
	//c.SaveAs(savedir+"/charmYields.pdf");
	//hJ.Draw("colz text45");
	//c.SaveAs(savedir+"/totalYields.pdf");
	//hRc.Draw("colz text45");
	//c.SaveAs(savedir+"/charmFractions.pdf");
	////hZ.Draw("colz text45");
	////c.SaveAs(savedir+"/fracZ.pdf");
	//hCSS.Draw("colz text45");
	//c.SaveAs(savedir+"/fracCharmBkg.pdf");
	//hBSS.Draw("colz text45");
	//c.SaveAs(savedir+"/fracBeautyBkg.pdf");

	for(uint i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			//std::cout << ptBinBoundaries[i-1] << "-" << ptBinBoundaries[i] << "," << zyBounds[j-1] << "-" << zyBounds[j] << ":\t" << hJ.GetBinContent(i,j) << "+/-" << hJ.GetBinError(i,j) << "\t" << hC.GetBinContent(i,j) << "+/-" << hC.GetBinError(i,j) << "\t" << hEc.GetBinContent(i,j) << "+/-" << hEc.GetBinError(i,j) << "\t" << hRc.GetBinContent(i,j) << "+/-" << hRc.GetBinError(i,j) << std::endl;
			printf("%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f\\pm%.3f\\pm%.3f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hJ.GetBinContent(i,j), hJ.GetBinError(i,j), hC.GetBinContent(i,j), hC.GetBinError(i,j), hEcFull.GetBinContent(i,j), hEcFull.GetBinError(i,j), hRc.GetBinContent(i,j), hRc_statErr.GetBinError(i,j), hRc_effErr.GetBinError(i,j), hRc_bfErr.GetBinError(i,j));
		}
	}
	std::cout << "unfolded" << std::endl;
	for(uint i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			//std::cout << ptBinBoundaries[i-1] << "-" << ptBinBoundaries[i] << "," << zyBounds[j-1] << "-" << zyBounds[j] << ":\t" << hUJ->GetBinContent(i,j) << "+/-" << hUJ->GetBinError(i,j) << "\t" << hUC->GetBinContent(i,j) << "+/-" << hUC->GetBinError(i,j) << "\t" << hEc.GetBinContent(i,j) << "+/-" << hEc.GetBinError(i,j) << "\t" << hURc.GetBinContent(i,j) << "+/-" << hURc.GetBinError(i,j) << std::endl;
			printf("%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f\\pm%.3f\\pm%.3f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hUJ->GetBinContent(i,j), hUJ->GetBinError(i,j), hUC->GetBinContent(i,j), hUC->GetBinError(i,j), hEcFull.GetBinContent(i,j), hEcFull.GetBinError(i,j), hURc.GetBinContent(i,j), hURc_statErr.GetBinError(i,j), hURc_effErr.GetBinError(i,j), hURc_bfErr.GetBinError(i,j));
		}
	}

	for(uint i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			printf("%3.0f--%3.0f & %.2f--%.2f & $%.4f\\pm%.4f\\pm%.4f\\pm%.4f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hURc.GetBinContent(i,j), hURc_statErr.GetBinError(i,j), hURc_effErr.GetBinError(i,j), hURc_bfErr.GetBinError(i,j));
		}
	}
	for(uint i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			printf("%3.0f--%3.0f & %.2f--%.2f & $%.4f\\pm%.4f\\pm%.4f\\pm%.4f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hURb.GetBinContent(i,j), hURb_statErr.GetBinError(i,j), hURb_effErr.GetBinError(i,j), hURb_bfErr.GetBinError(i,j));
		}
	}

	FILE* rFile = fopen(savedir+"/results.log","w");
	for(uint i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			fprintf(rFile,"%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hJ.GetBinContent(i,j), hJ.GetBinError(i,j), hC.GetBinContent(i,j), hC.GetBinError(i,j), hEcFull.GetBinContent(i,j), hEcFull.GetBinError(i,j), 1e2*hRc.GetBinContent(i,j), 1e2*hRc_statErr.GetBinError(i,j), 1e2*hRc_effErr.GetBinError(i,j), 1e2*hRc_bfErr.GetBinError(i,j));
		}
	}
	fprintf(rFile,"unfolded\n");
	for(uint i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			fprintf(rFile,"%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hUJ->GetBinContent(i,j), hUJ->GetBinError(i,j), hUC->GetBinContent(i,j), hUC->GetBinError(i,j), hEcFull.GetBinContent(i,j), hEcFull.GetBinError(i,j), 1e2*hURc.GetBinContent(i,j), 1e2*hURc_statErr.GetBinError(i,j), 1e2*hURc_effErr.GetBinError(i,j), 1e2*hURc_bfErr.GetBinError(i,j));
		}
	}
	fprintf(rFile,"beauty\n");
	for(uint i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			fprintf(rFile,"%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hJ.GetBinContent(i,j), hJ.GetBinError(i,j), hB.GetBinContent(i,j), hB.GetBinError(i,j), hEbFull.GetBinContent(i,j), hEbFull.GetBinError(i,j), 1e2*hRb.GetBinContent(i,j), 1e2*hRb_statErr.GetBinError(i,j), 1e2*hRb_effErr.GetBinError(i,j), 1e2*hRb_bfErr.GetBinError(i,j));
		}
	}
	fprintf(rFile,"unfolded\n");
	for(uint i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(uint j=1; j<=nzy; ++j) {
			fprintf(rFile,"%3.0f--%3.0f & %.2f--%.2f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f\\pm%.2f\\pm%.2f$ \\\\\n", ptBinBoundaries[i-1]/1e3, ptBinBoundaries[i]/1e3, yBinBoundaries[j-1], yBinBoundaries[j], hUJ->GetBinContent(i,j), hUJ->GetBinError(i,j), hUB->GetBinContent(i,j), hUB->GetBinError(i,j), hEbFull.GetBinContent(i,j), hEbFull.GetBinError(i,j), 1e2*hURb.GetBinContent(i,j), 1e2*hURb_statErr.GetBinError(i,j), 1e2*hURb_effErr.GetBinError(i,j), 1e2*hURb_bfErr.GetBinError(i,j));
		}
	}
	fclose(rFile);

	std::cout << "pT>20GeV, 2.0<y<4.5" << std::endl;
	printf("c %.4f\\pm%.4f\\pm%.4f\\pm%.4f\n", Rc20up, Rc20up_statErr, Rc20up_effErr, Rc20up_bfErr);
	printf("b %.4f\\pm%.4f\\pm%.4f\\pm%.4f\n", Rb20up, Rb20up_statErr, Rb20up_effErr, Rb20up_bfErr);
	std::cout << "unfolded" << std::endl;
	printf("c %.4f\\pm%.4f\\pm%.4f\\pm%.4f\n", URc20up, URc20up_statErr, URc20up_effErr, URc20up_bfErr);
	printf("b %.4f\\pm%.4f\\pm%.4f\\pm%.4f\n", URb20up, URb20up_statErr, URb20up_effErr, URb20up_bfErr);

	//inputs for python script
	std::ofstream write;
	write.open(savedir+"/resultsOutput.log");
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUC->GetBinContent(i,j) << " ";
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUB->GetBinContent(i,j) << " ";
	write << std::endl;
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUC->GetBinError(i,j) << " ";
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUB->GetBinError(i,j) << " ";
	write << std::endl;
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUJ->GetBinContent(i,j) << " ";
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUJ->GetBinContent(i,j) << " ";
	write << std::endl;
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUJ->GetBinError(i,j) << " ";
	for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j) write << hUJ->GetBinError(i,j) << " ";
	write << std::endl;
	for (unsigned int i=1; i<=npt; ++i)  /*for(uint j=1; j<=nzy; ++j)*/ write << hEc.GetBinContent(i,1/*j*/) << " ";//shared between y bins
	for (unsigned int i=1; i<=npt; ++i)  /*for(uint j=1; j<=nzy; ++j)*/ write << hEb.GetBinContent(i,1/*j*/) << " ";//shared between y bins
	write << std::endl;
	for (unsigned int i=1; i<=npt; ++i)  /*for(uint j=1; j<=nzy; ++j)*/ write << hEc.GetBinError(i,1/*j*/) << " ";//shared between y bins
	for (unsigned int i=1; i<=npt; ++i)  /*for(uint j=1; j<=nzy; ++j)*/ write << hEb.GetBinError(i,1/*j*/) << " ";//shared between y bins
	write << std::endl;
	/*for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j)*/ write << hBFcErr.GetBinError(1,1/*i,j*/) << " ";//shared between pT and y bins
	/*for (unsigned int i=1; i<=npt; ++i)  for(uint j=1; j<=nzy; ++j)*/ write << hBFbErr.GetBinError(1,1/*i,j*/) << " ";//shared between pT and y bins
	write << std::endl;
	write.close();

	TFile* fout = TFile::Open(savedir+"/hists.root","RECREATE");
	hC.Write();
	//hCSS.Write();
	hB.Write();
	//hBSS.Write();
	hJ.Write();
	hEc.Write();
	hEb.Write();
	hRc.Write();
	hRb.Write();
	hUC->Write();
	hUB->Write();
	hUJ->Write();
	hURc.Write();
	hURb.Write();
	hURc1520Full->Write();
	hURc2030Full->Write();
	hURc3050Full->Write();
	hURc50upFull->Write();
	hURc20upFull->Write();
	hURb1520Full->Write();
	hURb2030Full->Write();
	hURb3050Full->Write();
	hURb50upFull->Write();
	hURb20upFull->Write();
	hT.Write();
	hT4.Write();
	hT5.Write();
	hTSV4.Write();
	hTSV5.Write();
	hTU.Write();
	hTU4.Write();
	hTU5.Write();
	hTUSV4.Write();
	hTUSV5.Write();
	fout->Close();

	return 0;
}
