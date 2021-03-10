#include <vector>
#include <fstream>

#include <boost/progress.hpp>
#include <boost/program_options.hpp>

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

#include "DatasetManager.h"
//#include "DFitter.h"
#include "SimDFitter.h"
#include "MCJets.h"
#include "SVFitter.h"

namespace po = boost::program_options;

//globals to save passing these around
TString savedir("output");
//bool useSimpleEff(false);
bool useRhoZEffCor(true);
bool weightRecPt(false);

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
bool scaledSim = false;

TString charmSimFile  = "/data/dijets/dijets_sim4_reco.root";
TString beautySimFile = "/data/dijets/dijets_sim5_reco.root";
TString scaledCharmSimFile  = "/data/dijets/dijets_sim4_scaled_reco.root";//TODO
TString scaledBeautySimFile = "/data/dijets/dijets_sim5_scaled_reco.root";
//TString scaledCharmSimFile  = "/data/dijets/dijets_sim4_scaled2_reco.root";
//TString scaledBeautySimFile = "/data/dijets/dijets_sim5_scaled2_reco.root";
TString lightDataFile  = "/data/dijets/dijets_back_2016.root";
TString lightSimFile  = "/data/dijets/dijets_sim0_reco.root";
TString dataFile      = "/data/dijets/dijets_2016.root";

//efficiency inputs - update these with latest version numbers
//TString simpleEffFile = "../efficiencies/SimpleEffs_2XX_16x8bins_up190216.root";
TString oldEffFileD0 = "../efficiencies/D0Effs_2XX_16x8bins_up190213.root";
//TString effFileD0 = "../efficiencies/D0Effs_45_16x8bins_up201129.root";
TString effFileD0 = "../efficiencies/D0Effs_45_16x8bins_up210222.root";
TString accFileD0 = "../efficiencies/D0AccEffNewUp200518.root";
//TString accFileD0 = "../efficiencies/D0AccEffNewUp201118_fromD0.root";
//TString accFileD0 = "../efficiencies/D0AccEffNewUp201118_fromC.root";
//TString accFileD0 = "../efficiencies/D0AccEffNewUp201118_fromAll.root";
TString effFileDp = "../efficiencies/DpEffs_2XX_16x8bins_up200506.root";
TString accFileDp = "../efficiencies/DpAccEffNewUp200513.root";

//unfolding
TString unfoldingMode("bayes");
double unfoldingReg(4.);
double jetEnergyScale(0.95);
double jetEnergySmear(0.14);

bool dataIsMC(false);
bool dataIsResampledMC(false);

enum TaggingStage {
	stageEffs,
	stageWeight,
	stageSVTemplates,
	stageFittingD,
	stageFittingSV,
	stageUnfolding,
	stageEnd
};

int main(int argc, char** argv) {
	//print arguments used so we have them in the log file
	for(int iarg=0; iarg<argc; ++iarg) {
		printf("%s ",argv[iarg]);
	}
	std::cout << std::endl;
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kError;
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
//TODO	RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);//TODO TODO

	//options to pass to the D/SV fitters
	std::vector<double> ptBinBoundaries;
	bool evtByEvtWeighting;
	std::set<std::string> dfits;
	DFitterOptions dfitOpts;
	bool useOldEff(false);
	SVFitterOptions svOpts;
	bool useBackDataForLightShape;
	uint svEnhancePasses;
	std::string ptCorrFile;
	//START OF COMMANDLINE CONFIGURATION
	// Declare the supported options.
	po::options_description general_options("General");
	general_options.add_options()
	    ("help", "produce help message")
	    ("data", po::value<std::string>(), "input dataset (1-999 for toys) [2016]")
	    ("dir", po::value<std::string>(), "save directory [output]")
	    ("skip-to", po::value<int>(), "stage to skip to (if cached) [5]:\n  0: D0 efficiencies\n  1: MC weighting\n  2: SV templates\n  3: D0 fit\n  4: SV fit\n  5: Unfolding")
	    ("force-stages", po::value<std::vector<int>>()->multitoken()->zero_tokens()->composing(), "Stages to be forced to run")
	    ("pt-bins", po::value<std::vector<double>>(&ptBinBoundaries)->multitoken()->zero_tokens()->composing()->default_value(std::vector<double>{15e3,20e3,30e3,50e3,100e3}, "15e3, 20e3, 30e3, 50e3, 100e3"), "boundaries for pT bins [MeV]")
	    ("weight-sim", po::value<int>(), "weight true jet pT spectrum of simulated samples to be continuous (-1: default, 0: off, 1: on). Default to off for simulated datasets and on for real datasets")
	    ("scale-sim", po::value<bool>(&scaledSim)->default_value(false)->implicit_value(true), "used simulation scaled to match data in jet(pT) distribution")
	    ("use-evt-by-evt-weighting", po::value<bool>(&evtByEvtWeighting)->default_value(false)->implicit_value(true), "calculate denominator form -event-by-event weighted D0 candidates instead of average-efficiency weighted")
	    ("pt-corr-file",po::value<std::string>(&ptCorrFile)->default_value("")->implicit_value("ptCorrFactors.root"), "File to take hadron-pT correction factors from")
	;
	po::options_description dfit_options("D fit options");
	dfit_options.add_options()
	    ("dfits", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "D fits to run [D0]:\n  D0: D0 -> Kpi")
	    ("dfit-minpt", po::value<double>(&dfitOpts.minpt)->default_value(5000.), "minimum D0 pT in MeV")
	    ("dfit-maxpt", po::value<double>(&dfitOpts.maxpt)->default_value(-1.), "maximum D0 pT in MeV (-1 to turn off)")
	    ("dfit-skip-sumw2-fits", po::value<bool>(&dfitOpts.skipSumw2Fits)->default_value(false)->implicit_value(true), "skip the second stage of D0 fits used to correctly scale uncertainties for fits to weighted datasets")
	    ("dfit-truth-match",po::value<uint>(&dfitOpts.truthMatch)->default_value(0), "truth match data (if simulated) 0: off, 4: prompt, 5: displaced")
	    ("dfit-n-pt-bins",po::value<uint>(&dfitOpts.nDPtBins)->default_value(5u), "number of bins in pT(D0)")
	    ("dfit-use-ptfrac-bins",po::value<bool>(&dfitOpts.usePtFracBins)->default_value(false)->implicit_value(true), "bin in pT(D0)/pT(j) instead of pT(D0)")
	    ("dfit-comb-shape-syst",po::value<bool>(&dfitOpts.combShapeSyst)->default_value(false)->implicit_value(true), "use a smaller region of the D sidebands to evaluate the background chi^2_IP shape")
	    ("dfit-enhanced-fits",po::value<bool>(&dfitOpts.doEnhancedFits)->default_value(false)->implicit_value(true), "perform initial fits to datasets enhanced in beauty and charm to correct MC shape parameters")
	    ("dfit-allowMassWidthScale",po::value<bool>(&dfitOpts.allowMassWidthScale)->default_value(true)->implicit_value(true), "Allow a scale factor in sigma(D) between data and simulation")
	    ("dfit-allowPromptMeanShift",po::value<bool>(&dfitOpts.allowPromptMeanShift)->default_value(true)->implicit_value(true), "Allow a shift in mu_prompt between data and simulation")
	    ("dfit-allowPromptWidthScale",po::value<bool>(&dfitOpts.allowPromptWidthScale)->default_value(false)->implicit_value(true), "Allow a scale factor in sigma_prompt between data and simulation")
	    ("dfit-allowDisplcMeanShift",po::value<bool>(&dfitOpts.allowDisplcMeanShift)->default_value(false)->implicit_value(true), "Allow a shift in mu_displc between data and simulation")
	    ("dfit-allowDisplcWidthScale",po::value<bool>(&dfitOpts.allowDisplcWidthScale)->default_value(false)->implicit_value(true), "Allow a scale factor in sigma_displc between data and simulation")
	    ("dfit-enhanceMassWidthScale",po::value<bool>(&dfitOpts.enhanceMassWidthScale)->default_value(false)->implicit_value(true), "Allow scale factor in sigma(D) between data and simulation for enhanced data fits")
	    ("dfit-enhancePromptMeanShift",po::value<bool>(&dfitOpts.enhancePromptMeanShift)->default_value(false)->implicit_value(true), "Allow shift in mu_prompt between data and simulation for enhanced data fits")
	    ("dfit-enhancePromptWidthScale",po::value<bool>(&dfitOpts.enhancePromptWidthScale)->default_value(false)->implicit_value(true), "Allow scale factor in sigma_prompt between data and simulation for enhanced data fits")
	    ("dfit-enhanceDisplcMeanShift",po::value<bool>(&dfitOpts.enhanceDisplcMeanShift)->default_value(false)->implicit_value(true), "Allow shift in mu_displc between data and simulation for enhanced data fits")
	    ("dfit-enhanceDisplcWidthScale",po::value<bool>(&dfitOpts.enhanceDisplcWidthScale)->default_value(false)->implicit_value(true), "Allow scale factor in sigma_displc between data and simulation for enhanced data fits")
	    ("dfit-splitMassWidthScale",po::value<bool>(&dfitOpts.splitMassWidthScale)->default_value(false)->implicit_value(true), "Allow a bin-by-bin scale factor in sigma(D) between data and simulation")
	    ("dfit-splitPromptMeanShift",po::value<bool>(&dfitOpts.splitPromptMeanShift)->default_value(false)->implicit_value(true), "Allow a bin-by-bin shift in mu_prompt between data and simulation")
	    ("dfit-splitPromptWidthScale",po::value<bool>(&dfitOpts.splitPromptWidthScale)->default_value(false)->implicit_value(true), "Allow a bin-by-bin scale factor in sigma_prompt between data and simulation")
	    ("dfit-splitDisplcMeanShift",po::value<bool>(&dfitOpts.splitDisplcMeanShift)->default_value(false)->implicit_value(true), "Allow a bin-by-bin shift in mu_displc between data and simulation")
	    ("dfit-splitDisplcWidthScale",po::value<bool>(&dfitOpts.splitDisplcWidthScale)->default_value(false)->implicit_value(true), "Allow a bin-by-bin scale factor in sigma_displc between data and simulation")
	    ("dfit-splitBkgMassShape",po::value<bool>(&dfitOpts.splitBkgMassShape)->default_value(false)->implicit_value(true), "Allow a bin-by-bin gradient for the background mass shape")
	    ("dfit-splitBkgIPShape",po::value<bool>(&dfitOpts.splitBkgIPShape)->default_value(false)->implicit_value(true), "Allow a bin-by-bin KDE for the background IP shape")
	    ("dfit-allowLinearBkgMassShift",po::value<bool>(&dfitOpts.allowLinearBkgMassShift)->default_value(false)->implicit_value(true), "Allow a bin-by-bin gradient for the background mass shape but require the gradient to be a linear function of the average value of pT(D0) or pT(D0)/pT(j) in each bin")
	    ("dfit-setMassWidthScale",po::value<double>(&dfitOpts.setMassWidthScale)->default_value(-999.), "Set a scale factor in sigma(D) between data and simulation")
	    ("dfit-setPromptMeanShift",po::value<double>(&dfitOpts.setPromptMeanShift)->default_value(-999.), "Set a shift in mu_prompt between data and simulation")
	    ("dfit-setPromptWidthScale",po::value<double>(&dfitOpts.setPromptWidthScale)->default_value(-999.), "Set a scale factor in sigma_prompt between data and simulation")
	    ("dfit-setDisplcMeanShift",po::value<double>(&dfitOpts.setDisplcMeanShift)->default_value(-999.), "Set a shift in mu_displc between data and simulation")
	    ("dfit-setDisplcWidthScale",po::value<double>(&dfitOpts.setDisplcWidthScale)->default_value(-999.), "Set a scale factor in sigma_displc between data and simulation")
	    ("dfit-setDisplcWidthScale",po::value<double>(&dfitOpts.setDisplcWidthScale)->default_value(-999.), "Set a scale factor in sigma_displc between data and simulation")
	    ("dfit-simpleEffs",po::value<bool>(&dfitOpts.useSimpleEffs)->default_value(false)->implicit_value(true), "Use simplified model for D efficiencies")
	    ("dfit-effCorrFile4",po::value<std::string>(), "File to take prompt D0 efficiency corrections from (based on simulation closure test)")
	    ("dfit-effCorrFile5",po::value<std::string>(), "File to take displaced D0 efficiency corrections from (based on simulation closure test)")
	    ("dfit-toyStudy",po::value<bool>(&dfitOpts.runToyFits)->default_value(false)->implicit_value(true), "Run a toy study after each fit to data")
	    ("dfit-oldEff",po::value<bool>(&useOldEff)->default_value(false)->implicit_value(true), "Use old D0 efficiency files")
	;
	po::options_description svfit_options("SV fit options");
	svfit_options.add_options()
	    ("sv-nmcorbins", po::value<int>(&svOpts.nMCorBins)->default_value(94), "number of bins to use for MCor in SV fit (for binned fits and plots)")
	    ("sv-minmcor", po::value<double>(&svOpts.minMCor)->default_value(600), "minimum value of corrected mass to use in SV fits (in MeV)")
	    ("sv-maxmcor", po::value<double>(&svOpts.maxMCor)->default_value(10000), "maximum value of corrected mass to use in SV fits (in MeV)")
	    ("sv-nntrkbins", po::value<int>(&svOpts.nNTrkBins)->default_value(3), "number of bins to use for NTrk in SV fit")
	    ("sv-binnedtemplates", po::value<bool>(&svOpts.useBinnedTemplates)->default_value(false)->implicit_value(true), "use SV templates binned in MCor")
	    ("sv-smear-mcor", po::value<bool>(&svOpts.smearMCorShapes)->default_value(false)->implicit_value(true), "correct SV templates by smearing in MCor")
	    ("sv-additive-mcor-corr", po::value<bool>(&svOpts.additiveMCorCorrectionFactor)->default_value(false)->implicit_value(true), "When correcting the mcor shapes, perform polynomial fits to the residual rather than the correction scale factor")
	    ("sv-ptbintemplates", po::value<bool>(&svOpts.usePtBinnedTemplates)->default_value(false)->implicit_value(true), "bin simulated samples in jet pT to obtain SV fit templates")
	    ("sv-mistag-shape-from-back", po::value<bool>(&useBackDataForLightShape)->default_value(false)->implicit_value(true), "use the back-tagged data file for the mis-tag SV shape (instead of simulation)")
	    ("sv-light-yield-float", po::value<bool>(&svOpts.lightYieldFloat)->default_value(false)->implicit_value(true), "float the mistag yields")
	    ("sv-enhancement-passes", po::value<uint>(&svEnhancePasses)->default_value(2u), "number of times to run SV template enhancement")
	    ("sv-toyStudy",po::value<bool>(&svOpts.runToyFits)->default_value(false)->implicit_value(true), "Run a toy study after each fit to data")
	;
	po::options_description unfold_options("Unfolding options");
	unfold_options.add_options()
	    ("unfold-mode", po::value<std::string>(), "unfolding method (invert, [bayes], ids, svd, tunfold}")
	    ("unfold-reg", po::value<double>(), "unfolding regularisation parameter (-1 for method default) [-1]")
	    ("jet-energy-scale", po::value<double>(&jetEnergyScale), "scale factor to correct jet energy in simulation [0.95]")
	    ("jet-energy-smear", po::value<double>(&jetEnergySmear), "factor to smear simulated jet energy by [0.14]")
	;

	po::options_description desc("Options");
	desc.add(general_options).add(dfit_options).add(svfit_options).add(unfold_options);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	//print help text
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

	//set input data and output directory
	//0-999 for MC closure tests
	TString file="2016";
	if (vm.count("data")) {
		file = vm["data"].as<std::string>();
		std::cout << "Running over dataset:" << file << std::endl;
	}
	if (vm.count("dir")) {
		savedir = vm["dir"].as<std::string>();
		std::cout << "Saving to directory:" << savedir << std::endl;
	}

	//setup stages to run
	std::set<TaggingStage> runStages;
	TaggingStage skipTo = stageUnfolding;
	if (vm.count("skip-to")) {
		skipTo = static_cast<TaggingStage>(vm["skip-to"].as<int>());
		std::cout << "Will use cached results up to stage " << skipTo << std::endl;
	}
	//add non-skipped stages to the set to run
	for(int stage = static_cast<int>(skipTo); stage<static_cast<int>(stageEnd); ++stage) {
		runStages.insert(static_cast<TaggingStage>(stage));
	}
	//add any forced rerun stages
	if (vm.count("force-stages")) {
		auto stages = vm["force-stages"].as<std::vector<int> >();
		for(auto stage : stages) {
			runStages.insert(static_cast<TaggingStage>(stage));
			std::cout << "Will force run stage " << stage << std::endl;
		}
	}

	//simulation configuration
	int weightSimTruePt(-1);
	if (vm.count("weight-sim")) {
		weightSimTruePt = vm["weight-sim"].as<int>();
		if(weightSimTruePt>0) {
			std::cout << "Weighting simulated samples to obtain continuous true jet pT spectrum" << std::endl;
		} else if(weightSimTruePt==0) {
			std::cout << "Not weighting simulated samples in true jet pT" << std::endl;
		}
	}

	//D configuration
	if (vm.count("dfits")) {
		auto stages = vm["dfits"].as<std::vector<std::string> >();
		for(auto stage : stages) {
			dfits.insert(stage);
		}
	}

	if (vm.count("dfit-minpt")) {
		std::cout << "Minimum D0 pT set to " << dfitOpts.minpt << std::endl;
	}
	if (vm.count("dfit-maxpt")) {
		std::cout << "Maximum D0 pT set to " << dfitOpts.maxpt << std::endl;
	}
	if (vm.count("dfit-skip-sumw2-fits")) {
		if(dfitOpts.skipSumw2Fits) std::cout << "Will not run sumW2 fits to D0 datasets" << std::endl;
	}
	if (vm.count("dfit-truth-match")) {
		if(dfitOpts.truthMatch == 4) {
			std::cout << "Truth matching for prompt D0 decays in data" << std::endl;
		}
		if(dfitOpts.truthMatch == 5) {
			std::cout << "Truth matching for displaced D0 decays in data" << std::endl;
		}
	}
	TString d0EffCorrFile4="";
	TString d0EffCorrFile5="";
	if (vm.count("dfit-effCorrFile4")) {
		d0EffCorrFile4 = vm["dfit-effCorrFile4"].as<std::string>();
	}
	if (vm.count("dfit-effCorrFile5")) {
		d0EffCorrFile5 = vm["dfit-effCorrFile5"].as<std::string>();
	}
	if(useOldEff) effFileD0 = oldEffFileD0;

	//unfolding configuration
	if (vm.count("unfold-mode")) {
		unfoldingMode = vm["unfold-mode"].as<std::string>();
		std::cout << "Using unfolding mode: " << unfoldingMode << std::endl;
	}
	if (vm.count("unfold-reg")) {
		unfoldingReg = vm["unfold-reg"].as<double>();
		std::cout << "Unfolding regularisation parameter set to " << unfoldingReg << std::endl;
	}
	//END OF COMMANDLINE CONFIGURATION

	gSystem->Exec("mkdir -p "+savedir+"/fig");
	gSystem->Exec("mkdir -p "+savedir+"/svfig");
	gSystem->Exec("mkdir -p "+savedir+"/dat");
	gSystem->Exec("mkdir -p "+savedir+"/log");
	gSystem->Exec("mkdir -p "+savedir+"/dtoys");
	gSystem->Exec("mkdir -p "+savedir+"/svtoys");
	gSaveDir = savedir;

	if(scaledSim) {
		charmSimFile = scaledCharmSimFile;
		beautySimFile = scaledBeautySimFile;
	}

	DatasetManager* dm = DatasetManager::getInstance();
	if(file.IsDec() && atoi(file)<1000) {
		std::cout << "Running MC closure test " << file << " (" << atoi(file) << ")" << std::endl;
		gRandom->SetSeed(atoi(file));
		dataIsMC=true;
		dataIsResampledMC=true;

		dm->loadDataset("light",   "T",std::vector<TString>{lightSimFile});
		//dm->loadDataset("!charm", "T",std::vector<TString>{charmSimFile},200000);
		//dm->loadDataset("!beauty","T",std::vector<TString>{beautySimFile},200000);
		dm->loadDataset("charm", "T",std::vector<TString>{charmSimFile});
		dm->loadDataset("beauty","T",std::vector<TString>{beautySimFile});
		dm->loadDataset("datacharm", "T",std::vector<TString>{scaledCharmSimFile},200000);
		dm->loadDataset("databeauty","T",std::vector<TString>{scaledBeautySimFile},200000);
		dm->loadDataset("datalight", "T",std::vector<TString>{lightDataFile},300000);
		dm->buildDataset("data",   "T",std::vector<TString>{"datacharm","databeauty","datalight"});
		dm->loadDataset("back", "T",std::vector<TString>{lightDataFile});
	} else {
		if(file.BeginsWith("sim")) {
			dataIsMC=true;
		}
		dataFile = "/data/dijets/dijets_"+file+".root";
		dm->loadDataset("light", "T",std::vector<TString>{lightSimFile});
		dm->loadDataset("charm", "T",std::vector<TString>{charmSimFile});
		dm->loadDataset("beauty","T",std::vector<TString>{beautySimFile});
		dm->loadDataset("data",  "T",std::vector<TString>{dataFile});
		dm->loadDataset("back", "T",std::vector<TString>{lightDataFile});
	}

	//these histograms determined by hand from simulated samples such that JetTruePt is roughly continuous
	//used to weight MC for continuous true jet pT when performing fits to data
	unsigned int nmcpt=4;
	double* ptInputWeightBins  = new double[nmcpt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* jetTruePtWeights4  = new TH1D("jetTruePtWeights4", "",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights4d = new TH1D("jetTruePtWeights4d","",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights4s = new TH1D("jetTruePtWeights4s","",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights5  = new TH1D("jetTruePtWeights5", "",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights5d = new TH1D("jetTruePtWeights5d","",nmcpt,ptInputWeightBins);
	TH1D* jetTruePtWeights5s = new TH1D("jetTruePtWeights5s","",nmcpt,ptInputWeightBins);

	jetTruePtWeights4->SetBinContent(1,1.);
	jetTruePtWeights4->SetBinContent(2,0.6);
	jetTruePtWeights4->SetBinContent(3,0.6);//1.2);
	jetTruePtWeights4->SetBinContent(4,0.03);//0.07);
	jetTruePtWeights5->SetBinContent(1,1.);
	jetTruePtWeights5->SetBinContent(2,0.5);//0.6);
	jetTruePtWeights5->SetBinContent(3,0.4);//0.5);
	jetTruePtWeights5->SetBinContent(4,0.02);
	//D0/SV-sample-specific numbers
	jetTruePtWeights4d->SetBinContent(1,1.  *0.8);
	jetTruePtWeights4d->SetBinContent(2,0.6 *1.1);
	jetTruePtWeights4d->SetBinContent(3,0.6 *0.7);
	jetTruePtWeights4d->SetBinContent(4,0.03*1.0);
	jetTruePtWeights5d->SetBinContent(1,1.  *1.0);
	jetTruePtWeights5d->SetBinContent(2,0.5 *1.0);
	jetTruePtWeights5d->SetBinContent(3,0.4 *1.0);
	jetTruePtWeights5d->SetBinContent(4,0.02*1.0);
	jetTruePtWeights4s->SetBinContent(1,1.  *1.1);
	jetTruePtWeights4s->SetBinContent(2,0.6 *0.9);
	jetTruePtWeights4s->SetBinContent(3,0.6 *0.9);
	jetTruePtWeights4s->SetBinContent(4,0.03*1.0);
	jetTruePtWeights5s->SetBinContent(1,1.  *1.0);
	jetTruePtWeights5s->SetBinContent(2,0.5 *1.1);
	jetTruePtWeights5s->SetBinContent(3,0.4 *0.9);
	jetTruePtWeights5s->SetBinContent(4,0.02*1.0);
	//TODO tweaks to get better MC/data agreement in each unfolding

	//bining scheme in reco pT(jet)
	//unsigned int npt=5;
	//double* binsPt  = new double[npt +1]{17500.,20000.,30000.,50000.};
	//double* binsPt  = new double[npt +1]{10000.,15000.,20000.,30000.,50000.,100000.};
	//double* binsPt  = new double[npt +1]{15000.,20000.,30000.,100000.};
	if(ptBinBoundaries.size()<2) {
		std::cout << "Must define at least 2 pT bin boundaries" << std::endl;
		return 1;
	}
	unsigned int npt = ptBinBoundaries.size()-1;

	//histograms for reconstructed yields
	//D0
	TH1D recoD04_evtByEvtWeight("recoD04_evtByEvtWeight","",npt,ptBinBoundaries.data());
	TH1D recoD05_evtByEvtWeight("recoD05_evtByEvtWeight","",npt,ptBinBoundaries.data());
	TH1D recoD04_aveWeight("recoD04_aveWeight","",npt,ptBinBoundaries.data());
	TH1D recoD05_aveWeight("recoD05_aveWeight","",npt,ptBinBoundaries.data());
	TH1D recoD0Sel4("recoD0Sel4","",npt,ptBinBoundaries.data());
	TH1D recoD0Sel5("recoD0Sel5","",npt,ptBinBoundaries.data());
	TH1D recoD0Weight4("recoD0Weight4","",npt,ptBinBoundaries.data());
	TH1D recoD0Weight5("recoD0Weight5","",npt,ptBinBoundaries.data());
	//D+
	TH1D recoDp4_evtByEvtWeight("recoDp4_evtByEvtWeight","",npt,ptBinBoundaries.data());
	TH1D recoDp5_evtByEvtWeight("recoDp5_evtByEvtWeight","",npt,ptBinBoundaries.data());
	TH1D recoDp4_aveWeight("recoDp4_aveWeight","",npt,ptBinBoundaries.data());
	TH1D recoDp5_aveWeight("recoDp5_aveWeight","",npt,ptBinBoundaries.data());
	TH1D recoDpSel4("recoDpSel4","",npt,ptBinBoundaries.data());
	TH1D recoDpSel5("recoDpSel5","",npt,ptBinBoundaries.data());
	TH1D recoDpWeight4("recoDpWeight4","",npt,ptBinBoundaries.data());
	TH1D recoDpWeight5("recoDpWeight5","",npt,ptBinBoundaries.data());
	//SV
	TH1D recoSV4("recoSV4","",npt,ptBinBoundaries.data());
	TH1D recoSV5("recoSV5","",npt,ptBinBoundaries.data());
	TH1D reco4("reco4","",npt,ptBinBoundaries.data());
	TH1D reco5("reco5","",npt,ptBinBoundaries.data());
	TH1D unfolded4("unfolded4","",npt,ptBinBoundaries.data());
	TH1D unfolded5("unfolded5","",npt,ptBinBoundaries.data());

	//truth histograms for MC studies
	TH1D trueD04("trueD04","",npt,ptBinBoundaries.data());
	TH1D trueD05("trueD05","",npt,ptBinBoundaries.data());
	TH1D trueD0Sel4("trueD0Sel4","",npt,ptBinBoundaries.data());
	TH1D trueD0Sel5("trueD0Sel5","",npt,ptBinBoundaries.data());
	TH1D trueDp4("trueDp4","",npt,ptBinBoundaries.data());
	TH1D trueDp5("trueDp5","",npt,ptBinBoundaries.data());
	TH1D trueDpSel4("trueDpSel4","",npt,ptBinBoundaries.data());
	TH1D trueDpSel5("trueDpSel5","",npt,ptBinBoundaries.data());
	TH1D trueSV4("trueSV4","",npt,ptBinBoundaries.data());
	TH1D trueSV5("trueSV5","",npt,ptBinBoundaries.data());
	TH1D true4("true4","",npt,ptBinBoundaries.data());
	TH1D true5("true5","",npt,ptBinBoundaries.data());

	//truth histograms for unfolded
	TH1D trueUnfldD04("trueUnfldD04","",npt,ptBinBoundaries.data());
	TH1D trueUnfldD05("trueUnfldD05","",npt,ptBinBoundaries.data());
	TH1D trueUnfldD0Sel4("trueUnfldD0Sel4","",npt,ptBinBoundaries.data());
	TH1D trueUnfldD0Sel5("trueUnfldD0Sel5","",npt,ptBinBoundaries.data());
	TH1D trueUnfldDp4("trueUnfldDp4","",npt,ptBinBoundaries.data());
	TH1D trueUnfldDp5("trueUnfldDp5","",npt,ptBinBoundaries.data());
	TH1D trueUnfldDpSel4("trueUnfldDpSel4","",npt,ptBinBoundaries.data());
	TH1D trueUnfldDpSel5("trueUnfldDpSel5","",npt,ptBinBoundaries.data());
	TH1D trueUnfldSV4("trueUnfldSV4","",npt,ptBinBoundaries.data());
	TH1D trueUnfldSV5("trueUnfldSV5","",npt,ptBinBoundaries.data());
	TH1D trueUnfld4("trueUnfld4","",npt,ptBinBoundaries.data());
	TH1D trueUnfld5("trueUnfld5","",npt,ptBinBoundaries.data());

	//corrections for D0 efficiency non-closure
	TH1D effClosure4("effClosure4","",npt,ptBinBoundaries.data());
	TH1D effClosure5("effClosure5","",npt,ptBinBoundaries.data());
	for(int i=1; i<=npt; ++i) {
		effClosure4.SetBinContent(i, 1.);
		effClosure5.SetBinContent(i, 1.);
		effClosure4.SetBinError(i, 0.);
		effClosure5.SetBinError(i, 0.);
	}

	if(d0EffCorrFile4) {
		TFile* ftmp = TFile::Open(d0EffCorrFile4);
		if(ftmp) {
			TH1D* htmp = static_cast<TH1D*>(ftmp->Get("closure"));
			if(htmp) {
				for(int i=1; i<=npt; ++i) {
					effClosure4.SetBinContent(i, htmp->GetBinContent(i));
				}
			}
		}
	}
	//for(int i=1; i<=npt; ++i) {
	//	std::cout << "!" << effClosure4.GetBinContent(i) << std::endl;//TODO
	//}

	if(d0EffCorrFile5) {
		TFile* ftmp = TFile::Open(d0EffCorrFile5);
		if(ftmp) {
			TH1D* htmp = static_cast<TH1D*>(ftmp->Get("closure"));
			if(htmp) {
				for(int i=1; i<=npt; ++i) {
					effClosure5.SetBinContent(i, htmp->GetBinContent(i));
				}
			}
		}
	}

	double* binsY  = new double[2]{0.,10.};
	TH2D hE4("hE4","",npt,ptBinBoundaries.data(),1,binsY);
	TH2D hE5("hE5","",npt,ptBinBoundaries.data(),1,binsY);
	TH2D hUE4("hUE4","",npt,ptBinBoundaries.data(),1,binsY);
	TH2D hUE5("hUE5","",npt,ptBinBoundaries.data(),1,binsY);

	//Setup the classes to do the work
	MCJets mcj("jj");
	mcj.setDPtRange(dfitOpts.minpt,dfitOpts.maxpt);
	SimDFitter d0fit("jjDFit"+file, &trueD04);
	SimDFitter dpfit("jjDFit"+file, &trueD04);
	SVFitter svfit("jjFit"+file, &trueD04);

	//the order of the following steps is slightly complicated
	//for fits to data, we want to weight simulation to match pT(D)/pT(j) from efficiency corrected data
	//this means DFitter::addEffs needs to be called BEFORE MCJets::weightMC
	//but for resampled MC, we want to perform the D0 fits to the resampled tuple produced by MCJets::setupTrainTestSamples
	//and we don't want to reweight in pT(D)/pT(j)
	//
	//So the correct order is:
	//  -setup MCJets
	//  -(for MC) run MCJets::setupTrainTestSamples
	//  -setup DFitter (datafile input now available for MC)
	//  -run DFitter::addEffs
	//  -run MCJets::weightMC

	if(dataIsMC) {
		//for MC don't set any InputTruePtWeights and don't pass a D0 data file
		//mcj.setInputs("",charmSimFile,beautySimFile,"","");
		mcj.setInputs("","charm","beauty","","");
		mcj.setPtCorrFactorsFile("ptCorrFactors.root");//"toyPtCorrFactors.root");
		//mcj.getTruth("data",&true4,&true5,&trueD04,&trueD05,&trueD0Sel4,&trueD0Sel5,&trueSV4,&trueSV5);
		//mcj.getTruth("data",&trueUnfld4,&trueUnfld5,&trueUnfldD04,&trueUnfldD05,&trueUnfldD0Sel4,&trueUnfldD0Sel5,&trueUnfldSV4,&trueUnfldSV5,true);
		if(weightSimTruePt<0) weightSimTruePt=0;
	} else {
		mcj.setInputs("","charm","beauty","",d0fit.dFileName());
		if(weightSimTruePt<0) weightSimTruePt=1;
	}
	if(ptCorrFile!="") mcj.setPtCorrFactorsFile(ptCorrFile);

	if(weightSimTruePt) {// !scaledSim) {
		mcj.setInputTruePtWeights(MCJets::jetRecoD04,jetTruePtWeights4);//d);
		mcj.setInputTruePtWeights(MCJets::jetRecoD05,jetTruePtWeights5);//d);
		mcj.setInputTruePtWeights(MCJets::jetRecoSV4,jetTruePtWeights4);//s);
		mcj.setInputTruePtWeights(MCJets::jetRecoSV5,jetTruePtWeights5);//s);
	}

	//run the D0 efficiencies first because in data we need these to weight the MC
	//if(useSimpleEff) {
	//	return 1;
	//	//d0fit.setInputs("data","charm","beauty",simpleEffFile,"",true,dataIsMC,useRhoZEffCor);
	//} else {
		d0fit.setInputs("data","charm","beauty",effFileD0,accFileD0,dataIsMC);//,false,dataIsMC,useRhoZEffCor);
	//}
	d0fit.setOptions(dfitOpts);
	if(runStages.count(stageEffs)) d0fit.setRerunEffs();
	//d0fit.setTruthMatchType(dTruthMatch);
	//if(skipSumW2Fits) d0fit.skipSumW2Fits();
	d0fit.addEffs();

	//now weight the MC
	if(runStages.count(stageWeight)) mcj.setRerunWeights();
	if(!mcj.weightMC(MCJets::jetRecoSV4,&hUE4,&hE4)) return 1;
	if(!mcj.weightMC(MCJets::jetRecoSV5,&hUE5,&hE5)) return 1;

	if(!mcj.weightMC(MCJets::jetRecoD04)) return 1;
	if(!mcj.weightMC(MCJets::jetRecoD05)) return 1;

	if(dataIsMC) {
		//mcj.getTruth(mcj.resampledName(),&true4,&true5,&trueD04,&trueD05,&trueD0Sel4,&trueD0Sel5,&trueSV4,&trueSV5);
		//mcj.getTruth(mcj.resampledName(),&trueUnfld4,&trueUnfld5,&trueUnfldD04,&trueUnfldD05,&trueUnfldD0Sel4,&trueUnfldD0Sel5,&trueUnfldSV4,&trueUnfldSV5,true);
		mcj.getTruth("data",&true4,&true5,&trueD04,&trueD05,&trueD0Sel4,&trueD0Sel5,&trueSV4,&trueSV5);
		mcj.getTruth("data",&trueUnfld4,&trueUnfld5,&trueUnfldD04,&trueUnfldD05,&trueUnfldD0Sel4,&trueUnfldD0Sel5,&trueUnfldSV4,&trueUnfldSV5,true);
		//for the real MC samples, test D0 efficiencies and return
		//for permutations, compare truth results to extracted yields
		if(!dataIsResampledMC) {
			if(file.BeginsWith("sim4")) {
				d0fit.testEffs(4);
			}
			if(file.BeginsWith("sim5")) {
				d0fit.testEffs(5);
			}
			//return 0;//TODO
		}
	}

	dm->loadDataset("charmSV", "T",std::vector<TString>{mcj.outputName(MCJets::jetRecoSV4)});
	dm->loadDataset("beautySV","T",std::vector<TString>{mcj.outputName(MCJets::jetRecoSV5)});
	if(useBackDataForLightShape) {
		svfit.setInputs("back","charmSV","beautySV","data",false);
	} else {
		svfit.setInputs("light","charmSV","beautySV","data",true);
	}
	//svfit.setInputs(lightSimFile,mcj.outputName(MCJets::jetRecoSV4),mcj.outputName(MCJets::jetRecoSV5),dataFile);
	//svfit.setSVBinning(94,600.,10000.,3);
	svfit.setOptions(svOpts);
	if(runStages.count(stageSVTemplates)) svfit.setRerunTemplates();
	svfit.makeSVFitHists();

	if(!gSystem->AccessPathName(savedir+"/histsOut.root")) {
		TFile* fin = TFile::Open(savedir+"/histsOut.root");
		TH1D* hds4 = static_cast<TH1D*>(fin->Get("recoD0Sel4"));
		TH1D* hds5 = static_cast<TH1D*>(fin->Get("recoD0Sel5"));
		TH1D* hd4 = static_cast<TH1D*>(fin->Get("recoD04_evtByEvtWeight"));
		TH1D* hd5 = static_cast<TH1D*>(fin->Get("recoD05_evtByEvtWeight"));
		TH1D* hda4 = static_cast<TH1D*>(fin->Get("recoD04_aveWeight"));
		TH1D* hda5 = static_cast<TH1D*>(fin->Get("recoD05_aveWeight"));
		TH1D* hdw4 = static_cast<TH1D*>(fin->Get("recoD0Weight4"));
		TH1D* hdw5 = static_cast<TH1D*>(fin->Get("recoD0Weight5"));
		TH1D* hs4 = static_cast<TH1D*>(fin->Get("recoSV4"));
		TH1D* hs5 = static_cast<TH1D*>(fin->Get("recoSV5"));

		if(!hds4 || !hds5 || !hd4 || !hd5 || !hda4 || !hda5 || !hdw4 || !hdw5) {
			std::cout << "Existing D0 reco histograms not found - rerunning fits" << std::endl;
			runStages.insert(stageFittingD);
		} else if(static_cast<uint>(hds4->GetNbinsX())!=npt ||
		          static_cast<uint>(hds5->GetNbinsX())!=npt ||
		          static_cast<uint>(hd4->GetNbinsX())!=npt ||
		          static_cast<uint>(hd5->GetNbinsX())!=npt ||
		          static_cast<uint>(hda4->GetNbinsX())!=npt ||
		          static_cast<uint>(hda5->GetNbinsX())!=npt ||
		          static_cast<uint>(hdw4->GetNbinsX())!=npt ||
		          static_cast<uint>(hdw5->GetNbinsX())!=npt) {
			std::cout << "Existing D0 reco histograms have wrong number of bins - rerunning fits" << std::endl;
			runStages.insert(stageFittingD);
		} else {
			for(uint i=0; i<npt; ++i) {
				recoD0Sel4.SetBinContent(i+1, hds4->GetBinContent(i+1));
				recoD0Sel4.SetBinError(i+1, hds4->GetBinError(i+1));
				recoD0Sel5.SetBinContent(i+1, hds5->GetBinContent(i+1));
				recoD0Sel5.SetBinError(i+1, hds5->GetBinError(i+1));
				recoD04_evtByEvtWeight.SetBinContent(i+1, hd4->GetBinContent(i+1));
				recoD04_evtByEvtWeight.SetBinError(i+1, hd4->GetBinError(i+1));
				recoD05_evtByEvtWeight.SetBinContent(i+1, hd5->GetBinContent(i+1));
				recoD05_evtByEvtWeight.SetBinError(i+1, hd5->GetBinError(i+1));
				recoD04_aveWeight.SetBinContent(i+1, hda4->GetBinContent(i+1));
				recoD04_aveWeight.SetBinError(i+1, hda4->GetBinError(i+1));
				recoD05_aveWeight.SetBinContent(i+1, hda5->GetBinContent(i+1));
				recoD05_aveWeight.SetBinError(i+1, hda5->GetBinError(i+1));
				recoD0Weight4.SetBinContent(i+1, hdw4->GetBinContent(i+1));
				recoD0Weight4.SetBinError(i+1, hdw4->GetBinError(i+1));
				recoD0Weight5.SetBinContent(i+1, hdw5->GetBinContent(i+1));
				recoD0Weight5.SetBinError(i+1, hdw5->GetBinError(i+1));
			}
		}

		if(!hs4 || !hs5) {
			std::cout << "Existing SV reco histograms not found - running fits" << std::endl;
			runStages.insert(stageFittingSV);
		} else if(static_cast<uint>(hs4->GetNbinsX())!=npt ||
		          static_cast<uint>(hs5->GetNbinsX())!=npt) {
			std::cout << "Existing SV reco histograms have wrong number of bins - rerunning fits" << std::endl;
			runStages.insert(stageFittingSV);
		} else {
			for(uint i=0; i<npt; ++i) {
				recoSV4.SetBinContent(i+1, hs4->GetBinContent(i+1));
				recoSV4.SetBinError(i+1, hs4->GetBinError(i+1));
				recoSV5.SetBinContent(i+1, hs5->GetBinContent(i+1));
				recoSV5.SetBinError(i+1, hs5->GetBinError(i+1));
			}
		}

		fin->Close();

	} else {
		runStages.insert(stageFittingD);
		runStages.insert(stageFittingSV);
	}
	if(runStages.count(stageFittingD)) {
		//first do D0 fits for denominators
		double yield4(0.), yield5(0.), error4(0.), error5(0.), corr4(1.), corr5(1.), aveWeight4(1.), aveWeight5(1.), aveWeightErr4(0.), aveWeightErr5(0.);
		//fit all pT bins together
		d0fit.fitD(yield4,error4,yield5,error5,-1,0,0);
		for(int i=1; i<=recoD04_evtByEvtWeight.GetNbinsX(); ++i) {
			corr4 = mcj.getPtCorrFactor(MCJets::jetRecoD04,recoD04_evtByEvtWeight.GetBinLowEdge(i),recoD04_evtByEvtWeight.GetBinLowEdge(i+1));
			corr5 = mcj.getPtCorrFactor(MCJets::jetRecoD05,recoD05_evtByEvtWeight.GetBinLowEdge(i),recoD05_evtByEvtWeight.GetBinLowEdge(i+1));
			if(d0fit.fitD(yield4,error4,yield5,error5,i-1,0,0)) {
				aveWeight4 = d0fit.getAveWeight(4);
				aveWeight5 = d0fit.getAveWeight(5);
				aveWeightErr4 = d0fit.getAveWeightError(4);
				aveWeightErr5 = d0fit.getAveWeightError(5);
				printf("eff weights: %.3f+/-%.3f %.3f+/-%.3f\n",aveWeight4,aveWeightErr4,aveWeight5,aveWeightErr5);
				recoD0Sel4.SetBinContent(i,yield4*corr4);
				recoD0Sel4.SetBinError(  i,error4*corr4);
				recoD0Sel5.SetBinContent(i,yield5*corr5);
				recoD0Sel5.SetBinError(  i,error5*corr5);
				recoD04_aveWeight.SetBinContent(i,yield4*corr4*aveWeight4);
				recoD04_aveWeight.SetBinError(  i,error4*corr4*aveWeight4);
				recoD05_aveWeight.SetBinContent(i,yield5*corr5*aveWeight5);
				recoD05_aveWeight.SetBinError(  i,error5*corr5*aveWeight5);
				recoD0Weight4.SetBinContent(i,aveWeight4);
				recoD0Weight4.SetBinError(  i,aveWeightErr4);
				recoD0Weight5.SetBinContent(i,aveWeight5);
				recoD0Weight5.SetBinError(  i,aveWeightErr5);
			} else return 1;
			if(d0fit.fitD(yield4,error4,yield5,error5,i-1,0,4)) {
				recoD04_evtByEvtWeight.SetBinContent(i,yield4*corr4);
				recoD04_evtByEvtWeight.SetBinError(  i,error4*corr4);
			} //TODO//else return 1;
			if(d0fit.fitD(yield4,error4,yield5,error5,i-1,0,5)) {
				recoD05_evtByEvtWeight.SetBinContent(i,yield5*corr5);
				recoD05_evtByEvtWeight.SetBinError(  i,error5*corr5);
			} //TODO//else return 1;
		}
	}

	if(runStages.count(stageFittingSV)) {
		//make templates using all pT bins for comparison
		double nB(0.), eB(0.), nC(0.), eC(0.), nQ(0.), eQ(0.);
		svfit.enhanceTemplatesAndFit(nB,eB,nC,eC,nQ,eQ,0,0,svEnhancePasses);

		//do SV fits for numerator
		for(int i=1; i<=recoSV4.GetNbinsX(); ++i) {
			if(!svfit.enhanceTemplatesAndFit(nB,eB,nC,eC,nQ,eQ,i,0,svEnhancePasses)) continue;
			double corrC = mcj.getPtCorrFactor(MCJets::jetRecoSV4,recoSV4.GetBinLowEdge(i),recoSV4.GetBinLowEdge(i+1));
			double corrB = mcj.getPtCorrFactor(MCJets::jetRecoSV5,recoSV5.GetBinLowEdge(i),recoSV5.GetBinLowEdge(i+1));
			recoSV4.SetBinContent(i,nC*corrC);
			recoSV4.SetBinError(  i,eC*corrC);
			recoSV5.SetBinContent(i,nB*corrB);
			recoSV5.SetBinError(  i,eB*corrB);
		}
	}

	//do unfolding
	mcj.setUnfoldingMethod(unfoldingMode,unfoldingReg);
	mcj.setUnfoldingScaleSmear(jetEnergyScale,jetEnergySmear);
	TH1D* unfoldedD04_evtByEvtWeight = mcj.unfold(&recoD04_evtByEvtWeight, MCJets::jetRecoD04);
	TH1D* unfoldedD05_evtByEvtWeight = mcj.unfold(&recoD05_evtByEvtWeight, MCJets::jetRecoD05);
	TH1D* unfoldedD04_aveWeight = mcj.unfold(&recoD04_aveWeight, MCJets::jetRecoD04);
	TH1D* unfoldedD05_aveWeight = mcj.unfold(&recoD05_aveWeight, MCJets::jetRecoD05);
	TH1D* unfoldedD0Sel4 = mcj.unfold(&recoD0Sel4, MCJets::jetRecoD04);
	TH1D* unfoldedD0Sel5 = mcj.unfold(&recoD0Sel5, MCJets::jetRecoD05);
	TH1D* unfoldedSV4 = mcj.unfold(&recoSV4, MCJets::jetRecoSV4);
	TH1D* unfoldedSV5 = mcj.unfold(&recoSV5, MCJets::jetRecoSV5);

	//correct for D0 efficiency closure
	unfoldedD04_aveWeight->Divide(&effClosure4);
	unfoldedD04_evtByEvtWeight->Divide(&effClosure4);
	unfoldedD05_aveWeight->Divide(&effClosure5);
	unfoldedD05_evtByEvtWeight->Divide(&effClosure5);
	

	if(!unfoldedD04_evtByEvtWeight || !unfoldedD05_evtByEvtWeight) return 1;
	if(!unfoldedD04_aveWeight || !unfoldedD05_aveWeight) return 1;
	if(!unfoldedD0Sel4 || !unfoldedD0Sel5) return 1;
	if(!unfoldedSV4 || !unfoldedSV5) return 1;

	unfoldedD04_evtByEvtWeight->SetName("unfoldedD04_evtByEvtWeight");
	unfoldedD05_evtByEvtWeight->SetName("unfoldedD05_evtByEvtWeight");
	unfoldedD04_aveWeight->SetName("unfoldedD04_aveWeight");
	unfoldedD05_aveWeight->SetName("unfoldedD05_aveWeight");
	unfoldedD0Sel4->SetName("unfoldedD0Sel4");
	unfoldedD0Sel5->SetName("unfoldedD0Sel5");
	unfoldedSV4->SetName("unfoldedSV4");
	unfoldedSV5->SetName("unfoldedSV5");
	unfoldedD04_evtByEvtWeight->SetTitle("");
	unfoldedD05_evtByEvtWeight->SetTitle("");
	unfoldedD04_aveWeight->SetTitle("");
	unfoldedD05_aveWeight->SetTitle("");
	unfoldedD0Sel4->SetTitle("");
	unfoldedD0Sel5->SetTitle("");
	unfoldedSV4->SetTitle("");
	unfoldedSV5->SetTitle("");

	//For plots
	std::vector<TH1D*> plotD04;
	plotD04.push_back(&trueD04);
	plotD04.push_back(&recoD04_evtByEvtWeight);
	plotD04.push_back(&trueUnfldD04);
	plotD04.push_back(unfoldedD04_evtByEvtWeight);
	std::vector<TH1D*> plotD05;
	plotD05.push_back(&trueD05);
	plotD05.push_back(&recoD05_evtByEvtWeight);
	plotD05.push_back(&trueUnfldD05);
	plotD05.push_back(unfoldedD05_evtByEvtWeight);
	std::vector<TH1D*> plotD0Sel4;
	plotD0Sel4.push_back(&trueD0Sel4);
	plotD0Sel4.push_back(&recoD0Sel4);
	plotD0Sel4.push_back(&trueUnfldD0Sel4);
	plotD0Sel4.push_back(unfoldedD0Sel4);
	std::vector<TH1D*> plotD0Sel5;
	plotD0Sel5.push_back(&trueD0Sel5);
	plotD0Sel5.push_back(&recoD0Sel5);
	plotD0Sel5.push_back(&trueUnfldD0Sel5);
	plotD0Sel5.push_back(unfoldedD0Sel5);
	std::vector<TH1D*> plotSV4;
	plotSV4.push_back(&trueSV4);
	plotSV4.push_back(&recoSV4);
	plotSV4.push_back(&trueUnfldSV4);
	plotSV4.push_back(unfoldedSV4);
	std::vector<TH1D*> plotSV5;
	plotSV5.push_back(&trueSV5);
	plotSV5.push_back(&recoSV5);
	plotSV5.push_back(&trueUnfldSV5);
	plotSV5.push_back(unfoldedSV5);
	std::vector<TH1D*> plot4;
	plot4.push_back(&true4);
	plot4.push_back(&reco4);
	plot4.push_back(&trueUnfld4);
	plot4.push_back(&unfolded4);
	std::vector<TH1D*> plot5;
	plot5.push_back(&true5);
	plot5.push_back(&reco5);
	plot5.push_back(&trueUnfld5);
	plot5.push_back(&unfolded5);

	//print results
	double bfD0 = 0.03950;
	double errBFD0 = 0.00031;
	double ffc2D0 = 0.6017;//0.532;//0.542;
	double errFFc2D0 = 0.0075;//0.024;//TMath::Sqrt(.024*.024 + .007*.007);
	double ffb2D0 = 0.577;//0.587;
	double errFFb2D0 = 0.021;//TMath::Sqrt(.021*.021 + .008*.008);

	double bfDp = 0.0938;
	double errBFDp = 0.0016;
	double ffc2Dp = 0.2409;//0.219;//0.225;
	double errFFc2Dp = 0.0067;//0.010;//TMath::Sqrt(.010*.010 + .005*.005);
	double ffb2Dp = 0.218;//0.223;
	double errFFb2Dp = 0.011;//TMath::Sqrt(.011*.011 + .005*.005);

	//correction factors on FF*BF for simulation
	//true BFs and FFs in simulated sampels don't match world averages
	//vary by jet pT bin but take rough average
	//put full correction in FF
	if(dataIsMC) {
		ffc2D0 = 0.633;//*= 1.168;
		ffb2D0 = 0.607;//*= 1.034;
	}

	double bfffErr4 = TMath::Sqrt( (errFFc2D0/ffc2D0)*(errFFc2D0/ffc2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );
	double bfffErr5 = TMath::Sqrt( (errFFb2D0/ffb2D0)*(errFFb2D0/ffb2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );

	std::cout << "MC efficiency reco/true" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << hE4.GetBinContent(i,1) << " +/- " << hE4.GetBinError(i,1) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << hE5.GetBinContent(i,1) << " +/- " << hE5.GetBinError(i,1) << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << hUE4.GetBinContent(i,1) << " +/- " << hUE4.GetBinError(i,1) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << hUE5.GetBinContent(i,1) << " +/- " << hUE5.GetBinError(i,1) << std::endl;

	//D0 results
	std::cout << "D0 reco no eff4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD0Sel4.GetBinContent(i) << " +/- " << recoD0Sel4.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD0Sel5.GetBinContent(i) << " +/- " << recoD0Sel5.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD0Sel4->GetBinContent(i) << " +/- " << unfoldedD0Sel4->GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD0Sel5->GetBinContent(i) << " +/- " << unfoldedD0Sel5->GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD0Sel4.GetBinContent(i) << " +/- " << trueD0Sel4.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD0Sel5.GetBinContent(i) << " +/- " << trueD0Sel5.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD0Sel4.GetBinContent(i) << " +/- " << trueUnfldD0Sel4.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD0Sel5.GetBinContent(i) << " +/- " << trueUnfldD0Sel5.GetBinError(i) << std::endl;
	std::cout << std::endl;

	//D0 results
	std::cout << "D0 reco eff corr4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04_evtByEvtWeight.GetBinError(i)/recoD04_evtByEvtWeight.GetBinContent(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05_evtByEvtWeight.GetBinError(i)/recoD05_evtByEvtWeight.GetBinContent(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04_evtByEvtWeight.GetBinContent(i) << " +/- " << recoD04_evtByEvtWeight.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05_evtByEvtWeight.GetBinContent(i) << " +/- " << recoD05_evtByEvtWeight.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04_evtByEvtWeight->GetBinContent(i) << " +/- " << unfoldedD04_evtByEvtWeight->GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05_evtByEvtWeight->GetBinContent(i) << " +/- " << unfoldedD05_evtByEvtWeight->GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) << " +/- " << trueD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) << " +/- " << trueD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) << " +/- " << trueUnfldD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) << " +/- " << trueUnfldD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	std::cout << "D0 reco eff ave weight" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD04_aveWeight.GetBinContent(i), recoD04_aveWeight.GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD05_aveWeight.GetBinContent(i), recoD05_aveWeight.GetBinError(i));
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",unfoldedD04_aveWeight->GetBinContent(i), unfoldedD04_aveWeight->GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",unfoldedD05_aveWeight->GetBinContent(i), unfoldedD05_aveWeight->GetBinError(i));
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "D0 reco eff weights" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.3f\\pm%6.3f$ ",recoD0Weight4.GetBinContent(i), recoD0Weight4.GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.3f\\pm%6.3f$ ",recoD0Weight5.GetBinContent(i), recoD0Weight5.GetBinError(i));
	std::cout << std::endl;
	std::cout << std::endl;

	//D0 results scaled
	for (unsigned int i=1; i<=npt; ++i) {
		reco4.SetBinContent(i, recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0));
		reco4.SetBinError(i, TMath::Power(TMath::Power(recoD04_evtByEvtWeight.GetBinError(i) / (bfD0 * ffc2D0),2) + TMath::Power(bfffErr4*recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0),2),0.5));
		reco5.SetBinContent(i, recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0));
		reco5.SetBinError(i, TMath::Power(TMath::Power(recoD05_evtByEvtWeight.GetBinError(i) / (bfD0 * ffb2D0),2) + TMath::Power(bfffErr5*recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0),2),0.5));
		unfolded4.SetBinContent(i, unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0));
		unfolded4.SetBinError(i, TMath::Power(TMath::Power(unfoldedD04_evtByEvtWeight->GetBinError(i) / (bfD0 * ffc2D0),2) + TMath::Power(bfffErr4*unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0),2),0.5));
		unfolded5.SetBinContent(i, unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0));
		unfolded5.SetBinError(i, TMath::Power(TMath::Power(unfoldedD05_evtByEvtWeight->GetBinError(i) / (bfD0 * ffb2D0),2) + TMath::Power(bfffErr5*unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0),2),0.5));
	}

	std::cout << "jet reco4/5/unfld4/5/true(fromD0)4/5/true(fromD0)unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << recoD04_evtByEvtWeight.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << recoD05_evtByEvtWeight.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << unfoldedD04_evtByEvtWeight->GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << unfoldedD05_evtByEvtWeight->GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueUnfldD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueUnfldD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << true4.GetBinContent(i) << " +/- " << true4.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << true5.GetBinContent(i) << " +/- " << true5.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfld4.GetBinContent(i) << " +/- " << trueUnfld4.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfld5.GetBinContent(i) << " +/- " << trueUnfld5.GetBinError(i) << std::endl;
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

	//make comp plots
	plotComparison("nc2D0","n(#it{c}#rightarrow#it{D}^{0})",plotD04);
	plotComparison("nb2D0","n(#it{b}#rightarrow#it{D}^{0})",plotD05);
	plotComparison("nc2D0sel","n(#it{c}#rightarrow#it{D}^{0}) selected",plotD0Sel4);
	plotComparison("nb2D0sel","n(#it{b}#rightarrow#it{D}^{0}) selected",plotD0Sel5);
	plotComparison("nc2SV","n(#it{c}#rightarrow SV)",plotSV4);
	plotComparison("nb2SV","n(#it{b}#rightarrow SV)",plotSV5);
	plotComparison("nc","n(#it{c})",plot4);
	plotComparison("nb","n(#it{b})",plot5);

	std::ofstream fout("toysResults.log",std::ofstream::app);
	fout << file << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV4.GetBinContent(i) << "\t" << trueSV4.GetBinContent(i) << "\t" << recoSV4.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV5.GetBinContent(i) << "\t" << trueSV5.GetBinContent(i) << "\t" << recoSV5.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD0Sel4.GetBinContent(i) << "\t" << trueD0Sel4.GetBinContent(i) << "\t" << recoD0Sel4.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD0Sel5.GetBinContent(i) << "\t" << trueD0Sel5.GetBinContent(i) << "\t" << recoD0Sel5.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD04_evtByEvtWeight.GetBinContent(i) << "\t" << trueD04.GetBinContent(i) << "\t" << recoD04_evtByEvtWeight.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD05_evtByEvtWeight.GetBinContent(i) << "\t" << trueD05.GetBinContent(i) << "\t" << recoD05_evtByEvtWeight.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0) << "\t" << true4.GetBinContent(i) << "\t" << recoD04_evtByEvtWeight.GetBinError(i) / (bfD0 * ffc2D0) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0) << "\t" << true5.GetBinContent(i) << "\t" << recoD05_evtByEvtWeight.GetBinError(i) / (bfD0 * ffb2D0) << "\t";
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
		if(evtByEvtWeighting) {
			ratioRec4.push_back(recoSV4.GetBinContent(i) / (recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0)));
			ratioRec5.push_back(recoSV5.GetBinContent(i) / (recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0)));
		} else {
			ratioRec4.push_back(recoSV4.GetBinContent(i) / (recoD04_aveWeight.GetBinContent(i) / (bfD0 * ffc2D0)));
			ratioRec5.push_back(recoSV5.GetBinContent(i) / (recoD05_aveWeight.GetBinContent(i) / (bfD0 * ffb2D0)));
		}
		if(evtByEvtWeighting) {
			ratioUnfld4.push_back(unfoldedSV4->GetBinContent(i) / (unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0)));
			ratioUnfld5.push_back(unfoldedSV5->GetBinContent(i) / (unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0)));
		} else {
			ratioUnfld4.push_back(unfoldedSV4->GetBinContent(i) / (unfoldedD04_aveWeight->GetBinContent(i) / (bfD0 * ffc2D0)));
			ratioUnfld5.push_back(unfoldedSV5->GetBinContent(i) / (unfoldedD05_aveWeight->GetBinContent(i) / (bfD0 * ffb2D0)));
		}
		ratioTrue4.push_back(trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrue5.push_back(trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioTrueUnfld4.push_back(trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrueUnfld5.push_back(trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		if(evtByEvtWeighting) {
			errorRec4.push_back(ratioRec4[i-1]*TMath::Sqrt(TMath::Power(recoSV4.GetBinError(i)/recoSV4.GetBinContent(i),2)+TMath::Power(recoD04_evtByEvtWeight.GetBinError(i)/recoD04_evtByEvtWeight.GetBinContent(i),2)));
			errorRec5.push_back(ratioRec5[i-1]*TMath::Sqrt(TMath::Power(recoSV5.GetBinError(i)/recoSV5.GetBinContent(i),2)+TMath::Power(recoD05_evtByEvtWeight.GetBinError(i)/recoD05_evtByEvtWeight.GetBinContent(i),2)));
			errorUnfld4.push_back(ratioUnfld4[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV4->GetBinError(i)/unfoldedSV4->GetBinContent(i),2)+TMath::Power(unfoldedD04_evtByEvtWeight->GetBinError(i)/unfoldedD04_evtByEvtWeight->GetBinContent(i),2)));
			errorUnfld5.push_back(ratioUnfld5[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV5->GetBinError(i)/unfoldedSV5->GetBinContent(i),2)+TMath::Power(unfoldedD05_evtByEvtWeight->GetBinError(i)/unfoldedD05_evtByEvtWeight->GetBinContent(i),2)));
		} else {
			errorRec4.push_back(ratioRec4[i-1]*TMath::Sqrt(TMath::Power(recoSV4.GetBinError(i)/recoSV4.GetBinContent(i),2)+TMath::Power(recoD04_aveWeight.GetBinError(i)/recoD04_aveWeight.GetBinContent(i),2)));
			errorRec5.push_back(ratioRec5[i-1]*TMath::Sqrt(TMath::Power(recoSV5.GetBinError(i)/recoSV5.GetBinContent(i),2)+TMath::Power(recoD05_aveWeight.GetBinError(i)/recoD05_aveWeight.GetBinContent(i),2)));
			errorUnfld4.push_back(ratioUnfld4[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV4->GetBinError(i)/unfoldedSV4->GetBinContent(i),2)+TMath::Power(unfoldedD04_aveWeight->GetBinError(i)/unfoldedD04_aveWeight->GetBinContent(i),2)));
			errorUnfld5.push_back(ratioUnfld5[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV5->GetBinError(i)/unfoldedSV5->GetBinContent(i),2)+TMath::Power(unfoldedD05_aveWeight->GetBinError(i)/unfoldedD05_aveWeight->GetBinContent(i),2)));
		}
	}
	std::cout << "ratios reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("%.3f+/-%.3f+/-%.3f   ", ratioRec4[i], errorRec4[i], ratioRec4[i]*bfffErr4);
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("%.3f+/-%.3f+/-%.3f   ", ratioRec5[i], errorRec5[i], ratioRec5[i]*bfffErr5);
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("%.3f+/-%.3f+/-%.3f   ", ratioUnfld4[i], errorUnfld4[i], ratioUnfld4[i]*bfffErr4);
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("%.3f+/-%.3f+/-%.3f   ", ratioUnfld5[i], errorUnfld5[i], ratioUnfld5[i]*bfffErr5);
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("%.3f                   ",  trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("%.3f                   ",  trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)));
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("%.3f                   ", trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("%.3f                   ", trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)));
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Unfolding table:" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%4.0f$ ",unfoldedD04_aveWeight->GetBinContent(i), unfoldedD04_aveWeight->GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%4.0f$ ",unfoldedD05_aveWeight->GetBinContent(i), unfoldedD05_aveWeight->GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%4.0f$ ",unfoldedSV4->GetBinContent(i), unfoldedSV4->GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%4.0f$ ",unfoldedSV5->GetBinContent(i), unfoldedSV5->GetBinError(i));
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("& $%.3f\\pm%.3f\\pm%.3f$ ", ratioUnfld4[i], errorUnfld4[i], ratioUnfld4[i]*bfffErr4);
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  printf("& $%.3f\\pm%.3f\\pm%.3f$ ", ratioUnfld5[i], errorUnfld5[i], ratioUnfld5[i]*bfffErr5);
	std::cout << std::endl;

	if(dataIsMC) {
	//ratios and pulls from true values
		if(!file.BeginsWith("sim5")) {
			std::cout << "charm reco pT bins" << std::endl;
			std::cout << "D0 no efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD0Sel4.GetBinContent(i)/trueD0Sel4.GetBinContent(i), (recoD0Sel4.GetBinContent(i)-trueD0Sel4.GetBinContent(i))/recoD0Sel4.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD04_evtByEvtWeight.GetBinContent(i)/trueD04.GetBinContent(i), (recoD04_evtByEvtWeight.GetBinContent(i)-trueD04.GetBinContent(i))/recoD04_evtByEvtWeight.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with alt efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD04_aveWeight.GetBinContent(i)/trueD04.GetBinContent(i), (recoD04_aveWeight.GetBinContent(i)-trueD04.GetBinContent(i))/recoD04_aveWeight.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with BF/FF" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD04_evtByEvtWeight.GetBinContent(i)/(bfD0 * ffc2D0)/true4.GetBinContent(i), (recoD04_evtByEvtWeight.GetBinContent(i)/(bfD0 * ffc2D0)-true4.GetBinContent(i))/TMath::Sqrt(TMath::Power(recoD04_evtByEvtWeight.GetBinError(i) / (bfD0 * ffc2D0),2.) + TMath::Power(bfffErr4*recoD04_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffc2D0),2.)));
			std::cout << std::endl;
			std::cout << "SV" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoSV4.GetBinContent(i)/trueSV4.GetBinContent(i), (recoSV4.GetBinContent(i)-trueSV4.GetBinContent(i))/recoSV4.GetBinError(i));
			std::cout << std::endl;
			std::cout << "charm true pT bins" << std::endl;
			std::cout << "D0 no efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD0Sel4->GetBinContent(i)/trueUnfldD0Sel4.GetBinContent(i), (unfoldedD0Sel4->GetBinContent(i)-trueUnfldD0Sel4.GetBinContent(i))/unfoldedD0Sel4->GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with alt efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD04_evtByEvtWeight->GetBinContent(i)/trueUnfldD04.GetBinContent(i), (unfoldedD04_evtByEvtWeight->GetBinContent(i)-trueUnfldD04.GetBinContent(i))/unfoldedD04_evtByEvtWeight->GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with BF/FF" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD04_evtByEvtWeight->GetBinContent(i)/(bfD0 * ffc2D0)/trueUnfld4.GetBinContent(i), (unfoldedD04_evtByEvtWeight->GetBinContent(i)/(bfD0 * ffc2D0)-trueUnfld4.GetBinContent(i))/TMath::Sqrt(TMath::Power(unfoldedD04_evtByEvtWeight->GetBinError(i) / (bfD0 * ffc2D0),2.) + TMath::Power(bfffErr4*unfoldedD04_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffc2D0),2.)));
			std::cout << std::endl;
			std::cout << "SV" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedSV4->GetBinContent(i)/trueUnfldSV4.GetBinContent(i), (unfoldedSV4->GetBinContent(i)-trueUnfldSV4.GetBinContent(i))/unfoldedSV4->GetBinError(i));
			std::cout << std::endl;
		}
		if(!file.BeginsWith("sim4")) {
			std::cout << "beauty reco pT bins" << std::endl;
			std::cout << "D0 no efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD0Sel5.GetBinContent(i)/trueD0Sel5.GetBinContent(i), (recoD0Sel5.GetBinContent(i)-trueD0Sel5.GetBinContent(i))/recoD0Sel5.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD05_evtByEvtWeight.GetBinContent(i)/trueD05.GetBinContent(i), (recoD05_evtByEvtWeight.GetBinContent(i)-trueD05.GetBinContent(i))/recoD05_evtByEvtWeight.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with alt efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD05_aveWeight.GetBinContent(i)/trueD05.GetBinContent(i), (recoD05_aveWeight.GetBinContent(i)-trueD05.GetBinContent(i))/recoD05_aveWeight.GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with BF/FF" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoD05_evtByEvtWeight.GetBinContent(i)/(bfD0 * ffb2D0)/true5.GetBinContent(i), (recoD05_evtByEvtWeight.GetBinContent(i)/(bfD0 * ffb2D0)-true5.GetBinContent(i))/TMath::Sqrt(TMath::Power(recoD05_evtByEvtWeight.GetBinError(i) / (bfD0 * ffb2D0),2.) + TMath::Power(bfffErr5*recoD05_evtByEvtWeight.GetBinContent(i) / (bfD0 * ffb2D0),2.)));
			std::cout << std::endl;
			std::cout << "SV" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", recoSV5.GetBinContent(i)/trueSV5.GetBinContent(i), (recoSV5.GetBinContent(i)-trueSV5.GetBinContent(i))/recoSV5.GetBinError(i));
			std::cout << std::endl;
			std::cout << "beauty true pT bins" << std::endl;
			std::cout << "D0 no efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD0Sel5->GetBinContent(i)/trueUnfldD0Sel5.GetBinContent(i), (unfoldedD0Sel5->GetBinContent(i)-trueUnfldD0Sel5.GetBinContent(i))/unfoldedD0Sel5->GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with alt efficiencies" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD05_evtByEvtWeight->GetBinContent(i)/trueUnfldD05.GetBinContent(i), (unfoldedD05_evtByEvtWeight->GetBinContent(i)-trueUnfldD05.GetBinContent(i))/unfoldedD05_evtByEvtWeight->GetBinError(i));
			std::cout << std::endl;
			std::cout << "D0 with BF/FF" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedD05_evtByEvtWeight->GetBinContent(i)/(bfD0 * ffb2D0)/trueUnfld5.GetBinContent(i), (unfoldedD05_evtByEvtWeight->GetBinContent(i)/(bfD0 * ffb2D0)-trueUnfld5.GetBinContent(i))/TMath::Sqrt(TMath::Power(unfoldedD05_evtByEvtWeight->GetBinError(i) / (bfD0 * ffb2D0),2.) + TMath::Power(bfffErr5*unfoldedD05_evtByEvtWeight->GetBinContent(i) / (bfD0 * ffb2D0),2.)));
			std::cout << std::endl;
			std::cout << "SV" << std::endl;
			for (unsigned int i=1; i<=npt; ++i) printf("% 6.2f (% 6.2f)   ", unfoldedSV5->GetBinContent(i)/trueUnfldSV5.GetBinContent(i), (unfoldedSV5->GetBinContent(i)-trueUnfldSV5.GetBinContent(i))/unfoldedSV5->GetBinError(i));
			std::cout << std::endl;
		}
	}

	//ANA D0 YIELDS TABLE 
	std::cout << std::endl;
	std::cout << "ANA TABLES" << std::endl;
	std::cout << std::endl;
	printf("$N_{\\rm prmpt}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD0Sel4.GetBinContent(i), recoD0Sel4.GetBinError(i));
	printf("\\\\\n$N_{\\rm displ}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD0Sel5.GetBinContent(i), recoD0Sel5.GetBinError(i));
	printf("\\\\\n\\midrule\n$w_{\\rm prmpt}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.3f\\pm%6.3f$ ",recoD0Weight4.GetBinContent(i), recoD0Weight4.GetBinError(i));
	printf("\\\\\n$w_{\\rm displ}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.3f\\pm%6.3f$ ",recoD0Weight5.GetBinContent(i), recoD0Weight5.GetBinError(i));
	printf("\\\\\n\\midrule\n$N_{\\rm prmpt}^{\\rm corr}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD04_aveWeight.GetBinContent(i), recoD04_aveWeight.GetBinError(i));
	printf("\\\\\n$N_{\\rm displ}^{\\rm corr}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%5.0f\\pm%3.0f$ ",recoD05_aveWeight.GetBinContent(i), recoD05_aveWeight.GetBinError(i));
	printf("\\\\\n");
	std::cout << std::endl;
	std::cout << std::endl;

	//ANA SV YIELDS TABLE
	printf("$N_{\\cquark}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%3.0f$ ",recoSV4.GetBinContent(i), recoSV4.GetBinError(i));
	printf("\\\\\n$N_{\\bquark}$ ");
	for (unsigned int i=1; i<=npt; ++i) printf("& $%6.0f\\pm%3.0f$ ",recoSV5.GetBinContent(i), recoSV5.GetBinError(i));
	printf("\\\\\n");
	std::cout << std::endl;
	std::cout << std::endl;

	//ANA naive tagging effs

	//inputs for python script
	fout.open(savedir+"/resultsOutput.log");
	for (unsigned int i=0; i<npt; ++i)  fout << ratioUnfld4[i] << " ";
	for (unsigned int i=0; i<npt; ++i)  fout << ratioUnfld5[i] << " ";
	fout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  fout << errorUnfld4[i] << " ";
	for (unsigned int i=0; i<npt; ++i)  fout << errorUnfld5[i] << " ";
	fout << std::endl;
	for (unsigned int i=0; i<npt; ++i)  fout << ratioUnfld4[i]*bfffErr4 << " ";
	for (unsigned int i=0; i<npt; ++i)  fout << ratioUnfld5[i]*bfffErr5 << " ";
	fout << std::endl;
	fout.close();

	fout.open(savedir+"/systMCStats.log");
	for (unsigned int i=1; i<=npt; ++i)  fout << recoD0Weight4.GetBinError(i)/recoD0Weight4.GetBinContent(i) << " ";
	for (unsigned int i=1; i<=npt; ++i)  fout << recoD0Weight5.GetBinError(i)/recoD0Weight5.GetBinContent(i) << " ";
	fout << std::endl;
	fout.close();
	fout.open(savedir+"/systWeighting.log");
	for (unsigned int i=1; i<=npt; ++i)  fout << (recoD04_evtByEvtWeight.GetBinContent(i)-recoD04_aveWeight.GetBinContent(i))/recoD04_aveWeight.GetBinContent(i) << " ";
	for (unsigned int i=1; i<=npt; ++i)  fout << (recoD05_evtByEvtWeight.GetBinContent(i)-recoD05_aveWeight.GetBinContent(i))/recoD05_aveWeight.GetBinContent(i) << " ";
	fout << std::endl;
	fout.close();

	//plots
	recoD04_evtByEvtWeight.SetMinimum(0.);
	recoD04_evtByEvtWeight.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04_evtByEvtWeight.GetMaximum(),recoD05_evtByEvtWeight.GetMaximum()),TMath::Max(unfoldedD04_evtByEvtWeight->GetMaximum(),unfoldedD05_evtByEvtWeight->GetMaximum())));
	recoD04_evtByEvtWeight.SetLineColor(kBlue);
	recoD05_evtByEvtWeight.SetLineColor(kRed);
	unfoldedD04_evtByEvtWeight->SetLineColor(kBlue);
	unfoldedD05_evtByEvtWeight->SetLineColor(kRed);
	unfoldedD04_evtByEvtWeight->SetLineStyle(kDashed);
	unfoldedD05_evtByEvtWeight->SetLineStyle(kDashed);
	trueD04.SetLineColor(kBlue);
	trueD05.SetLineColor(kRed);
	trueD04.SetLineStyle(kDotted);
	trueD05.SetLineStyle(kDotted);

	recoSV4.SetMinimum(0.);
	recoSV4.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04_evtByEvtWeight.GetMaximum(),recoD05_evtByEvtWeight.GetMaximum()),TMath::Max(unfoldedD04_evtByEvtWeight->GetMaximum(),unfoldedD05_evtByEvtWeight->GetMaximum())));
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
	recoD04_evtByEvtWeight.Draw();
	unfoldedD04_evtByEvtWeight->Draw("same");
	recoD05_evtByEvtWeight.Draw("same");
	unfoldedD05_evtByEvtWeight->Draw("same");
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

	// Save histograms
	TFile* fHistsOut = TFile::Open(savedir+"/histsOut.root","RECREATE");
	trueD04.Write();
	recoD04_evtByEvtWeight.Write();
	recoD04_aveWeight.Write();
	recoD0Weight4.Write();
	trueUnfldD04.Write();
	unfoldedD04_evtByEvtWeight->Write();
	unfoldedD04_aveWeight->Write();
	trueD05.Write();
	recoD05_evtByEvtWeight.Write();
	recoD05_aveWeight.Write();
	recoD0Weight5.Write();
	trueUnfldD05.Write();
	unfoldedD05_evtByEvtWeight->Write();
	unfoldedD05_aveWeight->Write();
	trueD0Sel4.Write();
	recoD0Sel4.Write();
	trueUnfldD0Sel4.Write();
	unfoldedD0Sel4->Write();
	trueD0Sel5.Write();
	recoD0Sel5.Write();
	trueUnfldD0Sel5.Write();
	unfoldedD0Sel5->Write();
	trueSV4.Write();
	recoSV4.Write();
	trueUnfldSV4.Write();
	unfoldedSV4->Write();
	trueSV5.Write();
	recoSV5.Write();
	trueUnfldSV5.Write();
	unfoldedSV5->Write();
	true4.Write();
	reco4.Write();
	trueUnfld4.Write();
	unfolded4.Write();
	true5.Write();
	reco5.Write();
	trueUnfld5.Write();
	unfolded5.Write();
	fHistsOut->Close();

	//output tagging efficiencies for use by Z+j fit
	TFile* fTagEffOut = TFile::Open(savedir+"/taggingEffs.root","RECREATE");
	TH1D tagEff4("tagEff4","",npt,ptBinBoundaries.data());
	TH1D tagEff5("tagEff5","",npt,ptBinBoundaries.data());
	TH1D hist_bfffErr4("bfffErr4","",npt,ptBinBoundaries.data());
	TH1D hist_bfffErr5("bfffErr5","",npt,ptBinBoundaries.data());
	TH1D tagEffFull4("tagEffFull4","",npt,ptBinBoundaries.data());
	TH1D tagEffFull5("tagEffFull5","",npt,ptBinBoundaries.data());
	for(uint i=0; i<npt; ++i) {
		tagEff4.SetBinContent(i+1, ratioUnfld4[i]);
		tagEff4.SetBinError(i+1, errorUnfld4[i]);
		tagEff5.SetBinContent(i+1, ratioUnfld5[i]);
		tagEff5.SetBinError(i+1, errorUnfld5[i]);
		hist_bfffErr4.SetBinContent(i+1, 1.);
		hist_bfffErr4.SetBinError(i+1, bfffErr4);
		hist_bfffErr5.SetBinContent(i+1, 1.);
		hist_bfffErr5.SetBinError(i+1, bfffErr5);
		tagEffFull4.SetBinContent(i+1, ratioUnfld4[i]);
		tagEffFull4.SetBinError(i+1, TMath::Sqrt(TMath::Power(errorUnfld4[i],2.) + TMath::Power(ratioUnfld4[i]*bfffErr4,2.)));
		tagEffFull5.SetBinContent(i+1, ratioUnfld5[i]);
		tagEffFull5.SetBinError(i+1, TMath::Sqrt(TMath::Power(errorUnfld5[i],2.) + TMath::Power(ratioUnfld5[i]*bfffErr5,2.)));
	}
	tagEff4.Write();
	tagEff5.Write();
	hist_bfffErr4.Write();
	hist_bfffErr5.Write();
	tagEffFull4.Write();
	tagEffFull5.Write();
	fTagEffOut->Close();

	return 0;
}
