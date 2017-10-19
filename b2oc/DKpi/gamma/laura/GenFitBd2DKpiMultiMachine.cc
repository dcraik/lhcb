
#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"

#include "LauConstants.hh"
#include "LauCPFitModel.hh"
#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauWeightedSumEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauMagPhaseCoeffSet.hh"
#include "LauRealImagCoeffSet.hh"
#include "LauRealImagCPCoeffSet.hh"
#include "LauCartesianCPCoeffSet.hh"
#include "LauCartesianGammaCPCoeffSet.hh"
#include "LauPolarGammaCPCoeffSet.hh"
#include "LauCartesianGammaCPCoeffSet.hh"
#include "LauPolarGammaCPCoeffSet.hh"
#include "LauRealImagGammaCPCoeffSet.hh"
#include "LauCleoCPCoeffSet.hh"
#include "LauVetoes.hh"

#include "LauKinematics.hh"

#include "LauLASSNRRes.hh"
#include "LauLASSRes.hh"
#include "LauModIndPartWaveMagPhase.hh"
#include "LauModIndPartWaveRealImag.hh"
#include "LauEFKLLMRes.hh"

void usage( std::ostream& out, const TString& progName )
{
	out<<"Usage:\n";
	out<<progName<<" gen [nExpt = 1] [firstExpt = 0]\n";
	out<<"or\n";
	out<<progName<<" fit <iFit> <port> [model = \"nominal\"] [ddec = \"Kpi\"] [nnbin = \"all\"] [blind = 0] [nExpt = 1] [firstExpt = 0]"<<std::endl;
}

int main( int argc, char** argv )
{
	// Process command-line arguments
	// Usage:
	// ./GenFitDKpi gen [nExpt = 1] [firstExpt = 0]
	// or
	// ./GenFitDKpi fit <iFit> [nExpt = 1] [firstExpt = 0]
	if ( argc < 2 ) {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

	TString command = argv[1];
	command.ToLower();
	Int_t port(0);
	Int_t iFit(0);
	Int_t nExpt(1);
	Int_t firstExpt(0);
	TString dDecay("Kpi");
	TString model("nominal");
	TString NNbin("all");
	Bool_t blind(false);
	TString host("localhost");
	if ( command == "gen" ) {
		if ( argc > 2 ) {
			nExpt = atoi( argv[2] );
			if ( argc > 3 ) {
				firstExpt = atoi( argv[3] );
			}
		}
	} else if ( command == "fit" ) {
		if ( argc < 4 ) {
			usage( std::cerr, argv[0] );
			return EXIT_FAILURE;
		}
		iFit = atoi( argv[2] );
		port = atoi( argv[3] );
		host = argv[4];
		if ( argc > 5 ) {
			model = argv[5];
			if ( argc > 6 ) {
				dDecay = argv[6];
				if ( argc > 7 ) {
					NNbin = argv[7];
					if ( argc > 8 ) {
						blind = atoi( argv[8] );
						if ( argc > 9 ) {
							nExpt = atoi( argv[9] );
							if ( argc > 10 ) {
								firstExpt = atoi( argv[10] );
							}
						}
					}
				}
			}
		}
	} else {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

	Bool_t Ds(kTRUE), DsCPV(kTRUE), LASSCPV(kFALSE), floatShapes(kFALSE), oldLASSParams(kFALSE), kappa(kFALSE), mipw(kFALSE), DsSWave(kFALSE), DsSWaveCPV(kFALSE), DsPWave(kFALSE), DsPWaveCPV(kFALSE), extraLASSCPV(kFALSE), oldEff(kFALSE), newLASSParams(kFALSE), newnewLASSParams(kFALSE), altDstKpi(kFALSE), floatMGs(kFALSE), constrainMGs(kFALSE), glass(kTRUE), fixDK(kFALSE), DsPWaveCPVAlt(kFALSE), efkllm(kFALSE), lassCutoff(kFALSE), akima(kFALSE), prodAsym(kFALSE), gamma892(kFALSE), fixedDeltas(kFALSE), singleNNShape(kFALSE), altDpiSWave(kFALSE), combAsym(kFALSE), no1410(kFALSE);//, efkllmFixPhase(kFALSE); //, constrainDK(kFALSE);
	Int_t newMGs(0);
	Int_t shapeParams(0);
	Int_t mipwV(0);

	dDecay.ToLower();
	std::cout << dDecay << std::endl;
	NNbin.ToLower();
	std::cout << NNbin << std::endl;

	switch(iFit/100) {
		case 20:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			oldLASSParams = kTRUE;
			kappa = kFALSE;
			std::cout << "running 2000-series fit..." << std::endl;
			break;
		case 21:
			LASSCPV = kFALSE;
			floatShapes = kFALSE;
			oldLASSParams = kTRUE;
			kappa = kFALSE;
			std::cout << "running 2100-series fit..." << std::endl;
			break;
		case 22:
			LASSCPV = kFALSE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kFALSE;
			std::cout << "running 2200-series fit..." << std::endl;
			break;
		case 23:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kFALSE;
			std::cout << "running 2300-series fit..." << std::endl;
			break;
		case 24:
			LASSCPV = kFALSE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kTRUE;
			std::cout << "running 2400-series fit..." << std::endl;
			break;
		case 25:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kTRUE;
			std::cout << "running 2500-series fit..." << std::endl;
			break;
		case 26:
			LASSCPV = kFALSE;
			floatShapes = kTRUE;
			mipw = kTRUE;
			std::cout << "running 2600-series fit..." << std::endl;
			break;
		case 27:
			LASSCPV = kFALSE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			std::cout << "running 2700-series fit..." << std::endl;
			break;
		case 28:
			LASSCPV = kFALSE;
			floatShapes = kFALSE;
			oldLASSParams = kFALSE;
			kappa = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kTRUE;
			std::cout << "running 2800-series fit..." << std::endl;
			break;
		case 29:
			LASSCPV = kFALSE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			std::cout << "running 2900-series fit..." << std::endl;
			break;
		case 30:
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kTRUE;
			Ds = kTRUE;
			DsCPV = kTRUE;
			std::cout << "running 3000-series fit..." << std::endl;
			break;
		case 31:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			std::cout << "running 3100-series fit..." << std::endl;
			break;
		case 32:
			kappa = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			std::cout << "running 3200-series fit..." << std::endl;
			break;
		case 33:
			LASSCPV = kTRUE;
			extraLASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			std::cout << "running 3300-series fit..." << std::endl;
			break;
		case 34:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			std::cout << "running 3400-series fit..." << std::endl;
			break;
		case 35:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			newLASSParams = kTRUE;
			std::cout << "running 3500-series fit..." << std::endl;
			break;
		case 36:
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kTRUE;
			Ds = kTRUE;
			DsCPV = kTRUE;
			newLASSParams = kTRUE;
			std::cout << "running 3600-series fit..." << std::endl;
			break;
		case 37:
			LASSCPV = kTRUE;
			extraLASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			newLASSParams = kTRUE;
			std::cout << "running 3700-series fit..." << std::endl;
			break;
		case 38:
			kappa = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			newLASSParams = kTRUE;
			std::cout << "running 3800-series fit..." << std::endl;
			break;
		case 39://3DP
		case 40://1DP 0-19 c-e, 20-39 b-e
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			oldLASSParams = kTRUE;
			std::cout << "running 3900-series fit..." << std::endl;
			break;
		case 41://3DP
		case 42://1DP 0-19 c-e, 20-39 b-e
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			oldLASSParams = kTRUE;
			std::cout << "running 4100-series fit..." << std::endl;
			break;
		case 43:
			LASSCPV = kFALSE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newLASSParams = kTRUE;
			std::cout << "running 4300-series fit..." << std::endl;
			break;
		case 44:
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newLASSParams = kTRUE;
			std::cout << "running 4400-series fit..." << std::endl;
			break;
		case 45://new baseline
		case 49://11DP
		case 50://4DP
		case 51://14DP - new DCP geom efficiencies
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 4500-series fit..." << std::endl;
			break;
		case 46://new add DK CPV
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kTRUE;
			Ds = kTRUE;
			DsCPV = kTRUE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 4600-series fit..." << std::endl;
			break;
		case 47://new extra LASS CPV
			extraLASSCPV = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 4700-series fit..." << std::endl;
			break;
		case 48://new kappa
			kappa = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 4800-series fit..." << std::endl;
			break;
		case 52://14DP - no Ds2
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 5200-series fit..." << std::endl;
			break;
		case 53://14DP - no Ds NR
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 5300-series fit..." << std::endl;
			break;
		case 54://14DP - slight modification of baseline - shapeParams=5 - go back to newnew - these seem less stable
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=5;
			std::cout << "running 5400-series fit..." << std::endl;
			break;
		case 55://14DP - add 2700
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 4500-series fit..." << std::endl;
			break;
		case 56://14DP - alternative D*Kpi bkg
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			altDstKpi = kTRUE;
			std::cout << "running 5600-series fit..." << std::endl;
			break;
		case 57://new kappa float
			kappa = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 5700-series fit..." << std::endl;
			break;
		case 58://new kappa float no CPV
			kappa = kTRUE;
			LASSCPV = kFALSE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			std::cout << "running 5800-series fit..." << std::endl;
			break;
		case 59://14DP - new masses and widths
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			newMGs = 1;
			std::cout << "running 5900-series fit..." << std::endl;
			break;
		case 60://14DP - new masses and widths
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			newMGs = 2;
			std::cout << "running 6000-series fit..." << std::endl;
			break;
		case 61://14DP - float masses and widths
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			floatMGs = kTRUE;
			std::cout << "running 6100-series fit..." << std::endl;
			break;
		case 62://14DP - float masses, widths and shapes
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			std::cout << "running 6200-series fit..." << std::endl;
			break;
		case 63://14DP - new masses, widths and shapes
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			//newnewLASSParams = kTRUE;
			shapeParams=6;
			newMGs = 3;
			std::cout << "running 6300-series fit..." << std::endl;
			break;
		case 64://14DP - kappa, float masses, widths and shapes
			kappa = kTRUE;
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			newnewLASSParams = kTRUE;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			std::cout << "running 6400-series fit..." << std::endl;
			break;
		case 65://14DP - add 2700
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			oldEff = kTRUE;
			//newnewLASSParams = kTRUE;
			shapeParams=6;
			newMGs = 3;
			std::cout << "running 6500-series fit..." << std::endl;
			break;
		case 66://14DP - no DK D-wave
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			oldEff = kTRUE;
			//newnewLASSParams = kTRUE;
			shapeParams=6;
			newMGs = 3;
			std::cout << "running 6600-series fit..." << std::endl;
			break;
		case 67://14DP - no DK S-wave or D-wave
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			oldEff = kTRUE;
			//newnewLASSParams = kTRUE;
			shapeParams=6;
			newMGs = 3;
			std::cout << "running 6700-series fit..." << std::endl;
			break;
		case 68://14DP - new masses, widths and shapes, simple LASS
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kTRUE;
			DsSWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			//newnewLASSParams = kTRUE;
			shapeParams=6;
			newMGs = 3;
			glass = kFALSE;
			std::cout << "running 6800-series fit..." << std::endl;
			break;
		case 69://14DP - float masses, widths and shapes, constr DK
		case 77://14DP - float masses, widths and shapes, constr DK, new constr
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 6900-series fit..." << std::endl;
			break;
		case 70://14DP - float masses, widths and shapes, fix DK
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			fixDK = kTRUE;
			std::cout << "running 7000-series fit..." << std::endl;
			break;
		case 71://14DP - float masses, widths and shapes, constr DK add Ds2
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7100-series fit..." << std::endl;
			break;
		case 72://14DP - float masses, widths and shapes, constr DK add Ds CPV
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7200-series fit..." << std::endl;
			break;
		case 73://14DP - float masses, widths and shapes, constr DK, sLASS
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			glass = kFALSE;
			//constrainDK = kTRUE;
			std::cout << "running 7300-series fit..." << std::endl;
			break;
		case 74://14DP - float masses, widths and shapes, constr DK add Ds CPV alt
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			DsPWaveCPVAlt = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7400-series fit..." << std::endl;
			break;
		case 75://14DP - float masses, widths and shapes, unconstr DK add Ds CPV alt
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			DsPWaveCPVAlt = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7500-series fit..." << std::endl;
			break;
		case 76://14DP - float masses, widths and shapes, fix DK add Ds CPV alt
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			DsPWaveCPVAlt = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			fixDK = kTRUE;
			std::cout << "running 7600-series fit..." << std::endl;
			break;
		case 78://14DP - float masses, widths and shapes, constrain DK, EFKLLM
			LASSCPV = kTRUE;
			efkllm = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7800-series fit..." << std::endl;
			break;
		case 79://14DP - new masses, widths and shapes, constr DK, new constr
		case 92://14DP - xChecks - 00-19 no pipi, 20-39 no KK, 40-59 no A bins
		case 93://14DP - new data tuples
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 7900-series fit..." << std::endl;
			break;
		case 80://14DP - new masses, widths and shapes, constr DK, new constr, add Ds CPV
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			DsPWaveCPVAlt = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8000-series fit..." << std::endl;
			break;
		case 81://14DP - new masses, widths and shapes, constr DK, new constr, add Ds2*
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kTRUE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8100-series fit..." << std::endl;
			break;
		case 82://14DP - new masses, widths and shapes, constr DK, new constr, sLASS
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			glass = kFALSE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8200-series fit..." << std::endl;
			break;
		case 83://14DP - new masses, widths and shapes, constr DK, new constr, efkllm
			LASSCPV = kTRUE;
			efkllm = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8300-series fit..." << std::endl;
			break;
		case 84://14DP - float masses, widths and shapes, constr DK, new constr, kappa
			LASSCPV = kTRUE;
			kappa = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=6;
			newMGs = 3;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8400-series fit..." << std::endl;
			break;
		case 85://14DP - new masses, widths and shapes, constr DK, new constr, no SWaveCPV
			LASSCPV = kFALSE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8500-series fit..." << std::endl;
			break;
		case 86://14DP - new masses, widths and shapes, constr DK, new constr, extra SWaveCPV
			LASSCPV = kTRUE;
			extraLASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8600-series fit..." << std::endl;
			break;
		case 87://14DP - new masses, widths and shapes, constr DK, new constr, sLASS cutoff
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			glass = kFALSE;
			lassCutoff = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8700-series fit..." << std::endl;
			break;
		case 88://14DP - new masses, widths and shapes, constr DK, new constr, float MIPW
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			mipw = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8800-series fit..." << std::endl;
			break;
		case 89://14DP - new masses, widths and shapes, constr DK, new constr, float MIPW Akima
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			mipw = kTRUE;
			akima = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 8900-series fit..." << std::endl;
			break;
		case 90://14DP - new masses, widths and shapes, constr DK, new constr, float MIPW V1 Akima
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			mipw = kTRUE;
			akima = kTRUE;
			mipwV = 1;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9000-series fit..." << std::endl;
			break;
		case 91://14DP - new masses, widths and shapes, constr DK, new constr, float MIPW V2 Akima
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			mipw = kTRUE;
			akima = kTRUE;
			mipwV = 2;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9100-series fit..." << std::endl;
			break;
		case 94://14DP - add production asymmetry
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			prodAsym = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9400-series fit..." << std::endl;
			break;
		case 95://14DP - gamma892 
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			gamma892 = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9400-series fit..." << std::endl;
			break;
		case 96://14DP - fixed Deltas
			LASSCPV = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			fixedDeltas = kTRUE;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9600-series fit..." << std::endl;
			break;
		case 97://14DP - new data tuples, kappa
			LASSCPV = kTRUE;
			kappa = kTRUE;
			//floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			//floatMGs = kTRUE;
			//constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9700-series fit..." << std::endl;
			break;
		case 98://14DP - new data tuples, kappa, float shapes
			LASSCPV = kTRUE;
			kappa = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9800-series fit..." << std::endl;
			break;
		case 99://14DP - new data tuples, float shapes
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 9900-series fit..." << std::endl;
			break;
		case 100://14DP - new data tuples, single NNshape
			LASSCPV = kTRUE;
//			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			singleNNShape = kTRUE;
//			floatMGs = kTRUE;
//			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10000-series fit..." << std::endl;
			break;
		case 101://14DP - new data tuples, Ds CPV
			LASSCPV = kTRUE;
//			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kTRUE;
			DsPWaveCPVAlt = kTRUE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
//			floatMGs = kTRUE;
//			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10100-series fit..." << std::endl;
			break;
		case 102://14DP - new data tuples, alt Dpi S-wave
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			altDpiSWave = kTRUE;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10200-series fit..." << std::endl;
			break;
		case 103://14DP - new data tuples, comb CPV
			LASSCPV = kTRUE;
//			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			combAsym = kTRUE;
//			floatMGs = kTRUE;
//			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10300-series fit..." << std::endl;
			break;
		case 104://14DP - new data tuples, no 1410
			LASSCPV = kTRUE;
//			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			no1410 = kTRUE;
//			floatMGs = kTRUE;
//			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10400-series fit..." << std::endl;
			break;
		case 105://14DP - new data tuples, alt Dpi S-wave
			LASSCPV = kTRUE;
			floatShapes = kFALSE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			altDpiSWave = kTRUE;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10500-series fit..." << std::endl;
			break;
		case 106://14DP - new data tuples, no 1410, floatShapes
			LASSCPV = kTRUE;
			floatShapes = kTRUE;
			DsSWave = kFALSE;
			DsSWaveCPV = kFALSE;
			DsPWave = kTRUE;
			DsPWaveCPV = kFALSE;
			Ds = kFALSE;
			DsCPV = kFALSE;
			oldEff = kTRUE;
			shapeParams=7;
			newMGs = 4;
			no1410 = kTRUE;
			floatMGs = kTRUE;
			constrainMGs = kTRUE;
			//constrainDK = kTRUE;
			std::cout << "running 10600-series fit..." << std::endl;
			break;
			//case 85://14DP - new masses, widths and shapes, constr DK, new constr, efkllm, fix phase
			//	LASSCPV = kTRUE;
			//	efkllm = kTRUE;
			//	efkllmFixPhase = kTRUE;
			//	//floatShapes = kTRUE;
			//	DsSWave = kFALSE;
			//	DsSWaveCPV = kFALSE;
			//	DsPWave = kTRUE;
			//	DsPWaveCPV = kFALSE;
			//	Ds = kFALSE;
			//	DsCPV = kFALSE;
			//	oldEff = kTRUE;
			//	shapeParams=7;
			//	newMGs = 4;
			//	//floatMGs = kTRUE;
			//	//constrainMGs = kTRUE;
			//	//constrainDK = kTRUE;
			//	std::cout << "running 8500-series fit..." << std::endl;
		default:
			std::cout << "Unknown config..." << std::endl;
			return 0;

	}

	TString kpiNRName("LASSNR0");
	if(kappa) kpiNRName="kappa0";

	// If you want to use square DP histograms for efficiency,
	// backgrounds or you just want the square DP co-ordinates
	// stored in the toy MC ntuple then set this to kTRUE
	//Bool_t squareDP = kFALSE;
	Bool_t squareDP = kTRUE;

	// This defines the DP => decay is B0 -> D0_bar pi- K+
	// Particle 1 = D0_bar
	// Particle 2 = pi-
	// Particle 3 = K+
	// The DP is defined in terms of m13Sq and m23Sq
	LauDaughters* posDaughters = new LauDaughters( "B0", "D0_bar", "pi-", "K+", squareDP );
	LauDaughters* negDaughters = new LauDaughters( "B0_bar", "D0", "pi+", "K-", squareDP );

	// Optionally apply some vetoes to the DP
	LauVetoes* vetoes = new LauVetoes();
	Double_t DstMid = 2.01028;
	Double_t DstMin = DstMid-0.0025;
	Double_t DstMax = DstMid+0.0025;
	//Double_t DstMin = 2.00028;
	//Double_t DstMax = 2.02028;
	Double_t D0Min = 1.835;
	Double_t D0Max = 1.880;
	vetoes->addMassVeto(3, DstMin, DstMax); // Dst veto, m12
	vetoes->addMassVeto(1, D0Min, D0Max); // D0 veto, m23

	// Define the efficiency model (defaults to unity everywhere)
	// Can optionally provide a histogram to model variation over DP
	// (example syntax given in commented-out section)
	LauAbsEffModel* effModel;
	Bool_t useInterpolation = kTRUE;
	Bool_t fluctuateBins = kFALSE;
	Bool_t useUpperHalf = kFALSE;

	LauEffModel* effModelTOS    = new LauEffModel(posDaughters, vetoes);
	LauEffModel* effModelNotTOS = new LauEffModel(posDaughters, vetoes);
	effModel = new LauWeightedSumEffModel(posDaughters);

	if(dDecay=="kpi" && NNbin=="old") {
		TFile *effHistFile0 = TFile::Open("eff/d2kpi/geom.root", "read");
		TFile *effHistFile1 = TFile::Open("eff/d2kpi/seltrigTOS.root", "read");
		TFile *effHistFile2 = TFile::Open("eff/d2kpi/seltrigNotTOS.root", "read");
		TFile *effHistFile3 = TFile::Open("eff/d2kpi/pid.root", "read");
		TFile *effHistFile4 = TFile::Open("eff/d2kpi/trk.root", "read");
		TFile *effHistFile5 = TFile::Open("eff/d2kpi/l0h.root", "read");
		TFile *effHistFile6 = TFile::Open("eff/d2kpi/l0noth.root", "read");
		TH2* effHist0 = dynamic_cast<TH2*>(effHistFile0->Get("efficiency"));
		TH2* effHist1 = dynamic_cast<TH2*>(effHistFile1->Get("efficiency"));
		TH2* effHist2 = dynamic_cast<TH2*>(effHistFile2->Get("efficiency"));
		TH2* effHist3 = dynamic_cast<TH2*>(effHistFile3->Get("efficiency"));
		TH2* effHist4 = dynamic_cast<TH2*>(effHistFile4->Get("all"));
		TH2* effHist5 = dynamic_cast<TH2*>(effHistFile5->Get("ratio"));
		TH2* effHist6 = dynamic_cast<TH2*>(effHistFile6->Get("ratio"));

		effModelTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist1, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist5, 0.0, 0.0, useUpperHalf, squareDP);

		effModelNotTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist2, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist6, 0.0, 0.0, useUpperHalf, squareDP);

		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelTOS,1.006);
		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelNotTOS,0.988);
	} else if(oldEff) {
		//TFile *effHistFile0 = TFile::Open("eff-DpiK/geom.root", "read");
		//TFile *effHistFile1 = TFile::Open("eff-DpiK/seltrigTOS.root", "read");
		//TFile *effHistFile2 = TFile::Open("eff-DpiK/seltrigNotTOS.root", "read");
		//TFile *effHistFile3 = TFile::Open("eff-DpiK/pid.root", "read");
		//TFile *effHistFile4 = TFile::Open("eff-DpiK/trk.root", "read");
		//TFile *effHistFile5 = TFile::Open("eff-DpiK/l0h.root", "read");
		//TFile *effHistFile6 = TFile::Open("eff-DpiK/l0noth.root", "read");
		TFile *effHistFile0 = TFile::Open("eff/d2"+dDecay+"/geom.root", "read");
		TFile *effHistFile1 = TFile::Open("eff/d2"+dDecay+"/seltrigTOS.root", "read");
		TFile *effHistFile2 = TFile::Open("eff/d2"+dDecay+"/seltrigNotTOS.root", "read");
		TFile *effHistFile3 = TFile::Open("eff/d2"+dDecay+"/pid.root", "read");
		TFile *effHistFile4 = TFile::Open("eff/d2"+dDecay+"/trk.root", "read");
		TFile *effHistFile5 = TFile::Open("eff/d2"+dDecay+"/l0h.root", "read");
		TFile *effHistFile6 = TFile::Open("eff/d2"+dDecay+"/l0noth.root", "read");
		TH2* effHist0 = dynamic_cast<TH2*>(effHistFile0->Get("efficiency"));
		TH2* effHist1 = dynamic_cast<TH2*>(effHistFile1->Get("efficiency"));
		TH2* effHist2 = dynamic_cast<TH2*>(effHistFile2->Get("efficiency"));
		TH2* effHist3 = dynamic_cast<TH2*>(effHistFile3->Get("efficiency"));
		TH2* effHist4 = dynamic_cast<TH2*>(effHistFile4->Get("all"));
		TH2* effHist5 = dynamic_cast<TH2*>(effHistFile5->Get("ratio"));
		TH2* effHist6 = dynamic_cast<TH2*>(effHistFile6->Get("ratio"));

		effModelTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist1, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist5, 0.0, 0.0, useUpperHalf, squareDP);

		effModelNotTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist2, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist6, 0.0, 0.0, useUpperHalf, squareDP);

		//OLD//		all     TIS     TOS     TISandTOS       TISnotTOS
		//OLD//		3790    2133    2508    851     	1282
		//		all     TIS     TOS     TISandTOS       TISnotTOS
		//		3108    1738    2085    715     	1023
		//
		//MC		all		TOS			!TOS
		//		42346		28240			14106
		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelTOS,1.006);
		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelNotTOS,0.988);
	} else {
		TFile *effHistFile0 = TFile::Open("eff/d2"+dDecay+"/geom.root", "read");
		TFile *effHistFile1 = TFile::Open("eff/d2"+dDecay+"/seltrigTOS_"+NNbin+".root", "read");
		TFile *effHistFile2 = TFile::Open("eff/d2"+dDecay+"/seltrigNotTOS_"+NNbin+".root", "read");
		TFile *effHistFile3 = TFile::Open("eff/d2"+dDecay+"/pid_"+NNbin+".root", "read");
		TFile *effHistFile4 = TFile::Open("eff/d2"+dDecay+"/trk.root", "read");
		TFile *effHistFile5 = TFile::Open("eff/d2"+dDecay+"/l0h.root", "read");
		TFile *effHistFile6 = TFile::Open("eff/d2"+dDecay+"/l0noth.root", "read");
		//TFile *effHistFile = TFile::Open("Bs2DKpi_eff_all_combined_50.root", "read");
		//TFile *effHistFile = TFile::Open("toysCombined/eff19.root", "read");
		TH2* effHist0 = dynamic_cast<TH2*>(effHistFile0->Get("efficiency"));
		TH2* effHist1 = dynamic_cast<TH2*>(effHistFile1->Get("efficiency"));
		TH2* effHist2 = dynamic_cast<TH2*>(effHistFile2->Get("efficiency"));
		TH2* effHist3 = dynamic_cast<TH2*>(effHistFile3->Get("efficiency"));
		TH2* effHist4 = dynamic_cast<TH2*>(effHistFile4->Get("all"));
		TH2* effHist5 = dynamic_cast<TH2*>(effHistFile5->Get("ratio"));
		TH2* effHist6 = dynamic_cast<TH2*>(effHistFile6->Get("ratio"));

		//	effHist0->Print();
		//	effHist1->Print();
		//	effHist2->Print();
		//	effHist3->Print();
		//	effHist4->Print();
		//	effHist5->Print();
		//	effHist6->Print();

		effModelTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist1, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelTOS->addEffSpline(effHist5, 0.0, 0.0, useUpperHalf, squareDP);

		effModelNotTOS->setEffSpline(effHist0, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist2, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist3, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist4, 0.0, 0.0, useUpperHalf, squareDP);
		effModelNotTOS->addEffSpline(effHist6, 0.0, 0.0, useUpperHalf, squareDP);
		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelTOS,0.992);
		dynamic_cast<LauWeightedSumEffModel*>(effModel)->addEffModel(effModelNotTOS,1.015);
	}

	if(efkllm) LauEFKLLMRes::setupFormFactor("usfactor.dat");

	LauAbsResonance* res(0);

	// Create the isobar model
	LauIsobarDynamics* posSigModel = new LauIsobarDynamics(posDaughters, effModel);
	res = posSigModel->addResonance("K*0(892)",     1, LauAbsResonance::RelBW);
	res->changeResonance(0.89581, 0.0474,-1);
	if(!no1410) res = posSigModel->addResonance("K*0(1410)",    1, LauAbsResonance::RelBW);
	if(kappa) {
		res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::Flatte);
		res->changeResonance(1.513, -1, -1);
		res->setResonanceParameter("g1", 0.304);
		res->setResonanceParameter("g2", 0.380);


		res = posSigModel->addResonance("kappa0",       1, LauAbsResonance::Kappa);
		res->setResonanceParameter("b1", 5.0);
		res->setResonanceParameter("b2", 0.0);
		res->setResonanceParameter("A",  4.0);
		res->setResonanceParameter("m0", 4.2);
		if(floatShapes) {
			res->floatResonanceParameter("A");
			res->floatResonanceParameter("m0");
		}
	} else if(efkllm) {
		res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::EFKLLM);
		res->setResonanceParameter("massFactor", -2);
		res = posSigModel->addResonance(kpiNRName,      1, LauAbsResonance::EFKLLM);
		res->setResonanceParameter("massFactor",  0);
	} else if(mipw) {
		//res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::MIPW_RealImag);
		//std::set<Double_t> knots;
		//knots.insert(0.8);
		//knots.insert(1.2);
		//knots.insert(1.4);
		//knots.insert(1.6);
		//knots.insert(1.9);
		//knots.insert(3.2);
		//((LauModIndPartWaveRealImag*)res)->defineKnots(knots);
		//if(floatShapes) {
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(0, 0.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(1, 0.3,         -0.6,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(2, 0.5,         -0.4,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(3, 1.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(4, 0.4,          0.4,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(5, 0.2,          0.2,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(6, 0.1,          0.24,       kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(7, 0.0,          0.0,        kTRUE, kTRUE);
		//} else {
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(0, 0.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(1, 0.3,         -0.6,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(2, 0.5,         -0.4,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(3, 1.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(4, 0.4,          0.4,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(5, 0.2,          0.2,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(6, 0.1,          0.24,       kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(7, 0.0,          0.0,        kTRUE, kTRUE);
		//}
		res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::MIPW_MagPhase);
		if(akima) {
			dynamic_cast<LauModIndPartWaveMagPhase*>(res)->setType(Lau1DCubicSpline::AkimaSpline,Lau1DCubicSpline::AkimaSpline);
		}
		std::set<Double_t> knots;
		Double_t _a0(0.), _d0(0.);
		Double_t _a1(0.), _d1(0.);
		Double_t _a2(0.), _d2(0.);
		Double_t _a3(0.), _d3(0.);
		Double_t _a4(0.), _d4(0.);
		Double_t _a5(0.), _d5(0.);
		Double_t _a6(0.), _d6(0.);
		Double_t _a7(0.), _d7(0.);
		Double_t _a8(0.), _d8(0.);
		Double_t _a9(0.), _d9(0.);
		switch(mipwV) {
			case 0:
			default:
				knots.insert(0.8);
				knots.insert(1.0);
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(2.4);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 0.90; _d1 = 2.00;//0.8
				_a2 = 0.90; _d2 = 2.50;//1.0
				_a3 = 0.94; _d3 = 3.00;//1.2
				_a4 = 1.00; _d4 = 3.70;//1.4
				_a5 = 0.40; _d5 = 6.50;//1.7
				_a6 = 0.43; _d6 = 7.80;//1.9
				_a7 = 0.45; _d7 = 8.20;//2.4
				_a8 = 0.50; _d8 = 8.40;//3.2
				_a9 = 0.50; _d9 = 8.60;
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(7, _a7, _d7, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(8, _a8, _d8, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(9, _a9, _d9, kTRUE, kTRUE);
				break;
			case 1:
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 0.94; _d1 = 3.00;//0.8
				_a2 = 1.00; _d2 = 3.70;//1.0
				_a3 = 0.40; _d3 = 6.50;//1.2
				_a4 = 0.43; _d4 = 7.80;//1.4
				_a5 = 0.50; _d5 = 8.40;//1.7
				_a6 = 0.50; _d6 = 8.60;//1.9
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kTRUE, kTRUE);
				break;
			case 2:
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 9.84348e-01; _d1 = 3.81292e+00;
				_a2 = 1.00; _d2 = 3.70;//1.0
				_a3 = 5.29594e-01; _d3 = 5.66127e+00;
				_a4 = 3.76638e-01; _d4 = 5.04171e+00;
				_a5 = 1.73238e-01; _d5 = 4.79568e+00;
				_a6 = 0.50; _d6 = 8.60;//1.9
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kTRUE, kTRUE);
				break;
		}

	} else if(!glass) {
		res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::LASS);
		if(shapeParams==7) {
			res->changeResonance(1.479, 0.321, -1);
			res->setResonanceParameter("a", 3.6);
			res->setResonanceParameter("r", 0.0);
		} else if(shapeParams==6) {
			res->changeResonance(1.487, 0.431, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.0);//0.0);
		} else if(shapeParams==5) {
			res->changeResonance(1.446, 0.411, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.6);//4.9);
			res->setResonanceParameter("r", 1.3);//0.0);
		} else if(newnewLASSParams) {
			res->changeResonance(1.452, 0.418, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.2);//0.0);
		} else if(newLASSParams) {
			res->changeResonance(1.445, 0.391, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.0);//0.0);
		} else if(oldLASSParams) {
			res->changeResonance(1.451, 0.404, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.9);//0.0);
		} else {
			res->changeResonance(1.471, 0.241, -1); //1.467, 0.233, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 4.40); //4.48);//4.9);
			res->setResonanceParameter("r", 0.00); //0.06);//0.0);
		}
		if(floatShapes) {
			res->fixMass(kFALSE);
			res->fixWidth(kFALSE);
			res->floatResonanceParameter("a");
			res->floatResonanceParameter("r");
		}
		if(!lassCutoff) ((LauLASSRes*)res)->setCutOff(10.);
	} else {
		res = posSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::LASS_BW);
		if(shapeParams==7) {
			res->changeResonance(1.479, 0.321, -1);
			res->setResonanceParameter("a", 3.6);
			res->setResonanceParameter("r", 0.0);
		} else if(shapeParams==6) {
			res->changeResonance(1.487, 0.431, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.0);//0.0);
		} else if(shapeParams==5) {
			res->changeResonance(1.446, 0.411, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.6);//4.9);
			res->setResonanceParameter("r", 1.3);//0.0);
		} else if(newnewLASSParams) {
			res->changeResonance(1.452, 0.418, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.2);//0.0);
		} else if(newLASSParams) {
			res->changeResonance(1.445, 0.391, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.0);//0.0);
		} else if(oldLASSParams) {
			res->changeResonance(1.451, 0.404, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.9);//0.0);
		} else {
			res->changeResonance(1.471, 0.241, -1); //1.467, 0.233, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 4.40); //4.48);//4.9);
			res->setResonanceParameter("r", 0.00); //0.06);//0.0);
		}
		if(floatShapes) {
			res->fixMass(kFALSE);
			res->fixWidth(kFALSE);
			res->floatResonanceParameter("a");
			res->floatResonanceParameter("r");
		}
		res = posSigModel->addResonance(kpiNRName,      1, LauAbsResonance::LASS_NR);
		if(!lassCutoff) ((LauLASSNRRes*)res)->setCutOff(10.);
	}
	res = posSigModel->addResonance("K*0_2(1430)",  1, LauAbsResonance::RelBW);
	res = posSigModel->addResonance("D*-_0",        3, LauAbsResonance::RelBW);
	res->changeResonance(2.403,  0.283, -1);
	if(newMGs==1) res->changeResonance(2.349, 0.217, -1);
	if(newMGs==2) res->changeResonance(2.360, 0.255, -1);
	if(newMGs==3) res->changeResonance(2.353, 0.241, -1);
	if(newMGs==4) res->changeResonance(2.353, 0.248, -1);
	if(floatMGs) {
		res->fixMass(kFALSE);
		res->fixWidth(kFALSE);
		if(constrainMGs) {
			res->getMassPar()->addGaussianConstraint(2.3510,  0.0068);
			res->getWidthPar()->addGaussianConstraint(0.2302, 0.0158);
		}
	}
	res = posSigModel->addResonance("D*-_2",        3, LauAbsResonance::RelBW);
	res->changeResonance(2.4643,  0.037, -1);
	if(newMGs==1) res->changeResonance(2.4686, 0.0473, -1);
	if(newMGs==2) res->changeResonance(2.4656, 0.0460, -1);
	if(newMGs==3) res->changeResonance(2.4655, 0.0469, -1);
	if(newMGs==4) res->changeResonance(2.4656, 0.0471, -1);
	if(floatMGs) {
		res->fixMass(kFALSE);
		res->fixWidth(kFALSE);
		if(constrainMGs) {
			res->getMassPar()->addGaussianConstraint(2.4654,  0.0004);
			res->getWidthPar()->addGaussianConstraint(0.0471,   0.0012);
		}
	}
	//res = posSigModel->addResonance("D-(2760)",     3, LauAbsResonance::RelBW);
	//res->changeResonance(2.798, 0.105, 3);
	if(altDpiSWave) {
		res = posSigModel->addResonance("BelleNR_Swave-",     3, LauAbsResonance::BelleNR);
		res->setResonanceParameter("alpha", 0.48);
		res->floatResonanceParameter("alpha");
	} else {
		res = posSigModel->addResonance("dabba-",     3, LauAbsResonance::Dabba);
	}
	res = posSigModel->addResonance("BelleNR_Pwave-",     3, LauAbsResonance::BelleNR);
	if(shapeParams==7) {
		res->setResonanceParameter("alpha", 0.50);
	} else if(shapeParams==6) {
		res->setResonanceParameter("alpha", 0.48);
	} else if(shapeParams==5) {
		res->setResonanceParameter("alpha", 0.52);
	} else if(newnewLASSParams) {
		res->setResonanceParameter("alpha", 0.45); //0.39);
	} else if(newLASSParams) {
		res->setResonanceParameter("alpha", 0.55); //0.39);
	} else if(oldLASSParams) {
		res->setResonanceParameter("alpha", 0.88); //0.39);
	} else {
		res->setResonanceParameter("alpha", 0.40);//0.42); //0.39);
	}
	if(floatShapes) {
		res->floatResonanceParameter("alpha");
	}
	if(Ds) {
		res = posSigModel->addResonance("Ds*+_2(2573)",   2, LauAbsResonance::RelBW);
		if(newMGs==1) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==2) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==3) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==4) res->changeResonance(2.5684, 0.0169, -1);
	}
	if(DsSWave) {
		res = posSigModel->addResonance("BelleNR+",  2, LauAbsResonance::BelleNR);
		if(shapeParams==7) {
			res->setResonanceParameter("alpha", 0.14);
		} else if(shapeParams==6) {
			res->setResonanceParameter("alpha", 0.14);
		} else if(shapeParams==5) {
			res->setResonanceParameter("alpha", 0.24);
		} else if(newnewLASSParams) {
			res->setResonanceParameter("alpha", 0.24);
		} else if(newLASSParams) {
			res->setResonanceParameter("alpha", 0.091);
		} else if(oldLASSParams) {
			res->setResonanceParameter("alpha", 0.412);
		} else {
			res->setResonanceParameter("alpha", 0.12);//0.412);
		}
		if(floatShapes) {
			res->floatResonanceParameter("alpha");
		}
	}
	if(DsPWave) res = posSigModel->addResonance("Ds*+_1(2700)",  2, LauAbsResonance::RelBW);
	res = posSigModel->addIncoherentResonance("D*-",  3, LauAbsResonance::GaussIncoh);
	res->changeResonance(2.01026,  0.0049, -1);
	//if(floatShapes) res->fixWidth(kFALSE);

	// Reset the maximum signal DP ASq value
	// This will be automatically adjusted to avoid bias or extreme
	// inefficiency if you get the value wrong but best to set this by
	// hand once you've found the right value through some trial and
	// error.
	posSigModel->setASqMaxValue(0.055);  
	posSigModel->setIntegralBinWidths(0.002,0.002);

	posSigModel->setIntegralBinningFactor(10);
	posSigModel->setNarrowResonanceThreshold(0.001);

	LauIsobarDynamics* negSigModel = new LauIsobarDynamics(negDaughters, effModel);
	res = negSigModel->addResonance("K*0(892)",     1, LauAbsResonance::RelBW);
	res->changeResonance(0.89581, 0.0474,-1);
	if(!no1410) res = negSigModel->addResonance("K*0(1410)",    1, LauAbsResonance::RelBW);
	if(kappa) {
		res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::Flatte);
		res->changeResonance(1.513, -1, -1);
		res->setResonanceParameter("g1", 0.304);
		res->setResonanceParameter("g2", 0.380);

		res = negSigModel->addResonance("kappa0",       1, LauAbsResonance::Kappa);
		res->setResonanceParameter("b1", 5.0);
		res->setResonanceParameter("b2", 0.0);
		res->setResonanceParameter("A",  4.0);
		res->setResonanceParameter("m0", 4.2);
	} else if(efkllm) {
		res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::EFKLLM);
		res->setResonanceParameter("massFactor", -2);
		res = negSigModel->addResonance(kpiNRName,      1, LauAbsResonance::EFKLLM);
		res->setResonanceParameter("massFactor",  0);
	} else if(mipw) {
		//res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::MIPW_RealImag);
		//std::set<Double_t> knots;
		//knots.insert(0.8);
		//knots.insert(1.2);
		//knots.insert(1.4);
		//knots.insert(1.6);
		//knots.insert(1.9);
		//knots.insert(3.2);
		//((LauModIndPartWaveRealImag*)res)->defineKnots(knots);
		//if(floatShapes) {
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(0, 0.0,         -0.1,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(1, 0.3,         -0.6,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(2, 0.5,         -0.4,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(3, 1.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(4, 0.4,          0.4,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(5, 0.2,          0.2,        kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(6, 0.1,          0.24,       kFALSE,kFALSE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(7, 0.0,          0.1,        kTRUE, kTRUE);
		//} else {
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(0, 0.0,         -0.1,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(1, 0.3,         -0.6,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(2, 0.5,         -0.4,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(3, 1.0,          0.0,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(4, 0.4,          0.4,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(5, 0.2,          0.2,        kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(6, 0.1,          0.24,       kTRUE, kTRUE);
		//	((LauModIndPartWaveRealImag*)res)->setKnotAmp(7, 0.0,          0.1,        kTRUE, kTRUE);
		//}
		res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::MIPW_MagPhase);
		if(akima) {
			dynamic_cast<LauModIndPartWaveMagPhase*>(res)->setType(Lau1DCubicSpline::AkimaSpline,Lau1DCubicSpline::AkimaSpline);
		}
		std::set<Double_t> knots;
		Double_t _a0(0.), _d0(0.);
		Double_t _a1(0.), _d1(0.);
		Double_t _a2(0.), _d2(0.);
		Double_t _a3(0.), _d3(0.);
		Double_t _a4(0.), _d4(0.);
		Double_t _a5(0.), _d5(0.);
		Double_t _a6(0.), _d6(0.);
		Double_t _a7(0.), _d7(0.);
		Double_t _a8(0.), _d8(0.);
		Double_t _a9(0.), _d9(0.);
		switch(mipwV) {
			case 0:
			default:
				knots.insert(0.8);
				knots.insert(1.0);
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(2.4);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 0.90; _d1 = 2.00;//0.8
				_a2 = 0.90; _d2 = 2.50;//1.0
				_a3 = 0.94; _d3 = 3.00;//1.2
				_a4 = 1.00; _d4 = 3.70;//1.4
				_a5 = 0.40; _d5 = 6.50;//1.7
				_a6 = 0.43; _d6 = 7.80;//1.9
				_a7 = 0.45; _d7 = 8.20;//2.4
				_a8 = 0.50; _d8 = 8.40;//3.2
				_a9 = 0.50; _d9 = 8.60;
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(7, _a7, _d7, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(8, _a8, _d8, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(9, _a9, _d9, kTRUE, kTRUE);
				break;
			case 1:
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 0.94; _d1 = 3.00;//0.8
				_a2 = 1.00; _d2 = 3.70;//1.0
				_a3 = 0.40; _d3 = 6.50;//1.2
				_a4 = 0.43; _d4 = 7.80;//1.4
				_a5 = 0.50; _d5 = 8.40;//1.7
				_a6 = 0.50; _d6 = 8.60;//1.9
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kFALSE,kFALSE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kTRUE, kTRUE);
				break;
			case 2:
				knots.insert(1.2);
				knots.insert(1.4);
				knots.insert(1.7);
				knots.insert(1.9);
				knots.insert(3.2);
				_a0 = 0.94; _d0 = 1.20;
				_a1 = 9.84348e-01; _d1 = 3.81292e+00;
				_a2 = 1.00; _d2 = 3.70;//1.0
				_a3 = 5.29594e-01; _d3 = 5.66127e+00;
				_a4 = 3.76638e-01; _d4 = 5.04171e+00;
				_a5 = 1.73238e-01; _d5 = 4.79568e+00;
				_a6 = 0.50; _d6 = 8.60;//1.9
				((LauModIndPartWaveMagPhase*)res)->defineKnots(knots);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(0, _a0, _d0, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(1, _a1, _d1, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(2, _a2, _d2, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(3, _a3, _d3, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(4, _a4, _d4, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(5, _a5, _d5, kTRUE, kTRUE);
				((LauModIndPartWaveMagPhase*)res)->setKnotAmp(6, _a6, _d6, kTRUE, kTRUE);
				break;
		}

	} else if(!glass) {
		res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::LASS);
		if(shapeParams==7) {
			res->changeResonance(1.479, 0.321, -1);
			res->setResonanceParameter("a", 3.6);
			res->setResonanceParameter("r", 0.0);
		} else if(shapeParams==6) {
			res->changeResonance(1.487, 0.431, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.0);//0.0);
		} else if(shapeParams==5) {
			res->changeResonance(1.446, 0.411, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.6);//4.9);
			res->setResonanceParameter("r", 1.3);//0.0);
		} else if(newnewLASSParams) {
			res->changeResonance(1.452, 0.418, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.2);//0.0);
		} else if(newLASSParams) {
			res->changeResonance(1.445, 0.391, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.0);//0.0);
		} else if(oldLASSParams) {
			res->changeResonance(1.451, 0.404, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.9);//0.0);
		} else {
			res->changeResonance(1.467, 0.233, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 4.48);//4.9);
			res->setResonanceParameter("r", 0.06);//0.0);
		}
		if(!lassCutoff) ((LauLASSRes*)res)->setCutOff(10.);
	} else {
		res = negSigModel->addResonance("K*0_0(1430)",  1, LauAbsResonance::LASS_BW);
		if(shapeParams==7) {
			res->changeResonance(1.479, 0.321, -1);
			res->setResonanceParameter("a", 3.6);
			res->setResonanceParameter("r", 0.0);
		} else if(shapeParams==6) {
			res->changeResonance(1.487, 0.431, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.0);//0.0);
		} else if(shapeParams==5) {
			res->changeResonance(1.446, 0.411, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.6);//4.9);
			res->setResonanceParameter("r", 1.3);//0.0);
		} else if(newnewLASSParams) {
			res->changeResonance(1.452, 0.418, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.2);//0.0);
		} else if(newLASSParams) {
			res->changeResonance(1.445, 0.391, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 2.7);//4.9);
			res->setResonanceParameter("r", 1.0);//0.0);
		} else if(oldLASSParams) {
			res->changeResonance(1.451, 0.404, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 3.2);//4.9);
			res->setResonanceParameter("r", 0.9);//0.0);
		} else {
			res->changeResonance(1.467, 0.233, -1); //1.552, 0.195,-1);
			res->setResonanceParameter("a", 4.48);//4.9);
			res->setResonanceParameter("r", 0.06);//0.0);
		}
		res = negSigModel->addResonance(kpiNRName,      1, LauAbsResonance::LASS_NR);
		if(!lassCutoff) ((LauLASSNRRes*)res)->setCutOff(10.);
	}
	res = negSigModel->addResonance("K*0_2(1430)",  1, LauAbsResonance::RelBW);
	res = negSigModel->addResonance("D*+_0",        3, LauAbsResonance::RelBW);
	res->changeResonance(2.403,  0.283, -1);
	if(newMGs==1) res->changeResonance(2.349, 0.217, -1);
	if(newMGs==2) res->changeResonance(2.360, 0.255, -1);
	if(newMGs==3) res->changeResonance(2.353, 0.241, -1);
	if(newMGs==4) res->changeResonance(2.353, 0.248, -1);
	res = negSigModel->addResonance("D*+_2",        3, LauAbsResonance::RelBW);
	res->changeResonance(2.4643,  0.037, -1);
	if(newMGs==1) res->changeResonance(2.4686, 0.0473, -1);
	if(newMGs==2) res->changeResonance(2.4656, 0.0460, -1);
	if(newMGs==3) res->changeResonance(2.4655, 0.0469, -1);
	if(newMGs==4) res->changeResonance(2.4656, 0.0471, -1);
	//res = negSigModel->addResonance("D+(2760)",     3, LauAbsResonance::RelBW);
	//res->changeResonance(2.798, 0.105, 3);
	if(altDpiSWave) {
		res = negSigModel->addResonance("BelleNR_Swave+",     3, LauAbsResonance::BelleNR);
		res->setResonanceParameter("alpha", 0.48);
		res->floatResonanceParameter("alpha");
	} else {
		res = negSigModel->addResonance("dabba+",     3, LauAbsResonance::Dabba);
	}
	res = negSigModel->addResonance("BelleNR_Pwave+",     3, LauAbsResonance::BelleNR);
	if(shapeParams==7) {
		res->setResonanceParameter("alpha", 0.50);
	} else if(shapeParams==6) {
		res->setResonanceParameter("alpha", 0.48);
	} else if(shapeParams==5) {
		res->setResonanceParameter("alpha", 0.52);
	} else if(newnewLASSParams) {
		res->setResonanceParameter("alpha", 0.45); //0.39);
	} else if(newLASSParams) {
		res->setResonanceParameter("alpha", 0.55); //0.39);
	} else if(oldLASSParams) {
		res->setResonanceParameter("alpha", 0.88); //0.39);
	} else {
		res->setResonanceParameter("alpha", 0.42); //0.39);
	}
	if(Ds) {
		res = negSigModel->addResonance("Ds*-_2(2573)",   2, LauAbsResonance::RelBW);
		if(newMGs==1) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==2) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==3) res->changeResonance(2.5684, 0.0169, -1);
		if(newMGs==4) res->changeResonance(2.5684, 0.0169, -1);
	}
	if(DsSWave) {
		res = negSigModel->addResonance("BelleNR-",  2, LauAbsResonance::BelleNR);
		if(shapeParams==7) {
			res->setResonanceParameter("alpha", 0.14);
		} else if(shapeParams==6) {
			res->setResonanceParameter("alpha", 0.14);
		} else if(shapeParams==5) {
			res->setResonanceParameter("alpha", 0.24);
		} else if(newnewLASSParams) {
			res->setResonanceParameter("alpha", 0.24);
		} else if(newLASSParams) {
			res->setResonanceParameter("alpha", 0.091);
		} else if(oldLASSParams) {
			res->setResonanceParameter("alpha", 0.412);
		} else {
			res->setResonanceParameter("alpha", 0.412);
		}
	}
	if(DsPWave) res = negSigModel->addResonance("Ds*-_1(2700)",  2, LauAbsResonance::RelBW);
	res->changeResonance(2.7082,  0.1202, -1);
	res = negSigModel->addIncoherentResonance("D*+",  3, LauAbsResonance::GaussIncoh);
	res->changeResonance(2.01026,  0.0049, -1);

	// Reset the maximum signal DP ASq value
	// This will be automatically adjusted to avoid bias or extreme
	// inefficiency if you get the value wrong but best to set this by
	// hand once you've found the right value through some trial and
	// error.
	negSigModel->setASqMaxValue(5.);  
	negSigModel->setIntegralBinWidths(0.002,0.002);

	negSigModel->setIntegralBinningFactor(10);
	negSigModel->setNarrowResonanceThreshold(0.001);

	// Create the fit model
	LauCPFitModel* fitModel = new LauCPFitModel(posSigModel,negSigModel);

	// Create the complex coefficients for the isobar model
	// Use Real and Imaginary parts not Mag and Phase
	// Set rough starting values 
	std::vector<LauAbsCoeffSet*> coeffset;
	LauAbsCoeffSet* coeff(0);
	if(dDecay=="kpi") {
		if(model=="gamma" || gamma892) {
			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0(892)",    0.03, -1.32, 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE) );
		} else if(model=="gamma2") {
			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0(892)",    0.03, -1.32, 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE) );
		} else if(model=="alt") {
			coeffset.push_back( new LauCartesianGammaCPCoeffSet("K*0(892)",    0.03, -1.32, 0., 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
		} else {
			coeffset.push_back( new LauRealImagGammaCPCoeffSet("K*0(892)",    0.03, -1.32,  0., 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
		}
		if(!no1410) coeffset.push_back( new LauRealImagCoeffSet("K*0(1410)",  -0.24, -0.02,  kFALSE, kFALSE) );
		//coeff = new LauRealImagCoeffSet("K*0_0(1430)",-0.43,  0.54,  kFALSE, kFALSE);
		//coeffset.push_back( coeff );
		coeffset.push_back( new LauRealImagCoeffSet("K*0_0(1430)",-0.43,  0.54,  kFALSE, kFALSE) );
		if(!mipw && glass) {
			//if(efkllmFixPhase) {
			//	coeffset.push_back(coeff->createClone(kpiNRName, LauAbsCoeffSet::TiePhase)); 
			//} else {
			coeffset.push_back( new LauRealImagCoeffSet(kpiNRName,    -0.40,  0.58,  kFALSE, kFALSE) );
			//}
		}
		if(model=="gamma") {
			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE) );
		} else if(model=="gamma2") {
			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE) );
		} else if(model=="alt") {
			coeffset.push_back( new LauCartesianGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0., 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
		} else {
			coeffset.push_back( new LauRealImagGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61,  0., 0., 0., 0., kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
		}
	} else {
		if(model=="gamma" || gamma892) {
			//			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0(892)",  0.03, -1.32, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
			coeff = new LauPolarGammaCPCoeffSet("K*0(892)",  0.03, -1.32, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE);
			if(blind) {
				coeff->blindParameter("r","r ela etkahan oli. Sai",                2.);
				coeff->blindParameter("delta","delta sahvoorin nuo semmoiset",     LauConstants::pi*0.6);
				coeff->blindParameter("gamma","gamma Enempaa tavaraa pidatte kai", LauConstants::pi*0.6);
			}
			coeffset.push_back( coeff );
		} else if(model=="gamma2") {
			coeff = new LauPolarGammaCPCoeffSet("K*0(892)",  0.03, -1.32, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE);
			if(blind) {
				coeff->blindParameter("r","r sellainen asiakseen en ai",     2.);
				coeff->blindParameter("delta", "delta kaikkiaan se. Nuo",    LauConstants::pi*0.6);
				coeff->blindParameter("gamma", "gamma merirosvo. Koettelen", LauConstants::pi*0.6);
			}
			coeffset.push_back( coeff );
		} else if(model=="alt") {
			if(fixedDeltas) {
				coeff = new LauCartesianGammaCPCoeffSet("K*0(892)",  0.03, -1.32, 0.07, 0.03, 0.00, 0.00, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE);
			} else {
				coeff = new LauCartesianGammaCPCoeffSet("K*0(892)",  0.03, -1.32, 0.07, 0.03, -0.00, 0.00, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
			}
			if(blind) {
				coeff->blindParameter("XCP", "XCP jaa tuo ulapan emanta jai/",            2.);
				coeff->blindParameter("YCP", "YCP virkkoi jos ole",                       2.);
				coeff->blindParameter("DeltaXCP", "DeltaXCP Vie iso semmoinen ajelehtii", 0.4);
				coeff->blindParameter("DeltaYCP", "DeltaYCP avaimen. Sisaan paasta",      0.4);
			}
			coeffset.push_back( coeff );
		} else {
			coeff = new LauRealImagGammaCPCoeffSet("K*0(892)",   0.03, -1.32, -0.01, 0.25, 0.16, -0.20, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
			if(blind) {
				coeff->blindParameter("XCP", "XCP varmuuden osa. Tai",                      2.);
				coeff->blindParameter("YCP", "YCP luo sittenkuin kitupiikki",               2.);
				coeff->blindParameter("XbarCP", "XbarCP tulisikaan enemmankin. Puhisee",    2.);
				coeff->blindParameter("YbarCP", "YbarCP oli vanhoiksi tarvitaan kovinkaan", 2.);
			}
			coeffset.push_back( coeff );
		}
		if(!no1410) coeffset.push_back( new LauRealImagCoeffSet("K*0(1410)",  -0.24, -0.02,  kFALSE, kFALSE) );

		if(LASSCPV) {
			if(model=="gamma") {
				coeff = new LauPolarGammaCPCoeffSet("K*0_0(1430)",  0.03, -1.32, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE);
			} else {
				coeff = new LauRealImagGammaCPCoeffSet("K*0_0(1430)",   -0.43,  0.54, -0.01, 0.25, 0.16, -0.20, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
				if(blind) {
					coeff->blindParameter("XCP", "XCP On nalkainen teistakin kenenkaan on et",                2.);
					coeff->blindParameter("YCP", "YCP Tuvassa tapahdu oma kas monissa kynansa nyt",           2.);
					coeff->blindParameter("XbarCP", "XbarCP Saa vakituiset poikimatta rukoukseni tietamatta", 2.);
					coeff->blindParameter("YbarCP", "YbarCP tai kaksisataa puheenaihe",                       2.);
				}
			}
			coeffset.push_back( coeff );
			if(!mipw && glass) {
				if(extraLASSCPV) {
					if(model=="gamma") {
						coeff = new LauPolarGammaCPCoeffSet(kpiNRName,  0.03, -1.32, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE);
					} else {
						coeff = new LauRealImagGammaCPCoeffSet(kpiNRName,    -0.40,  0.58, -0.01, 0.25, 0.16, -0.20, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
						if(blind) {
							coeff->blindParameter("XCP", "XCP Kymmeneksi nyt ela vatvotusta liikkeella",  2.);
							coeff->blindParameter("YCP", "YCP rukoiltiin kitupiikki oli",                 2.);
							coeff->blindParameter("XbarCP", "XbarCP Rikki he vetta no ja silla revon",    2.);
							coeff->blindParameter("YbarCP", "YbarCP Hankkii et menevan ei se soittaa",    2.);
						}
					}
				} else {
					coeff = coeff->createClone(kpiNRName, LauAbsCoeffSet::TieCPPars); 
				}
				coeffset.push_back( coeff );
			}
		} else {
			coeffset.push_back( new LauRealImagCoeffSet("K*0_0(1430)",-0.43,  0.54,  kFALSE, kFALSE) );
			if(!mipw && glass) coeffset.push_back( new LauRealImagCoeffSet(kpiNRName,    -0.40,  0.58,  kFALSE, kFALSE) );
		}

		if(model=="gamma") {
			//			coeffset.push_back( new LauPolarGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE) );
			coeff = new LauPolarGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE);
			if(blind) {
				coeff->blindParameter("r","r Ole laskea toista nae omille",         2.);
				coeff->blindParameter("delta","delta milla ota pyori joita huone.", LauConstants::pi*0.6);
			}
			coeffset.push_back( coeff );
		} else if(model=="gamma2") {
			coeff = new LauPolarGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0.25, 0.35, 1.25, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE);
			if(blind) {
				coeff->blindParameter("r","r ai kesat esiin asken",       2.);
				coeff->blindParameter("delta", "delta isa jaa merkit leipaa", LauConstants::pi*0.6);
				coeff->blindParameter("gamma","gamma Osaat ai en usvaa nehan", LauConstants::pi*0.6);
			}
			coeffset.push_back( coeff );
		} else if(model=="alt") {
			coeff = new LauCartesianGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, 0.07, 0.03, -0.08, 0.22, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
			if(blind) {
				coeff->blindParameter("XCP", "XCP ne en saisit. Saa viittiloi",         2.);
				coeff->blindParameter("YCP", "YCP vai tai viinojaan unelmanne sen",     2.);
				coeff->blindParameter("DeltaXCP", "DeltaXCP Tuoma kalaa ole kay sanot", 0.4);
				coeff->blindParameter("DeltaYCP", "DeltaYCP Ja ruuhen sijaan me",       0.4);
			}
			coeffset.push_back( coeff );
		} else {
			coeff = new LauRealImagGammaCPCoeffSet("K*0_2(1430)",-0.16, -0.61, -0.01, 0.25, 0.16, -0.20, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE);
			if(blind) {
				coeff->blindParameter("XCP", "XCP Venhekin moottori porstuan",             2.);
				coeff->blindParameter("YCP", "YCP jattaisi kuuluvat voi tuo saarella loi", 2.);
				coeff->blindParameter("XbarCP", "XbarCP Et ai ryit kesy on ilta",          2.);
				coeff->blindParameter("YbarCP", "YbarCP ela paksu luo viina venhe",        2.);
			}
			coeffset.push_back( coeff );
		}
	}
	coeffset.push_back( new LauRealImagCoeffSet("D*-_0",        2.16909e-01,  2.00294e-01, kFALSE, kFALSE) );
	if(prodAsym) {
		coeffset.push_back( new LauCartesianCPCoeffSet("D*-_2", 1.00, 0.00, 0.00, 0.00, kTRUE, kTRUE, kFALSE, kTRUE) );
	} else {
		coeffset.push_back( new LauRealImagCoeffSet("D*-_2",        1.00,         0.00,        kTRUE,  kTRUE) );
	}
	//coeffset.push_back( new LauRealImagCoeffSet("D-(2760)",     2.43779e-01, -1.38701e-01, kFALSE, kFALSE) );
	if(altDpiSWave) {
		coeffset.push_back( new LauRealImagCoeffSet("BelleNR_Swave-",     4.51147e-01, -4.22586e-01, kFALSE, kFALSE) );
	} else {
		coeffset.push_back( new LauRealImagCoeffSet("dabba-",     4.51147e-01, -4.22586e-01, kFALSE, kFALSE) );
	}
	coeffset.push_back( new LauRealImagCoeffSet("BelleNR_Pwave-",     4.51147e-01, -4.22586e-01, kFALSE, kFALSE) );
	if(Ds) {
		if(dDecay=="kpi") {
			coeffset.push_back( new LauRealImagCoeffSet("Ds*+_2(2573)", 0.00, 0.00, kTRUE, kTRUE) );
		} else {
			if(DsCPV) {
				coeff = new LauRealImagCPCoeffSet("Ds*+_2(2573)", 0.01, 0.01, 0.01, 0.01, kFALSE, kFALSE, kFALSE, kFALSE);
				if(blind) {
					coeff->blindParameter("X", "X Pilkkanaan valahtivat pyyhkimaan kaupunkien",     2.);
					coeff->blindParameter("Y", "Y kaupunkiin han jos jaa ero",                      2.);
					coeff->blindParameter("Xbar", "Xbar Han pahaa hyvin menet isa",                 2.);
					coeff->blindParameter("Ybar", "Ybar Ai et lasianne oletteko helgahan tuhansia", 2.);
				}
				//coeff = new LauCartesianCPCoeffSet("Ds*+_2(2573)", 0.01, 0.01, 0.00, 0.00, kFALSE, kFALSE, kTRUE, kTRUE);
				//if(blind) {
				//	coeff->blindParameter("DeltaX", "DeltaX Ai et lasianne oletteko helgahan tuhansia", 0.4);
				//	coeff->blindParameter("DeltaY", "DeltaY Han pahaa hyvin menet isa",                 0.4);
				//}
			} else {
				coeff = new LauRealImagCoeffSet("Ds*+_2(2573)", 0.01, 0.01, kFALSE, kFALSE);
			}
			coeffset.push_back( coeff );
		}
	}
	if(DsSWave) {
		if(dDecay=="kpi") {
			coeffset.push_back( new LauRealImagCoeffSet("BelleNR+", 0.00, 0.00, kTRUE, kTRUE) );
		} else {
			if(DsSWaveCPV) {
				coeff = new LauRealImagCPCoeffSet("BelleNR+", 0.01, 0.01, 0.01, 0.01, kFALSE, kFALSE, kFALSE, kFALSE);
				if(blind) {
					coeff->blindParameter("X", "X kaupunkiin han jos jaa ero",                             2.);
					coeff->blindParameter("Y", "Y Pilkkanaan valahtivat pyyhkimaan kaupunkien",            2.);
					coeff->blindParameter("Xbar", "Xbar Ai et lasianne oletteko helgahan tuhansia",        2.);
					coeff->blindParameter("Ybar", "Ybar Han pahaa hyvin menet isa",                        2.);
				}
			} else {
				coeff = new LauRealImagCoeffSet("BelleNR+", 0.01, 0.01, kFALSE, kFALSE);
			}
			coeffset.push_back( coeff );
		}
	}
	if(DsPWave) {
		if(dDecay=="kpi") {
			coeffset.push_back( new LauRealImagCoeffSet("Ds*+_1(2700)", 0.00, 0.00, kTRUE, kTRUE) );
		} else {
			if(DsPWaveCPV) {
				if(DsPWaveCPVAlt) {
					coeff = new LauCartesianCPCoeffSet("Ds*+_1(2700)", 0.1, 0.1, 0.1, 0.1, kFALSE, kFALSE, kFALSE, kFALSE);
					if(blind) {
						coeff->blindParameter("DeltaX", "DeltaX Ai et lasianne oletteko helgahan tuhansia", 0.4);
						coeff->blindParameter("DeltaY", "DeltaY Han pahaa hyvin menet isa",      0.4);
					}
				} else {
					coeff = new LauCleoCPCoeffSet("Ds*+_1(2700)", 0.1, 0.0, 0.0, 0.0, kFALSE, kFALSE, kFALSE, kFALSE);
					//					coeff = new LauRealImagCPCoeffSet("Ds*+_1(2700)", 0.01, 0.01, 0.01, 0.01, kFALSE, kFALSE, kFALSE, kFALSE);
					if(blind) {
						//coeff->blindParameter("X", "X Ai et lasianne oletteko helgahan tuhansia",         2.);
						//coeff->blindParameter("Y", "Y Han pahaa hyvin menet isa",                         2.);
						//coeff->blindParameter("Xbar", "Xbar kaupunkiin han jos jaa ero",                  2.);
						//coeff->blindParameter("Ybar", "Ybar Pilkkanaan valahtivat pyyhkimaan kaupunkien", 2.);
						//	coeff->blindParameter("A",     "A Ai et lasianne oletteko helgahan tuhansia",     2.);
						coeff->blindParameter("Delta", "delta Han pahaa hyvin menet isa",                 LauConstants::pi*0.6);
						coeff->blindParameter("B",     "B kaupunkiin han jos jaa ero",                    0.4);
						coeff->blindParameter("Phi",   "phi Pilkkanaan valahtivat pyyhkimaan kaupunkien", LauConstants::pi*0.6);
					}
				}
			} else {
				//				if(constrainDK) coeff->addGaussianConstraint("A", TMath::Sqrt(0.217), TMath::Sqrt(0.055));
				if(fixDK) {
					coeff = new LauMagPhaseCoeffSet("Ds*+_1(2700)", TMath::Sqrt(0.217), 0., kTRUE, kFALSE);
				} else {
					coeff = new LauRealImagCoeffSet("Ds*+_1(2700)", 0.1, 0.1, kFALSE, kFALSE);
					//coeff = new LauMagPhaseCoeffSet("Ds*+_1(2700)", 0.12, 0.0, kFALSE, kFALSE);
					//if(constrainDK) coeff->addGaussianConstraint("A", TMath::Sqrt(0.118), TMath::Sqrt(0.029));
				}
			}
			coeffset.push_back( coeff );
		}
	}

	coeffset.push_back( new LauMagPhaseCoeffSet("D*-",          3.35713e-01,  0.0,         kFALSE, kTRUE) );
	for (std::vector<LauAbsCoeffSet*>::iterator iter=coeffset.begin(); iter!=coeffset.end(); ++iter) {
		fitModel->setAmpCoeffSet(*iter);
	}

	//if(constrainDK) {
	//	if(DsPWave && !Ds && !DsSWave) {
	//		if(!DsPWaveCPV) {
	//			std::vector<TString> vs;
	//			vs.push_back("A9_X");
	//			vs.push_back("A9_Y");
	//			fitModel->addConstraint("A9_X*A9_X + A9_Y*A9_Y", vs, 0.118, 0.029);
	//				if(constrainDK) coeff->addGaussianConstraint("A", TMath::Sqrt(0.217), TMath::Sqrt(0.055));
	//		}
	//	}
	//}

	Int_t nSigEvents   = 500;
	Int_t nCombEvents  = 0;
	Int_t nDstkEvents  = 0;
	Int_t nDpipiEvents = 0;
	Int_t nDppiEvents  = 0;
	Int_t nDpkEvents   = 0;
	Int_t nDstkpiEvents= 0;
	Int_t nDkpiEvents  = 0;
	Bool_t fixNEvents = kTRUE;
	if(dDecay=="kpi") {
		if(NNbin=="a") {
			nSigEvents   = 597;
			nCombEvents  = 540;
			nDstkEvents  = 287+18;
			nDpipiEvents = 20;
			nDpkEvents   = 21;
		} else if(NNbin=="b") {
			nSigEvents   = 546;
			nCombEvents  = 58;
			nDstkEvents  = 31+2;
			nDpipiEvents = 18;
			nDpkEvents   = 19;
		} else if(NNbin=="c") {
			nSigEvents   = 585;
			nCombEvents  = 16;
			nDstkEvents  = 8+1;
			nDpipiEvents = 20;
			nDpkEvents   = 21;
		} else if(NNbin=="d") {
			nSigEvents   = 571;
			nCombEvents  = 6;
			nDstkEvents  = 3+0;
			nDpipiEvents = 19;
			nDpkEvents   = 20;
		} else if(NNbin=="e") {
			nSigEvents   = 540;
			nCombEvents  = 1;
			nDstkEvents  = 1+0;
			nDpipiEvents = 18;
			nDpkEvents   = 19;
		} else if(NNbin=="b2e") {
			nSigEvents   = 546  + 585 + 571 + 540;
			nCombEvents  = 58   + 16  + 6   + 1;
			nDstkEvents  = 31+2 + 8+1 + 3+0 + 1+0;
			nDpipiEvents = 18   + 20  + 19  + 18;
			nDpkEvents   = 19   + 21  + 20  + 19;
		} else if(NNbin=="c2e") {
			nSigEvents   = 585 + 571 + 540;
			nCombEvents  = 16  + 6   + 1;
			nDstkEvents  = 8+1 + 3+0 + 1+0;
			nDpipiEvents = 20  + 19  + 18;
			nDpkEvents   = 21  + 20  + 19;
		} else if(NNbin == "old") {
			nSigEvents = 2344.;
			nCombEvents = 684.;
			nDpipiEvents = 51.;
		}
	} else if(dDecay=="kk") {
		if(NNbin=="a") {
			nSigEvents   = 70;
			nCombEvents  = 173;
			nDpipiEvents = 4;
			nDppiEvents  = 11;
			nDstkpiEvents= 19;
			nDkpiEvents  = 5;
		} else if(NNbin=="b") {
			nSigEvents   = 63;
			nCombEvents  = 19;
			nDpipiEvents = 3;
			nDppiEvents  = 10;
			nDstkpiEvents= 28;
			nDkpiEvents  = 5;
		} else if(NNbin=="c") {
			nSigEvents   = 68;
			nCombEvents  = 9;
			nDpipiEvents = 4;
			nDppiEvents  = 10;
			nDstkpiEvents= 34;
			nDkpiEvents  = 5;
		} else if(NNbin=="d") {
			nSigEvents   = 73;
			nCombEvents  = 3;
			nDpipiEvents = 4;
			nDppiEvents  = 11;
			nDstkpiEvents= 28;
			nDkpiEvents  = 6;
		} else if(NNbin=="e") {
			nSigEvents   = 65;
			nCombEvents  = 0;
			nDpipiEvents = 3;
			nDppiEvents  = 10;
			nDstkpiEvents= 20;
			nDkpiEvents  = 5;
		}
	} else if(dDecay=="pipi") {
		if(NNbin=="a") {
			nSigEvents   = 36;
			nCombEvents  = 119;
			nDpipiEvents = 2;
			nDppiEvents  = 6;
			nDstkpiEvents= 9;
			nDkpiEvents  = 3;
		} else if(NNbin=="b") {
			nSigEvents   = 31;
			nCombEvents  = 17;
			nDpipiEvents = 2;
			nDppiEvents  = 5;
			nDstkpiEvents= 16;
			nDkpiEvents  = 2;
		} else if(NNbin=="c") {
			nSigEvents   = 38;
			nCombEvents  = 4;
			nDpipiEvents = 2;
			nDppiEvents  = 6;
			nDstkpiEvents= 15;
			nDkpiEvents  = 3;
		} else if(NNbin=="d") {
			nSigEvents   = 32;
			nCombEvents  = 3;
			nDpipiEvents = 2;
			nDppiEvents  = 5;
			nDstkpiEvents= 12;
			nDkpiEvents  = 3;
		} else if(NNbin=="e") {
			nSigEvents   = 31;
			nCombEvents  = 2;
			nDpipiEvents = 2;
			nDppiEvents  = 5;
			nDstkpiEvents= 10;
			nDkpiEvents  = 2;
		}
	}

	TString sigEventsName = "signalEvents" + dDecay + NNbin;
	LauParameter* nSig = new LauParameter(sigEventsName,nSigEvents,-2*nSigEvents,2*nSigEvents,fixNEvents);
	fitModel->setNSigEvents(nSig);

	// Set the number of experiments to generate or fit and which
	// experiment to start with
	fitModel->setNExpts( nExpt, firstExpt );

	// Optionally load in continuum background DP model histogram
	// (example syntax given in commented-out section)
	std::vector<TString> bkgndNames(2);
	if(NNbin!="old") { 
		bkgndNames.push_back("");
		bkgndNames.push_back("");
	}
	if(dDecay!="kpi") bkgndNames.push_back("");

	bkgndNames[0] = "comb"  + dDecay + NNbin;
	bkgndNames[1] = "dpipi" + dDecay + NNbin;

	if(dDecay=="kpi") {
		if(NNbin!="old") {
			bkgndNames[2] = "dpk"  + dDecay + NNbin;
			bkgndNames[3] = "dstk"  + dDecay + NNbin;
		}
	} else {
		bkgndNames[2] = "dppi"  + dDecay + NNbin;
		bkgndNames[3] = "dstkpi"  + dDecay + NNbin;
		bkgndNames[4] = "bs2dkpi"  + dDecay + NNbin;
	}

	fitModel->setBkgndClassNames(bkgndNames);

	// Name of the background must match an entry in bkgndNames
	LauParameter *nComb(0), *nDstk(0), *nDpipi(0), *nDppi(0), *nDpk(0), *nDstkpi(0), *nDkpi(0);

	nComb = new LauParameter(bkgndNames[0],nCombEvents,-2*nCombEvents,2*nCombEvents,fixNEvents);
	fitModel->setNBkgndEvents(nComb);
	nDpipi = new LauParameter(bkgndNames[1],nDpipiEvents,-2*nDpipiEvents,2*nDpipiEvents,fixNEvents);
	fitModel->setNBkgndEvents(nDpipi);
	if(dDecay=="kpi") {
		if(NNbin!="old") {
			nDpk = new LauParameter(bkgndNames[2],nDpkEvents,-2*nDpkEvents,2*nDpkEvents,fixNEvents);
			fitModel->setNBkgndEvents(nDpk);
			nDstk = new LauParameter(bkgndNames[3],nDstkEvents,-2*nDstkEvents,2*nDstkEvents,fixNEvents);
			fitModel->setNBkgndEvents(nDstk);
		}
	} else {
		nDppi = new LauParameter(bkgndNames[2],nDppiEvents,-2*nDppiEvents,2*nDppiEvents,fixNEvents);
		fitModel->setNBkgndEvents(nDppi);
		nDstkpi = new LauParameter(bkgndNames[3],nDstkpiEvents,-2*nDstkpiEvents,2*nDstkpiEvents,fixNEvents);
		fitModel->setNBkgndEvents(nDstkpi);
		nDkpi = new LauParameter(bkgndNames[4],nDkpiEvents,-2*nDkpiEvents,2*nDkpiEvents,fixNEvents);
		fitModel->setNBkgndEvents(nDkpi);
	}

	TString combFileName("");
	TFile* combFile(0);
	TString combFileName2("");
	TFile* combFile2(0);
	TString combHistName("");
	//	if(dDecay=="kpi") {
	//		combFileName = "bkg/d2"+dDecay+"_comb_WS.root";
	//		combFile = TFile::Open(combFileName.Data(), "read");
	//		if(NNbin=="a") combHistName="combA_SDP";
	//		else combHistName="combB2E_SDP";
	TH2 *combDP(0), *combDP2(0);
	if(combAsym) {
		combFileName = "bkg/d2"+dDecay+"_comb_Bz.root";
		combFile = TFile::Open(combFileName.Data(), "read");
		combFileName2 = "bkg/d2"+dDecay+"_comb_Bzb.root";
		combFile2 = TFile::Open(combFileName2.Data(), "read");
		if(NNbin=="a") combHistName="combA_SDP_sub";
		else combHistName="combB2E_SDP_sub";
		combDP = dynamic_cast<TH2*>(combFile->Get(combHistName)); // m', theta'
		combDP2 = dynamic_cast<TH2*>(combFile2->Get(combHistName)); // m', theta'
	} else {
		if(NNbin=="old") {
			combFileName = "bkg-DpiK/comb_SDP_bkg_bkgrndSub_Dst25.root";
			combFile = TFile::Open(combFileName.Data(), "read");
			combHistName = "comb_SDP";
		} else {
			combFileName = "bkg/d2"+dDecay+"_comb.root";
			combFile = TFile::Open(combFileName.Data(), "read");
			if(NNbin=="a" || singleNNShape) combHistName="combA_SDP_sub";
			else combHistName="combB2E_SDP_sub";
		}
		combDP = dynamic_cast<TH2*>(combFile->Get(combHistName)); // m', theta'
		combDP2 = dynamic_cast<TH2*>(combFile->Get(combHistName)); // m', theta'
	}
	combDP->Print();
	combDP2->Print();
	LauBkgndDPModel* combModelPos = new LauBkgndDPModel(posDaughters, vetoes);
	LauBkgndDPModel* combModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
	combModelPos->setBkgndHisto(combDP,  useInterpolation, fluctuateBins, useUpperHalf, squareDP);
	combModelNeg->setBkgndHisto(combDP2, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
	fitModel->setBkgndDPModels( bkgndNames[0], combModelPos, combModelNeg );

	TString dpipiFileName("");
	TFile* dpipiFile(0);
	TString dpipiHistName("");
	if(NNbin=="old") {
		dpipiFileName = "bkg-DpiK/Dpipi_Dstpipi_SDP_bkg_Dst25_60.root";
		dpipiFile = TFile::Open(dpipiFileName.Data(), "read");
		dpipiHistName = "Dpipi_SDP";
	} else {
		dpipiFileName = "bkg/d2"+dDecay+"_Bd2D0pipi.root";
		dpipiFile = TFile::Open(dpipiFileName.Data(), "read");
		dpipiHistName = "Bd2D0pipi_SDP";
	}
	TH2* dpipiDP = dynamic_cast<TH2*>(dpipiFile->Get(dpipiHistName)); //bkg")); // m', theta'
	dpipiDP->Print();
	LauBkgndDPModel* dpipiModelPos = new LauBkgndDPModel(posDaughters, vetoes);
	LauBkgndDPModel* dpipiModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
	dpipiModelPos->setBkgndHisto(dpipiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
	dpipiModelNeg->setBkgndHisto(dpipiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
	fitModel->setBkgndDPModels( bkgndNames[1], dpipiModelPos, dpipiModelNeg );

	if(dDecay=="kpi") {
		if(NNbin!="old") {
			TString dpkFileName("bkg/d2"+dDecay+"_Lb2D0pK.root");
			TFile* dpkFile = TFile::Open(dpkFileName.Data(), "read");
			TH2* dpkDP = dynamic_cast<TH2*>(dpkFile->Get("Lb2D0pK_SDP")); // m', theta'
			dpkDP->Print();
			LauBkgndDPModel* dpkModelPos = new LauBkgndDPModel(posDaughters, vetoes);
			LauBkgndDPModel* dpkModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
			dpkModelPos->setBkgndHisto(dpkDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
			dpkModelNeg->setBkgndHisto(dpkDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
			fitModel->setBkgndDPModels( bkgndNames[2], dpkModelPos, dpkModelNeg );

			TString dstkFileName("bkg/d2"+dDecay+"_Bu2DstK.root");
			TFile* dstkFile = TFile::Open(dstkFileName.Data(), "read");
			TString dstkHistName("");
			if(NNbin=="a" || singleNNShape) dstkHistName="partcombA_SDP";
			else dstkHistName="partcombB2E_SDP";
			//if(NNbin=="a") dstkHistName="Bu2DstKA_SDP";
			//else dstkHistName="Bu2DstKB2E_SDP";
			TH2* dstkDP = dynamic_cast<TH2*>(dstkFile->Get(dstkHistName)); // m', theta'
			dstkDP->Print();
			LauBkgndDPModel* dstkModelPos = new LauBkgndDPModel(posDaughters, vetoes);
			LauBkgndDPModel* dstkModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
			dstkModelPos->setBkgndHisto(dstkDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
			dstkModelNeg->setBkgndHisto(dstkDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
			fitModel->setBkgndDPModels( bkgndNames[3], dstkModelPos, dstkModelNeg );
		}
	} else {
		TString dppiFileName("bkg/d2"+dDecay+"_Lb2D0ppi.root");
		TFile* dppiFile = TFile::Open(dppiFileName.Data(), "read");
		TH2* dppiDP = dynamic_cast<TH2*>(dppiFile->Get("Lb2D0ppi_SDP")); // m', theta'
		dppiDP->Print();
		LauBkgndDPModel* dppiModelPos = new LauBkgndDPModel(posDaughters, vetoes);
		LauBkgndDPModel* dppiModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
		dppiModelPos->setBkgndHisto(dppiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		dppiModelNeg->setBkgndHisto(dppiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		fitModel->setBkgndDPModels( bkgndNames[2], dppiModelPos, dppiModelNeg );

		TString dstkpiFileName("bkg/d2"+dDecay+"_Bs2Dst0Kpi.root");
		if(altDstKpi) dstkpiFileName = "bkg/d2kk_Bs2Dst0KpiAlt.root";
		TFile* dstkpiFile = TFile::Open(dstkpiFileName.Data(), "read");
		TH2* dstkpiDP = dynamic_cast<TH2*>(dstkpiFile->Get("Bs2DstKpi_SDP")); // m', theta'
		dstkpiDP->Print();
		LauBkgndDPModel* dstkpiModelPos = new LauBkgndDPModel(posDaughters, vetoes);
		LauBkgndDPModel* dstkpiModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
		dstkpiModelPos->setBkgndHisto(dstkpiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		dstkpiModelNeg->setBkgndHisto(dstkpiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		fitModel->setBkgndDPModels( bkgndNames[3], dstkpiModelPos, dstkpiModelNeg );

		TString dkpiFileName("bkg/d2"+dDecay+"_Bs2D0Kpi.root");
		TFile* dkpiFile = TFile::Open(dkpiFileName.Data(), "read");
		TH2* dkpiDP = dynamic_cast<TH2*>(dkpiFile->Get("Bs2D0Kpi_SDP")); // m', theta'
		dkpiDP->Print();
		LauBkgndDPModel* dkpiModelPos = new LauBkgndDPModel(posDaughters, vetoes);
		LauBkgndDPModel* dkpiModelNeg = new LauBkgndDPModel(negDaughters, vetoes);
		dkpiModelPos->setBkgndHisto(dkpiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		dkpiModelNeg->setBkgndHisto(dkpiDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
		fitModel->setBkgndDPModels( bkgndNames[4], dkpiModelPos, dkpiModelNeg );
	}

	// Switch on/off calculation of asymmetric errors.
	fitModel->useAsymmFitErrors(kFALSE);

	// Randomise initial fit values for the signal mode
	fitModel->useRandomInitFitPars(kTRUE);

	// Switch on/off Poissonian smearing of total number of events
	fitModel->doPoissonSmearing(kTRUE);

	// Switch on/off Extended ML Fit option
	Bool_t emlFit = ( fitModel->nBkgndClasses() > 0 );
	fitModel->doEMLFit(emlFit);

	//fitModel->twoStageFit(kTRUE);

	//Generate Toy to check the fit
	TString toyName = "Bdfits_sim_NNbins/"; toyName += model;
	toyName += "/"; toyName += dDecay; toyName += "/";  toyName += NNbin;
	toyName += "/fitToy"; toyName += iFit; toyName += ".root";
	fitModel->compareFitData(100,toyName,"",kTRUE);

	// Set the names of the files to read/write
	TString dataFile("");
	if(dDecay=="kpi") {
		if(NNbin=="old") {
			dataFile="/data/lhcb/phrkbf/B2DKpi/Job2012_DKpi_Scaled/TSS_baseline_DKpi_selBd_allVetoes_Dst25_DD55_PID3p_NNf1_25_DpiKSDP.root";
		} else {
			dataFile="/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2Kpi_Scaled/B2D0Kpi_D02Kpi_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bd25_"+NNbin+".root";
		}
	} else if(dDecay=="kk") {
		dataFile="/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2KK_Scaled/B2D0Kpi_D02KK_selBd_Dsignal_Dst25_vetoes_PID3_NND2KK_addIMs_Bd25_"+NNbin+".root";
	} else if(dDecay=="pipi") {
		dataFile="/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2pipi_Scaled/B2D0Kpi_D02pipi_selBd_Dsignal_Dst25_vetoes_PID3_NND2pipi_addIMs_Bd25_"+NNbin+".root";
	}
	TString treeName("DecayTree");
	TString rootFileName("");
	TString tableFileName("");
	TString splotFileName("splot_");
	if (command == "fit") {
		rootFileName = "Bdfits_sim_NNbins/"; rootFileName += model;
		rootFileName += "/"; rootFileName += dDecay; 
		rootFileName += "/";  rootFileName += NNbin;
		rootFileName += "/fit"; rootFileName += iFit;
		rootFileName += ".root";
		tableFileName = "fitResults_"; tableFileName += iFit;
		splotFileName += iFit;
		splotFileName += ".root";
	} else {
		rootFileName = "dummy.root";
		tableFileName = "genResults";
	}

	// Execute the generation/fit
	if ( command == "fit" ) {
		fitModel->runSlave( dataFile, treeName, rootFileName, tableFileName, host, port );
	} else {
		fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );
	}

	return EXIT_SUCCESS;
}
