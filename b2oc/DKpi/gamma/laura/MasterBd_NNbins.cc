
#include <cstdlib>
#include <iostream>

#include <vector>

#include "TFile.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"

#include "LauSimFitMaster.hh"

void usage( std::ostream& out, const TString& progName )
{
	out<<"Usage:\n";
	out<<progName<<" <iFit> <nExpt> [model = \"nominal\"] [firstExpt = 0] [numSlaves = 15] [port = 0]\n";
}

int main(const int argc, const  char ** argv)
{
	if ( argc < 3 ) {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

	UInt_t iFit = atoi( argv[1] );
	UInt_t nExpt = atoi( argv[2] );
	TString model = "nominal";
	UInt_t firstExpt = 0;
	UInt_t nSlaves = 15;
	UInt_t port = 0;

	if ( argc > 3 ) {
		model = argv[3];

		if ( argc > 4 ) {
			firstExpt = atoi( argv[4] );

			if ( argc > 5 ) {
				nSlaves = atoi( argv[5] );

				if ( argc > 6 ) {
					port = atoi( argv[6] );
				}
			}
		}
	}

	Bool_t constrainDK(kFALSE), DKCPV(kFALSE), DsDWave(kFALSE), sLASS(kFALSE), DKCPVAlt(kFALSE), newConstraint(kFALSE), constrainDWave(kFALSE), no1410(kFALSE);

	switch(iFit/100) {
		case 69:
			constrainDK = kTRUE;
			std::cout << "running 6900-series fit..." << std::endl;
			break;
		case 71:
			constrainDK = kTRUE;
			DsDWave = kTRUE;
			std::cout << "running 7100-series fit..." << std::endl;
			break;
		case 72:
			constrainDK = kTRUE;
			DKCPV = kTRUE;
			std::cout << "running 7200-series fit..." << std::endl;
			break;
		case 73:
			constrainDK = kTRUE;
			sLASS = kTRUE;
			std::cout << "running 7300-series fit..." << std::endl;
			break;
		case 74:
			constrainDK = kTRUE;
			DKCPV = kTRUE;
			DKCPVAlt = kTRUE;
			std::cout << "running 7400-series fit..." << std::endl;
			break;
		case 77:
		case 78:
			constrainDK = kTRUE;
			newConstraint = kTRUE;
			std::cout << "running 7700-series fit..." << std::endl;
			break;
		case 79:
		case 83:
		case 84:
		case 85:
		case 86:
		case 87:
		case 92:
		case 93:
		case 94:
		case 95:
		case 97:
		case 98:
		case 99:
		case 100:
		case 101:
		case 102:
		case 103:
		case 105:
			constrainDK = kTRUE;
			newConstraint = kTRUE;
			std::cout << "running 7900-series fit..." << std::endl;
			break;
		case 80:
			constrainDK = kTRUE;
			DKCPV = kTRUE;
			DKCPVAlt = kTRUE;
			newConstraint = kTRUE;
			std::cout << "running 8000-series fit..." << std::endl;
			break;
		case 81:
			constrainDK = kTRUE;
			DsDWave = kTRUE;
			constrainDWave = kTRUE;
			newConstraint = kTRUE;
			std::cout << "running 8100-series fit..." << std::endl;
			break;
		case 82:
		case 88:
		case 89:
		case 90:
		case 91:
			constrainDK = kTRUE;
			newConstraint = kTRUE;
			sLASS = kTRUE;
			std::cout << "running 8200-series fit..." << std::endl;
			break;
		case 104:
		case 106:
			constrainDK = kTRUE;
			newConstraint = kTRUE;
			no1410 = kTRUE;
			std::cout << "running 10400-series fit..." << std::endl;
			break;
	}

	TString ntupleName = "Bdfits_sim_NNbins/";
	ntupleName += model;
	ntupleName += "/master-ntuple-";
	ntupleName += iFit;
	ntupleName += ".root";

	LauSimFitMaster master( nSlaves, port );

	if(constrainDK) {
		std::vector<TString> vs;
		if(DsDWave) {
			vs.push_back("A10_X");
			vs.push_back("A10_Y");
		} else {
			if(sLASS || no1410) {
				vs.push_back("A8_X");
				vs.push_back("A8_Y");
			} else {
				if(DKCPV && !DKCPVAlt) {
					vs.push_back("A9_A");
				} else {
					vs.push_back("A9_X");
					vs.push_back("A9_Y");
				}
			}
		}
		if(DKCPV && !DKCPVAlt) {
			master.addConstraint("[0]*[0]", vs, 0.217, 0.055);
		} else {
			if(newConstraint) {
				master.addConstraint("[0]*[0] + [1]*[1]", vs, 0.217, 0.085);
			} else {
				master.addConstraint("[0]*[0] + [1]*[1]", vs, 0.217, 0.055);
			}
		}
	}

	if(constrainDWave) {
		std::vector<TString> vs;
		vs.push_back("A9_X");
		vs.push_back("A9_Y");
		vs.push_back("A10_X");
		vs.push_back("A10_Y");
		master.addConstraint("([0]*[0] + [1]*[1])/([2]*[2] + [3]*[3])", vs, 0.048, 0.026);
	}
	
	master.runSimFit( ntupleName, nExpt, firstExpt );

	return EXIT_SUCCESS;
}
