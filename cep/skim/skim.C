#define skim_cxx
#include "skim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

void skim::Loop()
{
	//first make files for phi tuples
	for(unsigned int iChain=0; iChain<phiChains.size(); ++iChain) {
		TChain* fChain = phiChains.at(iChain);
		TString which = phiNames.at(iChain);
		bool flipPID(false);
		if(which=="PiPi" || which=="PiPiWS") {
			flipPID=true;
		}
		if (fChain == 0) continue;

		TString outName=outPath;
		outName+="/skimmed";
		if(which!="") outName+="-"+which;
		outName+=".root";

		TFile* fout = TFile::Open(outName,"RECREATE");
		TTree *tout = fChain->CloneTree(0);
		TTree *lumiout = lumi->CloneTree();

		double pPb_l0_flags,  pPb_hlt1_flags,  pPb_hlt2_flags;
		double PbPb_l0_flags, PbPb_hlt1_flags, PbPb_hlt2_flags;

		tout->Branch("pPb_l0_flags",   &pPb_l0_flags);
		tout->Branch("pPb_hlt1_flags", &pPb_hlt1_flags);
		tout->Branch("pPb_hlt2_flags", &pPb_hlt2_flags);
		tout->Branch("PbPb_l0_flags",   &PbPb_l0_flags);
		tout->Branch("PbPb_hlt1_flags", &PbPb_hlt1_flags);
		tout->Branch("PbPb_hlt2_flags", &PbPb_hlt2_flags);

		Long64_t nentries = fChain->GetEntries();

		int counts[10] = {0,0,0,0,0,0,0,0,0,0};

		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			fChain->GetEntry(jentry);

			++counts[0];
			if(nLongTracks>4 || nVeloTracks>4) continue;//TODO keep up to 4 to perform "background" fits
			++counts[1];
			//if(LoKi_nNeutrals>0) continue;
			++counts[2];
			if(phi_M>1200.) continue;
			++counts[3];
			if(!flipPID && Kplus_MC15TuneV1_ProbNNk<=0.05) continue;
			if(!flipPID && Kminus_MC15TuneV1_ProbNNk<=0.05) continue;
			//for pipi, flip PID cuts and remove KS
			if(flipPID && Kplus_MC15TuneV1_ProbNNk>0.05 && Kminus_MC15TuneV1_ProbNNk>0.05) continue;
			if(flipPID && phi_ENDVERTEX_Z>100.) continue;
			++counts[4];

			pPb_l0_flags=0;
			pPb_hlt1_flags=0;
			pPb_hlt2_flags=0;
			PbPb_l0_flags=0;
			PbPb_hlt1_flags=0;
			PbPb_hlt2_flags=0;

			if(phi_L0DiHadron_lowMultDecision_Dec) pPb_l0_flags += 0x0001;
			if(phi_L0Muon_lowMultDecision_Dec) pPb_l0_flags += 0x0002;
			if(phi_L0DiMuon_lowMultDecision_Dec) pPb_l0_flags += 0x0004;
			if(phi_L0Electron_lowMultDecision_Dec) pPb_l0_flags += 0x0008;
			if(phi_L0Photon_lowMultDecision_Dec) pPb_l0_flags += 0x0010;
			if(phi_L0DiEM_lowMultDecision_Dec) pPb_l0_flags += 0x0020;
			if(phi_L0SPDDecision_Dec) pPb_l0_flags += 0x0040;
			if(phi_L0CALODecision_Dec) pPb_l0_flags += 0x0080;
			if(phi_L0PUDecision_Dec) pPb_l0_flags += 0x0100;
			if(phi_L0MuonDecision_Dec) pPb_l0_flags += 0x0200;
			if(phi_Hlt1LowMultPassThroughDecision_Dec) pPb_hlt1_flags += 0x0001;
			if(phi_Hlt1LowMultVeloCut_HadronsDecision_Dec) pPb_hlt1_flags += 0x0002;
			if(phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec) pPb_hlt1_flags += 0x0004;
			if(phi_Hlt1BBMicroBiasVeloDecision_Dec) pPb_hlt1_flags += 0x0008;
			if(phi_Hlt1BBHighMultDecision_Dec) pPb_hlt1_flags += 0x0010;
			if(phi_Hlt1MBNoBiasDecision_Dec) pPb_hlt1_flags += 0x0020;
			if(phi_Hlt2LowMultLMR2HHDecision_Dec) pPb_hlt2_flags += 0x0001;
			if(phi_Hlt2LowMultLMR2HHWSDecision_Dec) pPb_hlt2_flags += 0x0002;
			if(phi_Hlt2MBMicroBiasVeloDecision_Dec) pPb_hlt2_flags += 0x0004;
			if(phi_Hlt2MBHighMultDecision_Dec) pPb_hlt2_flags += 0x0008;
			if(phi_Hlt2MBNoBiasDecision_Dec) pPb_hlt2_flags += 0x0010;
			if(phi_Hlt2PassThroughDecision_Dec) pPb_hlt2_flags += 0x0020;

			if(phi_L0HadronDecision_Dec) PbPb_l0_flags += 0x0001;
			if(phi_L0MuonDecision_Dec) PbPb_l0_flags += 0x0002;
			if(phi_L0SPDDecision_Dec) PbPb_l0_flags += 0x0004;
			if(phi_L0HadronLowMultDecision_Dec) PbPb_l0_flags += 0x0008;
			if(phi_L0MuonLowMultDecision_Dec) PbPb_l0_flags += 0x0010;
			if(phi_L0ElectronLowMultDecision_Dec) PbPb_l0_flags += 0x0020;
			if(phi_L0PhotonLowMultDecision_Dec) PbPb_l0_flags += 0x0040;
			if(phi_L0DiEMLowMultDecision_Dec) PbPb_l0_flags += 0x0080;
			if(phi_L0SPDLowMultDecision_Dec) PbPb_l0_flags += 0x0100;
			if(phi_L0SoftCEPDecision_Dec) PbPb_l0_flags += 0x0200;
			if(phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec) PbPb_hlt1_flags += 0x0001;
			if(phi_Hlt1BBMicroBiasSoftCEPDecision_Dec) PbPb_hlt1_flags += 0x0002;
			if(phi_Hlt1BBHasTrackDecision_Dec) PbPb_hlt1_flags += 0x0004;
			if(phi_Hlt2BBPassThroughDecision_Dec) PbPb_hlt2_flags += 0x0001;
			if(phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec) PbPb_hlt2_flags += 0x0002;
			if(phi_Hlt2BBLongTrackDecision_Dec) PbPb_hlt2_flags += 0x0004;
			if(phi_Hlt2SingleTrackDecision_Dec) PbPb_hlt2_flags += 0x0008;
			if(phi_Hlt2MBMicroBiasVeloDecision_Dec) PbPb_hlt2_flags += 0x0010;

			++counts[5];
		
			tout->Fill();
		}

		for(int i=0; i<6; ++i) {
			std::cout << counts[i] << "\t";
		}
		std::cout << std::endl;

		tout->AutoSave();
		lumiout->AutoSave();
		fout->Close();
	}
	////make file for Z tuple
	//if (zChain == 0) return;

	//TString outName=outPath;
	//outName+="/skimmed-Z.root";

	//TFile* fout = TFile::Open(outName,"RECREATE");
	//TTree *tout = zChain->CloneTree(0);
	//TTree *lumiout = lumi->CloneTree();

	//Long64_t nentries = zChain->GetEntries();

	//int counts[10] = {0,0,0,0,0,0,0,0,0,0};

	//for (Long64_t jentry=0; jentry<nentries;jentry++) {
	//	//selection from https://cds.cern.ch/record/2017816/files/LHCb-ANA-2015-029.pdf
	//	zChain->GetEntry(jentry);
	//	++counts[0];
	//	if(Z_M<60e3 || Z_M>120e3) continue;
	//	++counts[1];
	//	if(muplus_PT<20e3 || muminus_PT<20e3) continue;
	//	++counts[2];
	//	TLorentzVector p4p(muplus_PX,muplus_PY,muplus_PZ,muplus_PE);
	//	TLorentzVector p4m(muminus_PX,muminus_PY,muminus_PZ,muminus_PE);
	//	if(p4p.Eta()<2. || p4p.Eta()>4.5 || p4m.Eta()<2. || p4m.Eta()>4.5) continue;
	//	++counts[3];
	//	//TODO add dp/p<10%
	//	++counts[4];
	//	if(muplus_TRACK_PCHI2<0.01 || muplus_TRACK_PCHI2<0.01) continue;
	//	++counts[5];
	//	if(Z_L0MuonEWDecision_TOS==0) continue;
	//	++counts[6];
	//	if(Z_Hlt1SingleMuonHighPTDecision_TOS==0) continue;
	//	++counts[7];
	//	if(Z_Hlt2SingleMuonHighPTDecision_TOS==0) continue;
	//	++counts[8];
	//	tout->Fill();
	//}

	//for(int i=0; i<9; ++i) {
	//	std::cout << counts[i] << "\t";
	//}
	//std::cout << std::endl;

	//tout->AutoSave();
	//lumiout->AutoSave();
	//fout->Close();

}

int main(int argc, char** argv) {
	int job(0), sjob(-1);
	TString dir="/tmp/dcraik";
	if(argc>1) job = atoi(argv[1]);
	if(argc>2) sjob = atoi(argv[2]);
	if(argc>3) dir = argv[3];
	skim a(job,sjob,dir);
	a.Loop();
	return 0;
}
