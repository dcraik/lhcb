#define skim_cxx
#include "skim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void skim::Loop()
{
	if (fChain == 0) return;

	TFile* fout = TFile::Open(outName,"RECREATE");
	//TTree *tout = new TTree("T","");
	TTree *tout = fChain->CloneTree(0);
	TTree *lumiout = lumi->CloneTree();

	double phi_l0_flags, phi_hlt1_flags, phi_hlt2_flags;

	//tout->Branch("phi_m", &phi_M);
	//tout->Branch("phi_p", &phi_P);
	//tout->Branch("phi_pt", &phi_PT);
	//tout->Branch("phi_ip_chi2", &phi_IPCHI2_OWNPV);
	//tout->Branch("phi_vtx_chi2", &phi_ENDVERTEX_CHI2);
	//tout->Branch("phi_fd_chi2", &phi_FDCHI2_OWNPV);

	//tout->Branch("phi_l0_flags", &phi_l0_flags);
	//tout->Branch("phi_hlt1_flags", &phi_hlt1_flags);
	//tout->Branch("phi_hlt2_flags", &phi_hlt2_flags);

	Long64_t nentries = fChain->GetEntries();

	int counts[10] = {0,0,0,0,0,0,0,0,0,0};

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		++counts[0];
		if(nLongTracks>2 || nVeloTracks>2) continue;
		++counts[1];
		//if(LoKi_nNeutrals>0) continue;
		++counts[2];
		if(phi_M>1200.) continue;
		++counts[3];
		if(Kplus_MC15TuneV1_ProbNNk<=0.1) continue;
		if(Kminus_MC15TuneV1_ProbNNk<=0.1) continue;
		++counts[4];

		phi_l0_flags=0;
		phi_hlt1_flags=0;
		phi_hlt2_flags=0;

		if(phi_L0DiHadron_lowMultDecision_TOS) phi_l0_flags+=1;
		
		if(phi_Hlt1LowMultVeloCut_HadronsDecision_Dec) phi_hlt1_flags+=1;
		if(phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec) phi_hlt1_flags+=2;

		if(phi_Hlt2LowMultLMR2HHDecision_TOS) phi_hlt2_flags+=1;

		if(phi_l0_flags==0 || phi_hlt1_flags==0 || phi_hlt2_flags==0) continue;
		++counts[5];
	
		tout->Fill();
	}

	for(int i=0; i<6; ++i) {
		std::cout << counts[i] << std::endl;
	}

	tout->AutoSave();
	lumiout->AutoSave();
	fout->Close();
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
