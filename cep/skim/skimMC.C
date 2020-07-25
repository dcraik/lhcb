#define skimMC_cxx
#include "skimMC.h"

#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void skimMC::fillOutput(int gp, int gkp, int gkm, int rp, int rkp, int rkm)
{
	Kplus_TRUEP  = gen_p ->at(gkp);
	Kplus_TRUEPT = gen_pt->at(gkp);
	Kplus_TRUEPX = gen_px->at(gkp);
	Kplus_TRUEPY = gen_py->at(gkp);
	Kplus_TRUEPZ = gen_pz->at(gkp);
	Kplus_TRUEE  = gen_e-> at(gkp);

	Kminus_TRUEP  = gen_p ->at(gkm);
	Kminus_TRUEPT = gen_pt->at(gkm);
	Kminus_TRUEPX = gen_px->at(gkm);
	Kminus_TRUEPY = gen_py->at(gkm);
	Kminus_TRUEPZ = gen_pz->at(gkm);
	Kminus_TRUEE  = gen_e-> at(gkm);

	phi_TRUEP  = gen_p ->at(gp);
	phi_TRUEPT = gen_pt->at(gp);
	phi_TRUEPX = gen_px->at(gp);
	phi_TRUEPY = gen_py->at(gp);
	phi_TRUEPZ = gen_pz->at(gp);
	phi_TRUEE  = gen_e-> at(gp);

	if(rkp>-1) {
		Kplus_P  = trk_p ->at(rkp);
		Kplus_PT = trk_pt->at(rkp);
		Kplus_PX = trk_px->at(rkp);
		Kplus_PY = trk_py->at(rkp);
		Kplus_PZ = trk_pz->at(rkp);
		Kplus_E  = trk_e-> at(rkp);
	} else {
		Kplus_P  = 0.;
		Kplus_PT = 0.;
		Kplus_PX = 0.;
		Kplus_PY = 0.;
		Kplus_PZ = 0.;
		Kplus_E  = 0.;
	}

	if(rkm>-1) {
		Kminus_P  = trk_p ->at(rkm);
		Kminus_PT = trk_pt->at(rkm);
		Kminus_PX = trk_px->at(rkm);
		Kminus_PY = trk_py->at(rkm);
		Kminus_PZ = trk_pz->at(rkm);
		Kminus_E  = trk_e-> at(rkm);
	} else {
		Kminus_P  = 0.;
		Kminus_PT = 0.;
		Kminus_PX = 0.;
		Kminus_PY = 0.;
		Kminus_PZ = 0.;
		Kminus_E  = 0.;
	}

	if(rp>-1) {
		phi_P  = phi_p ->at(rp);
		phi_PT = phi_pt->at(rp);
		phi_PX = phi_px->at(rp);
		phi_PY = phi_py->at(rp);
		phi_PZ = phi_pz->at(rp);
		phi_E  = phi_e-> at(rp);
		phi_M  = phi_m-> at(rp);
	} else {
		phi_P  = 0.;
		phi_PT = 0.;
		phi_PX = 0.;
		phi_PY = 0.;
		phi_PZ = 0.;
		phi_E  = 0.;
		phi_M  = 0.;
	}
	tout->Fill();
}

void skimMC::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		for(uint gp=0; gp<gen_pid->size(); ++gp) {
			if(gen_pid->at(gp)!=333) continue;
			for(uint gkp=0; gkp<gen_pid->size(); ++gkp) {
				if(gen_pid->at(gkp)!=321) continue;
				if(gen_idx_prnt->at(gkp)!=gp) continue;
				for(uint gkm=0; gkm<gen_pid->size(); ++gkm) {
					if(gen_pid->at(gkm)!=-321) continue;
					if(gen_idx_prnt->at(gkm)!=gp) continue;
					//found a true phi->KK decay
					//now see if any are reconstructed
					int rkp(-1), rkm(-1), rp(-1);
					for(uint t=0; t<trk_idx_gen->size(); ++t) {
						if(trk_idx_gen->at(t)==gkp) {
							rkp = t;
							break;
						}
					}
					for(uint t=0; t<trk_idx_gen->size(); ++t) {
						if(trk_idx_gen->at(t)==gkm) {
							rkm = t;
							break;
						}
					}
					if(rkp>-1 && rkm>-1) {
						for(uint p=0; p<phi_m->size(); ++p) {
							if((phi_idx_trk0->at(p)==rkp && phi_idx_trk1->at(p)==rkm) ||
							   (phi_idx_trk0->at(p)==rkm && phi_idx_trk1->at(p)==rkp)) {
								rp = p;
								break;
							}
						}
					}
					fillOutput(gp,gkp,gkm,rp,rkp,rkm);

				}
			}
		}
	}

	tout->AutoSave();
	fout->Close();
}

int main(int argc, char** argv) {
	int job(0), sjob(-1);
	TString dir="/tmp/dcraik";
	if(argc>1) job = atoi(argv[1]);
	if(argc>2) sjob = atoi(argv[2]);
	if(argc>3) dir = argv[3];
	skimMC a(job,sjob,dir);
	a.Loop();
	return 0;
}
