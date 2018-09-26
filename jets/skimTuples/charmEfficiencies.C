#define charmEfficiencies_cxx
#include "charmEfficiencies.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TVector3.h>

bool charmEfficiencies::checkForD0(int idxD, int idx0, int idx1, int idxNxt, int origin, bool ordered) {
	//first check indices are in range
	if(idxD<0 || idxD>=static_cast<int>(gen_pid->size())) return false;
	if(idx0<0 || idx0>=static_cast<int>(gen_pid->size())) return false;
	if(idx1<0 || idx1>=static_cast<int>(gen_pid->size())) return false;

	//now check we have a D0
	if(TMath::Abs(gen_pid->at(idxD))!=421.) return false;

	//now check the hierarchy
	if(gen_idx_prnt->at(idx0) != idxD) return false;
	if(gen_idx_prnt->at(idx1) != idxD) return false;

	//check which type of quark produced the D
	if(origin>0) {
		int idxPrnt = idxD;
		while(gen_idx_prnt->at(idxPrnt)>=0 && gen_idx_prnt->at(idxPrnt)<static_cast<int>(gen_pid->size())) {
			idxPrnt = gen_idx_prnt->at(idxPrnt);
		}
		if(TMath::Abs(gen_prnt_pid->at(idxPrnt)) != origin) return false;
	}


	//if we have a next particle then check it is either a photon (assume radiative) or unrelated
	if(idxNxt>=0 && idxNxt<static_cast<int>(gen_pid->size())) {
		if(gen_pid->at(idxNxt) != 22. && gen_idx_prnt->at(idxNxt) == idxD) return false;
	}

	//check for the correct number of pions and kaons
	int npi(0), nk(0);
	if(TMath::Abs(gen_pid->at(idx0)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx0)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx1)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx1)) == 321.) ++nk;

	if(nk!=1 || npi!=1) return false;

	//if ordered then require the kaon to be in idx0
	if(ordered && TMath::Abs(gen_pid->at(idx0)) != 321.) return false;

	return true;
}

bool charmEfficiencies::checkForDp(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered) {
	//first check indices are in range
	if(idxD<0 || idxD>=static_cast<int>(gen_pid->size())) return false;
	if(idx0<0 || idx0>=static_cast<int>(gen_pid->size())) return false;
	if(idx1<0 || idx1>=static_cast<int>(gen_pid->size())) return false;
	if(idx2<0 || idx2>=static_cast<int>(gen_pid->size())) return false;

	//now check we have a Dp
	if(TMath::Abs(gen_pid->at(idxD))!=411.) return false;

	//now check the hierarchy
	if(gen_idx_prnt->at(idx0) != idxD) return false;
	if(gen_idx_prnt->at(idx1) != idxD) return false;
	if(gen_idx_prnt->at(idx2) != idxD) return false;

	//check which type of quark produced the D
	if(origin>0) {
		int idxPrnt = idxD;
		while(gen_idx_prnt->at(idxPrnt)>=0 && gen_idx_prnt->at(idxPrnt)<static_cast<int>(gen_pid->size())) {
			idxPrnt = gen_idx_prnt->at(idxPrnt);
		}
		if(TMath::Abs(gen_prnt_pid->at(idxPrnt)) != origin) return false;
	}


	//if we have a next particle then check it is either a photon (assume radiative) or unrelated
	if(idxNxt>=0 && idxNxt<static_cast<int>(gen_pid->size())) {
		if(gen_pid->at(idxNxt) != 22. && gen_idx_prnt->at(idxNxt) == idxD) return false;
	}

	//check for the correct number of pions and kaons
	int npi(0), nk(0);
	if(TMath::Abs(gen_pid->at(idx0)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx0)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx1)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx1)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx2)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx2)) == 321.) ++nk;

	if(nk!=1 || npi!=2) return false;

	//if ordered then require the kaon to be in idx0
	if(ordered && TMath::Abs(gen_pid->at(idx0)) != 321.) return false;

	return true;
}

bool charmEfficiencies::checkForDs(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered) {
	//first check indices are in range
	if(idxD<0 || idxD>=static_cast<int>(gen_pid->size())) return false;
	if(idx0<0 || idx0>=static_cast<int>(gen_pid->size())) return false;
	if(idx1<0 || idx1>=static_cast<int>(gen_pid->size())) return false;
	if(idx2<0 || idx2>=static_cast<int>(gen_pid->size())) return false;

	//now check we have a Ds
	if(TMath::Abs(gen_pid->at(idxD))!=431.) return false;

	//now check the hierarchy
	if(gen_idx_prnt->at(idx0) != idxD) return false;
	if(gen_idx_prnt->at(idx1) != idxD) return false;
	if(gen_idx_prnt->at(idx2) != idxD) return false;

	//check which type of quark produced the D
	if(origin>0) {
		int idxPrnt = idxD;
		while(gen_idx_prnt->at(idxPrnt)>=0 && gen_idx_prnt->at(idxPrnt)<static_cast<int>(gen_pid->size())) {
			idxPrnt = gen_idx_prnt->at(idxPrnt);
		}
		if(TMath::Abs(gen_prnt_pid->at(idxPrnt)) != origin) return false;
	}

	//if we have a next particle then check it is either a photon (assume radiative) or unrelated
	if(idxNxt>=0 && idxNxt<static_cast<int>(gen_pid->size())) {
		if(gen_pid->at(idxNxt) != 22. && gen_idx_prnt->at(idxNxt) == idxD) return false;
	}

	//check for the correct number of pions and kaons (from phi)
	int npi(0), nk(0);
	int k1(-1), k2(-1);
	if(TMath::Abs(gen_pid->at(idx0)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx0)) == 321. /*&& gen_res_pid->at(idx0) == 333.*/) {
		++nk;
		k1=idx0;
	}
	if(TMath::Abs(gen_pid->at(idx1)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx1)) == 321. /*&& gen_res_pid->at(idx1) == 333.*/) {
		++nk;
		if(k1==-1) k1=idx1;
		else k2=idx1;
	}
	if(TMath::Abs(gen_pid->at(idx2)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx2)) == 321. /*&& gen_res_pid->at(idx2) == 333.*/) {
		++nk;
		k2=idx2;
	}

	if(nk!=2 || npi!=1) return false;

	//resonance not set for events created by DDALITZ model so instead use a mass window
	TLorentzVector phi(gen_px->at(k1)+gen_px->at(k2),
			   gen_py->at(k1)+gen_py->at(k2),
			   gen_pz->at(k1)+gen_pz->at(k2),
			   gen_e->at( k1)+gen_e->at( k2));
	if(phi.M()<990.||phi.M()>1050.) return false;

	//if ordered then require the pion to be in idx2
	if(ordered && TMath::Abs(gen_pid->at(idx2)) != 211.) return false;

	return true;
}

bool charmEfficiencies::checkForLc(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered) {
	//first check indices are in range
	if(idxD<0 || idxD>=static_cast<int>(gen_pid->size())) return false;
	if(idx0<0 || idx0>=static_cast<int>(gen_pid->size())) return false;
	if(idx1<0 || idx1>=static_cast<int>(gen_pid->size())) return false;
	if(idx2<0 || idx2>=static_cast<int>(gen_pid->size())) return false;

	//now check we have a Lc
	if(TMath::Abs(gen_pid->at(idxD))!=4122.) return false;

	//now check the hierarchy
	if(gen_idx_prnt->at(idx0) != idxD) return false;
	if(gen_idx_prnt->at(idx1) != idxD) return false;
	if(gen_idx_prnt->at(idx2) != idxD) return false;

	//check which type of quark produced the D
	if(origin>0) {
		int idxPrnt = idxD;
		while(gen_idx_prnt->at(idxPrnt)>=0 && gen_idx_prnt->at(idxPrnt)<static_cast<int>(gen_pid->size())) {
			idxPrnt = gen_idx_prnt->at(idxPrnt);
		}
		if(TMath::Abs(gen_prnt_pid->at(idxPrnt)) != origin) return false;
	}


	//if we have a next particle then check it is either a photon (assume radiative) or unrelated
	if(idxNxt>=0 && idxNxt<static_cast<int>(gen_pid->size())) {
		if(gen_pid->at(idxNxt) != 22. && gen_idx_prnt->at(idxNxt) == idxD) return false;
	}

	//check for the correct number of pions and kaons
	int npi(0), nk(0), np(0);
	if(TMath::Abs(gen_pid->at(idx0)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx0)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx0)) == 2212.) ++np;
	if(TMath::Abs(gen_pid->at(idx1)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx1)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx1)) == 2212.) ++np;
	if(TMath::Abs(gen_pid->at(idx2)) == 211.) ++npi;
	if(TMath::Abs(gen_pid->at(idx2)) == 321.) ++nk;
	if(TMath::Abs(gen_pid->at(idx2)) == 2212.) ++np;

	if(nk!=1 || npi!=1 || np!=1) return false;

	//if ordered then require the pion to be in idx2 and kaon in idx1
	if(ordered && TMath::Abs(gen_pid->at(idx2)) != 211.) return false;
	if(ordered && TMath::Abs(gen_pid->at(idx1)) != 321.) return false;

	return true;
}

bool charmEfficiencies::inLHCb(int idx) {//check "DaughtersInLHCb"
	TVector3 p(gen_px->at(idx),
			gen_py->at(idx),
			gen_pz->at(idx));

	return (p.Theta()>0.010 && p.Theta()<0.400);
}

double charmEfficiencies::getProbPi(int idx) {
	if(!pionPIDHist) return 0.;
	double PX = trk_px->at(idx);
	double PY = trk_py->at(idx);
	double PZ = trk_pz->at(idx);
	double P  = TMath::Sqrt(PX*PX + PY*PY + PZ*PZ);
	double PT = TMath::Sqrt(PX*PX + PY*PY);
	int bin = pionPIDHist->FindBin(P,PT);
	return pionPIDHist->GetBinContent(bin);
}

double charmEfficiencies::getProbK(int idx) {
	if(!kaonPIDHist) return 0.;
	double PX = trk_px->at(idx);
	double PY = trk_py->at(idx);
	double PZ = trk_pz->at(idx);
	double P  = TMath::Sqrt(PX*PX + PY*PY + PZ*PZ);
	double PT = TMath::Sqrt(PX*PX + PY*PY);
	int bin = kaonPIDHist->FindBin(P,PT);
	return kaonPIDHist->GetBinContent(bin);
}

double charmEfficiencies::getProbP(int idx) {
	if(!protonPIDHist) return 0.;
	double PX = trk_px->at(idx);
	double PY = trk_py->at(idx);
	double PZ = trk_pz->at(idx);
	double P  = TMath::Sqrt(PX*PX + PY*PY + PZ*PZ);
	double PT = TMath::Sqrt(PX*PX + PY*PY);
	int bin = protonPIDHist->FindBin(P,PT);
	return protonPIDHist->GetBinContent(bin);
}

void charmEfficiencies::Loop()
{
	gStyle->SetOptStat(0);

	if (fChain == 0) return;

	TH1D etaD("etaD","",100,0.,6.);
	TH1D ptD("ptD","",100,0.,10000.);

	TH1D hPhi("hPhi","",100,0.,2000.);
	
	int npt(25), neta(6);
	double ptBins[] = {0.,1000.,1500.,2000.,2250.,2500.,2750.,3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,7000.,8000.,9000.,10000.,11000.,12000.,13000.,14000.,15000.,17500.,20000.};
	double etaBins[] = {2.,2.5,3.,3.5,4.,4.5,5.};

	TH2D pidD0("pidD0","",npt,ptBins,neta,etaBins);
	TH2D numD0("numD0","",npt,ptBins,neta,etaBins);
	TH2D denomD0("denomD0","",npt,ptBins,neta,etaBins);

	TH2D pidDp("pidDp","",npt,ptBins,neta,etaBins);
	TH2D numDp("numDp","",npt,ptBins,neta,etaBins);
	TH2D denomDp("denomDp","",npt,ptBins,neta,etaBins);

	TH2D pidDs("pidDs","",npt,ptBins,neta,etaBins);
	TH2D numDs("numDs","",npt,ptBins,neta,etaBins);
	TH2D denomDs("denomDs","",npt,ptBins,neta,etaBins);

	TH2D pidLc("pidLc","",npt,ptBins,neta,etaBins);
	TH2D numLc("numLc","",npt,ptBins,neta,etaBins);
	TH2D denomLc("denomLc","",npt,ptBins,neta,etaBins);

	pidD0.Sumw2();
	numD0.Sumw2();
	denomD0.Sumw2();

	pidDp.Sumw2();
	numDp.Sumw2();
	denomDp.Sumw2();

	pidDs.Sumw2();
	numDs.Sumw2();
	denomDs.Sumw2();

	pidLc.Sumw2();
	numLc.Sumw2();
	denomLc.Sumw2();

	int genD0(0), genDp(0), genDs(0), genLc(0);
	double recoD0(0), recoDp(0), recoDs(0), recoLc(0);
	double recotrueD0(0), recotrueDp(0), recotrueDs(0), recotrueLc(0);
	double selD0(0), selDp(0), selDs(0), selLc(0);

	Long64_t nentries = fChain->GetEntries();
	if(nmax_>=0 && nmax_<nentries) nentries=nmax_;

	Long64_t nbytes = 0, nb = 0;
	boost::progress_display progress( nentries );
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		++progress;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		//loop over generated D's
		for(uint iGen=0; iGen<gen_pid->size()-2; ++iGen) {
			if(checkForD0(iGen,iGen+1,iGen+2,iGen+3,origin_)) {
				if(inLHCb(iGen+1) && inLHCb(iGen+2)) {
					TLorentzVector p4D(gen_px->at(iGen),gen_py->at(iGen),gen_pz->at(iGen),gen_e->at(iGen));
					ptD.Fill(p4D.Pt());
					etaD.Fill(p4D.Rapidity());
					denomD0.Fill(p4D.Pt(),p4D.Rapidity());
					++genD0;
				}
			}
			if(checkForDp(iGen,iGen+1,iGen+2,iGen+3,iGen+4,origin_)) {
				if(inLHCb(iGen+1) && inLHCb(iGen+2) && inLHCb(iGen+3)) {
					TLorentzVector p4D(gen_px->at(iGen),gen_py->at(iGen),gen_pz->at(iGen),gen_e->at(iGen));
					ptD.Fill(p4D.Pt());
					etaD.Fill(p4D.Rapidity());
					denomDp.Fill(p4D.Pt(),p4D.Rapidity());
					++genDp;
				}
			}
			if(checkForDs(iGen,iGen+1,iGen+2,iGen+3,iGen+4,origin_)) {
				if(inLHCb(iGen+1) && inLHCb(iGen+2) && inLHCb(iGen+3)) {
					TLorentzVector p4D(gen_px->at(iGen),gen_py->at(iGen),gen_pz->at(iGen),gen_e->at(iGen));
					ptD.Fill(p4D.Pt());
					etaD.Fill(p4D.Rapidity());
					denomDs.Fill(p4D.Pt(),p4D.Rapidity());
					++genDs;
				}
			}
			if(checkForLc(iGen,iGen+1,iGen+2,iGen+3,iGen+4,origin_)) {
				if(inLHCb(iGen+1) && inLHCb(iGen+2) && inLHCb(iGen+3)) {
					TLorentzVector p4D(gen_px->at(iGen),gen_py->at(iGen),gen_pz->at(iGen),gen_e->at(iGen));
					ptD.Fill(p4D.Pt());
					etaD.Fill(p4D.Rapidity());
					denomLc.Fill(p4D.Pt(),p4D.Rapidity());
					++genLc;
				}
			}
		}//end loop over generated D's

		//loop over reconstructed D0's
		recoD0+=d0_m->size();
		for(uint iD=0; iD<d0_m->size(); ++iD) {
			int trk0 = d0_idx_trk0->at(iD);
			int trk1 = d0_idx_trk1->at(iD);
			int gen0 = trk_idx_gen->at(trk0);
			int gen1 = trk_idx_gen->at(trk1);
			if(gen0<0 || gen1<0) continue;
			int genD = TMath::Min(gen0,gen1)-1;
			if(genD<0) continue;
			int genNxt = genD+3;
			if(checkForD0(genD,gen0,gen1,genNxt,origin_,true)) {
				if(inLHCb(gen0) && inLHCb(gen1)) {
					double probK  = getProbK(trk0);
					double probPi = getProbPi(trk1);
					TLorentzVector p4D(gen_px->at(genD),gen_py->at(genD),gen_pz->at(genD),gen_e->at(genD));
					numD0.Fill(p4D.Pt(),p4D.Rapidity());
					pidD0.Fill(p4D.Pt(),p4D.Rapidity(),probPi*probK);
					selD0 += probPi*probK;
					++recotrueD0;
				}
			}
		}//end loop over reconstructed D0's

		//loop over reconstructed Dp's
		recoDp+=dp_m->size();
		for(uint iD=0; iD<dp_m->size(); ++iD) {
			int trk0 = dp_idx_trk0->at(iD);
			int trk1 = dp_idx_trk1->at(iD);
			int trk2 = dp_idx_trk2->at(iD);
			int gen0 = trk_idx_gen->at(trk0);
			int gen1 = trk_idx_gen->at(trk1);
			int gen2 = trk_idx_gen->at(trk2);
			if(gen0<0 || gen1<0 || gen2<0) continue;
			int genD = TMath::Min(gen0,TMath::Min(gen1,gen2))-1;
			if(genD<0) continue;
			int genNxt = genD+4;
			if(checkForDp(genD,gen0,gen1,gen2,genNxt,origin_,true)) {
				if(inLHCb(gen0) && inLHCb(gen1) && inLHCb(gen2)) {
					double probK   = getProbK(trk0);
					double probPi1 = getProbPi(trk1);
					double probPi2 = getProbPi(trk2);
					TLorentzVector p4D(gen_px->at(genD),gen_py->at(genD),gen_pz->at(genD),gen_e->at(genD));
					numDp.Fill(p4D.Pt(),p4D.Rapidity());
					pidDp.Fill(p4D.Pt(),p4D.Rapidity(),probPi1*probPi2*probK);
					selDp += probPi1*probPi2*probK;
					++recotrueDp;
				}
			}
		}//end loop over reconstructed Dp's

		//loop over reconstructed Ds's
		recoDs+=ds_m->size();
		for(uint iD=0; iD<ds_m->size(); ++iD) {
			int trk0 = ds_idx_trk0->at(iD);
			int trk1 = ds_idx_trk1->at(iD);
			int trk2 = ds_idx_trk2->at(iD);
			int gen0 = trk_idx_gen->at(trk0);
			int gen1 = trk_idx_gen->at(trk1);
			int gen2 = trk_idx_gen->at(trk2);
			if(gen0<0 || gen1<0 || gen2<0) continue;
			int genD = TMath::Min(gen0,TMath::Min(gen1,gen2))-1;
			if(genD<0) continue;
			int genNxt = genD+4;
			if(checkForDs(genD,gen0,gen1,gen2,genNxt,origin_,true)) {
				if(inLHCb(gen0) && inLHCb(gen1) && inLHCb(gen2)) {
					TLorentzVector k1, k2;
					k1.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),496.677);
					k2.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),496.677);
					TLorentzVector phi = k1+k2;
					hPhi.Fill(phi.M());
					if(phi.M()<990.||phi.M()>1050.) continue;
					double probK1 = getProbK(trk0);
					double probK2 = getProbK(trk1);
					double probPi = getProbPi(trk2);
					TLorentzVector p4D(gen_px->at(genD),gen_py->at(genD),gen_pz->at(genD),gen_e->at(genD));
					numDs.Fill(p4D.Pt(),p4D.Rapidity());
					pidDs.Fill(p4D.Pt(),p4D.Rapidity(),probPi*probK1*probK2);
					selDs += probPi*probK1*probK2;
					++recotrueDs;
				}
			}

		}//end loop over reconstructed Lc's
		//loop over reconstructed Lc's
		recoLc+=lc_m->size();
		for(uint iD=0; iD<lc_m->size(); ++iD) {
			int trk0 = lc_idx_trk0->at(iD);
			int trk1 = lc_idx_trk1->at(iD);
			int trk2 = lc_idx_trk2->at(iD);
			int gen0 = trk_idx_gen->at(trk0);
			int gen1 = trk_idx_gen->at(trk1);
			int gen2 = trk_idx_gen->at(trk2);
			if(gen0<0 || gen1<0 || gen2<0) continue;
			int genD = TMath::Min(gen0,TMath::Min(gen1,gen2))-1;
			if(genD<0) continue;
			int genNxt = genD+4;
			if(checkForLc(genD,gen0,gen1,gen2,genNxt,origin_,true)) {
				if(inLHCb(gen0) && inLHCb(gen1) && inLHCb(gen2)) {
					double probP  = getProbP(trk0);
					double probK  = getProbK(trk1);
					double probPi = getProbPi(trk2);
					TLorentzVector p4D(gen_px->at(genD),gen_py->at(genD),gen_pz->at(genD),gen_e->at(genD));
					numLc.Fill(p4D.Pt(),p4D.Rapidity());
					pidLc.Fill(p4D.Pt(),p4D.Rapidity(),probPi*probK*probP);
					selLc += probPi*probK*probP;
					++recotrueLc;
				}
			}
		}//end loop over reconstructed Ds's
	}

	TEfficiency effD0(numD0,denomD0);
	TEfficiency effDp(numDp,denomDp);
	TEfficiency effDs(numDs,denomDs);
	TEfficiency effLc(numLc,denomLc);

	TString typeStr; typeStr+=type_; typeStr+="_"; typeStr+=origin_;
	TCanvas c;
	ptD.Draw();
	c.SaveAs("ptD"+typeStr+".pdf");
	etaD.Draw();
	c.SaveAs("etaD"+typeStr+".pdf");

	hPhi.Draw();
	c.SaveAs("m_phi.pdf");

	pidD0.Divide(&numD0);
	pidD0.Draw("colz");
	c.SaveAs("pidD0"+typeStr+".pdf");

	pidDp.Divide(&numDp);
	pidDp.Draw("colz");
	c.SaveAs("pidDp"+typeStr+".pdf");

	pidDs.Divide(&numDs);
	pidDs.Draw("colz");
	c.SaveAs("pidDs"+typeStr+".pdf");

	pidLc.Divide(&numLc);
	pidLc.Draw("colz");
	c.SaveAs("pidLc"+typeStr+".pdf");

	numD0.Divide(&denomD0);
	numD0.Draw("colz");
	c.SaveAs("effD0"+typeStr+".pdf");

	numDp.Divide(&denomDp);
	numDp.Draw("colz");
	c.SaveAs("effDp"+typeStr+".pdf");

	numDs.Divide(&denomDs);
	numDs.Draw("colz");
	c.SaveAs("effDs"+typeStr+".pdf");

	numLc.Divide(&denomLc);
	numLc.Draw("colz");
	c.SaveAs("effLc"+typeStr+".pdf");

	denomD0.Draw("colz");
	c.SaveAs("statsD0"+typeStr+".pdf");

	denomDp.Draw("colz");
	c.SaveAs("statsDp"+typeStr+".pdf");

	denomDs.Draw("colz");
	c.SaveAs("statsDs"+typeStr+".pdf");

	denomLc.Draw("colz");
	c.SaveAs("statsLc"+typeStr+".pdf");

	std::cout << genD0 << "\t" << genDp << "\t" << genDs << "\t" << genLc << std::endl;
	std::cout << recoD0 << "\t" << recoDp << "\t" << recoDs <<"\t" << recoLc << std::endl;
	std::cout << recotrueD0 << "\t" << recotrueDp << "\t" << recotrueDs <<"\t" << recotrueLc << std::endl;
	std::cout << selD0 << "\t" << selDp << "\t" << selDs <<"\t" << selLc << std::endl;

	TFile* fout = TFile::Open("efficiencies"+typeStr+".root","RECREATE");
	ptD.SetName("pt");
	ptD.Write();
	etaD.SetName("eta");
	etaD.Write();
	effD0.SetName("effD0");
	effD0.Write();
	effDp.SetName("effDp");
	effDp.Write();
	effDs.SetName("effDs");
	effDs.Write();
	effLc.SetName("effLc");
	effLc.Write();
	numD0.SetName("efficiencyD0");
	numD0.Write();
	numDp.SetName("efficiencyDp");
	numDp.Write();
	numDs.SetName("efficiencyDs");
	numDs.Write();
	numLc.SetName("efficiencyLc");
	numLc.Write();
	pidD0.SetName("pidD0");
	pidD0.Write();
	pidDp.SetName("pidDp");
	pidDp.Write();
	pidDs.SetName("pidDs");
	pidDs.Write();
	pidLc.SetName("pidLc");
	pidLc.Write();
	denomD0.SetName("statsD0");
	denomD0.Write();
	denomDp.SetName("statsDp");
	denomDp.Write();
	denomDs.SetName("statsDs");
	denomDs.Write();
	denomLc.SetName("statsLc");
	denomLc.Write();
	fout->Close();
}

int main(int argc, char** argv) {
	int type(0), origin(0), nmax(-1);
	if(argc>1) type = atoi(argv[1]);
	if(argc>2) origin = atoi(argv[2]);
	if(argc>3) nmax = atoi(argv[3]);

	charmEfficiencies c(type,origin,nmax);

	if(origin>5) {
		c.setOrigin(4);
		c.Loop();
		c.setOrigin(5);
	}
		c.Loop();
}
