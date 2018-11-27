#define makeCombShape_cxx
#include "makeCombShape.h"

#include <iostream>

#include <TCanvas.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "TMatrixD.h"
#include <TRandom3.h>
#include <TStyle.h>
#include "TVectorT.h"

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
// Use this one! Verified against the second method and documented here:
// http://geomalgorithms.com/a07-_distance.html
double makeCombShape::calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
		const TVector3 &pb, const TVector3 &db) {

	//lines defined by pa + s*da and pb + t*db
	//vector connecting pa and pb
	TVector3 w0 = pa-pb;

	//define some dot products to simplify the maths
	double a = da.Mag2();
	double b = da.Dot(db);
	double c = db.Mag2();
	double d = da.Dot(w0);
	double e = db.Dot(w0);

	//calculate coefficients of closest point
	double sc(0), tc(0);

	if(a*c - b*b == 0) {
		//lines are parallel - set one coefficient to zero and solve for t'other
		sc = 0;
		tc = e/c;
	} else {
		//general case - see http://geomalgorithms.com/a07-_distance.html
		sc=(b*e - c*d)/(a*c - b*b);
		tc=(a*e - b*d)/(a*c - b*b);
	}

	//points on lines with shortest distance
	TVector3 Pc = pa + sc*da;
	TVector3 Qc = pb + tc*db;

	//give vertex at centre point of the connecting line
	v = Pc+Qc;
	v *= 0.5;

	//return separation at closest point
	return (Pc-Qc).Mag();
}

//bool makeCombShape::checkTruth(int idxD, int needPID, int needPi, int needK, int needP, int needMu, int needOppPi, int needOppK, int needOppP) {
//	int ng = gen_pid->size();
//	int pid = gen_pid->at(idxD);
//	if(TMath::Abs(pid)!=needPID) return false;
//
//	int foundPi(0), foundK(0), foundP(0), foundMu(0), foundOppPi(0), foundOppK(0), foundOppP(0);
//	if(idxD<0 || idxD>=ng) return false;
//	for(int g=idxD+1; g<ng; ++g) {
//		if(gen_idx_prnt->at(g) != idxD) continue;
//		int childPID = gen_pid->at(g);
//		switch(TMath::Abs(childPID)) {
//			case 22:
//				if(gen_res_pid->at(g)!=-1) return false;
//				break;
//			case 13:
//				++foundMu;
//				break;
//			case 211:
//				++foundPi;
//				if(pid*childPID<0) ++ foundOppPi;
//				break;
//			case 321:
//				++foundK;
//				if(pid*childPID<0) ++ foundOppK;
//				break;
//			case 2212:
//				++foundP;
//				if(pid*childPID<0) ++ foundOppP;
//				break;
//			default:
//				return false;
//		}
//	}
//	if(foundPi!=needPi || foundK!=needK || foundP!=needP || foundMu!=needMu) return false;
//	if(needOppPi>-1 && foundOppPi!=needOppPi) return false;
//	if(needOppK>-1 && foundOppK!=needOppK) return false;
//	if(needOppP>-1 && foundOppP!=needOppP) return false;
//
//	return true;
//}
bool makeCombShape::checkTruth(int idxD, std::vector<int>& idcs, int needPID, int needPi, int needK, int needP, int needMu, int needOppPi, int needOppK, int needOppP, bool allowPart) {
	TString debugStr;
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(needPID!=0 && TMath::Abs(pid)!=needPID) return false;
	debugStr+=pid;
	debugStr+=" ->";

	int foundPi(0), foundK(0), foundP(0), foundMu(0), foundOppPi(0), foundOppK(0), foundOppP(0);
	if(idxD<0 || idxD>=ng) return false;
	for(int g=idxD+1; g<ng; ++g) {
		if(gen_idx_prnt->at(g) != idxD) continue;
		int childPID = gen_pid->at(g);
		switch(TMath::Abs(childPID)) {
			case 22:
				if(allowPart) continue;
				if(gen_res_pid->at(g)!=-1) return false;
				continue;
			case 13:
				++foundMu;
				break;
			case 211:
				++foundPi;
				if(pid*childPID<0) ++ foundOppPi;
				break;
			case 321:
				++foundK;
				if(pid*childPID<0) ++ foundOppK;
				break;
			case 2212:
				++foundP;
				if(pid*childPID<0) ++ foundOppP;
				break;
			default:
				if(allowPart) continue;
				return false;
		}
		idcs.push_back(g);
		debugStr+=" ";
		debugStr+=childPID;
	}
	//std::cout << std::endl << debugStr << std::endl; //TODO debug
	if(needPi>-1 && foundPi!=needPi) return false;
	if(needK>-1  && foundK!=needK) return false;
	if(needP>-1  && foundP!=needP) return false;
	if(needMu>-1 && foundMu!=needMu) return false;
	if(needOppPi>-1 && foundOppPi!=needOppPi) return false;
	if(needOppK>-1 && foundOppK!=needOppK) return false;
	if(needOppP>-1 && foundOppP!=needOppP) return false;
	//std::cout << debugStr << std::endl; //TODO debug

	return true;
}

bool makeCombShape::checkD0Truth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,421,1,1,0,0);
}

bool makeCombShape::checkD0Truth(int idxK, int idxPi, int& idxD, bool ordered) {
	//first check PIDs
	if(idxK<0 || idxPi<0) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxK))!=321) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxPi))!=211) return false;

	//if not ordered then apply a less stringent requirements
	if(!ordered && (TMath::Abs(gen_pid->at(idxK))==321)+(TMath::Abs(gen_pid->at(idxPi))==321)!=1) return false;
	if(!ordered && (TMath::Abs(gen_pid->at(idxK))==211)+(TMath::Abs(gen_pid->at(idxPi))==211)!=1) return false;

	//check shared parent
	if(gen_idx_prnt->at(idxK)<0 || gen_idx_prnt->at(idxK) != gen_idx_prnt->at(idxPi)) return false;
	idxD = gen_idx_prnt->at(idxK);

	//now check parent is truth
	return checkD0Truth(idxD);
}

bool makeCombShape::checkDpTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,411,2,1,0,0);
}

bool makeCombShape::checkDpTruth(int idxK, int idxPi1, int idxPi2, int& idxD, bool ordered) {
	//first check PIDs
	if(idxK<0 || idxPi1<0 || idxPi2<0) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxK))!=321) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxPi1))!=211) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxPi2))!=211) return false;

	//if not ordered then apply a less stringent requirements
	if(!ordered && (TMath::Abs(gen_pid->at(idxK))==321)+(TMath::Abs(gen_pid->at(idxPi1))==321)+(TMath::Abs(gen_pid->at(idxPi2))==321)!=1) return false;
	if(!ordered && (TMath::Abs(gen_pid->at(idxK))==211)+(TMath::Abs(gen_pid->at(idxPi1))==211)+(TMath::Abs(gen_pid->at(idxPi2))==211)!=2) return false;

	//check shared parent
	if(gen_idx_prnt->at(idxK)<0 || gen_idx_prnt->at(idxK) != gen_idx_prnt->at(idxPi1) || gen_idx_prnt->at(idxK) != gen_idx_prnt->at(idxPi2)) return false;
	idxD = gen_idx_prnt->at(idxK);

	//now check parent is truth
	return checkDpTruth(idxD);
}

bool makeCombShape::checkDsTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,431,1,2,0,0);
}

bool makeCombShape::checkDsTruth(int idxK1, int idxK2, int idxPi, int& idxD, bool ordered) {
	//first check PIDs
	if(idxK1<0 || idxK2<0 || idxPi<0) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxK1))!=321) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxK2))!=321) return false;
	if(ordered && TMath::Abs(gen_pid->at(idxPi))!=211) return false;

	//if not ordered then apply a less stringent requirements
	if(!ordered && (TMath::Abs(gen_pid->at(idxK1))==321)+(TMath::Abs(gen_pid->at(idxK2))==321)+(TMath::Abs(gen_pid->at(idxPi))==321)!=2) return false;
	if(!ordered && (TMath::Abs(gen_pid->at(idxK1))==211)+(TMath::Abs(gen_pid->at(idxK2))==211)+(TMath::Abs(gen_pid->at(idxPi))==211)!=1) return false;

	//check shared parent
	if(gen_idx_prnt->at(idxK1)<0 || gen_idx_prnt->at(idxK1) != gen_idx_prnt->at(idxK2) || gen_idx_prnt->at(idxK1) != gen_idx_prnt->at(idxPi)) return false;
	idxD = gen_idx_prnt->at(idxK1);

	//now check parent is truth
	return checkDsTruth(idxD);
}

bool makeCombShape::checkJpsiTruth(int idxJpsi) {
	std::vector<int> idcs;
	return checkTruth(idxJpsi,idcs,443,0,0,0,2);
}

bool makeCombShape::checkJpsiTruth(int idxMu1, int idxMu2, int& idxJpsi) {
	//first check PIDs
	if(idxMu1<0 || idxMu2<0) return false;
	if(TMath::Abs(gen_pid->at(idxMu1))!=13) return false;
	if(TMath::Abs(gen_pid->at(idxMu2))!=13) return false;

	//check shared parent
	if(gen_idx_prnt->at(idxMu1)<0 || gen_idx_prnt->at(idxMu1) != gen_idx_prnt->at(idxMu2)) return false;
	idxJpsi = gen_idx_prnt->at(idxMu1);

	//now check parent is truth
	return checkJpsiTruth(idxJpsi);
}

//check for any long-lived hadron with at least two charged track children
bool makeCombShape::checkSVTruth(int idxD) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLived(pid)) return false;

	int foundTrk(0);
	if(idxD<0 || idxD>=ng) return false;
	for(int g=idxD+1; g<ng; ++g) {
		//search up the decay tree to find either our hadron or no parent
		int idxPrnt = gen_idx_prnt->at(g);
		while(idxPrnt>-1 && idxPrnt!= idxD) idxPrnt = gen_idx_prnt->at(idxPrnt);
		if(idxPrnt != idxD) continue;
		int childPID = gen_pid->at(g);
		switch(TMath::Abs(childPID)) {
			case 11:
			case 13:
			case 211:
			case 321:
			case 2212:
				++foundTrk;
				break;
			default:
				break;
		}
	}
	if(foundTrk<2) return false;

	return true;
}

//check for any long-lived c-hadron with at least two charged track children
bool makeCombShape::checkDSVTruth(int idxD) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLivedC(pid)) return false;

	int foundTrk(0);
	if(idxD<0 || idxD>=ng) return false;
	for(int g=idxD+1; g<ng; ++g) {
		//search up the decay tree to find either our c-hadron or no parent
		int idxPrnt = gen_idx_prnt->at(g);
		while(idxPrnt>-1 && idxPrnt!= idxD) idxPrnt = gen_idx_prnt->at(idxPrnt);
		if(idxPrnt != idxD) continue;
		int childPID = gen_pid->at(g);
		switch(TMath::Abs(childPID)) {
			case 11:
			case 13:
			case 211:
			case 321:
			case 2212:
				++foundTrk;
				break;
			default:
				break;
		}
	}
	if(foundTrk<2) return false;

	return true;
}

//check for any long-lived b-hadron with at least two charged track children
bool makeCombShape::checkBSVTruth(int idxD) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLivedB(pid)) return false;

	int foundTrk(0);
	if(idxD<0 || idxD>=ng) return false;
	for(int g=idxD+1; g<ng; ++g) {
		//search up the decay tree to find either our b-hadron or no parent
		int idxPrnt = gen_idx_prnt->at(g);
		while(idxPrnt>-1 && idxPrnt!= idxD) idxPrnt = gen_idx_prnt->at(idxPrnt);
		if(idxPrnt != idxD) continue;
		int childPID = gen_pid->at(g);
		switch(TMath::Abs(childPID)) {
			case 11:
			case 13:
			case 211:
			case 321:
			case 2212:
				++foundTrk;
				break;
			default:
				break;
		}
	}
	if(foundTrk<2) return false;

	return true;
}

bool makeCombShape::checkFromB(int idxD) {
	int ng = gen_pid->size();
	if(idxD<0 || idxD>=ng) return false;

	int idxPrnt = idxD;
	while(gen_idx_prnt->at(idxPrnt)>-1) {
//		std::cout << gen_pid->at(idxPrnt) << "<-";//TODO
//		if(longLivedB(gen_pid->at(idxPrnt))) return true;
		idxPrnt = gen_idx_prnt->at(idxPrnt);
	}

	int prntPID = gen_prnt_pid->at(idxPrnt);
	if(prntPID==-1) prntPID=gen_pid->at(idxPrnt);
//	std::cout << prntPID << std::endl;//TODO

	int origin = getFlavour(prntPID);
//	std::cout << origin << std::endl;//TODO
	if(origin==5) return true;
	return false;
}

int makeCombShape::getFlavour(int pid) {
	pid = TMath::Abs(pid);
	
	if(pid<100) return pid;

	int pid1=pid%10000/1000;
	int pid2=pid%1000/100;
	int pid3=pid%100/10;

	return TMath::Max(pid1,TMath::Max(pid2,pid3));
}

bool makeCombShape::longLivedB(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==511 || pid==521 || pid==531 || pid==5122 || pid==5132 || pid==5232 || pid==5332 );
}

bool makeCombShape::longLivedC(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==411 || pid==421 || pid==431 || pid==4122 || pid==4132 || pid==4232 || pid==4332 );
}

bool makeCombShape::longLivedS(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==130 || pid==310 || pid==3122 || pid==3112 || pid==3222 || pid==3312 || pid==3322 );
}

bool makeCombShape::longLived(int pid) {
	return (longLivedB(pid) || longLivedC(pid) || longLivedS(pid));
}

void makeCombShape::fillSVCands(int j, int t)
{
	TLorentzVector p4j, p4t;
	if(j>-1) p4j.SetPxPyPzE(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	if(t>-1) p4t.SetPxPyPzE(jet_px->at(t),jet_py->at(t),jet_pz->at(t),jet_e->at(t));

	int ipv(0);// currently veto any events with more than one PV - if we unveto then at least this won't break for jet data/MC
	if(j>-1) ipv = jet_idx_pvr->at(j);
	TVector3 pv(pvr_x->at(ipv),pvr_y->at(ipv),pvr_z->at(ipv)); 

	for(uint iD=0; iD<d0_m->size(); ++iD) {
		TLorentzVector p4D(d0_px->at(iD),
				d0_py->at(iD),
				d0_pz->at(iD),
				d0_e->at( iD));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk0, p4Dtrk1;

		int trk0=d0_idx_trk0->at(iD);
		int trk1=d0_idx_trk1->at(iD);

		nj = (trk_idx_jet->at(trk0)==j) + (trk_idx_jet->at(trk1)==j);
		p4Dtrk0.SetPxPyPzE(trk_px->at(trk0),
				   trk_py->at(trk0),
				   trk_pz->at(trk0),
				   trk_e ->at(trk0));
		if(j>-1) maxdr = TMath::Max(p4Dtrk0.DeltaR(p4j), maxdr);
		p4Dtrk1.SetPxPyPzE(trk_px->at(trk1),
				   trk_py->at(trk1),
				   trk_pz->at(trk1),
				   trk_e ->at(trk1));
		if(j>-1) maxdr = TMath::Max(p4Dtrk1.DeltaR(p4j), maxdr);

		//if(nj == 0) continue;//require at least one track in probe jet

		TVector3 pd(d0_px->at(iD),d0_py->at(iD),d0_pz->at(iD));
		TVector3 xd(d0_x ->at(iD) - pvr_x->at(d0_idx_pvr->at(iD)),
		            d0_y ->at(iD) - pvr_y->at(d0_idx_pvr->at(iD)),
		            d0_z ->at(iD) - pvr_z->at(d0_idx_pvr->at(iD)));
		
		//Tighter PID, pT and ghost requirements 
		//kaon PID reversed for background sample
		if(trk_pnn_k->at(trk1)>0.2 || trk_pnn_k->at( trk0)>0.2) continue;
		if(!(pd.Pt()>2000)) continue;
		if(!(p4Dtrk0.Pt()>500)) continue;
		if(!(p4Dtrk1.Pt()>500)) continue;
		if(!(trk_prb_ghost->at(trk0)<0.2 && trk_prb_ghost->at(trk1)<0.2)) continue;

		hist_->Fill(d0_m->at(iD),TMath::Log(d0_ip_chi2->at(iD)));
		//std::cout << d0_m->at(iD) << "\t" << TMath::Log(d0_ip_chi2->at(iD)) << std::endl;
	}
}

//function to identify the probe (j) and tag (t) jets in a dijet data event
bool makeCombShape::tagSVJet(int& j, int& t) {
	std::vector<int> trigs;
	std::vector<int> tags;
	std::vector<int> probes;
	std::vector<double> dRs;
	std::vector<double> pTAs;

	bool trigFound(false);

	for(unsigned int trig=0; trig<evt_dec->size(); ++trig) {
		if(evt_dec->at(trig)!=1 &&
		   evt_dec->at(trig)!=6 &&
		   evt_dec->at(trig)!=11) continue;
		
		if(evt_j1_idx->at(trig)<0 || evt_j2_idx->at(trig)<0) continue;
		if(evt_j1_idx->at(trig)==evt_j2_idx->at(trig)) {
			++tagCounts[4];
			continue;
		}

		++tagCounts[2];
		if(trigFound) ++tagCounts[3];
		trigFound=true;

		bool tag1 = (evt_j1_idx->at(trig)>=0 && (evt_j1_nsv->at(trig)>0));
		bool tag2 = (evt_j2_idx->at(trig)>=0 && (evt_j2_nsv->at(trig)>0));

		//tag jet must pass L0 and HLT1 triggers
		if(tag1) {
			tag1 &= (jet_l0_hadron_tos  ->at(evt_j1_idx->at(trig)) ||
				 jet_l0_photon_tos  ->at(evt_j1_idx->at(trig)) ||
				 jet_l0_electron_tos->at(evt_j1_idx->at(trig)) ||
				 jet_l0_muon_tos    ->at(evt_j1_idx->at(trig)) ||
				 jet_l0_dimuon_tos  ->at(evt_j1_idx->at(trig)));

			tag1 &= (jet_hlt1_track_tos  ->at(evt_j1_idx->at(trig)) ||
				 jet_hlt1_ditrack_tos->at(evt_j1_idx->at(trig)));
		}
		if(tag2) {
			tag2 &= (jet_l0_hadron_tos  ->at(evt_j2_idx->at(trig)) ||
				 jet_l0_photon_tos  ->at(evt_j2_idx->at(trig)) ||
				 jet_l0_electron_tos->at(evt_j2_idx->at(trig)) ||
				 jet_l0_muon_tos    ->at(evt_j2_idx->at(trig)) ||
				 jet_l0_dimuon_tos  ->at(evt_j2_idx->at(trig)));

			tag2 &= (jet_hlt1_track_tos  ->at(evt_j2_idx->at(trig)) ||
				 jet_hlt1_ditrack_tos->at(evt_j2_idx->at(trig)));
		}

		if(!tag1 && !tag2) continue;

		unsigned int jtag;
		++tagCounts[5];
		
		if(tag1 && tag2) {
			++tagCounts[6];
			if(gRandom->Integer(2)==0) {
				tag2=false;
				jtag=evt_j1_idx->at(trig);
			} else {
				tag1=false;
				jtag=evt_j2_idx->at(trig);
			}
		} else if(tag1) {
			jtag=evt_j1_idx->at(trig);
		} else {//tag2
			jtag=evt_j2_idx->at(trig);
		}

		TLorentzVector p4tag(jet_px->at(jtag),
				     jet_py->at(jtag),
				     jet_pz->at(jtag),
				     jet_e->at( jtag));

		bool probeFound(false);

		for(unsigned int jprobe=0; jprobe<jet_pt->size(); ++jprobe) {
			if(jtag==jprobe) continue;
			if(jprobe!=evt_j1_idx->at(trig) && jprobe!=evt_j2_idx->at(trig)) continue;//check if the probe jet was part of the trigger //TODO off for now
			++tagCounts[7];

			TLorentzVector p4probe(jet_px->at(jprobe),
					       jet_py->at(jprobe),
					       jet_pz->at(jprobe),
					       jet_e->at( jprobe));
			if(p4probe.Pt()<10000. || p4probe.Eta()<2.5 || p4probe.Eta()>4.0) continue;

			double dR = p4probe.DeltaPhi(p4tag);
			double ptAsym = (jet_pt->at(jprobe)-jet_pt->at(jtag)) /
				(jet_pt->at(jprobe)+jet_pt->at(jtag));
	
			if(TMath::Abs(dR)<2.0) continue;
			if(TMath::Abs(ptAsym)>0.25) continue;

			++tagCounts[8];
			if(probeFound) ++tagCounts[9];
			probeFound=true;
			//std::cout << jtag << "\t" << jprobe << "\t" << evt_j1_idx->at(trig) << "\t" << evt_j2_idx->at(trig) << std::endl;//TODO

			trigs.push_back(trig);
			tags.push_back(jtag);
			probes.push_back(jprobe);
			dRs.push_back(dR);
			pTAs.push_back(ptAsym);

		}//end loop over probe jet
	}//end loop over trigger dec

	if(trigs.size()==0) return false;

	++tagCounts[10];
	if(trigs.size()>1) {
		bool same(false);
		bool switched(false);
		bool different(false);

		for(unsigned int i=1; i<trigs.size(); ++i) {
			if(tags[i]==tags[0] && probes[i]==probes[0]) same=true;
			else if(tags[i]==probes[0] && probes[i]==tags[0]) switched=true;
			else different=true;
		}

		if(switched&&different) ++tagCounts[13];
		else if(switched) ++tagCounts[11];
		else if(different) ++tagCounts[12];
		else if(same) ++tagCounts[14];
	}

	int whichTrig = gRandom->Integer(trigs.size());
	j = probes[whichTrig];
	t = tags[whichTrig];

	return true;
}

void makeCombShape::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	int firstEntry(0), lastEntry(nentries);

	boost::progress_display progress( lastEntry-firstEntry );
	for (Long64_t jentry=firstEntry; jentry<lastEntry;jentry++) {
		++progress;
		fChain->GetEntry(jentry);

		++tagCounts[0];
		if(evt_pvr_n!=1) continue;//keep to deal with topo bug
		++tagCounts[1];

		//select the jets to use
		int jet1(-1), jet2(-1);
		
		switch(tagType_) {
			case NoJets: //if non-jet sample then just fill the SV, D, J/psi containers
				fillSVCands(-1,-1);
				break;
			case JetTag: //if dijet data then tag on one jet and keep the other
				if(tagSVJet(jet1, jet2)) fillSVCands(jet1,jet2);
				break;
			default:
				break;
		}


	}//end loop over entries

	for(int i=0; i<15; ++i) {
		std::cout << tagCounts[i] << "\t";
	}
	std::cout << std::endl;

	fout_->cd();
	hist_->Write();
	fout_->Close();

}

int main(int /*argc*/, char** /*argv*/) {
	makeCombShape a;
	a.Loop();
	return 0;
}
