#define makeNewCalibTuples_cxx
#include "makeNewCalibTuples.h"

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
double makeNewCalibTuples::calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
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

//bool makeNewCalibTuples::checkTruth(int idxD, int needPID, int needPi, int needK, int needP, int needMu, int needOppPi, int needOppK, int needOppP) {
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
bool makeNewCalibTuples::checkTruth(int idxD, std::vector<int>& idcs, int needPID, int needPi, int needK, int needP, int needMu, int needOppPi, int needOppK, int needOppP, bool allowPart) {
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

bool makeNewCalibTuples::checkD0Truth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,421,1,1,0,0);
}

bool makeNewCalibTuples::checkD0Truth(int idxK, int idxPi, int& idxD, bool ordered) {
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

bool makeNewCalibTuples::checkDpTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,411,2,1,0,0);
}

bool makeNewCalibTuples::checkDpTruth(int idxK, int idxPi1, int idxPi2, int& idxD, bool ordered) {
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

bool makeNewCalibTuples::checkDsTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,431,1,2,0,0);
}

bool makeNewCalibTuples::checkDsTruth(int idxK1, int idxK2, int idxPi, int& idxD, bool ordered) {
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

bool makeNewCalibTuples::checkJpsiTruth(int idxJpsi) {
	std::vector<int> idcs;
	return checkTruth(idxJpsi,idcs,443,0,0,0,2);
}

bool makeNewCalibTuples::checkJpsiTruth(int idxMu1, int idxMu2, int& idxJpsi) {
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
bool makeNewCalibTuples::checkSVTruth(int idxD) {
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
bool makeNewCalibTuples::checkDSVTruth(int idxD) {
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
bool makeNewCalibTuples::checkBSVTruth(int idxD) {
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

bool makeNewCalibTuples::checkFromB(int idxD) {
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

int makeNewCalibTuples::getFlavour(int pid) {
	pid = TMath::Abs(pid);
	
	if(pid<100) return pid;

	int pid1=pid%10000/1000;
	int pid2=pid%1000/100;
	int pid3=pid%100/10;

	return TMath::Max(pid1,TMath::Max(pid2,pid3));
}

bool makeNewCalibTuples::longLivedB(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==511 || pid==521 || pid==531 || pid==5122 || pid==5132 || pid==5232 || pid==5332 );
}

bool makeNewCalibTuples::longLivedC(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==411 || pid==421 || pid==431 || pid==4122 || pid==4132 || pid==4232 || pid==4332 );
}

bool makeNewCalibTuples::longLivedS(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==130 || pid==310 || pid==3122 || pid==3112 || pid==3222 || pid==3312 || pid==3322 );
}

bool makeNewCalibTuples::longLived(int pid) {
	return (longLivedB(pid) || longLivedC(pid) || longLivedS(pid));
}

int makeNewCalibTuples::checkRealSV(std::vector<int> indices) {
	//first count the number of our particles descended from each ancestor particle
	std::map<int,int> ancestorCounts;
	for(unsigned int i=0; i<indices.size(); ++i) {
		int idx = indices[i];
		if(idx < 0 || idx > static_cast<int>(gen_pid->size())) continue;
		while(gen_idx_prnt->at(idx)!=-1) {
			idx = gen_idx_prnt->at(idx);
			++ancestorCounts[idx];
		}

	}

	//now find the best according to the following criteria
	// - if flavour_ is 4 or 5 then only accept long lived c- or b-hadrons, respectively
	// - prioritise ancestors with more matched descendents
	// - if equal then prioritise ancestors with higher indices (likely more recent)
	int bestAn(-1);
	int nInBestAn(2); //ignore counts <2
	std::map<int,int>::iterator it = ancestorCounts.begin();

	for( ; it!=ancestorCounts.end(); ++it) {
		int idxAn = (*it).first;
		int nInAn = (*it).second;

		//ignore anything with fewer descendents than best found (or fewer than 2)
		if(nInAn<nInBestAn) continue;

		int apid = TMath::Abs(gen_pid->at(idxAn));
		//std::cout << idxAn << "\t" << apid << "\t" << nInAn << std::endl;//TODO

		switch(flavour_) {
			case 4:
				if(longLivedC(apid)) {
					bestAn = idxAn;
					nInBestAn = nInAn;
				}
				break;
			case 5:
				if(longLivedB(apid)) {
					bestAn = idxAn;
					nInBestAn = nInAn;
				}
				break;
			default:
				if(longLived(apid)) {
					bestAn = idxAn;
					nInBestAn = nInAn;
				}
		}
	}

	//if we didn't find anything then return 0
	if(bestAn<0) return 0;
	return nInBestAn;
}

int makeNewCalibTuples::checkBestD0Cand(int dA, int dB) {
	if(dA==dB) return 0;

	int trkA0 = d0_idx_trk0->at(dA);
	int trkA1 = d0_idx_trk1->at(dA);
	int trkB0 = d0_idx_trk0->at(dB);
	int trkB1 = d0_idx_trk1->at(dB);

	int genA0 = trk_idx_gen->at(trkA0);
	int genA1 = trk_idx_gen->at(trkA1);
	int genB0 = trk_idx_gen->at(trkB0);
	int genB1 = trk_idx_gen->at(trkB1);

	//if only one candidate is fully truth-matched, keep that one
	if(genA0>-1 && genA1>-1 && !(genB0>-1 && genB1>-1)) return 0;
	if(genB0>-1 && genB1>-1 && !(genA0>-1 && genA1>-1)) return 1;

	//if truth-matched and only one has correct child PID assignment, keep that one
	if(genA0>-1 && genA1>-1 && genB0>-1 && genB1>-1) {
		if(TMath::Abs(gen_pid->at(genA0))==321 && TMath::Abs(gen_pid->at(genA1))==211 && !(TMath::Abs(gen_pid->at(genB0))==321 && TMath::Abs(gen_pid->at(genB1))==211)) return 0;
		if(TMath::Abs(gen_pid->at(genB0))==321 && TMath::Abs(gen_pid->at(genB1))==211 && !(TMath::Abs(gen_pid->at(genA0))==321 && TMath::Abs(gen_pid->at(genA1))==211)) return 1;
	}

	//keep candidate with higher pT
	if(d0_pt->at(dA)>d0_pt->at(dB)) return 0;
	if(d0_pt->at(dB)>d0_pt->at(dA)) return 1;

	//fall back to keeping the one with the most kaon-like kaon
	if(trk_pnn_k->at( trkA0) < trk_pnn_k->at(trkB0)) return 1;
	else return 0;
}

void makeNewCalibTuples::fillOutput(int j, int t)
{
	//fill the output tuples
	TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));

	int ipv = jet_idx_pvr->at(j);
	TVector3 pv(pvr_x->at(ipv),pvr_y->at(ipv),pvr_z->at(ipv)); 

	int ihard = -1, imu = -1, nmu = 0, jnchr = 0, jnneu = 0, ndispl6 = 0, ndispl9 = 0, ndispl16 = 0;
	double jptd = 0, jetq = 0, ry = 0, rp = 0, m11 = 0, m12 = 0, m22 = 0, sumpt2 = 0;//, pnnmu_best = 0;
	TLorentzVector p4mu,p4hard;
	for(unsigned int i=0; i<trk_p->size(); i++){
		if(trk_idx_jet->at(i) != j) continue;
		jnchr++;
		TLorentzVector p4trk(trk_px->at(i),trk_py->at(i),trk_pz->at(i),trk_e->at(i));
		int q = trk_q->at(i);
		jetq += q*p4trk.Pt();
		jptd += pow(p4trk.Pt(),2);
		if(p4trk.Pt() > p4hard.Pt()) {p4hard = p4trk; ihard = i;}
		if(trk_is_mu->at(i) > 0 && trk_pnn_mu->at(i) > 0.5 && p4trk.Pt() > 500){
			nmu++;
			//if(trk_pnnmu->at(i) > pnnmu_best) 
			if(p4trk.Pt() > p4mu.Pt()) {p4mu = p4trk; imu=i;}// pnnmu_best = trk_pnn_mu->at(i);
		}
		double dy = p4trk.Rapidity()-p4j.Rapidity();
		double dp = p4trk.DeltaPhi(p4j);
		double r = sqrt(dy*dy+dp*dp);
		ry += r*dy*p4trk.Pt();
		rp += r*dp*p4trk.Pt();
		m11 += pow(p4trk.Pt()*dy,2);
		m22 += pow(p4trk.Pt()*dp,2);
		m12 += pow(p4trk.Pt(),2)*dy*dp;
		sumpt2 += pow(p4trk.Pt(),2);
		if(p4trk.Pt() > 500){
			if(trk_ip_chi2->at(i) > 6) ndispl6++;
			if(trk_ip_chi2->at(i) > 9) ndispl9++;
			if(trk_ip_chi2->at(i) > 16) ndispl16++;
		}
	}

	jetq /= p4j.Pt();
	jptd = sqrt(jptd) / p4j.Pt();      
	for(unsigned int i=0; i<neu_p->size(); i++){
		if(neu_idx_jet->at(i) != j) continue; 
		jnneu++;
		TLorentzVector p4neu(neu_px->at(i),neu_py->at(i),neu_pz->at(i),neu_e->at(i));
		double dy = p4neu.Rapidity()-p4j.Rapidity();
		double dp = p4neu.DeltaPhi(p4j);
		double r = sqrt(dy*dy+dp*dp);
		ry += r*dy*p4neu.Pt();
		rp += r*dp*p4neu.Pt();
		m11 += pow(p4neu.Pt()*dy,2);
		m22 += pow(p4neu.Pt()*dp,2);
		m12 += pow(p4neu.Pt(),2)*dy*dp;
		sumpt2 += pow(p4neu.Pt(),2);
	}

	ry /= p4j.Pt();
	rp /= p4j.Pt();
	TMatrixD m(2,2);
	m(0,0)=m11;
	m(1,1)=m22;
	m(1,0)=m12;
	m(0,1)=m12;
	TVectorD eig;
	m.EigenVectors(eig);
	double js1=sqrt(eig[0]/sumpt2);
	double js2=sqrt(eig[1]/sumpt2);

	TLorentzVector p4mcj;
	unsigned int bestTrueJet(-1);
	// match to true jet
	for(unsigned int g=0; g<gen_pid->size(); ++g) {
		if(fabs(gen_pid->at(g)) != 98) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(p4g.DeltaR(p4j) < 0.5 && p4g.Pt() > p4mcj.Pt()) {
			p4mcj = p4g;
			bestTrueJet=g;
		}
	}

	JPX = p4j.Px();
	JPY = p4j.Py();
	JPZ = p4j.Pz();
	JE = p4j.E();
	JPT = p4j.Pt();
	JETA = p4j.Eta();
	if(p4mcj.Pt()>0.) {
		JTRUEPX = p4mcj.Px();
		JTRUEPY = p4mcj.Py();
		JTRUEPZ = p4mcj.Pz();
		JTRUEE = p4mcj.E();
		JTRUEPT = p4mcj.Pt();
		JTRUEETA = p4mcj.Eta();
		JTRUEDR = p4j.DeltaR(p4mcj);
		usedTruthJets_.insert(bestTrueJet);
	} else {
		JTRUEPX = 0.;
		JTRUEPY = 0.;
		JTRUEPZ = 0.;
		JTRUEE = 0.;
		JTRUEPT = 0.;
		JTRUEETA = 0.;
		JTRUEDR = 0.;
	}
	JS1 = js1; 
	JS2 = js2; 
	JQ = jetq; 
	JN = jnchr + jnneu; 
	JNQ = jnchr;
	JNN = jnneu;
	JPTD = jptd;
	JTRUEb=false;
	JTRUEc=false;
	JTRUED0=false;
	JTRUEDP=false;
	JTRUEDS=false;
	JTRUEJPSI=false;
	JTRUEDSV=false;
	JTRUEBSV=false;
	JTRUESV=false;
	JTRUEDPX = 0.;
	JTRUEDPY = 0.;
	JTRUEDPZ = 0.;
	JTRUEDE  = 0.;
	JTRUEDDR = 10.;
	int bestGen(-1);
	double bestGenDr(10.);
	for(unsigned int g=0; g<gen_pid->size(); ++g) {
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(p4g.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
			//std::cout << g << "\t" << gen_pid->at(g) << "\t" << longLivedC(gen_pid->at(g)) << std::endl;
			if(TMath::Abs(gen_pid->at(g))==5) JTRUEb=true;
			if(TMath::Abs(gen_pid->at(g))==4) JTRUEc=true;
			if(checkD0Truth(g)) JTRUED0=true;
			if(checkDpTruth(g)) JTRUEDP=true;
			if(checkDsTruth(g)) JTRUEDS=true;
			if(checkJpsiTruth(g)) JTRUEJPSI=true;
			if(checkDSVTruth(g)) JTRUEDSV=true;
			if(checkBSVTruth(g)) JTRUEBSV=true;
			if(checkSVTruth(g)) JTRUESV=true;
			if(longLivedC(gen_pid->at(g)) && p4j.DeltaR(p4g)<bestGenDr) {
				bestGenDr = p4j.DeltaR(p4g);
				bestGen = g;
			}
		}
	}

	if(bestGen>-1) {
		JTRUEDPX = gen_px->at(bestGen);
		JTRUEDPY = gen_py->at(bestGen);
		JTRUEDPZ = gen_pz->at(bestGen);
		JTRUEDE  = gen_e->at(bestGen);
		JTRUEDDR = bestGenDr;
	}

	if(t>-1) {
		TLorentzVector p4t(jet_px->at(t),jet_py->at(t),jet_pz->at(t),jet_e->at(t));
		DR = p4j.DeltaR(p4t);
		DPHI = p4j.DeltaPhi(p4t);

		TPX = p4t.Px();
		TPY = p4t.Py();
		TPZ = p4t.Pz();
		TE = p4t.E();
		TPT = p4t.Pt();
		TETA = p4t.Eta();
	}

	PVX = pv.X();
	PVY = pv.Y();
	PVZ = pv.Z();
	NDISPL6 = ndispl6;
	NDISPL9 = ndispl9;
	NDISPL16 = ndispl16;

	NMU = nmu;
	if(NMU > 0){	
		MUPT = p4mu.Pt();
		MUIPCHI2 = trk_ip_chi2->at(imu);
		MUDR = p4j.DeltaR(p4mu);
		MUPNN = trk_pnn_mu->at(imu);
	}
	else{
		MUPT = -1000;
		MUIPCHI2 = -1000;
		MUDR = -1000;
		MUPNN = -1000;
	}
	if(p4hard.Pt() > 500){
		HPT = p4hard.Pt();
		HIPCHI2 = trk_ip_chi2->at(ihard);
		HDR = p4hard.DeltaR(p4j);
	}
	else{
		HPT = -1000;
		HIPCHI2 = -1000;
		HDR = -1000;
	}


	fillSVCands(j,t);
	tout->Fill();
}

void makeNewCalibTuples::fillTruthOutput() {
	//std::cout << std::endl;//TODO
	for(unsigned int j=0; j<gen_pid->size(); ++j) {
		if(gen_pid->at(j)!=98) continue;
		//std::cout << std::endl;//TODO
		if(usedTruthJets_.count(j)!=0) continue;

		clearOutputs();

		TLorentzVector p4(gen_px->at(j),gen_py->at(j),gen_pz->at(j),gen_e->at(j));
		if(p4.Pt()<minTruePt_ || p4.Pt()>maxTruePt_) continue;
		////TODO start
		//for(unsigned int j2=0; j2<gen_pid->size(); ++j2) {
		//	if(gen_pid->at(j2)!=98) continue;
		//	TLorentzVector p42(gen_px->at(j2),gen_py->at(j2),gen_pz->at(j2),gen_e->at(j2));
		//	printf("\t% .4f", p4.DeltaPhi(p42));
		//}
		////TODO end
		if(gen_pid->at(j)!=98) continue;
		JTRUEPX  = p4.Px();
		JTRUEPY  = p4.Py();
		JTRUEPZ  = p4.Pz();
		JTRUEE   = p4.E();
		JTRUEPT  = p4.Pt();
		JTRUEETA = p4.Eta();
		JTRUEDR = 10.;
		JTRUEb  = 0.;
		JTRUEc  = 0.;
		JTRUED0  = 0.;
		JTRUEDP  = 0.;
		JTRUEDS  = 0.;
		JTRUEJPSI = 0.;
		JTRUEDSV = 0.;
		JTRUEBSV = 0.;
		JTRUESV = 0.;
		JTRUEDPX = 0.;
		JTRUEDPY = 0.;
		JTRUEDPZ = 0.;
		JTRUEDE  = 0.;
		JTRUEDDR = 10.;
		int bestGen(-1);
		double bestGenDr(10.);

		for(unsigned int g=0; g<gen_pid->size(); ++g) {
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
			if(p4g.Pt()>0.&&p4.Pt()>0.&&p4.DeltaR(p4g)<0.5) {
				if(TMath::Abs(gen_pid->at(g))==5) JTRUEb=true;
				if(TMath::Abs(gen_pid->at(g))==4) JTRUEc=true;
				if(checkD0Truth(g)) JTRUED0=true;
				if(checkDpTruth(g)) JTRUEDP=true;
				if(checkDsTruth(g)) JTRUEDS=true;
				if(checkJpsiTruth(g)) JTRUEJPSI=true;
				if(checkDSVTruth(g)) JTRUEDSV=true;
				if(checkBSVTruth(g)) JTRUEBSV=true;
				if(checkSVTruth(g)) JTRUESV=true;
				if(longLivedC(gen_pid->at(g)) && p4.DeltaR(p4g)<bestGenDr) {
					bestGenDr = p4.DeltaR(p4g);
					bestGen = g;
				}
			}
		}

		for(unsigned int r=0; r<jet_px->size(); ++r) {
			TLorentzVector p4r(jet_px->at(r),jet_py->at(r),jet_pz->at(r),jet_e->at(r));
			if(p4.DeltaR(p4r)<JTRUEDR) JTRUEDR=p4.DeltaR(p4r);
			//if(p4.DeltaR(p4r)<0.5) std::cout << std::endl << EVT << "\t" << r << "\t" << p4r.Pt() << "\t" << p4.Pt() << "\t" << JTRUEc << "\t" << JTRUEb;//TODO
		}

		if(bestGen>-1) {
			JTRUEDPX = gen_px->at(bestGen);
			JTRUEDPY = gen_py->at(bestGen);
			JTRUEDPZ = gen_pz->at(bestGen);
			JTRUEDE  = gen_e->at(bestGen);
			JTRUEDDR = bestGenDr;
		}


		tout->Fill();
	}
}

void makeNewCalibTuples::fillSVCands(int j, int t)
{
	NSV = 0; 
	NTSV = 0; 
	clearOutputVectors();

	TLorentzVector p4j, p4t;
	if(j>-1) p4j.SetPxPyPzE(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	if(t>-1) p4t.SetPxPyPzE(jet_px->at(t),jet_py->at(t),jet_pz->at(t),jet_e->at(t));

	int ipv(0);// currently veto any events with more than one PV - if we unveto then at least this won't break for jet data/MC
	if(j>-1) ipv = jet_idx_pvr->at(j);
	TVector3 pv(pvr_x->at(ipv),pvr_y->at(ipv),pvr_z->at(ipv)); 

	//keep generated D's
	std::map<int,int> foundDs;
	for(uint g=0; g<gen_pid->size(); ++g) {
		//if(EVT==6202) std::cout << g << "\t" << gen_pid->at(g) << "\t" << gen_idx_prnt->at(g) << std::endl;//TODO
		if(!longLivedC(gen_pid->at(g))) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(j>-1 && p4j.DeltaR(p4g)>=0.5) continue;

		std::vector<int> idcs;

		int pid = gen_pid->at(g);
		//find indices
		switch(TMath::Abs(pid)) {
			case 411:
				if(checkTruth(g, idcs,411,2,1)) break;
				//continue;//now keep ID, PX, PY, PZ, E, FROMB for all D's, not just our decay
			case 421:
				//for(int gg=g; gg<gen_pid->size(); ++gg) {//TODO
				//	std::cout << gg << "\t" << gen_pid->at(gg) << "\t" << gen_idx_prnt->at(gg) << std::endl;//TODO
				//}//TODO
				//std::cout << std::endl;//TODO
				if(checkTruth(g, idcs,421,1,1)) break;
				idcs.clear();
				if(checkTruth(g, idcs,421,3,1)) break;
				idcs.clear();
				break;
			case 431:
				if(checkTruth(g, idcs,431,1,2)) break;
				idcs.clear();
				break;
			case 443:
				if(checkTruth(g, idcs,443,0,0,0,2)) break;
				idcs.clear();
				break;
			case 4122:
				if(checkTruth(g, idcs,4122,1,1,1)) break;
				idcs.clear();
				break;
			default:
				break;
		}


		foundDs.insert(std::pair<int,int>(g, TRUEDID->size()));

		TRUEDID     ->push_back(pid);
		TRUEDPX     ->push_back(gen_px->at(g));
		TRUEDPY     ->push_back(gen_py->at(g));
		TRUEDPZ     ->push_back(gen_pz->at(g));
		TRUEDE      ->push_back(gen_e->at(g));
		TRUEDFROMB  ->push_back(checkFromB(g));
		TRUEDSEL    ->push_back(-1);
		//if(TMath::Abs(pid)==421) {
		//	std::cout << pid << " ->";
		//	for(unsigned int i=0; i<idcs.size(); ++i) {
		//		std::cout << gen_pid->at(idcs[i]) << " ";
		//	}
		//}
		//std::cout << std::endl;
		TVector3 p3;
		//double theta;
		bool inacc, reco;
		int rectrk;
		double recpnnk, recpnnpi;
		if(idcs.size()>1) {
			TRUEDTRK0IDX ->push_back(idcs[0]);
			TRUEDTRK0ID  ->push_back(gen_pid->at(idcs[0]));
			TRUEDTRK0P   ->push_back(gen_p->at(idcs[0]));
			TRUEDTRK0PT  ->push_back(gen_pt->at(idcs[0]));
			TRUEDTRK0PX  ->push_back(gen_px->at(idcs[0]));
			TRUEDTRK0PY  ->push_back(gen_py->at(idcs[0]));
			TRUEDTRK0PZ  ->push_back(gen_pz->at(idcs[0]));
			//old acceptance code
			//theta = TVector3(gen_px->at(idcs[0]),gen_py->at(idcs[0]),gen_pz->at(idcs[0])).Theta();
			//if(theta>0.01 && theta<0.4) inacc=true;
			//else inacc=false;
			//new acceptance code
			inacc = false;
			p3.SetXYZ(gen_px->at(idcs[0]),gen_py->at(idcs[0]),gen_pz->at(idcs[0]));
			if(p3.Eta()>2. && p3.Eta()<4.5 && p3.Pt()>500. && p3.Mag()>5000.) inacc=true;
			TRUEDTRK0INACC  ->push_back(inacc);
			reco=false;
			rectrk=-1;
			recpnnk=-1.;
			recpnnpi=-1.;
			for(unsigned int t=0; t<trk_idx_gen->size(); ++t) {
				if(trk_idx_gen->at(t)==idcs[0]) {
					reco=true;
					rectrk=t;
					recpnnk  = trk_pnn_k->at(rectrk);
					recpnnpi = trk_pnn_pi->at(rectrk);
					break;
				}
			}
			TRUEDTRK0RECO   ->push_back(reco);
			TRUEDTRK0PNNK   ->push_back(recpnnk);
			TRUEDTRK0PNNPI  ->push_back(recpnnpi);
			TRUEDTRK1IDX ->push_back(idcs[1]);
			TRUEDTRK1ID  ->push_back(gen_pid->at(idcs[1]));
			TRUEDTRK1P   ->push_back(gen_p->at(idcs[1]));
			TRUEDTRK1PT  ->push_back(gen_pt->at(idcs[1]));
			TRUEDTRK1PX  ->push_back(gen_px->at(idcs[1]));
			TRUEDTRK1PY  ->push_back(gen_py->at(idcs[1]));
			TRUEDTRK1PZ  ->push_back(gen_pz->at(idcs[1]));
			//old acceptance code
			//theta = TVector3(gen_px->at(idcs[1]),gen_py->at(idcs[1]),gen_pz->at(idcs[1])).Theta();
			//if(theta>0.01 && theta<0.4) inacc=true;
			//else inacc=false;
			//new acceptance code
			inacc = false;
			p3.SetXYZ(gen_px->at(idcs[1]),gen_py->at(idcs[1]),gen_pz->at(idcs[1]));
			if(p3.Eta()>2. && p3.Eta()<4.5 && p3.Pt()>500. && p3.Mag()>5000.) inacc=true;
			TRUEDTRK1INACC  ->push_back(inacc);
			reco=false;
			rectrk=-1;
			recpnnk=-1.;
			recpnnpi=-1.;
			for(unsigned int t=0; t<trk_idx_gen->size(); ++t) {
				if(trk_idx_gen->at(t)==idcs[1]) {
					reco=true;
					rectrk=t;
					recpnnk  = trk_pnn_k->at(rectrk);
					recpnnpi = trk_pnn_pi->at(rectrk);
					break;
				}
			}
			TRUEDTRK1RECO   ->push_back(reco);
			TRUEDTRK1PNNK   ->push_back(recpnnk);
			TRUEDTRK1PNNPI  ->push_back(recpnnpi);
		} else {
			TRUEDTRK0IDX ->push_back(-1);
			TRUEDTRK0ID  ->push_back(-1);
			TRUEDTRK0P   ->push_back(-1);
			TRUEDTRK0PT  ->push_back(-1);
			TRUEDTRK0PX  ->push_back(-1);
			TRUEDTRK0PY  ->push_back(-1);
			TRUEDTRK0PZ  ->push_back(-1);
			TRUEDTRK0INACC  ->push_back(-1);
			TRUEDTRK0RECO   ->push_back(-1);
			TRUEDTRK0PNNK   ->push_back(-1);
			TRUEDTRK0PNNPI  ->push_back(-1);
			TRUEDTRK1IDX ->push_back(-1);
			TRUEDTRK1ID  ->push_back(-1);
			TRUEDTRK1P   ->push_back(-1);
			TRUEDTRK1PT  ->push_back(-1);
			TRUEDTRK1PX  ->push_back(-1);
			TRUEDTRK1PY  ->push_back(-1);
			TRUEDTRK1PZ  ->push_back(-1);
			TRUEDTRK1INACC  ->push_back(-1);
			TRUEDTRK1RECO   ->push_back(-1);
			TRUEDTRK1PNNK   ->push_back(-1);
			TRUEDTRK1PNNPI  ->push_back(-1);
		}
		if(idcs.size()>2) {
			TRUEDTRK2IDX ->push_back(idcs[2]);
			TRUEDTRK2ID  ->push_back(gen_pid->at(idcs[2]));
			TRUEDTRK2P   ->push_back(gen_p->at(idcs[2]));
			TRUEDTRK2PT  ->push_back(gen_pt->at(idcs[2]));
			TRUEDTRK2PX  ->push_back(gen_px->at(idcs[2]));
			TRUEDTRK2PY  ->push_back(gen_py->at(idcs[2]));
			TRUEDTRK2PZ  ->push_back(gen_pz->at(idcs[2]));
			//old acceptance code
			//theta = TVector3(gen_px->at(idcs[2]),gen_py->at(idcs[2]),gen_pz->at(idcs[2])).Theta();
			//if(theta>0.01 && theta<0.4) inacc=true;
			//else inacc=false;
			//new acceptance code
			inacc = false;
			p3.SetXYZ(gen_px->at(idcs[2]),gen_py->at(idcs[2]),gen_pz->at(idcs[2]));
			if(p3.Eta()>2. && p3.Eta()<4.5 && p3.Pt()>500. && p3.Mag()>5000.) inacc=true;
			TRUEDTRK2INACC  ->push_back(inacc);
			reco=false;
			rectrk=-1;
			recpnnk=-1.;
			recpnnpi=-1.;
			for(unsigned int t=0; t<trk_idx_gen->size(); ++t) {
				if(trk_idx_gen->at(t)==idcs[2]) {
					reco=true;
					rectrk=t;
					recpnnk  = trk_pnn_k->at(rectrk);
					recpnnpi = trk_pnn_pi->at(rectrk);
					break;
				}
			}
			TRUEDTRK2RECO   ->push_back(reco);
			TRUEDTRK2PNNK   ->push_back(recpnnk);
			TRUEDTRK2PNNPI  ->push_back(recpnnpi);
		} else {
			TRUEDTRK2IDX ->push_back(-1);
			TRUEDTRK2ID  ->push_back(-1);
			TRUEDTRK2P   ->push_back(-1);
			TRUEDTRK2PT  ->push_back(-1);
			TRUEDTRK2PX  ->push_back(-1);
			TRUEDTRK2PY  ->push_back(-1);
			TRUEDTRK2PZ  ->push_back(-1);
			TRUEDTRK2INACC  ->push_back(-1);
			TRUEDTRK2RECO   ->push_back(-1);
			TRUEDTRK2PNNK   ->push_back(-1);
			TRUEDTRK2PNNPI  ->push_back(-1);
		}
		if(idcs.size()>3) {
			TRUEDTRK3IDX ->push_back(idcs[3]);
			TRUEDTRK3ID  ->push_back(gen_pid->at(idcs[3]));
			TRUEDTRK3P   ->push_back(gen_p->at(idcs[3]));
			TRUEDTRK3PT  ->push_back(gen_pt->at(idcs[3]));
			TRUEDTRK3PX  ->push_back(gen_px->at(idcs[3]));
			TRUEDTRK3PY  ->push_back(gen_py->at(idcs[3]));
			TRUEDTRK3PZ  ->push_back(gen_pz->at(idcs[3]));
			//old acceptance code
			//theta = TVector3(gen_px->at(idcs[3]),gen_py->at(idcs[3]),gen_pz->at(idcs[3])).Theta();
			//if(theta>0.01 && theta<0.4) inacc=true;
			//else inacc=false;
			//new acceptance code
			inacc = false;
			p3.SetXYZ(gen_px->at(idcs[3]),gen_py->at(idcs[3]),gen_pz->at(idcs[3]));
			if(p3.Eta()>2. && p3.Eta()<4.5 && p3.Pt()>500. && p3.Mag()>5000.) inacc=true;
			TRUEDTRK3INACC  ->push_back(inacc);
			reco=false;
			rectrk=-1;
			recpnnk=-1.;
			recpnnpi=-1.;
			for(unsigned int t=0; t<trk_idx_gen->size(); ++t) {
				if(trk_idx_gen->at(t)==idcs[3]) {
					reco=true;
					rectrk=t;
					recpnnk  = trk_pnn_k->at(rectrk);
					recpnnpi = trk_pnn_pi->at(rectrk);
					break;
				}
			}
			TRUEDTRK3RECO   ->push_back(reco);
			TRUEDTRK3PNNK   ->push_back(recpnnk);
			TRUEDTRK3PNNPI  ->push_back(recpnnpi);
		} else {
			TRUEDTRK3IDX ->push_back(-1);
			TRUEDTRK3ID  ->push_back(-1);
			TRUEDTRK3P   ->push_back(-1);
			TRUEDTRK3PT  ->push_back(-1);
			TRUEDTRK3PX  ->push_back(-1);
			TRUEDTRK3PY  ->push_back(-1);
			TRUEDTRK3PZ  ->push_back(-1);
			TRUEDTRK3INACC  ->push_back(-1);
			TRUEDTRK3RECO   ->push_back(-1);
			TRUEDTRK3PNNK   ->push_back(-1);
			TRUEDTRK3PNNPI  ->push_back(-1);
		}
	}

	std::set<long> foundSVs;

	for(unsigned int s = 0; s < svr_p->size(); s++){
		if(svr_z->at(s) != svr_z->at(s)) continue;//remove NaN entries
		TVector3 sv(svr_x->at(s),svr_y->at(s),svr_z->at(s));
		TVector3 fly = sv-pv;
		if(backwards_) fly = -fly;

		//Now we want to keep tag SVs too (for further tagging)
		//First check if we pass the deltaR requirement for either jet
		//Then check which is closer once we know if we have a good SV
		//If j/t==-1 then we're not using jets so we keep everything
		if(j>-1 && p4j.Vect().DeltaR(fly) > 0.5 && t>-1 && p4t.Vect().DeltaR(fly) > 0.5) continue;//reject if too far from jet
//OLD		if(j>-1 && p4j.Vect().DeltaR(fly) > 0.5) continue;//reject if too far from jet
//OLD		if(j>-1 && t>-1 && p4j.Vect().DeltaR(fly) >= p4t.Vect().DeltaR(fly)) continue;//reject if closer to tag jet
//OLD//fixed 180825	if(p4j.Vect().DeltaR(fly) > 0.5 && p4j.Vect().DeltaR(fly) >= p4t.Vect().DeltaR(fly)) continue;//reject if closer to tag jet

		double svn       = 0; 
		double svnj      = 0;
		double svnt      = 0;
		double svq       = 0; 
		double svsumipchi2 = 0;
		double svminipchi2 = 1e9;
		double svmaxghost = 0; 
		bool svisd0(false), svisdp(false);
		double svd0m(0.), svdpm(0.);
		int svtrued(-1);
		std::vector<int> svtrueindices;
		std::vector<int> svrecoindices;
		std::vector<double> svrecop;
		std::vector<double> svrecopt;
		std::vector<double> svrecopx;
		std::vector<double> svrecopy;
		std::vector<double> svrecopz;
		std::vector<double> svrecopnnpi;
		std::vector<double> svrecopnnk;
		//	TLorentzVector p4sv(svpx->at(s),svpy->at(s),svpz->at(s),sve->at(s));
		TLorentzVector p4sv;
		p4sv.SetPxPyPzE(0.,0.,0.,0.);
		bool veto=false;
		for(int i=0; i<10; i++){
			if(svtrk[i]->at(s) < 0) break;
			svn++;
			int ii = svtrk[i]->at(s);
//			std::cout << i << "\t" << ii << "\t" << trk_idx_gen->at(ii) << "\t" << trk_pid->at(ii) << "\t" << ((trk_idx_gen->at(ii)>=0)?gen_pid->at(trk_idx_gen->at(ii)):0) << std::endl;//TODO
			if(trk_idx_gen->at(ii) >= 0) svtrueindices.push_back(trk_idx_gen->at(ii));
			svrecoindices.push_back(ii);
			svrecop.push_back(trk_p->at(ii));
			svrecopt.push_back(trk_pt->at(ii));
			svrecopx.push_back(trk_px->at(ii));
			svrecopy.push_back(trk_py->at(ii));
			svrecopz.push_back(trk_pz->at(ii));
			svrecopnnpi.push_back(trk_pnn_pi->at(ii));
			svrecopnnk.push_back(trk_pnn_k->at(ii));
			TVector3 hit = TVector3(trk_x->at(ii),trk_y->at(ii),trk_z->at(ii));
			if(hit.Z() < sv.Z()) veto=true;
			TLorentzVector p4trk(trk_px->at(ii),trk_py->at(ii),trk_pz->at(ii),trk_e->at(ii));
			p4sv += p4trk;
			if(j>-1 && p4trk.DeltaR(p4j) < 0.5) svnj++;
			if(t>-1 && p4trk.DeltaR(p4t) < 0.5) svnt++;
			svq += trk_q->at(ii);
			double ipchi2 = trk_ip_chi2->at(ii);
			//cout << "sv track IP chi2: " << ipchi2 << endl;
			svsumipchi2 += ipchi2;
			if(ipchi2 < svminipchi2) svminipchi2 = ipchi2;
			double ghost = trk_prb_ghost->at(ii);
			if(ghost > svmaxghost) svmaxghost = ghost;
		}
		if(svmaxghost>0.2) continue;
		if(fly.Z()*p4sv.M()/p4sv.Pz()/(3e11)*(1e12) > 10) continue;	  
		if(svr_ip_chi2_min_trk->at(s)<minipchi2cut_ || svr_m_cor_err_full->at(s)>maxmcorerrcut_) continue;
		if(p4sv.M()<400.) continue;
		if (!(svr_fd_min->at(s) < 15)) continue;
		if (!(svr_z->at(s) < 200)) continue;
		if (!(svr_m_cor->at(s) > 600)) continue;
		if (!(p4sv.Pt() > 2000)) continue;
		if (!(svr_fd_chi2->at(s) > 32)) continue;
		if(detector->distance(svr_x->at(s)+detOffsetX,svr_y->at(s)+detOffsetY,svr_z->at(s),svr_dx->at(s),svr_dy->at(s),svr_dz->at(s))<0.5) continue;
		if(/*svnj < 1 || */veto) continue;

		if(svn>4) continue;//throw away if more than 4 tracks
		//create a unique fingerprint for each set of reco tracks (ignore permutations)
		std::sort(svrecoindices.begin(), svrecoindices.end());
		int fingerprint = 0;
		for(unsigned int itrk=0; itrk<svrecoindices.size(); ++itrk) {
			fingerprint += svrecoindices[itrk]*TMath::Power(trk_pt->size(),static_cast<int>(itrk));
		}
		if(!foundSVs.insert(fingerprint).second) continue;//if insert into set fails then we already have this SV

//		for(int i=0; i<svtrueindices.size(); ++i) {//TODO
//			std::cout << "---" << i << "\t" << svtrueindices[i] << std::endl;//TODO
//		}//TODO

		//try to make D0 candidates from the SVs
		if(svn==2) {
			int trk0 = svtrk[0]->at(s);
			int trk1 = svtrk[1]->at(s);
			TLorentzVector p0;
			TLorentzVector p1;

			int dg(-1);
			if(checkD0Truth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),dg,false)) {
				if(foundDs.find(dg)!=foundDs.end()) svtrued = foundDs.at(dg);
			}

			if(isMC_) {//for MC don't apply PID, use highest PIDk as kaon
				//std::cout << (trk_pnn_k->at(trk0)>trk_pnn_k->at(trk1)) << "\t" << trk_pnn_k->at(trk0) << "\t" << trk_pnn_k->at(trk1) << std::endl;
				if(trk_pnn_k->at(trk0)>trk_pnn_k->at(trk1)) {
					//std::cout << "A" << std::endl;
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					svisd0 = true;
					svd0m = (p0+p1).M();
				} else {
					//std::cout << "B" << std::endl;
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					svisd0 = true;
					svd0m = (p0+p1).M();
				}
				//TLorentzVector pA; pA.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
				//TLorentzVector pB; pB.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
				//TLorentzVector pC; pC.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
				//TLorentzVector pD; pD.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
				//std::cout << svd0m << "\t" << (pA+pB).M() << "\t" << (pC+pD).M() << std::endl;
			} else {
				if(trk_pnn_pi->at(trk1)>0.2 && trk_pnn_k->at( trk0)>0.3) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					svisd0 = true;
					svd0m = (p0+p1).M();
					//SV2D0M->push_back((p0+p1).M());
				} else if(trk_pnn_pi->at(trk0)>0.2 && trk_pnn_k->at( trk1)>0.3) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					svisd0 = true;
					svd0m = (p0+p1).M();
					//SV2D0M->push_back((p0+p1).M());
				}
			}
		}
		//std::cout << svd0m << std::endl;//TODO

		//try to make D+ candidates from the SVs
		if(svn==3) {
			int trk0 = svtrk[0]->at(s);
			int trk1 = svtrk[1]->at(s);
			int trk2 = svtrk[2]->at(s);
			TLorentzVector p0;
			TLorentzVector p1;
			TLorentzVector p2;

			int dg(-1);
			if(checkDpTruth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),trk_idx_gen->at(trk2),dg,false)) {
				if(foundDs.find(dg)!=foundDs.end()) svtrued = foundDs.at(dg);
			}

			if(isMC_) {//for MC don't apply PID, use highest PIDk as kaon
				if(trk_pnn_k->at(trk0)>trk_pnn_k->at(trk1) && trk_pnn_k->at(trk0)>trk_pnn_k->at(trk2)) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),139.6);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				} else if(trk_pnn_k->at(trk1)>trk_pnn_k->at(trk2)) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),139.6);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				} else {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),493.7);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				}
			} else {
				if(trk_pnn_pi->at(trk2)>0.2 && trk_pnn_pi->at(trk1)>0.2 && trk_pnn_k->at( trk0)>0.3) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),139.6);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				} else if(trk_pnn_pi->at(trk2)>0.2 && trk_pnn_k->at( trk1)>0.3 && trk_pnn_pi->at(trk0)>0.2) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),139.6);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				} else if(trk_pnn_k->at(trk2)>0.3 && trk_pnn_pi->at( trk1)>0.2 && trk_pnn_pi->at(trk0)>0.2) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					p2.SetXYZM(trk_px->at(trk2),trk_py->at(trk2),trk_pz->at(trk2),493.7);
					svisdp = true;
					svdpm = (p0+p1+p2).M();
				}
			}
		}

		//now check where to store this SV
		bool inTag(false);
		//if j not set then we're not using jets, if t not set then there is no tag - either way, put it in SV...
		if(j==-1 || t==-1) inTag=false;
		//else if j set, dRjet<0.5 and dRjet<dRtag - put it in SV...
		else if(p4j.Vect().DeltaR(fly) < 0.5 && p4j.Vect().DeltaR(fly) < p4t.Vect().DeltaR(fly)) inTag=false;
		//else if dRtag<0.5 - put it in TSV...
		else if(p4t.Vect().DeltaR(fly) < 0.5) inTag=true;

		if(inTag) {
			NTSV++;

			TSVX          ->push_back(sv.X());
			TSVY          ->push_back(sv.Y());
			TSVZ          ->push_back(sv.Z());
			TSVPERP       ->push_back(sv.Perp());
			TSVMCOR       ->push_back(svr_m_cor->at(s)); 
			TSVMCORERR    ->push_back(svr_m_cor_err_full->at(s)); 
			TSVMINPERP    ->push_back(svr_fd_min->at(s));
			TSVDRJ        ->push_back(p4j.Vect().DeltaR(fly));
			TSVDRT        ->push_back(p4t.Vect().DeltaR(fly));
			TSVN          ->push_back(svn);
			TSVNTRUE      ->push_back(checkRealSV(svtrueindices));
			TSVNJ         ->push_back(svnj);
			TSVNT         ->push_back(svnt);
			TSVQ          ->push_back(svq);
			TSVSUMIPCHI2  ->push_back(svsumipchi2);
			TSVIPCHI2     ->push_back(svr_ip_chi2->at(s));
			TSVMINIPCHI2  ->push_back(svminipchi2);
			TSVMAXGHOST   ->push_back(svmaxghost);
			TSVPX         ->push_back(p4sv.Px());
			TSVPY         ->push_back(p4sv.Py());
			TSVPZ         ->push_back(p4sv.Pz());
			TSVE          ->push_back(p4sv.E());
			TSVPT         ->push_back(p4sv.Pt());
			TSVETA        ->push_back(p4sv.Eta());
			TSVM          ->push_back(p4sv.M());
			TSVTZ         ->push_back(fly.Z()*p4sv.M()/p4sv.Pz()/(3e11)*(1e12));
			TSVISD0       ->push_back(svisd0);
			TSVISDP       ->push_back(svisdp);
			TSVD0M        ->push_back(svd0m);
			TSVDPM        ->push_back(svdpm);
			TSVTRUEIDX    ->push_back(svtrued);
		} else {
			NSV++;

			SVX          ->push_back(sv.X());
			SVY          ->push_back(sv.Y());
			SVZ          ->push_back(sv.Z());
			SVPERP       ->push_back(sv.Perp());
			SVMCOR       ->push_back(svr_m_cor->at(s)); 
			SVMCORERR    ->push_back(svr_m_cor_err_full->at(s)); 
			SVMINPERP    ->push_back(svr_fd_min->at(s));
			if(j>-1)SVDRJ->push_back(p4j.Vect().DeltaR(fly));
			if(t>-1)SVDRT->push_back(p4t.Vect().DeltaR(fly));
			SVN          ->push_back(svn);
			SVNTRUE      ->push_back(checkRealSV(svtrueindices));
			SVNJ         ->push_back(svnj);
			SVNT         ->push_back(svnt);
			SVQ          ->push_back(svq);
			SVSUMIPCHI2  ->push_back(svsumipchi2);
			SVIPCHI2     ->push_back(svr_ip_chi2->at(s));
			SVMINIPCHI2  ->push_back(svminipchi2);
			SVMAXGHOST   ->push_back(svmaxghost);
			SVPX         ->push_back(p4sv.Px());
			SVPY         ->push_back(p4sv.Py());
			SVPZ         ->push_back(p4sv.Pz());
			SVE          ->push_back(p4sv.E());
			SVPT         ->push_back(p4sv.Pt());
			SVETA        ->push_back(p4sv.Eta());
			SVM          ->push_back(p4sv.M());
			SVTZ         ->push_back(fly.Z()*p4sv.M()/p4sv.Pz()/(3e11)*(1e12));
			SVISD0       ->push_back(svisd0);
			SVISDP       ->push_back(svisdp);
			SVD0M        ->push_back(svd0m);
			SVDPM        ->push_back(svdpm);
			SVTRUEIDX    ->push_back(svtrued);
			SVTRUETRK0IDX    ->push_back((svtrueindices.size()>0?svtrueindices[0]:-1));
			SVTRUETRK1IDX    ->push_back((svtrueindices.size()>1?svtrueindices[1]:-1));
			SVTRUETRK2IDX    ->push_back((svtrueindices.size()>2?svtrueindices[2]:-1));
			SVTRUETRK3IDX    ->push_back((svtrueindices.size()>3?svtrueindices[3]:-1));
			SVTRK0P          ->push_back((svrecop.size()>0?svrecop[0]:-1));
			SVTRK1P          ->push_back((svrecop.size()>1?svrecop[1]:-1));
			SVTRK2P          ->push_back((svrecop.size()>2?svrecop[2]:-1));
			SVTRK3P          ->push_back((svrecop.size()>3?svrecop[3]:-1));
			SVTRK0PT         ->push_back((svrecopt.size()>0?svrecopt[0]:-1));
			SVTRK1PT         ->push_back((svrecopt.size()>1?svrecopt[1]:-1));
			SVTRK2PT         ->push_back((svrecopt.size()>2?svrecopt[2]:-1));
			SVTRK3PT         ->push_back((svrecopt.size()>3?svrecopt[3]:-1));
			SVTRK0PX         ->push_back((svrecopx.size()>0?svrecopx[0]:-1));
			SVTRK1PX         ->push_back((svrecopx.size()>1?svrecopx[1]:-1));
			SVTRK2PX         ->push_back((svrecopx.size()>2?svrecopx[2]:-1));
			SVTRK3PX         ->push_back((svrecopx.size()>3?svrecopx[3]:-1));
			SVTRK0PY         ->push_back((svrecopy.size()>0?svrecopy[0]:-1));
			SVTRK1PY         ->push_back((svrecopy.size()>1?svrecopy[1]:-1));
			SVTRK2PY         ->push_back((svrecopy.size()>2?svrecopy[2]:-1));
			SVTRK3PY         ->push_back((svrecopy.size()>3?svrecopy[3]:-1));
			SVTRK0PZ         ->push_back((svrecopz.size()>0?svrecopz[0]:-1));
			SVTRK1PZ         ->push_back((svrecopz.size()>1?svrecopz[1]:-1));
			SVTRK2PZ         ->push_back((svrecopz.size()>2?svrecopz[2]:-1));
			SVTRK3PZ         ->push_back((svrecopz.size()>3?svrecopz[3]:-1));
			SVTRK0PNNPI      ->push_back((svrecopnnpi.size()>0?svrecopnnpi[0]:-1));
			SVTRK1PNNPI      ->push_back((svrecopnnpi.size()>1?svrecopnnpi[1]:-1));
			SVTRK2PNNPI      ->push_back((svrecopnnpi.size()>2?svrecopnnpi[2]:-1));
			SVTRK3PNNPI      ->push_back((svrecopnnpi.size()>3?svrecopnnpi[3]:-1));
			SVTRK0PNNK       ->push_back((svrecopnnk.size()>0?svrecopnnk[0]:-1));
			SVTRK1PNNK       ->push_back((svrecopnnk.size()>1?svrecopnnk[1]:-1));
			SVTRK2PNNK       ->push_back((svrecopnnk.size()>2?svrecopnnk[2]:-1));
			SVTRK3PNNK       ->push_back((svrecopnnk.size()>3?svrecopnnk[3]:-1));
		}
	}

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
		//Do PID after checking for truth match to keep PID efficiency separate
		if(!(pd.Pt()>2000)) continue;
		if(!(p4Dtrk0.Pt()>500)) continue;
		if(!(p4Dtrk1.Pt()>500)) continue;
		if(!(trk_prb_ghost->at(trk0)<0.2 && trk_prb_ghost->at(trk1)<0.2)) continue;

		//try to match to a true D0 and update the TRUED if we find it
		int dg(-1);
		bool dfoundtrue = checkD0Truth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),dg);
		int dtrueid(-1);
		int dfromb(-1);
		int ddr(-1);
		if(dg>-1) {
			if(foundDs.find(dg)!=foundDs.end()) {
				dtrueid = foundDs.at(dg);
				if(TRUEDSEL->at(dtrueid)<0.) TRUEDSEL->at(dtrueid) = 999.;//999 signifies found but not yet matched to a kept D0 (which may fail PID cuts)
			}
			dfromb = checkFromB(dg);
			TLorentzVector p4g(gen_px->at(dg),gen_py->at(dg),gen_pz->at(dg),gen_e->at(dg));
			if(j>-1) ddr = p4j.DeltaR(p4g);
		}

		//PID off for MC (simply require kaon to be the most kaon-like)
		if(!isMC_) {
			if(trk_pnn_k->at( trk0)<0.1) {//PID cut off for pion
				if(trk_pt->at(trk0)<25000. && trk_p->at(trk0)<500000.) continue;//PID turned off for kaons with PT>25GeV or P>500GeV
			}
		}
		//for MC only use PID to distinguish double mis-ID
		if( isMC_ && d0_m->size()>1 ) {

			bool foundBetter(false);
			for(uint jD=0; jD<d0_m->size(); ++jD) {
				if(jD==iD) continue;
				//first check if there's a second D0 with the same tracks
				if((trk0==d0_idx_trk0->at(jD) && trk1==d0_idx_trk1->at(jD)) ||
				   (trk0==d0_idx_trk1->at(jD) && trk1==d0_idx_trk0->at(jD))) {
					if(checkBestD0Cand(iD,jD)==1) {
						foundBetter=true;
						break;
					}
				}
				//also check for cases where multiple tracks have same truth match
				int gen0 = trk_idx_gen->at(trk0);
				int gen1 = trk_idx_gen->at(trk1);
				if(gen0>-1 && gen1>-1) {
					if((gen0==trk_idx_gen->at(d0_idx_trk0->at(jD)) && gen1==trk_idx_gen->at(d0_idx_trk1->at(jD))) ||
					   (gen0==trk_idx_gen->at(d0_idx_trk1->at(jD)) && gen1==trk_idx_gen->at(d0_idx_trk0->at(jD)))) {
						if(checkBestD0Cand(iD,jD)==1) {
							foundBetter=true;
							break;
						}
					}
				} else if(gen0>-1) {
					if((gen0==trk_idx_gen->at(d0_idx_trk0->at(jD)) && trk1==d0_idx_trk1->at(jD)) ||
					   (gen0==trk_idx_gen->at(d0_idx_trk1->at(jD)) && trk1==d0_idx_trk0->at(jD))) {
						if(checkBestD0Cand(iD,jD)==1) {
							foundBetter=true;
							break;
						}
					}
				} else if(gen1>-1) {
					if((gen1==trk_idx_gen->at(d0_idx_trk0->at(jD)) && trk0==d0_idx_trk1->at(jD)) ||
					   (gen1==trk_idx_gen->at(d0_idx_trk1->at(jD)) && trk0==d0_idx_trk0->at(jD))) {
						if(checkBestD0Cand(iD,jD)==1) {
							foundBetter=true;
							break;
						}
					}
				}
			}
			if(foundBetter) continue;
		}

		D0M       ->push_back(d0_m->at(iD));
		D0PX      ->push_back(d0_px->at(iD));
		D0PY      ->push_back(d0_py->at(iD));
		D0PZ      ->push_back(d0_pz->at(iD));
		D0E       ->push_back(d0_e->at(iD));
		D0X       ->push_back(d0_x->at(iD));
		D0Y       ->push_back(d0_y->at(iD));
		D0Z       ->push_back(d0_z->at(iD));
		D0IP      ->push_back(d0_ip->at(iD));
		D0IPCHI2  ->push_back(d0_ip_chi2->at(iD));
		D0FD      ->push_back(d0_fd->at(iD));
		D0FDCHI2  ->push_back(d0_fd_chi2->at(iD));
		D0TAU     ->push_back(d0_tau->at(iD));
		D0DIRA    ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		D0VTXCHI2 ->push_back(d0_vtx_chi2->at(iD));
		D0VTXNDOF ->push_back(d0_vtx_ndof->at(iD));
		D0DRJET   ->push_back(dDRj);
		D0DRTAG   ->push_back(dDRt);

		D0KPNNK   ->push_back(trk_pnn_k->at( trk0));
		D0KPNNPI  ->push_back(trk_pnn_pi->at(trk0));
		D0PIPNNK  ->push_back(trk_pnn_k->at( trk1));
		D0PIPNNPI ->push_back(trk_pnn_pi->at(trk1));

		D0P       ->push_back(pd.Mag());
		D0PT      ->push_back(pd.Pt());
		D0ETA     ->push_back(pd.Eta());
		D0KP      ->push_back(p4Dtrk0.P());
		D0PIP     ->push_back(p4Dtrk1.P());
		D0KPT     ->push_back(p4Dtrk0.Pt());
		D0PIPT    ->push_back(p4Dtrk1.Pt());
		D0KETA    ->push_back(p4Dtrk0.Eta());
		D0PIETA   ->push_back(p4Dtrk1.Eta());
		D0KPX     ->push_back(p4Dtrk0.Px());
		D0PIPX    ->push_back(p4Dtrk1.Px());
		D0KPY     ->push_back(p4Dtrk0.Py());
		D0PIPY    ->push_back(p4Dtrk1.Py());
		D0KPZ     ->push_back(p4Dtrk0.Pz());
		D0PIPZ    ->push_back(p4Dtrk1.Pz());
		D0KIPCHI2 ->push_back(trk_ip_chi2->at(trk0));
		D0PIIPCHI2->push_back(trk_ip_chi2->at(trk1));

		double pk   = p4Dtrk0.P();
		double ptk  = p4Dtrk0.Pt();
		double ppi  = p4Dtrk1.P();
		double ptpi = p4Dtrk1.Pt();

		if(pk  >500000.) pk  =499000.;
		if(ppi >500000.) ppi =499000.;
		if(ptk > 50000.) ptk = 49000.;
		if(ptpi> 25000.) ptpi= 24000.;

		D0KWEIGHT ->push_back(pidK->GetBinContent(pidK->FindBin(pk,ptk)));
		D0PIWEIGHT->push_back(1.);//pidPi->GetBinContent(pidPi->FindBin(ppi,ptpi)));

		D0TRK0    ->push_back(trk0);
		D0TRK1    ->push_back(trk1);
		D0TRUETRK0->push_back(trk_idx_gen->at(trk0));
		D0TRUETRK1->push_back(trk_idx_gen->at(trk1));

		D0NJ      ->push_back(nj);
		D0MAXDR   ->push_back(maxdr);

		D0TRUE    ->push_back(dfoundtrue);
		D0TRUEIDX->push_back(dtrueid);
		D0FROMB->push_back(dfromb);
		D0TRUEDR->push_back(ddr);

		if(dtrueid>-1) {
			if(TRUEDSEL->at(dtrueid)<0. || TRUEDSEL->at(dtrueid)==999.) {//if multiple match and pass PID then keep the first
				TRUEDSEL->at(dtrueid) = D0TRUE->size()-1;//also link the true D back to the reconstructed D0
			}
		}
	}

	for(uint iD=0; iD<dp_m->size(); ++iD) {
		TLorentzVector p4D(dp_px->at(iD),
				dp_py->at(iD),
				dp_pz->at(iD),
				dp_e->at( iD));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk;

		int trk0=dp_idx_trk0->at(iD);
		int trk1=dp_idx_trk1->at(iD);
		int trk2=dp_idx_trk2->at(iD);

		nj = (trk_idx_jet->at(trk0)==j) + 
		     (trk_idx_jet->at(trk1)==j) + 
		     (trk_idx_jet->at(trk2)==j);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk0),
				  trk_py->at(trk0),
				  trk_pz->at(trk0),
				  trk_e ->at(trk0));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk1),
				  trk_py->at(trk1),
				  trk_pz->at(trk1),
				  trk_e ->at(trk1));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk2),
				  trk_py->at(trk2),
				  trk_pz->at(trk2),
				  trk_e ->at(trk2));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);

		//std::cout << nj << "\t" << dDRj << "\t" << dr0 << "\t" << dr1 << "\t" << dr2 << "\t" << dDRt << std::endl;

		//if(nj == 0) continue;//require at least one track in probe jet
//TODO		if((trk_ip_chi2->at(trk0)>50)+(trk_ip_chi2->at(trk1)>50)+(trk_ip_chi2->at(trk2)>50) == 0 ||
//TODO		   (trk_ip_chi2->at(trk0)>10)+(trk_ip_chi2->at(trk1)>10)+(trk_ip_chi2->at(trk2)>10) <= 1 ||
//TODO		   dp_tau->at(iD)<0.00015 || dp_vtx_chi2->at(iD)/dp_vtx_ndof->at(iD)>6.) continue;

		TVector3 pd(dp_px->at(iD),dp_py->at(iD),dp_pz->at(iD));
		TVector3 xd(dp_x ->at(iD) - pvr_x->at(dp_idx_pvr->at(iD)),
		            dp_y ->at(iD) - pvr_y->at(dp_idx_pvr->at(iD)),
		            dp_z ->at(iD) - pvr_z->at(dp_idx_pvr->at(iD)));
		
//TODO		if(pd.Dot(xd)/pd.Mag()/xd.Mag() < 0.9999) continue;

		//PID and pT requirements 
		//PID off for MC
		if(!isMC_ && 
		   !(trk_pnn_pi->at(trk2)>0.2 &&
		     trk_pnn_pi->at(trk1)>0.2 && 
		     trk_pnn_k->at( trk0)>0.3)) continue;
		if(!(pd.Pt()>2000)) continue;

		DPM      ->push_back(dp_m->at(iD));
		DPPX     ->push_back(dp_px->at(iD));
		DPPY     ->push_back(dp_py->at(iD));
		DPPZ     ->push_back(dp_pz->at(iD));
		DPE      ->push_back(dp_e->at(iD));
		DPX      ->push_back(dp_x->at(iD));
		DPY      ->push_back(dp_y->at(iD));
		DPZ      ->push_back(dp_z->at(iD));
		DPIP     ->push_back(dp_ip->at(iD));
		DPIPCHI2 ->push_back(dp_ip_chi2->at(iD));
		DPFD     ->push_back(dp_fd->at(iD));
		DPFDCHI2 ->push_back(dp_fd_chi2->at(iD));
		DPTAU    ->push_back(dp_tau->at(iD));
		DPDIRA   ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		DPVTXCHI2->push_back(dp_vtx_chi2->at(iD));
		DPVTXNDOF->push_back(dp_vtx_ndof->at(iD));
		DPDRJET  ->push_back(dDRj);
		DPDRTAG  ->push_back(dDRt);
		DPKPNN   ->push_back(trk_pnn_k->at( trk0));
		DPPI1PNN ->push_back(trk_pnn_pi->at(trk1));
		DPPI2PNN ->push_back(trk_pnn_pi->at(trk2));

		DPP       ->push_back(dp_p->at(iD));
		DPPT      ->push_back(dp_pt->at(iD));
		DPKPT     ->push_back(trk_pt->at(trk0));
		DPPI1PT   ->push_back(trk_pt->at(trk1));
		DPPI2PT   ->push_back(trk_pt->at(trk2));

		DPTRK0   ->push_back(trk0);
		DPTRK1   ->push_back(trk1);
		DPTRK2   ->push_back(trk2);

		DPNJ      ->push_back(nj);
		DPMAXDR   ->push_back(maxdr);

		int dg(-1);
		DPTRUE    ->push_back(checkDpTruth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),trk_idx_gen->at(trk2),dg));
		int dtrueid(-1);
		int dfromb(-1);
		int ddr(-1);
		if(dg>-1) {
			if(foundDs.find(dg)!=foundDs.end()) dtrueid = foundDs.at(dg);
			dfromb = checkFromB(dg);
			TLorentzVector p4g(gen_px->at(dg),gen_py->at(dg),gen_pz->at(dg),gen_e->at(dg));
			if(j>-1) ddr = p4j.DeltaR(p4g);
		}
		DPTRUEIDX->push_back(dtrueid);
		DPFROMB->push_back(dfromb);
		DPTRUEDR->push_back(ddr);
	}

	for(uint iD=0; iD<ds_m->size(); ++iD) {
		TLorentzVector p4D(ds_px->at(iD),
				   ds_py->at(iD),
				   ds_pz->at(iD),
				   ds_e->at( iD));

		TLorentzVector phi(trk_px->at(ds_idx_trk0->at(iD))+trk_px->at(ds_idx_trk1->at(iD)),
				   trk_py->at(ds_idx_trk0->at(iD))+trk_py->at(ds_idx_trk1->at(iD)),
				   trk_pz->at(ds_idx_trk0->at(iD))+trk_pz->at(ds_idx_trk1->at(iD)),
				   trk_e ->at(ds_idx_trk0->at(iD))+trk_e ->at(ds_idx_trk1->at(iD)));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk;

		int trk0=ds_idx_trk0->at(iD);
		int trk1=ds_idx_trk1->at(iD);
		int trk2=ds_idx_trk2->at(iD);

		nj = (trk_idx_jet->at(trk0)==j) + 
		     (trk_idx_jet->at(trk1)==j) + 
		     (trk_idx_jet->at(trk2)==j);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk0),
				  trk_py->at(trk0),
				  trk_pz->at(trk0),
				  trk_e ->at(trk0));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk1),
				  trk_py->at(trk1),
				  trk_pz->at(trk1),
				  trk_e ->at(trk1));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk2),
				  trk_py->at(trk2),
				  trk_pz->at(trk2),
				  trk_e ->at(trk2));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);

		//if(nj == 0) continue;//require at least one track in probe jet
//TODO		if((trk_ip_chi2->at(trk0)>50)+(trk_ip_chi2->at(trk1)>50)+(trk_ip_chi2->at(trk2)>50) == 0 ||
//TODO		   (trk_ip_chi2->at(trk0)>10)+(trk_ip_chi2->at(trk1)>10)+(trk_ip_chi2->at(trk2)>10) <= 1 ||
//TODO		   ds_tau->at(iD)<0.00015 || ds_vtx_chi2->at(iD)/ds_vtx_ndof->at(iD)>6.) continue;

		TVector3 pd(ds_px->at(iD),ds_py->at(iD),ds_pz->at(iD));
		TVector3 xd(ds_x ->at(iD) - pvr_x->at(ds_idx_pvr->at(iD)),
		            ds_y ->at(iD) - pvr_y->at(ds_idx_pvr->at(iD)),
		            ds_z ->at(iD) - pvr_z->at(ds_idx_pvr->at(iD)));
		
//TODO		if(pd.Dot(xd)/pd.Mag()/xd.Mag() < 0.9999) continue;

		if(phi.M()<990.||phi.M()>1050.) continue;

		//PID and pT requirements 
		//PID off for MC
		if(!isMC_ &&
		   !(trk_pnn_pi->at(trk2)>0.2 &&
		     trk_pnn_k->at( trk1)>0.3 && 
		     trk_pnn_k->at( trk0)>0.3)) continue;
		if(!(pd.Pt()>2000)) continue;

		DSM      ->push_back(ds_m->at(iD));
		DSPX     ->push_back(ds_px->at(iD));
		DSPY     ->push_back(ds_py->at(iD));
		DSPZ     ->push_back(ds_pz->at(iD));
		DSE      ->push_back(ds_e->at(iD));
		DSX      ->push_back(ds_x->at(iD));
		DSY      ->push_back(ds_y->at(iD));
		DSZ      ->push_back(ds_z->at(iD));
		DSIP     ->push_back(ds_ip->at(iD));
		DSIPCHI2 ->push_back(ds_ip_chi2->at(iD));
		DSFD     ->push_back(ds_fd->at(iD));
		DSFDCHI2 ->push_back(ds_fd_chi2->at(iD));
		DSTAU    ->push_back(ds_tau->at(iD));
		DSDIRA   ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		DSVTXCHI2->push_back(ds_vtx_chi2->at(iD));
		DSVTXNDOF->push_back(ds_vtx_ndof->at(iD));
		DSDRJET  ->push_back(dDRj);
		DSDRTAG  ->push_back(dDRt);
		DSK1PNN  ->push_back(trk_pnn_k->at( trk1));
		DSK2PNN  ->push_back(trk_pnn_k->at( trk0));
		DSPIPNN  ->push_back(trk_pnn_pi->at(trk2));
		DSPHIM   ->push_back(phi.M());

		DSP       ->push_back(ds_p->at(iD));
		DSPT      ->push_back(ds_pt->at(iD));
		DSK1PT    ->push_back(trk_pt->at(trk0));
		DSK2PT    ->push_back(trk_pt->at(trk1));
		DSPIPT    ->push_back(trk_pt->at(trk2));
		DSPHIPT   ->push_back(phi.Pt());

		DSTRK0   ->push_back(trk0);
		DSTRK1   ->push_back(trk1);
		DSTRK2   ->push_back(trk2);

		DSNJ      ->push_back(nj);
		DSMAXDR   ->push_back(maxdr);

		int dg(-1);
		DSTRUE    ->push_back(checkDsTruth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),trk_idx_gen->at(trk2),dg));
		int dtrueid(-1);
		int dfromb(-1);
		int ddr(-1);
		if(dg>-1) {
			if(foundDs.find(dg)!=foundDs.end()) dtrueid = foundDs.at(dg);
			dfromb = checkFromB(dg);
			TLorentzVector p4g(gen_px->at(dg),gen_py->at(dg),gen_pz->at(dg),gen_e->at(dg));
			if(j>-1) ddr = p4j.DeltaR(p4g);
		}
		DSTRUEIDX->push_back(dtrueid);
		DSFROMB->push_back(dfromb);
		DSTRUEDR->push_back(ddr);
	}

	for(uint iD=0; iD<lc_m->size(); ++iD) {
		TLorentzVector p4D(lc_px->at(iD),
				lc_py->at(iD),
				lc_pz->at(iD),
				lc_e->at( iD));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk;

		int trk0=lc_idx_trk0->at(iD);
		int trk1=lc_idx_trk1->at(iD);
		int trk2=lc_idx_trk2->at(iD);

		nj = (trk_idx_jet->at(trk0)==j) + 
		     (trk_idx_jet->at(trk1)==j) + 
		     (trk_idx_jet->at(trk2)==j);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk0),
				  trk_py->at(trk0),
				  trk_pz->at(trk0),
				  trk_e ->at(trk0));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk1),
				  trk_py->at(trk1),
				  trk_pz->at(trk1),
				  trk_e ->at(trk1));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk2),
				  trk_py->at(trk2),
				  trk_pz->at(trk2),
				  trk_e ->at(trk2));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);

		double maxdoca(0.);
		TVector3 v;
		TVector3 p0(trk_px->at(trk0),
			    trk_py->at(trk0),
			    trk_pz->at(trk0));
		TVector3 x0(trk_x ->at(trk0),
			    trk_y ->at(trk0),
			    trk_z ->at(trk0));
		TVector3 p1(trk_px->at(trk1),
			    trk_py->at(trk1),
			    trk_pz->at(trk1));
		TVector3 x1(trk_x ->at(trk1),
			    trk_y ->at(trk1),
			    trk_z ->at(trk1));
		TVector3 p2(trk_px->at(trk2),
			    trk_py->at(trk2),
			    trk_pz->at(trk2));
		TVector3 x2(trk_x ->at(trk2),
			    trk_y ->at(trk2),
			    trk_z ->at(trk2));

		maxdoca = TMath::Max(maxdoca, calcDoca(v, p0, x0, p1, x1));
		maxdoca = TMath::Max(maxdoca, calcDoca(v, p0, x0, p2, x2));
		maxdoca = TMath::Max(maxdoca, calcDoca(v, p1, x1, p2, x2));

//		if(maxdoca<1.) std::cout << maxdoca << std::endl;

		//if(nj == 0) continue;//require at least one track in probe jet
//TODO		if((trk_ip_chi2->at(trk0)>6)+(trk_ip_chi2->at(trk1)>6)+(trk_ip_chi2->at(trk2)>6) == 0 ||
//TODO		   (trk_ip_chi2->at(trk0)>4)+(trk_ip_chi2->at(trk1)>4)+(trk_ip_chi2->at(trk2)>4) <= 1) continue;// || 
		   //lc_vtx_chi2->at(iD)/lc_vtx_ndof->at(iD)>6.) continue;
//TODO		if(trk_pt->at(trk0)<3000 ||
//TODO		   trk_pt->at(trk1)<3000 ||
//TODO		   trk_pt->at(trk2)<3000) continue;

		TVector3 pd(lc_px->at(iD),lc_py->at(iD),lc_pz->at(iD));
		TVector3 xd(lc_x ->at(iD) - pvr_x->at(lc_idx_pvr->at(iD)),
		            lc_y ->at(iD) - pvr_y->at(lc_idx_pvr->at(iD)),
		            lc_z ->at(iD) - pvr_z->at(lc_idx_pvr->at(iD)));
		
		//if(pd.Dot(xd)/pd.Mag()/xd.Mag() < 0.9999) continue;

		//PID and pT requirements 
		//PID off for MC
		if(!isMC_ &&
		   !(trk_pnn_pi->at(trk2)>0.2 &&
		     trk_pnn_k->at( trk1)>0.3 && 
		     trk_pnn_p->at( trk0)>0.3)) continue;
		if(!(pd.Pt()>2000)) continue;

		LCM      ->push_back(lc_m->at(iD));
		LCPX     ->push_back(lc_px->at(iD));
		LCPY     ->push_back(lc_py->at(iD));
		LCPZ     ->push_back(lc_pz->at(iD));
		LCE      ->push_back(lc_e->at(iD));
		LCX      ->push_back(lc_x->at(iD));
		LCY      ->push_back(lc_y->at(iD));
		LCZ      ->push_back(lc_z->at(iD));
		LCIP     ->push_back(lc_ip->at(iD));
		LCIPCHI2 ->push_back(lc_ip_chi2->at(iD));
		LCFD     ->push_back(lc_fd->at(iD));
		LCFDCHI2 ->push_back(lc_fd_chi2->at(iD));
		LCTAU    ->push_back(lc_tau->at(iD));
		LCDIRA   ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		LCVTXCHI2->push_back(lc_vtx_chi2->at(iD));
		LCVTXNDOF->push_back(lc_vtx_ndof->at(iD));
		LCDRJET  ->push_back(dDRj);
		LCDRTAG  ->push_back(dDRt);
		LCPPNN   ->push_back(trk_pnn_p->at( trk0));
		LCKPNN   ->push_back(trk_pnn_k->at( trk1));
		LCPIPNN  ->push_back(trk_pnn_pi->at(trk2));

		LCP      ->push_back(lc_p->at(iD));
		LCPT     ->push_back(lc_pt->at(iD));
		LCPPT    ->push_back(trk_pt->at(trk0));
		LCKPT    ->push_back(trk_pt->at(trk1));
		LCPIPT   ->push_back(trk_pt->at(trk2));

		LCTRK0   ->push_back(trk0);
		LCTRK1   ->push_back(trk1);
		LCTRK2   ->push_back(trk2);

		LCNJ      ->push_back(nj);
		LCMAXDR   ->push_back(maxdr);
	}

	for(uint iD=0; iD<k3pi_m->size(); ++iD) {
		TLorentzVector p4D(k3pi_px->at(iD),
				   k3pi_py->at(iD),
				   k3pi_pz->at(iD),
				   k3pi_e->at( iD));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk;

		int trk0=k3pi_idx_trk0->at(iD);
		int trk1=k3pi_idx_trk1->at(iD);
		int trk2=k3pi_idx_trk2->at(iD);
		int trk3=k3pi_idx_trk3->at(iD);

		nj = (trk_idx_jet->at(trk0)==j) +
		     (trk_idx_jet->at(trk1)==j) +
		     (trk_idx_jet->at(trk2)==j) +
		     (trk_idx_jet->at(trk3)==j);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk0),
				  trk_py->at(trk0),
				  trk_pz->at(trk0),
				  trk_e ->at(trk0));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk1),
				  trk_py->at(trk1),
				  trk_pz->at(trk1),
				  trk_e ->at(trk1));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk2),
				  trk_py->at(trk2),
				  trk_pz->at(trk2),
				  trk_e ->at(trk2));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(trk3),
				  trk_py->at(trk3),
				  trk_pz->at(trk3),
				  trk_e ->at(trk3));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);

		//if(nj == 0) continue;//require at least one track in probe jet
//TODO		if((trk_ip_chi2->at(trk0)>50)+(trk_ip_chi2->at(trk1)>50)+(trk_ip_chi2->at(trk2)>50)+(trk_ip_chi2->at(trk3)>50) == 0 ||
//TODO		   (trk_ip_chi2->at(trk0)>10)+(trk_ip_chi2->at(trk1)>10)+(trk_ip_chi2->at(trk2)>10)+(trk_ip_chi2->at(trk3)>10) <= 1 ||
//TODO		   k3pi_tau->at(iD)<0.00010 || k3pi_vtx_chi2->at(iD)/k3pi_vtx_ndof->at(iD)>6.) continue;

		TVector3 pd(k3pi_px->at(iD),k3pi_py->at(iD),k3pi_pz->at(iD));
		TVector3 xd(k3pi_x ->at(iD) - pvr_x->at(k3pi_idx_pvr->at(iD)),
		            k3pi_y ->at(iD) - pvr_y->at(k3pi_idx_pvr->at(iD)),
		            k3pi_z ->at(iD) - pvr_z->at(k3pi_idx_pvr->at(iD)));
		
		//if(pd.Dot(xd)/pd.Mag()/xd.Mag() < 0.9999) continue;

		//PID and pT requirements 
		//PID off for MC
		if(!isMC_ &&
		   !(trk_pnn_pi->at(trk3)>0.2 &&
		     trk_pnn_pi->at(trk2)>0.2 && 
		     trk_pnn_pi->at(trk1)>0.2 && 
		     trk_pnn_k->at( trk0)>0.3)) continue;
		if(!(pd.Pt()>2000)) continue;

		K3PIM      ->push_back(k3pi_m->at(iD));
		K3PIPX     ->push_back(k3pi_px->at(iD));
		K3PIPY     ->push_back(k3pi_py->at(iD));
		K3PIPZ     ->push_back(k3pi_pz->at(iD));
		K3PIE      ->push_back(k3pi_e->at(iD));
		K3PIX      ->push_back(k3pi_x->at(iD));
		K3PIY      ->push_back(k3pi_y->at(iD));
		K3PIZ      ->push_back(k3pi_z->at(iD));
		K3PIIP     ->push_back(k3pi_ip->at(iD));
		K3PIIPCHI2 ->push_back(k3pi_ip_chi2->at(iD));
		K3PIFD     ->push_back(k3pi_fd->at(iD));
		K3PIFDCHI2 ->push_back(k3pi_fd_chi2->at(iD));
		K3PITAU    ->push_back(k3pi_tau->at(iD));
		K3PIDIRA   ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		K3PIVTXCHI2->push_back(k3pi_vtx_chi2->at(iD));
		K3PIVTXNDOF->push_back(k3pi_vtx_ndof->at(iD));
		K3PIDRJET  ->push_back(dDRj);
		K3PIDRTAG  ->push_back(dDRt);
		K3PIKPNN   ->push_back(trk_pnn_k->at( trk0));
		K3PIPI1PNN ->push_back(trk_pnn_pi->at(trk1));
		K3PIPI2PNN ->push_back(trk_pnn_pi->at(trk2));
		K3PIPI3PNN ->push_back(trk_pnn_pi->at(trk3));

		K3PIP       ->push_back(k3pi_p->at(iD));
		K3PIPT      ->push_back(k3pi_pt->at(iD));
		K3PIKPT     ->push_back(trk_pt->at(trk0));
		K3PIPI1PT   ->push_back(trk_pt->at(trk1));
		K3PIPI2PT   ->push_back(trk_pt->at(trk2));
		K3PIPI3PT   ->push_back(trk_pt->at(trk3));

		K3PITRK0   ->push_back(trk0);
		K3PITRK1   ->push_back(trk1);
		K3PITRK2   ->push_back(trk2);
		K3PITRK3   ->push_back(trk3);

		K3PINJ      ->push_back(nj);
		K3PIMAXDR   ->push_back(maxdr);
	}

	for(uint iD=0; iD<jpsi_m->size(); ++iD) {
		TLorentzVector p4D(jpsi_px->at(iD),
				jpsi_py->at(iD),
				jpsi_pz->at(iD),
				jpsi_e->at( iD));

		double dDRj(0.), dDRt(1.);
		if(j>-1) dDRj = p4D.DeltaR(p4j);
		if(t>-1) dDRt = p4D.DeltaR(p4t);

		if(dDRj>0.5 || dDRt<dDRj) continue;//require D to be close to probe and far from tag

		int nj(0);
		double maxdr(0.);
		TLorentzVector p4Dtrk;

		nj = (trk_idx_jet->at(jpsi_idx_trk0->at(iD))==j) + (trk_idx_jet->at(jpsi_idx_trk1->at(iD))==j);
		p4Dtrk.SetPxPyPzE(trk_px->at(jpsi_idx_trk0->at(iD)),
				  trk_py->at(jpsi_idx_trk0->at(iD)),
				  trk_pz->at(jpsi_idx_trk0->at(iD)),
				  trk_e ->at(jpsi_idx_trk0->at(iD)));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);
		p4Dtrk.SetPxPyPzE(trk_px->at(jpsi_idx_trk1->at(iD)),
				  trk_py->at(jpsi_idx_trk1->at(iD)),
				  trk_pz->at(jpsi_idx_trk1->at(iD)),
				  trk_e ->at(jpsi_idx_trk1->at(iD)));
		if(j>-1) maxdr = TMath::Max(p4Dtrk.DeltaR(p4j), maxdr);

		int trk0=jpsi_idx_trk0->at(iD);
		int trk1=jpsi_idx_trk1->at(iD);

		//if(nj == 0) continue;//require at least one track in probe jet

		TVector3 pd(jpsi_px->at(iD),jpsi_py->at(iD),jpsi_pz->at(iD));
		TVector3 xd(jpsi_x ->at(iD) - pvr_x->at(jpsi_idx_pvr->at(iD)),
		            jpsi_y ->at(iD) - pvr_y->at(jpsi_idx_pvr->at(iD)),
		            jpsi_z ->at(iD) - pvr_z->at(jpsi_idx_pvr->at(iD)));
		
		JPSIM      ->push_back(jpsi_m->at(iD));
		JPSIPX     ->push_back(jpsi_px->at(iD));
		JPSIPY     ->push_back(jpsi_py->at(iD));
		JPSIPZ     ->push_back(jpsi_pz->at(iD));
		JPSIE      ->push_back(jpsi_e->at(iD));
		JPSIX      ->push_back(jpsi_x->at(iD));
		JPSIY      ->push_back(jpsi_y->at(iD));
		JPSIZ      ->push_back(jpsi_z->at(iD));
		JPSIIP     ->push_back(jpsi_ip->at(iD));
		JPSIIPCHI2 ->push_back(jpsi_ip_chi2->at(iD));
		JPSIFD     ->push_back(jpsi_fd->at(iD));
		JPSIFDCHI2 ->push_back(jpsi_fd_chi2->at(iD));
		JPSITAU    ->push_back(jpsi_tau->at(iD));
		JPSIDIRA   ->push_back(pd.Dot(xd)/pd.Mag()/xd.Mag());
		JPSIVTXCHI2->push_back(jpsi_vtx_chi2->at(iD));
		JPSIVTXNDOF->push_back(jpsi_vtx_ndof->at(iD));
		JPSIDRJET  ->push_back(dDRj);
		JPSIDRTAG  ->push_back(dDRt);

		JPSIP      ->push_back(jpsi_p->at(iD));
		JPSIPT     ->push_back(jpsi_pt->at(iD));
		JPSIKPT    ->push_back(trk_pt->at(jpsi_idx_trk0->at(iD)));
		JPSIPIPT   ->push_back(trk_pt->at(jpsi_idx_trk1->at(iD)));

		JPSITRK0   ->push_back(jpsi_idx_trk0->at(iD));
		JPSITRK1   ->push_back(jpsi_idx_trk1->at(iD));

		JPSINJ      ->push_back(nj);
		JPSIMAXDR   ->push_back(maxdr);

		int dg(-1);
		JPSITRUE    ->push_back(checkJpsiTruth(trk_idx_gen->at(trk0),trk_idx_gen->at(trk1),dg));
		int dtrueid(-1);
		int dfromb(-1);
		int ddr(-1);
		if(dg>-1) {
			if(foundDs.find(dg)!=foundDs.end()) dtrueid = foundDs.at(dg);
			dfromb = checkFromB(dg);
			TLorentzVector p4g(gen_px->at(dg),gen_py->at(dg),gen_pz->at(dg),gen_e->at(dg));
			if(j>-1) ddr = p4j.DeltaR(p4g);
		}
		JPSITRUEIDX->push_back(dtrueid);
		JPSIFROMB->push_back(dfromb);
		JPSITRUEDR->push_back(ddr);
	}
	//if no jets then this is the main fill function
	//only fill if we found something
	if(tagType_==NoJets) {
		if(TRUEDID->size()>0||SVM->size()>0||D0M->size()>0||DPM->size()>0||DSM->size()>0||LCM->size()>0||K3PIM->size()>0||JPSIM->size()>0)
			tout->Fill();
	}
}

//function to identify the two jets (j1 & j2) in a dijet MC event
bool makeNewCalibTuples::tagTruthJet(int& j1, int& j2) {
	int ng = gen_pid->size();
	int nj = jet_pz->size();

	int c1(-1), c2(-1), b1(-1), b2(-1), q1(-1), q2(-1);
	TLorentzVector p4c1, p4c2;
	TLorentzVector p4b1, p4b2;
	TLorentzVector p4q1, p4q2;

	int countMatched(false);

	for(int j=0; j<nj; ++j) {
		TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
		bool matched(false);

		TLorentzVector p4jmc;
		//truth match jet
		for(int g=0; g<ng; ++g) {
			if(gen_pid->at(g)!=98) continue;
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));

			if(p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
				////TODO start
				//bool foundC(false);
				//for(int q=0; q<ng; ++q) {
				//	TLorentzVector p4q(gen_px->at(q),gen_py->at(q),gen_pz->at(q),gen_e->at(q));

				//	if(p4q.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4q)<0.5) {
				//		if(TMath::Abs(gen_pid->at(q)) == 4) {
				//			foundC=true;
				//		}
				//	}
				//}
				//std::cout << EVT << "\t" << j << "\t" << p4g.Pt() << "\t" << p4j.Pt() << "\t" << foundC << std::endl;//TODO
				////TODO end
				if(p4g.Pt()>minTruePt_ && p4g.Pt()<maxTruePt_) matched=true;
				p4jmc = p4g;
			}
		}
		if(!matched) continue;
		++countMatched;

		//for light sample - pick jet with highest pT
		if(p4j.Pt()>p4q1.Pt()) {
			p4q1 = p4j;
			q1 = j;
		}

		//for heavy samples - match jets to ccbar or bbbar
		for(int g=0; g<ng; ++g) {
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));

			//require quark to be within DeltaR<0.5 of either true or reco jet
			if((p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) ||
			   (p4g.Pt()>0.&&p4jmc.Pt()>0.&&p4jmc.DeltaR(p4g)<0.5)) {
				if(gen_pid->at(g) == 4 && p4j.Pt() > p4c1.Pt()) {
					p4c1 = p4j;
					c1 = j;
				}
				if(gen_pid->at(g) ==-4 && p4j.Pt() > p4c2.Pt()) {
					p4c2 = p4j;
					c2 = j;
				}
				if(gen_pid->at(g) == 5 && p4j.Pt() > p4b1.Pt()) {
					p4b1 = p4j;
					b1 = j;
				}
				if(gen_pid->at(g) ==-5 && p4j.Pt() > p4b2.Pt()) {
					p4b2 = p4j;
					b2 = j;
				}
			}
		}
	}

	if(countMatched==1) ++tagCounts[2];
	if(countMatched==2) ++tagCounts[3];
	if(countMatched>2)  ++tagCounts[4];

	//veto heavier flavours
	if(flavour_<5 && (b1>-1 || b2>-1)) return false;
	if(flavour_<4 && (c1>-1 || c2>-1)) return false;

	if(flavour_==5) {
		if(b1==-1 && b2==-1) return false;//TODO
		//if(b1==-1 || b2==-1) return false;
		//if(TMath::Abs(p4b1.DeltaPhi(p4b2))>2.0) {
			j1 = b1;
			j2 = b2;
			return true;
		//}
	} else if(flavour_==4) {
		if(c1==-1 && c2==-1) return false;//TODO
		//TODO if(c1==-1 || c2==-1) return false;
		//if(TMath::Abs(p4c1.DeltaPhi(p4c2))>2.0) {
			j1 = c1;
			j2 = c2;
			return true;
		//}
	} else {
		if(q1==-1) return false;
		for(int j=0; j<nj; ++j) {
			TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
			if(TMath::Abs(p4j.DeltaPhi(p4q1))>2.0 && p4j.Pt() > p4q2.Pt()) {
				p4q2 = p4j;
				q2 = j;
			}
		}
		if(q2==-1) return false;
		j1=q1;
		j2=q2;
		return true;
	}


	return false;
}

bool makeNewCalibTuples::tagTruthNoMuJet(int& j1) {
	int ng = gen_pid->size();
	int nj = jet_pz->size();

	int c1(-1), b1(-1), q1(-1);
	TLorentzVector p4c1;
	TLorentzVector p4b1;
	TLorentzVector p4q1;

	int countMatched(false);

	for(int j=0; j<nj; ++j) {
		TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
		bool matched(false);

		TLorentzVector p4jmc;
		//truth match jet
		for(int g=0; g<ng; ++g) {
			if(gen_pid->at(g)!=98) continue;
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));

			if(p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
				if(p4g.Pt()>minTruePt_ && p4g.Pt()<maxTruePt_) matched=true;
				p4jmc = p4g;
			}
		}
		if(!matched) continue;

		//veto if there's a muon from either a Z or W decay within 0.5 of the jet
		for(int g=0; g<ng; ++g) {
			if(TMath::Abs(gen_pid->at(g))!=13) continue;
			if(TMath::Abs(gen_prnt_pid->at(g))!=23 && TMath::Abs(gen_prnt_pid->at(g))!=24) continue;
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));

			//unmatch if we're too close to the muon
			if(p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
				matched = false;
			}
		}
		if(!matched) continue;
		++countMatched;

		//for light sample - pick jet with highest pT
		if(p4j.Pt()>p4q1.Pt()) {
			p4q1 = p4j;
			q1 = j;
		}

		//for heavy samples - match jets to c or b
		for(int g=0; g<ng; ++g) {
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));

			//require quark to be within DeltaR<0.5 of either true or reco jet
			if((p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) ||
			   (p4g.Pt()>0.&&p4jmc.Pt()>0.&&p4jmc.DeltaR(p4g)<0.5)) {
				if(TMath::Abs(gen_pid->at(g)) == 4 && p4j.Pt() > p4c1.Pt()) {
					p4c1 = p4j;
					c1 = j;
				}
				if(TMath::Abs(gen_pid->at(g)) == 5 && p4j.Pt() > p4b1.Pt()) {
					p4b1 = p4j;
					b1 = j;
				}
			}
		}
	}

	if(countMatched==1) ++tagCounts[2];
	if(countMatched==2) ++tagCounts[3];
	if(countMatched>2)  ++tagCounts[4];

	//veto heavier flavours
	if(flavour_<5 && b1>-1) return false;
	if(flavour_<4 && c1>-1) return false;

	if(flavour_==5) {
		if(b1==-1) return false;
		j1 = b1;
		return true;
	} else if(flavour_==4) {
		if(c1==-1) return false;
		j1 = c1;
		return true;
	} else {
		if(q1==-1) return false;
		j1=q1;
		return true;
	}


	return false;
}

//function to identify the probe (j) and tag (t) jets in a dijet data event
bool makeNewCalibTuples::tagSVJet(int& j, int& t) {
	std::vector<int> trigs;
	std::vector<int> tags;
	std::vector<int> probes;
	std::vector<double> dRs;
	std::vector<double> pTAs;

	bool trigFound(false);

	for(unsigned int trig=0; trig<evt_dec->size(); ++trig) {
		//backwards samples do not need to be flavour-enriched so don't require SV trigger
		if(!backwards_ &&
		   evt_dec->at(trig)!=1 &&
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

		//backwards samples do not need to be flavour-enriched so don't require tag SV
		bool tag1 = (evt_j1_idx->at(trig)>=0 && (backwards_ || evt_j1_nsv->at(trig)>0));
		bool tag2 = (evt_j2_idx->at(trig)>=0 && (backwards_ || evt_j2_nsv->at(trig)>0));

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
	DEC=trigs[whichTrig];
	j = probes[whichTrig];
	t = tags[whichTrig];

	return true;
}

//function to identify the probe (j) and tag (t) jets in a J/psi + jet data event
//if the J/psi has a nearby jet then return that as the tag jet, otherwise 0
bool makeNewCalibTuples::tagJpsi(int& j, int& t) {
	std::vector<int> jpsi;
	std::vector<int> tags;
	std::vector<int> probes;

	for(unsigned int jp=0; jp<jpsi_pt->size(); ++jp) {
		if(TMath::Abs(jpsi_m->at(jp) - 3096) > 50.) continue;

		int iTag(-1);
		int iJet(-1);
		double probePt(0.), tagPt(0.);

		TLorentzVector p4Jpsi(jpsi_px->at(jp),
				      jpsi_py->at(jp),
				      jpsi_pz->at(jp),
				      jpsi_e->at( jp));

		if(!(p4Jpsi.Pt()>0.)) continue;//protect against odd entries

		for(unsigned int j=0; j<jet_pt->size(); ++j) {
			TLorentzVector p4Jet(jet_px->at(j),
					     jet_py->at(j),
					     jet_pz->at(j),
					     jet_e->at( j));

			if(!(p4Jet.Pt()>0.)) continue;//protect against odd entries

			//require dPhi>2.0 to associate probe jet
			if(p4Jet.DeltaPhi(p4Jpsi) > 2.0) {
				if(p4Jet.Pt() > probePt) {
					iJet = j;
					probePt = p4Jet.Pt();
				}

			//require dR<0.5 to associate tag jet
			} else if(p4Jet.DeltaR(p4Jpsi) < 0.5) {
				if(p4Jet.Pt() > tagPt) {
					iTag = j;
					tagPt = p4Jet.Pt();
				}
			}
		}

		if(iJet>-1) {
			jpsi.push_back(jp);
			tags.push_back(iTag);
			probes.push_back(iJet);
		}
	}

	if(jpsi.size()==0) return false;

	int whichTrig = gRandom->Integer(jpsi.size());
	j = probes[whichTrig];
	t = tags[whichTrig];

	return true;
}

void makeNewCalibTuples::Loop(int nmax)
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	int firstEntry(0), lastEntry(nentries);
	if(part_>=0 && part_<10) {
		firstEntry = part_*nentries/10;
		lastEntry  = (part_+1)*nentries/10;
	}
	if(nmax>=0 && nmax<(lastEntry-firstEntry)) {
		lastEntry = firstEntry+nmax;
	}

	boost::progress_display progress( lastEntry-firstEntry );
	for (Long64_t jentry=firstEntry; jentry<lastEntry;jentry++) {
		++progress;
		fChain->GetEntry(jentry);
		clearEventOutputs();
		EVT=jentry;

		++tagCounts[0];
		if(evt_pvr_n!=1) continue;//keep to deal with topo bug
		NPV = evt_pvr_n;
		++tagCounts[1];

		//select the jets to use
		int jet1(-1), jet2(-1);
		
		switch(tagType_) {
			case NoTag: //fill all jets without a tag
				for(unsigned int j=0; j<jet_pt->size(); ++j) {
					fillOutput(j,-1);
				}
				break;
			case NoJets: //if non-jet sample then just fill the SV, D, J/psi containers
				fillSVCands(-1,-1);
				break;
			case TruthTag: //if MC then identify the dijet and keep both
				if(tagTruthJet(jet1,jet2)) {
					if(jet1>-1) fillOutput(jet1,jet2);
					if(jet2>-1) fillOutput(jet2,jet1);
				}
				fillTruthOutput();
				break;
			case TruthNoMuTag: //if W,Z + jet MC then identify jet and keep those without a muon from the boson decay
				if(tagTruthNoMuJet(jet1)) {
					if(jet1>-1) fillOutput(jet1,-1);
				}
//				fillTruthOutput();
				break;
			case JetTag: //if dijet data then tag on one jet and keep the other
				if(tagSVJet(jet1, jet2)) fillOutput(jet1,jet2);
				break;
			case JpsiTag: //if J/psi+jet data the tag on J/psi and keep jet
				if(tagJpsi(jet1, jet2)) fillOutput(jet1,jet2);
				break;
		}


	}//end loop over entries

	for(int i=0; i<15; ++i) {
		std::cout << tagCounts[i] << "\t";
	}
	std::cout << std::endl;

	tout->AutoSave();
	fout->Close();
}

int main(int argc, char** argv) {
	int year(0), part(-1), nmax(-1);
	if(argc>1) year = atoi(argv[1]);
	if(argc>2) part = atoi(argv[2]);
	if(argc>3) nmax = atoi(argv[3]);

	makeNewCalibTuples a(year,part);
	a.Loop(nmax);
	return 0;
}
