#define skimTuples_cxx
#include "skimTuples.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>

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
double skimTuples::calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
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

void skimTuples::pairPVs() {
	if(pvsPaired) return;

	lookupTruePV.clear();
	lookupRecoPV.clear();

	std::set<int> genPVset;
	std::set<int> recPVset;

	//first identify true and reco PVs based on which ones are associated with gen and trk objects
	//could also do this by checking if pvr_dz==-1 (gen) or !=-1 (rec)
	for(unsigned int g=0; g<gen_idx_pvr->size(); ++g) {
		genPVset.insert(gen_idx_pvr->at(g));
	}
	for(unsigned int t=0; t<trk_idx_pvr->size(); ++t) {
		recPVset.insert(trk_idx_pvr->at(t));
	}

	//if we have different numbers of true and reco PVs then we need to pad the containers with -1's
	int nPairs = TMath::Max(genPVset.size(),recPVset.size());
	int padGen = nPairs - genPVset.size();
	int padRec = nPairs - recPVset.size();

	//make C-arrays so that we can permute with next_permutation
	int* genPVs = new int[nPairs]();
	int* recPVs = new int[nPairs]();

	std::set<int>::iterator itG = genPVset.begin();
	std::set<int>::iterator itR = recPVset.begin();

	//sets are already ordered and all indices are >-1 so get "lexographically first" permutation by left-padding the sets
	for(int i=0; i<nPairs; ++i) {
		if(i<padGen) {
			genPVs[i] = -1;
		} else {
			genPVs[i] = (*itG);
			++itG;
		}
		if(i<padRec) {
			recPVs[i] = -1;
		} else {
			recPVs[i] = (*itR);
			++itR;
		}
	}

	//now permute the reco array and keep track of the permutation that minimises 
	int* bestMatch = new int[nPairs]();
	double minDist(999999.);

	do {
		//minimise the total significance of the separation of true and reco PVs
		//dist = SUM_i^nPairs [ (truX-recX)**2/deltaX**2 + (truY-recY)**2/deltaY**2 + (truZ-recZ)**2/deltaZ**2 ]**0.5
		double dist(0.);
		for(int i=0; i<nPairs; ++i) {
			//exclude the "fake" or "missed" PVs
			int g = genPVs[i];
			int r = recPVs[i];
			if(g>-1 && r>-1) {
				dist += TMath::Sqrt( TMath::Power((pvr_x->at(g)-pvr_x->at(r))/pvr_dx->at(r),2.)
						   + TMath::Power((pvr_y->at(g)-pvr_y->at(r))/pvr_dy->at(r),2.)
						   + TMath::Power((pvr_z->at(g)-pvr_z->at(r))/pvr_dz->at(r),2.));
			}
		}
		if(dist<minDist) {
			//found a new "best" permutation
			//save the "distance" and the order of recPVs
			minDist = dist;
			std::copy(recPVs,recPVs+nPairs,bestMatch);
		}
	} while(std::next_permutation(recPVs,recPVs+nPairs));

	//once we have a best permutation, fill the lookup maps
	//also map true->true and reco->reco so functions work on all PVs
	for(int i=0; i<nPairs; ++i) {
		int g = genPVs[i];
		int r = bestMatch[i];

		lookupTruePV[g] = g;
		lookupTruePV[r] = g;
		lookupRecoPV[g] = r;
		lookupRecoPV[r] = r;
	}

	pvsPaired=true;
}

int skimTuples::getTruePV(int pv) {
	if(!pvsPaired) pairPVs();

	if(lookupTruePV.find(pv)==lookupTruePV.end()) return -1;
	return lookupTruePV[pv];
}

int skimTuples::getRecoPV(int pv) {
	if(!pvsPaired) pairPVs();

	if(lookupRecoPV.find(pv)==lookupRecoPV.end()) return -1;
	return lookupRecoPV[pv];
}

bool skimTuples::checkTruth(int idxD, std::vector<int>& idcs, int needPID, int needPi, int needK, int needP, int needMu, int needOppPi, int needOppK, int needOppP, bool allowPart) {
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
	if(needPi>-1 && foundPi!=needPi) return false;
	if(needK>-1  && foundK!=needK) return false;
	if(needP>-1  && foundP!=needP) return false;
	if(needMu>-1 && foundMu!=needMu) return false;
	if(needOppPi>-1 && foundOppPi!=needOppPi) return false;
	if(needOppK>-1 && foundOppK!=needOppK) return false;
	if(needOppP>-1 && foundOppP!=needOppP) return false;

	return true;
}

bool skimTuples::checkD0Truth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,421,1,1,0,0);
}

bool skimTuples::checkD0Truth(int idxK, int idxPi, int& idxD, bool ordered) {
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

bool skimTuples::checkDpTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,411,2,1,0,0);
}

bool skimTuples::checkDpTruth(int idxK, int idxPi1, int idxPi2, int& idxD, bool ordered) {
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

bool skimTuples::checkDsTruth(int idxD) {
	std::vector<int> idcs;
	return checkTruth(idxD,idcs,431,1,2,0,0);
}

bool skimTuples::checkDsTruth(int idxK1, int idxK2, int idxPi, int& idxD, bool ordered) {
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

bool skimTuples::checkJpsiTruth(int idxJpsi) {
	std::vector<int> idcs;
	return checkTruth(idxJpsi,idcs,443,0,0,0,2);
}

bool skimTuples::checkJpsiTruth(int idxMu1, int idxMu2, int& idxJpsi) {
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
bool skimTuples::checkSVTruth(int idxD, double minPT) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLived(pid)) return false;
	if(gen_pt->at(idxD)<minPT) return false;

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
bool skimTuples::checkDSVTruth(int idxD, double minPT) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLivedC(pid)) return false;
	if(gen_pt->at(idxD)<minPT) return false;

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
bool skimTuples::checkBSVTruth(int idxD, double minPT) {
	int ng = gen_pid->size();
	int pid = gen_pid->at(idxD);
	if(!longLivedB(pid)) return false;
	if(gen_pt->at(idxD)<minPT) return false;

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

bool skimTuples::checkFromB(int idxD) {
	int ng = gen_pid->size();
	if(idxD<0 || idxD>=ng) return false;

	int idxPrnt = idxD;
	while(idxPrnt!=-1) {
		if(getFlavour(gen_pid->at(idxPrnt))==5) return true;
		if(getFlavour(gen_prnt_pid->at(idxPrnt))==5) return true;
		idxPrnt = gen_idx_prnt->at(idxPrnt);
	}

	return false;
}

int skimTuples::getFlavour(int pid) {
	pid = TMath::Abs(pid);
	
	if(pid<100) return pid;

	int pid1=pid%10000/1000;
	int pid2=pid%1000/100;
	int pid3=pid%100/10;

	return TMath::Max(pid1,TMath::Max(pid2,pid3));
}

bool skimTuples::longLivedB(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==511 || pid==521 || pid==531 || pid==5122 || pid==5132 || pid==5232 || pid==5332 );
}

bool skimTuples::longLivedC(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==411 || pid==421 || pid==431 || pid==4122 || pid==4132 || pid==4232 || pid==4332 );
}

bool skimTuples::longLivedS(int pid) {
	pid = TMath::Abs(pid);
	return ( pid==130 || pid==310 || pid==3122 || pid==3112 || pid==3222 || pid==3312 || pid==3322 );
}

bool skimTuples::longLived(int pid) {
	return (longLivedB(pid) || longLivedC(pid) || longLivedS(pid));
}

int skimTuples::checkRealSV(std::vector<int> indices, int* common) {
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

	if(common) {
		*common = bestAn;
	}

	//if we didn't find anything then return 0
	if(bestAn<0) return 0;
	return nInBestAn;
}

int skimTuples::checkBestD0Cand(int dA, int dB) {
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

void skimTuples::reattachLostParticles() {
	std::set<int> missingParents;
	std::vector<unsigned int> detachedChildren;

	for(unsigned int g=0; g<gen_pid->size(); ++g) {
		if(gen_prnt_pid->at(g)!=-1 && gen_idx_prnt->at(g)==-1) {
			if(getFlavour(gen_pid->at(g)) !=4) continue;//for now only deal with lost Ds

			missingParents.insert(gen_prnt_pid->at(g));
			detachedChildren.push_back(g);
		}
	}

	if(missingParents.size()==0) return;

	for(unsigned int gp=0; gp<gen_pid->size(); ++gp) {
		if(missingParents.find(gen_pid->at(gp)) == missingParents.end()) continue;

		//now find all already matched children and add together 4-mom
		TLorentzVector pp(gen_px->at(gp),gen_py->at(gp),gen_pz->at(gp),gen_e->at(gp));
		TLorentzVector sumpc;

		for(unsigned int gc=0; gc<gen_pid->size(); ++gc) {
			if(gen_idx_prnt->at(gc)==gp) {
				sumpc += TLorentzVector(gen_px->at(gc),gen_py->at(gc),gen_pz->at(gc),gen_e->at(gc));
			}
		}
		//if the 4-momenta add up then nothing is missing
		if(TMath::Abs((sumpc-pp).M())<1.) {
			continue;
		}

		//now try to add one or two of the detached children to see if the 4-momenta add up
		bool foundMatch(false);
		for(unsigned int i=0; i<detachedChildren.size(); ++i) {
			unsigned int gi = detachedChildren[i];
			TLorentzVector pi(gen_px->at(gi),gen_py->at(gi),gen_pz->at(gi),gen_e->at(gi));
			if(TMath::Abs((sumpc+pi-pp).M())<1.) {
				gen_idx_prnt->at(gi) = gp;
				detachedChildren.erase(detachedChildren.begin()+i);
				foundMatch=true;
				break;
			}
			for(unsigned int j=0; j<detachedChildren.size(); ++j) {
				unsigned int gj = detachedChildren[j];
				TLorentzVector pj(gen_px->at(gj),gen_py->at(gj),gen_pz->at(gj),gen_e->at(gj));
				if(TMath::Abs((sumpc+pi+pj-pp).M())<1.) {
					gen_idx_prnt->at(gi) = gp;
					gen_idx_prnt->at(gj) = gp;
					detachedChildren.erase(detachedChildren.begin()+j);
					detachedChildren.erase(detachedChildren.begin()+i);
					foundMatch=true;
					break;
				}
			}
			if(foundMatch) break;
		}
	}
}

//fill counters for PV-matching efficiencies
// 0 - N total true PVs
// 1 - N no truth information
// 2 - N matched PVs of jet
// 3 - N matched PVs of SV
// 4 - N matched PVs of majority of jet tracks
// 5 - N matched PVs of majority of SV tracks
// 6 - N matched PVs of pT-weighted majority of jet tracks
// the true PV is identified by the c-quark
void skimTuples::fillPVEffs(int s, int j) {
	if(j<0 || j>=static_cast<int>(jet_idx_pvr->size()) || s<0 || s>=static_cast<int>(svr_idx_pvr->size())) return;

	//first get truth-level quark for jet
	int gc = matchJetTruth(j, 4, false);

	pvEffsFilled=true;

	if(gc<0) {
		++nPVCount[1];
		return;
	}
	++nPVCount[0];

	int truePV = gen_idx_pvr->at(gc);

	if(getTruePV(jet_idx_pvr->at(j)) == truePV) ++nPVCount[2];
	if(getTruePV(svr_idx_pvr->at(s)) == truePV) ++nPVCount[3];

	std::map<int,int> jetTrkPVs;
	std::map<int,int> svrTrkPVs;
	std::map<int,double> jetTrkPVsWeight;

	for(unsigned int t=0; t<trk_p->size(); t++){
		if(trk_idx_jet->at(t) != j) continue;
		jetTrkPVs[trk_idx_pvr->at(t)]+=1;
		jetTrkPVsWeight[trk_idx_pvr->at(t)]+=trk_pt->at(t);
	}
	for(int i=0; i<10; ++i) {
		int t=svtrk[i]->at(s);
		if(t<0) break;
		svrTrkPVs[trk_idx_pvr->at(t)]+=1;
	}

	int bestPVJ(-1), bestPVS(-1), bestPVJW(-1);
	int mostTrksJ(0), mostTrksS(0), mostTrksJW(0);

	for(std::map<int,int>::iterator it = jetTrkPVs.begin(); it!=jetTrkPVs.end(); ++it) {
		if((*it).second > mostTrksJ) {
			mostTrksJ = (*it).second;
			bestPVJ = (*it).first;
		}
	}
	if(bestPVJ>-1 && getTruePV(bestPVJ) == truePV) ++nPVCount[4];

	for(std::map<int,int>::iterator it = svrTrkPVs.begin(); it!=svrTrkPVs.end(); ++it) {
		if((*it).second > mostTrksS) {
			mostTrksS = (*it).second;
			bestPVS = (*it).first;
		}
	}
	if(bestPVS>-1 && getTruePV(bestPVS) == truePV) ++nPVCount[5];

	for(std::map<int,double>::iterator it = jetTrkPVsWeight.begin(); it!=jetTrkPVsWeight.end(); ++it) {
		if((*it).second > mostTrksJW) {
			mostTrksJW = (*it).second;
			bestPVJW = (*it).first;
		}
	}
	if(bestPVJW>-1 && getTruePV(bestPVJW) == truePV) ++nPVCount[6];

	return;
}

int skimTuples::getSVCategory(int s, int j, std::vector<int> indices) {
	// 1 - no two tracks from same PV
	// 2 - two tracks from same PV but different PV to jet
	// 3 - two tracks from same long-lived decay but different PV to jet
	// 4 - two tracks from same PV as jet but no two from same long-lived decay
	// 5 - two tracks from same long-lived decay but unrelated to jet
	// 6 - two tracks related to jet but not from same long-lived decay
	// 7 - two tracks from same long-lived decay in jet

	//protect against cases where we have no true particles
	if(indices.size()<1) return 0;
	int idxCommon(-1);
	//get true c-quark for jet
	int gc = matchJetTruth(j, 4, false);
	int nPVMatched(0), nJetMatched(0);

	if(gc>-1) {
		for(unsigned int i=0; i<indices.size(); ++i) {
			if(gen_idx_pvr->at(indices[i]) == gen_idx_pvr->at(gc)) ++nPVMatched;
			int prnt = indices[i];
			while(prnt!=-1) {
				if(prnt==gc) {
					++nJetMatched;
					break;
				}
				prnt = gen_idx_prnt->at(prnt);
			}
		}
	}

	//match to same PV either if the reconstructed PV of the SV matches that of the jet or if the true PV of more than one track matches that of the quark
	if(svr_idx_pvr->at(s)!=jet_idx_pvr->at(j) && nPVMatched<2) {
		//possible options: 1,2,3
		if(checkRealSV(indices)<2) {
			//possible options: 1,2
			std::set<int> foundPVs;
			for(int i=0; i<10; ++i) {
				if(svtrk[i]->at(s)<0) break;
				if(!foundPVs.insert(trk_idx_pvr->at(svtrk[i]->at(s))).second) {
					//two tracks with same PV
					return 2;
				}
			}
			//no two tracks with same PV
			return 1;
		} else {
			//possible options: 3
			return 3;
		}

	} else {
		//possible options 4,5,6,7
		if(checkRealSV(indices,&idxCommon)<2) {
			//possible options 4,6
			if(nJetMatched<2) return 4;
			else return 6;
		} else {
			//possible options 5,7
			int prnt = idxCommon;
			while(prnt!=-1) {
				if(prnt==gc) {
					//if we have a real SV and it's the first one in this event then fill counters to determine which method of identifying the PV is best
					if(!pvEffsFilled) fillPVEffs(s,j);
					return 7;
				}
				prnt = gen_idx_prnt->at(prnt);
			}
			return 5;
		}
	}

	//should be unreachable
	return 0;
}

//match a jet to a truth candidate based on dR
//pid defaults to 98 (truth jet) but can be set to +/-4,5 to match to heavy quarks
////matchCharge (default true) require same sign
//useHiPt (default true) prioritises matches with higher pt
//if false, then prioritise lower deltaR
int skimTuples::matchJetTruth(int j, int pid, bool matchCharge, bool useHiPt) {
	TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	int match(-1);
	double ptMatch(0.), drMatch(10.);
	double jpt(0.), gpt(0.), jgdr(0.);
	jpt = p4j.Pt();
	if(jpt<=0.) return -1;
	if(!matchCharge) pid = TMath::Abs(pid);

	int ng = gen_pid->size();
	int gpid(0);
	for(int g=0; g<ng; ++g) {
		gpid = gen_pid->at(g);
		if(!matchCharge) gpid = TMath::Abs(gpid);
		if(gpid!=pid) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		gpt = p4g.Pt();
		jgdr = p4j.DeltaR(p4g);

		if(gpt>0. && jgdr<0.5) {
			if(( useHiPt && gpt>ptMatch) ||
			   (!useHiPt && jgdr<drMatch)) {
				match=g;
				ptMatch = gpt;
				drMatch = jgdr;
			}
		}
	}

	return match;
}

void skimTuples::fillZ(int z, int j) {
	TLorentzVector p4z( z0_px->at(z), z0_py->at(z), z0_pz->at(z), z0_e->at(z));
	TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	TLorentzVector p4mu0( trk_px->at(z0_idx_trk0->at(z)), trk_py->at(z0_idx_trk0->at(z)), trk_pz->at(z0_idx_trk0->at(z)), trk_e->at(z0_idx_trk0->at(z)));
	TLorentzVector p4mu1( trk_px->at(z0_idx_trk1->at(z)), trk_py->at(z0_idx_trk1->at(z)), trk_pz->at(z0_idx_trk1->at(z)), trk_e->at(z0_idx_trk1->at(z)));

	ZM  = p4z.M();
	ZP  = p4z.P();
	ZPX = p4z.Px();
	ZPY = p4z.Py();
	ZPZ = p4z.Pz();
	ZPT = p4z.Pt();
	ZE  = p4z.E();
	ZETA= p4z.Eta();
	ZY  = p4z.Rapidity();
	ZDR = p4z.DeltaR(p4j);

	MU0PX = p4mu0.Px();
	MU0PY = p4mu0.Py();
	MU0PZ = p4mu0.Pz();
	MU0PT = p4mu0.Pt();
	MU0IP = trk_ip->at(z0_idx_trk0->at(z));
	MU0IPCHI2 = trk_ip_chi2->at(z0_idx_trk0->at(z));
	double mu0jet = z0_idx_jet_trk0->at(z);
	if(mu0jet>-1) {
		MU0FPT = p4mu0.Pt()/jet_pt->at(z0_idx_jet_trk0->at(z));
	} else {
		MU0FPT = -1.;
	}
	MU0DR = p4mu0.DeltaR(p4j);

	MU1PX = p4mu1.Px();
	MU1PY = p4mu1.Py();
	MU1PZ = p4mu1.Pz();
	MU1PT = p4mu1.Pt();
	MU1IP = trk_ip->at(z0_idx_trk1->at(z));
	MU1IPCHI2 = trk_ip_chi2->at(z0_idx_trk1->at(z));
	double mu1jet = z0_idx_jet_trk1->at(z);
	if(mu1jet>-1) {
		MU1FPT = p4mu1.Pt()/jet_pt->at(z0_idx_jet_trk1->at(z));
	} else {
		MU1FPT = -1.;
	}
	MU1DR = p4mu1.DeltaR(p4j);
}


void skimTuples::fillTrueZ(int j) {
	int gz(-1), gmu0(-1), gmu1(-1);
	double zpt(0.);

	for(uint g=0; g<gen_pid->size(); ++g) {
		if(gen_pid->at(g)!=23) continue;
		if(gen_pt->at(g)<zpt) continue; //if we find more than one true Z then keep hardest

		gz=g;
		zpt = gen_pt->at(g);

		for(uint gg=0; gg<gen_pid->size(); ++gg) {
			if(gen_pid->at(gg)== 13 && gen_idx_prnt->at(gg)==g) gmu0=gg;
			if(gen_pid->at(gg)==-13 && gen_idx_prnt->at(gg)==g) gmu1=gg;
		}
	}

	if(gz<0) return;

	TLorentzVector p4z(gen_px->at(gz),gen_py->at(gz),gen_pz->at(gz),gen_e->at(gz));
	TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	TLorentzVector p4mu0, p4mu1;

	if(gmu0>-1) p4mu0.SetXYZT( gen_px->at(gmu0), gen_py->at(gmu0), gen_pz->at(gmu0), gen_e->at(gmu0));
	if(gmu1>-1) p4mu1.SetXYZT( gen_px->at(gmu1), gen_py->at(gmu1), gen_pz->at(gmu1), gen_e->at(gmu1));

	ZTRUEM  = p4z.M();
	ZTRUEP  = p4z.P();
	ZTRUEPX = p4z.Px();
	ZTRUEPY = p4z.Py();
	ZTRUEPZ = p4z.Pz();
	ZTRUEPT = p4z.Pt();
	ZTRUEE  = p4z.E();
	ZTRUEETA= p4z.Eta();
	ZTRUEDR = p4z.DeltaR(p4j);

	MU0TRUEPX = p4mu0.Px();
	MU0TRUEPY = p4mu0.Py();
	MU0TRUEPZ = p4mu0.Pz();
	MU0TRUEPT = p4mu0.Pt();

	MU1TRUEPX = p4mu1.Px();
	MU1TRUEPY = p4mu1.Py();
	MU1TRUEPZ = p4mu1.Pz();
	MU1TRUEPT = p4mu1.Pt();
}

void skimTuples::fillOutput(int j, int t)
{
	//fill the output tuples
	TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	TLorentzVector p4sum, p4sumAll;

	std::map<int,int> trkPVs;
	for(unsigned int i=0; i<trk_p->size(); i++){
		if(trk_idx_jet->at(i) != j) continue;
		trkPVs[trk_idx_pvr->at(i)]+=1;
	}
	int mostTrks(0);
	for(std::map<int,int>::iterator it = trkPVs.begin(); it!=trkPVs.end(); ++it) {
		if((*it).second > mostTrks) {
			mostTrks = (*it).second;
			JPV = (*it).first;
		}
	}
	TVector3 pv(pvr_x->at(JPV),pvr_y->at(JPV),pvr_z->at(JPV)); 

	int ihard = -1, imu = -1, nmu = 0, jnchr = 0, jnneu = 0, ndispl6 = 0, ndispl9 = 0, ndispl16 = 0;
	double jptd = 0, jetq = 0, ry = 0, rp = 0, m11 = 0, m12 = 0, m22 = 0, sumpt2 = 0;
	TLorentzVector p4mu,p4hard;
	for(unsigned int i=0; i<trk_p->size(); i++){
		if(trk_idx_jet->at(i) != j) continue;
		jnchr++;
		TLorentzVector p4trk(trk_px->at(i),trk_py->at(i),trk_pz->at(i),trk_e->at(i));
		if(trk_idx_pvr->at(i) == JPV) p4sum+=p4trk;
		p4sumAll+=p4trk;
		int q = trk_q->at(i);
		jetq += q*p4trk.Pt();
		jptd += pow(p4trk.Pt(),2);
		if(p4trk.Pt() > p4hard.Pt()) {p4hard = p4trk; ihard = i;}
		if(trk_is_mu->at(i) > 0 && trk_pnn_mu->at(i) > 0.5 && p4trk.Pt() > 500){
			nmu++;
			if(p4trk.Pt() > p4mu.Pt()) {p4mu = p4trk; imu=i;}
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

	for(unsigned int i=0; i<neu_p->size(); i++){
		if(neu_idx_jet->at(i) != j) continue;
	}

	jetq /= p4j.Pt();
	jptd = sqrt(jptd) / p4j.Pt();      
	for(unsigned int i=0; i<neu_p->size(); i++){
		if(neu_idx_jet->at(i) != j) continue; 
		jnneu++;
		TLorentzVector p4neu(neu_px->at(i),neu_py->at(i),neu_pz->at(i),neu_e->at(i));
		p4sum+=p4neu;
		p4sumAll+=p4neu;
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

	TLorentzVector p4mcsum, p4mcsumAll;
	for(unsigned int g=0; g<gen_pid->size(); ++g) {
		int gpid = fabs(gen_pid->at(g));
		if(gpid==12 || gpid==14 || gpid==16) continue;
		bool decayFound(false);
		for(unsigned int gg=0; gg<gen_pid->size(); ++gg) {
			if(gg==g) continue;
			if(gen_idx_prnt->at(gg) == g) decayFound=true;
		}
		if(decayFound) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(p4g.DeltaR(p4mcj) < 0.5) {
			p4mcsumAll+=p4g;
			if(gen_idx_pvr->at(g) == gen_idx_pvr->at(bestTrueJet)) {
				p4mcsum+=p4g;
			}
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
	JTRUEDPT = 0.;
	JTRUEDDR = 10.;
	int bestGen(-1);
	double bestGenDr(10.);
	for(unsigned int g=0; g<gen_pid->size(); ++g) {
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(p4g.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
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
		JTRUEDPT = gen_pt->at(bestGen);
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

		TLorentzVector p4mct;
		unsigned int bestTrueTag(-1);
		// match to true jet
		for(unsigned int g=0; g<gen_pid->size(); ++g) {
			if(fabs(gen_pid->at(g)) != 98) continue;
			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
			if(p4g.DeltaR(p4t) < 0.5 && p4g.Pt() > p4mct.Pt()) {
				p4mct = p4g;
				bestTrueTag=g;
			}
		}
		if(p4mct.Pt()>0.) {
			TTRUEPX = p4mct.Px();
			TTRUEPY = p4mct.Py();
			TTRUEPZ = p4mct.Pz();
			TTRUEE = p4mct.E();
			TTRUEPT = p4mct.Pt();
			TTRUEETA = p4mct.Eta();
			TTRUEDR = p4t.DeltaR(p4mct);
			usedTruthJets_.insert(bestTrueTag);
		}
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

void skimTuples::fillTruthOutput() {
	for(unsigned int j=0; j<gen_pid->size(); ++j) {
		if(gen_pid->at(j)!=98) continue;
		if(usedTruthJets_.count(j)!=0) continue;

		clearOutputs();

		TLorentzVector p4(gen_px->at(j),gen_py->at(j),gen_pz->at(j),gen_e->at(j));
		if(p4.Pt()<minTruePt_ || p4.Pt()>maxTruePt_) continue;
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
		JTRUEDPT = 0.;
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
		}

		if(bestGen>-1) {
			JTRUEDPX = gen_px->at(bestGen);
			JTRUEDPY = gen_py->at(bestGen);
			JTRUEDPZ = gen_pz->at(bestGen);
			JTRUEDE  = gen_e->at(bestGen);
			JTRUEDPT = gen_pt->at(bestGen);
			JTRUEDDR = bestGenDr;
		}


		tout->Fill();
	}
}

void skimTuples::fillSVCands(int j, int t)
{
	//boolean flags to track SV tagging progress
	bool foundTrueC(false), foundTrueHc(false), foundTrueHc5(false), foundTrueSV(false), 
	     foundRecoSV(false), foundSVPassDR(false), foundSVPassSel(false), foundSVPass4(false), foundSVPassNJ(false);

	if(JTRUEc>0.) foundTrueC=true;
	if(JTRUEDDR<10.) foundTrueHc=true;
	if(JTRUEDPT>5000.) foundTrueHc5=true;
	if(JTRUEDSV>0.) foundTrueSV=true;

	NSV = 0; 
	NTSV = 0; 
	NSVTRK = 0;
	NSVTRUETRK = 0;
	clearOutputVectors();

	TLorentzVector p4j, p4t;
	if(j>-1) p4j.SetPxPyPzE(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
	if(t>-1) p4t.SetPxPyPzE(jet_px->at(t),jet_py->at(t),jet_pz->at(t),jet_e->at(t));

	//map track indices in input tuple to indices in output
	std::map<int,int> jetTrks;

	for(uint t=0; t<trk_pt->size(); ++t) {
		if(j==-1 || trk_idx_jet->at(t) != j) continue;
		jetTrks.insert(std::pair<int,int>(t, TRKPT->size()));
		TRKPT->push_back(trk_pt->at(t));
		TRKIPCHI2->push_back(trk_ip_chi2->at(t));
		TRKINJET->push_back(1.);
	}

	//keep generated B's
	std::map<int,int> foundBs;
	for(uint g=0; g<gen_pid->size(); ++g) {
		if(!longLivedB(gen_pid->at(g))) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(j>-1 && p4j.DeltaR(p4g)>=0.5) continue;
		int pid = gen_pid->at(g);

		foundBs.insert(std::pair<int,int>(g, TRUEBID->size()));

		TRUEBID     ->push_back(pid);
		TRUEBPX     ->push_back(gen_px->at(g));
		TRUEBPY     ->push_back(gen_py->at(g));
		TRUEBPZ     ->push_back(gen_pz->at(g));
		TRUEBPT     ->push_back(gen_pt->at(g));
		TRUEBE      ->push_back(gen_e->at(g));
		//gen_{x,y,z} store ORIGIN coordinates - we want decay coordinates of Hb so find a child
		double gx(-999.), gy(-999.), gz(-999.);
		for(uint gd=0; gd<gen_pid->size(); ++gd) {
			if(gen_idx_prnt->at(gd)!=g) continue;
			gx = gen_x->at(gd);
			gy = gen_y->at(gd);
			gz = gen_z->at(gd);
			break;
		}
		TRUEBX      ->push_back(gx);
		TRUEBY      ->push_back(gy);
		TRUEBZ      ->push_back(gz);
	}

	//keep generated D's
	std::map<int,int> foundDs;
	for(uint g=0; g<gen_pid->size(); ++g) {
		if(!longLivedC(gen_pid->at(g))) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(j>-1 && p4j.DeltaR(p4g)>=0.5) continue;

		std::vector<int> idcs;

		int pid = gen_pid->at(g);
		//find indices
		switch(TMath::Abs(pid)) {
			case 411:
				if(checkTruth(g, idcs,411,2,1)) break;
				idcs.clear();
				break;
			case 421:
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

		//check foundBs for an ancestor
		int bMatch(-1);
		if(foundBs.size()>0) {
			int prnt = g;
			while(prnt!=-1) {
				if(foundBs.find(prnt)!=foundBs.end()) {
					bMatch=foundBs.at(prnt);
					break;
				}
				prnt = gen_idx_prnt->at(prnt);
			}
		}

		TRUEDID     ->push_back(pid);
		TRUEDPX     ->push_back(gen_px->at(g));
		TRUEDPY     ->push_back(gen_py->at(g));
		TRUEDPZ     ->push_back(gen_pz->at(g));
		TRUEDPT     ->push_back(gen_pt->at(g));
		TRUEDE      ->push_back(gen_e->at(g));
		//gen_{x,y,z} store ORIGIN coordinates - we want decay coordinates of Hc so find a child
		double gx(-999.), gy(-999.), gz(-999.);
		for(uint gd=0; gd<gen_pid->size(); ++gd) {
			if(gen_idx_prnt->at(gd)!=g) continue;
			gx = gen_x->at(gd);
			gy = gen_y->at(gd);
			gz = gen_z->at(gd);
			break;
		}
		TRUEDX      ->push_back(gx);
		TRUEDY      ->push_back(gy);
		TRUEDZ      ->push_back(gz);
		TRUEDFROMB  ->push_back(checkFromB(g));
		TRUEDTRUEB  ->push_back(bMatch);
		TRUEDSEL    ->push_back(-1);
		TVector3 p3;
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
	std::set<long> foundSVTrks;
	std::set<long> foundSVTrueTrks;

	int bestSVCat(0);

	for(unsigned int s = 0; s < svr_p->size(); s++){
		if(svr_z->at(s) != svr_z->at(s)) continue;//remove NaN entries

		//From MC studies, the PV assignment of the jet is often wrong O(30%)
		//Calculate flight direction using PV assigned to the SV (correct for O(99%) of simulated SVs)
		TVector3 pv(pvr_x->at(svr_idx_pvr->at(s)),pvr_y->at(svr_idx_pvr->at(s)),pvr_z->at(svr_idx_pvr->at(s))); 

		TVector3 sv(svr_x->at(s),svr_y->at(s),svr_z->at(s));
		TVector3 fly = sv-pv;
		if(backwards_) fly = -fly;

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

		TLorentzVector p4sv;
		p4sv.SetPxPyPzE(0.,0.,0.,0.);
		bool veto=false;
		for(int i=0; i<10; i++){
			if(svtrk[i]->at(s) < 0) break;
			svn++;
			int ii = svtrk[i]->at(s);

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
			svsumipchi2 += ipchi2;
			if(ipchi2 < svminipchi2) svminipchi2 = ipchi2;
			double ghost = trk_prb_ghost->at(ii);
			if(ghost > svmaxghost) svmaxghost = ghost;
		}

		int svCat = getSVCategory(s,j,svtrueindices);
		if(svCat>bestSVCat) bestSVCat = svCat;

		//Now we want to keep tag SVs too (for further tagging)
		//First check if we pass the deltaR requirement for either jet
		//Then check which is closer once we know if we have a good SV
		//If j==-1 && t==-1 then we're not using jets so we keep everything
		foundRecoSV=true;
		if(!((j>-1 && p4j.Vect().DeltaR(fly) < 0.5) || 
		     (t>-1 && p4t.Vect().DeltaR(fly) < 0.5) ||
		     (j==-1 && t==-1))) {
			continue;//reject if too far from jet
		}
		foundSVPassDR=true;

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
		foundSVPassSel=true;

		if(svn>4) continue;//throw away if more than 4 tracks
		foundSVPass4=true;
		//apply 2,3-body additional requirements here too
		//TODO//if(p4sv.M()>5279.4) continue;
		//TODO//double ndof=2*svn-3;
		//TODO//if(svr_chi2->at(s)/ndof>10.) continue;
		//create a unique fingerprint for each set of reco tracks (ignore permutations)
		std::sort(svrecoindices.begin(), svrecoindices.end());
		int fingerprint = 0;
		for(unsigned int itrk=0; itrk<svrecoindices.size(); ++itrk) {
			fingerprint += svrecoindices[itrk]*TMath::Power(trk_pt->size(),static_cast<int>(itrk));
		}
		if(!foundSVs.insert(fingerprint).second) continue;//if insert into set fails then we already have this SV

		//only count tracks in SVs that pass selection
		for(int i=0; i<10; i++){
			if(svtrk[i]->at(s) < 0) break;
			int ii = svtrk[i]->at(s);
			if(trk_idx_gen->at(ii) >= 0) {
				foundSVTrueTrks.insert(trk_idx_gen->at(ii));
			}
			foundSVTrks.insert(ii);
		}

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
				if(trk_pnn_k->at(trk0)>trk_pnn_k->at(trk1)) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					svisd0 = true;
					svd0m = (p0+p1).M();
				} else {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					svisd0 = true;
					svd0m = (p0+p1).M();
				}
			} else {
				if(trk_pnn_pi->at(trk1)>0.2 && trk_pnn_k->at( trk0)>0.3) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),493.7);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),139.6);
					svisd0 = true;
					svd0m = (p0+p1).M();
				} else if(trk_pnn_pi->at(trk0)>0.2 && trk_pnn_k->at( trk1)>0.3) {
					p0.SetXYZM(trk_px->at(trk0),trk_py->at(trk0),trk_pz->at(trk0),139.6);
					p1.SetXYZM(trk_px->at(trk1),trk_py->at(trk1),trk_pz->at(trk1),493.7);
					svisd0 = true;
					svd0m = (p0+p1).M();
				}
			}
		}

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

		int idx0(-1), idx1(-1), idx2(-1), idx3(-1);
		//lookup new indices for tracks
		if(svn>0) {
			if(jetTrks.find(svtrk[0]->at(s))==jetTrks.end()) {
				uint t = svtrk[0]->at(s);
				jetTrks.insert(std::pair<int,int>(t, TRKPT->size()));
				TRKPT->push_back(trk_pt->at(t));
				TRKIPCHI2->push_back(trk_ip_chi2->at(t));
				TRKINJET->push_back(0.);
			}
			idx0 = jetTrks.at(svtrk[0]->at(s));
		}
		if(svn>1) {
			if(jetTrks.find(svtrk[1]->at(s))==jetTrks.end()) {
				uint t = svtrk[1]->at(s);
				jetTrks.insert(std::pair<int,int>(t, TRKPT->size()));
				TRKPT->push_back(trk_pt->at(t));
				TRKIPCHI2->push_back(trk_ip_chi2->at(t));
				TRKINJET->push_back(0.);
			}
			idx1 = jetTrks.at(svtrk[1]->at(s));
		}
		if(svn>2) {
			if(jetTrks.find(svtrk[2]->at(s))==jetTrks.end()) {
				uint t = svtrk[2]->at(s);
				jetTrks.insert(std::pair<int,int>(t, TRKPT->size()));
				TRKPT->push_back(trk_pt->at(t));
				TRKIPCHI2->push_back(trk_ip_chi2->at(t));
				TRKINJET->push_back(0.);
			}
			idx2 = jetTrks.at(svtrk[2]->at(s));
		}
		if(svn>3) {
			if(jetTrks.find(svtrk[3]->at(s))==jetTrks.end()) {
				uint t = svtrk[3]->at(s);
				jetTrks.insert(std::pair<int,int>(t, TRKPT->size()));
				TRKPT->push_back(trk_pt->at(t));
				TRKIPCHI2->push_back(trk_ip_chi2->at(t));
				TRKINJET->push_back(0.);
			}
			idx3 = jetTrks.at(svtrk[3]->at(s));
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
			if(svnt<1) continue;
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
			if(svnj<1) continue;
			foundSVPassNJ=true;
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
			SVVXCHI2     ->push_back(svr_chi2->at(s));
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
			SVTRK0IDX    ->push_back(idx0);
			SVTRK1IDX    ->push_back(idx1);
			SVTRK2IDX    ->push_back(idx2);
			SVTRK3IDX    ->push_back(idx3);
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
	NSVTRK = foundSVTrks.size();
	NSVTRUETRK = foundSVTrueTrks.size();
	//now find number of linked tracks for each SV
	//for each saved SV, loop over the other SVs adding tracks from any with a shared track
	//repeat until the number of tracks remains the same over an iteration
	//this gives a measure of how high the SV track multiplicity would be if "vertex" requirements are ignored
	for(uint s=0; s<NSV; ++s) {
		std::set<int> linkedTrks;
		std::set<int> usedSVs;
		uint nLinked = 0;

		usedSVs.insert(s);
		if(SVN->at(s)>0) linkedTrks.insert(SVTRUETRK0IDX->at(s));
		if(SVN->at(s)>1) linkedTrks.insert(SVTRUETRK1IDX->at(s));
		if(SVN->at(s)>2) linkedTrks.insert(SVTRUETRK2IDX->at(s));
		if(SVN->at(s)>3) linkedTrks.insert(SVTRUETRK3IDX->at(s));

		while(nLinked<linkedTrks.size()) {
			nLinked = linkedTrks.size();

			for(uint ss=0; ss<NSV; ++ss) {
				if(usedSVs.count(ss)>0) continue;
				if((SVN->at(ss)>0 && linkedTrks.count(SVTRUETRK0IDX->at(ss))>0) ||
				   (SVN->at(ss)>1 && linkedTrks.count(SVTRUETRK1IDX->at(ss))>0) ||
				   (SVN->at(ss)>2 && linkedTrks.count(SVTRUETRK2IDX->at(ss))>0) ||
				   (SVN->at(ss)>3 && linkedTrks.count(SVTRUETRK3IDX->at(ss))>0)) {
					if(SVN->at(ss)>0) linkedTrks.insert(SVTRUETRK0IDX->at(ss));
					if(SVN->at(ss)>1) linkedTrks.insert(SVTRUETRK1IDX->at(ss));
					if(SVN->at(ss)>2) linkedTrks.insert(SVTRUETRK2IDX->at(ss));
					if(SVN->at(ss)>3) linkedTrks.insert(SVTRUETRK3IDX->at(ss));
					usedSVs.insert(ss);
				}
			}
		}

		SVNLINKED->push_back(nLinked);
	}
	int bestSVPassLevel(0);
	if(foundTrueC) {
		++svTagCounts[0];
		++bestSVPassLevel;
		if(foundTrueHc) {
			++svTagCounts[1];
			++bestSVPassLevel;
			if(foundTrueHc5) {
				++svTagCounts[2];
				++bestSVPassLevel;
				if(foundTrueSV) {
					++svTagCounts[3];
					++bestSVPassLevel;
					if(foundRecoSV) {
						++svTagCounts[4];
						++bestSVPassLevel;
						if(foundSVPassDR) {
							++svTagCounts[5];
							++bestSVPassLevel;
							if(foundSVPassSel) {
								++svTagCounts[6];
								++bestSVPassLevel;
								if(foundSVPass4) {
									++svTagCounts[7];
									++bestSVPassLevel;
									if(foundSVPassNJ) {
										++svTagCounts[8];
										++bestSVPassLevel;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	svCatHist_->Fill(bestSVCat,bestSVPassLevel);
	

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

		TVector3 pd(dp_px->at(iD),dp_py->at(iD),dp_pz->at(iD));
		TVector3 xd(dp_x ->at(iD) - pvr_x->at(dp_idx_pvr->at(iD)),
		            dp_y ->at(iD) - pvr_y->at(dp_idx_pvr->at(iD)),
		            dp_z ->at(iD) - pvr_z->at(dp_idx_pvr->at(iD)));
		
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

		TVector3 pd(ds_px->at(iD),ds_py->at(iD),ds_pz->at(iD));
		TVector3 xd(ds_x ->at(iD) - pvr_x->at(ds_idx_pvr->at(iD)),
		            ds_y ->at(iD) - pvr_y->at(ds_idx_pvr->at(iD)),
		            ds_z ->at(iD) - pvr_z->at(ds_idx_pvr->at(iD)));
		
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

		TVector3 pd(lc_px->at(iD),lc_py->at(iD),lc_pz->at(iD));
		TVector3 xd(lc_x ->at(iD) - pvr_x->at(lc_idx_pvr->at(iD)),
		            lc_y ->at(iD) - pvr_y->at(lc_idx_pvr->at(iD)),
		            lc_z ->at(iD) - pvr_z->at(lc_idx_pvr->at(iD)));
		
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

		TVector3 pd(k3pi_px->at(iD),k3pi_py->at(iD),k3pi_pz->at(iD));
		TVector3 xd(k3pi_x ->at(iD) - pvr_x->at(k3pi_idx_pvr->at(iD)),
		            k3pi_y ->at(iD) - pvr_y->at(k3pi_idx_pvr->at(iD)),
		            k3pi_z ->at(iD) - pvr_z->at(k3pi_idx_pvr->at(iD)));
		
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
bool skimTuples::tagTruthJet(int& j1, int& j2) {
	int nj = jet_pz->size();

	//vectors to store pairings of light, charm and beauty dijets
	//dijets must share flavour and PV
	//if a jet has no matches save in a pair with -1
	//first in pair is the quark jet, second is anti-quark
	std::vector<std::pair<int,int> > qqDijets;
	std::vector<std::pair<int,int> > ccDijets;
	std::vector<std::pair<int,int> > bbDijets;

	//a jet may appear in multiple pairs but protect against adding single jet if it's already in a pair
	std::set<int> usedJets;

	//number of truth-matched jets
	int countMatched(0);

	for(int j=0; j<nj; ++j) {
		if(jet_idx_pvr->at(j)<0) continue;//protect against weird events where jets aren't matched to PVs
		TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
		if(p4j.Pt()<=0. || p4j.Eta()<2.2 || p4j.Eta()>4.2) continue;

		//truth match jet
		int g = matchJetTruth(j);
//		printf("jet%d %d %f\n",j,g,p4j.Pt());//TODO
		if(g<0) continue;
		TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
		if(p4g.Pt()<minTruePt_ || p4g.Pt()>maxTruePt_) continue;
		++countMatched;

		//try to- match jet to heavy quark
		int b  = matchJetTruth(j, 5);
		int bb = matchJetTruth(j,-5);
		int c  = matchJetTruth(j, 4);
		int cb = matchJetTruth(j,-4);
//		printf("jet%d %d %d %d %d\n",j,b,bb,c,cb);//TODO

		int jetFlav(0);
		if(b>-1) {
			if(bb>-1) {
//				std::cout << "ambiguous flavour b and bbar found" << std::endl;//TODO
				continue;
			}
			jetFlav=5;
		} else if(bb>-1) {
			jetFlav=-5;
		} else if(c>-1) {
			if(cb>-1) {
//				std::cout << "ambiguous flavour c and cbar found" << std::endl;//TODO
				continue;
			}
			jetFlav=4;
		} else if(cb>-1) {
			jetFlav=-4;
		}

		//try to create di-jets
		for(int jj=j+1; jj<nj; ++jj) {
			TLorentzVector p4jj(jet_px->at(jj),jet_py->at(jj),jet_pz->at(jj),jet_e->at(jj));

			int gg = matchJetTruth(jj);
			if(gg<0) continue;
			TLorentzVector p4gg(gen_px->at(gg),gen_py->at(gg),gen_pz->at(gg),gen_e->at(gg));
			if(p4gg.Pt()<minTruePt_ || p4gg.Pt()>maxTruePt_) continue;

			//check PVs match
			if(gen_idx_pvr->at(g) != gen_idx_pvr->at(gg)) continue;
			if(jet_idx_pvr->at(j) != jet_idx_pvr->at(jj)) {
//				std::cout << "true PVs match but not reco PVs" << std::endl;//TODO
				continue;
			}

			//now ensure jet flavours correct
			//if we have a flavour-matched jet pair
			//add to the correct vector and add both to list of used jets
			switch(jetFlav) {
				case -5:
					//require a b-jet
					//veto bbar
					if(matchJetTruth(jj,5)==-1 || matchJetTruth(jj,-5)>-1) continue;
					bbDijets.push_back(std::make_pair(jj,j));
					break;
				case -4:
					//require a c-jet
					//veto cbar, b and bbar
					if(matchJetTruth(jj,4)==-1 || matchJetTruth(jj,-4)>-1 || matchJetTruth(jj,5)>-1 || matchJetTruth(jj,-5)>-1) continue;
					ccDijets.push_back(std::make_pair(jj,j));
					break;
				case 0:
					//require another light jet
					if(matchJetTruth(jj,4)>-1 || matchJetTruth(jj,-4)>-1 || matchJetTruth(jj,5)>-1 || matchJetTruth(jj,-5)>-1) continue;
					qqDijets.push_back(std::make_pair(j,jj));
					break;
				case 4:
					//require a cbar-jet
					//veto c, b and bbar
					if(matchJetTruth(jj,-4)==-1 || matchJetTruth(jj,4)>-1 || matchJetTruth(jj,5)>-1 || matchJetTruth(jj,-5)>-1) continue;
					ccDijets.push_back(std::make_pair(j,jj));
					break;
				case 5:
					//require a bbar-jet
					//veto b
					if(matchJetTruth(jj,-5)==-1 || matchJetTruth(jj,5)>-1) continue;
					bbDijets.push_back(std::make_pair(j,jj));
					break;
				default:
					std::cout << "SHOULD NOT REACH HERE! " << jetFlav << std::endl;
			}

			usedJets.insert(j);
			usedJets.insert(jj);
		}

		//if we didn't find a pair then add single jet
		if(usedJets.count(j)==0) {
			switch(jetFlav) {
				case -5:
					bbDijets.push_back(std::make_pair(-1,j));
					break;
				case -4:
					ccDijets.push_back(std::make_pair(-1,j));
					break;
				case 0:
					qqDijets.push_back(std::make_pair(j,-1));
					break;
				case 4:
					ccDijets.push_back(std::make_pair(j,-1));
					break;
				case 5:
					bbDijets.push_back(std::make_pair(j,-1));
					break;
				default:
					std::cout << "SHOULD NOT REACH HERE! " << jetFlav << std::endl;
			}
			usedJets.insert(j);
		}
	}

	//point dijets at the list we're using
	std::vector<std::pair<int,int> >* dijets(0);
	switch(flavour_) {
		case 0:
			dijets = &qqDijets;
			break;
		case 4:
			dijets = &ccDijets;
			break;
		case 5:
			dijets = &bbDijets;
			break;
		default:
			std::cout << "SHOULD NOT REACH HERE! " << flavour_ << std::endl;
	}
	if(!dijets || dijets->empty()) return false;

	int best(-1);
	double bestPtSum(0.);

	for(unsigned int i=0; i<dijets->size(); ++i) {
		int jet1 = dijets->at(i).first;
		int jet2 = dijets->at(i).second;
		double ptSum(0.);
		if(jet1>-1) ptSum += jet_pt->at(jet1);
		if(jet2>-1) ptSum += jet_pt->at(jet2);
		if(ptSum>bestPtSum) {
			best = i;
			bestPtSum = ptSum;
		}
	}

	j1 = dijets->at(best).first;
	j2 = dijets->at(best).second;
//	printf("debug %d %d %d %d %d %d\n", flavour_, j1, j2, static_cast<int>(qqDijets.size()), static_cast<int>(ccDijets.size()), static_cast<int>(bbDijets.size()));//TODO

	if(countMatched==1) ++tagCounts[2];
	if(countMatched==2) ++tagCounts[3];
	if(countMatched>2)  ++tagCounts[4];
	if(dijets->size()==1) ++tagCounts[5];
	if(dijets->size()>1) ++tagCounts[6];
	if(j1<0 || j2<0) ++tagCounts[7];
	if(j1>-1 && j2>-1) ++tagCounts[8];

	//add trigger information
	for(unsigned int trig=0; trig<evt_dec->size(); ++trig) {
		if(j1>-1 && (evt_j1_idx->at(trig)==j1 || evt_j2_idx->at(trig)==j1)) {
			if(evt_j1_idx->at(trig)==j1) JTRIGPT = TMath::Sqrt(evt_j1_px->at(trig)*evt_j1_px->at(trig)+evt_j1_py->at(trig)*evt_j1_py->at(trig));
			if(evt_j2_idx->at(trig)==j1) JTRIGPT = TMath::Sqrt(evt_j2_px->at(trig)*evt_j2_px->at(trig)+evt_j2_py->at(trig)*evt_j2_py->at(trig));
			if(evt_dec->at(trig)==5) JTRIG10=1.;
			if(evt_dec->at(trig)==0) JTRIG17=1.;
			if(evt_dec->at(trig)==10) JTRIG60=1.;
		}
		if(j2>-1 && (evt_j1_idx->at(trig)==j2 || evt_j2_idx->at(trig)==j2)) {
			if(evt_j1_idx->at(trig)==j2) TTRIGPT = TMath::Sqrt(evt_j1_px->at(trig)*evt_j1_px->at(trig)+evt_j1_py->at(trig)*evt_j1_py->at(trig));
			if(evt_j2_idx->at(trig)==j2) TTRIGPT = TMath::Sqrt(evt_j2_px->at(trig)*evt_j2_px->at(trig)+evt_j2_py->at(trig)*evt_j2_py->at(trig));
			if(evt_dec->at(trig)==5) TTRIG10=1.;
			if(evt_dec->at(trig)==0) TTRIG17=1.;
			if(evt_dec->at(trig)==10) TTRIG60=1.;
		}
	}

	return true;
}
//bool skimTuples::tagTruthJet(int& j1, int& j2) {
//	int ng = gen_pid->size();
//	int nj = jet_pz->size();
//
//	int c1(-1), c2(-1), b1(-1), b2(-1), q1(-1), q2(-1);
//	TLorentzVector p4c1, p4c2;
//	TLorentzVector p4b1, p4b2;
//	TLorentzVector p4q1, p4q2;
//
//	int countMatched(false);
//
//	for(int j=0; j<nj; ++j) {
//		if(jet_idx_pvr->at(j)<0) continue;//protect against weird events where jets aren't matched to PVs
//		TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
//		bool matched(false);
//
//		TLorentzVector p4jmc;
//		//truth match jet
//		for(int g=0; g<ng; ++g) {
//			if(gen_pid->at(g)!=98) continue;
//			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
//
//			if(p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) {
//				if(p4g.Pt()>minTruePt_ && p4g.Pt()<maxTruePt_) matched=true;
//				p4jmc = p4g;
//			}
//		}
//		if(!matched) continue;
//		++countMatched;
//
//		//for light sample - pick jet with highest pT
//		if(p4j.Pt()>p4q1.Pt()) {
//			p4q1 = p4j;
//			q1 = j;
//		}
//
//		//for heavy samples - match jets to ccbar or bbbar
//		for(int g=0; g<ng; ++g) {
//			TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
//
//			//require quark to be within DeltaR<0.5 of either true or reco jet
//			if((p4g.Pt()>0.&&p4j.Pt()>0.&&p4j.DeltaR(p4g)<0.5) ||
//			   (p4g.Pt()>0.&&p4jmc.Pt()>0.&&p4jmc.DeltaR(p4g)<0.5)) {
//				if(gen_pid->at(g) == 4 && p4j.Pt() > p4c1.Pt()) {
//					p4c1 = p4j;
//					c1 = j;
//				}
//				if(gen_pid->at(g) ==-4 && p4j.Pt() > p4c2.Pt()) {
//					p4c2 = p4j;
//					c2 = j;
//				}
//				if(gen_pid->at(g) == 5 && p4j.Pt() > p4b1.Pt()) {
//					p4b1 = p4j;
//					b1 = j;
//				}
//				if(gen_pid->at(g) ==-5 && p4j.Pt() > p4b2.Pt()) {
//					p4b2 = p4j;
//					b2 = j;
//				}
//			}
//		}
//	}
//
//	if(countMatched==1) ++tagCounts[2];
//	if(countMatched==2) ++tagCounts[3];
//	if(countMatched>2)  ++tagCounts[4];
//
//	//veto heavier flavours
//	if(flavour_<5 && (b1>-1 || b2>-1)) return false;
//	if(flavour_<4 && (c1>-1 || c2>-1)) return false;
//
//	if(flavour_==5) {
//		if(b1==-1 && b2==-1) return false;
//			j1 = b1;
//			j2 = b2;
//			return true;
//		//}
//	} else if(flavour_==4) {
//		if(c1==-1 && c2==-1) return false;
//			j1 = c1;
//			j2 = c2;
//			return true;
//		//}
//	} else {
//		if(q1==-1) return false;
//		for(int j=0; j<nj; ++j) {
//			TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
//			if(TMath::Abs(p4j.DeltaPhi(p4q1))>2.0 && p4j.Pt() > p4q2.Pt()) {
//				p4q2 = p4j;
//				q2 = j;
//			}
//		}
//		if(q2==-1) return false;
//		j1=q1;
//		j2=q2;
//		return true;
//	}
//
//
//	return false;
//}

//selection slightly different from Z+jet data (we keep leading jet even if Z not found)
bool skimTuples::tagTruthNoMuJet(int& j1, int& z) {
	int ng = gen_pid->size();
	int nj = jet_pz->size();

	int c1(-1), b1(-1), q1(-1);
	TLorentzVector p4c1;
	TLorentzVector p4b1;
	TLorentzVector p4q1;

	int countMatched(false);

	double ptZ(0.);

	//first get Z if we have it (in the unlikely even we find two candidates, keep high-pT)
	for(unsigned int iz=0; iz<z0_pt->size(); ++iz) {
		int mu0 = z0_idx_trk0->at(iz);
		int mu1 = z0_idx_trk1->at(iz);

		TLorentzVector p4Z(z0_px->at(iz),
				   z0_py->at(iz),
				   z0_pz->at(iz),
				   z0_e->at( iz));
		TLorentzVector p4mu0( trk_px->at(mu0),
				      trk_py->at(mu0),
				      trk_pz->at(mu0),
				      trk_e->at( mu0));
		TLorentzVector p4mu1( trk_px->at(mu1),
				      trk_py->at(mu1),
				      trk_pz->at(mu1),
				      trk_e->at( mu1));

		if(!(p4Z.Pt()>0.)) continue;//protect against odd entries

		//Z selection start
		if(TMath::Abs(z0_m->at(iz) - 90000.) > 30000.) continue;
		if(p4mu0.Pt()<20000. || p4mu1.Pt()<20000.) continue;
		if(p4mu0.Eta()<2. || p4mu0.Eta()>4.5) continue;
		if(p4mu1.Eta()<2. || p4mu1.Eta()>4.5) continue;
		if(TMath::Prob(trk_chi2->at(mu0),trk_ndof->at(mu0))<0.01) continue;
		if(TMath::Prob(trk_chi2->at(mu1),trk_ndof->at(mu1))<0.01) continue;
		if(trk_ip->at(mu0)>0.04 || trk_ip->at(mu1)>0.04) continue;
		//if(trk_dp->at(mu0)/trk_p->at(mu0)<0.1) continue;//TODO off until rerun of DV
		//if(trk_dp->at(mu1)/trk_p->at(mu1)<0.1) continue;//TODO off until rerun of DV
		//if(!(z0_l0_ewmuon_tos->at(iz))) continue; //TODO off until rerun of DV
		//if(!(z0_hlt1_hiptmu_tos->at(iz) && z0_hlt2_hiptmu_tos->at(iz))) continue;
		//Z selection end
		if(p4Z.Pt()>ptZ) {
			z=iz;
			ptZ = p4Z.Pt();
		}
	}


	//select the jet independently
	for(int j=0; j<nj; ++j) {
		if(jet_idx_pvr->at(j)<0) continue;//protect against weird events where jets aren't matched to PVs
		TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),jet_e->at(j));
		bool matched(false);

		//require eta within accepted range
		if(p4j.Pt()<=0. || p4j.Eta()<2.2 || p4j.Eta()>4.2) continue;

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
		//as an extra sanity check - if we have a Z candidate, check for it's muon jets
		if(z>-1 && (j == z0_idx_jet_trk0->at(z) ||
		            j == z0_idx_jet_trk1->at(z))) continue;
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
bool skimTuples::tagSVJet(int& j, int& t) {
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
		if(evt_j1_idx->at(trig)>=jet_pt->size() || evt_j2_idx->at(trig)>=jet_pt->size()) continue;//TODO seems to be required for a few pathological events
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
		if(p4tag.Pt()<10000. || p4tag.Eta()<2.2 || p4tag.Eta()>4.2) continue;

		bool probeFound(false);

		for(unsigned int jprobe=0; jprobe<jet_pt->size(); ++jprobe) {
			if(jtag==jprobe) continue;
			if(jprobe!=evt_j1_idx->at(trig) && jprobe!=evt_j2_idx->at(trig)) continue;//check if the probe jet was part of the trigger
			++tagCounts[7];

			TLorentzVector p4probe(jet_px->at(jprobe),
					       jet_py->at(jprobe),
					       jet_pz->at(jprobe),
					       jet_e->at( jprobe));
			if(p4probe.Pt()<10000. || p4probe.Eta()<2.2 || p4probe.Eta()>4.2) continue;

			double dR = p4probe.DeltaPhi(p4tag);
			double ptAsym = (jet_pt->at(jprobe)-jet_pt->at(jtag)) /
				(jet_pt->at(jprobe)+jet_pt->at(jtag));
	
			if(TMath::Abs(dR)<2.0) continue;
			if(TMath::Abs(ptAsym)>0.25) continue;

			++tagCounts[8];
			if(probeFound) ++tagCounts[9];
			probeFound=true;

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

//function to identify the probe (j) jet and Z in a Z + jet data event
//don't return a tag jet (instead return Z index)
bool skimTuples::tagZJet(int& j, int& z) {
	std::vector<int> zs;
	std::vector<int> probes;

	for(unsigned int iz=0; iz<z0_pt->size(); ++iz) {
		int iJet(-1);
		double probePt(0.);
		int mu0 = z0_idx_trk0->at(iz);
		int mu1 = z0_idx_trk1->at(iz);

		TLorentzVector p4Z(z0_px->at(iz),
				   z0_py->at(iz),
				   z0_pz->at(iz),
				   z0_e->at( iz));
		TLorentzVector p4mu0( trk_px->at(mu0),
				      trk_py->at(mu0),
				      trk_pz->at(mu0),
				      trk_e->at( mu0));
		TLorentzVector p4mu1( trk_px->at(mu1),
				      trk_py->at(mu1),
				      trk_pz->at(mu1),
				      trk_e->at( mu1));

		if(!(p4Z.Pt()>0.)) continue;//protect against odd entries

		//Z selection start
		if(TMath::Abs(z0_m->at(iz) - 90000.) > 30000.) continue;
		if(p4mu0.Pt()<20000. || p4mu1.Pt()<20000.) continue;
		if(p4mu0.Eta()<2. || p4mu0.Eta()>4.5) continue;
		if(p4mu1.Eta()<2. || p4mu1.Eta()>4.5) continue;
		if(TMath::Prob(trk_chi2->at(mu0),trk_ndof->at(mu0))<0.01) continue;
		if(TMath::Prob(trk_chi2->at(mu1),trk_ndof->at(mu1))<0.01) continue;
		if(trk_ip->at(mu0)>0.04 || trk_ip->at(mu1)>0.04) continue;
		//if(trk_dp->at(mu0)/trk_p->at(mu0)<0.1) continue;//TODO off until rerun of DV
		//if(trk_dp->at(mu1)/trk_p->at(mu1)<0.1) continue;//TODO off until rerun of DV
		//if(!(z0_l0_ewmuon_tos->at(iz))) continue; //TODO off until rerun of DV
		if(!(z0_hlt1_hiptmu_tos->at(iz) && z0_hlt2_hiptmu_tos->at(iz))) continue;
		//Z selection end

		for(unsigned int j=0; j<jet_pt->size(); ++j) {
			//skip Z muon "jet"s
			if(j == z0_idx_jet_trk0->at(iz) ||
			   j == z0_idx_jet_trk1->at(iz)) continue;

			TLorentzVector p4Jet(jet_px->at(j),
					     jet_py->at(j),
					     jet_pz->at(j),
					     jet_e->at( j));

			//if(!(p4Jet.Pt()>0.)) continue;//protect against odd entries
			if(p4Jet.Pt()<=0. || p4Jet.Eta()<2.2 || p4Jet.Eta()>4.2) continue;

			//require same PV and separation between jet and Z muons
			if(jet_idx_pvr->at(j) != trk_idx_pvr->at(mu0) || jet_idx_pvr->at(j) != trk_idx_pvr->at(mu1)) continue;
			if(p4Jet.DeltaR(p4mu0)<0.5 || p4Jet.DeltaR(p4mu1)<0.5) continue;
			if(p4Jet.Pt() > probePt) {
				iJet = j;
				probePt = p4Jet.Pt();
			}
		}

		if(iJet>-1) {
			zs.push_back(iz);
			probes.push_back(iJet);
		}
	}

	if(zs.size()==0) return false;

	int whichTrig = gRandom->Integer(zs.size());
	j = probes[whichTrig];
	z = zs[whichTrig];

	return true;
}

//function to identify the probe (j) and tag (t) jets in a J/psi + jet data event
//if the J/psi has a nearby jet then return that as the tag jet, otherwise -1
bool skimTuples::tagJpsi(int& j, int& t) {
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

void skimTuples::Loop(int nmax)
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	int firstEntry(0), lastEntry(nentries);

	if(nmax>=0 && nmax<(lastEntry-firstEntry)) {
		lastEntry = firstEntry+nmax;
	}

//TODO	boost::progress_display progress( lastEntry-firstEntry );
	for (Long64_t jentry=firstEntry; jentry<lastEntry;jentry++) {
//TODO		++progress;
		fChain->GetEntry(jentry);
		clearEventOutputs();
		EVT=jentry;
		pvEffsFilled=false;
		pvsPaired=false;

		//first fix generated particles tree
		reattachLostParticles();

		++tagCounts[0];
		if((tagType_==JetTag || tagType_==TruthTag) && evt_pvr_n!=1) continue;//keep only sinlge PV for dijet
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
				if(tagTruthNoMuJet(jet1,jet2)) {
					if(jet2>-1) fillZ(jet2,jet1);//fill Z branches in output
					fillTrueZ(jet1);//add truth-level information for any Z candidates
					fillOutput(jet1,-1);
				}
//				fillTruthOutput();
				break;
			case JetTag: //if dijet data then tag on one jet and keep the other
				if(tagSVJet(jet1, jet2)) fillOutput(jet1,jet2);
				break;
			case ZTag: //if Z+jet data the tag on Z and keep jet
				if(tagZJet(jet1, jet2)) {
					fillZ(jet2,jet1);//fill Z branches in output
					fillOutput(jet1,-1);
				}
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
	for(int i=0; i<9; ++i) {
		std::cout << svTagCounts[i] << "\t";
	}
	std::cout << std::endl;

	std::ofstream fs(svTagLog,std::ofstream::out);
	for(int i=0; i<9; ++i) {
		fs << svTagCounts[i] << "\t";
	}
	fs << std::endl;
	fs.close();

	fs.open(pvEffLog,std::ofstream::out);
	for(int i=0; i<7; ++i) {
		fs << nPVCount[i] << "\t";
	}
	fs << std::endl;
	fs.close();

	TFile* fhist = TFile::Open(histName,"RECREATE");
	svCatHist_->Write();
	fhist->Close();

	tout->AutoSave();
	if(lumiout) lumiout->AutoSave();
	fout->Close();
}

int main(int argc, char** argv) {
	int year(0), part(-1), nmax(-1);
	TString dir="/tmp/dcraik";
	if(argc>1) year = atoi(argv[1]);
	if(argc>2) part = atoi(argv[2]);
	if(argc>3) dir = argv[3];
	if(argc>4) nmax = atoi(argv[4]);

	skimTuples a(year,part,dir);
	a.Loop(nmax);
	return 0;
}
