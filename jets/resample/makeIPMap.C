#define makeIPMap_cxx
#include "makeIPMap.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//Calculate distance of closest approach of line given by point A and direction A to point B
//Equation from Wolfram Alpha
double calcDocaPoint(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	TVector3 x0 = pb;
	TVector3 x1 = pa;
	TVector3 x2 = pa + da;

	return (x0-x1).Cross(x0-x2).Mag()/(x2-x1).Mag();
}

//Calculate the offset in x and y (from point B) when point A is propagated along direction A to the Z-plane of point B
TVector3 calcDxDy(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	double tx = da.X() / da.Z();
	double ty = da.Y() / da.Z();

	double dx = pa.X() + (pb.Z() - pa.Z())*tx - pb.X();
	double dy = pa.Y() + (pb.Z() - pa.Z())*ty - pb.Y();

	return TVector3(dx,dy,0.);
}

int makeIPMap::getNSharedVeloHits(int idxi, int idxj) {
	int nshared(0);
	if(trk_vhit0->at(idxi) != -1 && trk_vhit0->at(idxi) == trk_vhit0->at(idxj)) ++nshared;
	if(trk_vhit1->at(idxi) != -1 && trk_vhit1->at(idxi) == trk_vhit1->at(idxj)) ++nshared;
	if(trk_vhit2->at(idxi) != -1 && trk_vhit2->at(idxi) == trk_vhit2->at(idxj)) ++nshared;
	if(trk_vhit3->at(idxi) != -1 && trk_vhit3->at(idxi) == trk_vhit3->at(idxj)) ++nshared;
	if(trk_vhit4->at(idxi) != -1 && trk_vhit4->at(idxi) == trk_vhit4->at(idxj)) ++nshared;
	if(trk_vhit5->at(idxi) != -1 && trk_vhit5->at(idxi) == trk_vhit5->at(idxj)) ++nshared;
	if(trk_vhit6->at(idxi) != -1 && trk_vhit6->at(idxi) == trk_vhit6->at(idxj)) ++nshared;
	if(trk_vhit7->at(idxi) != -1 && trk_vhit7->at(idxi) == trk_vhit7->at(idxj)) ++nshared;
	if(trk_vhit8->at(idxi) != -1 && trk_vhit8->at(idxi) == trk_vhit8->at(idxj)) ++nshared;
	if(trk_vhit9->at(idxi) != -1 && trk_vhit9->at(idxi) == trk_vhit9->at(idxj)) ++nshared;
	if(trk_vhit10->at(idxi) != -1 && trk_vhit10->at(idxi) == trk_vhit10->at(idxj)) ++nshared;
	if(trk_vhit11->at(idxi) != -1 && trk_vhit11->at(idxi) == trk_vhit11->at(idxj)) ++nshared;
	if(trk_vhit12->at(idxi) != -1 && trk_vhit12->at(idxi) == trk_vhit12->at(idxj)) ++nshared;
	if(trk_vhit13->at(idxi) != -1 && trk_vhit13->at(idxi) == trk_vhit13->at(idxj)) ++nshared;
	if(trk_vhit14->at(idxi) != -1 && trk_vhit14->at(idxi) == trk_vhit14->at(idxj)) ++nshared;
	if(trk_vhit15->at(idxi) != -1 && trk_vhit15->at(idxi) == trk_vhit15->at(idxj)) ++nshared;
	if(trk_vhit16->at(idxi) != -1 && trk_vhit16->at(idxi) == trk_vhit16->at(idxj)) ++nshared;
	if(trk_vhit17->at(idxi) != -1 && trk_vhit17->at(idxi) == trk_vhit17->at(idxj)) ++nshared;
	if(trk_vhit18->at(idxi) != -1 && trk_vhit18->at(idxi) == trk_vhit18->at(idxj)) ++nshared;
	if(trk_vhit19->at(idxi) != -1 && trk_vhit19->at(idxi) == trk_vhit19->at(idxj)) ++nshared;
	if(trk_vhit20->at(idxi) != -1 && trk_vhit20->at(idxi) == trk_vhit20->at(idxj)) ++nshared;
	if(trk_vhit21->at(idxi) != -1 && trk_vhit21->at(idxi) == trk_vhit21->at(idxj)) ++nshared;
	if(trk_vhit22->at(idxi) != -1 && trk_vhit22->at(idxi) == trk_vhit22->at(idxj)) ++nshared;
	if(trk_vhit23->at(idxi) != -1 && trk_vhit23->at(idxi) == trk_vhit23->at(idxj)) ++nshared;
	if(trk_vhit24->at(idxi) != -1 && trk_vhit24->at(idxi) == trk_vhit24->at(idxj)) ++nshared;
	if(trk_vhit25->at(idxi) != -1 && trk_vhit25->at(idxi) == trk_vhit25->at(idxj)) ++nshared;
	if(trk_vhit26->at(idxi) != -1 && trk_vhit26->at(idxi) == trk_vhit26->at(idxj)) ++nshared;
	if(trk_vhit27->at(idxi) != -1 && trk_vhit27->at(idxi) == trk_vhit27->at(idxj)) ++nshared;
	if(trk_vhit28->at(idxi) != -1 && trk_vhit28->at(idxi) == trk_vhit28->at(idxj)) ++nshared;
	if(trk_vhit29->at(idxi) != -1 && trk_vhit29->at(idxi) == trk_vhit29->at(idxj)) ++nshared;
	if(trk_vhit30->at(idxi) != -1 && trk_vhit30->at(idxi) == trk_vhit30->at(idxj)) ++nshared;
	if(trk_vhit31->at(idxi) != -1 && trk_vhit31->at(idxi) == trk_vhit31->at(idxj)) ++nshared;
	if(trk_vhit32->at(idxi) != -1 && trk_vhit32->at(idxi) == trk_vhit32->at(idxj)) ++nshared;
	if(trk_vhit33->at(idxi) != -1 && trk_vhit33->at(idxi) == trk_vhit33->at(idxj)) ++nshared;
	if(trk_vhit34->at(idxi) != -1 && trk_vhit34->at(idxi) == trk_vhit34->at(idxj)) ++nshared;
	if(trk_vhit35->at(idxi) != -1 && trk_vhit35->at(idxi) == trk_vhit35->at(idxj)) ++nshared;
	if(trk_vhit36->at(idxi) != -1 && trk_vhit36->at(idxi) == trk_vhit36->at(idxj)) ++nshared;
	if(trk_vhit37->at(idxi) != -1 && trk_vhit37->at(idxi) == trk_vhit37->at(idxj)) ++nshared;
	if(trk_vhit38->at(idxi) != -1 && trk_vhit38->at(idxi) == trk_vhit38->at(idxj)) ++nshared;
	if(trk_vhit39->at(idxi) != -1 && trk_vhit39->at(idxi) == trk_vhit39->at(idxj)) ++nshared;
	if(trk_vhit40->at(idxi) != -1 && trk_vhit40->at(idxi) == trk_vhit40->at(idxj)) ++nshared;
	if(trk_vhit41->at(idxi) != -1 && trk_vhit41->at(idxi) == trk_vhit41->at(idxj)) ++nshared;
	if(trk_vhit42->at(idxi) != -1 && trk_vhit42->at(idxi) == trk_vhit42->at(idxj)) ++nshared;
	if(trk_vhit43->at(idxi) != -1 && trk_vhit43->at(idxi) == trk_vhit43->at(idxj)) ++nshared;
	if(trk_vhit44->at(idxi) != -1 && trk_vhit44->at(idxi) == trk_vhit44->at(idxj)) ++nshared;
	if(trk_vhit45->at(idxi) != -1 && trk_vhit45->at(idxi) == trk_vhit45->at(idxj)) ++nshared;
	if(trk_vhit46->at(idxi) != -1 && trk_vhit46->at(idxi) == trk_vhit46->at(idxj)) ++nshared;
	if(trk_vhit47->at(idxi) != -1 && trk_vhit47->at(idxi) == trk_vhit47->at(idxj)) ++nshared;
	if(trk_vhit48->at(idxi) != -1 && trk_vhit48->at(idxi) == trk_vhit48->at(idxj)) ++nshared;
	if(trk_vhit49->at(idxi) != -1 && trk_vhit49->at(idxi) == trk_vhit49->at(idxj)) ++nshared;
	if(trk_vhit50->at(idxi) != -1 && trk_vhit50->at(idxi) == trk_vhit50->at(idxj)) ++nshared;
	if(trk_vhit51->at(idxi) != -1 && trk_vhit51->at(idxi) == trk_vhit51->at(idxj)) ++nshared;
	if(trk_vhit52->at(idxi) != -1 && trk_vhit52->at(idxi) == trk_vhit52->at(idxj)) ++nshared;
	if(trk_vhit53->at(idxi) != -1 && trk_vhit53->at(idxi) == trk_vhit53->at(idxj)) ++nshared;
	if(trk_vhit54->at(idxi) != -1 && trk_vhit54->at(idxi) == trk_vhit54->at(idxj)) ++nshared;
	if(trk_vhit55->at(idxi) != -1 && trk_vhit55->at(idxi) == trk_vhit55->at(idxj)) ++nshared;
	if(trk_vhit56->at(idxi) != -1 && trk_vhit56->at(idxi) == trk_vhit56->at(idxj)) ++nshared;
	if(trk_vhit57->at(idxi) != -1 && trk_vhit57->at(idxi) == trk_vhit57->at(idxj)) ++nshared;
	if(trk_vhit58->at(idxi) != -1 && trk_vhit58->at(idxi) == trk_vhit58->at(idxj)) ++nshared;
	if(trk_vhit59->at(idxi) != -1 && trk_vhit59->at(idxi) == trk_vhit59->at(idxj)) ++nshared;

	return nshared;
}

void makeIPMap::Loop()
{
	gStyle->SetOptStat(0);
	//   In a ROOT session, you can do:
	//      root> .L makeIPMap.C
	//      root> makeIPMap t
	//      root> t.GetEntry(12); // Fill t data members with entry number 12
	//      root> t.Show();       // Show values of entry 12
	//      root> t.Show(16);     // Read and show values of entry 16
	//      root> t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	double maxX(0.), minX(1.), maxY(0.), minY(1.);
	double binsip[61] = {0.,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,6.,8.,10.,20.,30.,75.,150.,500.,5000.};
	//   double binsinvpt[5] = {0.000,0.0005,0.001,0.0015,0.002};
	double binsinvpt[4] = {0.000,0.0005,0.001,0.002};
	//double binstheta[11] = {0.00,0.01,0.02,0.03,0.04,0.05,0.10,0.20,0.30,0.50,3.142};
	//double binsdip[11] = {0.00,0.003,0.01,0.02,0.03,0.05,0.10,0.20,0.30,0.50,1.00};
	double binscosphi[11] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	double binsratio[11] = {0.003,0.01,0.03,0.1,0.3,1.0,3.0,10.0,30.0,100.0,300.0};
	double binsnshare[4] = {0.,2.,4.,20.};
	TH2D prompt(   "prompt",   "", 60, binsip, 3, binsinvpt);
	TH2D nonprompt("nonprompt","", 60, binsip, 3, binsinvpt);
	TH2D truth(    "truth",    "", 60, binsip, 3, binsinvpt);
	TH2D nontruth( "nontruth", "", 60, binsip, 3, binsinvpt);
	TH2D all(      "all",      "", 60, binsip, 3, binsinvpt);
	TH2D promptandnontruth("promptandnontruth","", 60, binsip, 3, binsinvpt);
	TH2D subnonprompt("subnonprompt","", 60, binsip, 3, binsinvpt);

	TH2D deltaDR("deltaDR","", 10, binscosphi, 3, binsnshare);
	TH2D deltaDR2("deltaDR2","", 10, binsratio, 3, binsnshare);

	TH2D phiprompt(   "phiprompt",   "", 10, binscosphi, 3, binsinvpt);

	int pvr_idx(0);
	double tx(0.), ty(0.), dx(0.), dy(0.), sigma(0.), invpt(0.);
	int nSkipped(0), nTruth(0), nPrompt(0), nTotal(0);

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		//store list of prompt tracks for each pvr to perform opening angle study
		std::map<int, std::vector<int> > promptTrksByPV;

		for(int idx=0; idx!= trk_z->size(); ++idx) {
			if(trk_type->at(idx) != 3) continue;
			if(trk_pt->at(idx)<500 || trk_prb_ghost->at(idx)>0.3) continue;
			++nTotal;

			bool isPrompt(false), isTruth(false);
			int gen_idx = trk_idx_gen->at(idx);

			//if we have truth level information we can separate prompt and non-prompt tracks
			if(gen_idx>=0) isTruth=true;

			int pvr_idx(-1);
			//if we have truth-level information try to get the PV from this.
			if(gen_idx>=0) {
				pvr_idx=gen_idx_pvr->at(gen_idx);
			}
			//next use the fact that true PVs are stored before the corresponding reconstructed PV and are close in Z
			if(pvr_idx<0) {
				int rec_pvr_idx = trk_idx_pvr->at(idx);
				double minDist=2.;
				for(int ipvr=0; ipvr<rec_pvr_idx; ++ipvr) {
					double d = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(rec_pvr_idx));
					if(d < minDist) {
						minDist=d;
						pvr_idx = ipvr;
					}
				}
			}
			//if we still haven't found it (no reco pvr associated with track or reco pvr is far from true pvr) maybe there's only one true pvr?
			//if so there are at most 2 pvrs and the true one is at index 0
			if(pvr_idx<0) {
				if (pvr_z->size()<3) pvr_idx=0;
			}
			//otherwise assign the PV with the smallest DOCA
			if(pvr_idx<0) {
				int bestpvr(-1);
				double minDOCA=100.;
				TVector3 pp(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));
				TVector3 dp(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
				for(int ipvr=0; ipvr<pvr_z->size(); ++ipvr) {
					TVector3 pv(pvr_x->at(ipvr),pvr_y->at(ipvr),pvr_z->at(ipvr));
					double doca = calcDocaPoint(pp,dp,pv);
					if(doca<minDOCA) {
						bestpvr = ipvr;
						minDOCA = doca;
//						std::cout << "best doca so far is " << minDOCA << std::endl;
					}
				}
				pvr_idx=bestpvr;
				//this might be reconstructed so look for close pv with lower index
				double minDist(2.);
				for(int ipvr=0; ipvr<bestpvr; ++ipvr) {
					double d = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(bestpvr));
					if(d < minDist) {
						minDist=d;
						pvr_idx = ipvr;
					}
				}
			}
			if(pvr_idx<0) {
				//otherwise, skip this track - should no longer be able to reach here
				++nSkipped;
				continue;
			}
			if(isTruth && TMath::Abs(gen_x->at(gen_idx) - pvr_x->at(pvr_idx)) < 1e-5 &&
				      TMath::Abs(gen_y->at(gen_idx) - pvr_y->at(pvr_idx)) < 1e-5 &&
				      TMath::Abs(gen_z->at(gen_idx) - pvr_z->at(pvr_idx)) < 1e-5) isPrompt=true;

			if(isPrompt || !isTruth) {
				TVector3 pi(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
				TVector3 xi(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));
				TVector3 pv(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx));

				TVector3 dxdyi = calcDxDy(xi, pi, pv);

				for(int jtrk=0; jtrk<promptTrksByPV[pvr_idx].size(); ++jtrk) {
					int idxj = promptTrksByPV[pvr_idx][jtrk];

					int nshared = getNSharedVeloHits(idx,idxj);

					TVector3 pj(trk_px->at(idxj), trk_py->at(idxj), trk_pz->at(idxj));
					TVector3 xj(trk_x->at(idxj), trk_y->at(idxj), trk_z->at(idxj));

					TVector3 dxdyj = calcDxDy(xj, pj, pv);
					double dij = dxdyi.Dot(dxdyj)/(dxdyi.Mag()*dxdyj.Mag());
					double dij2 = dxdyi.Mag()/dxdyj.Mag();

					//double cosTheta(-2.);
					//Double_t ptot2 = pi.Mag2()*pj.Mag2();
					//if(ptot2 > 0) {
					//	cosTheta= pi.Dot(pj)/TMath::Sqrt(ptot2);
					//	if(cosTheta>1.0) cosTheta=1.0;
					//	if(cosTheta<-1.0) cosTheta=-1.0;
					//}
					double theta = pi.Angle(pj);

					deltaDR.Fill(dij,nshared);//theta);
					deltaDR2.Fill(dij2,nshared);//theta);
					deltaDR2.Fill(1./dij2,nshared);//theta);
				}
				promptTrksByPV[pvr_idx].push_back(idx);
			}

			invpt = 1./TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx));

			tx = trk_px->at(idx) / trk_pz->at(idx);
			ty = trk_py->at(idx) / trk_pz->at(idx);

			dx = TMath::Abs(trk_x->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*tx - pvr_x->at(pvr_idx));
			dy = TMath::Abs(trk_y->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*ty - pvr_y->at(pvr_idx));

			sigma = TMath::Sqrt(trk_ip->at(idx)*trk_ip->at(idx) / trk_ip_chi2->at(idx)); //assume unchanged by shift

			if(TMath::Sqrt(dx*dx+dy*dy)/sigma>maxX) maxX=TMath::Sqrt(dx*dx+dy*dy)/sigma;
			if(dx/sigma<minX) minX=dx/sigma;
			if(dy/sigma<minX) minX=dy/sigma;
			if(invpt>maxY) maxY=invpt;
			if(invpt<minY) minY=invpt;

			all.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
			if(isTruth) {
				++nTruth;
				truth.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
				if(isPrompt) {
					++nPrompt;
					prompt.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
					phiprompt.Fill(TMath::ATan2(dy,dx), invpt);
				}
				else nonprompt.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
			}
			else nontruth.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
		}

	}

	promptandnontruth.Add(&prompt,&nontruth);
	subnonprompt.Add(&promptandnontruth,&nonprompt,1.,-((nTotal-nSkipped-nTruth)/static_cast<double>(nTruth)));

	TFile* f = TFile::Open("ipmap-new.root","RECREATE");
	prompt.Write();
	nonprompt.Write();
	truth.Write();
	nontruth.Write();
	all.Write();
	promptandnontruth.Write();
	subnonprompt.Write();
	deltaDR.Write();
	deltaDR2.Write();
	phiprompt.Write();
	f->Close();

	std::cout << maxX << "\t" << minX << "\t" << maxY << "\t" << minY << std::endl;
	std::cout << nSkipped << "\t" << nTruth << "\t" << nPrompt << "\t" << nTotal << std::endl;
	//138.68  5.6146e-06      0.00995349      2.36823e-05
}
