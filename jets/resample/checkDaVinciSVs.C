#define checkDaVinciSVs_cxx
#include "checkDaVinciSVs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include <iostream>

//interpolate value of histogram as linear function between bin centres
double interpolate(TH1* hist, double x) {
	Int_t xbin = hist->FindBin(x);
	Double_t x0,x1,y0,y1;

	if(x<=hist->GetBinCenter(1)) {
		//below first bin centre - extrapolate from bin centre down to low edge
		//force function to zero at low edge of first bin
		y0 = 0;
		x0 = hist->GetBinLowEdge(1);
		y1 = hist->GetBinContent(1);
		x1 = hist->GetBinCenter(1);
	} else if(x>=hist->GetBinCenter(hist->GetNbinsX())) {
		//above final bin centre - extrapolate up from bin centre to high edge
		//force function to zero at high edge of final bin
		y0 = hist->GetBinContent(hist->GetNbinsX());
		x0 = hist->GetBinCenter(hist->GetNbinsX());
		y1 = 0;
		x1 = hist->GetBinLowEdge(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX());
	} else {
		//for all other points interpolate between nearest two bin centres
		//Note the factor hist->GetBinWidth(xbin)/hist->GetBinWidth(xbin-1) accounts for different bin widths
		if(x<=hist->GetBinCenter(xbin)) {
			y0 = hist->GetBinContent(xbin-1)*hist->GetBinWidth(xbin)/hist->GetBinWidth(xbin-1);
			x0 = hist->GetBinCenter(xbin-1);
			y1 = hist->GetBinContent(xbin);
			x1 = hist->GetBinCenter(xbin);
		} else {
			y0 = hist->GetBinContent(xbin);
			x0 = hist->GetBinCenter(xbin);
			y1 = hist->GetBinContent(xbin+1)*hist->GetBinWidth(xbin)/hist->GetBinWidth(xbin+1);
			x1 = hist->GetBinCenter(xbin+1);
		}
	}
	return y0 + (x-x0)*((y1-y0)/(x1-x0));
}

//get random number from histogram with interpolation between bin centres
double getRandom(TH1* hist, TRandom3* r) {
	//first pick a bin based on integrated rate of each bin
	//Note we could just do
	//  double rho=hist->GetRandom();
	//  int bin = hist->FindBin(rho);
	//but prefer to ensure we're using the same TRandom throughout
	double rho(0);
	//get the maximum for MC method
	//at this point we're just choosing a bin so no need to normalise to bin width
	double max = hist->GetMaximum();
	int bin(0);
	while(true) {
		//throw bin numbers until we win
		bin = r->Rndm()*hist->GetNbinsX();
		if(hist->GetBinContent(bin) > r->Rndm()*max) break;
	}

	//now interpolate rho within the bin
	double binLo = hist->GetBinLowEdge(bin);
	double binHi = binLo + hist->GetBinWidth(bin);
	//get the maximum for MC method
	//linear between bin centres so just check 3 centres
	max = hist->GetBinContent(bin);
	if(interpolate(hist,binLo) > max) max = interpolate(hist,binLo);
	if(interpolate(hist,binHi) > max) max = interpolate(hist,binHi);
	while(true) {
		//throw numbers within the bin until we win
		rho=binLo + r->Rndm()*(binHi-binLo);
		if(interpolate(hist,rho) > r->Rndm()*max) break;
	}
	return rho;
}

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
// Use this one! Verified against the second method and documented here:
// http://geomalgorithms.com/a07-_distance.html
double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
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

//Calculate the offset in x and y (from point B) when point A is propagated along direction A to the Z-plane of point B
TVector3 calcDxDy(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	double tx = da.X() / da.Z();
	double ty = da.Y() / da.Z();

	double dx = pa.X() + (pb.Z() - pa.Z())*tx - pb.X();
	double dy = pa.Y() + (pb.Z() - pa.Z())*ty - pb.Y();

	return TVector3(dx,dy,0.);
}

//Calculate distance of closest approach of line given by point A and direction A to point B
//Equation from Wolfram Alpha
double calcDocaPoint(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	TVector3 x0 = pb;
	TVector3 x1 = pa;
	TVector3 x2 = pa + da;

	return (x0-x1).Cross(x0-x2).Mag()/(x2-x1).Mag();
}

int checkDaVinciSVs::getNSharedVeloHits(int idxi, int idxj) {
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

int checkDaVinciSVs::getNSharedEarlyVeloHits(int idxi, int idxj) {
	int nshared(0);
	if(trk_vhit0->at(idxi) != -1 && trk_vhit0->at(idxi) == trk_vhit0->at(idxj)) ++nshared;
	if(trk_vhit1->at(idxi) != -1 && trk_vhit1->at(idxi) == trk_vhit1->at(idxj)) ++nshared;
	if(trk_vhit2->at(idxi) != -1 && trk_vhit2->at(idxi) == trk_vhit2->at(idxj)) ++nshared;
	if(trk_vhit3->at(idxi) != -1 && trk_vhit3->at(idxi) == trk_vhit3->at(idxj)) ++nshared;

	return nshared;
}

int checkDaVinciSVs::findBestGenPvrForTrk(int idx) {
	int pvr_idx(-1);

     	//if we have truth-level information try to get the PV from this.
	int gen_idx = trk_idx_gen->at(idx);

     	if(gen_idx>=0) {
     		pvr_idx=gen_idx_pvr->at(gen_idx);
     	}
     	//next use the fact that true PVs are stored before the corresponding reconstructed PV and are close in Z
     	if(pvr_idx<0) {
     		int rec_pvr_idx = trk_idx_pvr->at(idx);
     		double minDist=2.;
     		for(int ipvr=0; ipvr<rec_pvr_idx; ++ipvr) {
     			double dist = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(rec_pvr_idx));
     			if(dist < minDist) {
     				minDist=dist;
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
     		double minDOCA=9999.;
     		TVector3 pp(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));
     		TVector3 dp(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
     		for(uint ipvr=0; ipvr<pvr_z->size(); ++ipvr) {
     			TVector3 pv(pvr_x->at(ipvr),pvr_y->at(ipvr),pvr_z->at(ipvr));
     			double doca = calcDocaPoint(pp,dp,pv);
     			if(doca<minDOCA) {
     				bestpvr = ipvr;
     				minDOCA = doca;
     			}
     		}
     		pvr_idx=bestpvr;
     		//this might be reconstructed so look for close pv with lower index
     		double minDist(2.);
     		for(int ipvr=0; ipvr<bestpvr; ++ipvr) {
     			double dist = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(bestpvr));
     			if(dist < minDist) {
     				minDist=dist;
     				pvr_idx = ipvr;
     			}
     		}
     	}
     	if(pvr_idx<0) {
     		//if we get here then the initial value of minDOCA might be too small
     		std::cout << "WARNING - shouldn't reach here!" << std::endl;
     	}
	return pvr_idx;
}

int checkDaVinciSVs::findBestPvrForSvr(TVector3 x, TVector3 p, bool gen) {
	//find best pvr for our svr
	int pvr_idx(-1);
	int bestpvr(-1);
	double minDOCA(9999.);
	for(uint ipvr=0; ipvr<pvr_z->size(); ++ipvr) {
		TVector3 pv(pvr_x->at(ipvr),pvr_y->at(ipvr),pvr_z->at(ipvr));
		double doca = calcDocaPoint(x,p,pv);
		if(doca<minDOCA) {
			bestpvr = ipvr;
			minDOCA = doca;
		}
	}
	pvr_idx=bestpvr;

	double minDist(2.);
	if(gen) {
		//this might be reconstructed so look for close pv with lower index
		for(int ipvr=0; ipvr<bestpvr; ++ipvr) {
			double dist = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(bestpvr));
			if(dist < minDist) {
				minDist=dist;
				pvr_idx = ipvr;
			}
		}
	} else {
		//this might be generated so look for close pv with higher index
		for(uint ipvr=bestpvr+1; ipvr<pvr_z->size(); ++ipvr) {
			double dist = TMath::Abs(pvr_z->at(ipvr) - pvr_z->at(bestpvr));
			if(dist < minDist) {
				minDist=dist;
				pvr_idx = ipvr;
			}
		}
	}
	if(pvr_idx<0) {
		std::cout << "WARNING - shouldn't reach here! (2)" << std::endl;
	}

	return pvr_idx;
}


//main loop
void checkDaVinciSVs::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	int foundNBSV(0), foundSV(0), foundTrks(0), totalSVs(0);

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//		if(_irep<0 && jentry%1000==0) std::cout << jentry << " of " << nentries << std::endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

//		svsBefore += svr_z->size();
		for(uint isvr=0; isvr<svr_z->size(); ++isvr) {
			std::vector<int> trks;
			std::vector<TVector3> trkhit, trkdir;
			int nExpected(0), nFound(0);
//			int pvr_idx = svr_idx_pvr->at(isvr);
			
			//int jet_idx = svr_idx_jet->at(isvr);
			//int pvr_idx = svr_idx_pvr->at(isvr);
			//if(jet_idx>=0 && pvr_idx>=0) {
			//	TLorentzVector pjet(jet_px->at(jet_idx),jet_py->at(jet_idx),jet_pz->at(jet_idx), jet_e->at(jet_idx));
			//	TLorentzVector xsvr(svr_x->at(isvr),svr_y->at(isvr),svr_z->at(isvr), 0);
			//	TLorentzVector xpvr(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx), 0);
			//	TLorentzVector fv = xsvr-xpvr;
			//	double phi1(pjet.Phi()), phi2(fv.Phi());
			//	double dPhi(phi1 - phi2);
			//	while (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
			//	while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
			//	double dEta(pjet.Eta()- fv.Eta());
			//	double dR = TMath::Sqrt(dPhi * dPhi + dEta * dEta);
			//	std::cout << dR <<"\t"<< pjet.DeltaR(fv) <<"\t"<< svr_jet_dr->at(isvr) << std::endl;
			//}
			//continue;

			if(svr_idx_trk0->at(isvr)!=-1) {
				trks.push_back(svr_idx_trk0->at(isvr));
				if(svr_idx_trk1->at(isvr)!=-1) {
					trks.push_back(svr_idx_trk1->at(isvr));
					if(svr_idx_trk2->at(isvr)!=-1) {
						trks.push_back(svr_idx_trk2->at(isvr));
						if(svr_idx_trk3->at(isvr)!=-1) {
							trks.push_back(svr_idx_trk3->at(isvr));
							if(svr_idx_trk4->at(isvr)!=-1) {
								trks.push_back(svr_idx_trk4->at(isvr));
								if(svr_idx_trk5->at(isvr)!=-1) {
									trks.push_back(svr_idx_trk5->at(isvr));
									if(svr_idx_trk6->at(isvr)!=-1) {
										trks.push_back(svr_idx_trk6->at(isvr));
										if(svr_idx_trk7->at(isvr)!=-1) {
											trks.push_back(svr_idx_trk7->at(isvr));
											if(svr_idx_trk8->at(isvr)!=-1) {
												trks.push_back(svr_idx_trk8->at(isvr));
												if(svr_idx_trk9->at(isvr)!=-1) {
													trks.push_back(svr_idx_trk9->at(isvr));
													nExpected=10;
												} else nExpected=9;
											} else nExpected=8;
										} else nExpected=7;
									} else nExpected=6;
								} else nExpected=5;
							} else nExpected=4;
						} else nExpected=3;
					} else nExpected=2;
				} else nExpected=1;
			} else nExpected=0;

			for(uint itrk=0; itrk<trks.size(); ++itrk) {
				//first check the tracks meet our requirements
				//almost all should work by definition, however, some vertices may have non-prompt or non-truth tracks
				int idx = trks[itrk];
				if(trk_vid->at(idx) == 1) continue;
				if(trk_type->at(idx) != 3) continue;
				if(trk_pt->at(idx)<500 || trk_prb_ghost->at(idx)>0.3) continue;

				//if we have truth level information we can separate prompt and non-prompt tracks
				bool isPrompt(false), isTruth(false);
				int gen_idx = trk_idx_gen->at(idx);
				if(gen_idx>-1) isTruth=true;

     				int pvr_idx = findBestGenPvrForTrk(idx);

     				if(isTruth) {
     					if(TMath::Abs(gen_x->at(gen_idx) - pvr_x->at(pvr_idx)) < 1e-5 &&
     							TMath::Abs(gen_y->at(gen_idx) - pvr_y->at(pvr_idx)) < 1e-5 &&
     							TMath::Abs(gen_z->at(gen_idx) - pvr_z->at(pvr_idx)) < 1e-5) isPrompt=true;
     				}
				if(!isPrompt) break;

				//check chi2 requirement
     				//pvr_idx = svr_idx_pvr->at(isvr);

				TVector3 pv(pvr_x->at(pvr_idx), pvr_y->at(pvr_idx), pvr_z->at(pvr_idx));
				TVector3 p(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
				TVector3 x(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));

				TVector3 dxdy = calcDxDy(x, p, pv);
				//double doca = calcDocaPoint(x, p, pv);
				double sigma2 = trk_ip->at(idx)*trk_ip->at(idx) / trk_ip_chi2->at(idx);
				if(sigma2<0) break;
//				if(trk_ip_chi2->at(idx) < 13.) break;
//				if(dxdy.Perp2() / sigma2 < 13.) break;
				if(dxdy.Perp2() / sigma2 < 9.) break;
				//if(doca*doca / sigma2 < 16.) break;
				
				++nFound;

				//now get what we need to remake the vertex
				//get direction and point on track for vertexing
				TVector3 dir(trk_px->at(idx),trk_py->at(idx),trk_pz->at(idx));
				dir.SetMag(1.);
				TVector3 hit(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));

				double stepX = 10*TMath::Tan(dir.Theta())*TMath::Cos(dir.Phi());
				double stepY = 10*TMath::Tan(dir.Theta())*TMath::Sin(dir.Phi());

				while((hit.X()*hit.X()+hit.Y()*hit.Y()) < 6*6){
					hit.SetZ(hit.Z()+10);
					hit.SetX(hit.X()+stepX);
					hit.SetY(hit.Y()+stepY);
				}

				//trks.push_back(idx);
				trkhit.push_back(hit);
				trkdir.push_back(dir);
			}//end loop over tracks

			++totalSVs;
			if(nExpected!=nFound) {
				continue;
				//std::cout << nExpected << "\t" << nFound << std::endl;
			} else {
				++foundTrks;
			}

			//now try to make simple 2-body vertices
			int ntrk = trks.size();
			std::vector<std::pair<int,int> > sv2ij;
			std::vector<TVector3> sv2;
			std::vector<TVector3> sv2fv;

			for( int itrk=0; itrk<ntrk; ++itrk) {
				TLorentzVector p4i(trk_px->at(trks[itrk]), trk_py->at(trks[itrk]), trk_pz->at(trks[itrk]), 0);
				p4i.SetE(TMath::Sqrt(p4i.P()*p4i.P() + 140*140));
				for( int jtrk=itrk+1; jtrk<ntrk; ++jtrk) {
					//if(gen_idx_pvr->at(trk_idx_gen->at(trks[itrk])) != gen_idx_pvr->at(trk_idx_gen->at(trks[jtrk]))) continue;

					TLorentzVector p4j(trk_px->at(trks[jtrk]), trk_py->at(trks[jtrk]), trk_pz->at(trks[jtrk]), 0);
					p4j.SetE(TMath::Sqrt(p4j.P()*p4j.P() + 140*140));

					TLorentzVector p4 = p4i+p4j;
					double m = p4.M();
					TVector3 sv;
					double d = calcDoca(sv, trkhit[itrk], trkdir[itrk], trkhit[jtrk], trkdir[jtrk]);
					
					int pvr_idx = findBestPvrForSvr(sv,p4.Vect(),false);

					if(pvr_idx<0) {
						std::cout << "WARNING - shouldn't reach here! (2)" << std::endl;
						continue;
					}
					//if(pvr_idx != svr_idx_pvr->at(isvr)) continue;
					//pvr_idx = svr_idx_pvr->at(isvr);

					TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
							sv.Y() - pvr_y->at(pvr_idx),
							sv.Z() - pvr_z->at(pvr_idx));
					double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();
					double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);
					double tz = fv.Z()*m/p4.Z()/(3e11)*(1e12);

					//Apply 2-body SV requirements
					if(m < 0. || m > 5000.) continue;
					if(d > 0.2) continue;
					if(mcor<600.) continue;
					if(tz<0. || tz>=10.) continue;
					//TODO add VTXCHI2<10, FDCHI2>25, HITS>0 
					sv2ij.push_back(std::pair<int,int>(trks[itrk],trks[jtrk]));
					sv2.push_back(sv);
					sv2fv.push_back(fv);
				}
			}
			if(sv2.size() > 0) ++foundSV;
			//std::cout << ntrk << "\t" << sv2.size() << std::endl;

			//Now do linking
			std::vector<TVector3> svN; //vertices of linked SVs
			std::vector<std::vector<int> > svNij; //tracks associated with each linked SV

			for(uint ijet=0; ijet<jet_pz->size(); ++ijet) {
				std::set<int> used; // only use each SV once

				//more than one SV - now link
				for(uint s=0; s<sv2.size(); ++s) {
					if(used.count(s)>0) continue;
					used.insert(s);

					//check if this SV is within the jet
					TLorentzVector p4jet;
					p4jet.SetPxPyPzE(jet_px->at(ijet), jet_py->at(ijet), jet_pz->at(ijet), jet_e->at(ijet));
					if(p4jet.DeltaR(TLorentzVector(sv2fv[s],0.)) > 0.5) continue;


					std::vector<int> svTrks;
					std::vector<TVector3> svSv2s;

					TVector3 sv = sv2[s];
					svTrks.push_back(sv2ij[s].first);
					svTrks.push_back(sv2ij[s].second);
					svSv2s.push_back(sv2[s]);

					int nsvCombined(1);
					//int jet_idx(-1);

					//iterate through each track we add and check for additional SVs
					for(uint t=0; t<svTrks.size(); ++t) {
						//loop over remaining SVs
						for(uint ss=s+1; ss<sv2.size(); ++ss) {
							if(used.count(ss)>0) continue;
							if(p4jet.DeltaR(TLorentzVector(sv2fv[ss],0.)) > 0.5) continue;
							bool found(false);
							if(svTrks[t]==sv2ij[ss].first) {
								svTrks.push_back(sv2ij[ss].second);
								found=true;
							}
							if(svTrks[t]==sv2ij[ss].second) {
								svTrks.push_back(sv2ij[ss].first);
								found=true;
							}
							//if found then add this SV to be linked
							if(found) {
								used.insert(ss);
								svSv2s.push_back(sv2[ss]);
								sv += sv2[ss];
								++nsvCombined;
							}
						}
					}
					//sort tracks and remove repeats
					std::sort(svTrks.begin(), svTrks.end());
					svTrks.erase(std::unique(svTrks.begin(),svTrks.end()),svTrks.end());
					//normalise average vertex //TODO weight?
					sv *= (1./nsvCombined);

					TLorentzVector p4;
					std::vector<TLorentzVector> p4ij;

					for(uint itrk=0; itrk<svTrks.size(); ++itrk) {
						TLorentzVector p4i(trk_px->at(svTrks[itrk]), trk_py->at(svTrks[itrk]), trk_pz->at(svTrks[itrk]), 0);
						p4i.SetE(TMath::Sqrt(p4i.P()*p4i.P() + 140*140));
						p4ij.push_back(p4i);
						p4 += p4i;
					}

					//find best pvr for our svr
					int pvr_idx = findBestPvrForSvr(sv,p4.Vect(),false);
					
					if(pvr_idx<0) {
						std::cout << "WARNING - shouldn't reach here! (3)" << std::endl;
						continue;
					}

					//find minimum 2-body vertex radial FD from the assigned PV
					double fd_min(9999.);
					for(uint isv2=0; isv2<svSv2s.size(); ++isv2) {
						double fd = TMath::Sqrt(TMath::Power(svSv2s[isv2].X() - pvr_x->at(pvr_idx),2) +
								TMath::Power(svSv2s[isv2].Y() - pvr_y->at(pvr_idx),2));
						if(fd<fd_min) fd_min = fd;

					}

					//p4.SetPxPyPzE(svr_px->at(isvr), svr_py->at(isvr), svr_pz->at(isvr), svr_e->at(isvr));//TODO

					double m = p4.M();

					TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
							sv.Y() - pvr_y->at(pvr_idx),
							sv.Z() - pvr_z->at(pvr_idx));
					double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();//note this is transverse wrt fd NOT the z-axis
					double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);
					//double tz = fv.Z()*m/p4.Z()/(3e11)*(1e12);

					//TLorentzVector p4jet;
					////loop over jets to find a match
					//double minDeltaR=10.;
					//for(uint ijet=0; ijet<jet_pz->size(); ++ijet) {
					//	TLorentzVector p4jeti;
					//	p4jeti.SetPxPyPzE(jet_px->at(ijet), jet_py->at(ijet), jet_pz->at(ijet), 0);
					//	double deltaR = p4jeti.DeltaR(TLorentzVector(fv,0.));
					//	if(deltaR<minDeltaR) {
					//		if(minDeltaR<0.5) std::cout << "multiple jets found" << std::endl;
					//		minDeltaR = deltaR;
					//		jet_idx = ijet;
					//		p4jet = p4jeti;
					//	}
					//}
					//printf("%.4f\t%d\t%d\n", minDeltaR, jet_idx, static_cast<int>(svr_idx_jet->at(isvr)));

					//no matching jet
					//if(minDeltaR >= 0.5) {
					//	continue;
					//}
					

					//count tracks in jet
					int nTrkInJet(0);
					for(uint itrk=0; itrk<svTrks.size(); ++itrk) {
						double deltaR = p4jet.DeltaR(TLorentzVector(trk_px->at(svTrks[itrk]), trk_py->at(svTrks[itrk]), trk_pz->at(svTrks[itrk]), 0));
						if(deltaR<0.5) ++nTrkInJet;
					}

					//if(svr_pass->at(isvr) == 0) continue;

					//check if we pass n-body requirements
					if(svTrks.size() == 2 && TMath::Abs(m - 500.) < 20.) continue;
					if(fd_min > 15.) continue;
					if(sv.Z() > 200.) continue;
					if(mcor < 600.) continue;
					if(nTrkInJet<2) continue;
					if(p4.Perp() < 2000.) continue;

					//if (!(svr_tau->at(isvr) < 1.5)) continue;
					//if (!(svr_fd_chi2->at(isvr) > 32)) continue;
					//TODO////if (m_parent->m_nbvSelect && !pass) return false;

					//save vertex and list of tracks
					svN.push_back(sv);
					svNij.push_back(svTrks);
				}
			}

			if(svN.size() > 0) ++foundNBSV;
			std::cout << svN.size() << "\t" << jentry << std::endl;

		}

		
	}//end loop over events

	std::cout << foundNBSV << "\t" << foundSV << "\t" << foundTrks << "\t" << totalSVs << std::endl;
}

int main(int argc, char** argv) {
	checkDaVinciSVs a;
	a.Loop();
	return 0;
}
