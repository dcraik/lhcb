#define resampleTrackIPs_cxx
#include "resampleTrackIPs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TVector3.h>
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
// The first DOCA method has a bug which I haven't bothered to track down.
// The other two DOCA methods agree so ignore this one for now
//double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
//		               const TVector3 &pb, const TVector3 &db) {
//	if(da.Cross(db).Mag2() == 0) { //special case - lines are parallel
//		if(da.Cross(pa-pb).Mag2() == 0) { // colinear define v as midpoint of pa and pb
//			v = (pa + pb);
//			v *= 0.5;
//			return 0.;
//		} else {
//			//shortest distance from line At + a to At + b
//			//both lines lie in a plane with normal n = Ax(a-b)
//			//shortest distance lies in this plane and is perpendicular to A
//			//solution of the form (a-b) + At* = r(Axn)   (1)
//			//dot with (Axn) to remove At* term
//			//(a-b).(Axn) = r|Axn|^2   (2)
//			TVector3 n = da.Cross(pa-pb);
//			double r = (pa-pb).Dot(da.Cross(n)) / da.Cross(n).Mag2();
//
//			//define midpoint based on pa
//			v = pa - r/2. * da.Cross(n);
//
//			//			std::cout << TMath::Abs(r * da.Cross(n).Mag()) << std::endl;
//			//			v.Print();
//
//			//shortest distance is length of r(Axn)
//			return TMath::Abs(r * da.Cross(n).Mag());
//		}
//	} else {
//		//shortest distance from line At + a to line Bs + b
//		//equivalent to perpendicular distance from plane At - Bs to point (a-b)
//		//solution of form At* - Bs* + (a-b) = r(AxB)   (1)
//		//dot with (AxB) to remove t* and s* terms
//		//(a-b).(AxB) = r(AxB).(AxB)   (2)
//		double r = (pa-pb).Dot(da.Cross(db))/da.Cross(db).Mag2();
//
//		//now find midpoint
//		//dotting (1) with A and B gives
//		//|A|^2 t* -  A.B  s* + (a-b).A = 0   (3)
//		// A.B  t* - |B|^2 s* + (a-b).B = 0   (4)
//		//
//		// (3) - A.B/|B|^2 (4) gives
//		// ( |A|^2 - |A.B|^2/|B|^2 ) t* + (a-b).A - ((a-b).B)(A.B)/|B|^2
//		//
//		// - (4) + A.B/|A|^2 (3) gives
//		// ( |B|^2 - |A.B|^2/|A|^2 ) s* - (a-b).B + ((a-b).A)(A.B)/|A|^2
//		// TODO check fudge factor of -1
//		double tSt = -1*((pa-pb).Dot(da) - (pa-pb).Dot(db)*(da.Dot(db))/db.Mag2()) /
//			(da.Mag2() - da.Dot(db)/db.Mag2());
//
//		double sSt = -1*((pb-pa).Dot(db) - (pb-pa).Dot(da)*(db.Dot(da))/da.Mag2()) /
//			(db.Mag2() - db.Dot(da)/da.Mag2());
//
//		TVector3 va = pa + tSt*da;
//		TVector3 vb = pb + sSt*db;
//		v = (va+vb);
//		v *= 0.5;
//
//		//std::cout << TMath::Abs(r * da.Cross(db).Mag()) << "\t" << (va-vb).Mag() << std::endl;
//		//va.Print();
//		//vb.Print();
//		//v.Print();
//		//(va-vb).Print();
//
//		//shortest distance is length of r(AxB)
//		return (va-vb).Mag();//TMath::Abs(r * da.Cross(db).Mag());
//	}
//	return 0.;
//}

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
//This method comes from a Google search but I've lost the link so use the third instead (this one gives matching results)
//double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
//		const TVector3 &pb, const TVector3 &db) {
//	double x01 = pa.X();
//	double x02 = pb.X();
//	double y01 = pa.Y();
//	double y02 = pb.Y();
//	double z01 = pa.Z();
//	double z02 = pb.Z();
//
//	double l1 = da.X();
//	double l2 = db.X();
//	double m1 = da.Y();
//	double m2 = db.Y();
//	double n1 = da.Z();
//	double n2 = db.Z();
//
//
//
//	double t1 = (-l1*(pow(m2,2) + pow(n2,2))*(x01 - x02) - (m2*n1 - m1*n2)*(n2*(-y01 + y02) + 
//				m2*(z01 - z02)) + pow(l2,2)*(m1*(-y01 + y02) + n1*(-z01 + z02)) + 
//			l2*(m1*m2*(x01 - x02) + n1*n2*(x01 - x02) + l1*(m2*y01 - m2*y02 + n2*z01 - n2*z02))) / 
//		(pow(l2,2)*(pow(m1,2) + pow(n1,2)) + pow(m2*n1 - m1*n2, 2) - 2*l1*l2*(m1*m2 + n1*n2) + pow(l1,2) * (pow(m2,2) + pow(n2,2)));
//
//	double t2 = -((l1*l2 + m1*m2 + n1*n2)*(l1*(x01 - x02) + m1*(y01 - y02) + 
//				n1*(z01 - z02)) - (pow(l1,2) + pow(m1,2) + pow(n1,2))*(l2*(x01 - x02) + 
//				m2*(y01 - y02) + n2*(z01 - z02))) / 
//		(-pow(l1*l2 + m1*m2 + n1*n2, 2) + (pow(l1,2) + pow(m1,2) + pow(n1,2))*(pow(l2,2) + pow(m2,2) + pow(n2,2)));
//
//	TVector3 p(x01 + l1 * t1,
//			y01 + m1 * t1,
//			z01 + n1 * t1);
//
//	TVector3 q(x02 + l2 * t2,
//			y02 + m2 * t2,
//			z02 + n2 * t2);
//
//	v = 0.5*(p+q);
//
//	//pa.Print();
//	//pb.Print();
//	//da.Print();
//	//db.Print();
//	//p.Print();
//	//q.Print();
//	//v.Print();
//	//(p-q).Print();
//
//	return (p-q).Mag();
//}

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

//main loop
void resampleTrackIPs::Loop()
{
	if (fChain == 0) return;

	//setup histograms - to plot these run with _irep=-1 (runs generated vs resampled comparison)
	//IP
	TH1D hist0("hist0", "",   10000, 0, 10);
	TH1D hist1("hist1", "",   10000, 0, 10);
	TH1D hist1b("hist1b", "", 10000, 0, 10);
	TH1D hist2("hist2", "",   10000, 0, 10);
	//IP : sigma^2
	TH2D hist3("hist3", "", 10000, 0, 100, 3, 0, 1e-3);
	TH2D hist4("hist4", "", 10000, 0, 100, 3, 0, 1e-3);
	//IPchi2
	TH1D hist5("hist5", "", 10000, 0, 30);
	TH1D hist6("hist6", "", 10000, 0, 30);
	TH1D hist7("hist7", "", 10000, 0, 30);
	//IP (chi2<16)
	TH1D hist8("hist8", "", 100, 0, 1);
	TH1D hist9("hist9", "", 100, 0, 1);
	TH1D hist10("hist10", "", 100, 0, 1);
	//IP (chi2>16)
	TH1D hist11("hist11", "", 100, 0, 1);
	TH1D hist12("hist12", "", 100, 0, 1);
	TH1D hist13("hist13", "", 100, 0, 1);
	//IPx
	TH1D hist14("hist14", "", 100, -1., 1.);
	TH1D hist15("hist15", "", 100, -1., 1.);
	//IPy
	TH1D hist16("hist16", "", 100, -1., 1.);
	TH1D hist17("hist17", "", 100, -1., 1.);

	//IP over sigma
	double binsx[59] = {0.,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,6.,8.,10.,20.,30.,75.,150.};
	double binsy[4] = {0.000,0.0005,0.001,0.002};
	TH2D hist18("hist18","", 58, binsx, 3, binsy);
	TH2D hist19("hist19","", 58, binsx, 3, binsy);
	TH2D hist20("hist20","", 58, binsx, 3, binsy);

	//shift in PV
	TH1D hist21("hist21", "", 1000, -0.1, 0.1);
	TH1D hist22("hist22", "", 1000, -0.1, 0.1);
	TH1D hist23("hist23", "", 1000, -0.1, 0.1);

	//track separation x,y and z
	TH1D vhist0a("vhist0a","",20,0,5.);//40);
	TH1D vhist0b("vhist0b","",20,0,5.);//40);
	TH1D vhist1a("vhist1a","",20,0,5.);//40);
	TH1D vhist1b("vhist1b","",20,0,5.);//40);

	TH1D vhist0c("vhist0c","",20,0,5.);//40);
	TH1D vhist0d("vhist0d","",20,0,5.);//40);
	TH1D vhist1c("vhist1c","",20,0,5.);//40);
	TH1D vhist1d("vhist1d","",20,0,5.);//40);
	//TH1D vhist2c("vhist2c","",40,0,20.);//200);
	//TH1D vhist2d("vhist2d","",40,0,20.);//200);

	//Ntracks
	TH1D vhist3a("vhist3a","",50,0,50);
	TH1D vhist3b("vhist3b","",50,0,50);

	//tz
	TH1D vhist4a("vhist4a","",50,-1,1);
	TH1D vhist4b("vhist4b","",50,-1,1);
	//fvz
	TH1D vhist5z("vhist5z","",50,-10,10);
	TH1D vhist5a("vhist5a","",50,-10,10);
	TH1D vhist5b("vhist5b","",50,-10,10);
	//m
	TH1D vhist6a("vhist6a","",50,0,5000);
	TH1D vhist6b("vhist6b","",50,0,5000);
	//pz
	TH1D vhist7a("vhist7a","",50,0,100000);
	TH1D vhist7b("vhist7b","",50,0,100000);
	//mcor
	TH1D vhist8a("vhist8a","",50,0,5000);
	TH1D vhist8b("vhist8b","",50,0,5000);
	//doca
	TH1D vhist9a("vhist9a","",100,0,2.);
	TH1D vhist9b("vhist9b","",100,0,2.);

	//N svrs
	TH1D vhist10a("vhist10a","",20,0,20);
	TH1D vhist10b("vhist10b","",20,0,20);

	//post-selection histograms
	TH1D vhist4c("vhist4c","",50,-1,1);
	TH1D vhist4d("vhist4d","",50,-1,1);
	TH1D vhist5c("vhist5c","",50,-10,10);
	TH1D vhist5d("vhist5d","",50,-10,10);
	TH1D vhist6c("vhist6c","",50,0,5000);
	TH1D vhist6d("vhist6d","",50,0,5000);
	TH1D vhist7c("vhist7c","",50,0,100000);
	TH1D vhist7d("vhist7d","",50,0,100000);
	TH1D vhist8c("vhist8c","",50,0,5000);
	TH1D vhist8d("vhist8d","",50,0,5000);
	TH1D vhist9c("vhist9c","",100,0,2.);
	TH1D vhist9d("vhist9d","",100,0,2.);

	//load the IP/sigma map for each 1/pT bin
	TFile* fh = TFile::Open("ipmap.root");
	TH2D* h = dynamic_cast<TH2D*>(fh->Get("ipr"));
	std::vector<TH1D*> hs;
	for(int i=1; i<=h->GetNbinsY(); ++i) {
		TString name = "_bin"; name+=i;
		hs.push_back(h->ProjectionX(name,i,i));
	}

	//keep track of various stats
	int countBefore(0), countAfter(0), total(0), nonPrompt(0);
	int svsBefore(0), svsAfter1(0), svsAfter2(0), svNsAfter1(0), svNsAfter2(0), sv2sAfter1(0), sv2sAfter2(0), sv3sAfter1(0), sv3sAfter2(0), sv4sAfter1(0), sv4sAfter2(0), svMoresAfter1(0), svMoresAfter2(0);

	int seed(0);
	if(_irep>=0) seed=1000+_irep;
	TRandom3 r(seed);

	int sumntrk1(0), sumncomb1(0);
	int sumntrk2(0), sumncomb2(0);

	int pvr_idx(0), binInvPt(0);
	double tx(0.), ty(0.), dx(0.), dy(0.), newDx(0.), newDy(0.), invpt(0.), sigma(0.), sigma2(0.);

	int firstrep(0), endrep(101);
	if(_irep>=0) {
		firstrep=_irep;
		endrep=_irep+1;
	}

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		if(_irep<0 && jentry%1000==0) std::cout << jentry << " of " << nentries << std::endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		svsBefore += svr_z->size();
		for(uint isvr=0; isvr<svr_z->size(); ++isvr) {
			pvr_idx = svr_idx_pvr->at(isvr);

			vhist5z.Fill(svr_z->at(isvr) - pvr_z->at(pvr_idx));
		}

		std::vector<double>* origtrk_x = new std::vector<double>();
		std::vector<double>* origtrk_y = new std::vector<double>();
		for(uint i=0; i<trk_x->size(); ++i) {
			origtrk_x->push_back(trk_x->at(i));
			origtrk_y->push_back(trk_y->at(i));
		}

		//resample each event as required
		//loop zero is unaltered
		//when filling TTrees we only run this loop once to fill a single TTree
		for(int irep=firstrep; irep<endrep; ++irep) {
			//clear out the old SVR variables to be overwritten
			svr_idx_pvr->clear();
			svr_idx_jet->clear();
			svr_idx_trk0->clear();
			svr_idx_trk1->clear();
			svr_idx_trk2->clear();
			svr_idx_trk3->clear();
			svr_idx_trk4->clear();
			svr_idx_trk5->clear();
			svr_idx_trk6->clear();
			svr_idx_trk7->clear();
			svr_idx_trk8->clear();
			svr_idx_trk9->clear();
			svr_px->clear();
			svr_py->clear();
			svr_pz->clear();
			svr_e->clear();
			svr_x->clear();
			svr_y->clear();
			svr_z->clear();
			svr_fd_min->clear();
			svr_m_cor->clear();

			std::vector<int> trks;
			std::vector<TVector3> trkhit, trkdir;

			for(uint idx=0; idx<trk_idx_gen->size(); ++idx) {
				int gen_idx = trk_idx_gen->at(idx);
				if(gen_idx<0) continue;
				if(trk_idx_pvr->at(idx)<0) continue;//TODO
				pvr_idx = gen_idx_pvr->at(gen_idx);
				if(pvr_idx<0) continue;
				if(trk_vid->at(idx) == 1) continue; //ignore upstream tracks TODO can we ignore VELO tracks here too?
				if(trk_pt->at(idx)<500 || trk_prb_ghost->at(idx)>0.3) continue; //check track meets requirements that won't change (Phys/JetTagging/src/LoKiBDTTag.cpp:262)
				if(TMath::Abs(gen_x->at(gen_idx) - pvr_x->at(pvr_idx)) > 1e-5 ||
						TMath::Abs(gen_y->at(gen_idx) - pvr_y->at(pvr_idx)) > 1e-5 ||
						TMath::Abs(gen_z->at(gen_idx) - pvr_z->at(pvr_idx)) > 1e-5) {++nonPrompt; continue;}//only use prompt tracks
				//pvr_idx = trk_idx_pvr->at(idx);//TODO
				//if(pvr_idx<0) continue;//TODO
				//TODO track type not saved offline

				//resample IP of track
				//if(pvr_idx != trk_idx_pvr->at(idx)) std::cout << "!!" << std::endl;
				//int pvr_idx = trk_idx_pvr->at(idx);

				//invpt = 1./TMath::Sqrt(gen_px->at(gen_idx)*gen_px->at(gen_idx) + gen_py->at(gen_idx)*gen_py->at(gen_idx));
				invpt = 1./TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx));

				tx = trk_px->at(idx) / trk_pz->at(idx);
				ty = trk_py->at(idx) / trk_pz->at(idx);

				dx = origtrk_x->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*tx - pvr_x->at(pvr_idx);
				dy = origtrk_y->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*ty - pvr_y->at(pvr_idx);
				//TVector3 tmp = calcDxDy(TVector3(origtrk_x->at(idx),origtrk_y->at(idx),trk_z->at(idx)),
				//		        TVector3(trk_px->at(idx),trk_py->at(idx),trk_pz->at(idx)),
				//		        TVector3(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx)));
				//tmp.Print();
				//std::cout << dx << "\t" << dy << std::endl;

				binInvPt = h->GetYaxis()->FindBin(invpt);
				sigma2 = trk_ip->at(idx)*trk_ip->at(idx) / trk_ip_chi2->at(idx); //assume unchanged by shift
				sigma = TMath::Sqrt(sigma2);

				++total;
				//      	        if(trk_ip_chi2->at(idx) >= 16.) ++countBefore;
				//std::cout << gen_idx_pvr->at(gen_idx) << "\t" << trk_idx_pvr->at(idx) << std::endl;//TODO
				//std::cout << pvr_x->at(gen_idx_pvr->at(gen_idx)) << "\t" << pvr_x->at(trk_idx_pvr->at(idx)) << std::endl;//TODO
				//std::cout << pvr_y->at(gen_idx_pvr->at(gen_idx)) << "\t" << pvr_y->at(trk_idx_pvr->at(idx)) << std::endl;//TODO
				//std::cout << pvr_z->at(gen_idx_pvr->at(gen_idx)) << "\t" << pvr_z->at(trk_idx_pvr->at(idx)) << std::endl;//TODO

				newDx = dx;
				newDy = dy;

				if(irep>0) {
					//newDx = hs[binInvPt-1]->GetRandom();
					//newDy = hs[binInvPt-1]->GetRandom();
					//if(r.Integer(2)==0) newDx=-newDx;
					//if(r.Integer(2)==0) newDy=-newDy;

					double phi=r.Rndm()*2*TMath::Pi();
					//std::cout << trk_ip->at(idx) << "\t" << trk_ip_chi2->at(idx) << "\t" << sigma2 << "\t" << sigma << std::endl;
					double rho = sigma*getRandom(hs[binInvPt-1],&r); //interpolate rho within the bin
					//double rho=sigma*hs[binInvPt-1]->GetRandom();
					newDx = rho*TMath::Sin(phi);
					newDy = rho*TMath::Cos(phi);

					trk_x->at(idx) = origtrk_x->at(idx) + newDx - dx;
					trk_y->at(idx) = origtrk_y->at(idx) + newDy - dy;
				}

				//		std::cout << dx << "\t" << newDx << std::endl;

				//h1.Fill(TMath::Sqrt(dx*dx + dy*dy),invpt);
				//h2.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy),invpt);

				//now check track meets IP chi2 requirement (Tab 2 ANA-2014-074)
				//IP
				if(irep==0) {
					hist0.Fill(trk_ip->at(idx));
					hist1.Fill(TMath::Sqrt(dx*dx + dy*dy));
					double tmpX = origtrk_x->at(idx) + (pvr_z->at(trk_idx_pvr->at(idx)) - trk_z->at(idx))*tx - pvr_x->at(trk_idx_pvr->at(idx));
					double tmpY = origtrk_y->at(idx) + (pvr_z->at(trk_idx_pvr->at(idx)) - trk_z->at(idx))*ty - pvr_y->at(trk_idx_pvr->at(idx));
					hist1b.Fill(TMath::Sqrt(tmpX*tmpX + tmpY*tmpY));
				} else {
					hist2.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy));
				}
				//with chi2 cuts
				if(irep==0) {
					if(trk_ip_chi2->at(idx) < 16) hist8.Fill(trk_ip->at(idx));
					else hist11.Fill(trk_ip->at(idx));
					if((dx*dx + dy*dy)/sigma2 < 16) hist9.Fill(TMath::Sqrt(dx*dx + dy*dy));
					else hist12.Fill(TMath::Sqrt(dx*dx + dy*dy));
				} else {
					if((newDx*newDx + newDy*newDy)/sigma2 < 16) hist10.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy));
					else hist13.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy));
				}
				if(irep==0) {
					//IPx
					hist14.Fill(dx);
					//IPy
					hist16.Fill(dy);
					//IP vs resolution^2
					hist3.Fill(TMath::Sqrt(dx*dx + dy*dy), sigma2);
				} else {
					//IPx
					hist15.Fill(newDx);
					//IPy
					hist17.Fill(newDy);
					//IP vs resolution^2
					hist4.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy), sigma2);
				}
				//IP/sigma
				if(irep==0) {
					hist18.Fill(TMath::Sqrt(trk_ip_chi2->at(idx)),invpt);
					hist19.Fill(TMath::Sqrt(dx*dx + dy*dy)/sigma,invpt);

					//std::cout << TMath::Sqrt(trk_ip->at(idx)) << "\t" << TMath::Sqrt(dx*dx + dy*dy) << std::endl;
				} else {
					hist20.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy)/sigma,invpt);
				}
				//IPchi2
				if(irep==0) {
					hist5.Fill(trk_ip_chi2->at(idx));
					hist6.Fill((dx*dx + dy*dy)/sigma2);
				} else {
					hist7.Fill((newDx*newDx + newDy*newDy)/sigma2);
				}
				if(irep>0) {
					hist21.Fill(pvr_x->at(trk_idx_pvr->at(idx)) - pvr_x->at(pvr_idx));
					hist22.Fill(pvr_y->at(trk_idx_pvr->at(idx)) - pvr_y->at(pvr_idx));
					hist23.Fill(pvr_z->at(trk_idx_pvr->at(idx)) - pvr_z->at(pvr_idx));
				}

				//apply the IPchi2 cut to the new IPchi2
				if((newDx*newDx + newDy*newDy) / sigma2 < 16.) continue;
				if(irep>0) ++countAfter;
				else ++countBefore;

				//get direction and point on track for vertexing
				TVector3 dir(trk_px->at(idx),trk_py->at(idx),trk_pz->at(idx));
				dir.SetMag(1.);
				TVector3 hit(origtrk_x->at(idx) + newDx - dx,
						origtrk_y->at(idx) + newDy - dy,
						trk_z->at(idx));

				double stepX = 10*TMath::Tan(dir.Theta())*TMath::Cos(dir.Phi());
				double stepY = 10*TMath::Tan(dir.Theta())*TMath::Sin(dir.Phi());

				while((hit.X()*hit.X()+hit.Y()*hit.Y()) < 6*6){
					hit.SetZ(hit.Z()+10);
					hit.SetX(hit.X()+stepX);
					hit.SetY(hit.Y()+stepY);
				}

				trks.push_back(idx);
				trkhit.push_back(hit);
				trkdir.push_back(dir);
			}

			int ntrk = trks.size();
			if(irep==0) vhist3a.Fill(ntrk);
			else vhist3b.Fill(ntrk);
			if(irep==0) {
				sumntrk1+=ntrk;
				sumncomb1+=ntrk*(ntrk-1)/2.;
			} else {
				sumntrk2+=ntrk;
				sumncomb2+=ntrk*(ntrk-1)/2.;
			}

			//check we have at least two tracks to vertex
			if(ntrk>1) {
				//now make 2body SVs
				std::vector<std::pair<int,int> > sv2ij;
				std::vector<TVector3> sv2;

				for( int itrk=0; itrk<ntrk; ++itrk) {
					TLorentzVector p4i(trk_px->at(trks[itrk]), trk_py->at(trks[itrk]), trk_pz->at(trks[itrk]), 0);
					p4i.SetE(TMath::Sqrt(p4i.P()*p4i.P() + 140*140));
					for( int jtrk=itrk+1; jtrk<ntrk; ++jtrk) {
						if(gen_idx_pvr->at(trk_idx_gen->at(trks[itrk])) != gen_idx_pvr->at(trk_idx_gen->at(trks[jtrk]))) continue;

						TVector3 dxdyi = calcDxDy(trkhit[itrk], trkdir[itrk], TVector3(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx)));
						TVector3 dxdyj = calcDxDy(trkhit[jtrk], trkdir[jtrk], TVector3(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx)));
						if(irep==0) {
							vhist0a.Fill(dxdyi.X() - dxdyj.X());
							vhist1a.Fill(dxdyi.Y() - dxdyj.Y());
						} else {
							vhist0b.Fill(dxdyi.X() - dxdyj.X());
							vhist1b.Fill(dxdyi.Y() - dxdyj.Y());
						}
						//if(TMath::Abs(trkhit[itrk].X() - trkhit[jtrk].X())<0.01 && TMath::Abs(trkhit[itrk].Y() - trkhit[jtrk].Y())<0.01 && TMath::Abs(trkhit[itrk].Z() - trkhit[jtrk].Z())==0.) {
						//	std::cout << trk_p->at(trks[itrk]) << "\t" << trk_p->at(trks[jtrk]) << "\t" << trk_pid->at(trks[itrk]) << "\t" << trk_pid->at(trks[jtrk]) << std::endl;
						//	std::cout << origtrk_x->at(trks[itrk]) << "\t" << origtrk_x->at(trks[jtrk]) << "\t" << origtrk_y->at(trks[itrk]) << "\t" << origtrk_y->at(trks[jtrk]) << std::endl;

						//	trkhit[itrk].Print();
						//	trkdir[itrk].Print();
						//	trkhit[jtrk].Print();
						//	trkdir[jtrk].Print();
						//}
						TLorentzVector p4j(trk_px->at(trks[jtrk]), trk_py->at(trks[jtrk]), trk_pz->at(trks[jtrk]), 0);
						p4j.SetE(TMath::Sqrt(p4j.P()*p4j.P() + 140*140));
						TLorentzVector p4 = p4i+p4j;
						double m = p4.M();
						TVector3 sv;
						double d = calcDoca(sv, trkhit[itrk], trkdir[itrk], trkhit[jtrk], trkdir[jtrk]);
						pvr_idx = gen_idx_pvr->at(trk_idx_gen->at(trks[itrk]));
						//pvr_idx = trk_idx_pvr->at(trks[itrk]);
						TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
								sv.Y() - pvr_y->at(pvr_idx),
								sv.Z() - pvr_z->at(pvr_idx));
						double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();
						double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);
						double tz = fv.Z()*m/p4.Z()/(3e11)*(1e12);
						if(irep==0) {
							vhist4a.Fill(tz);
							vhist5a.Fill(fv.Z());
							vhist6a.Fill(m);
							vhist7a.Fill(p4.Z());
							vhist8a.Fill(mcor);
							vhist9a.Fill(d);
						} else {
							vhist4b.Fill(tz);
							vhist5b.Fill(fv.Z());
							vhist6b.Fill(m);
							vhist7b.Fill(p4.Z());
							vhist8b.Fill(mcor);
							vhist9b.Fill(d);
						}
						//			std::cout << fv.Z() << "\t" << m << "\t" << p4.Z() << "\t" << tz << std::endl;
						//Apply 2-body SV requirements
						if(m < 0. || m > 5000.) continue;
						if(d > 0.2) continue;
						if(mcor<600.) continue;
						if(tz<0. || tz>=10.) continue;
						//TODO add VTXCHI2<10, FDCHI2>25, HITS>0 
						sv2ij.push_back(std::pair<int,int>(trks[itrk],trks[jtrk]));
						sv2.push_back(sv);
						if(irep==0) {
							vhist0c.Fill(dxdyi.X() - dxdyj.X());//TMath::Abs(trkhit[itrk].X() - trkhit[jtrk].X()));
							vhist1c.Fill(dxdyi.Y() - dxdyj.Y());//TMath::Abs(trkhit[itrk].Y() - trkhit[jtrk].Y()));
							//vhist2c.Fill(TMath::Abs(trkhit[itrk].Z() - trkhit[jtrk].Z()));
							vhist4c.Fill(tz);
							vhist5c.Fill(fv.Z());
							vhist6c.Fill(m);
							vhist7c.Fill(p4.Z());
							vhist8c.Fill(mcor);
							vhist9c.Fill(d);
						} else {
							vhist0d.Fill(dxdyi.X() - dxdyj.X());//TMath::Abs(trkhit[itrk].X() - trkhit[jtrk].X()));
							vhist1d.Fill(dxdyi.Y() - dxdyj.Y());//TMath::Abs(trkhit[itrk].Y() - trkhit[jtrk].Y()));
							//vhist2d.Fill(TMath::Abs(trkhit[itrk].Z() - trkhit[jtrk].Z()));
							vhist4d.Fill(tz);
							vhist5d.Fill(fv.Z());
							vhist6d.Fill(m);
							vhist7d.Fill(p4.Z());
							vhist8d.Fill(mcor);
							vhist9d.Fill(d);
						}
					}
				}


				if(irep==0) {
					svsAfter1+=sv2.size();
					vhist10a.Fill(sv2.size());
				} else {
					svsAfter2+=sv2.size();
					vhist10b.Fill(sv2.size());
				}

				std::vector<TVector3> svN; //vertices of linked SVs
				std::vector<std::vector<int> > svNij; //tracks associated with each linked SV
				std::set<int> used; // only use each SV once

				//check we have at least one 2body SV
				if(sv2.size()==1) {
					//only one SV - no need to link
					svN.push_back(sv2[1]);
					if(irep==0) ++sv2sAfter1;
					else ++sv2sAfter2;
				} else if(sv2.size()>1) {
					//more than one SV - now link
					for(uint s=0; s<sv2.size(); ++s) {
						if(used.count(s)>0) continue;
						used.insert(s);
						std::vector<int> svTrks;
						TVector3 sv = sv2[s];
						svTrks.push_back(sv2ij[s].first);
						svTrks.push_back(sv2ij[s].second);
						int nsvCombined(1);
						int jet_idx(-1);
						pvr_idx = gen_idx_pvr->at(trk_idx_gen->at(svTrks[0]));

						double fd_min = TMath::Sqrt(TMath::Power(sv2[s].X() - pvr_x->at(pvr_idx),2) +
							                    TMath::Power(sv2[s].Y() - pvr_y->at(pvr_idx),2));

						//iterate through each track we add and check for additional SVs
						for(uint t=0; t<svTrks.size(); ++t) {
							if(jet_idx==-1 && trk_idx_jet->at(svTrks[t])!=-1) jet_idx=trk_idx_jet->at(svTrks[t]);

							//loop over remaining SVs
							for(uint ss=s+1; ss<sv2.size(); ++ss) {
								if(used.count(ss)>0) continue;
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
									sv += sv2[ss];
									++nsvCombined;
									double fd = TMath::Sqrt(TMath::Power(sv2[s].X() - pvr_x->at(pvr_idx),2) +
											        TMath::Power(sv2[s].Y() - pvr_y->at(pvr_idx),2));
									if(fd<fd_min) fd_min = fd;
								}
							}
						}
						//sort tracks and remove repeats
						std::sort(svTrks.begin(), svTrks.end());
						svTrks.erase(std::unique(svTrks.begin(),svTrks.end()),svTrks.end());
						//normalise average vertex //TODO weight?
						sv *= (1./nsvCombined);
						//save vertex and list of tracks
						svN.push_back(sv);
						svNij.push_back(svTrks);
						//count SVs
						if(svTrks.size()==2) {
							if(irep==0) ++sv2sAfter1;
							else ++sv2sAfter2;
						} else if(svTrks.size()==3) {
							if(irep==0) ++sv3sAfter1;
							else ++sv3sAfter2;
						} else if(svTrks.size()==4) {
							if(irep==0) ++sv4sAfter1;
							else ++sv4sAfter2;
						} else {
							if(irep==0) ++svMoresAfter1;
							else ++svMoresAfter2;
						}

						TLorentzVector p4;

						for(uint itrk=0; itrk<svTrks.size(); ++itrk) {
							TLorentzVector p4i(trk_px->at(svTrks[itrk]), trk_py->at(svTrks[itrk]), trk_pz->at(svTrks[itrk]), 0);
							p4i.SetE(TMath::Sqrt(p4i.P()*p4i.P() + 140*140));
							p4 += p4i;
						}

						double m = p4.M();
						//pvr_idx = trk_idx_pvr->at(trks[itrk]);
						TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
							    sv.Y() - pvr_y->at(pvr_idx),
							    sv.Z() - pvr_z->at(pvr_idx));
						double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();
						double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);

						//pad with -1 for filling tuple - size of svTrks no longer meaningful
						while(svTrks.size()<10) {
							svTrks.push_back(-1);
						}

						svr_idx_pvr->push_back(trk_idx_pvr->at(svTrks[0]));
						svr_idx_jet->push_back(jet_idx);
						svr_idx_trk0->push_back(svTrks[0]);
						svr_idx_trk1->push_back(svTrks[1]);
						svr_idx_trk2->push_back(svTrks[2]);
						svr_idx_trk3->push_back(svTrks[3]);
						svr_idx_trk4->push_back(svTrks[4]);
						svr_idx_trk5->push_back(svTrks[5]);
						svr_idx_trk6->push_back(svTrks[6]);
						svr_idx_trk7->push_back(svTrks[7]);
						svr_idx_trk8->push_back(svTrks[8]);
						svr_idx_trk9->push_back(svTrks[9]);
						svr_px->push_back(p4.Px());
						svr_py->push_back(p4.Py());
						svr_pz->push_back(p4.Pz());
						svr_e->push_back(p4.E());
						svr_x->push_back(sv.X());
						svr_y->push_back(sv.Y());
						svr_z->push_back(sv.Z());
						svr_fd_min->push_back(fd_min);
						svr_m_cor->push_back(mcor);

						//std::cout << svTrks.size() << "-track vertex found ( ";
						//for(int i=0; i<svTrks.size(); ++i) std::cout << svTrks[i] << " ";
						//std::cout << ")" << std::endl;
						//if(used.size()!=sv2.size()) {
						//	std::cout << "linked " << used.size() << " of " << sv2.size() << " vertices" << std::endl;
						//	std::cout << svTrks.size() << "-track vertex found ( ";
						//	for(uint i=0; i<svTrks.size(); ++i) std::cout << svTrks[i] << " ";
						//	std::cout << ") from " << trks.size() << " tracks" << std::endl;
						//}
					}
				}

				if(irep==0) svNsAfter1+=svN.size();
				else svNsAfter2+=svN.size();
			}

			if(newtree) newtree->Fill();
		}

		origtrk_x->clear();
		origtrk_y->clear();

		delete origtrk_x;
		delete origtrk_y;

	}

	if(newtree) {
		newtree->AutoSave();
		newfile->Close();
	}

	if(_irep>=0) {
		std::cout << countBefore+countAfter << " qualifying tracks used and " << svNsAfter1+svNsAfter2 << " SVs generated..." << std::endl;
		return;
	}
	//below here it's just stats and histograms

	std::cout << "Ntrk/Ntrkin\t" << static_cast<double>(countBefore)/total << " +/- " << TMath::Sqrt(countBefore)/total << "\t" << static_cast<double>(countAfter)/total/100. << " +/- " << TMath::Sqrt(countAfter)/total/100. << "\t" << total << std::endl;
	std::cout << "Nisv\t" << "---" << "\t" << static_cast<double>(svsAfter1)/1. << "\t" << static_cast<double>(svsAfter2)/100. << std::endl;
	std::cout << "Nsv2\t" << "---" << "\t" << static_cast<double>(sv2sAfter1)/1. << "\t" << static_cast<double>(sv2sAfter2)/100. << std::endl;
	std::cout << "Nsv3\t" << "---" << "\t" << static_cast<double>(sv3sAfter1)/1. << "\t" << static_cast<double>(sv3sAfter2)/100. << std::endl;
	std::cout << "Nsv4\t" << "---" << "\t" << static_cast<double>(sv4sAfter1)/1. << "\t" << static_cast<double>(sv4sAfter2)/100. << std::endl;
	std::cout << "Nsv5+\t" << "---" << "\t" << static_cast<double>(svMoresAfter1)/1. << "\t" << static_cast<double>(svMoresAfter2)/100. << std::endl;
	std::cout << "NsvN\t" << svsBefore << "\t" << static_cast<double>(svNsAfter1)/1. << "\t" << static_cast<double>(svNsAfter2)/100. << std::endl;
	std::cout << "Ntrk\t" << sumntrk1 << "\t" << static_cast<double>(sumntrk2)/100. << std::endl;
	std::cout << "Ncomb\t" << sumncomb1 << "\t" << static_cast<double>(sumncomb2)/100. << std::endl;
	std::cout << "Nisv/Ncomb\t" << static_cast<double>(svsAfter1)/sumncomb1 << " +/- " << TMath::Sqrt(svsAfter1)/sumncomb1 << "\t" << static_cast<double>(svsAfter2)/sumncomb2 << " +/- " << TMath::Sqrt(svsAfter2)/sumncomb2 << std::endl;
	std::cout << "prompt:" << static_cast<double>(total)/(total+nonPrompt) << std::endl;

	gStyle->SetOptStat(0);
	TH1D* htmp(0);

	hist0.Scale(1.   );//,"width");
	hist1.Scale(1    );//,"width");
	hist1b.Scale(1    );//,"width");
	hist2.Scale(0.01 );//,"width");
	hist3.Scale(1    );//,"width");
	hist4.Scale(0.01 );//,"width");
	hist5.Scale(1.   );//,"width");
	hist6.Scale(1    );//,"width");
	hist7.Scale(0.01 );//,"width");
	hist8.Scale(1.   );//,"width");
	hist9.Scale(1    );//,"width");
	hist10.Scale(0.01);//,"width");
	hist11.Scale(1.  );//,"width");
	hist12.Scale(1   );//,"width");
	hist13.Scale(0.01);//,"width");
	hist14.Scale(1   );//,"width");
	hist15.Scale(0.01);//,"width");
	hist16.Scale(1   );//,"width");
	hist17.Scale(0.01);//,"width");
	hist18.Scale(1.  );//,"width");
	hist19.Scale(1   );//,"width");
	hist20.Scale(0.01);//,"width");

	hist0.GetXaxis()->SetTitle("IP");
	hist3.GetXaxis()->SetTitle("IP");
	hist5.GetXaxis()->SetTitle("#chi^{2}_{IP}");
	hist8.GetXaxis()->SetTitle("IP (#chi^{2}_{IP} < 16)");
	hist11.GetXaxis()->SetTitle("IP (#chi^{2}_{IP} > 16)");
	hist14.GetXaxis()->SetTitle("IP_{X}");
	hist16.GetXaxis()->SetTitle("IP_{Y}");

	TCanvas c;
	c.SetLogx(1);
	c.SetLogy(0);
	if(hist1.GetMaximum()>hist0.GetMaximum()) hist0.SetMaximum(1.1*hist1.GetMaximum());
	if(hist2.GetMaximum()>hist0.GetMaximum()) hist0.SetMaximum(1.1*hist2.GetMaximum());
	hist0.Draw("EP");
	hist1.SetLineColor(kRed);
	hist1.Draw("same");
	hist1b.SetLineColor(kGreen);
	hist1b.Draw("same");
	hist2.SetLineColor(kMagenta);
	hist2.Draw("same");
	c.SaveAs("ip.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	if(hist9.GetMaximum()>hist8.GetMaximum()) hist8.SetMaximum(1.1*hist9.GetMaximum());
	if(hist10.GetMaximum()>hist8.GetMaximum()) hist8.SetMaximum(1.1*hist10.GetMaximum());
	hist8.Draw("EP");
	hist9.SetLineColor(kRed);
	hist9.Draw("same");
	hist10.SetLineColor(kMagenta);
	hist10.Draw("same");
	c.SaveAs("ip_ipchi2LT16.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	if(hist12.GetMaximum()>hist11.GetMaximum()) hist11.SetMaximum(1.1*hist12.GetMaximum());
	if(hist13.GetMaximum()>hist11.GetMaximum()) hist11.SetMaximum(1.1*hist13.GetMaximum());
	hist11.Draw("EP");
	hist12.SetLineColor(kRed);
	hist12.Draw("same");
	hist13.SetLineColor(kMagenta);
	hist13.Draw("same");
	c.SaveAs("ip_ipchi2GT16.pdf");

	c.SetLogx(0);
	c.SetLogy(0);
	if(hist15.GetMaximum()>hist14.GetMaximum()) hist14.SetMaximum(1.1*hist15.GetMaximum());
	hist14.SetLineColor(kRed);
	hist14.Draw();
	hist15.SetLineColor(kMagenta);
	hist15.Draw("same");
	c.SaveAs("ipx.pdf");

	c.SetLogx(0);
	c.SetLogy(0);
	if(hist17.GetMaximum()>hist16.GetMaximum()) hist16.SetMaximum(1.1*hist17.GetMaximum());
	hist16.SetLineColor(kRed);
	hist16.Draw();
	hist17.SetLineColor(kMagenta);
	hist17.Draw("same");
	c.SaveAs("ipy.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	if(hist6.GetMaximum()>hist5.GetMaximum()) hist5.SetMaximum(1.1*hist6.GetMaximum());
	if(hist7.GetMaximum()>hist5.GetMaximum()) hist5.SetMaximum(1.1*hist7.GetMaximum());
	hist5.Draw("EP");
	hist6.SetLineColor(kRed);
	hist6.Draw("same");
	hist7.SetLineColor(kMagenta);
	hist7.Draw("same");
	c.SaveAs("ipchi2.pdf");

	c.SetLogx(1);
	c.SetLogy(1);
	htmp = hist3.ProjectionX("hist1_1",1,1);
	htmp->SetLineColor(kBlue); htmp->Draw("");
	htmp = hist4.ProjectionX("hist2_1",1,1);
	htmp->SetLineColor(kBlue); htmp->SetLineStyle(kDashed); htmp->Draw("same");
	htmp = hist3.ProjectionX("hist1_2",2,2);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist4.ProjectionX("hist2_2",2,2);
	htmp->SetLineColor(kRed); htmp->SetLineStyle(kDashed); htmp->Draw("same");
	htmp = hist3.ProjectionX("hist1_3",3,3);
	htmp->SetLineColor(kMagenta); htmp->Draw("same");
	htmp = hist4.ProjectionX("hist2_3",3,3);
	htmp->SetLineColor(kMagenta); htmp->SetLineStyle(kDashed); htmp->Draw("same");
	c.SaveAs("ip_sigma2.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	htmp = hist18.ProjectionX("hist18_1",1,1);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_1",1,1);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_1",1,1);
	htmp->SetLineColor(kMagenta); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin1.pdf");
	htmp = hist18.ProjectionX("hist18_2",2,2);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_2",2,2);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_2",2,2);
	htmp->SetLineColor(kMagenta); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin2.pdf");
	htmp = hist18.ProjectionX("hist18_3",3,3);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_3",3,3);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_3",3,3);
	htmp->SetLineColor(kMagenta); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin3.pdf");

	vhist0c.Scale(1./vhist0a.Integral());
	vhist1c.Scale(1./vhist1a.Integral());
	vhist0d.Scale(1./vhist0b.Integral());
	vhist1d.Scale(1./vhist1b.Integral());
	vhist4c.Scale(1./vhist4a.Integral());
	vhist5c.Scale(1./vhist5a.Integral());
	vhist4d.Scale(1./vhist4b.Integral());
	vhist5d.Scale(1./vhist5b.Integral());
	vhist6c.Scale(1./vhist6a.Integral());
	vhist7c.Scale(1./vhist7a.Integral());
	vhist8c.Scale(1./vhist8a.Integral());
	vhist9c.Scale(1./vhist9a.Integral());
	vhist6d.Scale(1./vhist6b.Integral());
	vhist7d.Scale(1./vhist7b.Integral());
	vhist8d.Scale(1./vhist8b.Integral());
	vhist9d.Scale(1./vhist9b.Integral());

	vhist0a.Scale(1./vhist0a.Integral());
	vhist1a.Scale(1./vhist1a.Integral());
	vhist3a.Scale(1./vhist3a.Integral());
	vhist4a.Scale(1./vhist4a.Integral());
	vhist5z.Scale(1./vhist5z.Integral());
	vhist5a.Scale(1./vhist5a.Integral());
	vhist6a.Scale(1./vhist6a.Integral());
	vhist7a.Scale(1./vhist7a.Integral());
	vhist8a.Scale(1./vhist8a.Integral());
	vhist9a.Scale(1./vhist9a.Integral());

	vhist0b.Scale(1./vhist0b.Integral());
	vhist1b.Scale(1./vhist1b.Integral());
	vhist3b.Scale(1./vhist3b.Integral());
	vhist4b.Scale(1./vhist4b.Integral());
	vhist5b.Scale(1./vhist5b.Integral());
	vhist6b.Scale(1./vhist6b.Integral());
	vhist7b.Scale(1./vhist7b.Integral());
	vhist8b.Scale(1./vhist8b.Integral());
	vhist9b.Scale(1./vhist9b.Integral());
	vhist10b.Scale(0.01);

	vhist0a.GetXaxis()->SetTitle("#DeltaX");
	vhist1a.GetXaxis()->SetTitle("#DeltaY");
	vhist3a.GetXaxis()->SetTitle("N_{trk}");
	vhist4a.GetXaxis()->SetTitle("t_{z}");
	vhist5a.GetXaxis()->SetTitle("fv_{z}");
	vhist6a.GetXaxis()->SetTitle("m");
	vhist7a.GetXaxis()->SetTitle("p_{z}");
	vhist8a.GetXaxis()->SetTitle("m_{cor}");
	vhist9a.GetXaxis()->SetTitle("DOCA");
	vhist10a.GetXaxis()->SetTitle("nSV");

	c.SetLogx(0);
	c.SetLogy(0);
	vhist0a.SetLineColor(kRed);
	vhist0a.Draw();
	vhist0b.SetLineColor(kMagenta);
	vhist0b.Draw("same");
	vhist0c.SetLineColor(kRed);
	vhist0c.SetLineStyle(kDashed);
	vhist0c.Draw("same");
	vhist0d.SetLineColor(kMagenta);
	vhist0d.SetLineStyle(kDashed);
	vhist0d.Draw("same");
	c.SaveAs("deltaX.pdf");
	vhist1a.SetLineColor(kRed);
	vhist1a.Draw();
	vhist1b.SetLineColor(kMagenta);
	vhist1b.Draw("same");
	vhist1c.SetLineColor(kRed);
	vhist1c.SetLineStyle(kDashed);
	vhist1c.Draw("same");
	vhist1d.SetLineColor(kMagenta);
	vhist1d.SetLineStyle(kDashed);
	vhist1d.Draw("same");
	c.SaveAs("deltaY.pdf");
	
	vhist3a.SetLineColor(kRed);
	vhist3a.Draw();
	vhist3b.SetLineColor(kMagenta);
	vhist3b.Draw("same");
	c.SaveAs("ntrk.pdf");

	vhist4a.SetLineColor(kRed);
	vhist4a.Draw();
	vhist4b.SetLineColor(kMagenta);
	vhist4b.Draw("same");
	vhist4c.SetLineColor(kRed);
	vhist4c.SetLineStyle(kDashed);
	vhist4c.Draw("same");
	vhist4d.SetLineColor(kMagenta);
	vhist4d.SetLineStyle(kDashed);
	vhist4d.Draw("same");
	c.SaveAs("tz.pdf");

	if(vhist5b.GetMaximum()>vhist5a.GetMaximum()) vhist5a.SetMaximum(1.1*vhist5b.GetMaximum());
	if(vhist5z.GetMaximum()>vhist5a.GetMaximum()) vhist5a.SetMaximum(1.1*vhist5z.GetMaximum());
	vhist5a.SetLineColor(kRed);
	vhist5a.Draw();
	vhist5b.SetLineColor(kMagenta);
	vhist5b.Draw("same");
	vhist5z.SetLineColor(kBlue);
	vhist5z.Draw("same");
	vhist5c.SetLineColor(kRed);
	vhist5c.SetLineStyle(kDashed);
	vhist5c.Draw("same");
	vhist5d.SetLineColor(kMagenta);
	vhist5d.SetLineStyle(kDashed);
	vhist5d.Draw("same");
	c.SaveAs("fvz.pdf");

	vhist6a.SetLineColor(kRed);
	vhist6a.Draw();
	vhist6b.SetLineColor(kMagenta);
	vhist6b.Draw("same");
	vhist6c.SetLineColor(kRed);
	vhist6c.SetLineStyle(kDashed);
	vhist6c.Draw("same");
	vhist6d.SetLineColor(kMagenta);
	vhist6d.SetLineStyle(kDashed);
	vhist6d.Draw("same");
	c.SaveAs("m.pdf");

	vhist7a.SetLineColor(kRed);
	vhist7a.Draw();
	vhist7b.SetLineColor(kMagenta);
	vhist7b.Draw("same");
	vhist7c.SetLineColor(kRed);
	vhist7c.SetLineStyle(kDashed);
	vhist7c.Draw("same");
	vhist7d.SetLineColor(kMagenta);
	vhist7d.SetLineStyle(kDashed);
	vhist7d.Draw("same");
	c.SaveAs("pz.pdf");

	vhist8a.SetLineColor(kRed);
	vhist8a.Draw();
	vhist8b.SetLineColor(kMagenta);
	vhist8b.Draw("same");
	vhist8c.SetLineColor(kRed);
	vhist8c.SetLineStyle(kDashed);
	vhist8c.Draw("same");
	vhist8d.SetLineColor(kMagenta);
	vhist8d.SetLineStyle(kDashed);
	vhist8d.Draw("same");
	c.SaveAs("mcor.pdf");

	vhist9a.SetLineColor(kRed);
	vhist9a.Draw();
	vhist9b.SetLineColor(kMagenta);
	vhist9b.Draw("same");
	vhist9c.SetLineColor(kRed);
	vhist9c.SetLineStyle(kDashed);
	vhist9c.Draw("same");
	vhist9d.SetLineColor(kMagenta);
	vhist9d.SetLineStyle(kDashed);
	vhist9d.Draw("same");
	c.SaveAs("doca.pdf");

	vhist10a.SetLineColor(kRed);
	vhist10a.Draw();
	vhist10b.SetLineColor(kMagenta);
	vhist10b.Draw("same");
	c.SaveAs("nsvr.pdf");

	hist21.Draw();
	c.SaveAs("pvShiftX.pdf");
	hist22.Draw();
	c.SaveAs("pvShiftY.pdf");
	hist23.Draw();
	c.SaveAs("pvShiftZ.pdf");
}

//usage - run with no argument to compare generated to (100x) resampled tracks and produce histograms
//      - run with argument (n) to produce n resampled datasets
int main(int argc, char** argv) {
	int max(-1);
	if(argc>1) max=atoi(argv[1]);

	if(max<0) {
		resampleTrackIPs a;
		a.Loop();
	} else {
		for(int i=0; i<=max; ++i) {
			std::cout << i << " of " << max << std::endl;
			resampleTrackIPs a(i);
			a.Loop();
		}
	}
	return 0;
}
