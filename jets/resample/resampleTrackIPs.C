#define resampleTrackIPs_cxx
#include "resampleTrackIPs.h"
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

int resampleTrackIPs::getNSharedVeloHits(int idxi, int idxj) {
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

int resampleTrackIPs::getNSharedEarlyVeloHits(int idxi, int idxj) {
	int nshared(0);
	if(trk_vhit0->at(idxi) != -1 && trk_vhit0->at(idxi) == trk_vhit0->at(idxj)) ++nshared;
	if(trk_vhit1->at(idxi) != -1 && trk_vhit1->at(idxi) == trk_vhit1->at(idxj)) ++nshared;
	if(trk_vhit2->at(idxi) != -1 && trk_vhit2->at(idxi) == trk_vhit2->at(idxj)) ++nshared;
	if(trk_vhit3->at(idxi) != -1 && trk_vhit3->at(idxi) == trk_vhit3->at(idxj)) ++nshared;

	return nshared;
}

int resampleTrackIPs::findTrkInJet(int trk_idx, int jet_idx) {
	if( jet_idx_trk0->at(jet_idx) == trk_idx) return 0;
	if( jet_idx_trk1->at(jet_idx) == trk_idx) return 1;
	if( jet_idx_trk2->at(jet_idx) == trk_idx) return 2;
	if( jet_idx_trk3->at(jet_idx) == trk_idx) return 3;
	if( jet_idx_trk4->at(jet_idx) == trk_idx) return 4;
	if( jet_idx_trk5->at(jet_idx) == trk_idx) return 5;
	if( jet_idx_trk6->at(jet_idx) == trk_idx) return 6;
	if( jet_idx_trk7->at(jet_idx) == trk_idx) return 7;
	if( jet_idx_trk8->at(jet_idx) == trk_idx) return 8;
	if( jet_idx_trk9->at(jet_idx) == trk_idx) return 9;
	if(jet_idx_trk10->at(jet_idx) == trk_idx) return 10;
	if(jet_idx_trk11->at(jet_idx) == trk_idx) return 11;
	if(jet_idx_trk12->at(jet_idx) == trk_idx) return 12;
	if(jet_idx_trk13->at(jet_idx) == trk_idx) return 13;
	if(jet_idx_trk14->at(jet_idx) == trk_idx) return 14;
	if(jet_idx_trk15->at(jet_idx) == trk_idx) return 15;
	if(jet_idx_trk16->at(jet_idx) == trk_idx) return 16;
	if(jet_idx_trk17->at(jet_idx) == trk_idx) return 17;
	if(jet_idx_trk18->at(jet_idx) == trk_idx) return 18;
	if(jet_idx_trk19->at(jet_idx) == trk_idx) return 19;
	if(jet_idx_trk20->at(jet_idx) == trk_idx) return 20;
	if(jet_idx_trk21->at(jet_idx) == trk_idx) return 21;
	if(jet_idx_trk22->at(jet_idx) == trk_idx) return 22;
	if(jet_idx_trk23->at(jet_idx) == trk_idx) return 23;
	if(jet_idx_trk24->at(jet_idx) == trk_idx) return 24;
	if(jet_idx_trk25->at(jet_idx) == trk_idx) return 25;
	if(jet_idx_trk26->at(jet_idx) == trk_idx) return 26;
	if(jet_idx_trk27->at(jet_idx) == trk_idx) return 27;
	if(jet_idx_trk28->at(jet_idx) == trk_idx) return 28;
	if(jet_idx_trk29->at(jet_idx) == trk_idx) return 29;
	if(jet_idx_trk30->at(jet_idx) == trk_idx) return 30;
	if(jet_idx_trk31->at(jet_idx) == trk_idx) return 31;
	if(jet_idx_trk32->at(jet_idx) == trk_idx) return 32;
	if(jet_idx_trk33->at(jet_idx) == trk_idx) return 33;
	if(jet_idx_trk34->at(jet_idx) == trk_idx) return 34;
	if(jet_idx_trk35->at(jet_idx) == trk_idx) return 35;
	if(jet_idx_trk36->at(jet_idx) == trk_idx) return 36;
	if(jet_idx_trk37->at(jet_idx) == trk_idx) return 37;
	if(jet_idx_trk38->at(jet_idx) == trk_idx) return 38;
	if(jet_idx_trk39->at(jet_idx) == trk_idx) return 39;
	
	return -1;
}

int resampleTrackIPs::findBestGenPvrForTrk(int idx, double x, double y) {
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
		if(x!=-999.) pp.SetX(x);
		if(y!=-999.) pp.SetY(y);
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


int resampleTrackIPs::findBestRecPvrForTrk(int idx, double x, double y) {
	int pvr_idx(-1);

     	//assign the PV with the smallest DOCA
     	if(pvr_idx<0) {
     		int bestpvr(-1);
     		double minDOCA=9999.;
     		TVector3 pp(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));
		if(x!=-999.) pp.SetX(x);
		if(y!=-999.) pp.SetY(y);
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
     		//this might be generated so look for close pv with higher index
     		double minDist(2.);
     		for(uint ipvr=bestpvr+1; ipvr<pvr_z->size(); ++ipvr) {
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

int resampleTrackIPs::findBestPvrForSvr(TVector3 x, TVector3 p, bool gen) {
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

void resampleTrackIPs::printTrk(int itrk, int isvr) {
	int bestpvr = findBestGenPvrForTrk(itrk);
	int bestrecopvr = findBestRecPvrForTrk(itrk);
	int pvr = trk_idx_pvr->at(itrk);
	int svrpvr = svr_idx_pvr->at(isvr);

	TVector3 pi(trk_px->at(itrk), trk_py->at(itrk), trk_pz->at(itrk));
	TVector3 xi(trk_x->at(itrk), trk_y->at(itrk), trk_z->at(itrk));
	TVector3 owngenpv(pvr_x->at(bestpvr),pvr_y->at(bestpvr),pvr_z->at(bestpvr));
	TVector3 ownrecpv(pvr_x->at(bestrecopvr),pvr_y->at(bestrecopvr),pvr_z->at(bestrecopvr));
	TVector3 svrpv(pvr_x->at(svrpvr),pvr_y->at(svrpvr),pvr_z->at(svrpvr));
	TVector3 pv(pvr_x->at(pvr),pvr_y->at(pvr),pvr_z->at(pvr));

	double sigma2 = trk_ip->at(itrk)*trk_ip->at(itrk) / trk_ip_chi2->at(itrk);
	double sigma = TMath::Sqrt(sigma2);

	TVector3 dxdyrec = calcDxDy(xi, pi, ownrecpv);
	TVector3 dxdygen = calcDxDy(xi, pi, owngenpv);
	TVector3 dxdysvr = calcDxDy(xi, pi, svrpv);
	TVector3 dxdy    = calcDxDy(xi, pi, pv);

	double ipchi2rec = TMath::Power(dxdyrec.Mag() / sigma, 2);
	double ipchi2gen = TMath::Power(dxdygen.Mag() / sigma, 2);
	double ipchi2svr = TMath::Power(dxdysvr.Mag() / sigma, 2);
	double ipchi2    = TMath::Power(dxdy.Mag() / sigma, 2);

	printf("%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",bestpvr, bestrecopvr, pvr, svrpvr, TMath::Log(trk_ip_chi2->at(itrk)), TMath::Log(ipchi2rec), TMath::Log(ipchi2gen), TMath::Log(ipchi2svr), TMath::Log(ipchi2));

}

void resampleTrackIPs::printSvr(int isvr) {
	if(svr_idx_trk0->at(isvr)!=-1) printTrk(svr_idx_trk0->at(isvr), isvr);
	if(svr_idx_trk1->at(isvr)!=-1) printTrk(svr_idx_trk1->at(isvr), isvr);
	if(svr_idx_trk2->at(isvr)!=-1) printTrk(svr_idx_trk2->at(isvr), isvr);
	if(svr_idx_trk3->at(isvr)!=-1) printTrk(svr_idx_trk3->at(isvr), isvr);
	if(svr_idx_trk4->at(isvr)!=-1) printTrk(svr_idx_trk4->at(isvr), isvr);
	if(svr_idx_trk5->at(isvr)!=-1) printTrk(svr_idx_trk5->at(isvr), isvr);
	if(svr_idx_trk6->at(isvr)!=-1) printTrk(svr_idx_trk6->at(isvr), isvr);
	if(svr_idx_trk7->at(isvr)!=-1) printTrk(svr_idx_trk7->at(isvr), isvr);
	if(svr_idx_trk8->at(isvr)!=-1) printTrk(svr_idx_trk8->at(isvr), isvr);
	if(svr_idx_trk9->at(isvr)!=-1) printTrk(svr_idx_trk9->at(isvr), isvr);
	std::cout << std::endl;
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
	TH1D hist5("hist5", "", 100, 0, 20);
	TH1D hist6("hist6", "", 100, 0, 20);
	TH1D hist7("hist7", "", 100, 0, 20);
	TH1D hist5b("hist5b", "", 100, 0, 20);
	TH1D hist6b("hist6b", "", 100, 0, 20);
	TH1D hist7b("hist7b", "", 100, 0, 20);
	//IP (chi2<16)
	TH1D hist8("hist8", "", 100, 0, 1);
	TH1D hist9("hist9", "", 100, 0, 1);
	TH1D hist10("hist10", "", 100, 0, 1);
	//IP (chi2>16)
	TH1D hist11("hist11", "", 100, 0, 1);
	TH1D hist12("hist12", "", 100, 0, 1);
	TH1D hist13("hist13", "", 100, 0, 1);
	//IPx
	TH1D hist14("hist14", "", 100, -0.1, 0.1);
	TH1D hist15("hist15", "", 100, -0.1, 0.1);
	//IPy
	TH1D hist16("hist16", "", 100, -0.1, 0.1);
	TH1D hist17("hist17", "", 100, -0.1, 0.1);

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

	//Sum IPchi2
	TH1D hist24("hist24", "", 100, 0, 50);
	TH1D hist25("hist25", "", 100, 0, 50);
	TH1D hist26("hist26", "", 100, 0, 50);

	TH1D histA0("histA0", "", 10, 0, 10);
	TH1D histA1("histA1", "", 10, 0, 10);
	TH1D histA2("histA2", "", 10, 0, 10);
	TH1D histB0("histB0", "", 10, 0, 10);
	TH1D histB1("histB1", "", 10, 0, 10);
	TH1D histB2("histB2", "", 10, 0, 10);
	TH1D histC0("histC0", "", 10, 0, 10);
	TH1D histC1("histC1", "", 10, 0, 10);
	TH1D histC2("histC2", "", 10, 0, 10);
	TH1D histD0("histD0", "", 10, 0, 10);
	TH1D histD1("histD1", "", 10, 0, 10);
	TH1D histD2("histD2", "", 10, 0, 10);
	TH1D histE0("histE0", "", 10, 0, 10);
	TH1D histE1("histE1", "", 10, 0, 10);
	TH1D histE2("histE2", "", 10, 0, 10);

	//track separation x,y and angle
	TH1D vhist0a("vhist0a","",20,0,5.);//40);
	TH1D vhist0b("vhist0b","",20,0,5.);//40);
	TH1D vhist1a("vhist1a","",20,0,5.);//40);
	TH1D vhist1b("vhist1b","",20,0,5.);//40);
	TH1D vhist2a("vhist2a","",40,0,.5);
	TH1D vhist2b("vhist2b","",40,0,.5);

	TH1D vhist0c("vhist0c","",20,0,5.);//40);
	TH1D vhist0d("vhist0d","",20,0,5.);//40);
	TH1D vhist1c("vhist1c","",20,0,5.);//40);
	TH1D vhist1d("vhist1d","",20,0,5.);//40);
	TH1D vhist2c("vhist2c","",40,0,.5);
	TH1D vhist2d("vhist2d","",40,0,.5);

	TH2D vhist2e("vhist2e","",20,0,.5,10,0.,5.);
	TH2D vhist2f("vhist2f","",20,0,.5,10,0.,5.);
	TH2D vhist2g("vhist2g","",20,0,.5,10,0.,5.);
	TH2D vhist2h("vhist2h","",20,0,.5,10,0.,5.);

	TH2D vhist2i("vhist2i","",20,0,.5,15,0.,1500.);
	TH2D vhist2j("vhist2j","",20,0,.5,15,0.,1500.);
	TH2D vhist2k("vhist2k","",20,0,.5,15,0.,1500.);
	TH2D vhist2l("vhist2l","",20,0,.5,15,0.,1500.);

	//Ntracks
	TH1D vhist3a("vhist3a","",50,0,50);
	TH1D vhist3b("vhist3b","",50,0,50);

	//tz
	TH1D vhist4a("vhist4a","",50,-1,10);
	TH1D vhist4b("vhist4b","",50,-1,10);
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

	//N shared hits
	TH1D vhist11a("vhist11a","",30,0.,5.);
	TH1D vhist11b("vhist11b","",30,0.,5.);

	//fvz vs difference in shift
	TH2D vhist12a("vhist12a","",10,-2.,2.,10,0.,0.3);
	TH2D vhist12b("vhist12b","",10,-2.,2.,10,0.,0.3);

	//N in jet
	TH1D vhist13a("vhist13a","",6,-0.5,5.5);
	TH1D vhist13b("vhist13b","",6,-0.5,5.5);

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
	TFile* fh = TFile::Open("ipmap-new.root");
	TH2D* h = dynamic_cast<TH2D*>(fh->Get("prompt"));//andnontruth"));
	std::vector<TH1D*> hs;
	for(int i=1; i<=h->GetNbinsY(); ++i) {
		TString name = "_bin"; name+=i;
		hs.push_back(h->ProjectionX(name,i,i));
	}
	TH2D* hNT = dynamic_cast<TH2D*>(fh->Get("nontruth"));//promptandnontruth"));
	std::vector<TH1D*> hsNT;
	for(int i=1; i<=h->GetNbinsY(); ++i) {
		TString name = "_binNT"; name+=i;
		hsNT.push_back(hNT->ProjectionX(name,i,i));
	}
	TH2D* hNP = dynamic_cast<TH2D*>(fh->Get("nonprompt"));
	std::vector<TH1D*> hsNP;
	for(int i=1; i<=h->GetNbinsY(); ++i) {
		TString name = "_binNP"; name+=i;
		hsNP.push_back(hNP->ProjectionX(name,i,i));
	}
	TH1D* hphi0 = dynamic_cast<TH2D*>(fh->Get("deltaDR"))->ProjectionX("hphi0",1,1);
	TH1D* hphi1 = dynamic_cast<TH2D*>(fh->Get("deltaDR"))->ProjectionX("hphi1",2,2);
	TH1D* hphi2 = dynamic_cast<TH2D*>(fh->Get("deltaDR"))->ProjectionX("hphi2",3,3);
	hphi1->Sumw2();
	hphi2->Sumw2();
	hphi1->Divide(hphi0);
	hphi2->Divide(hphi0);
	TH1D* hratio0 = dynamic_cast<TH2D*>(fh->Get("deltaDR2"))->ProjectionX("hratio0",1,1);
	TH1D* hratio1 = dynamic_cast<TH2D*>(fh->Get("deltaDR2"))->ProjectionX("hratio1",2,2);
	TH1D* hratio2 = dynamic_cast<TH2D*>(fh->Get("deltaDR2"))->ProjectionX("hratio2",3,3);
	hratio1->Sumw2();
	hratio2->Sumw2();
	hratio1->Divide(hratio0);
	hratio2->Divide(hratio0);

	//keep track of various stats
	int countBefore(0), countAfter(0), total(0), nonPrompt(0), noTruth(0), realTotal(0);
	int svsBefore(0), svsAfter1(0), svsAfter2(0), svNsAfter1(0), svNsAfter2(0), sv2sAfter1(0), sv2sAfter2(0), sv3sAfter1(0), sv3sAfter2(0), sv4sAfter1(0), sv4sAfter2(0), svMoresAfter1(0), svMoresAfter2(0);
	int nLinkedSVs1(0), nLinkedSVs2(0);
	int nLinkedSVs11(0), nLinkedSVs21(0);
	int nLinkedSVs12(0), nLinkedSVs22(0);
	int nLinkedSVs13(0), nLinkedSVs23(0);
	int nLinkedSVs14(0), nLinkedSVs24(0);
	int nLinkedSVs15(0), nLinkedSVs25(0);
	int nLinkedSVs16(0), nLinkedSVs26(0);
	int n2BSVs10(0), n2BSVs20(0);
	int n2BSVs11(0), n2BSVs21(0);
	int n2BSVs12(0), n2BSVs22(0);
	int n2BSVs13(0), n2BSVs23(0);
	int n2BSVs14(0), n2BSVs24(0);
	int noTrks1(0), noTrks2(0), oneTrk1(0), oneTrk2(0), twoTrks1(0), twoTrks2(0);
	int nMatchGen1(0), nMatchGen2(0);
	int rerunsRequired(0);

	int nReruns[50] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int reruns(0);

	//events that pass online but fail offline
	std::set<int> inspect1;
	inspect1.insert(  608);	inspect1.insert( 4084);	inspect1.insert( 4146);	inspect1.insert( 8504);	inspect1.insert( 9439);	inspect1.insert( 9623);	inspect1.insert(10754);
	inspect1.insert(11926);	inspect1.insert(17877);	inspect1.insert(20383);	inspect1.insert(20815);	inspect1.insert(22037);	inspect1.insert(27006);	inspect1.insert(28301);
	inspect1.insert(29516);	inspect1.insert(29549);	inspect1.insert(32035);	inspect1.insert(36053);	inspect1.insert(40293);	inspect1.insert(40637);	inspect1.insert(43009);
	inspect1.insert(44522);	inspect1.insert(58061);	inspect1.insert(66410);

	//events that fail online but pass offline
	std::set<int> inspect2;
	inspect2.insert(  623);	inspect2.insert(  642);	inspect2.insert(  690);	inspect2.insert(  947);	inspect2.insert( 1587);	inspect2.insert( 1625);	inspect2.insert( 2057);
	inspect2.insert( 2074);	inspect2.insert( 2285);	inspect2.insert( 2607);	inspect2.insert( 3134);	inspect2.insert( 3210);	inspect2.insert( 3465);	inspect2.insert( 4193);
	inspect2.insert( 4422);	inspect2.insert( 5375);	inspect2.insert( 5402);	inspect2.insert( 5561);	inspect2.insert( 5764);	inspect2.insert( 5954);	inspect2.insert( 6003);
	inspect2.insert( 6129);	inspect2.insert( 6454);	inspect2.insert( 6660);	inspect2.insert( 7002);	inspect2.insert( 7158);	inspect2.insert( 7423);	inspect2.insert( 7546);
	inspect2.insert( 7979);	inspect2.insert( 9036);	inspect2.insert( 9437);	inspect2.insert( 9511);	inspect2.insert( 9985);	inspect2.insert(10539);	inspect2.insert(10698);
	inspect2.insert(11225);	inspect2.insert(11403);	inspect2.insert(11485);	inspect2.insert(11653);	inspect2.insert(12279);	inspect2.insert(13060);	inspect2.insert(13522);
	inspect2.insert(13727);	inspect2.insert(13903);	inspect2.insert(13969);	inspect2.insert(14223);	inspect2.insert(14873);	inspect2.insert(15096);	inspect2.insert(15424);
	inspect2.insert(15688);	inspect2.insert(15700);	inspect2.insert(15891);	inspect2.insert(16719);	inspect2.insert(16946);	inspect2.insert(17065);	inspect2.insert(17146);
	inspect2.insert(17760);	inspect2.insert(17909);	inspect2.insert(18017);	inspect2.insert(18621);	inspect2.insert(18831);	inspect2.insert(18953);	inspect2.insert(19049);
	inspect2.insert(19165);	inspect2.insert(19611);	inspect2.insert(19927);	inspect2.insert(20033);	inspect2.insert(20241);	inspect2.insert(20357);	inspect2.insert(20404);
	inspect2.insert(20577);	inspect2.insert(20903);	inspect2.insert(20913);	inspect2.insert(20928);	inspect2.insert(20971);	inspect2.insert(21133);	inspect2.insert(21170);
	inspect2.insert(21334);	inspect2.insert(21398);	inspect2.insert(21725);	inspect2.insert(21790);	inspect2.insert(21838);	inspect2.insert(23113);	inspect2.insert(23368);
	inspect2.insert(23444);	inspect2.insert(23499);	inspect2.insert(24243);	inspect2.insert(24512);	inspect2.insert(25022);	inspect2.insert(25100);	inspect2.insert(26256);
	inspect2.insert(26312);	inspect2.insert(26731);	inspect2.insert(27138);	inspect2.insert(27265);	inspect2.insert(27498);	inspect2.insert(27541);	inspect2.insert(28488);
	inspect2.insert(28596);	inspect2.insert(28737);	inspect2.insert(28763);	inspect2.insert(29015);	inspect2.insert(29582);	inspect2.insert(29709);	inspect2.insert(29821);
	inspect2.insert(30060);	inspect2.insert(30142);	inspect2.insert(30922);	inspect2.insert(31009);	inspect2.insert(31064);	inspect2.insert(31264);	inspect2.insert(31473);
	inspect2.insert(31479);	inspect2.insert(31979);	inspect2.insert(32659);	inspect2.insert(33008);	inspect2.insert(33192);	inspect2.insert(33242);	inspect2.insert(33727);
	inspect2.insert(33860);	inspect2.insert(34350);	inspect2.insert(34392);	inspect2.insert(34652);	inspect2.insert(35786);	inspect2.insert(36874);	inspect2.insert(37454);
	inspect2.insert(38996);	inspect2.insert(39021);	inspect2.insert(39196);	inspect2.insert(39243);	inspect2.insert(40012);	inspect2.insert(40182);	inspect2.insert(40547);
	inspect2.insert(40702);	inspect2.insert(41508);	inspect2.insert(41701);	inspect2.insert(41945);	inspect2.insert(42077);	inspect2.insert(42256);	inspect2.insert(42734);
	inspect2.insert(42766);	inspect2.insert(43094);	inspect2.insert(43214);	inspect2.insert(43523);	inspect2.insert(43808);	inspect2.insert(44128);	inspect2.insert(44211);
	inspect2.insert(44381);	inspect2.insert(44604);	inspect2.insert(44668);	inspect2.insert(44702);	inspect2.insert(45118);	inspect2.insert(45863);	inspect2.insert(45959);
	inspect2.insert(46207);	inspect2.insert(46390);	inspect2.insert(46404);	inspect2.insert(46417);	inspect2.insert(47141);	inspect2.insert(47161);	inspect2.insert(47611);
	inspect2.insert(47763);	inspect2.insert(47895);	inspect2.insert(48373);	inspect2.insert(49027);	inspect2.insert(49085);	inspect2.insert(49253);	inspect2.insert(50073);
	inspect2.insert(50239);	inspect2.insert(50505);	inspect2.insert(50793);	inspect2.insert(51242);	inspect2.insert(51297);	inspect2.insert(51337);	inspect2.insert(51698);
	inspect2.insert(51899);	inspect2.insert(52195);	inspect2.insert(52491);	inspect2.insert(52535);	inspect2.insert(52563);	inspect2.insert(52597);	inspect2.insert(53303);
	inspect2.insert(53547);	inspect2.insert(54053);	inspect2.insert(54127);	inspect2.insert(54308);	inspect2.insert(54541);	inspect2.insert(54747);	inspect2.insert(54803);
	inspect2.insert(55714);	inspect2.insert(55803);	inspect2.insert(55885);	inspect2.insert(56297);	inspect2.insert(56523);	inspect2.insert(57012);	inspect2.insert(57026);
	inspect2.insert(57453);	inspect2.insert(57708);	inspect2.insert(57744);	inspect2.insert(57825);	inspect2.insert(57845);	inspect2.insert(57894);	inspect2.insert(58495);
	inspect2.insert(58607);	inspect2.insert(58824);	inspect2.insert(58922);	inspect2.insert(58956);	inspect2.insert(58970);	inspect2.insert(59355);	inspect2.insert(60177);
	inspect2.insert(60402);	inspect2.insert(60465);	inspect2.insert(60488);	inspect2.insert(61137);	inspect2.insert(61277);	inspect2.insert(61735);	inspect2.insert(62103);
	inspect2.insert(62446);	inspect2.insert(62735);	inspect2.insert(62953);	inspect2.insert(63069);	inspect2.insert(63886);	inspect2.insert(64142);	inspect2.insert(64537);
	inspect2.insert(64879);	inspect2.insert(64883);	inspect2.insert(65117);	inspect2.insert(65350);	inspect2.insert(65536);	inspect2.insert(65990);	inspect2.insert(66403);

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

		std::set<int> trksToWatch;

		svsBefore += svr_z->size();
		for(uint isvr=0; isvr<svr_z->size(); ++isvr) {
			pvr_idx = svr_idx_pvr->at(isvr);

			vhist5z.Fill(svr_z->at(isvr) - pvr_z->at(pvr_idx));
			//printSvr(isvr);

			if(svr_idx_trk0->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk1->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk2->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk2->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk2->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk3->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk3->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk3->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk4->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk4->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk4->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk5->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk5->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk5->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk6->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk6->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk6->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk7->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk7->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk7->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk8->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk8->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk8->at(isvr)))>8) printSvr(isvr);}
			if(svr_idx_trk9->at(isvr)!=-1) {hist5b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk9->at(isvr)))); }//if(TMath::Log(trk_ip_chi2->at(svr_idx_trk9->at(isvr)))>8) printSvr(isvr);}

			//if(svr_idx_trk0->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk0->at(isvr));
			//if(svr_idx_trk1->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk1->at(isvr));
			//if(svr_idx_trk2->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk2->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk2->at(isvr));
			//if(svr_idx_trk3->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk3->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk3->at(isvr));
			//if(svr_idx_trk4->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk4->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk4->at(isvr));
			//if(svr_idx_trk5->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk5->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk5->at(isvr));
			//if(svr_idx_trk6->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk6->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk6->at(isvr));
			//if(svr_idx_trk7->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk7->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk7->at(isvr));
			//if(svr_idx_trk8->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk8->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk8->at(isvr));
			//if(svr_idx_trk9->at(isvr)!=-1 && TMath::Log(trk_ip_chi2->at(svr_idx_trk9->at(isvr)))>10.) trksToWatch.insert(svr_idx_trk9->at(isvr));
		}

		for(std::set<int>::iterator it=trksToWatch.begin(); it!=trksToWatch.end(); ++it) {
			std::cout << "DEBUG track " << (*it) << " is being watched (log(IPchi2)=" << TMath::Log(trk_ip_chi2->at(*it)) << ")" << std::endl;
		}

		std::vector<double>* origtrk_x = new std::vector<double>();
		std::vector<double>* origtrk_y = new std::vector<double>();
		std::vector<double>* origtrk_ip = new std::vector<double>();
		std::vector<double>* origtrk_ip_chi2 = new std::vector<double>();
		for(uint i=0; i<trk_x->size(); ++i) {
			origtrk_x->push_back(trk_x->at(i));
			origtrk_y->push_back(trk_y->at(i));
			origtrk_ip->push_back(trk_ip->at(i));
			origtrk_ip_chi2->push_back(trk_ip_chi2->at(i));
		}

		std::map<int,std::vector<int>> passedTracksByRep;
		//resample each event as required
		//loop zero is unaltered
		//when filling TTrees we only run this loop once to fill a single TTree
		for(int irep=firstrep; irep<endrep; ++irep) {
			bool rerunEvt(false);

			double sumIPchi2(0.);
			double trueSumIPchi2(0.);
			int nInIPchi2Bin[5] ={0,0,0,0,0};
			int nTrueInIPchi2Bin[5] ={0,0,0,0,0};

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
			svr_pass->clear();

			std::vector<int> trks;
			std::vector<TVector3> trkhit, trkdir;

			//store list of prompt tracks for each pvr to perform opening angle study
			std::map<int, std::vector<int> > promptTrksByPV;

			for(uint idx=0; idx<trk_idx_gen->size(); ++idx) {
				if(trksToWatch.count(idx)) std::cout << "DEBUG track " << idx << " START" << std::endl;
				if(trk_type->at(idx) != 3) continue;
				if(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx)<500*500 || trk_prb_ghost->at(idx)>0.3) continue;
				if(irep==0) ++realTotal;
				if(trksToWatch.count(idx)) std::cout << "DEBUG track " << idx << " PASS PRESEL" << std::endl;

				bool isPrompt(false), isTruth(false);
				int gen_idx = trk_idx_gen->at(idx);

				//if we have truth level information we can separate prompt and non-prompt tracks
				if(gen_idx>=0) isTruth=true;
				else ++noTruth;

				pvr_idx = -1;
				//if we have truth-level information try to get the PV from this.
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
					TVector3 pp(origtrk_x->at(idx), origtrk_y->at(idx), trk_z->at(idx));
					TVector3 dp(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
					for(uint ipvr=0; ipvr<pvr_z->size(); ++ipvr) {
						TVector3 pv(pvr_x->at(ipvr),pvr_y->at(ipvr),pvr_z->at(ipvr));
						double doca = calcDocaPoint(pp,dp,pv);
						if(doca<minDOCA) {
							bestpvr = ipvr;
							minDOCA = doca;
							//							std::cout << "best doca so far is " << minDOCA << std::endl;
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
					std::cout << "WARNING - shouldn't reach here! (1)" << std::endl;
					continue;
				}
				if(isTruth) {
					if(TMath::Abs(gen_x->at(gen_idx) - pvr_x->at(pvr_idx)) < 1e-5 &&
							TMath::Abs(gen_y->at(gen_idx) - pvr_y->at(pvr_idx)) < 1e-5 &&
							TMath::Abs(gen_z->at(gen_idx) - pvr_z->at(pvr_idx)) < 1e-5) isPrompt=true;
					else ++nonPrompt;
				}

				invpt = 1./TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx));

				tx = trk_px->at(idx) / trk_pz->at(idx);
				ty = trk_py->at(idx) / trk_pz->at(idx);

				dx = origtrk_x->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*tx - pvr_x->at(pvr_idx);
				dy = origtrk_y->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*ty - pvr_y->at(pvr_idx);

				binInvPt = h->GetYaxis()->FindBin(invpt);
				sigma2 = trk_ip->at(idx)*trk_ip->at(idx) / trk_ip_chi2->at(idx); //assume unchanged by shift
				sigma = TMath::Sqrt(sigma2);
				
				if(sigma!=sigma) continue;
				if(trksToWatch.count(idx)) std::cout << "DEBUG track " << idx << " PASS SIGMA" << std::endl;
				//if(!isPrompt) continue;//TODO
				if(trksToWatch.count(idx)) std::cout << "DEBUG track " << idx << " PASS PROMPT" << std::endl;

				++total;

				newDx = dx;
				newDy = dy;

				//resample tracks identified as prompt and tracks without truth information
				if(irep>0 && (isPrompt || !isTruth) ) {
					double phi=r.Rndm()*2*TMath::Pi();
					double rho(0);
					if(!isTruth) {
						rho = sigma*getRandom(hsNT[binInvPt-1],&r); //interpolate rho within the bin
					} else if(isPrompt) {
						rho = sigma*getRandom(hs[binInvPt-1],&r); //interpolate rho within the bin
					} else {
						rho = sigma*getRandom(hsNP[binInvPt-1],&r); //interpolate rho within the bin
					}

					newDx = rho*TMath::Sin(phi);
					newDy = rho*TMath::Cos(phi);

					//if(newDx!=newDx) std::cout << rho << "\t" << phi << "\t" << sigma << std::endl;//TODO

					trk_x->at(idx) = origtrk_x->at(idx) + newDx - dx;
					trk_y->at(idx) = origtrk_y->at(idx) + newDy - dy;
				}
				
				trk_ip->at(idx) = TMath::Sqrt(newDx*newDx + newDy*newDy);
				trk_ip_chi2->at(idx) = (newDx*newDx + newDy*newDy) / sigma2;

				//if this track shares any hits with earlier tracks then check the ratio of and angle between the shifts is good
				TVector3 pi(trk_px->at(idx), trk_py->at(idx), trk_pz->at(idx));
				TVector3 xi(trk_x->at(idx), trk_y->at(idx), trk_z->at(idx));
				TVector3 pv(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx));

				TVector3 dxdyi = calcDxDy(xi, pi, pv);

				//attempt to apply correlation corrections
				if(false) { //TODO isPrompt || !isTruth) {
					for(uint jtrk=0; jtrk<promptTrksByPV[pvr_idx].size(); ++jtrk) {
						int idxj = promptTrksByPV[pvr_idx][jtrk];

						TVector3 pj(trk_px->at(idxj), trk_py->at(idxj), trk_pz->at(idxj));
						TVector3 xj(trk_x->at(idxj), trk_y->at(idxj), trk_z->at(idxj));
						TVector3 dxdyj = calcDxDy(xj, pj, pv);
						double dij = dxdyi.Dot(dxdyj)/(dxdyi.Mag()*dxdyj.Mag());
						double dij2 = dxdyi.Mag()/dxdyj.Mag();

						int nShared = getNSharedVeloHits(idx,idxj);

						if(nShared>=2) {
							if(nShared>=4) {
								double randMax = hphi2->GetMaximum();
								if( hphi2->GetBinContent(hphi2->FindBin(dij)) < randMax*r.Rndm() ) {
									//std::cout << "failed phi in >=4 bin" << dij <<"\t"<< hphi2->GetBinContent(hphi2->FindBin(dij)) <<"\t"<< randMax << std::endl;
									rerunEvt=true;
									break;
								}
								randMax = hratio2->GetMaximum();
								if( hratio2->GetBinContent(hratio2->FindBin(dij2)) < randMax*r.Rndm() ) {
									//std::cout << "failed ratio in >=4 bin" << dij2 <<"\t"<< hratio2->GetBinContent(hratio2->FindBin(dij2)) <<"\t"<< randMax << std::endl;
									rerunEvt=true;
									break;
								}
							} else {
								double randMax = hphi1->GetMaximum();
								if( hphi1->GetBinContent(hphi1->FindBin(dij)) < randMax*r.Rndm() ) {
									//std::cout << "failed phi in >=2 bin" << dij <<"\t"<< hphi1->GetBinContent(hphi1->FindBin(dij)) <<"\t"<< randMax << std::endl;
									rerunEvt=true;
									break;
								}
								randMax = hratio1->GetMaximum();
								if( hratio1->GetBinContent(hratio1->FindBin(dij2)) < randMax*r.Rndm() ) {
									//std::cout << "failed ratio in >=2 bin" << dij2 <<"\t"<< hratio1->GetBinContent(hratio1->FindBin(dij2)) <<"\t"<< randMax << std::endl;
									rerunEvt=true;
									break;
								}
							}
						}
					}
					if(rerunEvt) {
						break;
					}
					promptTrksByPV[pvr_idx].push_back(idx);
				}

				//now check track meets IP chi2 requirement (Tab 2 ANA-2014-074)
				//IP
				if(irep==0) {
					hist0.Fill(origtrk_ip->at(idx));
					hist1.Fill(TMath::Sqrt(dx*dx + dy*dy));
					double tmpX = origtrk_x->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*tx - pvr_x->at(pvr_idx);
					double tmpY = origtrk_y->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*ty - pvr_y->at(pvr_idx);
					hist1b.Fill(TMath::Sqrt(tmpX*tmpX + tmpY*tmpY));
				} else {
					hist2.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy));
				}
				//with chi2 cuts
				if(irep==0) {
					if(origtrk_ip_chi2->at(idx) < 16) hist8.Fill(origtrk_ip->at(idx));
					else hist11.Fill(origtrk_ip->at(idx));
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
					hist18.Fill(TMath::Sqrt(origtrk_ip_chi2->at(idx)),invpt);
					hist19.Fill(TMath::Sqrt(dx*dx + dy*dy)/sigma,invpt);

					//std::cout << TMath::Sqrt(origtrk_ip->at(idx)) << "\t" << TMath::Sqrt(dx*dx + dy*dy) << std::endl;
				} else {
					hist20.Fill(TMath::Sqrt(newDx*newDx + newDy*newDy)/sigma,invpt);
				}
				//IPchi2
				if(irep==0) {
					hist5.Fill(TMath::Log(origtrk_ip_chi2->at(idx)));
					hist6.Fill(TMath::Log((dx*dx + dy*dy)/sigma2));
					
					trueSumIPchi2 += origtrk_ip_chi2->at(idx);
					sumIPchi2 += (dx*dx + dy*dy)/sigma2;

					if(origtrk_ip_chi2->at(idx) < 13.) ++nTrueInIPchi2Bin[0];
//					else if(origtrk_ip_chi2->at(idx) < 12.) ++nTrueInIPchi2Bin[1];
//					else if(origtrk_ip_chi2->at(idx) < 18.) ++nTrueInIPchi2Bin[2];
//					else if(origtrk_ip_chi2->at(idx) < 24.) ++nTrueInIPchi2Bin[3];
					else ++nTrueInIPchi2Bin[4];
					if((dx*dx + dy*dy)/sigma2 < 13.) ++nInIPchi2Bin[0];
//					else if((dx*dx + dy*dy)/sigma2 < 12.) ++nInIPchi2Bin[1];
//					else if((dx*dx + dy*dy)/sigma2 < 18.) ++nInIPchi2Bin[2];
//					else if((dx*dx + dy*dy)/sigma2 < 24.) ++nInIPchi2Bin[3];
					else ++nInIPchi2Bin[4];
				} else {
					hist7.Fill(TMath::Log((newDx*newDx + newDy*newDy)/sigma2));

					sumIPchi2 += (newDx*newDx + newDy*newDy)/sigma2;

					if((newDx*newDx + newDy*newDy)/sigma2 < 13.) ++nInIPchi2Bin[0];
//					else if((newDx*newDx + newDy*newDy)/sigma2 < 12.) ++nInIPchi2Bin[1];
//					else if((newDx*newDx + newDy*newDy)/sigma2 < 18.) ++nInIPchi2Bin[2];
//					else if((newDx*newDx + newDy*newDy)/sigma2 < 24.) ++nInIPchi2Bin[3];
					else ++nInIPchi2Bin[4];
				}
				if(irep==0 && trk_idx_pvr->at(idx)!=-1 && trk_idx_gen->at(idx)!=-1 && gen_idx_pvr->at(trk_idx_gen->at(idx))!=-1) {
					hist21.Fill(pvr_x->at(trk_idx_pvr->at(idx)) - pvr_x->at(gen_idx_pvr->at(trk_idx_gen->at(idx))));
					hist22.Fill(pvr_y->at(trk_idx_pvr->at(idx)) - pvr_y->at(gen_idx_pvr->at(trk_idx_gen->at(idx))));
					hist23.Fill(pvr_z->at(trk_idx_pvr->at(idx)) - pvr_z->at(gen_idx_pvr->at(trk_idx_gen->at(idx))));
				}

				//apply the IPchi2 cut to the new IPchi2
				//if(irep==0) {
				//	if(origtrk_ip_chi2->at(idx) < 16.) continue;
				//} else {
					if((newDx*newDx + newDy*newDy) / sigma2 < 13.) continue;
				//}
				if(irep>0) ++countAfter;
				else ++countBefore;
				if(trksToWatch.count(idx)) std::cout << "DEBUG track " << idx << " PASS IPCHI2" << std::endl;

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
			}//end loop over tracks

			//check if and pairs of tracks failed correlation tests
			if(rerunEvt) {
				--irep;
				++rerunsRequired;
				++reruns;
				//if(rerunsRequired % 1000 == 0) std::cout << "Now rerun an event " << rerunsRequired << " times." << std::endl;
				continue;
			} else {
				if(reruns>49) {
					std::cout << reruns << std::endl;
					reruns=49;
				}

				++nReruns[reruns];
				reruns=0;
			}

			if(trks.size()>0) {
				passedTracksByRep[irep] = trks;
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

			if(irep==0) {
				hist24.Fill(trueSumIPchi2);
				hist25.Fill(sumIPchi2);

				histA0.Fill(nTrueInIPchi2Bin[0]);
				histA1.Fill(nInIPchi2Bin[0]);
				histB0.Fill(nTrueInIPchi2Bin[1]);
				histB1.Fill(nInIPchi2Bin[1]);
				histC0.Fill(nTrueInIPchi2Bin[2]);
				histC1.Fill(nInIPchi2Bin[2]);
				histD0.Fill(nTrueInIPchi2Bin[3]);
				histD1.Fill(nInIPchi2Bin[3]);
				histE0.Fill(nTrueInIPchi2Bin[4]);
				histE1.Fill(nInIPchi2Bin[4]);
			} else {
				hist26.Fill(sumIPchi2);
				//std::cout << sumIPchi2<<std::endl;

				histA2.Fill(nInIPchi2Bin[0]);
				histB2.Fill(nInIPchi2Bin[1]);
				histC2.Fill(nInIPchi2Bin[2]);
				histD2.Fill(nInIPchi2Bin[3]);
				histE2.Fill(nInIPchi2Bin[4]);
			}

			//check we have at least two tracks to vertex
			if(ntrk==0) {
				if(irep==0) ++noTrks1;
				else ++noTrks2;
			} else if(ntrk==1) {
				if(irep==0) ++oneTrk1;
				else ++oneTrk2;
			} else if(ntrk>1) {
				if(irep==0) ++twoTrks1;
				else ++twoTrks2;

				//now make 2body SVs
				std::vector<std::pair<int,int> > sv2ij;
				std::vector<TVector3> sv2;
				std::vector<TVector3> sv2fv;

				for( int itrk=0; itrk<ntrk; ++itrk) {
					TLorentzVector p4i(trk_px->at(trks[itrk]), trk_py->at(trks[itrk]), trk_pz->at(trks[itrk]), 0);
					p4i.SetE(TMath::Sqrt(p4i.P()*p4i.P() + 140*140));
					for( int jtrk=itrk+1; jtrk<ntrk; ++jtrk) {
						if(trk_idx_gen->at(trks[itrk]) == trk_idx_gen->at(trks[jtrk])) {
							if(irep==0) ++nMatchGen1;
							else ++nMatchGen2;
						}
						//if(trk_vid->at(trks[itrk]) == trk_vid->at(trks[jtrk])) std::cout << "shared VELO track" << std::endl;
						int nSharedHits = getNSharedEarlyVeloHits(trks[itrk],trks[jtrk]);
						if(irep==0) vhist11a.Fill(nSharedHits);
						else vhist11b.Fill(nSharedHits);
						
						TLorentzVector p4j(trk_px->at(trks[jtrk]), trk_py->at(trks[jtrk]), trk_pz->at(trks[jtrk]), 0);
						p4j.SetE(TMath::Sqrt(p4j.P()*p4j.P() + 140*140));
						if(irep==0) {
							vhist2a.Fill(p4i.Angle(p4j.Vect()));
						} else {
							vhist2b.Fill(p4i.Angle(p4j.Vect()));
						}
						TLorentzVector p4 = p4i+p4j;
						double m = p4.M();
						TVector3 sv;
						double d = calcDoca(sv, trkhit[itrk], trkdir[itrk], trkhit[jtrk], trkdir[jtrk]);
						
						if(irep==0) {
							vhist2e.Fill(p4i.Angle(p4j.Vect()), sv.Perp());
							vhist2i.Fill(p4i.Angle(p4j.Vect()), p4.M());
						} else {
							vhist2f.Fill(p4i.Angle(p4j.Vect()), sv.Perp());
							vhist2j.Fill(p4i.Angle(p4j.Vect()), p4.M());
						}

						//find best pvr for our svr
						pvr_idx = findBestPvrForSvr(sv, p4.Vect(), false);
						if(pvr_idx<0) {
							std::cout << "WARNING - shouldn't reach here! (2)" << std::endl;
							continue;
						}

						TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
								sv.Y() - pvr_y->at(pvr_idx),
								sv.Z() - pvr_z->at(pvr_idx));
						double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();
						double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);
						double tz = fv.Z()*m/p4.Z()/(3e11)*(1e12);
						double sigma2i = origtrk_ip->at(trks[itrk])*origtrk_ip->at(trks[itrk]) / origtrk_ip_chi2->at(trks[itrk]);
						double sigma2j = origtrk_ip->at(trks[jtrk])*origtrk_ip->at(trks[jtrk]) / origtrk_ip_chi2->at(trks[jtrk]);
						double vtxchi2 = d*d/(sigma2i+sigma2j);

						TVector3 dxdyi = calcDxDy(trkhit[itrk], trkdir[itrk], TVector3(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx)));
						TVector3 dxdyj = calcDxDy(trkhit[jtrk], trkdir[jtrk], TVector3(pvr_x->at(pvr_idx),pvr_y->at(pvr_idx),pvr_z->at(pvr_idx)));

						if(irep==0) {
							vhist0a.Fill(dxdyi.X() - dxdyj.X());
							vhist1a.Fill(dxdyi.Y() - dxdyj.Y());
							vhist4a.Fill(tz);
							vhist5a.Fill(fv.Z());
							vhist6a.Fill(m);
							vhist7a.Fill(p4.Z());
							vhist8a.Fill(mcor);
							vhist9a.Fill(d);
							TVector3 xi(trk_x->at(trks[itrk]), trk_y->at(trks[itrk]), trk_z->at(trks[itrk]));
							TVector3 xj(trk_x->at(trks[jtrk]), trk_y->at(trks[jtrk]), trk_z->at(trks[jtrk]));
							int idx_pvi = findBestGenPvrForTrk(trks[itrk], origtrk_x->at(trks[itrk]), origtrk_y->at(trks[itrk]));
							int idx_pvj = findBestGenPvrForTrk(trks[jtrk], origtrk_x->at(trks[jtrk]), origtrk_y->at(trks[itrk]));
							TVector3 pvi(pvr_x->at(idx_pvi), pvr_y->at(idx_pvi), pvr_z->at(idx_pvi));
							TVector3 pvj(pvr_x->at(idx_pvj), pvr_y->at(idx_pvj), pvr_z->at(idx_pvj));
							TVector3 dxdyi = calcDxDy(xi,p4i.Vect(),pvi);
							TVector3 dxdyj = calcDxDy(xj,p4j.Vect(),pvj);
							vhist12a.Fill(fv.Z(),(dxdyi-dxdyj).Mag());
						} else {
							vhist0b.Fill(dxdyi.X() - dxdyj.X());
							vhist1b.Fill(dxdyi.Y() - dxdyj.Y());
							vhist4b.Fill(tz);
							vhist5b.Fill(fv.Z());
							vhist6b.Fill(m);
							vhist7b.Fill(p4.Z());
							vhist8b.Fill(mcor);
							vhist9b.Fill(d);
							TVector3 xi(trk_x->at(trks[itrk]), trk_y->at(trks[itrk]), trk_z->at(trks[itrk]));
							TVector3 xj(trk_x->at(trks[jtrk]), trk_y->at(trks[jtrk]), trk_z->at(trks[jtrk]));
							int idx_pvi = findBestGenPvrForTrk(trks[itrk], origtrk_x->at(trks[itrk]), origtrk_y->at(trks[itrk]));
							int idx_pvj = findBestGenPvrForTrk(trks[jtrk], origtrk_x->at(trks[jtrk]), origtrk_y->at(trks[itrk]));
							TVector3 pvi(pvr_x->at(idx_pvi),pvr_y->at(idx_pvi),pvr_z->at(idx_pvi));
							TVector3 pvj(pvr_x->at(idx_pvj),pvr_y->at(idx_pvj),pvr_z->at(idx_pvj));
							TVector3 dxdyi = calcDxDy(xi,p4i.Vect(),pvi);
							TVector3 dxdyj = calcDxDy(xj,p4j.Vect(),pvj);
							vhist12b.Fill(fv.Z(),(dxdyi-dxdyj).Mag());
						}

						//			std::cout << fv.Z() << "\t" << m << "\t" << p4.Z() << "\t" << tz << std::endl;
						//Apply 2-body SV requirements
						if(irep==0) ++n2BSVs10;
						else ++n2BSVs20;
						if(m < 0. || m > 5000.) continue;
						if(irep==0) ++n2BSVs11;
						else ++n2BSVs21;
						if(d > 0.2) continue;
						if(irep==0) ++n2BSVs12;
						else ++n2BSVs22;
						if(mcor<600.) continue;
						if(irep==0) ++n2BSVs13;
						else ++n2BSVs23;
						if(tz<0. || tz>=10.) continue;
						if(irep==0) ++n2BSVs14;
						else ++n2BSVs24;
						if(vtxchi2>=10) continue;
						//TODO add VTXCHI2<10, FDCHI2>25, HITS>0 
						sv2ij.push_back(std::pair<int,int>(trks[itrk],trks[jtrk]));
						sv2.push_back(sv);
						sv2fv.push_back(fv);
						if(irep==0) {
							vhist0c.Fill(dxdyi.X() - dxdyj.X());//TMath::Abs(trkhit[itrk].X() - trkhit[jtrk].X()));
							vhist1c.Fill(dxdyi.Y() - dxdyj.Y());//TMath::Abs(trkhit[itrk].Y() - trkhit[jtrk].Y()));
							vhist2c.Fill(p4i.Angle(p4j.Vect()));
							vhist2g.Fill(p4i.Angle(p4j.Vect()), sv.Perp());
							//vhist2k.Fill(p4i.Angle(p4j.Vect()), p4.M());
							vhist4c.Fill(tz);
							vhist5c.Fill(fv.Z());
							vhist6c.Fill(m);
							vhist7c.Fill(p4.Z());
							vhist8c.Fill(mcor);
							vhist9c.Fill(d);
						} else {
							vhist0d.Fill(dxdyi.X() - dxdyj.X());//TMath::Abs(trkhit[itrk].X() - trkhit[jtrk].X()));
							vhist1d.Fill(dxdyi.Y() - dxdyj.Y());//TMath::Abs(trkhit[itrk].Y() - trkhit[jtrk].Y()));
							vhist2d.Fill(p4i.Angle(p4j.Vect()));
							vhist2h.Fill(p4i.Angle(p4j.Vect()), sv.Perp());
							//vhist2l.Fill(p4i.Angle(p4j.Vect()), p4.M());
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

				for(uint ijet=0; ijet<jet_pz->size(); ++ijet) {
					std::set<int> used; // only use each SV once

				//check we have at least one 2body SV
					//more than one SV - now link
					for(uint s=0; s<sv2.size(); ++s) {
						//first check if we've used this SV already in this jet
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
						int bestpvr(-1);
						double minDOCA(9999.);
						for(uint ipvr=0; ipvr<pvr_z->size(); ++ipvr) {
							TVector3 pv(pvr_x->at(ipvr),pvr_y->at(ipvr),pvr_z->at(ipvr));
							double doca = calcDocaPoint(sv,p4.Vect(),pv);
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

						double m = p4.M();

						TVector3 fv(sv.X() - pvr_x->at(pvr_idx),
								sv.Y() - pvr_y->at(pvr_idx),
								sv.Z() - pvr_z->at(pvr_idx));
						double pt2 = p4.Vect().Cross(fv.Unit()).Mag2();//note this is transverse wrt fd NOT the z-axis
						double mcor = TMath::Sqrt(m*m + pt2) + TMath::Sqrt(pt2);

						if(irep==0) ++nLinkedSVs1;
						else ++nLinkedSVs2;

						//count tracks in jet
						int nTrkInJet(0);
						for(uint itrk=0; itrk<svTrks.size(); ++itrk) {
//							if(findTrkInJet(itrk, jet_idx) != -1) ++nTrkInJet;
							double deltaR = p4jet.DeltaR(TLorentzVector(trk_px->at(svTrks[itrk]), trk_py->at(svTrks[itrk]), trk_pz->at(svTrks[itrk]), 0));
							if(deltaR<0.5) ++nTrkInJet;
						}

						if(irep==0) vhist13a.Fill(nTrkInJet);
						else vhist13b.Fill(nTrkInJet);

						bool pass(true);

						//check if we pass n-body requirements
						if(svTrks.size() == 2 && TMath::Abs(m - 500.) < 20.) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs11;
							else ++nLinkedSVs21;
						}
						if(fd_min > 15.) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs12;
							else ++nLinkedSVs22;
						}
						if(sv.Z() > 200.) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs13;
							else ++nLinkedSVs23;
						}
						if(mcor < 600.) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs14;
							else ++nLinkedSVs24;
						}
						if(nTrkInJet<2) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs15;
							else ++nLinkedSVs25;
						}
						if(p4.Perp() < 2000.) pass=false;
						if(pass) {
							if(irep==0) ++nLinkedSVs16;
							else ++nLinkedSVs26;
						}
						//else if (!(m_tau < 1.5))                  pass = false;
						//else if (!(m_fdChi2 > 32))                pass = false;
						//else if (!(m_drSvrJet <= m_parent->m_dr)) pass = false;
						//if (m_parent->m_nbvSelect && !pass) return false;

						//save vertex and list of tracks
						svN.push_back(sv);
						svNij.push_back(svTrks);
						//count SVs
						if(svTrks.size()==2) {
							if(irep==0) ++sv2sAfter1;
							else ++sv2sAfter2;
							if(irep==0) vhist2k.Fill(p4ij[0].Angle(p4ij[1].Vect()), p4.M());
							else vhist2l.Fill(p4ij[0].Angle(p4ij[1].Vect()), p4.M());
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

						//pad with -1 for filling tuple - size of svTrks no longer meaningful
						while(svTrks.size()<10) {
							svTrks.push_back(-1);
						}

						svr_idx_pvr->push_back(trk_idx_pvr->at(svTrks[0]));
						svr_idx_jet->push_back(ijet);
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
						svr_pass->push_back(pass);
						//std::cout << jentry << std::endl;

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
				//}

				}

				if(irep==0) svNsAfter1+=svN.size();
				else svNsAfter2+=svN.size();
				
			}

			for(uint isvr=0; isvr<svr_z->size(); ++isvr) {
				pvr_idx = svr_idx_pvr->at(isvr);

				if(irep==0) {
					//if(svr_idx_trk0->at(isvr)!=-1 && TMath::Log(origtrk_ip_chi2->at(svr_idx_trk0->at(isvr)))>8.) std::cout << TMath::Log(origtrk_ip_chi2->at(svr_idx_trk0->at(isvr))) << "\t" << TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr))) << std::endl;
					//if(svr_idx_trk1->at(isvr)!=-1 && TMath::Log(origtrk_ip_chi2->at(svr_idx_trk1->at(isvr)))>8.) std::cout << TMath::Log(origtrk_ip_chi2->at(svr_idx_trk1->at(isvr))) << "\t" << TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr))) << std::endl;

					if(svr_idx_trk0->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr))));
					if(svr_idx_trk1->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr))));
					if(svr_idx_trk2->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk2->at(isvr))));
					if(svr_idx_trk3->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk3->at(isvr))));
					if(svr_idx_trk4->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk4->at(isvr))));
					if(svr_idx_trk5->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk5->at(isvr))));
					if(svr_idx_trk6->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk6->at(isvr))));
					if(svr_idx_trk7->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk7->at(isvr))));
					if(svr_idx_trk8->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk8->at(isvr))));
					if(svr_idx_trk9->at(isvr)!=-1) hist6b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk9->at(isvr))));
				} else {
					if(svr_idx_trk0->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk0->at(isvr))));
					if(svr_idx_trk1->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk1->at(isvr))));
					if(svr_idx_trk2->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk2->at(isvr))));
					if(svr_idx_trk3->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk3->at(isvr))));
					if(svr_idx_trk4->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk4->at(isvr))));
					if(svr_idx_trk5->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk5->at(isvr))));
					if(svr_idx_trk6->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk6->at(isvr))));
					if(svr_idx_trk7->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk7->at(isvr))));
					if(svr_idx_trk8->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk8->at(isvr))));
					if(svr_idx_trk9->at(isvr)!=-1) hist7b.Fill(TMath::Log(trk_ip_chi2->at(svr_idx_trk9->at(isvr))));
				}
			}


			if(newtree) newtree->Fill();
		}//end loop over reps

		//if(passedTracksByRep.size()>0 && passedTracksByRep.find(0) != passedTracksByRep.end() && passedTracksByRep[0].size()>1) {
		//	std::cout << jentry << std::endl;
		//	std::map<int,std::vector<int>>::iterator repIt = passedTracksByRep.begin();
		//	for( ; repIt!= passedTracksByRep.end(); ++repIt) {
		//		std::cout << (*repIt).first << ":\t";
		//		std::vector<int>::iterator trkIt = (*repIt).second.begin();
		//		for( ; trkIt!= (*repIt).second.end(); ++trkIt) {
		//			std::cout << (*trkIt) << "\t";
		//		}
		//		std::cout << std::endl;
		//	}
		//}

		origtrk_x->clear();
		origtrk_y->clear();

		delete origtrk_x;
		delete origtrk_y;

	}//end loop over events

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
//	std::cout << "Nisv\t" << "---" << "\t" << static_cast<double>(svsAfter1)/1. << "\t" << static_cast<double>(svsAfter2)/100. << std::endl;
	std::cout << "Nisv\t" << "---" << "\t" << static_cast<double>(n2BSVs10)/1. << "\t" << static_cast<double>(n2BSVs20)/100. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(n2BSVs11)/static_cast<double>(n2BSVs10)/1. << "\t" << static_cast<double>(n2BSVs21)/static_cast<double>(n2BSVs20)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(n2BSVs12)/static_cast<double>(n2BSVs11)/1. << "\t" << static_cast<double>(n2BSVs22)/static_cast<double>(n2BSVs21)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(n2BSVs13)/static_cast<double>(n2BSVs12)/1. << "\t" << static_cast<double>(n2BSVs23)/static_cast<double>(n2BSVs22)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(n2BSVs14)/static_cast<double>(n2BSVs13)/1. << "\t" << static_cast<double>(n2BSVs24)/static_cast<double>(n2BSVs23)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(svsAfter1)/static_cast<double>(n2BSVs14)/1. << "\t" << static_cast<double>(svsAfter2)/static_cast<double>(n2BSVs24)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(svsAfter1)/1. << "\t" << static_cast<double>(svsAfter2)/100. << std::endl;
	std::cout << std::endl;
	std::cout << "Nlinked\t" << "---" << "\t" << static_cast<double>(nLinkedSVs1)/1. << "\t" << static_cast<double>(nLinkedSVs2)/100. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs11)/static_cast<double>(nLinkedSVs1)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs21)/static_cast<double>(nLinkedSVs2)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs12)/static_cast<double>(nLinkedSVs11)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs22)/static_cast<double>(nLinkedSVs21)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs13)/static_cast<double>(nLinkedSVs12)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs23)/static_cast<double>(nLinkedSVs22)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs14)/static_cast<double>(nLinkedSVs13)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs24)/static_cast<double>(nLinkedSVs23)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs15)/static_cast<double>(nLinkedSVs14)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs25)/static_cast<double>(nLinkedSVs24)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs16)/static_cast<double>(nLinkedSVs15)/1. << "\t" 
		                                  << static_cast<double>(nLinkedSVs26)/static_cast<double>(nLinkedSVs25)/1. << std::endl;
	std::cout << ".......\t" << "---" << "\t" << static_cast<double>(nLinkedSVs16)/1. << "\t" << static_cast<double>(nLinkedSVs26)/100. << std::endl;
	std::cout << std::endl;
	std::cout << "Nsv2\t" << "---" << "\t" << static_cast<double>(sv2sAfter1)/1. << "\t" << static_cast<double>(sv2sAfter2)/100. << std::endl;
	std::cout << "Nsv3\t" << "---" << "\t" << static_cast<double>(sv3sAfter1)/1. << "\t" << static_cast<double>(sv3sAfter2)/100. << std::endl;
	std::cout << "Nsv4\t" << "---" << "\t" << static_cast<double>(sv4sAfter1)/1. << "\t" << static_cast<double>(sv4sAfter2)/100. << std::endl;
	std::cout << "Nsv5+\t" << "---" << "\t" << static_cast<double>(svMoresAfter1)/1. << "\t" << static_cast<double>(svMoresAfter2)/100. << std::endl;
	std::cout << "NsvN\t" << svsBefore << "\t" << static_cast<double>(svNsAfter1)/1. << "\t" << static_cast<double>(svNsAfter2)/100. << std::endl;
	std::cout << "Ntrk\t" << sumntrk1 << "\t" << static_cast<double>(sumntrk2)/100. << std::endl;
	std::cout << "Ncomb\t" << sumncomb1 << "\t" << static_cast<double>(sumncomb2)/100. << std::endl;
	std::cout << "N(same gen)\t" << nMatchGen1 << "\t" << static_cast<double>(nMatchGen2)/100. << std::endl;
	std::cout << "N(0 trks)\t" << noTrks1 << "\t" << static_cast<double>(noTrks2)/100. << std::endl;
	std::cout << "N(1 trk)\t" << oneTrk1 << "\t" << static_cast<double>(oneTrk2)/100. << std::endl;
	std::cout << "N(>2 trks)\t" << twoTrks1 << "\t" << static_cast<double>(twoTrks2)/100. << std::endl;
	std::cout << "Nisv/Ncomb\t" << static_cast<double>(svsAfter1)/sumncomb1 << " +/- " << TMath::Sqrt(svsAfter1)/sumncomb1 << "\t" << static_cast<double>(svsAfter2)/sumncomb2 << " +/- " << TMath::Sqrt(svsAfter2)/sumncomb2 << std::endl;
	std::cout << "prompt:" << static_cast<double>(total)/(total+nonPrompt) << std::endl;
	std::cout << "no truth:" << noTruth << "/" << realTotal << std::endl;
	std::cout << "reruns:" << rerunsRequired << std::endl;
	for(int i=0; i<50; ++i) {
		std::cout << nReruns[i] << "\t";
	}
	std::cout << std::endl;

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
	hist5b.Scale(1.   );//,"width");
	hist6b.Scale(1    );//,"width");
	hist7b.Scale(0.01 );//,"width");
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
	hist26.Scale(0.01);//,"width");
	histA2.Scale(0.01);//,"width");
	histB2.Scale(0.01);//,"width");
	histC2.Scale(0.01);//,"width");
	histD2.Scale(0.01);//,"width");
	histE2.Scale(0.01);//,"width");

	hist0.GetXaxis()->SetTitle("IP");
	hist3.GetXaxis()->SetTitle("IP");
	hist5.GetXaxis()->SetTitle("#chi^{2}_{IP}");
	hist8.GetXaxis()->SetTitle("IP (#chi^{2}_{IP} < 16)");
	hist11.GetXaxis()->SetTitle("IP (#chi^{2}_{IP} > 16)");
	hist14.GetXaxis()->SetTitle("IP_{X}");
	hist16.GetXaxis()->SetTitle("IP_{Y}");
	hist24.GetXaxis()->SetTitle("#Sigma#chi^{2}_{IP}");

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
	hist2.SetLineColor(kGreen+2);
	hist2.Draw("same");
	c.SaveAs("ip.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	if(hist9.GetMaximum()>hist8.GetMaximum()) hist8.SetMaximum(1.1*hist9.GetMaximum());
	if(hist10.GetMaximum()>hist8.GetMaximum()) hist8.SetMaximum(1.1*hist10.GetMaximum());
	hist8.Draw("EP");
	hist9.SetLineColor(kRed);
	hist9.Draw("same");
	hist10.SetLineColor(kGreen+2);
	hist10.Draw("same");
	c.SaveAs("ip_ipchi2LT16.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	if(hist12.GetMaximum()>hist11.GetMaximum()) hist11.SetMaximum(1.1*hist12.GetMaximum());
	if(hist13.GetMaximum()>hist11.GetMaximum()) hist11.SetMaximum(1.1*hist13.GetMaximum());
	hist11.Draw("EP");
	hist12.SetLineColor(kRed);
	hist12.Draw("same");
	hist13.SetLineColor(kGreen+2);
	hist13.Draw("same");
	c.SaveAs("ip_ipchi2GT16.pdf");

	c.SetLogx(0);
	c.SetLogy(0);
	if(hist15.GetMaximum()>hist14.GetMaximum()) hist14.SetMaximum(1.1*hist15.GetMaximum());
	hist14.SetLineColor(kRed);
	hist14.Draw();
	hist15.SetLineColor(kGreen+2);
	hist15.Draw("same");
	c.SaveAs("ipx.pdf");

	c.SetLogx(0);
	c.SetLogy(0);
	if(hist17.GetMaximum()>hist16.GetMaximum()) hist16.SetMaximum(1.1*hist17.GetMaximum());
	hist16.SetLineColor(kRed);
	hist16.Draw();
	hist17.SetLineColor(kGreen+2);
	hist17.Draw("same");
	c.SaveAs("ipy.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(hist6.GetMaximum()>hist5.GetMaximum()) hist5.SetMaximum(1.1*hist6.GetMaximum());
	if(hist7.GetMaximum()>hist5.GetMaximum()) hist5.SetMaximum(1.1*hist7.GetMaximum());
	hist5.Draw("EP");
	hist6.SetLineColor(kRed);
	hist6.Draw("same");
	hist7.SetLineColor(kGreen+2);
	hist7.Draw("same");
	hist5b.SetLineStyle(kDashed);
	hist5b.Draw("same");
	hist6b.SetLineStyle(kDashed);
	hist6b.SetLineColor(kRed);
	hist6b.Draw("same");
	hist7b.SetLineStyle(kDashed);
	hist7b.SetLineColor(kGreen+2);
	hist7b.Draw("same");
	c.SaveAs("ipchi2.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(hist25.GetMaximum()>hist24.GetMaximum()) hist24.SetMaximum(1.1*hist25.GetMaximum());
	if(hist26.GetMaximum()>hist24.GetMaximum()) hist24.SetMaximum(1.1*hist26.GetMaximum());
	hist24.Draw("EP");
	hist25.SetLineColor(kRed);
	hist25.Draw("same");
	hist26.SetLineColor(kGreen+2);
	hist26.Draw("same");
	c.SaveAs("sumipchi2.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(histA1.GetMaximum()>histA0.GetMaximum()) histA0.SetMaximum(1.1*histA1.GetMaximum());
	if(histA2.GetMaximum()>histA0.GetMaximum()) histA0.SetMaximum(1.1*histA2.GetMaximum());
	histA0.Draw("EP");
	histA1.SetLineColor(kRed);
	histA1.Draw("same");
	histA2.SetLineColor(kGreen+2);
	histA2.Draw("same");
	c.SaveAs("NinIPchi2BinA.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(histB1.GetMaximum()>histB0.GetMaximum()) histB0.SetMaximum(1.1*histB1.GetMaximum());
	if(histB2.GetMaximum()>histB0.GetMaximum()) histB0.SetMaximum(1.1*histB2.GetMaximum());
	histB0.Draw("EP");
	histB1.SetLineColor(kRed);
	histB1.Draw("same");
	histB2.SetLineColor(kGreen+2);
	histB2.Draw("same");
	c.SaveAs("NinIPchi2BinB.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(histC1.GetMaximum()>histC0.GetMaximum()) histC0.SetMaximum(1.1*histC1.GetMaximum());
	if(histC2.GetMaximum()>histC0.GetMaximum()) histC0.SetMaximum(1.1*histC2.GetMaximum());
	histC0.Draw("EP");
	histC1.SetLineColor(kRed);
	histC1.Draw("same");
	histC2.SetLineColor(kGreen+2);
	histC2.Draw("same");
	c.SaveAs("NinIPchi2BinC.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(histD1.GetMaximum()>histD0.GetMaximum()) histD0.SetMaximum(1.1*histD1.GetMaximum());
	if(histD2.GetMaximum()>histD0.GetMaximum()) histD0.SetMaximum(1.1*histD2.GetMaximum());
	histD0.Draw("EP");
	histD1.SetLineColor(kRed);
	histD1.Draw("same");
	histD2.SetLineColor(kGreen+2);
	histD2.Draw("same");
	c.SaveAs("NinIPchi2BinD.pdf");

	c.SetLogx(0);
	c.SetLogy(1);
	if(histE1.GetMaximum()>histE0.GetMaximum()) histE0.SetMaximum(1.1*histE1.GetMaximum());
	if(histE2.GetMaximum()>histE0.GetMaximum()) histE0.SetMaximum(1.1*histE2.GetMaximum());
	histE0.Draw("EP");
	histE1.SetLineColor(kRed);
	histE1.Draw("same");
	histE2.SetLineColor(kGreen+2);
	histE2.Draw("same");
	c.SaveAs("NinIPchi2BinE.pdf");

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
	htmp->SetLineColor(kGreen+2); htmp->Draw("same");
	htmp = hist4.ProjectionX("hist2_3",3,3);
	htmp->SetLineColor(kGreen+2); htmp->SetLineStyle(kDashed); htmp->Draw("same");
	c.SaveAs("ip_sigma2.pdf");

	c.SetLogx(1);
	c.SetLogy(0);
	htmp = hist18.ProjectionX("hist18_1",1,1);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_1",1,1);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_1",1,1);
	htmp->SetLineColor(kGreen+2); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin1.pdf");
	htmp = hist18.ProjectionX("hist18_2",2,2);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_2",2,2);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_2",2,2);
	htmp->SetLineColor(kGreen+2); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin2.pdf");
	htmp = hist18.ProjectionX("hist18_3",3,3);
	htmp->SetLineColor(kBlue); htmp->Draw("EP");
	htmp = hist19.ProjectionX("hist19_3",3,3);
	htmp->SetLineColor(kRed); htmp->Draw("same");
	htmp = hist20.ProjectionX("hist20_3",3,3);
	htmp->SetLineColor(kGreen+2); htmp->Draw("same");
	c.SaveAs("ipOverSigma_bin3.pdf");

	//vhist0c.Scale(1./vhist0a.Integral());
	//vhist1c.Scale(1./vhist1a.Integral());
	vhist0d.Scale(0.01);//1./vhist0b.Integral());
	vhist1d.Scale(0.01);//1./vhist1b.Integral());
	vhist2d.Scale(0.01);//1./vhist2b.Integral());
	vhist2f.Scale(0.01);//1./vhist2f.Integral());
	vhist2h.Scale(0.01);//1./vhist2f.Integral());
	vhist2j.Scale(0.01);//1./vhist2j.Integral());
	vhist2l.Scale(0.01);//1./vhist2j.Integral());
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

	//vhist0a.Scale(1./vhist0a.Integral());
	//vhist1a.Scale(1./vhist1a.Integral());
	vhist3a.Scale(1./vhist3a.Integral());
	vhist4a.Scale(1./vhist4a.Integral());
	vhist5z.Scale(1./vhist5z.Integral());
	vhist5a.Scale(1./vhist5a.Integral());
	vhist6a.Scale(1./vhist6a.Integral());
	vhist7a.Scale(1./vhist7a.Integral());
	vhist8a.Scale(1./vhist8a.Integral());
	vhist9a.Scale(1./vhist9a.Integral());

	vhist0b.Scale(0.01);//1./vhist0b.Integral());
	vhist1b.Scale(0.01);//1./vhist1b.Integral());
	vhist2b.Scale(0.01);//1./vhist2b.Integral());
	vhist3b.Scale(1./vhist3b.Integral());
	vhist4b.Scale(1./vhist4b.Integral());
	vhist5b.Scale(1./vhist5b.Integral());
	vhist6b.Scale(1./vhist6b.Integral());
	vhist7b.Scale(1./vhist7b.Integral());
	vhist8b.Scale(1./vhist8b.Integral());
	vhist9b.Scale(1./vhist9b.Integral());
	vhist10b.Scale(0.01);
	vhist11b.Scale(0.01);
	vhist13b.Scale(0.01);

	vhist0a.GetXaxis()->SetTitle("#DeltaX");
	vhist1a.GetXaxis()->SetTitle("#DeltaY");
	vhist2a.GetXaxis()->SetTitle("angle");
	vhist3a.GetXaxis()->SetTitle("N_{trk}");
	vhist4a.GetXaxis()->SetTitle("t_{z}");
	vhist5a.GetXaxis()->SetTitle("fv_{z}");
	vhist6a.GetXaxis()->SetTitle("m");
	vhist7a.GetXaxis()->SetTitle("p_{z}");
	vhist8a.GetXaxis()->SetTitle("m_{cor}");
	vhist9a.GetXaxis()->SetTitle("DOCA");
	vhist10a.GetXaxis()->SetTitle("nSV");
	vhist11a.GetXaxis()->SetTitle("shared hits");
	vhist13a.GetXaxis()->SetTitle("N in jet");

	c.SetLogx(0);
	c.SetLogy(0);
	vhist0a.SetLineColor(kRed);
	vhist0a.Draw();
	vhist0b.SetLineColor(kGreen+2);
	vhist0b.Draw("same");
	vhist0c.SetLineColor(kRed);
	vhist0c.SetLineStyle(kDashed);
	vhist0c.Draw("same");
	vhist0d.SetLineColor(kGreen+2);
	vhist0d.SetLineStyle(kDashed);
	vhist0d.Draw("same");
	c.SaveAs("deltaX.pdf");

	vhist1a.SetLineColor(kRed);
	vhist1a.Draw();
	vhist1b.SetLineColor(kGreen+2);
	vhist1b.Draw("same");
	vhist1c.SetLineColor(kRed);
	vhist1c.SetLineStyle(kDashed);
	vhist1c.Draw("same");
	vhist1d.SetLineColor(kGreen+2);
	vhist1d.SetLineStyle(kDashed);
	vhist1d.Draw("same");
	c.SaveAs("deltaY.pdf");

	vhist2a.SetLineColor(kRed);
	vhist2a.Draw();
	vhist2b.SetLineColor(kGreen+2);
	vhist2b.Draw("same");
	vhist2c.SetLineColor(kRed);
	vhist2c.SetLineStyle(kDashed);
	vhist2c.Draw("same");
	vhist2d.SetLineColor(kGreen+2);
	vhist2d.SetLineStyle(kDashed);
	vhist2d.Draw("same");
	c.SaveAs("angle.pdf");

	vhist3a.SetLineColor(kRed);
	vhist3a.Draw();
	vhist3b.SetLineColor(kGreen+2);
	vhist3b.Draw("same");
	c.SaveAs("ntrk.pdf");

	vhist4a.SetLineColor(kRed);
	vhist4a.Draw();
	vhist4b.SetLineColor(kGreen+2);
	vhist4b.Draw("same");
	vhist4c.SetLineColor(kRed);
	vhist4c.SetLineStyle(kDashed);
	vhist4c.Draw("same");
	vhist4d.SetLineColor(kGreen+2);
	vhist4d.SetLineStyle(kDashed);
	vhist4d.Draw("same");
	c.SaveAs("tz.pdf");

	if(vhist5b.GetMaximum()>vhist5a.GetMaximum()) vhist5a.SetMaximum(1.1*vhist5b.GetMaximum());
	if(vhist5z.GetMaximum()>vhist5a.GetMaximum()) vhist5a.SetMaximum(1.1*vhist5z.GetMaximum());
	vhist5a.SetLineColor(kRed);
	vhist5a.Draw();
	vhist5b.SetLineColor(kGreen+2);
	vhist5b.Draw("same");
	vhist5z.SetLineColor(kBlue);
	vhist5z.Draw("same");
	vhist5c.SetLineColor(kRed);
	vhist5c.SetLineStyle(kDashed);
	vhist5c.Draw("same");
	vhist5d.SetLineColor(kGreen+2);
	vhist5d.SetLineStyle(kDashed);
	vhist5d.Draw("same");
	c.SaveAs("fvz.pdf");

	vhist6a.SetLineColor(kRed);
	vhist6a.Draw();
	vhist6b.SetLineColor(kGreen+2);
	vhist6b.Draw("same");
	vhist6c.SetLineColor(kRed);
	vhist6c.SetLineStyle(kDashed);
	vhist6c.Draw("same");
	vhist6d.SetLineColor(kGreen+2);
	vhist6d.SetLineStyle(kDashed);
	vhist6d.Draw("same");
	c.SaveAs("m.pdf");

	vhist7a.SetLineColor(kRed);
	vhist7a.Draw();
	vhist7b.SetLineColor(kGreen+2);
	vhist7b.Draw("same");
	vhist7c.SetLineColor(kRed);
	vhist7c.SetLineStyle(kDashed);
	vhist7c.Draw("same");
	vhist7d.SetLineColor(kGreen+2);
	vhist7d.SetLineStyle(kDashed);
	vhist7d.Draw("same");
	c.SaveAs("pz.pdf");

	vhist8a.SetLineColor(kRed);
	vhist8a.Draw();
	vhist8b.SetLineColor(kGreen+2);
	vhist8b.Draw("same");
	vhist8c.SetLineColor(kRed);
	vhist8c.SetLineStyle(kDashed);
	vhist8c.Draw("same");
	vhist8d.SetLineColor(kGreen+2);
	vhist8d.SetLineStyle(kDashed);
	vhist8d.Draw("same");
	c.SaveAs("mcor.pdf");

	vhist9a.SetLineColor(kRed);
	vhist9a.Draw();
	vhist9b.SetLineColor(kGreen+2);
	vhist9b.Draw("same");
	vhist9c.SetLineColor(kRed);
	vhist9c.SetLineStyle(kDashed);
	vhist9c.Draw("same");
	vhist9d.SetLineColor(kGreen+2);
	vhist9d.SetLineStyle(kDashed);
	vhist9d.Draw("same");
	c.SaveAs("doca.pdf");

	vhist10a.SetLineColor(kRed);
	vhist10a.Draw();
	vhist10b.SetLineColor(kGreen+2);
	vhist10b.Draw("same");
	c.SaveAs("nsvr.pdf");

	vhist11a.SetLineColor(kRed);
	vhist11a.Draw();
	vhist11b.SetLineColor(kGreen+2);
	vhist11b.Draw("same");
	c.SaveAs("sharedhits.pdf");

	vhist13a.SetLineColor(kRed);
	vhist13a.Draw();
	vhist13b.SetLineColor(kGreen+2);
	vhist13b.Draw("same");
	c.SaveAs("ninjet.pdf");

	hist21.Draw();
	c.SaveAs("pvShiftX.pdf");
	hist22.Draw();
	c.SaveAs("pvShiftY.pdf");
	hist23.Draw();
	c.SaveAs("pvShiftZ.pdf");

	TCanvas c2;
	c2.Divide(2,2);
	c2.cd(1);
	vhist2e.Draw("colztext89");
	c2.cd(2);
	vhist2f.Draw("colztext89");
	c2.cd(3);
	vhist2g.Draw("colztext89");
	c2.cd(4);
	vhist2h.Draw("colztext89");
	c2.SaveAs("angleVsPerp.pdf");

	c2.cd(1);
	vhist2i.Draw("colztext89");
	c2.cd(2);
	vhist2j.Draw("colztext89");
	c2.cd(3);
	vhist2k.Draw("colztext89");
	c2.cd(4);
	vhist2l.Draw("colztext89");
	c2.SaveAs("angleVsMass.pdf");

	c2.Clear();
	c2.Divide(2,2);
	c2.cd(1);
	vhist12a.Draw("colztext89");
	c2.cd(2);
	vhist12b.Draw("colztext89");
	c2.SaveAs("fvzVsDeltaShift.pdf");
}

//usage - run with no argument to compare generated to (100x) resampled tracks and produce histograms
//      - run with argument (n) to produce n resampled datasets
int main(int argc, char** argv) {
	int max(-1);
	TString dir("/tmp/dcraik/");
	if(argc>1) max=atoi(argv[1]);
	if(argc>2) dir=argv[2];

	if(max<0) {
		resampleTrackIPs a(-1,dir);
		a.Loop();
	} else {
		for(int i=41; i<=max; ++i) {
			std::cout << i << " of " << max << std::endl;
			resampleTrackIPs a(i,dir);
			a.Loop();
		}
	}
	return 0;
}
