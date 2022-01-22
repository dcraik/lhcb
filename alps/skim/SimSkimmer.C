#define SimSkimmer_cxx
#include "SimSkimmer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <boost/progress.hpp>

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
// http://geomalgorithms.com/a07-_distance.html
double SimSkimmer::calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
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

//Calculate distance of closest approach of line given by point A and direction A to point B
//Equation from Wolfram Alpha
double SimSkimmer::calcDoca(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	TVector3 x0 = pb;
	TVector3 x1 = pa;
	TVector3 x2 = pa + da;

	return (x0-x1).Cross(x0-x2).Mag()/(x2-x1).Mag();
}

BPart* SimSkimmer::addB(TTree* tree, TString pre) {
	//create a struct to hold the new branches and add them to the tree
	BPart* part = new BPart;

	tree->Branch(pre+"_p"              , &(part->p));
	tree->Branch(pre+"_pt"             , &(part->pt));
	tree->Branch(pre+"_px"             , &(part->px));
	tree->Branch(pre+"_py"             , &(part->py));
	tree->Branch(pre+"_pz"             , &(part->pz));
	tree->Branch(pre+"_e"              , &(part->e));
	tree->Branch(pre+"_x"              , &(part->x));
	tree->Branch(pre+"_y"              , &(part->y));
	tree->Branch(pre+"_z"              , &(part->z));
	tree->Branch(pre+"_m"              , &(part->m));
	tree->Branch(pre+"_vtx_chi2"       , &(part->vtx_chi2));
	tree->Branch(pre+"_vtx_ndof"       , &(part->vtx_ndof));
	tree->Branch(pre+"_ip"             , &(part->ip));
	tree->Branch(pre+"_ip_chi2"        , &(part->ip_chi2));
	tree->Branch(pre+"_fd"             , &(part->fd));
	tree->Branch(pre+"_fd_chi2"        , &(part->fd_chi2));
	tree->Branch(pre+"_tau"            , &(part->tau));
	tree->Branch(pre+"_tau_err"        , &(part->tau_err));
	tree->Branch(pre+"_tau_chi2"       , &(part->tau_chi2));
	tree->Branch(pre+"_l0_hadron_dec"  , &(part->l0_hadron_dec));
	tree->Branch(pre+"_l0_hadron_tis"  , &(part->l0_hadron_tis));
	tree->Branch(pre+"_l0_hadron_tos"  , &(part->l0_hadron_tos));
	tree->Branch(pre+"_l0_photon_dec"  , &(part->l0_photon_dec));
	tree->Branch(pre+"_l0_photon_tis"  , &(part->l0_photon_tis));
	tree->Branch(pre+"_l0_photon_tos"  , &(part->l0_photon_tos));
	tree->Branch(pre+"_l0_electron_dec", &(part->l0_electron_dec));
	tree->Branch(pre+"_l0_electron_tis", &(part->l0_electron_tis));
	tree->Branch(pre+"_l0_electron_tos", &(part->l0_electron_tos));
	tree->Branch(pre+"_l0_muon_dec"    , &(part->l0_muon_dec));
	tree->Branch(pre+"_l0_muon_tis"    , &(part->l0_muon_tis));
	tree->Branch(pre+"_l0_muon_tos"    , &(part->l0_muon_tos));
	tree->Branch(pre+"_l0_dimuon_dec"  , &(part->l0_dimuon_dec));
	tree->Branch(pre+"_l0_dimuon_tis"  , &(part->l0_dimuon_tis));
	tree->Branch(pre+"_l0_dimuon_tos"  , &(part->l0_dimuon_tos));
	tree->Branch(pre+"_dira"           , &(part->dira));
	tree->Branch(pre+"_mcor"           , &(part->mcor));
	tree->Branch(pre+"_maxdoca"        , &(part->maxdoca));
	tree->Branch(pre+"_doca_kpi"       , &(part->doca_kpi));
	tree->Branch(pre+"_m_kpi"          , &(part->m_kpi));

	return part;
}

XPart* SimSkimmer::addX(TTree* tree, TString pre) {
	//create a struct to hold the new branches and add them to the tree
	XPart* part = new XPart;

	tree->Branch(pre+"_p"              , &(part->p));
	tree->Branch(pre+"_pt"             , &(part->pt));
	tree->Branch(pre+"_px"             , &(part->px));
	tree->Branch(pre+"_py"             , &(part->py));
	tree->Branch(pre+"_pz"             , &(part->pz));
	tree->Branch(pre+"_e"              , &(part->e));
	tree->Branch(pre+"_x"              , &(part->x));
	tree->Branch(pre+"_y"              , &(part->y));
	tree->Branch(pre+"_z"              , &(part->z));
	tree->Branch(pre+"_m"              , &(part->m));
	tree->Branch(pre+"_vtx_chi2"       , &(part->vtx_chi2));
	tree->Branch(pre+"_vtx_ndof"       , &(part->vtx_ndof));
	tree->Branch(pre+"_ip"             , &(part->ip));
	tree->Branch(pre+"_ip_best"        , &(part->ip_best));
	tree->Branch(pre+"_ip_bvtx"        , &(part->ip_bvtx));
	tree->Branch(pre+"_ip_chi2"        , &(part->ip_chi2));
	tree->Branch(pre+"_ip_chi2_best"   , &(part->ip_chi2_best));
	tree->Branch(pre+"_ip_chi2_bvtx"   , &(part->ip_chi2_bvtx));
	tree->Branch(pre+"_fd"             , &(part->fd));
	tree->Branch(pre+"_fd_best"        , &(part->fd_best));
	tree->Branch(pre+"_fd_bvtx"        , &(part->fd_bvtx));
	tree->Branch(pre+"_fd_chi2"        , &(part->fd_chi2));
	tree->Branch(pre+"_fd_chi2_best"   , &(part->fd_chi2_best));
	tree->Branch(pre+"_fd_chi2_bvtx"   , &(part->fd_chi2_bvtx));
	tree->Branch(pre+"_tau"            , &(part->tau));
	tree->Branch(pre+"_tau_err"        , &(part->tau_err));
	tree->Branch(pre+"_tau_chi2"       , &(part->tau_chi2));
	tree->Branch(pre+"_l0_hadron_dec"  , &(part->l0_hadron_dec));
	tree->Branch(pre+"_l0_hadron_tis"  , &(part->l0_hadron_tis));
	tree->Branch(pre+"_l0_hadron_tos"  , &(part->l0_hadron_tos));
	tree->Branch(pre+"_l0_photon_dec"  , &(part->l0_photon_dec));
	tree->Branch(pre+"_l0_photon_tis"  , &(part->l0_photon_tis));
	tree->Branch(pre+"_l0_photon_tos"  , &(part->l0_photon_tos));
	tree->Branch(pre+"_l0_electron_dec", &(part->l0_electron_dec));
	tree->Branch(pre+"_l0_electron_tis", &(part->l0_electron_tis));
	tree->Branch(pre+"_l0_electron_tos", &(part->l0_electron_tos));
	tree->Branch(pre+"_l0_muon_dec"    , &(part->l0_muon_dec));
	tree->Branch(pre+"_l0_muon_tis"    , &(part->l0_muon_tis));
	tree->Branch(pre+"_l0_muon_tos"    , &(part->l0_muon_tos));
	tree->Branch(pre+"_l0_dimuon_dec"  , &(part->l0_dimuon_dec));
	tree->Branch(pre+"_l0_dimuon_tis"  , &(part->l0_dimuon_tis));
	tree->Branch(pre+"_l0_dimuon_tos"  , &(part->l0_dimuon_tos));
	tree->Branch(pre+"_maxdoca"        , &(part->maxdoca));

	return part;
}

Track* SimSkimmer::addTrk(TTree* tree, TString pre) {
	//create a struct to hold the new branches and add them to the tree
	Track* part = new Track;

	tree->Branch(pre+"_p"              , &(part->p));
	tree->Branch(pre+"_pt"             , &(part->pt));
	tree->Branch(pre+"_px"             , &(part->px));
	tree->Branch(pre+"_py"             , &(part->py));
	tree->Branch(pre+"_pz"             , &(part->pz));
	tree->Branch(pre+"_e"              , &(part->e));
	tree->Branch(pre+"_pid"            , &(part->pid));
	tree->Branch(pre+"_q"              , &(part->q));
	tree->Branch(pre+"_ip"             , &(part->ip));
	tree->Branch(pre+"_ip_best"        , &(part->ip_best));
	tree->Branch(pre+"_ip_bvtx"        , &(part->ip_bvtx));
	tree->Branch(pre+"_ip_chi2"        , &(part->ip_chi2));
	tree->Branch(pre+"_ip_chi2_best"   , &(part->ip_chi2_best));
	tree->Branch(pre+"_ip_chi2_bvtx"   , &(part->ip_chi2_bvtx));
	tree->Branch(pre+"_pnn_e"          , &(part->pnn_e));
	tree->Branch(pre+"_pnn_mu"         , &(part->pnn_mu));
	tree->Branch(pre+"_pnn_pi"         , &(part->pnn_pi));
	tree->Branch(pre+"_pnn_k"          , &(part->pnn_k));
	tree->Branch(pre+"_pnn_p"          , &(part->pnn_p));
	tree->Branch(pre+"_pnn_ghost"      , &(part->pnn_ghost));
	tree->Branch(pre+"_ecal"           , &(part->ecal));
	tree->Branch(pre+"_hcal"           , &(part->hcal));
	tree->Branch(pre+"_prb_ghost"      , &(part->prb_ghost));
	tree->Branch(pre+"_type"           , &(part->type));
	tree->Branch(pre+"_is_mu"          , &(part->is_mu));
	tree->Branch(pre+"_vid"            , &(part->vid));
	tree->Branch(pre+"_x"              , &(part->x));
	tree->Branch(pre+"_y"              , &(part->y));
	tree->Branch(pre+"_z"              , &(part->z));
	tree->Branch(pre+"_l0_hadron_dec"  , &(part->l0_hadron_dec));
	tree->Branch(pre+"_l0_hadron_tis"  , &(part->l0_hadron_tis));
	tree->Branch(pre+"_l0_hadron_tos"  , &(part->l0_hadron_tos));
	tree->Branch(pre+"_l0_photon_dec"  , &(part->l0_photon_dec));
	tree->Branch(pre+"_l0_photon_tis"  , &(part->l0_photon_tis));
	tree->Branch(pre+"_l0_photon_tos"  , &(part->l0_photon_tos));
	tree->Branch(pre+"_l0_electron_dec", &(part->l0_electron_dec));
	tree->Branch(pre+"_l0_electron_tis", &(part->l0_electron_tis));
	tree->Branch(pre+"_l0_electron_tos", &(part->l0_electron_tos));
	tree->Branch(pre+"_l0_muon_dec"    , &(part->l0_muon_dec));
	tree->Branch(pre+"_l0_muon_tis"    , &(part->l0_muon_tis));
	tree->Branch(pre+"_l0_muon_tos"    , &(part->l0_muon_tos));
	tree->Branch(pre+"_l0_dimuon_dec"  , &(part->l0_dimuon_dec));
	tree->Branch(pre+"_l0_dimuon_tis"  , &(part->l0_dimuon_tis));
	tree->Branch(pre+"_l0_dimuon_tos"  , &(part->l0_dimuon_tos));

	return part;
}

Neutral* SimSkimmer::addNeu(TTree* tree, TString pre) {
	//create a struct to hold the new branches and add them to the tree
	Neutral* part = new Neutral;

	tree->Branch(pre+"_p"              , &(part->p));
	tree->Branch(pre+"_pt"             , &(part->pt));
	tree->Branch(pre+"_px"             , &(part->px));
	tree->Branch(pre+"_py"             , &(part->py));
	tree->Branch(pre+"_pz"             , &(part->pz));
	tree->Branch(pre+"_e"              , &(part->e));
	tree->Branch(pre+"_pid"            , &(part->pid));
	tree->Branch(pre+"_l0_hadron_dec"  , &(part->l0_hadron_dec));
	tree->Branch(pre+"_l0_hadron_tis"  , &(part->l0_hadron_tis));
	tree->Branch(pre+"_l0_hadron_tos"  , &(part->l0_hadron_tos));
	tree->Branch(pre+"_l0_photon_dec"  , &(part->l0_photon_dec));
	tree->Branch(pre+"_l0_photon_tis"  , &(part->l0_photon_tis));
	tree->Branch(pre+"_l0_photon_tos"  , &(part->l0_photon_tos));
	tree->Branch(pre+"_l0_electron_dec", &(part->l0_electron_dec));
	tree->Branch(pre+"_l0_electron_tis", &(part->l0_electron_tis));
	tree->Branch(pre+"_l0_electron_tos", &(part->l0_electron_tos));
	tree->Branch(pre+"_l0_muon_dec"    , &(part->l0_muon_dec));
	tree->Branch(pre+"_l0_muon_tis"    , &(part->l0_muon_tis));
	tree->Branch(pre+"_l0_muon_tos"    , &(part->l0_muon_tos));
	tree->Branch(pre+"_l0_dimuon_dec"  , &(part->l0_dimuon_dec));
	tree->Branch(pre+"_l0_dimuon_tis"  , &(part->l0_dimuon_tis));
	tree->Branch(pre+"_l0_dimuon_tos"  , &(part->l0_dimuon_tos));
	tree->Branch(pre+"_m"  , &(part->m));

	return part;
}


void SimSkimmer::fillB(BPart* part, int idx) {
	//check we're in range
	if(idx<0 || idx>=static_cast<int>(b_p->size())) {
		//Fill with -1
		part->p=-1.;
		part->pt=-1.;
		part->px=-1.;
		part->py=-1.;
		part->pz=-1.;
		part->e=-1.;
		part->x=-1.;
		part->y=-1.;
		part->z=-1.;
		part->m=-1.;
		part->vtx_chi2=-1.;
		part->vtx_ndof=-1.;
		part->ip=-1.;
		part->ip_chi2=-1.;
		part->fd=-1.;
		part->fd_chi2=-1.;
		part->tau=-1.;
		part->tau_err=-1.;
		part->tau_chi2=-1.;
		part->l0_hadron_dec=-1.;
		part->l0_hadron_tis=-1.;
		part->l0_hadron_tos=-1.;
		part->l0_photon_dec=-1.;
		part->l0_photon_tis=-1.;
		part->l0_photon_tos=-1.;
		part->l0_electron_dec=-1.;
		part->l0_electron_tis=-1.;
		part->l0_electron_tos=-1.;
		part->l0_muon_dec=-1.;
		part->l0_muon_tis=-1.;
		part->l0_muon_tos=-1.;
		part->l0_dimuon_dec=-1.;
		part->l0_dimuon_tis=-1.;
		part->l0_dimuon_tos=-1.;
		part->dira=-1.;
		part->mcor=-1.;
		part->maxdoca=-1.;
		part->doca_kpi = -1.;
		part->m_kpi = -1.;
	} else {
		part->p                = b_p->at(idx);
		part->pt               = b_pt->at(idx);
		part->px               = b_px->at(idx);
		part->py               = b_py->at(idx);
		part->pz               = b_pz->at(idx);
		part->e                = b_e->at(idx);
		part->x                = b_x->at(idx);
		part->y                = b_y->at(idx);
		part->z                = b_z->at(idx);
		part->m                = b_m->at(idx);
		part->vtx_chi2         = b_vtx_chi2->at(idx);
		part->vtx_ndof         = b_vtx_ndof->at(idx);
		part->ip               = b_ip->at(idx);
		part->ip_chi2          = b_ip_chi2->at(idx);
		part->fd               = b_fd->at(idx);
		part->fd_chi2          = b_fd_chi2->at(idx);
		part->tau              = b_tau->at(idx);
		part->tau_err          = b_tau_err->at(idx);
		part->tau_chi2         = b_tau_chi2->at(idx);
		part->l0_hadron_dec    = b_l0_hadron_dec->at(idx);
		part->l0_hadron_tis    = b_l0_hadron_tis->at(idx);
		part->l0_hadron_tos    = b_l0_hadron_tos->at(idx);
		part->l0_photon_dec    = b_l0_photon_dec->at(idx);
		part->l0_photon_tis    = b_l0_photon_tis->at(idx);
		part->l0_photon_tos    = b_l0_photon_tos->at(idx);
		part->l0_electron_dec  = b_l0_electron_dec->at(idx);
		part->l0_electron_tis  = b_l0_electron_tis->at(idx);
		part->l0_electron_tos  = b_l0_electron_tos->at(idx);
		part->l0_muon_dec      = b_l0_muon_dec->at(idx);
		part->l0_muon_tis      = b_l0_muon_tis->at(idx);
		part->l0_muon_tos      = b_l0_muon_tos->at(idx);
		part->l0_dimuon_dec    = b_l0_dimuon_dec->at(idx);
		part->l0_dimuon_tis    = b_l0_dimuon_tis->at(idx);
		part->l0_dimuon_tos    = b_l0_dimuon_tos->at(idx);
		//calculate dira and mCor
		int ipv = b_idx_pvr->at(idx);
		if(ipv<0 || ipv>=static_cast<int>(pvr_z->size())) {
			part->dira=-1.;
			part->mcor=-1.;
		} else {
			TVector3 pvr(pvr_x->at(ipv),pvr_y->at(ipv),pvr_z->at(ipv));
			TVector3 svr(part->x,part->y,part->z);
			TVector3 mom(part->px,part->py,part->pz);
			TVector3 fly = svr - pvr;
			part->dira = fly.Dot(mom)/(fly.Mag()*mom.Mag());

  			double pt2 = mom.Cross(fly.Unit()).Mag2();
  			part->mcor = TMath::Sqrt(part->m*part->m + pt2) + TMath::Sqrt(pt2);
		}
		//calculate max doca
		//first get list of charged tracks in decay
		std::vector<uint> trks;
		if(b_idx_k->at(idx)>-1) trks.push_back(b_idx_k->at(idx));
		if(b_idx_pi->at(idx)>-1) trks.push_back(b_idx_pi->at(idx));
		if(b_idx_x->at(idx)>-1) {
			uint idxX = b_idx_x->at(idx);
			if(x_idx_trk0->at(idxX)>-1) trks.push_back(x_idx_trk0->at(idxX));
			if(x_idx_trk1->at(idxX)>-1) trks.push_back(x_idx_trk1->at(idxX));
			if(x_idx_trk2->at(idxX)>-1) trks.push_back(x_idx_trk2->at(idxX));
			if(x_idx_trk3->at(idxX)>-1) trks.push_back(x_idx_trk3->at(idxX));
			if(x_idx_trk4->at(idxX)>-1) trks.push_back(x_idx_trk4->at(idxX));
			if(x_idx_trk5->at(idxX)>-1) trks.push_back(x_idx_trk5->at(idxX));
		}
		part->maxdoca = 0.;
		TVector3 v;
		for(uint itrk=0; itrk<trks.size(); ++itrk) {
			TVector3 xi(trk_x ->at(trks[itrk]),trk_y ->at(trks[itrk]),trk_z ->at(trks[itrk]));
			TVector3 pi(trk_px->at(trks[itrk]),trk_py->at(trks[itrk]),trk_pz->at(trks[itrk]));
			for(uint jtrk=itrk+1; jtrk<trks.size(); ++jtrk) {
				TVector3 xj(trk_x ->at(trks[jtrk]),trk_y ->at(trks[jtrk]),trk_z ->at(trks[jtrk]));
				TVector3 pj(trk_px->at(trks[jtrk]),trk_py->at(trks[jtrk]),trk_pz->at(trks[jtrk]));
				double doca = calcDoca(v,xi,pi,xj,pj);
				//if we have the K and pi from the B then save DOCA
				if(itrk==0 && jtrk==1) {
					if(b_idx_k->at(idx)>-1 && b_idx_pi->at(idx)>-1) {
						part->doca_kpi = doca;
					} else {
						part->doca_kpi = -1.;
					}
				}
				if(doca>part->maxdoca) part->maxdoca = doca;
			}
		}
		if(b_idx_k->at(idx)>-1 && b_idx_pi->at(idx)>-1) {
			uint idxK = b_idx_k->at(idx);
			uint idxPi = b_idx_pi->at(idx);
			TLorentzVector kst(trk_px->at(idxK)+trk_px->at(idxPi),
					   trk_py->at(idxK)+trk_py->at(idxPi),
					   trk_pz->at(idxK)+trk_pz->at(idxPi),
					   trk_e->at(idxK)+trk_e->at(idxPi));
			part->m_kpi = kst.M();
		} else {
			part->m_kpi = -1.;
		}
	}
}

void SimSkimmer::fillX(XPart* part, int idx) {
	//check we're in range
	if(idx<0 || idx>=static_cast<int>(x_p->size())) {
		//Fill with -1
		part->p=-1.;
		part->pt=-1.;
		part->px=-1.;
		part->py=-1.;
		part->pz=-1.;
		part->e=-1.;
		part->x=-1.;
		part->y=-1.;
		part->z=-1.;
		part->m=-1.;
		part->vtx_chi2=-1.;
		part->vtx_ndof=-1.;
		part->ip=-1.;
		part->ip_best=-1.;
		part->ip_bvtx=-1.;
		part->ip_chi2=-1.;
		part->ip_chi2_best=-1.;
		part->ip_chi2_bvtx=-1.;
		part->fd=-1.;
		part->fd_best=-1.;
		part->fd_bvtx=-1.;
		part->fd_chi2=-1.;
		part->fd_chi2_best=-1.;
		part->fd_chi2_bvtx=-1.;
		part->tau=-1.;
		part->tau_err=-1.;
		part->tau_chi2=-1.;
		part->l0_hadron_dec=-1.;
		part->l0_hadron_tis=-1.;
		part->l0_hadron_tos=-1.;
		part->l0_photon_dec=-1.;
		part->l0_photon_tis=-1.;
		part->l0_photon_tos=-1.;
		part->l0_electron_dec=-1.;
		part->l0_electron_tis=-1.;
		part->l0_electron_tos=-1.;
		part->l0_muon_dec=-1.;
		part->l0_muon_tis=-1.;
		part->l0_muon_tos=-1.;
		part->l0_dimuon_dec=-1.;
		part->l0_dimuon_tis=-1.;
		part->l0_dimuon_tos=-1.;
		part->maxdoca=-1.;
	} else {
		part->p                = x_p->at(idx);
		part->pt               = x_pt->at(idx);
		part->px               = x_px->at(idx);
		part->py               = x_py->at(idx);
		part->pz               = x_pz->at(idx);
		part->e                = x_e->at(idx);
		part->x                = x_x->at(idx);
		part->y                = x_y->at(idx);
		part->z                = x_z->at(idx);
		part->m                = x_m->at(idx);
		part->vtx_chi2         = x_vtx_chi2->at(idx);
		part->vtx_ndof         = x_vtx_ndof->at(idx);
		part->ip               = x_ip->at(idx);
		part->ip_best          = x_ip_best->at(idx);
		part->ip_bvtx          = x_ip_bvtx->at(idx);
		part->ip_chi2          = x_ip_chi2->at(idx);
		part->ip_chi2_best     = x_ip_chi2_best->at(idx);
		part->ip_chi2_bvtx     = x_ip_chi2_bvtx->at(idx);
		part->fd               = x_fd->at(idx);
		part->fd_best          = x_fd_best->at(idx);
		part->fd_bvtx          = x_fd_bvtx->at(idx);
		part->fd_chi2          = x_fd_chi2->at(idx);
		part->fd_chi2_best     = x_fd_chi2_best->at(idx);
		part->fd_chi2_bvtx     = x_fd_chi2_bvtx->at(idx);
		part->tau              = x_tau->at(idx);
		part->tau_err          = x_tau_err->at(idx);
		part->tau_chi2         = x_tau_chi2->at(idx);
		part->l0_hadron_dec    = x_l0_hadron_dec->at(idx);
		part->l0_hadron_tis    = x_l0_hadron_tis->at(idx);
		part->l0_hadron_tos    = x_l0_hadron_tos->at(idx);
		part->l0_photon_dec    = x_l0_photon_dec->at(idx);
		part->l0_photon_tis    = x_l0_photon_tis->at(idx);
		part->l0_photon_tos    = x_l0_photon_tos->at(idx);
		part->l0_electron_dec  = x_l0_electron_dec->at(idx);
		part->l0_electron_tis  = x_l0_electron_tis->at(idx);
		part->l0_electron_tos  = x_l0_electron_tos->at(idx);
		part->l0_muon_dec      = x_l0_muon_dec->at(idx);
		part->l0_muon_tis      = x_l0_muon_tis->at(idx);
		part->l0_muon_tos      = x_l0_muon_tos->at(idx);
		part->l0_dimuon_dec    = x_l0_dimuon_dec->at(idx);
		part->l0_dimuon_tis    = x_l0_dimuon_tis->at(idx);
		part->l0_dimuon_tos    = x_l0_dimuon_tos->at(idx);

		//calculate max doca
		//first get list of charged tracks in decay
		std::vector<uint> trks;
		if(x_idx_trk0->at(idx)>-1) trks.push_back(x_idx_trk0->at(idx));
		if(x_idx_trk1->at(idx)>-1) trks.push_back(x_idx_trk1->at(idx));
		if(x_idx_trk2->at(idx)>-1) trks.push_back(x_idx_trk2->at(idx));
		if(x_idx_trk3->at(idx)>-1) trks.push_back(x_idx_trk3->at(idx));
		if(x_idx_trk4->at(idx)>-1) trks.push_back(x_idx_trk4->at(idx));
		if(x_idx_trk5->at(idx)>-1) trks.push_back(x_idx_trk5->at(idx));

		part->maxdoca = 0.;
		TVector3 v;
		for(uint itrk=0; itrk<trks.size(); ++itrk) {
			TVector3 xi(trk_x ->at(trks[itrk]),trk_y ->at(trks[itrk]),trk_z ->at(trks[itrk]));
			TVector3 pi(trk_px->at(trks[itrk]),trk_py->at(trks[itrk]),trk_pz->at(trks[itrk]));
			for(uint jtrk=itrk+1; jtrk<trks.size(); ++jtrk) {
				TVector3 xj(trk_x ->at(trks[jtrk]),trk_y ->at(trks[jtrk]),trk_z ->at(trks[jtrk]));
				TVector3 pj(trk_px->at(trks[jtrk]),trk_py->at(trks[jtrk]),trk_pz->at(trks[jtrk]));
				double doca = calcDoca(v,xi,pi,xj,pj);
				if(doca>part->maxdoca) part->maxdoca = doca;
			}
		}
	}
}

void SimSkimmer::fillTrk(Track* part, int idx) {
	//check we're in range
	if(idx<0 || idx>=static_cast<int>(trk_p->size())) {
		//Fill with -1
		part->p=-1.;
		part->pt=-1.;
		part->px=-1.;
		part->py=-1.;
		part->pz=-1.;
		part->e=-1.;
		part->pid=-1.;
		part->q=-1.;
		part->ip=-1.;
		part->ip_best=-1.;
		part->ip_bvtx=-1.;
		part->ip_chi2=-1.;
		part->ip_chi2_best=-1.;
		part->ip_chi2_bvtx=-1.;
		part->pnn_e=-1.;
		part->pnn_mu=-1.;
		part->pnn_pi=-1.;
		part->pnn_k=-1.;
		part->pnn_p=-1.;
		part->pnn_ghost=-1.;
		part->ecal=-1.;
		part->hcal=-1.;
		part->prb_ghost=-1.;
		part->type=-1.;
		part->is_mu=-1.;
		part->vid=-1.;
		part->x=-1.;
		part->y=-1.;
		part->z=-1.;
		part->l0_hadron_dec=-1.;
		part->l0_hadron_tis=-1.;
		part->l0_hadron_tos=-1.;
		part->l0_photon_dec=-1.;
		part->l0_photon_tis=-1.;
		part->l0_photon_tos=-1.;
		part->l0_electron_dec=-1.;
		part->l0_electron_tis=-1.;
		part->l0_electron_tos=-1.;
		part->l0_muon_dec=-1.;
		part->l0_muon_tis=-1.;
		part->l0_muon_tos=-1.;
		part->l0_dimuon_dec=-1.;
		part->l0_dimuon_tis=-1.;
		part->l0_dimuon_tos=-1.;
	} else {
		part->p               = trk_p->at(idx);
		part->pt              = trk_pt->at(idx);
		part->px              = trk_px->at(idx);
		part->py              = trk_py->at(idx);
		part->pz              = trk_pz->at(idx);
		part->e               = trk_e->at(idx);
		part->pid             = trk_pid->at(idx);
		part->q               = trk_q->at(idx);
		part->ip              = trk_ip->at(idx);
		part->ip_best         = trk_ip_best->at(idx);
		part->ip_bvtx         = trk_ip_bvtx->at(idx);
		part->ip_chi2         = trk_ip_chi2->at(idx);
		part->ip_chi2_best    = trk_ip_chi2_best->at(idx);
		part->ip_chi2_bvtx    = trk_ip_chi2_bvtx->at(idx);
		part->pnn_e           = trk_pnn_e->at(idx);
		part->pnn_mu          = trk_pnn_mu->at(idx);
		part->pnn_pi          = trk_pnn_pi->at(idx);
		part->pnn_k           = trk_pnn_k->at(idx);
		part->pnn_p           = trk_pnn_p->at(idx);
		part->pnn_ghost       = trk_pnn_ghost->at(idx);
		part->ecal            = trk_ecal->at(idx);
		part->hcal            = trk_hcal->at(idx);
		part->prb_ghost       = trk_prb_ghost->at(idx);
		part->type            = trk_type->at(idx);
		part->is_mu           = trk_is_mu->at(idx);
		part->vid             = trk_vid->at(idx);
		part->x               = trk_x->at(idx);
		part->y               = trk_y->at(idx);
		part->z               = trk_z->at(idx);
		part->l0_hadron_dec   = trk_l0_hadron_dec->at(idx);
		part->l0_hadron_tis   = trk_l0_hadron_tis->at(idx);
		part->l0_hadron_tos   = trk_l0_hadron_tos->at(idx);
		part->l0_photon_dec   = trk_l0_photon_dec->at(idx);
		part->l0_photon_tis   = trk_l0_photon_tis->at(idx);
		part->l0_photon_tos   = trk_l0_photon_tos->at(idx);
		part->l0_electron_dec = trk_l0_electron_dec->at(idx);
		part->l0_electron_tis = trk_l0_electron_tis->at(idx);
		part->l0_electron_tos = trk_l0_electron_tos->at(idx);
		part->l0_muon_dec     = trk_l0_muon_dec->at(idx);
		part->l0_muon_tis     = trk_l0_muon_tis->at(idx);
		part->l0_muon_tos     = trk_l0_muon_tos->at(idx);
		part->l0_dimuon_dec   = trk_l0_dimuon_dec->at(idx);
		part->l0_dimuon_tis   = trk_l0_dimuon_tis->at(idx);
		part->l0_dimuon_tos   = trk_l0_dimuon_tos->at(idx);
	}
}

void SimSkimmer::fillNeu(Neutral* part, int idx) {
	//check we're in range
	if(idx<0 || idx>=static_cast<int>(neu_p->size())) {
		//Fill with -1
		part->p=-1.;
		part->pt=-1.;
		part->px=-1.;
		part->py=-1.;
		part->pz=-1.;
		part->e=-1.;
		part->pid=-1.;
		part->l0_hadron_dec=-1.;
		part->l0_hadron_tis=-1.;
		part->l0_hadron_tos=-1.;
		part->l0_photon_dec=-1.;
		part->l0_photon_tis=-1.;
		part->l0_photon_tos=-1.;
		part->l0_electron_dec=-1.;
		part->l0_electron_tis=-1.;
		part->l0_electron_tos=-1.;
		part->l0_muon_dec=-1.;
		part->l0_muon_tis=-1.;
		part->l0_muon_tos=-1.;
		part->l0_dimuon_dec=-1.;
		part->l0_dimuon_tis=-1.;
		part->l0_dimuon_tos=-1.;
		part->m=-1.;
	} else {
		part->p                = neu_p->at(idx);
		part->pt               = neu_pt->at(idx);
		part->px               = neu_px->at(idx);
		part->py               = neu_py->at(idx);
		part->pz               = neu_pz->at(idx);
		part->e                = neu_e->at(idx);
		part->pid              = neu_pid->at(idx);
		part->l0_hadron_dec    = neu_l0_hadron_dec->at(idx);
		part->l0_hadron_tis    = neu_l0_hadron_tis->at(idx);
		part->l0_hadron_tos    = neu_l0_hadron_tos->at(idx);
		part->l0_photon_dec    = neu_l0_photon_dec->at(idx);
		part->l0_photon_tis    = neu_l0_photon_tis->at(idx);
		part->l0_photon_tos    = neu_l0_photon_tos->at(idx);
		part->l0_electron_dec  = neu_l0_electron_dec->at(idx);
		part->l0_electron_tis  = neu_l0_electron_tis->at(idx);
		part->l0_electron_tos  = neu_l0_electron_tos->at(idx);
		part->l0_muon_dec      = neu_l0_muon_dec->at(idx);
		part->l0_muon_tis      = neu_l0_muon_tis->at(idx);
		part->l0_muon_tos      = neu_l0_muon_tos->at(idx);
		part->l0_dimuon_dec    = neu_l0_dimuon_dec->at(idx);
		part->l0_dimuon_tis    = neu_l0_dimuon_tis->at(idx);
		part->l0_dimuon_tos    = neu_l0_dimuon_tos->at(idx);
		TLorentzVector p4(neu_px->at(idx), neu_py->at(idx), neu_pz->at(idx), neu_e->at(idx));
		part->m=p4.M();
	}
}

bool SimSkimmer::checkTrue_Kpi_pipineu(int ib, int ix, int ik, int ipi, int ixpip, int ixpim, int ineu, int neuPID) {
	int genk    = trk_idx_gen->at(ik);
	int genpi   = trk_idx_gen->at(ipi);
	int genxpip = trk_idx_gen->at(ixpip);
	int genxpim = trk_idx_gen->at(ixpim);
	int genneu  = neu_idx_gen->at(ineu);

	//check all exist
	if(genk    < 0 || genk    > static_cast<int>(gen_pid->size())) return false;
	if(genpi   < 0 || genpi   > static_cast<int>(gen_pid->size())) return false;
	if(genxpip < 0 || genxpip > static_cast<int>(gen_pid->size())) return false;
	if(genxpim < 0 || genxpim > static_cast<int>(gen_pid->size())) return false;
	if(genneu  < 0 || genneu  > static_cast<int>(gen_pid->size())) return false;

	//check PIDs
	if(TMath::Abs(gen_pid->at(genk   )) != 321) return false;
	if(TMath::Abs(gen_pid->at(genpi  )) != 211) return false;
	if(TMath::Abs(gen_pid->at(genxpip)) != 211) return false;
	if(TMath::Abs(gen_pid->at(genxpim)) != 211) return false;
	if(TMath::Abs(gen_pid->at(genneu )) != neuPID) return false;

	//check topology (allow resonant or non-resonant Kpi)
	int genx = gen_idx_prnt->at(genxpip);

	//X decay
	if(genx < 0 || genx > static_cast<int>(gen_pid->size())) return false;
	if(genx != gen_idx_prnt->at(genxpim)) return false;
	if(genx != gen_idx_prnt->at(genneu)) return false;

	//K* decay
	int genkst = gen_idx_prnt->at(genk);

	if(genkst < 0 || genkst > static_cast<int>(gen_pid->size())) return false;
	if(genkst != gen_idx_prnt->at(genpi)) return false;

	//B decay
	int genb = gen_idx_prnt->at(genx);
	if(genb < 0 || genb > static_cast<int>(gen_pid->size())) return false;

	if(genb == genkst) {
		//non-resonant
	} else if(genb == gen_idx_prnt->at(genkst)) {
		//resonant
	} else {
		return false;
	}

	return true;

}

void SimSkimmer::Loop_Kpi_pipineu(int mode, int nmax)
{
	std::vector<uint> counts(19,0);
	std::vector<TString> countNames = {"all","mass(B)","mcor(B)","pT(B)","dira(B)","chi2_vtx(B)","doca(Kpi)","FD(B)","chi2_FD(B or X)","chi2_vtx(X)","forward or prompt X",//"chi2_IP(X)",
		                           "PID(neu)","pT(trk)","p(trk)","chi2_IP_PV(trk)","chi2_IP_SV(trk)","ghostProb(trk)","PNN_K(K)","sum(chi2_IP(trk))"};

	if (fChain == 0) return;

	TString fname=savedir;
	int neuPID(0);
	switch(mode) {
		case 5:
			fname += "kpi3pi-sim";
			neuPID=111;
			break;
		case 11:
			fname += "kpietapipi-sim";
			neuPID=221;
			break;
		default:
			std::cout << "Unknown mode" << std::endl;
			return;
	}
	if(nmax>-1) {
		fname+="_";
		fname+=nmax;
	}
	fname+=".root";

	TFile* fout = TFile::Open(fname,"RECREATE");
	TTree* tout = new TTree("T","");

	BPart* b = addB(tout, "b");
	XPart* x = addX(tout, "x");
	Track* k = addTrk(tout, "k");
	Track* pi = addTrk(tout, "pi");
	Track* xpip = addTrk(tout, "xpip");
	Track* xpim = addTrk(tout, "xpim");
	Neutral* neu(0);
	switch(mode) {
		case 5:
			neu = addNeu(tout, "xpiz");
			break;
		case 11:
			neu = addNeu(tout, "xeta");
			break;
	}

	Long64_t nentries = fChain->GetEntries();
	if(nmax>-1 && nmax<nentries) nentries=nmax;

	Long64_t nbytes = 0, nb = 0;

	//indices for B children
	int ix, ik, ipi, ixpip, ixpim, ineu;

	//track multiple candidates in event
	int iCand(0);
	tout->Branch("iCand", &iCand);

	boost::progress_display progress(nentries);
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		++progress;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		iCand=0;
		for(int ib=0; ib<static_cast<int>(b_mode->size()); ++ib) {
			if(b_mode->at(ib)!=mode) continue;

			//B children
			ix = b_idx_x->at(ib);
			ik = b_idx_k->at(ib);
			ipi = b_idx_pi->at(ib);

			//X children
			ineu = x_idx_neu0->at(ix);
			ixpip = x_idx_trk0->at(ix);
			ixpim = x_idx_trk1->at(ix);

			//flip pions if pim came first
			if(trk_pid->at(ixpip) < 0.) {
				ixpip = x_idx_trk1->at(ix);
				ixpim = x_idx_trk0->at(ix);
			}

			//truth match
			if(!checkTrue_Kpi_pipineu(ib,ix,ik,ipi,ixpip,ixpim,ineu,neuPID)) continue;
			++counts[0];

			//now fill all branches
			fillB(b, ib);
			fillX(x, ix);
			fillTrk(k, ik);
			fillTrk(pi, ipi);
			fillTrk(xpip, ixpip);
			fillTrk(xpim, ixpim);
			fillNeu(neu, ineu);

			//apply any cuts

			//B0 mass
			if(b->m<4800. || b->m>5800.) continue;
			++counts[1];
			//B0 corrected mass
			if(b->mcor<4000. || b->mcor>10000.) continue;
			++counts[2];
			//B0 pT
			if(b->pt<3000.) continue;
			++counts[3];
			//B cos(theta_dira)
			if(b->dira<0.9995) continue;
			++counts[4];
			//B0 vertex chi2
			if(b->vtx_chi2/b->vtx_ndof>4.) continue;
			++counts[5];
			//Kpi DOCA
			if(b->doca_kpi>0.2) continue;
			++counts[6];
			//B FD
			if(b->fd<3.) continue;
			++counts[7];
			//B&X FD chi2
			if(b->fd_chi2<500. && x->fd_chi2<500.) continue;
			++counts[8];

			//X vertex chi2
			if(x->vtx_chi2/x->vtx_ndof>4.) continue;
			++counts[9];
			//X forwards of B or prompt from B
			if(x->z - b->z < -20. || (b->maxdoca>0.2&&x->z-b->z<0.) || (b->maxdoca>1.0&&x->z-b->z<20.)) continue;
			++counts[10];
			//X IP chi2 from B decay vertex
			//if(x->ip_chi2_bvtx>10.) continue;
			//++counts[11];

			//eta/pi0
			if(neu->pid!=neuPID) continue;
			++counts[11];

			//track pT
			if(k->pt<100. || pi->pt<100. || xpip->pt<100. || xpim->pt<100.) continue;
			++counts[12];
			//track p
			if(k->p<1500. || pi->p<1500. || xpip->p<1500. || xpim->p<1500.) continue;
			++counts[13];
			//track IP chi2 from PV
			if(k->ip_chi2_best<10. || pi->ip_chi2_best<10. || xpip->ip_chi2_best<10. || xpim->ip_chi2_best<10.) continue;
			++counts[14];
			//track IP chi2 from SV
			if(k->ip_chi2_bvtx>10. || pi->ip_chi2_bvtx>10.) continue;
			++counts[15];
			//track ghost prob
			if(k->prb_ghost>0.1 || pi->prb_ghost>0.1 || xpip->prb_ghost>0.1 || xpim->prb_ghost>0.1) continue;
			++counts[16];
			//track PID
			if(k->pnn_k<0.3) continue;
			++counts[17];

			//sum of IP chi2
			if(k->ip_chi2_best + pi->ip_chi2_best + xpip->ip_chi2_best + xpim->ip_chi2_best < 200.) continue;
			++counts[18];

			tout->Fill();
			++iCand;
		}

	}
	for (uint i=0; i<counts.size(); ++i) {
		printf("%20s % 5d %.2f %.2f\n", countNames[i].Data(), counts[i], (i>0?(static_cast<double>(counts[i])/counts[i-1]):1.),static_cast<double>(counts[i])/counts[0]);
	}
	tout->Write();
	fout->Close();
}

int main(int argc, char** argv) {
	if(argc<2) {
		std::cout << "Usage: " << argv[0] << " <mode> [nmax]" << std::endl;
		return 0;
	}
	int mode = atoi(argv[1]);
	int nmax(-1), file(-1);
	if(argc>2) nmax = atoi(argv[2]);
	if(argc>3) file = atoi(argv[3]);

	SimSkimmer a(file);
	switch(mode) {
		case 5:
		case 11:
			a.Loop_Kpi_pipineu(mode,nmax);
			break;
		default:
			std::cout << "Unknown mode " << mode << std::endl;
	}
	return 0;
}
