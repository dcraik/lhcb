#define TreeSkimmer_cxx
#include "TreeSkimmer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <boost/progress.hpp>

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
// http://geomalgorithms.com/a07-_distance.html
double TreeSkimmer::calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
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
double TreeSkimmer::calcDoca(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	TVector3 x0 = pb;
	TVector3 x1 = pa;
	TVector3 x2 = pa + da;

	return (x0-x1).Cross(x0-x2).Mag()/(x2-x1).Mag();
}

BPart* TreeSkimmer::addB(TTree* tree, TString pre) {
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

XPart* TreeSkimmer::addX(TTree* tree, TString pre) {
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

Track* TreeSkimmer::addTrk(TTree* tree, TString pre) {
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

Neutral* TreeSkimmer::addNeu(TTree* tree, TString pre) {
	//create a struct to hold the new branches and add them to the tree
	Neutral* part = new Neutral;

	tree->Branch(pre+"_p"              , &(part->p));
	tree->Branch(pre+"_pt"             , &(part->pt));
	tree->Branch(pre+"_px"             , &(part->px));
	tree->Branch(pre+"_py"             , &(part->py));
	tree->Branch(pre+"_pz"             , &(part->pz));
	tree->Branch(pre+"_e"              , &(part->e));
	tree->Branch(pre+"_pid"            , &(part->pid));
	tree->Branch(pre+"_cl"             , &(part->cl));
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

DTFVars* TreeSkimmer::addDTF(TTree* tree, TString pre,
		             std::vector<TString> trks,
			     std::vector<TString> neus) {
	int ntrk = trks.size();
	int nneu = neus.size();
	if(ntrk>8) {
		std::cout << "WARNING in TreeSkimmer::addDTF : maximum of 8 tracks supported" << std::endl;
		ntrk=8;
	}
	if(nneu>8) {
		std::cout << "WARNING in TreeSkimmer::addDTF : maximum of 8 neutrals supported" << std::endl;
		nneu=8;
	}
	//create a struct to hold the new branches and add them to the tree
	DTFVars* dtf = new DTFVars(ntrk, nneu);

	tree->Branch(pre+"_b_px"  , &(dtf->px[0]));
	tree->Branch(pre+"_b_py"  , &(dtf->py[0]));
	tree->Branch(pre+"_b_pz"  , &(dtf->pz[0]));
	tree->Branch(pre+"_b_e"   , &(dtf->pe[0]));
	tree->Branch(pre+"_b_m"   , &(dtf->pm[0]));
	tree->Branch(pre+"_b_ctau", &(dtf->ctau[0]));
	tree->Branch(pre+"_x_m"   , &(dtf->pm[1]));
	tree->Branch(pre+"_x_ctau", &(dtf->ctau[1]));
	int ipart(2);
	for(int itrk=0; itrk<ntrk; ++itrk) {
		tree->Branch(pre+"_"+trks[itrk]+"_px"  , &(dtf->px[ipart]));
		tree->Branch(pre+"_"+trks[itrk]+"_py"  , &(dtf->py[ipart]));
		tree->Branch(pre+"_"+trks[itrk]+"_pz"  , &(dtf->pz[ipart]));
		tree->Branch(pre+"_"+trks[itrk]+"_e"   , &(dtf->pe[ipart]));
		tree->Branch(pre+"_"+trks[itrk]+"_m"   , &(dtf->pm[ipart]));
		++ipart;
	}
	for(int ineu=0; ineu<nneu; ++ineu) {
		tree->Branch(pre+"_"+neus[ineu]+"_px"  , &(dtf->px[ipart]));
		tree->Branch(pre+"_"+neus[ineu]+"_py"  , &(dtf->py[ipart]));
		tree->Branch(pre+"_"+neus[ineu]+"_pz"  , &(dtf->pz[ipart]));
		tree->Branch(pre+"_"+neus[ineu]+"_e"   , &(dtf->pe[ipart]));
		tree->Branch(pre+"_"+neus[ineu]+"_m"   , &(dtf->pm[ipart]));
		++ipart;
	}

	tree->Branch(pre+"_chi2"  , &(dtf->chi2));
	tree->Branch(pre+"_ndof"  , &(dtf->ndof));

	return dtf;
}


void TreeSkimmer::fillB(BPart* part, int idx) {
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

void TreeSkimmer::fillX(XPart* part, int idx) {
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

void TreeSkimmer::fillTrk(Track* part, int idx) {
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

void TreeSkimmer::fillNeu(Neutral* part, int idx) {
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
		part->cl=-1.;
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
		part->cl               = neu_cl->at(idx);
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

void TreeSkimmer::fillDTF(DTFVars* dtf, int idx) {
	//check we're in range
	if(idx<0 || idx>=static_cast<int>(b_p->size())) {
	} else {
		dtf->px[0]   = b_dtf_px->at(idx);
		dtf->py[0]   = b_dtf_py->at(idx);
		dtf->pz[0]   = b_dtf_pz->at(idx);
		dtf->pe[0]   = b_dtf_e->at(idx);
		dtf->pm[0]   = b_dtf_m->at(idx);
		dtf->ctau[0] = b_dtf_ctau->at(idx);
		dtf->pm[1]   = b_dtf_x_m->at(idx);
		dtf->ctau[1] = b_dtf_x_ctau->at(idx);

		int ipart(2);
		std::vector<std::vector<double>*> dtf_trk_px = {b_dtf_trk0_px,b_dtf_trk1_px,b_dtf_trk2_px,b_dtf_trk3_px,
						                b_dtf_trk4_px,b_dtf_trk5_px,b_dtf_trk6_px,b_dtf_trk7_px};
		std::vector<std::vector<double>*> dtf_trk_py = {b_dtf_trk0_py,b_dtf_trk1_py,b_dtf_trk2_py,b_dtf_trk3_py,
						                b_dtf_trk4_py,b_dtf_trk5_py,b_dtf_trk6_py,b_dtf_trk7_py};
		std::vector<std::vector<double>*> dtf_trk_pz = {b_dtf_trk0_pz,b_dtf_trk1_pz,b_dtf_trk2_pz,b_dtf_trk3_pz,
						                b_dtf_trk4_pz,b_dtf_trk5_pz,b_dtf_trk6_pz,b_dtf_trk7_pz};
		std::vector<std::vector<double>*> dtf_trk_e  = {b_dtf_trk0_e ,b_dtf_trk1_e ,b_dtf_trk2_e ,b_dtf_trk3_e ,
						                b_dtf_trk4_e ,b_dtf_trk5_e ,b_dtf_trk6_e ,b_dtf_trk7_e };
		std::vector<std::vector<double>*> dtf_trk_m  = {b_dtf_trk0_m ,b_dtf_trk1_m ,b_dtf_trk2_m ,b_dtf_trk3_m ,
							        b_dtf_trk4_m ,b_dtf_trk5_m ,b_dtf_trk6_m ,b_dtf_trk7_m };
		std::vector<std::vector<double>*> dtf_neu_px = {b_dtf_neu0_px,b_dtf_neu1_px,b_dtf_neu2_px,b_dtf_neu3_px,
						                b_dtf_neu4_px,b_dtf_neu5_px,b_dtf_neu6_px,b_dtf_neu7_px};
		std::vector<std::vector<double>*> dtf_neu_py = {b_dtf_neu0_py,b_dtf_neu1_py,b_dtf_neu2_py,b_dtf_neu3_py,
						                b_dtf_neu4_py,b_dtf_neu5_py,b_dtf_neu6_py,b_dtf_neu7_py};
		std::vector<std::vector<double>*> dtf_neu_pz = {b_dtf_neu0_pz,b_dtf_neu1_pz,b_dtf_neu2_pz,b_dtf_neu3_pz,
						                b_dtf_neu4_pz,b_dtf_neu5_pz,b_dtf_neu6_pz,b_dtf_neu7_pz};
		std::vector<std::vector<double>*> dtf_neu_e  = {b_dtf_neu0_e ,b_dtf_neu1_e ,b_dtf_neu2_e ,b_dtf_neu3_e ,
						                b_dtf_neu4_e ,b_dtf_neu5_e ,b_dtf_neu6_e ,b_dtf_neu7_e };
		std::vector<std::vector<double>*> dtf_neu_m  = {b_dtf_neu0_m ,b_dtf_neu1_m ,b_dtf_neu2_m ,b_dtf_neu3_m ,
							        b_dtf_neu4_m ,b_dtf_neu5_m ,b_dtf_neu6_m ,b_dtf_neu7_m };

		for(int itrk=0; itrk<dtf->_ntrk && itrk<8; ++itrk) {
			dtf->px[ipart] = dtf_trk_px[itrk]->at(idx);
			dtf->py[ipart] = dtf_trk_py[itrk]->at(idx);
			dtf->pz[ipart] = dtf_trk_pz[itrk]->at(idx);
			dtf->pe[ipart] = dtf_trk_e[itrk]->at(idx);
			dtf->pm[ipart] = dtf_trk_m[itrk]->at(idx);
			++ipart;
		}
		for(int ineu=0; ineu<dtf->_nneu && ineu<8; ++ineu) {
			dtf->px[ipart] = dtf_neu_px[ineu]->at(idx);
			dtf->py[ipart] = dtf_neu_py[ineu]->at(idx);
			dtf->pz[ipart] = dtf_neu_pz[ineu]->at(idx);
			dtf->pe[ipart] = dtf_neu_e[ineu]->at(idx);
			dtf->pm[ipart] = dtf_neu_m[ineu]->at(idx);
			++ipart;
		}

		dtf->chi2 = b_dtf_chi2->at(idx);
		dtf->ndof = b_dtf_ndof->at(idx);
	}
}


void TreeSkimmer::Loop_Kpi_pipineu(int mode, int nmax)
{
	std::vector<uint> counts(19,0);
	std::vector<TString> countNames = {"all","mass(B)","mcor(B)","pT(B)","dira(B)","chi2_vtx(B)","doca(Kpi)","FD(B)","chi2_FD(B or X)","chi2_vtx(X)","forward or prompt X",//"chi2_IP(X)",
		                           "PID(neu)","pT(trk)","p(trk)","chi2_IP_PV(trk)","chi2_IP_SV(trk)","ghostProb(trk)","PNN_K(K)","sum(chi2_IP(trk))"};
	std::vector<TH1D*> bmhists;
	std::vector<TH1D*> xmhists;
	for(int i=0; i<19; ++i) {
		bmhists.push_back(new TH1D(TString::Format("hb_%d",i),"",100,5000,6000));
		xmhists.push_back(new TH1D(TString::Format("hx_%d",i),"",100,0,5000));
	}

	if (fChain == 0) return;

	TString fname=savedir;
	int neuPID(0);
	switch(mode) {
		case 5:
			fname += "kpi3pi";
			neuPID=111;
			break;
		case 11:
			fname += "kpietapipi";
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

	std::vector<TString> trks{"k","pi","xpip","xpim"};
	std::vector<TString> neus;

	switch(mode) {
		case 5:
			neu = addNeu(tout, "xpiz");
			neus.push_back("xpiz");
			break;
		case 11:
			neu = addNeu(tout, "xeta");
			neus.push_back("xeta");
			break;
	}
	DTFVars* dtf = addDTF(tout, "dtf", trks, neus);

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

			//now fill all branches
			fillB(b, ib);
			fillX(x, ix);
			fillTrk(k, ik);
			fillTrk(pi, ipi);
			fillTrk(xpip, ixpip);
			fillTrk(xpim, ixpim);
			fillNeu(neu, ineu);
			//DTF branches are all stored in B candidate
			fillDTF(dtf, ib);

			//apply any cuts
			++counts[0];
			bmhists[0]->Fill(dtf->pm[0]);
			xmhists[0]->Fill(dtf->pm[1]);

			//B0 mass
			if(b->m<4800. || b->m>5800.) continue;
			++counts[1];
			bmhists[1]->Fill(dtf->pm[0]);
			xmhists[1]->Fill(dtf->pm[1]);
			//B0 corrected mass
			if(b->mcor<4000. || b->mcor>10000.) continue;
			++counts[2];
			bmhists[2]->Fill(dtf->pm[0]);
			xmhists[2]->Fill(dtf->pm[1]);
			//B0 pT
			if(b->pt<3000.) continue;
			++counts[3];
			bmhists[3]->Fill(dtf->pm[0]);
			xmhists[3]->Fill(dtf->pm[1]);
			//B cos(theta_dira)
			if(b->dira<0.9995) continue;
			++counts[4];
			bmhists[4]->Fill(dtf->pm[0]);
			xmhists[4]->Fill(dtf->pm[1]);
			//B0 vertex chi2
			if(b->vtx_chi2/b->vtx_ndof>4.) continue;
			++counts[5];
			bmhists[5]->Fill(dtf->pm[0]);
			xmhists[5]->Fill(dtf->pm[1]);
			//Kpi DOCA
			if(b->doca_kpi>0.2) continue;
			++counts[6];
			bmhists[6]->Fill(dtf->pm[0]);
			xmhists[6]->Fill(dtf->pm[1]);
			//B FD
			if(b->fd<3.) continue;
			++counts[7];
			bmhists[7]->Fill(dtf->pm[0]);
			xmhists[7]->Fill(dtf->pm[1]);
			//B&X FD chi2
			if(b->fd_chi2<500. && x->fd_chi2<500.) continue;
			++counts[8];
			bmhists[8]->Fill(dtf->pm[0]);
			xmhists[8]->Fill(dtf->pm[1]);

			//X vertex chi2
			if(x->vtx_chi2/x->vtx_ndof>4.) continue;
			++counts[9];
			bmhists[9]->Fill(dtf->pm[0]);
			xmhists[9]->Fill(dtf->pm[1]);
			//X forwards of B or prompt from B
			if(x->z - b->z < -20. || (b->maxdoca>0.2&&x->z-b->z<0.) || (b->maxdoca>1.0&&x->z-b->z<20.)) continue;
			++counts[10];
			bmhists[10]->Fill(dtf->pm[0]);
			xmhists[10]->Fill(dtf->pm[1]);
			//X IP chi2 from B decay vertex
			//if(x->ip_chi2_bvtx>10.) continue;
			//++counts[11];

			//eta/pi0
			if(neu->pid!=neuPID) continue;
			++counts[11];
			bmhists[11]->Fill(dtf->pm[0]);
			xmhists[11]->Fill(dtf->pm[1]);

			//track pT
			if(k->pt<100. || pi->pt<100. || xpip->pt<100. || xpim->pt<100.) continue;
			++counts[12];
			bmhists[12]->Fill(dtf->pm[0]);
			xmhists[12]->Fill(dtf->pm[1]);
			//track p
			if(k->p<1500. || pi->p<1500. || xpip->p<1500. || xpim->p<1500.) continue;
			++counts[13];
			bmhists[13]->Fill(dtf->pm[0]);
			xmhists[13]->Fill(dtf->pm[1]);
			//track IP chi2 from PV
			if(k->ip_chi2_best<10. || pi->ip_chi2_best<10. || xpip->ip_chi2_best<10. || xpim->ip_chi2_best<10.) continue;
			++counts[14];
			bmhists[14]->Fill(dtf->pm[0]);
			xmhists[14]->Fill(dtf->pm[1]);
			//track IP chi2 from SV
			if(k->ip_chi2_bvtx>10. || pi->ip_chi2_bvtx>10.) continue;
			++counts[15];
			bmhists[15]->Fill(dtf->pm[0]);
			xmhists[15]->Fill(dtf->pm[1]);
			//track ghost prob
			if(k->prb_ghost>0.1 || pi->prb_ghost>0.1 || xpip->prb_ghost>0.1 || xpim->prb_ghost>0.1) continue;
			++counts[16];
			bmhists[16]->Fill(dtf->pm[0]);
			xmhists[16]->Fill(dtf->pm[1]);
			//track PID
			if(k->pnn_k<0.3) continue;
			++counts[17];
			bmhists[17]->Fill(dtf->pm[0]);
			xmhists[17]->Fill(dtf->pm[1]);

			//sum of IP chi2
			if(k->ip_chi2_best + pi->ip_chi2_best + xpip->ip_chi2_best + xpim->ip_chi2_best < 200.) continue;
			++counts[18];
			bmhists[18]->Fill(dtf->pm[0]);
			xmhists[18]->Fill(dtf->pm[1]);

			tout->Fill();
			++iCand;
		}

	}
	for (uint i=0; i<counts.size(); ++i) {
		printf("%20s % 5d %.2f %.2f\n", countNames[i].Data(), counts[i], (i>0?(static_cast<double>(counts[i])/counts[i-1]):1.),static_cast<double>(counts[i])/counts[0]);
	}
	gStyle->SetOptStat(0);
	TCanvas c;
	bmhists[0]->Draw();
	for (uint i=1; i<bmhists.size(); ++i) {
		bmhists[i]->Draw("same");
	}
	c.SaveAs("mB_skim.pdf");
	xmhists[0]->Draw();
	for (uint i=0; i<xmhists.size(); ++i) {
		xmhists[i]->Draw("same");
	}
	c.SaveAs("mX_skim.pdf");
	tout->Write();
	fout->Close();
}

int main(int argc, char** argv) {
	if(argc<2) {
		std::cout << "Usage: " << argv[0] << " <mode> [nmax=-1] [file=-1]" << std::endl;
		return 0;
	}
	int mode = atoi(argv[1]);
	int nmax(-1), file(-1);
	if(argc>2) nmax = atoi(argv[2]);
	if(argc>3) file = atoi(argv[3]);

	TreeSkimmer a(file);
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
