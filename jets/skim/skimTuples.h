//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  4 23:02:48 2017 by ROOT version 5.34/36
// from TTree data/data
// found on file: /eos/lhcb/user/d/dcraik/jets/361/0/output.root
//////////////////////////////////////////////////////////

#ifndef skimTuples_h
#define skimTuples_h

#include "velo.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include <map>
#include <vector>

//TODO//#include <boost/progress.hpp>

// Fixed size dimensions of array or collections stored in the TTree if any.

class skimTuples {
	public :
		enum JetTagType {
			NoTag,    //keep all jets
			NoJets,   //only keep SVs etc
			JetTag,   //tag on back-to-back jet with SV
			ZTag,  //tag on Z/g* (require jet that isn't a Z muon)
			JpsiTag,  //tag on back-to-back J/psi (no jet needed but keep if found)
			TruthTag, //find a true dijet
			TruthNoMuTag //find true jets and veto jets including muons from Z or W
		};

		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Declaration of leaf types
		std::vector<double>  *gen_idx_pvr;
		std::vector<double>  *gen_idx_jet;
		std::vector<double>  *gen_idx_prnt;
		std::vector<double>  *gen_pid;
		std::vector<double>  *gen_q;
		std::vector<double>  *gen_p;
		std::vector<double>  *gen_pt;
		std::vector<double>  *gen_px;
		std::vector<double>  *gen_py;
		std::vector<double>  *gen_pz;
		std::vector<double>  *gen_e;
		std::vector<double>  *gen_x;
		std::vector<double>  *gen_y;
		std::vector<double>  *gen_z;
		std::vector<double>  *gen_prnt_pid;
		std::vector<double>  *gen_res_pid;
		std::vector<double>  *gen_from_sig;
		std::vector<double>  *pvr_x;
		std::vector<double>  *pvr_y;
		std::vector<double>  *pvr_z;
		std::vector<double>  *pvr_dx;
		std::vector<double>  *pvr_dy;
		std::vector<double>  *pvr_dz;
		std::vector<double>  *pvr_chi2;
		std::vector<double>  *pvr_ndof;
		std::vector<double>  *svr_idx_pvr;
		std::vector<double>  *svr_idx_jet;
		std::vector<double>  *svr_idx_trk0;
		std::vector<double>  *svr_idx_trk1;
		std::vector<double>  *svr_idx_trk2;
		std::vector<double>  *svr_idx_trk3;
		std::vector<double>  *svr_idx_trk4;
		std::vector<double>  *svr_idx_trk5;
		std::vector<double>  *svr_idx_trk6;
		std::vector<double>  *svr_idx_trk7;
		std::vector<double>  *svr_idx_trk8;
		std::vector<double>  *svr_idx_trk9;
		std::vector<double>  *svr_p;
		std::vector<double>  *svr_pt;
		std::vector<double>  *svr_px;
		std::vector<double>  *svr_py;
		std::vector<double>  *svr_pz;
		std::vector<double>  *svr_e;
		std::vector<double>  *svr_x;
		std::vector<double>  *svr_y;
		std::vector<double>  *svr_z;
		std::vector<double>  *svr_dx;
		std::vector<double>  *svr_dy;
		std::vector<double>  *svr_dz;
		std::vector<double>  *svr_m;
		std::vector<double>  *svr_m_cor;
		std::vector<double>  *svr_m_cor_err;
		std::vector<double>  *svr_m_cor_err_full;
		std::vector<double>  *svr_fd_min;
		std::vector<double>  *svr_fd_chi2;
		std::vector<double>  *svr_chi2;
		std::vector<double>  *svr_ip_chi2;
		std::vector<double>  *svr_ip_chi2_sum;
		std::vector<double>  *svr_ip_chi2_min_trk;
		std::vector<double>  *svr_abs_q_sum;
		std::vector<double>  *svr_tau;
		std::vector<double>  *svr_ntrk;
		std::vector<double>  *svr_ntrk_jet;
		std::vector<double>  *svr_jet_dr;
		std::vector<double>  *svr_jet_pt;
		std::vector<double>  *svr_pass;
		std::vector<double>  *svr_bdt0;
		std::vector<double>  *svr_bdt1;
		std::vector<double>  *jet_idx_pvr;
		std::vector<double>  *jet_ntrk;
		std::vector<double>  *jet_nneu;
		std::vector<double>  *jet_p;
		std::vector<double>  *jet_pt;
		std::vector<double>  *jet_px;
		std::vector<double>  *jet_py;
		std::vector<double>  *jet_pz;
		std::vector<double>  *jet_e;
		std::vector<double>  *jet_l0_hadron_tos;
		std::vector<double>  *jet_l0_photon_tos;
		std::vector<double>  *jet_l0_electron_tos;
		std::vector<double>  *jet_l0_muon_tos;
		std::vector<double>  *jet_l0_dimuon_tos;
		std::vector<double>  *jet_hlt1_track_tos;
		std::vector<double>  *jet_hlt1_ditrack_tos;
		std::vector<double>  *trk_idx_gen;
		std::vector<double>  *trk_idx_pvr;
		std::vector<double>  *trk_idx_jet;
		std::vector<double>  *trk_p;
		std::vector<double>  *trk_pt;
		std::vector<double>  *trk_px;
		std::vector<double>  *trk_py;
		std::vector<double>  *trk_pz;
		std::vector<double>  *trk_e;
		std::vector<double>  *trk_pid;
		std::vector<double>  *trk_q;
		std::vector<double>  *trk_ip;
		std::vector<double>  *trk_ip_chi2;
		std::vector<double>  *trk_pnn_e;
		std::vector<double>  *trk_pnn_mu;
		std::vector<double>  *trk_pnn_pi;
		std::vector<double>  *trk_pnn_k;
		std::vector<double>  *trk_pnn_p;
		std::vector<double>  *trk_pnn_ghost;
		std::vector<double>  *trk_ecal;
		std::vector<double>  *trk_hcal;
		std::vector<double>  *trk_prb_ghost;
		std::vector<double>  *trk_type;
		std::vector<double>  *trk_is_mu;
		std::vector<double>  *trk_vid;
		std::vector<double>  *trk_x;
		std::vector<double>  *trk_y;
		std::vector<double>  *trk_z;
		std::vector<double>  *trk_chi2;
		std::vector<double>  *trk_ndof;
		std::vector<double>  *neu_idx_gen;
		std::vector<double>  *neu_idx_jet;
		std::vector<double>  *neu_p;
		std::vector<double>  *neu_pt;
		std::vector<double>  *neu_px;
		std::vector<double>  *neu_py;
		std::vector<double>  *neu_pz;
		std::vector<double>  *neu_e;
		std::vector<double>  *neu_pid;
		std::vector<double>  *z0_idx_pvr;
		std::vector<double>  *z0_idx_jet_trk0;
		std::vector<double>  *z0_idx_jet_trk1;
		std::vector<double>  *z0_idx_jet_dr;
		std::vector<double>  *z0_p;
		std::vector<double>  *z0_pt;
		std::vector<double>  *z0_px;
		std::vector<double>  *z0_py;
		std::vector<double>  *z0_pz;
		std::vector<double>  *z0_e;
		std::vector<double>  *z0_x;
		std::vector<double>  *z0_y;
		std::vector<double>  *z0_z;
		std::vector<double>  *z0_m;
		std::vector<double>  *z0_ip;
		std::vector<double>  *z0_ip_chi2;
		std::vector<double>  *z0_dr_jet;
		std::vector<double>  *z0_idx_trk0;
		std::vector<double>  *z0_idx_trk1;
		std::vector<double>  *z0_l0_muon_tos;
		std::vector<double>  *z0_l0_dimuon_tos;
		std::vector<double>  *z0_hlt1_hiptmu_tos;
		std::vector<double>  *z0_hlt2_hiptmu_tos;
		std::vector<double>  *d0_idx_pvr;
		std::vector<double>  *d0_idx_jet;
		std::vector<double>  *d0_p;
		std::vector<double>  *d0_pt;
		std::vector<double>  *d0_px;
		std::vector<double>  *d0_py;
		std::vector<double>  *d0_pz;
		std::vector<double>  *d0_e;
		std::vector<double>  *d0_x;
		std::vector<double>  *d0_y;
		std::vector<double>  *d0_z;
		std::vector<double>  *d0_m;
		std::vector<double>  *d0_ip;
		std::vector<double>  *d0_ip_chi2;
		std::vector<double>  *d0_vtx_chi2;
		std::vector<double>  *d0_vtx_ndof;
		std::vector<double>  *d0_fd;
		std::vector<double>  *d0_fd_chi2;
		std::vector<double>  *d0_tau;
		std::vector<double>  *d0_tau_err;
		std::vector<double>  *d0_tau_chi2;
		std::vector<double>  *d0_ntrk_jet;
		std::vector<double>  *d0_idx_trk0;
		std::vector<double>  *d0_idx_trk1;
		std::vector<double>  *dp_idx_pvr;
		std::vector<double>  *dp_idx_jet;
		std::vector<double>  *dp_p;
		std::vector<double>  *dp_pt;
		std::vector<double>  *dp_px;
		std::vector<double>  *dp_py;
		std::vector<double>  *dp_pz;
		std::vector<double>  *dp_e;
		std::vector<double>  *dp_x;
		std::vector<double>  *dp_y;
		std::vector<double>  *dp_z;
		std::vector<double>  *dp_m;
		std::vector<double>  *dp_ip;
		std::vector<double>  *dp_ip_chi2;
		std::vector<double>  *dp_vtx_chi2;
		std::vector<double>  *dp_vtx_ndof;
		std::vector<double>  *dp_fd;
		std::vector<double>  *dp_fd_chi2;
		std::vector<double>  *dp_tau;
		std::vector<double>  *dp_tau_err;
		std::vector<double>  *dp_tau_chi2;
		std::vector<double>  *dp_ntrk_jet;
		std::vector<double>  *dp_idx_trk0;
		std::vector<double>  *dp_idx_trk1;
		std::vector<double>  *dp_idx_trk2;
		std::vector<double>  *ds_idx_pvr;
		std::vector<double>  *ds_idx_jet;
		std::vector<double>  *ds_p;
		std::vector<double>  *ds_pt;
		std::vector<double>  *ds_px;
		std::vector<double>  *ds_py;
		std::vector<double>  *ds_pz;
		std::vector<double>  *ds_e;
		std::vector<double>  *ds_x;
		std::vector<double>  *ds_y;
		std::vector<double>  *ds_z;
		std::vector<double>  *ds_m;
		std::vector<double>  *ds_ip;
		std::vector<double>  *ds_ip_chi2;
		std::vector<double>  *ds_vtx_chi2;
		std::vector<double>  *ds_vtx_ndof;
		std::vector<double>  *ds_fd;
		std::vector<double>  *ds_fd_chi2;
		std::vector<double>  *ds_tau;
		std::vector<double>  *ds_tau_err;
		std::vector<double>  *ds_tau_chi2;
		std::vector<double>  *ds_ntrk_jet;
		std::vector<double>  *ds_idx_trk0;
		std::vector<double>  *ds_idx_trk1;
		std::vector<double>  *ds_idx_trk2;
		std::vector<double>  *lc_idx_pvr;
		std::vector<double>  *lc_idx_jet;
		std::vector<double>  *lc_p;
		std::vector<double>  *lc_pt;
		std::vector<double>  *lc_px;
		std::vector<double>  *lc_py;
		std::vector<double>  *lc_pz;
		std::vector<double>  *lc_e;
		std::vector<double>  *lc_x;
		std::vector<double>  *lc_y;
		std::vector<double>  *lc_z;
		std::vector<double>  *lc_m;
		std::vector<double>  *lc_ip;
		std::vector<double>  *lc_ip_chi2;
		std::vector<double>  *lc_vtx_chi2;
		std::vector<double>  *lc_vtx_ndof;
		std::vector<double>  *lc_fd;
		std::vector<double>  *lc_fd_chi2;
		std::vector<double>  *lc_tau;
		std::vector<double>  *lc_tau_err;
		std::vector<double>  *lc_tau_chi2;
		std::vector<double>  *lc_ntrk_jet;
		std::vector<double>  *lc_idx_trk0;
		std::vector<double>  *lc_idx_trk1;
		std::vector<double>  *lc_idx_trk2;
		std::vector<double>  *k3pi_idx_pvr;
		std::vector<double>  *k3pi_idx_jet;
		std::vector<double>  *k3pi_p;
		std::vector<double>  *k3pi_pt;
		std::vector<double>  *k3pi_px;
		std::vector<double>  *k3pi_py;
		std::vector<double>  *k3pi_pz;
		std::vector<double>  *k3pi_e;
		std::vector<double>  *k3pi_x;
		std::vector<double>  *k3pi_y;
		std::vector<double>  *k3pi_z;
		std::vector<double>  *k3pi_m;
		std::vector<double>  *k3pi_ip;
		std::vector<double>  *k3pi_ip_chi2;
		std::vector<double>  *k3pi_vtx_chi2;
		std::vector<double>  *k3pi_vtx_ndof;
		std::vector<double>  *k3pi_fd;
		std::vector<double>  *k3pi_fd_chi2;
		std::vector<double>  *k3pi_tau;
		std::vector<double>  *k3pi_tau_err;
		std::vector<double>  *k3pi_tau_chi2;
		std::vector<double>  *k3pi_ntrk_jet;
		std::vector<double>  *k3pi_idx_trk0;
		std::vector<double>  *k3pi_idx_trk1;
		std::vector<double>  *k3pi_idx_trk2;
		std::vector<double>  *k3pi_idx_trk3;
		std::vector<double>  *jpsi_idx_pvr;
		std::vector<double>  *jpsi_idx_jet;
		std::vector<double>  *jpsi_p;
		std::vector<double>  *jpsi_pt;
		std::vector<double>  *jpsi_px;
		std::vector<double>  *jpsi_py;
		std::vector<double>  *jpsi_pz;
		std::vector<double>  *jpsi_e;
		std::vector<double>  *jpsi_x;
		std::vector<double>  *jpsi_y;
		std::vector<double>  *jpsi_z;
		std::vector<double>  *jpsi_m;
		std::vector<double>  *jpsi_ip;
		std::vector<double>  *jpsi_ip_chi2;
		std::vector<double>  *jpsi_vtx_chi2;
		std::vector<double>  *jpsi_vtx_ndof;
		std::vector<double>  *jpsi_fd;
		std::vector<double>  *jpsi_fd_chi2;
		std::vector<double>  *jpsi_tau;
		std::vector<double>  *jpsi_tau_err;
		std::vector<double>  *jpsi_tau_chi2;
		std::vector<double>  *jpsi_ntrk_jet;
		std::vector<double>  *jpsi_idx_trk0;
		std::vector<double>  *jpsi_idx_trk1;
		std::vector<double>  *evt_dec;
		std::vector<double>  *evt_j1_idx;
		std::vector<double>  *evt_j1_dR;
		std::vector<double>  *evt_j1_nsv;
		std::vector<double>  *evt_j1_nmu;
		std::vector<double>  *evt_j1_ntrk;
		std::vector<double>  *evt_j1_nneu;
		std::vector<double>  *evt_j1_p;
		std::vector<double>  *evt_j1_pt;
		std::vector<double>  *evt_j1_px;
		std::vector<double>  *evt_j1_py;
		std::vector<double>  *evt_j1_pz;
		std::vector<double>  *evt_j1_e;
		std::vector<double>  *evt_j2_idx;
		std::vector<double>  *evt_j2_dR;
		std::vector<double>  *evt_j2_nsv;
		std::vector<double>  *evt_j2_nmu;
		std::vector<double>  *evt_j2_ntrk;
		std::vector<double>  *evt_j2_nneu;
		std::vector<double>  *evt_j2_p;
		std::vector<double>  *evt_j2_pt;
		std::vector<double>  *evt_j2_px;
		std::vector<double>  *evt_j2_py;
		std::vector<double>  *evt_j2_pz;
		std::vector<double>  *evt_j2_e;
		Double_t        evt_pvr_n;
		Double_t        evt_trk_n;

		// List of branches
		TBranch        *b_gen_idx_pvr;   //!
		TBranch        *b_gen_idx_jet;   //!
		TBranch        *b_gen_idx_prnt;   //!
		TBranch        *b_gen_pid;   //!
		TBranch        *b_gen_q;   //!
		TBranch        *b_gen_p;   //!
		TBranch        *b_gen_pt;   //!
		TBranch        *b_gen_px;   //!
		TBranch        *b_gen_py;   //!
		TBranch        *b_gen_pz;   //!
		TBranch        *b_gen_e;   //!
		TBranch        *b_gen_x;   //!
		TBranch        *b_gen_y;   //!
		TBranch        *b_gen_z;   //!
		TBranch        *b_gen_prnt_pid;   //!
		TBranch        *b_gen_res_pid;   //!
		TBranch        *b_gen_from_sig;   //!
		TBranch        *b_pvr_x;   //!
		TBranch        *b_pvr_y;   //!
		TBranch        *b_pvr_z;   //!
		TBranch        *b_pvr_dx;   //!
		TBranch        *b_pvr_dy;   //!
		TBranch        *b_pvr_dz;   //!
		TBranch        *b_pvr_chi2;   //!
		TBranch        *b_pvr_ndof;   //!
		TBranch        *b_svr_idx_pvr;   //!
		TBranch        *b_svr_idx_jet;   //!
		TBranch        *b_svr_idx_trk0;   //!
		TBranch        *b_svr_idx_trk1;   //!
		TBranch        *b_svr_idx_trk2;   //!
		TBranch        *b_svr_idx_trk3;   //!
		TBranch        *b_svr_idx_trk4;   //!
		TBranch        *b_svr_idx_trk5;   //!
		TBranch        *b_svr_idx_trk6;   //!
		TBranch        *b_svr_idx_trk7;   //!
		TBranch        *b_svr_idx_trk8;   //!
		TBranch        *b_svr_idx_trk9;   //!
		TBranch        *b_svr_p;   //!
		TBranch        *b_svr_pt;   //!
		TBranch        *b_svr_px;   //!
		TBranch        *b_svr_py;   //!
		TBranch        *b_svr_pz;   //!
		TBranch        *b_svr_e;   //!
		TBranch        *b_svr_x;   //!
		TBranch        *b_svr_y;   //!
		TBranch        *b_svr_z;   //!
		TBranch        *b_svr_dx;   //!
		TBranch        *b_svr_dy;   //!
		TBranch        *b_svr_dz;   //!
		TBranch        *b_svr_m;   //!
		TBranch        *b_svr_m_cor;   //!
		TBranch        *b_svr_m_cor_err;   //!
		TBranch        *b_svr_m_cor_err_full;   //!
		TBranch        *b_svr_fd_min;   //!
		TBranch        *b_svr_fd_chi2;   //!
		TBranch        *b_svr_chi2;   //!
		TBranch        *b_svr_ip_chi2;   //!
		TBranch        *b_svr_ip_chi2_sum;   //!
		TBranch        *b_svr_ip_chi2_min_trk;
		TBranch        *b_svr_abs_q_sum;   //!
		TBranch        *b_svr_tau;   //!
		TBranch        *b_svr_ntrk;   //!
		TBranch        *b_svr_ntrk_jet;   //!
		TBranch        *b_svr_jet_dr;   //!
		TBranch        *b_svr_jet_pt;   //!
		TBranch        *b_svr_pass;   //!
		TBranch        *b_svr_bdt0;   //!
		TBranch        *b_svr_bdt1;   //!
		TBranch        *b_jet_idx_pvr;   //!
		TBranch        *b_jet_ntrk;   //!
		TBranch        *b_jet_nneu;   //!
		TBranch        *b_jet_p;   //!
		TBranch        *b_jet_pt;   //!
		TBranch        *b_jet_px;   //!
		TBranch        *b_jet_py;   //!
		TBranch        *b_jet_pz;   //!
		TBranch        *b_jet_e;   //!
		TBranch        *b_jet_l0_hadron_tos;
		TBranch        *b_jet_l0_photon_tos;
		TBranch        *b_jet_l0_electron_tos;
		TBranch        *b_jet_l0_muon_tos;
		TBranch        *b_jet_l0_dimuon_tos;
		TBranch        *b_jet_hlt1_track_tos;
		TBranch        *b_jet_hlt1_ditrack_tos;
		TBranch        *b_trk_idx_gen;   //!
		TBranch        *b_trk_idx_pvr;   //!
		TBranch        *b_trk_idx_jet;   //!
		TBranch        *b_trk_p;   //!
		TBranch        *b_trk_pt;   //!
		TBranch        *b_trk_px;   //!
		TBranch        *b_trk_py;   //!
		TBranch        *b_trk_pz;   //!
		TBranch        *b_trk_e;   //!
		TBranch        *b_trk_pid;   //!
		TBranch        *b_trk_q;   //!
		TBranch        *b_trk_ip;   //!
		TBranch        *b_trk_ip_chi2;   //!
		TBranch        *b_trk_pnn_e;   //!
		TBranch        *b_trk_pnn_mu;   //!
		TBranch        *b_trk_pnn_pi;   //!
		TBranch        *b_trk_pnn_k;   //!
		TBranch        *b_trk_pnn_p;   //!
		TBranch        *b_trk_pnn_ghost;   //!
		TBranch        *b_trk_ecal;   //!
		TBranch        *b_trk_hcal;   //!
		TBranch        *b_trk_prb_ghost;   //!
		TBranch        *b_trk_type;   //!
		TBranch        *b_trk_is_mu;   //!
		TBranch        *b_trk_vid;   //!
		TBranch        *b_trk_x;   //!
		TBranch        *b_trk_y;   //!
		TBranch        *b_trk_z;   //!
		TBranch        *b_trk_chi2;   //!
		TBranch        *b_trk_ndof;   //!
		TBranch        *b_neu_idx_gen;   //!
		TBranch        *b_neu_idx_jet;   //!
		TBranch        *b_neu_p;   //!
		TBranch        *b_neu_pt;   //!
		TBranch        *b_neu_px;   //!
		TBranch        *b_neu_py;   //!
		TBranch        *b_neu_pz;   //!
		TBranch        *b_neu_e;   //!
		TBranch        *b_neu_pid;   //!
		TBranch        *b_z0_idx_pvr;   //!
		TBranch        *b_z0_idx_jet_trk0;   //!
		TBranch        *b_z0_idx_jet_trk1;   //!
		TBranch        *b_z0_idx_jet_dr;   //!
		TBranch        *b_z0_p;   //!
		TBranch        *b_z0_pt;   //!
		TBranch        *b_z0_px;   //!
		TBranch        *b_z0_py;   //!
		TBranch        *b_z0_pz;   //!
		TBranch        *b_z0_e;   //!
		TBranch        *b_z0_x;   //!
		TBranch        *b_z0_y;   //!
		TBranch        *b_z0_z;   //!
		TBranch        *b_z0_m;   //!
		TBranch        *b_z0_ip;   //!
		TBranch        *b_z0_ip_chi2;   //!
		TBranch        *b_z0_dr_jet;   //!
		TBranch        *b_z0_idx_trk0;   //!
		TBranch        *b_z0_idx_trk1;   //!
		TBranch        *b_z0_l0_muon_tos;
		TBranch        *b_z0_l0_dimuon_tos;
		TBranch        *b_z0_hlt1_hiptmu_tos;
		TBranch        *b_z0_hlt2_hiptmu_tos;
		TBranch        *b_d0_idx_pvr;   //!
		TBranch        *b_d0_idx_jet;   //!
		TBranch        *b_d0_p;   //!
		TBranch        *b_d0_pt;   //!
		TBranch        *b_d0_px;   //!
		TBranch        *b_d0_py;   //!
		TBranch        *b_d0_pz;   //!
		TBranch        *b_d0_e;   //!
		TBranch        *b_d0_x;   //!
		TBranch        *b_d0_y;   //!
		TBranch        *b_d0_z;   //!
		TBranch        *b_d0_m;   //!
		TBranch        *b_d0_ip;   //!
		TBranch        *b_d0_ip_chi2;   //!
		TBranch        *b_d0_vtx_chi2;   //!
		TBranch        *b_d0_vtx_ndof;   //!
		TBranch        *b_d0_fd;   //!
		TBranch        *b_d0_fd_chi2;   //!
		TBranch        *b_d0_tau;   //!
		TBranch        *b_d0_tau_err;   //!
		TBranch        *b_d0_tau_chi2;   //!
		TBranch        *b_d0_ntrk_jet;   //!
		TBranch        *b_d0_idx_trk0;   //!
		TBranch        *b_d0_idx_trk1;   //!
		TBranch        *b_dp_idx_pvr;   //!
		TBranch        *b_dp_idx_jet;   //!
		TBranch        *b_dp_p;   //!
		TBranch        *b_dp_pt;   //!
		TBranch        *b_dp_px;   //!
		TBranch        *b_dp_py;   //!
		TBranch        *b_dp_pz;   //!
		TBranch        *b_dp_e;   //!
		TBranch        *b_dp_x;   //!
		TBranch        *b_dp_y;   //!
		TBranch        *b_dp_z;   //!
		TBranch        *b_dp_m;   //!
		TBranch        *b_dp_ip;   //!
		TBranch        *b_dp_ip_chi2;   //!
		TBranch        *b_dp_vtx_chi2;   //!
		TBranch        *b_dp_vtx_ndof;   //!
		TBranch        *b_dp_fd;   //!
		TBranch        *b_dp_fd_chi2;   //!
		TBranch        *b_dp_tau;   //!
		TBranch        *b_dp_tau_err;   //!
		TBranch        *b_dp_tau_chi2;   //!
		TBranch        *b_dp_ntrk_jet;   //!
		TBranch        *b_dp_idx_trk0;   //!
		TBranch        *b_dp_idx_trk1;   //!
		TBranch        *b_dp_idx_trk2;   //!
		TBranch        *b_ds_idx_pvr;   //!
		TBranch        *b_ds_idx_jet;   //!
		TBranch        *b_ds_p;   //!
		TBranch        *b_ds_pt;   //!
		TBranch        *b_ds_px;   //!
		TBranch        *b_ds_py;   //!
		TBranch        *b_ds_pz;   //!
		TBranch        *b_ds_e;   //!
		TBranch        *b_ds_x;   //!
		TBranch        *b_ds_y;   //!
		TBranch        *b_ds_z;   //!
		TBranch        *b_ds_m;   //!
		TBranch        *b_ds_ip;   //!
		TBranch        *b_ds_ip_chi2;   //!
		TBranch        *b_ds_vtx_chi2;   //!
		TBranch        *b_ds_vtx_ndof;   //!
		TBranch        *b_ds_fd;   //!
		TBranch        *b_ds_fd_chi2;   //!
		TBranch        *b_ds_tau;   //!
		TBranch        *b_ds_tau_err;   //!
		TBranch        *b_ds_tau_chi2;   //!
		TBranch        *b_ds_ntrk_jet;   //!
		TBranch        *b_ds_idx_trk0;   //!
		TBranch        *b_ds_idx_trk1;   //!
		TBranch        *b_ds_idx_trk2;   //!
		TBranch        *b_lc_idx_pvr;   //!
		TBranch        *b_lc_idx_jet;   //!
		TBranch        *b_lc_p;   //!
		TBranch        *b_lc_pt;   //!
		TBranch        *b_lc_px;   //!
		TBranch        *b_lc_py;   //!
		TBranch        *b_lc_pz;   //!
		TBranch        *b_lc_e;   //!
		TBranch        *b_lc_x;   //!
		TBranch        *b_lc_y;   //!
		TBranch        *b_lc_z;   //!
		TBranch        *b_lc_m;   //!
		TBranch        *b_lc_ip;   //!
		TBranch        *b_lc_ip_chi2;   //!
		TBranch        *b_lc_vtx_chi2;   //!
		TBranch        *b_lc_vtx_ndof;   //!
		TBranch        *b_lc_fd;   //!
		TBranch        *b_lc_fd_chi2;   //!
		TBranch        *b_lc_tau;   //!
		TBranch        *b_lc_tau_err;   //!
		TBranch        *b_lc_tau_chi2;   //!
		TBranch        *b_lc_ntrk_jet;   //!
		TBranch        *b_lc_idx_trk0;   //!
		TBranch        *b_lc_idx_trk1;   //!
		TBranch        *b_lc_idx_trk2;   //!
		TBranch        *b_k3pi_idx_pvr;   //!
		TBranch        *b_k3pi_idx_jet;   //!
		TBranch        *b_k3pi_p;   //!
		TBranch        *b_k3pi_pt;   //!
		TBranch        *b_k3pi_px;   //!
		TBranch        *b_k3pi_py;   //!
		TBranch        *b_k3pi_pz;   //!
		TBranch        *b_k3pi_e;   //!
		TBranch        *b_k3pi_x;   //!
		TBranch        *b_k3pi_y;   //!
		TBranch        *b_k3pi_z;   //!
		TBranch        *b_k3pi_m;   //!
		TBranch        *b_k3pi_ip;   //!
		TBranch        *b_k3pi_ip_chi2;   //!
		TBranch        *b_k3pi_vtx_chi2;   //!
		TBranch        *b_k3pi_vtx_ndof;   //!
		TBranch        *b_k3pi_fd;   //!
		TBranch        *b_k3pi_fd_chi2;   //!
		TBranch        *b_k3pi_tau;   //!
		TBranch        *b_k3pi_tau_err;   //!
		TBranch        *b_k3pi_tau_chi2;   //!
		TBranch        *b_k3pi_ntrk_jet;   //!
		TBranch        *b_k3pi_idx_trk0;   //!
		TBranch        *b_k3pi_idx_trk1;   //!
		TBranch        *b_k3pi_idx_trk2;   //!
		TBranch        *b_k3pi_idx_trk3;   //!
		TBranch        *b_jpsi_idx_pvr;   //!
		TBranch        *b_jpsi_idx_jet;   //!
		TBranch        *b_jpsi_p;   //!
		TBranch        *b_jpsi_pt;   //!
		TBranch        *b_jpsi_px;   //!
		TBranch        *b_jpsi_py;   //!
		TBranch        *b_jpsi_pz;   //!
		TBranch        *b_jpsi_e;   //!
		TBranch        *b_jpsi_x;   //!
		TBranch        *b_jpsi_y;   //!
		TBranch        *b_jpsi_z;   //!
		TBranch        *b_jpsi_m;   //!
		TBranch        *b_jpsi_ip;   //!
		TBranch        *b_jpsi_ip_chi2;   //!
		TBranch        *b_jpsi_vtx_chi2;   //!
		TBranch        *b_jpsi_vtx_ndof;   //!
		TBranch        *b_jpsi_fd;   //!
		TBranch        *b_jpsi_fd_chi2;   //!
		TBranch        *b_jpsi_tau;   //!
		TBranch        *b_jpsi_tau_err;   //!
		TBranch        *b_jpsi_tau_chi2;   //!
		TBranch        *b_jpsi_ntrk_jet;   //!
		TBranch        *b_jpsi_idx_trk0;   //!
		TBranch        *b_jpsi_idx_trk1;   //!
		TBranch        *b_evt_dec;   //!
		TBranch        *b_evt_j1_idx;   //!
		TBranch        *b_evt_j1_dR;   //!
		TBranch        *b_evt_j1_nsv;   //!
		TBranch        *b_evt_j1_nmu;   //!
		TBranch        *b_evt_j1_ntrk;   //!
		TBranch        *b_evt_j1_nneu;   //!
		TBranch        *b_evt_j1_p;   //!
		TBranch        *b_evt_j1_pt;   //!
		TBranch        *b_evt_j1_px;   //!
		TBranch        *b_evt_j1_py;   //!
		TBranch        *b_evt_j1_pz;   //!
		TBranch        *b_evt_j1_e;   //!
		TBranch        *b_evt_j2_idx;   //!
		TBranch        *b_evt_j2_dR;   //!
		TBranch        *b_evt_j2_nsv;   //!
		TBranch        *b_evt_j2_nmu;   //!
		TBranch        *b_evt_j2_ntrk;   //!
		TBranch        *b_evt_j2_nneu;   //!
		TBranch        *b_evt_j2_p;   //!
		TBranch        *b_evt_j2_pt;   //!
		TBranch        *b_evt_j2_px;   //!
		TBranch        *b_evt_j2_py;   //!
		TBranch        *b_evt_j2_pz;   //!
		TBranch        *b_evt_j2_e;   //!
		TBranch        *b_evt_pvr_n;   //!
		TBranch        *b_evt_trk_n;   //!

		std::vector<double> *svtrk[10];

		//output
		TFile* fout;
		TTree* tout;
		TTree* lumiout;

		int EVT, DEC, NPV;
		double JPX,JPY,JPZ,JE,JPT,JETA,JS1,JS2,JQ,JN,JNQ,JNN,JPTD,JTRIG10,JTRIG17,JTRIG60, JTRIGPT, JPV;
		double JTRUEPX,JTRUEPY,JTRUEPZ,JTRUEE,JTRUEPT,JTRUEETA, JTRUEDR;
		double JTRUEb, JTRUEc, JTRUED0, JTRUEDP, JTRUEDS, JTRUEJPSI, JTRUEDSV, JTRUEBSV, JTRUESV;
		double JTRUEDPX, JTRUEDPY, JTRUEDPZ, JTRUEDE, JTRUEDPT, JTRUEDDR;
		double TPX,TPY,TPZ,TE,TPT,TETA,TTRIG10,TTRIG17,TTRIG60, TTRIGPT;
		double TTRUEPX,TTRUEPY,TTRUEPZ,TTRUEE,TTRUEPT,TTRUEETA, TTRUEDR;
		double DR, DPHI;
		double ZM,ZP,ZPX,ZPY,ZPZ,ZPT,ZE,ZETA,ZY,ZDR;
		double MU0PX,MU0PY,MU0PZ,MU0PT,MU0IP,MU0IPCHI2,MU0FPT, MU0DR;
		double MU1PX,MU1PY,MU1PZ,MU1PT,MU1IP,MU1IPCHI2,MU1FPT, MU1DR;
		double ZTRUEM,ZTRUEP,ZTRUEPX,ZTRUEPY,ZTRUEPZ,ZTRUEPT,ZTRUEE,ZTRUEETA,ZTRUEDR;
		double MU0TRUEPX,MU0TRUEPY,MU0TRUEPZ,MU0TRUEPT;
		double MU1TRUEPX,MU1TRUEPY,MU1TRUEPZ,MU1TRUEPT;
		std::vector<double> *SVM;
		std::vector<double> *SVMCOR;
		std::vector<double> *SVMCORERR;
		std::vector<double> *SVMINPERP;
		std::vector<double> *SVPT;
		std::vector<double> *SVDRJ;
		std::vector<double> *SVDRT;
		std::vector<double> *SVN;
		std::vector<double> *SVNTRUE;
		std::vector<double> *SVNLINKED;
		std::vector<double> *SVNJ;
		std::vector<double> *SVNT;
		std::vector<double> *SVQ;
		std::vector<double> *SVPERP;
		std::vector<double> *SVETA;
		std::vector<double> *SVTZ;
		std::vector<double> *SVVXCHI2;
		std::vector<double> *SVIPCHI2;
		std::vector<double> *SVMINIPCHI2;
		std::vector<double> *SVPX;
		std::vector<double> *SVPY;
		std::vector<double> *SVPZ;
		std::vector<double> *SVE;
		std::vector<double> *SVMAXGHOST;
		std::vector<double> *SVX;
		std::vector<double> *SVY;
		std::vector<double> *SVZ;
		std::vector<double> *SVSUMIPCHI2;
		std::vector<double> *SVISD0;
		std::vector<double> *SVISDP;
		std::vector<double> *SVD0M;
		std::vector<double> *SVDPM;
		std::vector<double> *SVTRK0IDX;
		std::vector<double> *SVTRK1IDX;
		std::vector<double> *SVTRK2IDX;
		std::vector<double> *SVTRK3IDX;
		std::vector<double> *SVTRUEIDX;
		std::vector<double> *SVTRUETRK0IDX;
		std::vector<double> *SVTRUETRK1IDX;
		std::vector<double> *SVTRUETRK2IDX;
		std::vector<double> *SVTRUETRK3IDX;
		std::vector<double> *SVTRK0P;
		std::vector<double> *SVTRK1P;
		std::vector<double> *SVTRK2P;
		std::vector<double> *SVTRK3P;
		std::vector<double> *SVTRK0PT;
		std::vector<double> *SVTRK1PT;
		std::vector<double> *SVTRK2PT;
		std::vector<double> *SVTRK3PT;
		std::vector<double> *SVTRK0PX;
		std::vector<double> *SVTRK1PX;
		std::vector<double> *SVTRK2PX;
		std::vector<double> *SVTRK3PX;
		std::vector<double> *SVTRK0PY;
		std::vector<double> *SVTRK1PY;
		std::vector<double> *SVTRK2PY;
		std::vector<double> *SVTRK3PY;
		std::vector<double> *SVTRK0PZ;
		std::vector<double> *SVTRK1PZ;
		std::vector<double> *SVTRK2PZ;
		std::vector<double> *SVTRK3PZ;
		std::vector<double> *SVTRK0PNNPI;
		std::vector<double> *SVTRK1PNNPI;
		std::vector<double> *SVTRK2PNNPI;
		std::vector<double> *SVTRK3PNNPI;
		std::vector<double> *SVTRK0PNNK;
		std::vector<double> *SVTRK1PNNK;
		std::vector<double> *SVTRK2PNNK;
		std::vector<double> *SVTRK3PNNK;
		
		std::vector<double> *TSVM;
		std::vector<double> *TSVMCOR;
		std::vector<double> *TSVMCORERR;
		std::vector<double> *TSVMINPERP;
		std::vector<double> *TSVPT;
		std::vector<double> *TSVDRJ;
		std::vector<double> *TSVDRT;
		std::vector<double> *TSVN;
		std::vector<double> *TSVNTRUE;
		std::vector<double> *TSVNJ;
		std::vector<double> *TSVNT;
		std::vector<double> *TSVQ;
		std::vector<double> *TSVPERP;
		std::vector<double> *TSVETA;
		std::vector<double> *TSVTZ;
		std::vector<double> *TSVIPCHI2;
		std::vector<double> *TSVMINIPCHI2;
		std::vector<double> *TSVPX;
		std::vector<double> *TSVPY;
		std::vector<double> *TSVPZ;
		std::vector<double> *TSVE;
		std::vector<double> *TSVMAXGHOST;
		std::vector<double> *TSVX;
		std::vector<double> *TSVY;
		std::vector<double> *TSVZ;
		std::vector<double> *TSVSUMIPCHI2;
		std::vector<double> *TSVISD0;
		std::vector<double> *TSVISDP;
		std::vector<double> *TSVD0M;
		std::vector<double> *TSVDPM;
		std::vector<double> *TSVTRUEIDX;

		std::vector<double> *TRKPT;
		std::vector<double> *TRKIPCHI2;
		std::vector<double> *TRKINJET;

		double NSV,NSVTRK,NSVTRUETRK,NTSV,PVX,PVY,PVZ,NDISPL6,NDISPL9,NDISPL16;
		double MUPT,MUIPCHI2,MUDR,MUPNN,NMU;
		double HPT,HIPCHI2,HDR;
		double NTRK, NNEU;

		std::vector<double> *TRUEBID;
		std::vector<double> *TRUEBPX;
		std::vector<double> *TRUEBPY;
		std::vector<double> *TRUEBPZ;
		std::vector<double> *TRUEBPT;
		std::vector<double> *TRUEBE;
		std::vector<double> *TRUEBX;
		std::vector<double> *TRUEBY;
		std::vector<double> *TRUEBZ;
		std::vector<double> *TRUEDID;
		std::vector<double> *TRUEDPX;
		std::vector<double> *TRUEDPY;
		std::vector<double> *TRUEDPZ;
		std::vector<double> *TRUEDPT;
		std::vector<double> *TRUEDE;
		std::vector<double> *TRUEDX;
		std::vector<double> *TRUEDY;
		std::vector<double> *TRUEDZ;
		std::vector<double> *TRUEDFROMB;
		std::vector<double> *TRUEDTRUEB;
		std::vector<double> *TRUEDSEL;
		std::vector<double> *TRUEDTRK0IDX;
		std::vector<double> *TRUEDTRK0ID;
		std::vector<double> *TRUEDTRK0P;
		std::vector<double> *TRUEDTRK0PT;
		std::vector<double> *TRUEDTRK0INACC;
		std::vector<double> *TRUEDTRK0RECO;
		std::vector<double> *TRUEDTRK0PNNK;
		std::vector<double> *TRUEDTRK0PNNPI;
		std::vector<double> *TRUEDTRK1IDX;
		std::vector<double> *TRUEDTRK1ID;
		std::vector<double> *TRUEDTRK1P;
		std::vector<double> *TRUEDTRK1PT;
		std::vector<double> *TRUEDTRK1INACC;
		std::vector<double> *TRUEDTRK1RECO;
		std::vector<double> *TRUEDTRK1PNNK;
		std::vector<double> *TRUEDTRK1PNNPI;
		std::vector<double> *TRUEDTRK2IDX;
		std::vector<double> *TRUEDTRK2ID;
		std::vector<double> *TRUEDTRK2P;
		std::vector<double> *TRUEDTRK2PT;
		std::vector<double> *TRUEDTRK2INACC;
		std::vector<double> *TRUEDTRK2RECO;
		std::vector<double> *TRUEDTRK2PNNK;
		std::vector<double> *TRUEDTRK2PNNPI;
		std::vector<double> *TRUEDTRK3IDX;
		std::vector<double> *TRUEDTRK3ID;
		std::vector<double> *TRUEDTRK3P;
		std::vector<double> *TRUEDTRK3PT;
		std::vector<double> *TRUEDTRK3INACC;
		std::vector<double> *TRUEDTRK3RECO;
		std::vector<double> *TRUEDTRK3PNNK;
		std::vector<double> *TRUEDTRK3PNNPI;
		std::vector<double> *TRUEDTRK0PX;
		std::vector<double> *TRUEDTRK1PX;
		std::vector<double> *TRUEDTRK2PX;
		std::vector<double> *TRUEDTRK3PX;
		std::vector<double> *TRUEDTRK0PY;
		std::vector<double> *TRUEDTRK1PY;
		std::vector<double> *TRUEDTRK2PY;
		std::vector<double> *TRUEDTRK3PY;
		std::vector<double> *TRUEDTRK0PZ;
		std::vector<double> *TRUEDTRK1PZ;
		std::vector<double> *TRUEDTRK2PZ;
		std::vector<double> *TRUEDTRK3PZ;

		std::vector<double> *D0M;
		std::vector<double> *D0PX;
		std::vector<double> *D0PY;
		std::vector<double> *D0PZ;
		std::vector<double> *D0E;
		std::vector<double> *D0X;
		std::vector<double> *D0Y;
		std::vector<double> *D0Z;
		std::vector<double> *D0FD;
		std::vector<double> *D0FDCHI2;
		std::vector<double> *D0IP;
		std::vector<double> *D0IPCHI2;
		std::vector<double> *D0VTXCHI2;
		std::vector<double> *D0VTXNDOF;;
		std::vector<double> *D0DRJET;
		std::vector<double> *D0DRTAG;
		std::vector<double> *D0TAU;
		std::vector<double> *D0DIRA;
		std::vector<double> *D0PIPNNK;
		std::vector<double> *D0KPNNK;
		std::vector<double> *D0PIPNNPI;
		std::vector<double> *D0KPNNPI;
		std::vector<double> *D0PIWEIGHT;
		std::vector<double> *D0KWEIGHT;

		std::vector<double> *D0P;
		std::vector<double> *D0PT;
		std::vector<double> *D0ETA;
		std::vector<double> *D0PIP;
		std::vector<double> *D0KP;
		std::vector<double> *D0PIPT;
		std::vector<double> *D0KPT;
		std::vector<double> *D0PIETA;
		std::vector<double> *D0KETA;
		std::vector<double> *D0PIPX;
		std::vector<double> *D0KPX;
		std::vector<double> *D0PIPY;
		std::vector<double> *D0KPY;
		std::vector<double> *D0PIPZ;
		std::vector<double> *D0KPZ;
		std::vector<double> *D0PIIPCHI2;
		std::vector<double> *D0KIPCHI2;

		std::vector<double> *D0TRK0;
		std::vector<double> *D0TRK1;
		std::vector<double> *D0TRUETRK0;
		std::vector<double> *D0TRUETRK1;

		std::vector<double> *D0NJ;
		std::vector<double> *D0MAXDR;

		std::vector<double> *D0TRUE;
		std::vector<double> *D0TRUEDR;
		std::vector<double> *D0TRUEIDX;
		std::vector<double> *D0FROMB ;

		std::vector<double> *DPM;
		std::vector<double> *DPPX;
		std::vector<double> *DPPY;
		std::vector<double> *DPPZ;
		std::vector<double> *DPE;
		std::vector<double> *DPX;
		std::vector<double> *DPY;
		std::vector<double> *DPZ;
		std::vector<double> *DPFD;
		std::vector<double> *DPFDCHI2;
		std::vector<double> *DPIP;
		std::vector<double> *DPIPCHI2;
		std::vector<double> *DPVTXCHI2;
		std::vector<double> *DPVTXNDOF;;
		std::vector<double> *DPDRJET;
		std::vector<double> *DPDRTAG;
		std::vector<double> *DPTAU;
		std::vector<double> *DPDIRA;
		std::vector<double> *DPPI1PNN;
		std::vector<double> *DPPI2PNN;
		std::vector<double> *DPKPNN;

		std::vector<double> *DPP;
		std::vector<double> *DPPT;
		std::vector<double> *DPPI1PT;
		std::vector<double> *DPPI2PT;
		std::vector<double> *DPKPT;

		std::vector<double> *DPTRK0;
		std::vector<double> *DPTRK1;
		std::vector<double> *DPTRK2;

		std::vector<double> *DPNJ;
		std::vector<double> *DPMAXDR;

		std::vector<double> *DPTRUE;
		std::vector<double> *DPTRUEDR;
		std::vector<double> *DPTRUEIDX;
		std::vector<double> *DPFROMB ;

		std::vector<double> *DSM;
		std::vector<double> *DSPX;
		std::vector<double> *DSPY;
		std::vector<double> *DSPZ;
		std::vector<double> *DSE;
		std::vector<double> *DSX;
		std::vector<double> *DSY;
		std::vector<double> *DSZ;
		std::vector<double> *DSFD;
		std::vector<double> *DSFDCHI2;
		std::vector<double> *DSIP;
		std::vector<double> *DSIPCHI2;
		std::vector<double> *DSVTXCHI2;
		std::vector<double> *DSVTXNDOF;;
		std::vector<double> *DSDRJET;
		std::vector<double> *DSDRTAG;
		std::vector<double> *DSTAU;
		std::vector<double> *DSDIRA;
		std::vector<double> *DSPIPNN;
		std::vector<double> *DSK1PNN;
		std::vector<double> *DSK2PNN;
		std::vector<double> *DSPHIM;

		std::vector<double> *DSP;
		std::vector<double> *DSPT;
		std::vector<double> *DSPIPT;
		std::vector<double> *DSK1PT;
		std::vector<double> *DSK2PT;
		std::vector<double> *DSPHIPT;

		std::vector<double> *DSTRK0;
		std::vector<double> *DSTRK1;
		std::vector<double> *DSTRK2;

		std::vector<double> *DSNJ;
		std::vector<double> *DSMAXDR;

		std::vector<double> *DSTRUE;
		std::vector<double> *DSTRUEDR;
		std::vector<double> *DSTRUEIDX;
		std::vector<double> *DSFROMB ;

		std::vector<double> *LCM;
		std::vector<double> *LCPX;
		std::vector<double> *LCPY;
		std::vector<double> *LCPZ;
		std::vector<double> *LCE;
		std::vector<double> *LCX;
		std::vector<double> *LCY;
		std::vector<double> *LCZ;
		std::vector<double> *LCFD;
		std::vector<double> *LCFDCHI2;
		std::vector<double> *LCIP;
		std::vector<double> *LCIPCHI2;
		std::vector<double> *LCVTXCHI2;
		std::vector<double> *LCVTXNDOF;;
		std::vector<double> *LCDRJET;
		std::vector<double> *LCDRTAG;
		std::vector<double> *LCTAU;
		std::vector<double> *LCDIRA;
		std::vector<double> *LCPPNN;
		std::vector<double> *LCKPNN;
		std::vector<double> *LCPIPNN;

		std::vector<double> *LCP;
		std::vector<double> *LCPT;
		std::vector<double> *LCPPT;
		std::vector<double> *LCKPT;
		std::vector<double> *LCPIPT;

		std::vector<double> *LCTRK0;
		std::vector<double> *LCTRK1;
		std::vector<double> *LCTRK2;

		std::vector<double> *LCNJ;
		std::vector<double> *LCMAXDR;

		std::vector<double> *K3PIM;
		std::vector<double> *K3PIPX;
		std::vector<double> *K3PIPY;
		std::vector<double> *K3PIPZ;
		std::vector<double> *K3PIE;
		std::vector<double> *K3PIX;
		std::vector<double> *K3PIY;
		std::vector<double> *K3PIZ;
		std::vector<double> *K3PIFD;
		std::vector<double> *K3PIFDCHI2;
		std::vector<double> *K3PIIP;
		std::vector<double> *K3PIIPCHI2;
		std::vector<double> *K3PIVTXCHI2;
		std::vector<double> *K3PIVTXNDOF;;
		std::vector<double> *K3PIDRJET;
		std::vector<double> *K3PIDRTAG;
		std::vector<double> *K3PITAU;
		std::vector<double> *K3PIDIRA;
		std::vector<double> *K3PIKPNN;
		std::vector<double> *K3PIPI1PNN;
		std::vector<double> *K3PIPI2PNN;
		std::vector<double> *K3PIPI3PNN;

		std::vector<double> *K3PIP;
		std::vector<double> *K3PIPT;
		std::vector<double> *K3PIPI1PT;
		std::vector<double> *K3PIPI2PT;
		std::vector<double> *K3PIPI3PT;
		std::vector<double> *K3PIKPT;

		std::vector<double> *K3PITRK0;
		std::vector<double> *K3PITRK1;
		std::vector<double> *K3PITRK2;
		std::vector<double> *K3PITRK3;

		std::vector<double> *K3PINJ;
		std::vector<double> *K3PIMAXDR;

		std::vector<double> *JPSIM;
		std::vector<double> *JPSIPX;
		std::vector<double> *JPSIPY;
		std::vector<double> *JPSIPZ;
		std::vector<double> *JPSIE;
		std::vector<double> *JPSIX;
		std::vector<double> *JPSIY;
		std::vector<double> *JPSIZ;
		std::vector<double> *JPSIFD;
		std::vector<double> *JPSIFDCHI2;
		std::vector<double> *JPSIIP;
		std::vector<double> *JPSIIPCHI2;
		std::vector<double> *JPSIVTXCHI2;
		std::vector<double> *JPSIVTXNDOF;;
		std::vector<double> *JPSIDRJET;
		std::vector<double> *JPSIDRTAG;
		std::vector<double> *JPSITAU;
		std::vector<double> *JPSIDIRA;

		std::vector<double> *JPSIP;
		std::vector<double> *JPSIPT;
		std::vector<double> *JPSIPIPT;
		std::vector<double> *JPSIKPT;

		std::vector<double> *JPSITRK0;
		std::vector<double> *JPSITRK1;

		std::vector<double> *JPSINJ;
		std::vector<double> *JPSIMAXDR;

		std::vector<double> *JPSITRUE;
		std::vector<double> *JPSITRUEDR;
		std::vector<double> *JPSITRUEIDX;
		std::vector<double> *JPSIFROMB ;

		skimTuples(int year, int part=-1, TString dir=".");
		virtual ~skimTuples();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop(int nmax=-1);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		void clearEventOutputs();
		void clearOutputs();
		void clearOutputVectors();

		double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
				 	     const TVector3 &pb, const TVector3 &db);

		void pairPVs();
		int getTruePV(int pv);
		int getRecoPV(int pv);

		//check gen particle at idxD matches requirements and return indices of
		//tracks to idcs. Default values override requirements. Decays with 
		//final-state particles other than pi,K,p,mu always return false unless 
		//allowPart is set.
		bool checkTruth(int idxD, std::vector<int>& idcs, int needPID=0, int needPi=-1, int needK=-1, int needP=-1, int needMu=-1, int needOppPi=-1, int needOppK=-1, int needOppP=-1, bool allowPart=false);
		bool checkD0Truth(int idxD);
		bool checkD0Truth(int idxK,int idxPi, int& idxD, bool ordered=true);
		bool checkDpTruth(int idxD);
		bool checkDpTruth(int idxK,int idxPi1,int idxPi2, int& idxD, bool ordered=true);
		bool checkDsTruth(int idxD);
		bool checkDsTruth(int idxK1,int idxK2,int idxPi, int& idxD, bool ordered=true);
		bool checkJpsiTruth(int idxJpsi);
		bool checkJpsiTruth(int idxMu1,int idxMu2, int& idxJpsi);
		bool checkDSVTruth(int idxD, double minPT=5000.);
		bool checkBSVTruth(int idxD, double minPT=5000.);
		bool checkSVTruth(int idxD, double minPT=5000.);
		bool checkFromB(int idxD);

		int getFlavour(int pid);

		bool longLivedB(int pid);
		bool longLivedC(int pid);
		bool longLivedS(int pid);
		bool longLived(int pid);
		int checkRealSV(std::vector<int> indices, int* common=0);

		int checkBestD0Cand(int dA, int dB);

		void fillPVEffs(int s, int j);
		int getSVCategory(int s, int j, std::vector<int> indices);

		int matchJetTruth(int j, int pid=98, bool matchCharge=true, bool useHiPt=true);

		//tagging functions
		bool tagTruthJet(int& j1, int& j2);
		bool tagTruthNoMuJet(int& j1, int& z);
		bool tagSVJet(int& j, int& t);
		bool tagZJet(int& j, int& z);
		bool tagJpsi(int& j, int& t);

		//filling functions
		void fillZ(int z, int j);
		void fillTrueZ(int j);
		void fillOutput(int j, int t);
		void fillTruthOutput();
		void fillSVCands(int j, int t);//fill the output for SV, D and Jpsi containers

		//function to fix breaks in generated particle trees (caused when the D from a b->c decay is saved before the B)
		void reattachLostParticles();

		//determine which sample we're unning over
		int year_;
		int part_;

		//type of jet tagging to use
		JetTagType tagType_;

		//whether to use backwards SVs
		bool backwards_;

		//flavour and pt range for MC jets
		bool isMC_;
		int flavour_;
		double minTruePt_;
		double maxTruePt_;

		//adjust the SV requirements
		double minipchi2cut_;
		double maxmcorerrcut_;

		TH2F* pidPi;
		TH2F* pidK;

		VeloMaterial* detector;
		double detOffsetX;
		double detOffsetY;

		std::set<unsigned int> usedTruthJets_;

		std::map<int,int> lookupTruePV;
		std::map<int,int> lookupRecoPV;
		bool pvsPaired; //track whether PV matching has been performed for current event

		TString outName;
		TChain* lumi;

		int* tagCounts;
		//  0 - n events
		//  1 - n single PVs
		//  2 - single trigger / single good MC jet
		//  3 - multiple triggers / two good MC jets
		//  4 - same jet triggers / multiple good MC jets
		//  5 - SV tagged / single good MC dijet
		//  6 - double SV tag / multiple good MC dijets
		//  7 - probe in trigger / best dijet has jet missing
		//  8 - good probe found / best dijet has both jets
		//  9 - multiple probe found
		// 10 - dijet match
		// 11 - multiple dijets - switched
		// 12 - multiple dijets - different
		// 13 - multiple dijets - both switched and different
		// 14 - multiple dijets - actually the same
		int* svTagCounts;
		//  0 - N have true c
		//  1 - N have true Hc
		//  2 - N have true Hc (pT>5GeV)
		//  3 - N have true SV
		//  4 - N have reco SV
		//  5 - N pass dR cut
		//  6 - N pass SV selection
		//  7 - N pass SVN<=4
		//  8 - N pass SVNJ>0
		TString svTagLog;
		TString histName;
		TH2D* svCatHist_;

		int* nPVCount;
		// 0 - N total true PVs
		// 1 - N no truth information
		// 2 - N matched PVs of jet
		// 3 - N matched PVs of SV
		// 4 - N matched PVs of majority of jet tracks
		// 5 - N matched PVs of majority of SV tracks
		// 6 - N matched PVs of pT-weighted majority of jet tracks
		TString pvEffLog;//log file
		bool pvEffsFilled;//flag so only one SV/jet pair is used per event
		//(guard against clones or subsets of SVs)
};

#endif

#ifdef skimTuples_cxx
skimTuples::skimTuples(int year, int part, TString dir) : fChain(0), year_(year), part_(part), backwards_(false), isMC_(false), flavour_(0), pvsPaired(false), pvEffsFilled(false)
{
	outName=dir;
	outName+="/";
	outName+=year;
	outName+="/";
	outName+=part;
	outName+="/skimmed.root";

	svTagLog=dir;
	svTagLog+="/";
	svTagLog+=year;
	svTagLog+="/";
	svTagLog+=part;
	svTagLog+="/svTagCounts.log";

	pvEffLog=dir;
	pvEffLog+="/";
	pvEffLog+=year;
	pvEffLog+="/";
	pvEffLog+=part;
	pvEffLog+="/pvEffCounts.log";

	histName=dir;
	histName+="/";
	histName+=year;
	histName+="/";
	histName+=part;
	histName+="/svCats.root";
	svCatHist_= new TH2D("svCats","",8,-0.5,7.5,10,-0.5,9.5);

	if(part==-2) outName="testSkim.root";

	TFile* fpidk = TFile::Open("../efficiencies/pidcalib/PerfHists_K_Turbo16_MagDown_kaons4_Brunel_P_Brunel_PT.root");
	TFile* fpidpi= TFile::Open("../efficiencies/pidcalib/PerfHists_Pi_Turbo16_MagDown_pions2_Brunel_P_Brunel_PT.root");

	pidK = dynamic_cast<TH2F*>(fpidk->Get("K_Brunel_MC15TuneV1_ProbNNK > 0.2_All"));
	pidPi= dynamic_cast<TH2F*>(fpidpi->Get("Pi_Brunel_MC15TuneV1_ProbNNpi > 0.1_All"));

	detector = new VeloMaterial("run2.root");
	detOffsetX = 0.;//data  0.818);//MC
	detOffsetY = 0.;//data -0.173);//MC

	bool isMC(false);//if MC then reset the detector offsets

	minTruePt_=0.;
	maxTruePt_=100000.;

	minipchi2cut_=9.;
	maxmcorerrcut_=500.;

	tagCounts = new int[15];
	for(int i=0; i<15; ++i) tagCounts[i]=0;
	svTagCounts = new int[15];
	for(int i=0; i<15; ++i) svTagCounts[i]=0;
	nPVCount = new int[7];
	for(int i=0; i<7; ++i) nPVCount[i]=0;

	char str[256];
	TChain *t = new TChain("data");
	lumi = new TChain("GetIntegratedLuminosity/LumiTuple");
	//select tagging type based on dataset number
	switch(year) {
		case 100:
		case 101:
			tagType_ = JetTag;
			break;
		case 102:
			tagType_ = JpsiTag;
			break;
		case 15://TODO test run
		case 105:
		case 106:
		case 107:
		case 108:
			tagType_ = ZTag;
			break;
		case 116:
		case 117:
		case 118:
			tagType_ = JetTag;
			break;
		case 126:
		case 127:
		case 128:
			tagType_ = JpsiTag;
			break;
		case 135:
		case 136:
		case 137:
		case 138:
			tagType_ = ZTag;
			break;
		case 140:
		case 141:
		case 142:
		case 143:
		case 144:
			isMC_=true;
			tagType_ = TruthTag;
			flavour_=4;
			break;
		case 150:
		case 151:
		case 152:
		case 153:
		case 154:
			isMC_=true;
			tagType_ = TruthTag;
			flavour_=5;
			break;
		case 160:
			isMC_=true;
			tagType_ = TruthNoMuTag;
			flavour_=0;
			break;
		case 164:
			isMC_=true;
			tagType_ = TruthNoMuTag;
			flavour_=4;
			break;
		case 165:
			isMC_=true;
			tagType_ = TruthNoMuTag;
			flavour_=5;
			break;
		case 200:
		case 201:
		case 202:
		case 203:
		case 205:
		case 240:
		case 241:
		case 242:
		case 243:
		case 244:
		case 250:
		case 251:
		case 252:
		case 253:
		case 254:
			isMC_=true;
			tagType_ = NoTag;
			break;
		case 305:
		case 306:
		case 307:
		case 308:
			tagType_ = ZTag;
			break;
		case 316:
		case 317:
		case 318:
			tagType_ = JetTag;
			break;
		case 840:
		case 841:
		case 842:
		case 843:
		case 844:
		case 940:
		case 941:
		case 942:
		case 943:
		case 944:
		case 1040:
		case 1041:
		case 1042:
		case 1043:
		case 1044:
			isMC_=true;
			tagType_ = TruthTag;
			flavour_=4;
			break;
		case 850:
		case 851:
		case 852:
		case 853:
		case 854:
		case 950:
		case 951:
		case 952:
		case 953:
		case 954:
		case 1050:
		case 1051:
		case 1052:
		case 1053:
		case 1054:
			isMC_=true;
			tagType_ = TruthTag;
			flavour_=5;
			break;
		default:
			std::cout << "Unknown dataset: tagging not configured" << std::endl;
	}

	//any extra settings
	switch(year) {
		case 101:
			backwards_=true;
			break;
		case 140:
		case 150:
			minTruePt_=5000.;
			maxTruePt_=10000.;
			break;
		case 141:
		case 151:
			minTruePt_=10000.;
			maxTruePt_=15000.;
			break;
		case 142:
		case 152:
			minTruePt_=15000.;
			maxTruePt_=20000.;
			break;
		case 143:
		case 153:
			minTruePt_=20000.;
			maxTruePt_=50000.;
			break;
		case 144:
		case 154:
			minTruePt_=50000.;
			maxTruePt_=100000.;
			break;
		case 164:
		case 165:
			minTruePt_=0.;
			maxTruePt_=200000.;
			break;
		case 305:
		case 306:
		case 307:
		case 308:
		case 316:
		case 317:
		case 318:
			backwards_=true;
			break;
		case 840:
		case 850:
		case 940:
		case 950:
		case 1040:
		case 1050:
		     	minTruePt_=5000.;
		     	maxTruePt_=10000.;
		     	break;
		case 841:
		case 851:
		case 941:
		case 951:
		case 1041:
		case 1051:
		     	minTruePt_=10000.;
		     	maxTruePt_=15000.;
		     	break;
		case 842:
		case 852:
		case 942:
		case 952:
		case 1042:
		case 1052:
		     	minTruePt_=15000.;
		     	maxTruePt_=20000.;
		     	break;
		case 843:
		case 853:
		case 943:
		case 953:
		case 1043:
		case 1053:
		     	minTruePt_=20000.;
		     	maxTruePt_=50000.;
		     	break;
		case 844:
		case 854:
		case 944:
		case 954:
		case 1044:
		case 1054:
			minTruePt_=50000.;
			maxTruePt_=100000.;
			break;
	}
	if(part==-2) {
		t->Add("/tmp/dcraik/output.root");
		lumi->Add("/tmp/dcraik/LumiTuple.root");
	} else if(part>=0) {
		sprintf(str,"%s/%d/%d/output.root",dir.Data(),year,part);
		if(!gSystem->AccessPathName(str)) {
			t->Add(str);
		}
		sprintf(str,"%s/%d/%d/LumiTuple.root",dir.Data(),year,part);
		if(!gSystem->AccessPathName(str)) {
			lumi->Add(str);
		}
	} else {
//TODO		boost::progress_display show_addfile_progress( 3000 );
		for(int i=0; i<3000; ++i) {
//TODO			++show_addfile_progress;
			sprintf(str,"%s/%d/%d/output.root",dir.Data(),year,i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			sprintf(str,"%s/%d/%d/LumiTuple.root",dir.Data(),year,i);
			if(gSystem->AccessPathName(str)) continue;
			lumi->Add(str);
		}
	}
	//if(year==100) {
	//	tagType_ = JetTag;
	//	boost::progress_display show_addfile_progress( 1051 );
	//	for(int i=0; i<1051; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/577/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==101) {
	//	tagType_ = JetTag;
	//	backwards_=true;
	//	boost::progress_display show_addfile_progress( 1051 );
	//	for(int i=0; i<1051; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/583/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==102) {
	//	tagType_ = JpsiTag;
	//	boost::progress_display show_addfile_progress( 514 );
	//	for(int i=0; i<514; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/633/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==110) {//repeat with D0 PID selections turned off
	//	tagType_ = JetTag;
	//	boost::progress_display show_addfile_progress( 1051 );
	//	for(int i=0; i<1051; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/577/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==120) {//repeat with SV->{D0,D+} branches in
	//	tagType_ = JetTag;
	//	boost::progress_display show_addfile_progress( 1051 );
	//	for(int i=0; i<1051; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/577/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==122) {//repeat with SV->{D0,D+} branches in
	//	tagType_ = JpsiTag;
	//	boost::progress_display show_addfile_progress( 514 );
	//	for(int i=0; i<514; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/633/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==140) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=5000.;
	//	maxTruePt_=10000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/564/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/568/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==141) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=10000.;
	//	maxTruePt_=15000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/565/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/569/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==142) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=15000.;
	//	maxTruePt_=20000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/566/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/570/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==143) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=20000.;
	//	maxTruePt_=50000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/559/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/562/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==144) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=50000.;
	//	maxTruePt_=100000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/567/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/571/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==150) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=5;
	//	minTruePt_=5000.;
	//	maxTruePt_=10000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/596/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/600/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==151) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=5;
	//	minTruePt_=10000.;
	//	maxTruePt_=15000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/597/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/601/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==152) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=5;
	//	minTruePt_=15000.;
	//	maxTruePt_=20000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/598/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/602/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==153) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=5;
	//	minTruePt_=20000.;
	//	maxTruePt_=50000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/560/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/563/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==154) {
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=5;
	//	minTruePt_=50000.;
	//	maxTruePt_=100000.;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/599/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/603/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==164) {
	//	isMC_=true;
	//	tagType_ = TruthNoMuTag;
	//	flavour_=4;
	//	minTruePt_=0.;
	//	maxTruePt_=200000.;
	//	boost::progress_display show_addfile_progress( 488+452 );
	//	for(int i=0; i<488; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/5/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<452; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/7/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==165) {
	//	isMC_=true;
	//	tagType_ = TruthNoMuTag;
	//	flavour_=5;
	//	minTruePt_=0.;
	//	maxTruePt_=200000.;
	//	boost::progress_display show_addfile_progress( 417+431 );
	//	for(int i=0; i<417; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/6/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<431; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/8/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==200) {//mixed b->J/psi + incl. J/psi sample
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 61 );
	//	for(int i=0; i<30; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/608/%d/output.root",i);//214
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<31; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/609/%d/output.root",i);//214
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==201) {//D0 special
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 200 );
	//	for(int i=0; i<200; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/546/%d/output.root",i);//380
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==202) {//Dp special
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 200 );
	//	for(int i=0; i<200; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/547/%d/output.root",i);//381
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==203) {//Ds special
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 200 );
	//	for(int i=0; i<200; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/548/%d/output.root",i);//382
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==204) {//Lc special
	//	isMC_=true;
	//	tagType_ = NoJets;
	////TODO	boost::progress_display show_addfile_progress( 200 );
	////TODO	for(int i=0; i<200; ++i) {
	////TODO		++show_addfile_progress;
	////TODO		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/377/%d/output.root",i);//377//359
	////TODO		if(gSystem->AccessPathName(str)) continue;
	////TODO		t->Add(str);
	////TODO	}
	//} else if(year==205) {//K3pi special
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 220 );
	//	for(int i=0; i<220; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/549/%d/output.root",i);//383
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==240) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/564/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/568/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==241) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/565/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/569/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==242) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/566/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/570/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==243) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/559/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/562/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==244) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/567/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/571/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==250) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/596/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/600/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==251) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/597/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/601/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==252) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/598/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/602/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==253) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/560/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/563/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year==254) {
	//	isMC_=true;
	//	tagType_ = NoJets;
	//	boost::progress_display show_addfile_progress( 20 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/599/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/603/%d/output.root",i);
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//} else if(year>=900) {//test runs
	//	//tagType_ = JetTag;
	//	//boost::progress_display show_addfile_progress( 10 );
	//	//for(int i=0; i<10; ++i) {
	//	//	++show_addfile_progress;
	//	//	sprintf(str,"/eos/user/d/dcraik/jets/577/%d/output.root",i);
	//	//	if(gSystem->AccessPathName(str)) continue;
	//	//	t->Add(str);
	//	//}
	//	//tagType_ = NoJets;
	//	//t->Add("/eos/lhcb/user/d/dcraik/jets/546/17/output.root");
	//	isMC_=true;
	//	tagType_ = TruthTag;
	//	flavour_=4;
	//	minTruePt_=20000.;
	//	maxTruePt_=50000.;
	//	boost::progress_display show_addfile_progress( 10 );
	//	for(int i=0; i<10; ++i) {
	//		++show_addfile_progress;
	//		sprintf(str,"/eos/user/d/dcraik/jets/559/%d/output.root",i);//380
	//		//sprintf(str,"/eos/lhcb/user/d/dcraik/jets/546/%d/output.root",i);//380
	//		if(gSystem->AccessPathName(str)) continue;
	//		t->Add(str);
	//	}
	//}
	//558 RIIJ_JMC20180315MagDown49000003
	//559 RIIJ_JMC20180315MagDown49000043    --- 143
	//560 RIIJ_JMC20180315MagDown49000053    --- 153
	//561 RIIJ_JMC20180315MagUp49000003
	//562 RIIJ_JMC20180315MagUp49000043      --- 143
	//563 RIIJ_JMC20180315MagUp49000053      --- 153
	//564 RIIJ_JMC20180315MagDown49000040    --- 140
	//565 RIIJ_JMC20180315MagDown49000041    --- 141
	//566 RIIJ_JMC20180315MagDown49000042    --- 142
	//567 RIIJ_JMC20180315MagDown49000044    --- 144
	//568 RIIJ_JMC20180315MagUp49000040      --- 140
	//569 RIIJ_JMC20180315MagUp49000041      --- 141
	//570 RIIJ_JMC20180315MagUp49000042      --- 142
	//571 RIIJ_JMC20180315MagUp49000044      --- 144
	//577 RIIJ_dijet20180325MagDown2016      --- 100
	//578 RIIJ_dijet20180325MagDown2017A
	//579 RIIJ_dijet20180325MagDown2017B
	//583 RIIJ_dijet_back20180325MagDown2016 --- 101
	//584 RIIJ_JMC20180325MagDown49000050    --- 150
	//585 RIIJ_JMC20180325MagDown49000051    --- 151
	//586 RIIJ_JMC20180325MagDown49000052    --- 152
	//587 RIIJ_JMC20180325MagDown49000054    --- 154
	//588 RIIJ_JMC20180325MagUp49000050      --- 150
	//589 RIIJ_JMC20180325MagUp49000051      --- 151
	//590 RIIJ_JMC20180325MagUp49000052      --- 152
	//591 RIIJ_JMC20180325MagUp49000054      --- 154
	if(isMC) {
		detOffsetX = 0.818;//MC
		detOffsetY = 0.173;//MC
	}
	std::cout << t->GetEntries() << std::endl;
	Init(t);
}

skimTuples::~skimTuples()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t skimTuples::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t skimTuples::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void skimTuples::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	gen_idx_pvr = 0;
	gen_idx_jet = 0;
	gen_idx_prnt = 0;
	gen_pid = 0;
	gen_q = 0;
	gen_p = 0;
	gen_pt = 0;
	gen_px = 0;
	gen_py = 0;
	gen_pz = 0;
	gen_e = 0;
	gen_x = 0;
	gen_y = 0;
	gen_z = 0;
	gen_prnt_pid = 0;
	gen_res_pid = 0;
	gen_from_sig = 0;
	pvr_x = 0;
	pvr_y = 0;
	pvr_z = 0;
	pvr_dx = 0;
	pvr_dy = 0;
	pvr_dz = 0;
	pvr_chi2 = 0;
	pvr_ndof = 0;
	svr_idx_pvr = 0;
	svr_idx_jet = 0;
	svr_idx_trk0 = 0;
	svr_idx_trk1 = 0;
	svr_idx_trk2 = 0;
	svr_idx_trk3 = 0;
	svr_idx_trk4 = 0;
	svr_idx_trk5 = 0;
	svr_idx_trk6 = 0;
	svr_idx_trk7 = 0;
	svr_idx_trk8 = 0;
	svr_idx_trk9 = 0;
	svr_p = 0;
	svr_pt = 0;
	svr_px = 0;
	svr_py = 0;
	svr_pz = 0;
	svr_e = 0;
	svr_x = 0;
	svr_y = 0;
	svr_z = 0;
	svr_dx = 0;
	svr_dy = 0;
	svr_dz = 0;
	svr_m = 0;
	svr_m_cor = 0;
	svr_m_cor_err = 0;
	svr_m_cor_err_full = 0;
	svr_fd_min = 0;
	svr_fd_chi2 = 0;
	svr_chi2 = 0;
	svr_ip_chi2 = 0;
	svr_ip_chi2_sum = 0;
	svr_ip_chi2_min_trk = 0;
	svr_abs_q_sum = 0;
	svr_tau = 0;
	svr_ntrk = 0;
	svr_ntrk_jet = 0;
	svr_jet_dr = 0;
	svr_jet_pt = 0;
	svr_pass = 0;
	svr_bdt0 = 0;
	svr_bdt1 = 0;
	jet_idx_pvr = 0;
	jet_ntrk = 0;
	jet_nneu = 0;
	jet_p = 0;
	jet_pt = 0;
	jet_px = 0;
	jet_py = 0;
	jet_pz = 0;
	jet_e = 0;
	jet_l0_hadron_tos = 0;
	jet_l0_photon_tos = 0;
	jet_l0_electron_tos = 0;
	jet_l0_muon_tos = 0;
	jet_l0_dimuon_tos = 0;
	jet_hlt1_track_tos = 0;
	jet_hlt1_ditrack_tos = 0;
	trk_idx_gen = 0;
	trk_idx_pvr = 0;
	trk_idx_jet = 0;
	trk_p = 0;
	trk_pt = 0;
	trk_px = 0;
	trk_py = 0;
	trk_pz = 0;
	trk_e = 0;
	trk_pid = 0;
	trk_q = 0;
	trk_ip = 0;
	trk_ip_chi2 = 0;
	trk_pnn_e = 0;
	trk_pnn_mu = 0;
	trk_pnn_pi = 0;
	trk_pnn_k = 0;
	trk_pnn_p = 0;
	trk_pnn_ghost = 0;
	trk_ecal = 0;
	trk_hcal = 0;
	trk_prb_ghost = 0;
	trk_type = 0;
	trk_is_mu = 0;
	trk_vid = 0;
	trk_x = 0;
	trk_y = 0;
	trk_z = 0;
	trk_chi2 = 0;
	trk_ndof = 0;
	neu_idx_gen = 0;
	neu_idx_jet = 0;
	neu_p = 0;
	neu_pt = 0;
	neu_px = 0;
	neu_py = 0;
	neu_pz = 0;
	neu_e = 0;
	neu_pid = 0;
	z0_idx_pvr = 0;
	z0_idx_jet_trk0 = 0;
	z0_idx_jet_trk1 = 0;
	z0_idx_jet_dr = 0;
	z0_p = 0;
	z0_pt = 0;
	z0_px = 0;
	z0_py = 0;
	z0_pz = 0;
	z0_e = 0;
	z0_x = 0;
	z0_y = 0;
	z0_z = 0;
	z0_m = 0;
	z0_ip = 0;
	z0_ip_chi2 = 0;
	z0_dr_jet = 0;
	z0_idx_trk0 = 0;
	z0_idx_trk1 = 0;
	z0_l0_muon_tos = 0;
	z0_l0_dimuon_tos = 0;
	z0_hlt1_hiptmu_tos = 0;
	z0_hlt2_hiptmu_tos = 0;
	d0_idx_pvr = 0;
	d0_idx_jet = 0;
	d0_p = 0;
	d0_pt = 0;
	d0_px = 0;
	d0_py = 0;
	d0_pz = 0;
	d0_e = 0;
	d0_x = 0;
	d0_y = 0;
	d0_z = 0;
	d0_m = 0;
	d0_ip = 0;
	d0_ip_chi2 = 0;
	d0_vtx_chi2 = 0;
	d0_vtx_ndof = 0;
	d0_fd = 0;
	d0_fd_chi2 = 0;
	d0_tau = 0;
	d0_tau_err = 0;
	d0_tau_chi2 = 0;
	d0_ntrk_jet = 0;
	d0_idx_trk0 = 0;
	d0_idx_trk1 = 0;
	dp_idx_pvr = 0;
	dp_idx_jet = 0;
	dp_p = 0;
	dp_pt = 0;
	dp_px = 0;
	dp_py = 0;
	dp_pz = 0;
	dp_e = 0;
	dp_x = 0;
	dp_y = 0;
	dp_z = 0;
	dp_m = 0;
	dp_ip = 0;
	dp_ip_chi2 = 0;
	dp_vtx_chi2 = 0;
	dp_vtx_ndof = 0;
	dp_fd = 0;
	dp_fd_chi2 = 0;
	dp_tau = 0;
	dp_tau_err = 0;
	dp_tau_chi2 = 0;
	dp_ntrk_jet = 0;
	dp_idx_trk0 = 0;
	dp_idx_trk1 = 0;
	dp_idx_trk2 = 0;
	ds_idx_pvr = 0;
	ds_idx_jet = 0;
	ds_p = 0;
	ds_pt = 0;
	ds_px = 0;
	ds_py = 0;
	ds_pz = 0;
	ds_e = 0;
	ds_x = 0;
	ds_y = 0;
	ds_z = 0;
	ds_m = 0;
	ds_ip = 0;
	ds_ip_chi2 = 0;
	ds_vtx_chi2 = 0;
	ds_vtx_ndof = 0;
	ds_fd = 0;
	ds_fd_chi2 = 0;
	ds_tau = 0;
	ds_tau_err = 0;
	ds_tau_chi2 = 0;
	ds_ntrk_jet = 0;
	ds_idx_trk0 = 0;
	ds_idx_trk1 = 0;
	ds_idx_trk2 = 0;
	lc_idx_pvr = 0;
	lc_idx_jet = 0;
	lc_p = 0;
	lc_pt = 0;
	lc_px = 0;
	lc_py = 0;
	lc_pz = 0;
	lc_e = 0;
	lc_x = 0;
	lc_y = 0;
	lc_z = 0;
	lc_m = 0;
	lc_ip = 0;
	lc_ip_chi2 = 0;
	lc_vtx_chi2 = 0;
	lc_vtx_ndof = 0;
	lc_fd = 0;
	lc_fd_chi2 = 0;
	lc_tau = 0;
	lc_tau_err = 0;
	lc_tau_chi2 = 0;
	lc_ntrk_jet = 0;
	lc_idx_trk0 = 0;
	lc_idx_trk1 = 0;
	lc_idx_trk2 = 0;
	k3pi_idx_pvr = 0;
	k3pi_idx_jet = 0;
	k3pi_p = 0;
	k3pi_pt = 0;
	k3pi_px = 0;
	k3pi_py = 0;
	k3pi_pz = 0;
	k3pi_e = 0;
	k3pi_x = 0;
	k3pi_y = 0;
	k3pi_z = 0;
	k3pi_m = 0;
	k3pi_ip = 0;
	k3pi_ip_chi2 = 0;
	k3pi_vtx_chi2 = 0;
	k3pi_vtx_ndof = 0;
	k3pi_fd = 0;
	k3pi_fd_chi2 = 0;
	k3pi_tau = 0;
	k3pi_tau_err = 0;
	k3pi_tau_chi2 = 0;
	k3pi_ntrk_jet = 0;
	k3pi_idx_trk0 = 0;
	k3pi_idx_trk1 = 0;
	k3pi_idx_trk2 = 0;
	k3pi_idx_trk3 = 0;
	jpsi_idx_pvr = 0;
	jpsi_idx_jet = 0;
	jpsi_p = 0;
	jpsi_pt = 0;
	jpsi_px = 0;
	jpsi_py = 0;
	jpsi_pz = 0;
	jpsi_e = 0;
	jpsi_x = 0;
	jpsi_y = 0;
	jpsi_z = 0;
	jpsi_m = 0;
	jpsi_ip = 0;
	jpsi_ip_chi2 = 0;
	jpsi_vtx_chi2 = 0;
	jpsi_vtx_ndof = 0;
	jpsi_fd = 0;
	jpsi_fd_chi2 = 0;
	jpsi_tau = 0;
	jpsi_tau_err = 0;
	jpsi_tau_chi2 = 0;
	jpsi_ntrk_jet = 0;
	jpsi_idx_trk0 = 0;
	jpsi_idx_trk1 = 0;
	evt_dec = 0;
	evt_j1_idx = 0;
	evt_j1_dR = 0;
	evt_j1_nsv = 0;
	evt_j1_nmu = 0;
	evt_j1_ntrk = 0;
	evt_j1_nneu = 0;
	evt_j1_p = 0;
	evt_j1_pt = 0;
	evt_j1_px = 0;
	evt_j1_py = 0;
	evt_j1_pz = 0;
	evt_j1_e = 0;
	evt_j2_idx = 0;
	evt_j2_dR = 0;
	evt_j2_nsv = 0;
	evt_j2_nmu = 0;
	evt_j2_ntrk = 0;
	evt_j2_nneu = 0;
	evt_j2_p = 0;
	evt_j2_pt = 0;
	evt_j2_px = 0;
	evt_j2_py = 0;
	evt_j2_pz = 0;
	evt_j2_e = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("gen_idx_pvr", &gen_idx_pvr, &b_gen_idx_pvr);
	fChain->SetBranchAddress("gen_idx_jet", &gen_idx_jet, &b_gen_idx_jet);
	fChain->SetBranchAddress("gen_idx_prnt", &gen_idx_prnt, &b_gen_idx_prnt);
	fChain->SetBranchAddress("gen_pid", &gen_pid, &b_gen_pid);
	fChain->SetBranchAddress("gen_q", &gen_q, &b_gen_q);
	fChain->SetBranchAddress("gen_p", &gen_p, &b_gen_p);
	fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
	fChain->SetBranchAddress("gen_px", &gen_px, &b_gen_px);
	fChain->SetBranchAddress("gen_py", &gen_py, &b_gen_py);
	fChain->SetBranchAddress("gen_pz", &gen_pz, &b_gen_pz);
	fChain->SetBranchAddress("gen_e", &gen_e, &b_gen_e);
	fChain->SetBranchAddress("gen_x", &gen_x, &b_gen_x);
	fChain->SetBranchAddress("gen_y", &gen_y, &b_gen_y);
	fChain->SetBranchAddress("gen_z", &gen_z, &b_gen_z);
	fChain->SetBranchAddress("gen_prnt_pid", &gen_prnt_pid, &b_gen_prnt_pid);
	fChain->SetBranchAddress("gen_res_pid", &gen_res_pid, &b_gen_res_pid);
	fChain->SetBranchAddress("gen_from_sig", &gen_from_sig, &b_gen_from_sig);
	fChain->SetBranchAddress("pvr_x", &pvr_x, &b_pvr_x);
	fChain->SetBranchAddress("pvr_y", &pvr_y, &b_pvr_y);
	fChain->SetBranchAddress("pvr_z", &pvr_z, &b_pvr_z);
	fChain->SetBranchAddress("pvr_dx", &pvr_dx, &b_pvr_dx);
	fChain->SetBranchAddress("pvr_dy", &pvr_dy, &b_pvr_dy);
	fChain->SetBranchAddress("pvr_dz", &pvr_dz, &b_pvr_dz);
	fChain->SetBranchAddress("pvr_chi2", &pvr_chi2, &b_pvr_chi2);
	fChain->SetBranchAddress("pvr_ndof", &pvr_ndof, &b_pvr_ndof);
	fChain->SetBranchAddress("svr_idx_pvr", &svr_idx_pvr, &b_svr_idx_pvr);
	fChain->SetBranchAddress("svr_idx_jet", &svr_idx_jet, &b_svr_idx_jet);
	fChain->SetBranchAddress("svr_idx_trk0", &svr_idx_trk0, &b_svr_idx_trk0);
	fChain->SetBranchAddress("svr_idx_trk1", &svr_idx_trk1, &b_svr_idx_trk1);
	fChain->SetBranchAddress("svr_idx_trk2", &svr_idx_trk2, &b_svr_idx_trk2);
	fChain->SetBranchAddress("svr_idx_trk3", &svr_idx_trk3, &b_svr_idx_trk3);
	fChain->SetBranchAddress("svr_idx_trk4", &svr_idx_trk4, &b_svr_idx_trk4);
	fChain->SetBranchAddress("svr_idx_trk5", &svr_idx_trk5, &b_svr_idx_trk5);
	fChain->SetBranchAddress("svr_idx_trk6", &svr_idx_trk6, &b_svr_idx_trk6);
	fChain->SetBranchAddress("svr_idx_trk7", &svr_idx_trk7, &b_svr_idx_trk7);
	fChain->SetBranchAddress("svr_idx_trk8", &svr_idx_trk8, &b_svr_idx_trk8);
	fChain->SetBranchAddress("svr_idx_trk9", &svr_idx_trk9, &b_svr_idx_trk9);
	fChain->SetBranchAddress("svr_p", &svr_p, &b_svr_p);
	fChain->SetBranchAddress("svr_pt", &svr_pt, &b_svr_pt);
	fChain->SetBranchAddress("svr_px", &svr_px, &b_svr_px);
	fChain->SetBranchAddress("svr_py", &svr_py, &b_svr_py);
	fChain->SetBranchAddress("svr_pz", &svr_pz, &b_svr_pz);
	fChain->SetBranchAddress("svr_e", &svr_e, &b_svr_e);
	fChain->SetBranchAddress("svr_x", &svr_x, &b_svr_x);
	fChain->SetBranchAddress("svr_y", &svr_y, &b_svr_y);
	fChain->SetBranchAddress("svr_z", &svr_z, &b_svr_z);
	fChain->SetBranchAddress("svr_dx", &svr_dx, &b_svr_dx);
	fChain->SetBranchAddress("svr_dy", &svr_dy, &b_svr_dy);
	fChain->SetBranchAddress("svr_dz", &svr_dz, &b_svr_dz);
	fChain->SetBranchAddress("svr_m", &svr_m, &b_svr_m);
	fChain->SetBranchAddress("svr_m_cor", &svr_m_cor, &b_svr_m_cor);
	fChain->SetBranchAddress("svr_m_cor_err", &svr_m_cor_err, &b_svr_m_cor_err);
	fChain->SetBranchAddress("svr_m_cor_err_full", &svr_m_cor_err_full, &b_svr_m_cor_err_full);
	fChain->SetBranchAddress("svr_fd_min", &svr_fd_min, &b_svr_fd_min);
	fChain->SetBranchAddress("svr_fd_chi2", &svr_fd_chi2, &b_svr_fd_chi2);
	fChain->SetBranchAddress("svr_chi2", &svr_chi2, &b_svr_chi2);
	fChain->SetBranchAddress("svr_ip_chi2", &svr_ip_chi2, &b_svr_ip_chi2);
	fChain->SetBranchAddress("svr_ip_chi2_sum", &svr_ip_chi2_sum, &b_svr_ip_chi2_sum);
	fChain->SetBranchAddress("svr_ip_chi2_min_trk", &svr_ip_chi2_min_trk, &b_svr_ip_chi2_min_trk);
	fChain->SetBranchAddress("svr_abs_q_sum", &svr_abs_q_sum, &b_svr_abs_q_sum);
	fChain->SetBranchAddress("svr_tau", &svr_tau, &b_svr_tau);
	fChain->SetBranchAddress("svr_ntrk", &svr_ntrk, &b_svr_ntrk);
	fChain->SetBranchAddress("svr_ntrk_jet", &svr_ntrk_jet, &b_svr_ntrk_jet);
	fChain->SetBranchAddress("svr_jet_dr", &svr_jet_dr, &b_svr_jet_dr);
	fChain->SetBranchAddress("svr_jet_pt", &svr_jet_pt, &b_svr_jet_pt);
	fChain->SetBranchAddress("svr_pass", &svr_pass, &b_svr_pass);
	fChain->SetBranchAddress("svr_bdt0", &svr_bdt0, &b_svr_bdt0);
	fChain->SetBranchAddress("svr_bdt1", &svr_bdt1, &b_svr_bdt1);
	fChain->SetBranchAddress("jet_idx_pvr", &jet_idx_pvr, &b_jet_idx_pvr);
	fChain->SetBranchAddress("jet_ntrk", &jet_ntrk, &b_jet_ntrk);
	fChain->SetBranchAddress("jet_nneu", &jet_nneu, &b_jet_nneu);
	fChain->SetBranchAddress("jet_p", &jet_p, &b_jet_p);
	fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
	fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
	fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
	fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
	fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
	if(!isMC_) {
		fChain->SetBranchAddress("jet_l0_hadron_tos",    &jet_l0_hadron_tos,    &b_jet_l0_hadron_tos);
		fChain->SetBranchAddress("jet_l0_photon_tos",    &jet_l0_photon_tos,    &b_jet_l0_photon_tos);
		fChain->SetBranchAddress("jet_l0_electron_tos",  &jet_l0_electron_tos,  &b_jet_l0_electron_tos);
		fChain->SetBranchAddress("jet_l0_muon_tos",      &jet_l0_muon_tos,      &b_jet_l0_muon_tos);
		fChain->SetBranchAddress("jet_l0_dimuon_tos",    &jet_l0_dimuon_tos,    &b_jet_l0_dimuon_tos);
		fChain->SetBranchAddress("jet_hlt1_track_tos",   &jet_hlt1_track_tos,   &b_jet_hlt1_track_tos);
		fChain->SetBranchAddress("jet_hlt1_ditrack_tos", &jet_hlt1_ditrack_tos, &b_jet_hlt1_ditrack_tos);
	}
	fChain->SetBranchAddress("trk_idx_gen", &trk_idx_gen, &b_trk_idx_gen);
	fChain->SetBranchAddress("trk_idx_pvr", &trk_idx_pvr, &b_trk_idx_pvr);
	fChain->SetBranchAddress("trk_idx_jet", &trk_idx_jet, &b_trk_idx_jet);
	fChain->SetBranchAddress("trk_p", &trk_p, &b_trk_p);
	fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
	fChain->SetBranchAddress("trk_px", &trk_px, &b_trk_px);
	fChain->SetBranchAddress("trk_py", &trk_py, &b_trk_py);
	fChain->SetBranchAddress("trk_pz", &trk_pz, &b_trk_pz);
	fChain->SetBranchAddress("trk_e", &trk_e, &b_trk_e);
	fChain->SetBranchAddress("trk_pid", &trk_pid, &b_trk_pid);
	fChain->SetBranchAddress("trk_q", &trk_q, &b_trk_q);
	fChain->SetBranchAddress("trk_ip", &trk_ip, &b_trk_ip);
	fChain->SetBranchAddress("trk_ip_chi2", &trk_ip_chi2, &b_trk_ip_chi2);
	fChain->SetBranchAddress("trk_pnn_e", &trk_pnn_e, &b_trk_pnn_e);
	fChain->SetBranchAddress("trk_pnn_mu", &trk_pnn_mu, &b_trk_pnn_mu);
	fChain->SetBranchAddress("trk_pnn_pi", &trk_pnn_pi, &b_trk_pnn_pi);
	fChain->SetBranchAddress("trk_pnn_k", &trk_pnn_k, &b_trk_pnn_k);
	fChain->SetBranchAddress("trk_pnn_p", &trk_pnn_p, &b_trk_pnn_p);
	fChain->SetBranchAddress("trk_pnn_ghost", &trk_pnn_ghost, &b_trk_pnn_ghost);
	fChain->SetBranchAddress("trk_ecal", &trk_ecal, &b_trk_ecal);
	fChain->SetBranchAddress("trk_hcal", &trk_hcal, &b_trk_hcal);
	fChain->SetBranchAddress("trk_prb_ghost", &trk_prb_ghost, &b_trk_prb_ghost);
	fChain->SetBranchAddress("trk_type", &trk_type, &b_trk_type);
	fChain->SetBranchAddress("trk_is_mu", &trk_is_mu, &b_trk_is_mu);
	fChain->SetBranchAddress("trk_vid", &trk_vid, &b_trk_vid);
	fChain->SetBranchAddress("trk_x", &trk_x, &b_trk_x);
	fChain->SetBranchAddress("trk_y", &trk_y, &b_trk_y);
	fChain->SetBranchAddress("trk_z", &trk_z, &b_trk_z);
	fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
	fChain->SetBranchAddress("trk_ndof", &trk_ndof, &b_trk_ndof);
	fChain->SetBranchAddress("neu_idx_gen", &neu_idx_gen, &b_neu_idx_gen);
	fChain->SetBranchAddress("neu_idx_jet", &neu_idx_jet, &b_neu_idx_jet);
	fChain->SetBranchAddress("neu_p", &neu_p, &b_neu_p);
	fChain->SetBranchAddress("neu_pt", &neu_pt, &b_neu_pt);
	fChain->SetBranchAddress("neu_px", &neu_px, &b_neu_px);
	fChain->SetBranchAddress("neu_py", &neu_py, &b_neu_py);
	fChain->SetBranchAddress("neu_pz", &neu_pz, &b_neu_pz);
	fChain->SetBranchAddress("neu_e", &neu_e, &b_neu_e);
	fChain->SetBranchAddress("neu_pid", &neu_pid, &b_neu_pid);
	fChain->SetBranchAddress("z0_idx_pvr", &z0_idx_pvr, &b_z0_idx_pvr);
	fChain->SetBranchAddress("z0_idx_jet_trk0", &z0_idx_jet_trk0, &b_z0_idx_jet_trk0);
	fChain->SetBranchAddress("z0_idx_jet_trk1", &z0_idx_jet_trk1, &b_z0_idx_jet_trk1);
	fChain->SetBranchAddress("z0_idx_jet_dr", &z0_idx_jet_dr, &b_z0_idx_jet_dr);
	fChain->SetBranchAddress("z0_p", &z0_p, &b_z0_p);
	fChain->SetBranchAddress("z0_pt", &z0_pt, &b_z0_pt);
	fChain->SetBranchAddress("z0_px", &z0_px, &b_z0_px);
	fChain->SetBranchAddress("z0_py", &z0_py, &b_z0_py);
	fChain->SetBranchAddress("z0_pz", &z0_pz, &b_z0_pz);
	fChain->SetBranchAddress("z0_e", &z0_e, &b_z0_e);
	fChain->SetBranchAddress("z0_x", &z0_x, &b_z0_x);
	fChain->SetBranchAddress("z0_y", &z0_y, &b_z0_y);
	fChain->SetBranchAddress("z0_z", &z0_z, &b_z0_z);
	fChain->SetBranchAddress("z0_m", &z0_m, &b_z0_m);
	fChain->SetBranchAddress("z0_ip", &z0_ip, &b_z0_ip);
	fChain->SetBranchAddress("z0_ip_chi2", &z0_ip_chi2, &b_z0_ip_chi2);
	fChain->SetBranchAddress("z0_dr_jet", &z0_dr_jet, &b_z0_dr_jet);
	fChain->SetBranchAddress("z0_idx_trk0", &z0_idx_trk0, &b_z0_idx_trk0);
	fChain->SetBranchAddress("z0_idx_trk1", &z0_idx_trk1, &b_z0_idx_trk1);
	fChain->SetBranchAddress("z0_l0_muon_tos", &z0_l0_muon_tos, &b_z0_l0_muon_tos);
	fChain->SetBranchAddress("z0_l0_dimuon_tos", &z0_l0_dimuon_tos, &b_z0_l0_dimuon_tos);
	fChain->SetBranchAddress("z0_hlt1_hiptmu_tos", &z0_hlt1_hiptmu_tos, &b_z0_hlt1_hiptmu_tos);
	fChain->SetBranchAddress("z0_hlt2_hiptmu_tos", &z0_hlt2_hiptmu_tos, &b_z0_hlt2_hiptmu_tos);
	fChain->SetBranchAddress("d0_idx_pvr", &d0_idx_pvr, &b_d0_idx_pvr);
	fChain->SetBranchAddress("d0_idx_jet", &d0_idx_jet, &b_d0_idx_jet);
	fChain->SetBranchAddress("d0_p", &d0_p, &b_d0_p);
	fChain->SetBranchAddress("d0_pt", &d0_pt, &b_d0_pt);
	fChain->SetBranchAddress("d0_px", &d0_px, &b_d0_px);
	fChain->SetBranchAddress("d0_py", &d0_py, &b_d0_py);
	fChain->SetBranchAddress("d0_pz", &d0_pz, &b_d0_pz);
	fChain->SetBranchAddress("d0_e", &d0_e, &b_d0_e);
	fChain->SetBranchAddress("d0_x", &d0_x, &b_d0_x);
	fChain->SetBranchAddress("d0_y", &d0_y, &b_d0_y);
	fChain->SetBranchAddress("d0_z", &d0_z, &b_d0_z);
	fChain->SetBranchAddress("d0_m", &d0_m, &b_d0_m);
	fChain->SetBranchAddress("d0_ip", &d0_ip, &b_d0_ip);
	fChain->SetBranchAddress("d0_ip_chi2", &d0_ip_chi2, &b_d0_ip_chi2);
	fChain->SetBranchAddress("d0_vtx_chi2", &d0_vtx_chi2, &b_d0_vtx_chi2);
	fChain->SetBranchAddress("d0_vtx_ndof", &d0_vtx_ndof, &b_d0_vtx_ndof);
	fChain->SetBranchAddress("d0_fd", &d0_fd, &b_d0_fd);
	fChain->SetBranchAddress("d0_fd_chi2", &d0_fd_chi2, &b_d0_fd_chi2);
	fChain->SetBranchAddress("d0_tau", &d0_tau, &b_d0_tau);
	fChain->SetBranchAddress("d0_tau_err", &d0_tau_err, &b_d0_tau_err);
	fChain->SetBranchAddress("d0_tau_chi2", &d0_tau_chi2, &b_d0_tau_chi2);
	fChain->SetBranchAddress("d0_ntrk_jet", &d0_ntrk_jet, &b_d0_ntrk_jet);
	fChain->SetBranchAddress("d0_idx_trk0", &d0_idx_trk0, &b_d0_idx_trk0);
	fChain->SetBranchAddress("d0_idx_trk1", &d0_idx_trk1, &b_d0_idx_trk1);
	fChain->SetBranchAddress("dp_idx_pvr", &dp_idx_pvr, &b_dp_idx_pvr);
	fChain->SetBranchAddress("dp_idx_jet", &dp_idx_jet, &b_dp_idx_jet);
	fChain->SetBranchAddress("dp_p", &dp_p, &b_dp_p);
	fChain->SetBranchAddress("dp_pt", &dp_pt, &b_dp_pt);
	fChain->SetBranchAddress("dp_px", &dp_px, &b_dp_px);
	fChain->SetBranchAddress("dp_py", &dp_py, &b_dp_py);
	fChain->SetBranchAddress("dp_pz", &dp_pz, &b_dp_pz);
	fChain->SetBranchAddress("dp_e", &dp_e, &b_dp_e);
	fChain->SetBranchAddress("dp_x", &dp_x, &b_dp_x);
	fChain->SetBranchAddress("dp_y", &dp_y, &b_dp_y);
	fChain->SetBranchAddress("dp_z", &dp_z, &b_dp_z);
	fChain->SetBranchAddress("dp_m", &dp_m, &b_dp_m);
	fChain->SetBranchAddress("dp_ip", &dp_ip, &b_dp_ip);
	fChain->SetBranchAddress("dp_ip_chi2", &dp_ip_chi2, &b_dp_ip_chi2);
	fChain->SetBranchAddress("dp_vtx_chi2", &dp_vtx_chi2, &b_dp_vtx_chi2);
	fChain->SetBranchAddress("dp_vtx_ndof", &dp_vtx_ndof, &b_dp_vtx_ndof);
	fChain->SetBranchAddress("dp_fd", &dp_fd, &b_dp_fd);
	fChain->SetBranchAddress("dp_fd_chi2", &dp_fd_chi2, &b_dp_fd_chi2);
	fChain->SetBranchAddress("dp_tau", &dp_tau, &b_dp_tau);
	fChain->SetBranchAddress("dp_tau_err", &dp_tau_err, &b_dp_tau_err);
	fChain->SetBranchAddress("dp_tau_chi2", &dp_tau_chi2, &b_dp_tau_chi2);
	fChain->SetBranchAddress("dp_ntrk_jet", &dp_ntrk_jet, &b_dp_ntrk_jet);
	fChain->SetBranchAddress("dp_idx_trk0", &dp_idx_trk0, &b_dp_idx_trk0);
	fChain->SetBranchAddress("dp_idx_trk1", &dp_idx_trk1, &b_dp_idx_trk1);
	fChain->SetBranchAddress("dp_idx_trk2", &dp_idx_trk2, &b_dp_idx_trk2);
	fChain->SetBranchAddress("ds_idx_pvr", &ds_idx_pvr, &b_ds_idx_pvr);
	fChain->SetBranchAddress("ds_idx_jet", &ds_idx_jet, &b_ds_idx_jet);
	fChain->SetBranchAddress("ds_p", &ds_p, &b_ds_p);
	fChain->SetBranchAddress("ds_pt", &ds_pt, &b_ds_pt);
	fChain->SetBranchAddress("ds_px", &ds_px, &b_ds_px);
	fChain->SetBranchAddress("ds_py", &ds_py, &b_ds_py);
	fChain->SetBranchAddress("ds_pz", &ds_pz, &b_ds_pz);
	fChain->SetBranchAddress("ds_e", &ds_e, &b_ds_e);
	fChain->SetBranchAddress("ds_x", &ds_x, &b_ds_x);
	fChain->SetBranchAddress("ds_y", &ds_y, &b_ds_y);
	fChain->SetBranchAddress("ds_z", &ds_z, &b_ds_z);
	fChain->SetBranchAddress("ds_m", &ds_m, &b_ds_m);
	fChain->SetBranchAddress("ds_ip", &ds_ip, &b_ds_ip);
	fChain->SetBranchAddress("ds_ip_chi2", &ds_ip_chi2, &b_ds_ip_chi2);
	fChain->SetBranchAddress("ds_vtx_chi2", &ds_vtx_chi2, &b_ds_vtx_chi2);
	fChain->SetBranchAddress("ds_vtx_ndof", &ds_vtx_ndof, &b_ds_vtx_ndof);
	fChain->SetBranchAddress("ds_fd", &ds_fd, &b_ds_fd);
	fChain->SetBranchAddress("ds_fd_chi2", &ds_fd_chi2, &b_ds_fd_chi2);
	fChain->SetBranchAddress("ds_tau", &ds_tau, &b_ds_tau);
	fChain->SetBranchAddress("ds_tau_err", &ds_tau_err, &b_ds_tau_err);
	fChain->SetBranchAddress("ds_tau_chi2", &ds_tau_chi2, &b_ds_tau_chi2);
	fChain->SetBranchAddress("ds_ntrk_jet", &ds_ntrk_jet, &b_ds_ntrk_jet);
	fChain->SetBranchAddress("ds_idx_trk0", &ds_idx_trk0, &b_ds_idx_trk0);
	fChain->SetBranchAddress("ds_idx_trk1", &ds_idx_trk1, &b_ds_idx_trk1);
	fChain->SetBranchAddress("ds_idx_trk2", &ds_idx_trk2, &b_ds_idx_trk2);
	fChain->SetBranchAddress("lc_idx_pvr", &lc_idx_pvr, &b_lc_idx_pvr);
	fChain->SetBranchAddress("lc_idx_jet", &lc_idx_jet, &b_lc_idx_jet);
	fChain->SetBranchAddress("lc_p", &lc_p, &b_lc_p);
	fChain->SetBranchAddress("lc_pt", &lc_pt, &b_lc_pt);
	fChain->SetBranchAddress("lc_px", &lc_px, &b_lc_px);
	fChain->SetBranchAddress("lc_py", &lc_py, &b_lc_py);
	fChain->SetBranchAddress("lc_pz", &lc_pz, &b_lc_pz);
	fChain->SetBranchAddress("lc_e", &lc_e, &b_lc_e);
	fChain->SetBranchAddress("lc_x", &lc_x, &b_lc_x);
	fChain->SetBranchAddress("lc_y", &lc_y, &b_lc_y);
	fChain->SetBranchAddress("lc_z", &lc_z, &b_lc_z);
	fChain->SetBranchAddress("lc_m", &lc_m, &b_lc_m);
	fChain->SetBranchAddress("lc_ip", &lc_ip, &b_lc_ip);
	fChain->SetBranchAddress("lc_ip_chi2", &lc_ip_chi2, &b_lc_ip_chi2);
	fChain->SetBranchAddress("lc_vtx_chi2", &lc_vtx_chi2, &b_lc_vtx_chi2);
	fChain->SetBranchAddress("lc_vtx_ndof", &lc_vtx_ndof, &b_lc_vtx_ndof);
	fChain->SetBranchAddress("lc_fd", &lc_fd, &b_lc_fd);
	fChain->SetBranchAddress("lc_fd_chi2", &lc_fd_chi2, &b_lc_fd_chi2);
	fChain->SetBranchAddress("lc_tau", &lc_tau, &b_lc_tau);
	fChain->SetBranchAddress("lc_tau_err", &lc_tau_err, &b_lc_tau_err);
	fChain->SetBranchAddress("lc_tau_chi2", &lc_tau_chi2, &b_lc_tau_chi2);
	fChain->SetBranchAddress("lc_ntrk_jet", &lc_ntrk_jet, &b_lc_ntrk_jet);
	fChain->SetBranchAddress("lc_idx_trk0", &lc_idx_trk0, &b_lc_idx_trk0);
	fChain->SetBranchAddress("lc_idx_trk1", &lc_idx_trk1, &b_lc_idx_trk1);
	fChain->SetBranchAddress("lc_idx_trk2", &lc_idx_trk2, &b_lc_idx_trk2);
	fChain->SetBranchAddress("k3pi_idx_pvr", &k3pi_idx_pvr, &b_k3pi_idx_pvr);
	fChain->SetBranchAddress("k3pi_idx_jet", &k3pi_idx_jet, &b_k3pi_idx_jet);
	fChain->SetBranchAddress("k3pi_p", &k3pi_p, &b_k3pi_p);
	fChain->SetBranchAddress("k3pi_pt", &k3pi_pt, &b_k3pi_pt);
	fChain->SetBranchAddress("k3pi_px", &k3pi_px, &b_k3pi_px);
	fChain->SetBranchAddress("k3pi_py", &k3pi_py, &b_k3pi_py);
	fChain->SetBranchAddress("k3pi_pz", &k3pi_pz, &b_k3pi_pz);
	fChain->SetBranchAddress("k3pi_e", &k3pi_e, &b_k3pi_e);
	fChain->SetBranchAddress("k3pi_x", &k3pi_x, &b_k3pi_x);
	fChain->SetBranchAddress("k3pi_y", &k3pi_y, &b_k3pi_y);
	fChain->SetBranchAddress("k3pi_z", &k3pi_z, &b_k3pi_z);
	fChain->SetBranchAddress("k3pi_m", &k3pi_m, &b_k3pi_m);
	fChain->SetBranchAddress("k3pi_ip", &k3pi_ip, &b_k3pi_ip);
	fChain->SetBranchAddress("k3pi_ip_chi2", &k3pi_ip_chi2, &b_k3pi_ip_chi2);
	fChain->SetBranchAddress("k3pi_vtx_chi2", &k3pi_vtx_chi2, &b_k3pi_vtx_chi2);
	fChain->SetBranchAddress("k3pi_vtx_ndof", &k3pi_vtx_ndof, &b_k3pi_vtx_ndof);
	fChain->SetBranchAddress("k3pi_fd", &k3pi_fd, &b_k3pi_fd);
	fChain->SetBranchAddress("k3pi_fd_chi2", &k3pi_fd_chi2, &b_k3pi_fd_chi2);
	fChain->SetBranchAddress("k3pi_tau", &k3pi_tau, &b_k3pi_tau);
	fChain->SetBranchAddress("k3pi_tau_err", &k3pi_tau_err, &b_k3pi_tau_err);
	fChain->SetBranchAddress("k3pi_tau_chi2", &k3pi_tau_chi2, &b_k3pi_tau_chi2);
	fChain->SetBranchAddress("k3pi_ntrk_jet", &k3pi_ntrk_jet, &b_k3pi_ntrk_jet);
	fChain->SetBranchAddress("k3pi_idx_trk0", &k3pi_idx_trk0, &b_k3pi_idx_trk0);
	fChain->SetBranchAddress("k3pi_idx_trk1", &k3pi_idx_trk1, &b_k3pi_idx_trk1);
	fChain->SetBranchAddress("k3pi_idx_trk2", &k3pi_idx_trk2, &b_k3pi_idx_trk2);
	fChain->SetBranchAddress("k3pi_idx_trk3", &k3pi_idx_trk3, &b_k3pi_idx_trk3);
	fChain->SetBranchAddress("jpsi_idx_pvr", &jpsi_idx_pvr, &b_jpsi_idx_pvr);
	fChain->SetBranchAddress("jpsi_idx_jet", &jpsi_idx_jet, &b_jpsi_idx_jet);
	fChain->SetBranchAddress("jpsi_p", &jpsi_p, &b_jpsi_p);
	fChain->SetBranchAddress("jpsi_pt", &jpsi_pt, &b_jpsi_pt);
	fChain->SetBranchAddress("jpsi_px", &jpsi_px, &b_jpsi_px);
	fChain->SetBranchAddress("jpsi_py", &jpsi_py, &b_jpsi_py);
	fChain->SetBranchAddress("jpsi_pz", &jpsi_pz, &b_jpsi_pz);
	fChain->SetBranchAddress("jpsi_e", &jpsi_e, &b_jpsi_e);
	fChain->SetBranchAddress("jpsi_x", &jpsi_x, &b_jpsi_x);
	fChain->SetBranchAddress("jpsi_y", &jpsi_y, &b_jpsi_y);
	fChain->SetBranchAddress("jpsi_z", &jpsi_z, &b_jpsi_z);
	fChain->SetBranchAddress("jpsi_m", &jpsi_m, &b_jpsi_m);
	fChain->SetBranchAddress("jpsi_ip", &jpsi_ip, &b_jpsi_ip);
	fChain->SetBranchAddress("jpsi_ip_chi2", &jpsi_ip_chi2, &b_jpsi_ip_chi2);
	fChain->SetBranchAddress("jpsi_vtx_chi2", &jpsi_vtx_chi2, &b_jpsi_vtx_chi2);
	fChain->SetBranchAddress("jpsi_vtx_ndof", &jpsi_vtx_ndof, &b_jpsi_vtx_ndof);
	fChain->SetBranchAddress("jpsi_fd", &jpsi_fd, &b_jpsi_fd);
	fChain->SetBranchAddress("jpsi_fd_chi2", &jpsi_fd_chi2, &b_jpsi_fd_chi2);
	fChain->SetBranchAddress("jpsi_tau", &jpsi_tau, &b_jpsi_tau);
	fChain->SetBranchAddress("jpsi_tau_err", &jpsi_tau_err, &b_jpsi_tau_err);
	fChain->SetBranchAddress("jpsi_tau_chi2", &jpsi_tau_chi2, &b_jpsi_tau_chi2);
	fChain->SetBranchAddress("jpsi_ntrk_jet", &jpsi_ntrk_jet, &b_jpsi_ntrk_jet);
	fChain->SetBranchAddress("jpsi_idx_trk0", &jpsi_idx_trk0, &b_jpsi_idx_trk0);
	fChain->SetBranchAddress("jpsi_idx_trk1", &jpsi_idx_trk1, &b_jpsi_idx_trk1);
	fChain->SetBranchAddress("evt_dec", &evt_dec, &b_evt_dec);
	fChain->SetBranchAddress("evt_j1_idx", &evt_j1_idx, &b_evt_j1_idx);
	fChain->SetBranchAddress("evt_j1_dR", &evt_j1_dR, &b_evt_j1_dR);
	fChain->SetBranchAddress("evt_j1_nsv", &evt_j1_nsv, &b_evt_j1_nsv);
	fChain->SetBranchAddress("evt_j1_nmu", &evt_j1_nmu, &b_evt_j1_nmu);
	fChain->SetBranchAddress("evt_j1_ntrk", &evt_j1_ntrk, &b_evt_j1_ntrk);
	fChain->SetBranchAddress("evt_j1_nneu", &evt_j1_nneu, &b_evt_j1_nneu);
	fChain->SetBranchAddress("evt_j1_p", &evt_j1_p, &b_evt_j1_p);
	fChain->SetBranchAddress("evt_j1_pt", &evt_j1_pt, &b_evt_j1_pt);
	fChain->SetBranchAddress("evt_j1_px", &evt_j1_px, &b_evt_j1_px);
	fChain->SetBranchAddress("evt_j1_py", &evt_j1_py, &b_evt_j1_py);
	fChain->SetBranchAddress("evt_j1_pz", &evt_j1_pz, &b_evt_j1_pz);
	fChain->SetBranchAddress("evt_j1_e", &evt_j1_e, &b_evt_j1_e);
	fChain->SetBranchAddress("evt_j2_idx", &evt_j2_idx, &b_evt_j2_idx);
	fChain->SetBranchAddress("evt_j2_dR", &evt_j2_dR, &b_evt_j2_dR);
	fChain->SetBranchAddress("evt_j2_nsv", &evt_j2_nsv, &b_evt_j2_nsv);
	fChain->SetBranchAddress("evt_j2_nmu", &evt_j2_nmu, &b_evt_j2_nmu);
	fChain->SetBranchAddress("evt_j2_ntrk", &evt_j2_ntrk, &b_evt_j2_ntrk);
	fChain->SetBranchAddress("evt_j2_nneu", &evt_j2_nneu, &b_evt_j2_nneu);
	fChain->SetBranchAddress("evt_j2_p", &evt_j2_p, &b_evt_j2_p);
	fChain->SetBranchAddress("evt_j2_pt", &evt_j2_pt, &b_evt_j2_pt);
	fChain->SetBranchAddress("evt_j2_px", &evt_j2_px, &b_evt_j2_px);
	fChain->SetBranchAddress("evt_j2_py", &evt_j2_py, &b_evt_j2_py);
	fChain->SetBranchAddress("evt_j2_pz", &evt_j2_pz, &b_evt_j2_pz);
	fChain->SetBranchAddress("evt_j2_e", &evt_j2_e, &b_evt_j2_e);
	fChain->SetBranchAddress("evt_pvr_n", &evt_pvr_n, &b_evt_pvr_n);
	fChain->SetBranchAddress("evt_trk_n", &evt_trk_n, &b_evt_trk_n);

	for(int i=0; i<10; i++){
		svtrk[i] = new std::vector<double>();
		TString str="svr_idx_trk"; str+=i;
		fChain->SetBranchAddress(str,&svtrk[i]);
	}

	fout = TFile::Open(outName,"recreate");
	tout = new TTree("T","");
	if(lumi->GetNtrees()>0) lumiout = lumi->CloneTree();
	else lumiout = 0;

	tout->Branch("Evt",&EVT);
	tout->Branch("Dec",&DEC);
	tout->Branch("NPV",&NPV);
	tout->Branch("JetPx",&JPX);
	tout->Branch("JetPy",&JPY);
	tout->Branch("JetPz",&JPZ);
	tout->Branch("JetE",&JE);
	tout->Branch("JetPT",&JPT);
	tout->Branch("JetEta",&JETA);
	if(isMC_) {
		tout->Branch("JetTruePx",&JTRUEPX);
		tout->Branch("JetTruePy",&JTRUEPY);
		tout->Branch("JetTruePz",&JTRUEPZ);
		tout->Branch("JetTrueE",&JTRUEE);
		tout->Branch("JetTruePT",&JTRUEPT);
		tout->Branch("JetTrueEta",&JTRUEETA);
		tout->Branch("JetTrueDR",&JTRUEDR);
	}
	tout->Branch("JetSigma1",&JS1);
	tout->Branch("JetSigma2",&JS2);
	tout->Branch("JetQ",&JQ);
	tout->Branch("JetMult",&JN);
	tout->Branch("JetNChr",&JNQ);
	tout->Branch("JetNNeu",&JNN);
	tout->Branch("JetPTD",&JPTD);
	tout->Branch("JetTrigPT",&JTRIGPT);
	tout->Branch("JetTrig10",&JTRIG10);
	tout->Branch("JetTrig17",&JTRIG17);
	tout->Branch("JetTrig60",&JTRIG60);
	tout->Branch("JetTRUEb",&JTRUEb);
	tout->Branch("JetTRUEc",&JTRUEc);
	tout->Branch("JetTRUED0",&JTRUED0);
	tout->Branch("JetTRUEDP",&JTRUEDP);
	tout->Branch("JetTRUEDS",&JTRUEDS);
	tout->Branch("JetTRUEJPSI",&JTRUEJPSI);
	tout->Branch("JetTRUEDSV",&JTRUEDSV);
	tout->Branch("JetTRUEBSV",&JTRUEBSV);
	tout->Branch("JetTRUESV",&JTRUESV);
	tout->Branch("JetTRUEDPX",&JTRUEDPX);
	tout->Branch("JetTRUEDPY",&JTRUEDPY);
	tout->Branch("JetTRUEDPZ",&JTRUEDPZ);
	tout->Branch("JetTRUEDE", &JTRUEDE);
	tout->Branch("JetTRUEDPT",&JTRUEDPT);
	tout->Branch("JetTRUEDDR",&JTRUEDDR);
	tout->Branch("TagPx",&TPX);
	tout->Branch("TagPy",&TPY);
	tout->Branch("TagPz",&TPZ);
	tout->Branch("TagE",&TE);
	tout->Branch("TagPT",&TPT);
	if(isMC_) {
		tout->Branch("TagTruePx",&TTRUEPX);
		tout->Branch("TagTruePy",&TTRUEPY);
		tout->Branch("TagTruePz",&TTRUEPZ);
		tout->Branch("TagTrueE",&TTRUEE);
		tout->Branch("TagTruePT",&TTRUEPT);
		tout->Branch("TagTrueEta",&TTRUEETA);
		tout->Branch("TagTrueDR",&TTRUEDR);
	}
	tout->Branch("TagTrigPT",&TTRIGPT);
	tout->Branch("TagTrig10",&TTRIG10);
	tout->Branch("TagTrig17",&TTRIG17);
	tout->Branch("TagTrig60",&TTRIG60);
	tout->Branch("TagEta",&TETA);
	tout->Branch("DeltaR",&DR);
	tout->Branch("DeltaPhi",&DPHI);

	tout->Branch("ZM"  ,&ZM  );
	tout->Branch("ZP"  ,&ZP  );
	tout->Branch("ZPX" ,&ZPX );
	tout->Branch("ZPY" ,&ZPY );
	tout->Branch("ZPZ" ,&ZPZ );
	tout->Branch("ZPT" ,&ZPT );
	tout->Branch("ZE"  ,&ZE  );
	tout->Branch("ZETA",&ZETA);
	tout->Branch("ZY",&ZY);
	tout->Branch("ZDR" ,&ZDR );
	tout->Branch("MU0PX" ,&MU0PX );
	tout->Branch("MU0PY" ,&MU0PY );
	tout->Branch("MU0PZ" ,&MU0PZ );
	tout->Branch("MU0PT" ,&MU0PT );
	tout->Branch("MU0IP" ,&MU0IP );
	tout->Branch("MU0IPCHI2" ,&MU0IPCHI2 );
	tout->Branch("MU0FPT" ,&MU0FPT );
	tout->Branch("MU0DR" ,&MU0DR );
	tout->Branch("MU1PX" ,&MU1PX );
	tout->Branch("MU1PY" ,&MU1PY );
	tout->Branch("MU1PZ" ,&MU1PZ );
	tout->Branch("MU1PT" ,&MU1PT );
	tout->Branch("MU1IP" ,&MU1IP );
	tout->Branch("MU1IPCHI2" ,&MU1IPCHI2 );
	tout->Branch("MU1FPT" ,&MU1FPT );
	tout->Branch("MU1DR" ,&MU1DR );

	tout->Branch("ZTRUEM"  ,&ZTRUEM  );
	tout->Branch("ZTRUEP"  ,&ZTRUEP  );
	tout->Branch("ZTRUEPX" ,&ZTRUEPX );
	tout->Branch("ZTRUEPY" ,&ZTRUEPY );
	tout->Branch("ZTRUEPZ" ,&ZTRUEPZ );
	tout->Branch("ZTRUEPT" ,&ZTRUEPT );
	tout->Branch("ZTRUEE"  ,&ZTRUEE  );
	tout->Branch("ZTRUEETA",&ZTRUEETA);
	tout->Branch("ZTRUEDR" ,&ZTRUEDR );
	tout->Branch("MU0TRUEPX" ,&MU0TRUEPX );
	tout->Branch("MU0TRUEPY" ,&MU0TRUEPY );
	tout->Branch("MU0TRUEPZ" ,&MU0TRUEPZ );
	tout->Branch("MU0TRUEPT" ,&MU0TRUEPT );
	tout->Branch("MU1TRUEPX" ,&MU1TRUEPX );
	tout->Branch("MU1TRUEPY" ,&MU1TRUEPY );
	tout->Branch("MU1TRUEPZ" ,&MU1TRUEPZ );
	tout->Branch("MU1TRUEPT" ,&MU1TRUEPT );

	 SVM           = new std::vector<double>();
	 SVMCOR        = new std::vector<double>();
	 SVMCORERR     = new std::vector<double>();
	 SVMINPERP     = new std::vector<double>();
	 SVPT          = new std::vector<double>();
	 SVDRJ         = new std::vector<double>();
	 SVDRT         = new std::vector<double>();
	 SVN           = new std::vector<double>();
	 SVNTRUE       = new std::vector<double>();
	 SVNLINKED     = new std::vector<double>();
	 SVNJ          = new std::vector<double>();
	 SVNT          = new std::vector<double>();
	 SVQ           = new std::vector<double>();
	 SVPERP        = new std::vector<double>();
	 SVETA         = new std::vector<double>();
	 SVTZ          = new std::vector<double>();
	 SVVXCHI2      = new std::vector<double>();
	 SVIPCHI2      = new std::vector<double>();
	 SVMINIPCHI2   = new std::vector<double>();
	 SVPX          = new std::vector<double>();
	 SVPY          = new std::vector<double>();
	 SVPZ          = new std::vector<double>();
	 SVE           = new std::vector<double>();
	 SVMAXGHOST    = new std::vector<double>();
	 SVX           = new std::vector<double>();
	 SVY           = new std::vector<double>();
	 SVZ           = new std::vector<double>();
	 SVSUMIPCHI2   = new std::vector<double>();
	 SVISD0        = new std::vector<double>();
	 SVISDP        = new std::vector<double>();
	 SVD0M         = new std::vector<double>();
	 SVDPM         = new std::vector<double>();
	 SVTRK0IDX     = new std::vector<double>();
	 SVTRK1IDX     = new std::vector<double>();
	 SVTRK2IDX     = new std::vector<double>();
	 SVTRK3IDX     = new std::vector<double>();
	 SVTRUEIDX     = new std::vector<double>();
	 SVTRUETRK0IDX = new std::vector<double>();
	 SVTRUETRK1IDX = new std::vector<double>();
	 SVTRUETRK2IDX = new std::vector<double>();
	 SVTRUETRK3IDX = new std::vector<double>();
	 SVTRK0P       = new std::vector<double>();
	 SVTRK1P       = new std::vector<double>();
	 SVTRK2P       = new std::vector<double>();
	 SVTRK3P       = new std::vector<double>();
	 SVTRK0PT      = new std::vector<double>();
	 SVTRK1PT      = new std::vector<double>();
	 SVTRK2PT      = new std::vector<double>();
	 SVTRK3PT      = new std::vector<double>();
	 SVTRK0PX      = new std::vector<double>();
	 SVTRK1PX      = new std::vector<double>();
	 SVTRK2PX      = new std::vector<double>();
	 SVTRK3PX      = new std::vector<double>();
	 SVTRK0PY      = new std::vector<double>();
	 SVTRK1PY      = new std::vector<double>();
	 SVTRK2PY      = new std::vector<double>();
	 SVTRK3PY      = new std::vector<double>();
	 SVTRK0PZ      = new std::vector<double>();
	 SVTRK1PZ      = new std::vector<double>();
	 SVTRK2PZ      = new std::vector<double>();
	 SVTRK3PZ      = new std::vector<double>();
	 SVTRK0PNNPI   = new std::vector<double>();
	 SVTRK1PNNPI   = new std::vector<double>();
	 SVTRK2PNNPI   = new std::vector<double>();
	 SVTRK3PNNPI   = new std::vector<double>();
	 SVTRK0PNNK    = new std::vector<double>();
	 SVTRK1PNNK    = new std::vector<double>();
	 SVTRK2PNNK    = new std::vector<double>();
	 SVTRK3PNNK    = new std::vector<double>();

	 TSVM           = new std::vector<double>();
	 TSVMCOR        = new std::vector<double>();
	 TSVMCORERR     = new std::vector<double>();
	 TSVMINPERP     = new std::vector<double>();
	 TSVPT          = new std::vector<double>();
	 TSVDRJ         = new std::vector<double>();
	 TSVDRT         = new std::vector<double>();
	 TSVN           = new std::vector<double>();
	 TSVNTRUE       = new std::vector<double>();
	 TSVNJ          = new std::vector<double>();
	 TSVNT          = new std::vector<double>();
	 TSVQ           = new std::vector<double>();
	 TSVPERP        = new std::vector<double>();
	 TSVETA         = new std::vector<double>();
	 TSVTZ          = new std::vector<double>();
	 TSVIPCHI2      = new std::vector<double>();
	 TSVMINIPCHI2   = new std::vector<double>();
	 TSVPX          = new std::vector<double>();
	 TSVPY          = new std::vector<double>();
	 TSVPZ          = new std::vector<double>();
	 TSVE           = new std::vector<double>();
	 TSVMAXGHOST    = new std::vector<double>();
	 TSVX           = new std::vector<double>();
	 TSVY           = new std::vector<double>();
	 TSVZ           = new std::vector<double>();
	 TSVSUMIPCHI2   = new std::vector<double>();
	 TSVISD0        = new std::vector<double>();
	 TSVISDP        = new std::vector<double>();
	 TSVD0M         = new std::vector<double>();
	 TSVDPM         = new std::vector<double>();
	 TSVTRUEIDX     = new std::vector<double>();

	 TRKPT          = new std::vector<double>();
	 TRKIPCHI2      = new std::vector<double>();
	 TRKINJET       = new std::vector<double>();

	tout->Branch("SVX",&SVX);
	tout->Branch("SVY",&SVY);
	tout->Branch("SVZ",&SVZ);
	tout->Branch("SVPerp",&SVPERP);
	tout->Branch("SVPx",&SVPX);
	tout->Branch("SVPy",&SVPY);
	tout->Branch("SVPz",&SVPZ);
	tout->Branch("SVE",&SVE);
	tout->Branch("SVPT",&SVPT);
	tout->Branch("SVETA",&SVETA);
	tout->Branch("SVM",&SVM);
	tout->Branch("SVMCor",&SVMCOR);
	tout->Branch("SVMCorErr",&SVMCORERR);
	tout->Branch("SVMINPERP",&SVMINPERP);
	tout->Branch("SVDRJ",&SVDRJ);
	tout->Branch("SVDRT",&SVDRT);
	tout->Branch("SVN",&SVN);
	tout->Branch("SVNTRUE",&SVNTRUE);
	tout->Branch("SVNLINKED",&SVNLINKED);
	tout->Branch("SVNJ",&SVNJ);
	tout->Branch("SVNT",&SVNT);
	tout->Branch("SVQ",&SVQ);
	tout->Branch("SVSumIPChi2",&SVSUMIPCHI2);
	tout->Branch("SVTZ",&SVTZ);
	tout->Branch("SVVXCHI2",&SVVXCHI2);
	tout->Branch("SVIPCHI2",&SVIPCHI2);
	tout->Branch("SVMINIPCHI2",&SVMINIPCHI2);
	tout->Branch("SVGhostMax",&SVMAXGHOST);
	tout->Branch("SVISD0",&SVISD0);
	tout->Branch("SVISDP",&SVISDP);
	tout->Branch("SVD0M",&SVD0M);
	tout->Branch("SVDPM",&SVDPM);
	tout->Branch("SVTRK0IDX",&SVTRK0IDX);
	tout->Branch("SVTRK1IDX",&SVTRK1IDX);
	tout->Branch("SVTRK2IDX",&SVTRK2IDX);
	tout->Branch("SVTRK3IDX",&SVTRK3IDX);
	tout->Branch("SVTRUEIDX",&SVTRUEIDX);
	tout->Branch("SVTRUETRK0IDX",&SVTRUETRK0IDX);
	tout->Branch("SVTRUETRK1IDX",&SVTRUETRK1IDX);
	tout->Branch("SVTRUETRK2IDX",&SVTRUETRK2IDX);
	tout->Branch("SVTRUETRK3IDX",&SVTRUETRK3IDX);
	tout->Branch("SVTRK0P", &SVTRK0P);
	tout->Branch("SVTRK1P", &SVTRK1P);
	tout->Branch("SVTRK2P", &SVTRK2P);
	tout->Branch("SVTRK3P", &SVTRK3P);
	tout->Branch("SVTRK0PT",&SVTRK0PT);
	tout->Branch("SVTRK1PT",&SVTRK1PT);
	tout->Branch("SVTRK2PT",&SVTRK2PT);
	tout->Branch("SVTRK3PT",&SVTRK3PT);
	tout->Branch("SVTRK0PX",&SVTRK0PX);
	tout->Branch("SVTRK1PX",&SVTRK1PX);
	tout->Branch("SVTRK2PX",&SVTRK2PX);
	tout->Branch("SVTRK3PX",&SVTRK3PX);
	tout->Branch("SVTRK0PY",&SVTRK0PY);
	tout->Branch("SVTRK1PY",&SVTRK1PY);
	tout->Branch("SVTRK2PY",&SVTRK2PY);
	tout->Branch("SVTRK3PY",&SVTRK3PY);
	tout->Branch("SVTRK0PZ",&SVTRK0PZ);
	tout->Branch("SVTRK1PZ",&SVTRK1PZ);
	tout->Branch("SVTRK2PZ",&SVTRK2PZ);
	tout->Branch("SVTRK3PZ",&SVTRK3PZ);
	tout->Branch("SVTRK0PNNPI",&SVTRK0PNNPI);
	tout->Branch("SVTRK1PNNPI",&SVTRK1PNNPI);
	tout->Branch("SVTRK2PNNPI",&SVTRK2PNNPI);
	tout->Branch("SVTRK3PNNPI",&SVTRK3PNNPI);
	tout->Branch("SVTRK0PNNK",&SVTRK0PNNK);
	tout->Branch("SVTRK1PNNK",&SVTRK1PNNK);
	tout->Branch("SVTRK2PNNK",&SVTRK2PNNK);
	tout->Branch("SVTRK3PNNK",&SVTRK3PNNK);

	tout->Branch("TSVX",&TSVX);
	tout->Branch("TSVY",&TSVY);
	tout->Branch("TSVZ",&TSVZ);
	tout->Branch("TSVPerp",&TSVPERP);
	tout->Branch("TSVPx",&TSVPX);
	tout->Branch("TSVPy",&TSVPY);
	tout->Branch("TSVPz",&TSVPZ);
	tout->Branch("TSVE",&TSVE);
	tout->Branch("TSVPT",&TSVPT);
	tout->Branch("TSVETA",&TSVETA);
	tout->Branch("TSVM",&TSVM);
	tout->Branch("TSVMCor",&TSVMCOR);
	tout->Branch("TSVMCorErr",&TSVMCORERR);
	tout->Branch("TSVMINPERP",&TSVMINPERP);
	tout->Branch("TSVDRJ",&TSVDRJ);
	tout->Branch("TSVDRT",&TSVDRT);
	tout->Branch("TSVN",&TSVN);
	tout->Branch("TSVNTRUE",&TSVNTRUE);
	tout->Branch("TSVNJ",&TSVNJ);
	tout->Branch("TSVNT",&TSVNT);
	tout->Branch("TSVQ",&TSVQ);
	tout->Branch("TSVSumIPChi2",&TSVSUMIPCHI2);
	tout->Branch("TSVTZ",&TSVTZ);
	tout->Branch("TSVIPCHI2",&TSVIPCHI2);
	tout->Branch("TSVMINIPCHI2",&TSVMINIPCHI2);
	tout->Branch("TSVGhostMax",&TSVMAXGHOST);
	tout->Branch("TSVISD0",&TSVISD0);
	tout->Branch("TSVISDP",&TSVISDP);
	tout->Branch("TSVD0M",&TSVD0M);
	tout->Branch("TSVDPM",&TSVDPM);
	tout->Branch("TSVTRUEIDX",&TSVTRUEIDX);

	tout->Branch("TRKPT",&TRKPT);
	tout->Branch("TRKIPCHI2",&TRKIPCHI2);
	tout->Branch("TRKINJET",&TRKINJET);

	tout->Branch("NSV",&NSV);
	tout->Branch("NSVTRK",&NSVTRK);
	tout->Branch("NSVTRUETRK",&NSVTRUETRK);
	tout->Branch("NTSV",&NTSV);
	tout->Branch("NTRK",&NTRK);
	tout->Branch("NNEU",&NNEU);
	tout->Branch("PVX",&PVX);
	tout->Branch("PVY",&PVY);
	tout->Branch("PVZ",&PVZ);
	tout->Branch("NDispl6",&NDISPL6);
	tout->Branch("NDispl9",&NDISPL9);
	tout->Branch("NDispl16",&NDISPL16);
	tout->Branch("MuPT",&MUPT);
	tout->Branch("MuIPChi2",&MUIPCHI2);
	tout->Branch("MuDR",&MUDR);
	tout->Branch("MuPNN",&MUPNN);
	tout->Branch("NMu",&NMU);
	tout->Branch("HardPT",&HPT);
	tout->Branch("HardIPChi2",&HIPCHI2);
	tout->Branch("HardDR",&HDR);

	TRUEBID  = new std::vector<double>();
	TRUEBPX   = new std::vector<double>();
	TRUEBPY   = new std::vector<double>();
	TRUEBPZ   = new std::vector<double>();
	TRUEBPT   = new std::vector<double>();
	TRUEBE    = new std::vector<double>();
	TRUEBX   = new std::vector<double>();
	TRUEBY   = new std::vector<double>();
	TRUEBZ   = new std::vector<double>();
	TRUEDID  = new std::vector<double>();
	TRUEDPX   = new std::vector<double>();
	TRUEDPY   = new std::vector<double>();
	TRUEDPZ   = new std::vector<double>();
	TRUEDPT   = new std::vector<double>();
	TRUEDE    = new std::vector<double>();
	TRUEDX   = new std::vector<double>();
	TRUEDY   = new std::vector<double>();
	TRUEDZ   = new std::vector<double>();
	TRUEDFROMB= new std::vector<double>();
	TRUEDTRUEB= new std::vector<double>();
	TRUEDSEL  = new std::vector<double>();
	TRUEDTRK0IDX = new std::vector<double>();
	TRUEDTRK0ID  = new std::vector<double>();
	TRUEDTRK0P   = new std::vector<double>();
	TRUEDTRK0PT  = new std::vector<double>();
	TRUEDTRK0INACC = new std::vector<double>();
	TRUEDTRK0RECO  = new std::vector<double>();
	TRUEDTRK0PNNK  = new std::vector<double>();
	TRUEDTRK0PNNPI = new std::vector<double>();
	TRUEDTRK1IDX = new std::vector<double>();
	TRUEDTRK1ID  = new std::vector<double>();
	TRUEDTRK1P   = new std::vector<double>();
	TRUEDTRK1PT  = new std::vector<double>();
	TRUEDTRK1INACC = new std::vector<double>();
	TRUEDTRK1RECO  = new std::vector<double>();
	TRUEDTRK1PNNK  = new std::vector<double>();
	TRUEDTRK1PNNPI = new std::vector<double>();
	TRUEDTRK2IDX = new std::vector<double>();
	TRUEDTRK2ID  = new std::vector<double>();
	TRUEDTRK2P   = new std::vector<double>();
	TRUEDTRK2PT  = new std::vector<double>();
	TRUEDTRK2INACC = new std::vector<double>();
	TRUEDTRK2RECO  = new std::vector<double>();
	TRUEDTRK2PNNK  = new std::vector<double>();
	TRUEDTRK2PNNPI = new std::vector<double>();
	TRUEDTRK3IDX = new std::vector<double>();
	TRUEDTRK3ID  = new std::vector<double>();
	TRUEDTRK3P   = new std::vector<double>();
	TRUEDTRK3PT  = new std::vector<double>();
	TRUEDTRK3INACC = new std::vector<double>();
	TRUEDTRK3RECO  = new std::vector<double>();
	TRUEDTRK3PNNK  = new std::vector<double>();
	TRUEDTRK3PNNPI = new std::vector<double>();
	TRUEDTRK0PX   = new std::vector<double>();
	TRUEDTRK1PX   = new std::vector<double>();
	TRUEDTRK2PX   = new std::vector<double>();
	TRUEDTRK3PX   = new std::vector<double>();
	TRUEDTRK0PY   = new std::vector<double>();
	TRUEDTRK1PY   = new std::vector<double>();
	TRUEDTRK2PY   = new std::vector<double>();
	TRUEDTRK3PY   = new std::vector<double>();
	TRUEDTRK0PZ   = new std::vector<double>();
	TRUEDTRK1PZ   = new std::vector<double>();
	TRUEDTRK2PZ   = new std::vector<double>();
	TRUEDTRK3PZ   = new std::vector<double>();

	D0M       = new std::vector<double>();
	D0PX      = new std::vector<double>();
	D0PY      = new std::vector<double>();
	D0PZ      = new std::vector<double>();
	D0E       = new std::vector<double>();
	D0X       = new std::vector<double>();
	D0Y       = new std::vector<double>();
	D0Z       = new std::vector<double>();
	D0IP      = new std::vector<double>();
	D0IPCHI2  = new std::vector<double>();
	D0FD      = new std::vector<double>();
	D0FDCHI2  = new std::vector<double>();
	D0TAU     = new std::vector<double>();
	D0DIRA    = new std::vector<double>();
	D0VTXCHI2 = new std::vector<double>();
	D0VTXNDOF = new std::vector<double>();
	D0DRJET   = new std::vector<double>();
	D0DRTAG   = new std::vector<double>();
	D0PIPNNK  = new std::vector<double>();
	D0KPNNK   = new std::vector<double>();
	D0PIPNNPI = new std::vector<double>();
	D0KPNNPI  = new std::vector<double>();
	D0PIWEIGHT= new std::vector<double>();
	D0KWEIGHT = new std::vector<double>();

	D0P       = new std::vector<double>();
	D0PT      = new std::vector<double>();
	D0ETA     = new std::vector<double>();
	D0PIP     = new std::vector<double>();
	D0KP      = new std::vector<double>();
	D0PIPT    = new std::vector<double>();
	D0KPT     = new std::vector<double>();
	D0PIETA   = new std::vector<double>();
	D0KETA    = new std::vector<double>();
	D0PIPX    = new std::vector<double>();
	D0KPX     = new std::vector<double>();
	D0PIPY    = new std::vector<double>();
	D0KPY     = new std::vector<double>();
	D0PIPZ    = new std::vector<double>();
	D0KPZ     = new std::vector<double>();
	D0PIIPCHI2= new std::vector<double>();
	D0KIPCHI2 = new std::vector<double>();

	D0TRK0    = new std::vector<double>();
	D0TRK1    = new std::vector<double>();
	D0TRUETRK0= new std::vector<double>();
	D0TRUETRK1= new std::vector<double>();

	D0NJ      = new std::vector<double>();
	D0MAXDR   = new std::vector<double>();

	D0TRUE    = new std::vector<double>();
	D0TRUEDR  = new std::vector<double>();
	D0TRUEIDX = new std::vector<double>();
	D0FROMB   = new std::vector<double>();

	tout->Branch("TRUEBID",   &TRUEBID);
	tout->Branch("TRUEBPX",   &TRUEBPX);
	tout->Branch("TRUEBPY",   &TRUEBPY);
	tout->Branch("TRUEBPZ",   &TRUEBPZ);
	tout->Branch("TRUEBPT",   &TRUEBPT);
	tout->Branch("TRUEBE",    &TRUEBE);
	tout->Branch("TRUEBX",   &TRUEBX);
	tout->Branch("TRUEBY",   &TRUEBY);
	tout->Branch("TRUEBZ",   &TRUEBZ);
	tout->Branch("TRUEDID",   &TRUEDID);
	tout->Branch("TRUEDPX",   &TRUEDPX);
	tout->Branch("TRUEDPY",   &TRUEDPY);
	tout->Branch("TRUEDPZ",   &TRUEDPZ);
	tout->Branch("TRUEDPT",   &TRUEDPT);
	tout->Branch("TRUEDE",    &TRUEDE);
	tout->Branch("TRUEDX",   &TRUEDX);
	tout->Branch("TRUEDY",   &TRUEDY);
	tout->Branch("TRUEDZ",   &TRUEDZ);
	tout->Branch("TRUEDFROMB",&TRUEDFROMB);
	tout->Branch("TRUEDTRUEB",&TRUEDTRUEB);
	tout->Branch("TRUEDSEL",  &TRUEDSEL);
	tout->Branch("TRUEDTRK0IDX", &TRUEDTRK0IDX);
	tout->Branch("TRUEDTRK0ID",  &TRUEDTRK0ID);
	tout->Branch("TRUEDTRK0P",   &TRUEDTRK0P);
	tout->Branch("TRUEDTRK0PT",  &TRUEDTRK0PT);
	tout->Branch("TRUEDTRK0INACC", &TRUEDTRK0INACC);
	tout->Branch("TRUEDTRK0RECO",  &TRUEDTRK0RECO);
	tout->Branch("TRUEDTRK0PNNK",  &TRUEDTRK0PNNK);
	tout->Branch("TRUEDTRK0PNNPI", &TRUEDTRK0PNNPI);
	tout->Branch("TRUEDTRK1IDX", &TRUEDTRK1IDX);
	tout->Branch("TRUEDTRK1ID",  &TRUEDTRK1ID);
	tout->Branch("TRUEDTRK1P",   &TRUEDTRK1P);
	tout->Branch("TRUEDTRK1PT",  &TRUEDTRK1PT);
	tout->Branch("TRUEDTRK1INACC", &TRUEDTRK1INACC);
	tout->Branch("TRUEDTRK1RECO",  &TRUEDTRK1RECO);
	tout->Branch("TRUEDTRK1PNNK",  &TRUEDTRK1PNNK);
	tout->Branch("TRUEDTRK1PNNPI", &TRUEDTRK1PNNPI);
	tout->Branch("TRUEDTRK2IDX", &TRUEDTRK2IDX);
	tout->Branch("TRUEDTRK2ID",  &TRUEDTRK2ID);
	tout->Branch("TRUEDTRK2P",   &TRUEDTRK2P);
	tout->Branch("TRUEDTRK2PT",  &TRUEDTRK2PT);
	tout->Branch("TRUEDTRK2INACC", &TRUEDTRK2INACC);
	tout->Branch("TRUEDTRK2RECO",  &TRUEDTRK2RECO);
	tout->Branch("TRUEDTRK2PNNK",  &TRUEDTRK2PNNK);
	tout->Branch("TRUEDTRK2PNNPI", &TRUEDTRK2PNNPI);
	tout->Branch("TRUEDTRK3IDX", &TRUEDTRK3IDX);
	tout->Branch("TRUEDTRK3ID",  &TRUEDTRK3ID);
	tout->Branch("TRUEDTRK3P",   &TRUEDTRK3P);
	tout->Branch("TRUEDTRK3PT",  &TRUEDTRK3PT);
	tout->Branch("TRUEDTRK3INACC", &TRUEDTRK3INACC);
	tout->Branch("TRUEDTRK3RECO",  &TRUEDTRK3RECO);
	tout->Branch("TRUEDTRK3PNNK",  &TRUEDTRK3PNNK);
	tout->Branch("TRUEDTRK3PNNPI", &TRUEDTRK3PNNPI);
	tout->Branch("TRUEDTRK0PX",   &TRUEDTRK0PX);
	tout->Branch("TRUEDTRK1PX",   &TRUEDTRK1PX);
	tout->Branch("TRUEDTRK2PX",   &TRUEDTRK2PX);
	tout->Branch("TRUEDTRK3PX",   &TRUEDTRK3PX);
	tout->Branch("TRUEDTRK0PY",   &TRUEDTRK0PY);
	tout->Branch("TRUEDTRK1PY",   &TRUEDTRK1PY);
	tout->Branch("TRUEDTRK2PY",   &TRUEDTRK2PY);
	tout->Branch("TRUEDTRK3PY",   &TRUEDTRK3PY);
	tout->Branch("TRUEDTRK0PZ",   &TRUEDTRK0PZ);
	tout->Branch("TRUEDTRK1PZ",   &TRUEDTRK1PZ);
	tout->Branch("TRUEDTRK2PZ",   &TRUEDTRK2PZ);
	tout->Branch("TRUEDTRK3PZ",   &TRUEDTRK3PZ);

	tout->Branch("D0M",        &D0M);
	tout->Branch("D0PX",       &D0PX);
	tout->Branch("D0PY",       &D0PY);
	tout->Branch("D0PZ",       &D0PZ);
	tout->Branch("D0E",        &D0E);
	tout->Branch("D0X",        &D0X);
	tout->Branch("D0Y",        &D0Y);
	tout->Branch("D0Z",        &D0Z);
	tout->Branch("D0IP",       &D0IP);
	tout->Branch("D0IPCHI2",   &D0IPCHI2);
	tout->Branch("D0FD",       &D0FD);
	tout->Branch("D0FDCHI2",   &D0FDCHI2);
	tout->Branch("D0TAU",      &D0TAU);
	tout->Branch("D0DIRA",      &D0DIRA);
	tout->Branch("D0VTXCHI2",  &D0VTXCHI2);
	tout->Branch("D0VTXNDOF",  &D0VTXNDOF);
	tout->Branch("D0DRJET",    &D0DRJET);
	tout->Branch("D0DRTAG",    &D0DRTAG);
	tout->Branch("D0PIPNNK",   &D0PIPNNK);
	tout->Branch("D0KPNNK",    &D0KPNNK);
	tout->Branch("D0PIPNNPI",  &D0PIPNNPI);
	tout->Branch("D0KPNNPI",   &D0KPNNPI);
	tout->Branch("D0PIWEIGHT", &D0PIWEIGHT);
	tout->Branch("D0KWEIGHT",  &D0KWEIGHT);

	tout->Branch("D0P",        &D0P);
	tout->Branch("D0PT",       &D0PT);
	tout->Branch("D0ETA",      &D0ETA);
	tout->Branch("D0PIP",      &D0PIP);
	tout->Branch("D0KP",       &D0KP);
	tout->Branch("D0PIPT",     &D0PIPT);
	tout->Branch("D0KPT",      &D0KPT);
	tout->Branch("D0PIETA",    &D0PIETA);
	tout->Branch("D0KETA",     &D0KETA);
	tout->Branch("D0PIPX",     &D0PIPX);
	tout->Branch("D0KPX",      &D0KPX);
	tout->Branch("D0PIPY",     &D0PIPY);
	tout->Branch("D0KPY",      &D0KPY);
	tout->Branch("D0PIPZ",     &D0PIPZ);
	tout->Branch("D0KPZ",      &D0KPZ);
	tout->Branch("D0PIIPCHI2", &D0PIIPCHI2);
	tout->Branch("D0KIPCHI2",  &D0KIPCHI2);

	tout->Branch("D0TRK0",     &D0TRK0);
	tout->Branch("D0TRK1",     &D0TRK1);
	tout->Branch("D0TRUETRK0", &D0TRUETRK0);
	tout->Branch("D0TRUETRK1", &D0TRUETRK1);

	tout->Branch("D0NJ",       &D0NJ);
	tout->Branch("D0MAXDR",    &D0MAXDR);

	tout->Branch("D0TRUE",    &D0TRUE);
	tout->Branch("D0TRUEDR",  &D0TRUEDR);
	tout->Branch("D0TRUEIDX", &D0TRUEIDX);
	tout->Branch("D0FROMB",   &D0FROMB);

	DPM       = new std::vector<double>();
	DPPX      = new std::vector<double>();
	DPPY      = new std::vector<double>();
	DPPZ      = new std::vector<double>();
	DPE       = new std::vector<double>();
	DPX       = new std::vector<double>();
	DPY       = new std::vector<double>();
	DPZ       = new std::vector<double>();
	DPIP      = new std::vector<double>();
	DPIPCHI2  = new std::vector<double>();
	DPFD      = new std::vector<double>();
	DPFDCHI2  = new std::vector<double>();
	DPTAU     = new std::vector<double>();
	DPDIRA    = new std::vector<double>();
	DPVTXCHI2 = new std::vector<double>();
	DPVTXNDOF = new std::vector<double>();
	DPDRJET   = new std::vector<double>();
	DPDRTAG   = new std::vector<double>();
	DPPI1PNN  = new std::vector<double>();
	DPPI2PNN  = new std::vector<double>();
	DPKPNN    = new std::vector<double>();

	DPP       = new std::vector<double>();
	DPPT      = new std::vector<double>();
	DPPI1PT   = new std::vector<double>();
	DPPI2PT   = new std::vector<double>();
	DPKPT     = new std::vector<double>();

	DPTRK0    = new std::vector<double>();
	DPTRK1    = new std::vector<double>();
	DPTRK2    = new std::vector<double>();

	DPNJ      = new std::vector<double>();
	DPMAXDR   = new std::vector<double>();

	DPTRUE    = new std::vector<double>();
	DPTRUEDR  = new std::vector<double>();
	DPTRUEIDX = new std::vector<double>();
	DPFROMB   = new std::vector<double>();

	tout->Branch("DM",        &DPM);
	tout->Branch("DPX",       &DPPX);
	tout->Branch("DPY",       &DPPY);
	tout->Branch("DPZ",       &DPPZ);
	tout->Branch("DE",        &DPE);
	tout->Branch("DX",        &DPX);
	tout->Branch("DY",        &DPY);
	tout->Branch("DZ",        &DPZ);
	tout->Branch("DIP",       &DPIP);
	tout->Branch("DIPCHI2",   &DPIPCHI2);
	tout->Branch("DFD",       &DPFD);
	tout->Branch("DFDCHI2",   &DPFDCHI2);
	tout->Branch("DTAU",      &DPTAU);
	tout->Branch("DDIRA",      &DPDIRA);
	tout->Branch("DVTXCHI2",  &DPVTXCHI2);
	tout->Branch("DVTXNDOF",  &DPVTXNDOF);
	tout->Branch("DDRJET",    &DPDRJET);
	tout->Branch("DDRTAG",    &DPDRTAG);
	tout->Branch("DPI1PNN",   &DPPI1PNN);
	tout->Branch("DPI2PNN",   &DPPI2PNN);
	tout->Branch("DKPNN",     &DPKPNN);

	tout->Branch("DP",        &DPP);
	tout->Branch("DPT",       &DPPT);
	tout->Branch("DPI1PT",    &DPPI1PT);
	tout->Branch("DPI2PT",    &DPPI2PT);
	tout->Branch("DKPT",      &DPKPT);

	tout->Branch("DTRK0",     &DPTRK0);
	tout->Branch("DTRK1",     &DPTRK1);
	tout->Branch("DTRK2",     &DPTRK2);

	tout->Branch("DNJ",       &DPNJ);
	tout->Branch("DMAXDR",    &DPMAXDR);

	tout->Branch("DTRUE",    &DPTRUE);
	tout->Branch("DTRUEDR",  &DPTRUEDR);
	tout->Branch("DTRUEIDX", &DPTRUEIDX);
	tout->Branch("DFROMB",   &DPFROMB );

	DSM       = new std::vector<double>();
	DSPX      = new std::vector<double>();
	DSPY      = new std::vector<double>();
	DSPZ      = new std::vector<double>();
	DSE       = new std::vector<double>();
	DSX       = new std::vector<double>();
	DSY       = new std::vector<double>();
	DSZ       = new std::vector<double>();
	DSIP      = new std::vector<double>();
	DSIPCHI2  = new std::vector<double>();
	DSFD      = new std::vector<double>();
	DSFDCHI2  = new std::vector<double>();
	DSTAU     = new std::vector<double>();
	DSDIRA    = new std::vector<double>();
	DSVTXCHI2 = new std::vector<double>();
	DSVTXNDOF = new std::vector<double>();
	DSDRJET   = new std::vector<double>();
	DSDRTAG   = new std::vector<double>();
	DSPIPNN   = new std::vector<double>();
	DSK1PNN   = new std::vector<double>();
	DSK2PNN   = new std::vector<double>();
	DSPHIM    = new std::vector<double>();

	DSP       = new std::vector<double>();
	DSPT      = new std::vector<double>();
	DSPIPT    = new std::vector<double>();
	DSK1PT    = new std::vector<double>();
	DSK2PT    = new std::vector<double>();
	DSPHIPT   = new std::vector<double>();

	DSTRK0    = new std::vector<double>();
	DSTRK1    = new std::vector<double>();
	DSTRK2    = new std::vector<double>();

	DSNJ      = new std::vector<double>();
	DSMAXDR   = new std::vector<double>();

	DSTRUE    = new std::vector<double>();
	DSTRUEDR  = new std::vector<double>();
	DSTRUEIDX = new std::vector<double>();
	DSFROMB   = new std::vector<double>();

	tout->Branch("DSM",        &DSM);
	tout->Branch("DSPX",       &DSPX);
	tout->Branch("DSPY",       &DSPY);
	tout->Branch("DSPZ",       &DSPZ);
	tout->Branch("DSE",        &DSE);
	tout->Branch("DSX",        &DSX);
	tout->Branch("DSY",        &DSY);
	tout->Branch("DSZ",        &DSZ);
	tout->Branch("DSIP",       &DSIP);
	tout->Branch("DSIPCHI2",   &DSIPCHI2);
	tout->Branch("DSFD",       &DSFD);
	tout->Branch("DSFDCHI2",   &DSFDCHI2);
	tout->Branch("DSTAU",      &DSTAU);
	tout->Branch("DSDIRA",      &DSDIRA);
	tout->Branch("DSVTXCHI2",  &DSVTXCHI2);
	tout->Branch("DSVTXNDOF",  &DSVTXNDOF);
	tout->Branch("DSDRJET",    &DSDRJET);
	tout->Branch("DSDRTAG",    &DSDRTAG);
	tout->Branch("DSPIPNN",    &DSPIPNN);
	tout->Branch("DSK1PNN",    &DSK1PNN);
	tout->Branch("DSK2PNN",    &DSK2PNN);
	tout->Branch("DSPHIM",     &DSPHIM);

	tout->Branch("DSP",        &DSP);
	tout->Branch("DSPT",       &DSPT);
	tout->Branch("DSPIPT",     &DSPIPT);
	tout->Branch("DSK1PT",     &DSK1PT);
	tout->Branch("DSK2PT",     &DSK2PT);
	tout->Branch("DSPHIPT",    &DSPHIPT);

	tout->Branch("DSTRK0",     &DSTRK0);
	tout->Branch("DSTRK1",     &DSTRK1);
	tout->Branch("DSTRK2",     &DSTRK2);

	tout->Branch("DSNJ",       &DSNJ);
	tout->Branch("DSMAXDR",    &DSMAXDR);

	tout->Branch("DSTRUE",    &DSTRUE);
	tout->Branch("DSTRUEDR",  &DSTRUEDR);
	tout->Branch("DSTRUEIDX", &DSTRUEIDX);
	tout->Branch("DSFROMB",   &DSFROMB );

	LCM       = new std::vector<double>();
	LCPX      = new std::vector<double>();
	LCPY      = new std::vector<double>();
	LCPZ      = new std::vector<double>();
	LCE       = new std::vector<double>();
	LCX       = new std::vector<double>();
	LCY       = new std::vector<double>();
	LCZ       = new std::vector<double>();
	LCIP      = new std::vector<double>();
	LCIPCHI2  = new std::vector<double>();
	LCFD      = new std::vector<double>();
	LCFDCHI2  = new std::vector<double>();
	LCTAU     = new std::vector<double>();
	LCDIRA    = new std::vector<double>();
	LCVTXCHI2 = new std::vector<double>();
	LCVTXNDOF = new std::vector<double>();
	LCDRJET   = new std::vector<double>();
	LCDRTAG   = new std::vector<double>();
	LCPPNN    = new std::vector<double>();
	LCKPNN    = new std::vector<double>();
	LCPIPNN   = new std::vector<double>();

	LCP       = new std::vector<double>();
	LCPT      = new std::vector<double>();
	LCPPT     = new std::vector<double>();
	LCKPT     = new std::vector<double>();
	LCPIPT    = new std::vector<double>();

	LCTRK0    = new std::vector<double>();
	LCTRK1    = new std::vector<double>();
	LCTRK2    = new std::vector<double>();

	LCNJ      = new std::vector<double>();
	LCMAXDR   = new std::vector<double>();

	tout->Branch("LCM",        &LCM);
	tout->Branch("LCPX",       &LCPX);
	tout->Branch("LCPY",       &LCPY);
	tout->Branch("LCPZ",       &LCPZ);
	tout->Branch("LCE",        &LCE);
	tout->Branch("LCX",        &LCX);
	tout->Branch("LCY",        &LCY);
	tout->Branch("LCZ",        &LCZ);
	tout->Branch("LCIP",       &LCIP);
	tout->Branch("LCIPCHI2",   &LCIPCHI2);
	tout->Branch("LCFD",       &LCFD);
	tout->Branch("LCFDCHI2",   &LCFDCHI2);
	tout->Branch("LCTAU",      &LCTAU);
	tout->Branch("LCDIRA",      &LCDIRA);
	tout->Branch("LCVTXCHI2",  &LCVTXCHI2);
	tout->Branch("LCVTXNDOF",  &LCVTXNDOF);
	tout->Branch("LCDRJET",    &LCDRJET);
	tout->Branch("LCDRTAG",    &LCDRTAG);
	tout->Branch("LCPPNN",     &LCPPNN);
	tout->Branch("LCKPNN",     &LCKPNN);
	tout->Branch("LCPIPNN",    &LCPIPNN);

	tout->Branch("LCP",        &LCP);
	tout->Branch("LCPT",       &LCPT);
	tout->Branch("LCPPT",      &LCPPT);
	tout->Branch("LCKPT",      &LCKPT);
	tout->Branch("LCPIPT",     &LCPIPT);

	tout->Branch("LCTRK0",     &LCTRK0);
	tout->Branch("LCTRK1",     &LCTRK1);
	tout->Branch("LCTRK2",     &LCTRK2);

	tout->Branch("LCNJ",       &LCNJ);
	tout->Branch("LCMAXDR",    &LCMAXDR);

	K3PIM       = new std::vector<double>();
	K3PIPX      = new std::vector<double>();
	K3PIPY      = new std::vector<double>();
	K3PIPZ      = new std::vector<double>();
	K3PIE       = new std::vector<double>();
	K3PIX       = new std::vector<double>();
	K3PIY       = new std::vector<double>();
	K3PIZ       = new std::vector<double>();
	K3PIIP      = new std::vector<double>();
	K3PIIPCHI2  = new std::vector<double>();
	K3PIFD      = new std::vector<double>();
	K3PIFDCHI2  = new std::vector<double>();
	K3PITAU     = new std::vector<double>();
	K3PIDIRA    = new std::vector<double>();
	K3PIVTXCHI2 = new std::vector<double>();
	K3PIVTXNDOF = new std::vector<double>();
	K3PIDRJET   = new std::vector<double>();
	K3PIDRTAG   = new std::vector<double>();
	K3PIPI1PNN  = new std::vector<double>();
	K3PIPI2PNN  = new std::vector<double>();
	K3PIPI3PNN  = new std::vector<double>();
	K3PIKPNN    = new std::vector<double>();

	K3PIP       = new std::vector<double>();
	K3PIPT      = new std::vector<double>();
	K3PIPI1PT   = new std::vector<double>();
	K3PIPI2PT   = new std::vector<double>();
	K3PIPI3PT   = new std::vector<double>();
	K3PIKPT     = new std::vector<double>();

	K3PITRK0    = new std::vector<double>();
	K3PITRK1    = new std::vector<double>();
	K3PITRK2    = new std::vector<double>();
	K3PITRK3    = new std::vector<double>();

	K3PINJ      = new std::vector<double>();
	K3PIMAXDR   = new std::vector<double>();

	tout->Branch("K3PIM",        &K3PIM);
	tout->Branch("K3PIPX",       &K3PIPX);
	tout->Branch("K3PIPY",       &K3PIPY);
	tout->Branch("K3PIPZ",       &K3PIPZ);
	tout->Branch("K3PIE",        &K3PIE);
	tout->Branch("K3PIX",        &K3PIX);
	tout->Branch("K3PIY",        &K3PIY);
	tout->Branch("K3PIZ",        &K3PIZ);
	tout->Branch("K3PIIP",       &K3PIIP);
	tout->Branch("K3PIIPCHI2",   &K3PIIPCHI2);
	tout->Branch("K3PIFD",       &K3PIFD);
	tout->Branch("K3PIFDCHI2",   &K3PIFDCHI2);
	tout->Branch("K3PITAU",      &K3PITAU);
	tout->Branch("K3PIDIRA",      &K3PIDIRA);
	tout->Branch("K3PIVTXCHI2",  &K3PIVTXCHI2);
	tout->Branch("K3PIVTXNDOF",  &K3PIVTXNDOF);
	tout->Branch("K3PIDRJET",    &K3PIDRJET);
	tout->Branch("K3PIDRTAG",    &K3PIDRTAG);
	tout->Branch("K3PIPI1PNN",   &K3PIPI1PNN);
	tout->Branch("K3PIPI2PNN",   &K3PIPI2PNN);
	tout->Branch("K3PIPI3PNN",   &K3PIPI3PNN);
	tout->Branch("K3PIKPNN",     &K3PIKPNN);

	tout->Branch("K3PIP",        &K3PIP);
	tout->Branch("K3PIPT",       &K3PIPT);
	tout->Branch("K3PIPI1PT",    &K3PIPI1PT);
	tout->Branch("K3PIPI2PT",    &K3PIPI2PT);
	tout->Branch("K3PIPI3PT",    &K3PIPI3PT);
	tout->Branch("K3PIKPT",      &K3PIKPT);

	tout->Branch("K3PITRK0",     &K3PITRK0);
	tout->Branch("K3PITRK1",     &K3PITRK1);
	tout->Branch("K3PITRK2",     &K3PITRK2);
	tout->Branch("K3PITRK3",     &K3PITRK3);

	tout->Branch("K3PINJ",       &K3PINJ);
	tout->Branch("K3PIMAXDR",    &K3PIMAXDR);

	JPSIM       = new std::vector<double>();
	JPSIPX      = new std::vector<double>();
	JPSIPY      = new std::vector<double>();
	JPSIPZ      = new std::vector<double>();
	JPSIE       = new std::vector<double>();
	JPSIX       = new std::vector<double>();
	JPSIY       = new std::vector<double>();
	JPSIZ       = new std::vector<double>();
	JPSIIP      = new std::vector<double>();
	JPSIIPCHI2  = new std::vector<double>();
	JPSIFD      = new std::vector<double>();
	JPSIFDCHI2  = new std::vector<double>();
	JPSITAU     = new std::vector<double>();
	JPSIDIRA    = new std::vector<double>();
	JPSIVTXCHI2 = new std::vector<double>();
	JPSIVTXNDOF = new std::vector<double>();
	JPSIDRJET   = new std::vector<double>();
	JPSIDRTAG   = new std::vector<double>();

	JPSIP       = new std::vector<double>();
	JPSIPT      = new std::vector<double>();
	JPSIPIPT    = new std::vector<double>();
	JPSIKPT     = new std::vector<double>();

	JPSITRK0    = new std::vector<double>();
	JPSITRK1    = new std::vector<double>();

	JPSINJ      = new std::vector<double>();
	JPSIMAXDR   = new std::vector<double>();

	JPSITRUE    = new std::vector<double>();
	JPSITRUEDR  = new std::vector<double>();
	JPSITRUEIDX = new std::vector<double>();
	JPSIFROMB   = new std::vector<double>();

	tout->Branch("JPSIM",        &JPSIM);
	tout->Branch("JPSIPX",       &JPSIPX);
	tout->Branch("JPSIPY",       &JPSIPY);
	tout->Branch("JPSIPZ",       &JPSIPZ);
	tout->Branch("JPSIE",        &JPSIE);
	tout->Branch("JPSIX",        &JPSIX);
	tout->Branch("JPSIY",        &JPSIY);
	tout->Branch("JPSIZ",        &JPSIZ);
	tout->Branch("JPSIIP",       &JPSIIP);
	tout->Branch("JPSIIPCHI2",   &JPSIIPCHI2);
	tout->Branch("JPSIFD",       &JPSIFD);
	tout->Branch("JPSIFDCHI2",   &JPSIFDCHI2);
	tout->Branch("JPSITAU",      &JPSITAU);
	tout->Branch("JPSIDIRA",     &JPSIDIRA);
	tout->Branch("JPSIVTXCHI2",  &JPSIVTXCHI2);
	tout->Branch("JPSIVTXNDOF",  &JPSIVTXNDOF);
	tout->Branch("JPSIDRJET",    &JPSIDRJET);
	tout->Branch("JPSIDRTAG",    &JPSIDRTAG);

	tout->Branch("JPSIP",        &JPSIP);
	tout->Branch("JPSIPT",       &JPSIPT);
	tout->Branch("JPSIPIPT",     &JPSIPIPT);
	tout->Branch("JPSIKPT",      &JPSIKPT);

	tout->Branch("JPSITRK0",     &JPSITRK0);
	tout->Branch("JPSITRK1",     &JPSITRK1);

	tout->Branch("JPSINJ",       &JPSINJ);
	tout->Branch("JPSIMAXDR",    &JPSIMAXDR);

	tout->Branch("JPSITRUE",     &JPSITRUE);
	tout->Branch("JPSITRUEDR",   &JPSITRUEDR);
	tout->Branch("JPSITRUEIDX",  &JPSITRUEIDX);
	tout->Branch("JPSIFROMB",    &JPSIFROMB );

	Notify();
}

Bool_t skimTuples::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void skimTuples::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t skimTuples::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
void skimTuples::clearEventOutputs() {
	EVT = 0.;
	DEC = 0.;
	NPV = 0.;
	clearOutputs();

	usedTruthJets_.clear();
}
void skimTuples::clearOutputs() {
	JPX = 0.;
	JPY = 0.;
	JPZ = 0.;
	JE = 0.;
	JPT = 0.;
	JETA = 0.;
	JTRUEPX = 0.;
	JTRUEPY = 0.;
	JTRUEPZ = 0.;
	JTRUEE = 0.;
	JTRUEPT = 0.;
	JTRUEETA = 0.;
	JTRUEDR = 0.;
	JS1 = 0.;
	JS2 = 0.;
	JQ = 0.;
	JN = 0.;
	JNQ = 0.;
	JNN = 0.;
	JPTD = 0.;
	JPV = 0.;
	JTRIGPT = 0.;
	JTRIG10 = 0.;
	JTRIG17 = 0.;
	JTRIG60 = 0.;
	JTRUEb = 0.;
	JTRUEc = 0.;
	JTRUED0 = 0.;
	JTRUEDP = 0.;
	JTRUEDS = 0.;
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
	TPX = 0.;
	TPY = 0.;
	TPZ = 0.;
	TE = 0.;
	TPT = 0.;
	TETA = 0.;
	TTRUEPX = 0.;
	TTRUEPY = 0.;
	TTRUEPZ = 0.;
	TTRUEE = 0.;
	TTRUEPT = 0.;
	TTRUEETA = 0.;
	TTRUEDR = 0.;
	TTRIGPT = 0.;
	TTRIG10 = 0.;
	TTRIG17 = 0.;
	TTRIG60 = 0.;
	DR = 0.;
	DPHI = 0.;
	ZM  =0.;
	ZP  =0.;
	ZPX =0.;
	ZPY =0.;
	ZPZ =0.;
	ZPT =0.;
	ZE  =0.;
	ZETA=0.;
	ZY=0.;
	ZDR =0.;
	MU0PX =0.;
	MU0PY =0.;
	MU0PZ =0.;
	MU0PT =0.;
	MU0IP =0.;
	MU0IPCHI2 =0.;
	MU0FPT =0.;
	MU0DR =0.;
	MU1PX =0.;
	MU1PY =0.;
	MU1PZ =0.;
	MU1PT =0.;
	MU1IP =0.;
	MU1IPCHI2 =0.;
	MU1FPT =0.;
	MU1DR =0.;
	ZTRUEM  =0.;
	ZTRUEP  =0.;
	ZTRUEPX =0.;
	ZTRUEPY =0.;
	ZTRUEPZ =0.;
	ZTRUEPT =0.;
	ZTRUEE  =0.;
	ZTRUEETA=0.;
	ZTRUEDR =0.;
	MU0TRUEPX =0.;
	MU0TRUEPY =0.;
	MU0TRUEPZ =0.;
	MU0TRUEPT =0.;
	MU1TRUEPX =0.;
	MU1TRUEPY =0.;
	MU1TRUEPZ =0.;
	MU1TRUEPT =0.;
	clearOutputVectors();
}
void skimTuples::clearOutputVectors() {
	SVM        ->clear();
	SVMCOR     ->clear();
	SVMCORERR  ->clear();
	SVMINPERP  ->clear();
	SVPT       ->clear();
	SVDRJ      ->clear();
	SVDRT      ->clear();
	SVN        ->clear();
	SVNTRUE    ->clear();
	SVNLINKED  ->clear();
	SVNJ       ->clear();
	SVNT       ->clear();
	SVQ        ->clear();
	SVPERP     ->clear();
	SVETA      ->clear();
	SVTZ       ->clear();
	SVVXCHI2   ->clear();
	SVIPCHI2   ->clear();
	SVMINIPCHI2->clear();
	SVPX       ->clear();
	SVPY       ->clear();
	SVPZ       ->clear();
	SVE        ->clear();
	SVMAXGHOST ->clear();
	SVX        ->clear();
	SVY        ->clear();
	SVZ        ->clear();
	SVSUMIPCHI2->clear();
	SVISD0     ->clear();
	SVISDP     ->clear();
	SVD0M      ->clear();
	SVDPM      ->clear();
	SVTRK0IDX  ->clear();
	SVTRK1IDX  ->clear();
	SVTRK2IDX  ->clear();
	SVTRK3IDX  ->clear();
	SVTRUEIDX  ->clear();
	SVTRUETRK0IDX  ->clear();
	SVTRUETRK1IDX  ->clear();
	SVTRUETRK2IDX  ->clear();
	SVTRUETRK3IDX  ->clear();
	SVTRK0P  ->clear();
	SVTRK1P  ->clear();
	SVTRK2P  ->clear();
	SVTRK3P  ->clear();
	SVTRK0PT  ->clear();
	SVTRK1PT  ->clear();
	SVTRK2PT  ->clear();
	SVTRK3PT  ->clear();
	SVTRK0PX  ->clear();
	SVTRK1PX  ->clear();
	SVTRK2PX  ->clear();
	SVTRK3PX  ->clear();
	SVTRK0PY  ->clear();
	SVTRK1PY  ->clear();
	SVTRK2PY  ->clear();
	SVTRK3PY  ->clear();
	SVTRK0PZ  ->clear();
	SVTRK1PZ  ->clear();
	SVTRK2PZ  ->clear();
	SVTRK3PZ  ->clear();
	SVTRK0PNNPI->clear();
	SVTRK1PNNPI->clear();
	SVTRK2PNNPI->clear();
	SVTRK3PNNPI->clear();
	SVTRK0PNNK->clear();
	SVTRK1PNNK->clear();
	SVTRK2PNNK->clear();
	SVTRK3PNNK->clear();

	TSVM        ->clear();
	TSVMCOR     ->clear();
	TSVMCORERR  ->clear();
	TSVMINPERP  ->clear();
	TSVPT       ->clear();
	TSVDRJ      ->clear();
	TSVDRT      ->clear();
	TSVN        ->clear();
	TSVNTRUE    ->clear();
	TSVNJ       ->clear();
	TSVNT       ->clear();
	TSVQ        ->clear();
	TSVPERP     ->clear();
	TSVETA      ->clear();
	TSVTZ       ->clear();
	TSVIPCHI2   ->clear();
	TSVMINIPCHI2->clear();
	TSVPX       ->clear();
	TSVPY       ->clear();
	TSVPZ       ->clear();
	TSVE        ->clear();
	TSVMAXGHOST ->clear();
	TSVX        ->clear();
	TSVY        ->clear();
	TSVZ        ->clear();
	TSVSUMIPCHI2->clear();
	TSVISD0     ->clear();
	TSVISDP     ->clear();
	TSVD0M      ->clear();
	TSVDPM      ->clear();
	TSVTRUEIDX  ->clear();

	TRKPT      ->clear();
	TRKIPCHI2  ->clear();
	TRKINJET   ->clear();

	TRUEBID   ->clear();
	TRUEBPX   ->clear();
	TRUEBPY   ->clear();
	TRUEBPZ   ->clear();
	TRUEBPT   ->clear();
	TRUEBE    ->clear();
	TRUEBX   ->clear();
	TRUEBY   ->clear();
	TRUEBZ   ->clear();
	TRUEDID   ->clear();
	TRUEDPX   ->clear();
	TRUEDPY   ->clear();
	TRUEDPZ   ->clear();
	TRUEDPT   ->clear();
	TRUEDE    ->clear();
	TRUEDX   ->clear();
	TRUEDY   ->clear();
	TRUEDZ   ->clear();
	TRUEDFROMB->clear();
	TRUEDTRUEB->clear();
	TRUEDSEL  ->clear();
	TRUEDTRK0IDX ->clear();
	TRUEDTRK0ID  ->clear();
	TRUEDTRK0P   ->clear();
	TRUEDTRK0PT  ->clear();
	TRUEDTRK0INACC ->clear();
	TRUEDTRK0RECO  ->clear();
	TRUEDTRK0PNNK  ->clear();
	TRUEDTRK0PNNPI ->clear();
	TRUEDTRK1IDX ->clear();
	TRUEDTRK1ID  ->clear();
	TRUEDTRK1P   ->clear();
	TRUEDTRK1PT  ->clear();
	TRUEDTRK1INACC ->clear();
	TRUEDTRK1RECO  ->clear();
	TRUEDTRK1PNNK  ->clear();
	TRUEDTRK1PNNPI ->clear();
	TRUEDTRK2IDX ->clear();
	TRUEDTRK2ID  ->clear();
	TRUEDTRK2P   ->clear();
	TRUEDTRK2PT  ->clear();
	TRUEDTRK2INACC ->clear();
	TRUEDTRK2RECO  ->clear();
	TRUEDTRK2PNNK  ->clear();
	TRUEDTRK2PNNPI ->clear();
	TRUEDTRK3IDX ->clear();
	TRUEDTRK3ID  ->clear();
	TRUEDTRK3P   ->clear();
	TRUEDTRK3PT  ->clear();
	TRUEDTRK3INACC ->clear();
	TRUEDTRK3RECO  ->clear();
	TRUEDTRK3PNNK  ->clear();
	TRUEDTRK3PNNPI ->clear();
	TRUEDTRK0PX   ->clear();
	TRUEDTRK1PX   ->clear();
	TRUEDTRK2PX   ->clear();
	TRUEDTRK3PX   ->clear();
	TRUEDTRK0PY   ->clear();
	TRUEDTRK1PY   ->clear();
	TRUEDTRK2PY   ->clear();
	TRUEDTRK3PY   ->clear();
	TRUEDTRK0PZ   ->clear();
	TRUEDTRK1PZ   ->clear();
	TRUEDTRK2PZ   ->clear();
	TRUEDTRK3PZ   ->clear();

	D0M       ->clear();
	D0PX      ->clear();
	D0PY      ->clear();
	D0PZ      ->clear();
	D0E       ->clear();
	D0X       ->clear();
	D0Y       ->clear();
	D0Z       ->clear();
	D0IP      ->clear();
	D0IPCHI2  ->clear();
	D0FD      ->clear();
	D0FDCHI2  ->clear();
	D0TAU     ->clear();
	D0DIRA    ->clear();
	D0VTXCHI2 ->clear();
	D0VTXNDOF ->clear();
	D0DRJET   ->clear();
	D0DRTAG   ->clear();
	D0PIPNNK  ->clear();
	D0KPNNK   ->clear();
	D0PIPNNPI ->clear();
	D0KPNNPI  ->clear();
	D0PIWEIGHT->clear();
	D0KWEIGHT ->clear();

	D0P       ->clear();
	D0PT      ->clear();
	D0ETA     ->clear();
	D0PIP     ->clear();
	D0KP      ->clear();
	D0PIPT    ->clear();
	D0KPT     ->clear();
	D0PIETA   ->clear();
	D0KETA    ->clear();
	D0PIPX    ->clear();
	D0KPX     ->clear();
	D0PIPY    ->clear();
	D0KPY     ->clear();
	D0PIPZ    ->clear();
	D0KPZ     ->clear();
	D0PIIPCHI2->clear();
	D0KIPCHI2 ->clear();

	D0TRK0    ->clear();
	D0TRK1    ->clear();
	D0TRUETRK0->clear();
	D0TRUETRK1->clear();

	D0NJ      ->clear();
	D0MAXDR   ->clear();

	D0TRUE    ->clear();
	D0TRUEDR  ->clear();
	D0TRUEIDX ->clear();
	D0FROMB   ->clear();

	DPM       ->clear();
	DPPX      ->clear();
	DPPY      ->clear();
	DPPZ      ->clear();
	DPE       ->clear();
	DPX       ->clear();
	DPY       ->clear();
	DPZ       ->clear();
	DPIP      ->clear();
	DPIPCHI2  ->clear();
	DPFD      ->clear();
	DPFDCHI2  ->clear();
	DPTAU     ->clear();
	DPDIRA    ->clear();
	DPVTXCHI2 ->clear();
	DPVTXNDOF ->clear();
	DPDRJET   ->clear();
	DPDRTAG   ->clear();
	DPPI1PNN  ->clear();
	DPPI2PNN  ->clear();
	DPKPNN    ->clear();

	DPP       ->clear();
	DPPT      ->clear();
	DPPI1PT   ->clear();
	DPPI2PT   ->clear();
	DPKPT     ->clear();

	DPTRK0    ->clear();
	DPTRK1    ->clear();
	DPTRK2    ->clear();

	DPNJ      ->clear();
	DPMAXDR   ->clear();

	DPTRUE    ->clear();
	DPTRUEDR  ->clear();
	DPTRUEIDX ->clear();
	DPFROMB   ->clear();

	DSM       ->clear();
	DSPX      ->clear();
	DSPY      ->clear();
	DSPZ      ->clear();
	DSE       ->clear();
	DSX       ->clear();
	DSY       ->clear();
	DSZ       ->clear();
	DSIP      ->clear();
	DSIPCHI2  ->clear();
	DSFD      ->clear();
	DSFDCHI2  ->clear();
	DSTAU     ->clear();
	DSDIRA    ->clear();
	DSVTXCHI2 ->clear();
	DSVTXNDOF ->clear();
	DSDRJET   ->clear();
	DSDRTAG   ->clear();
	DSPIPNN   ->clear();
	DSK1PNN   ->clear();
	DSK2PNN   ->clear();
	DSPHIM    ->clear();

	DSP       ->clear();
	DSPT      ->clear();
	DSPIPT    ->clear();
	DSK1PT    ->clear();
	DSK2PT    ->clear();
	DSPHIPT   ->clear();

	DSTRK0    ->clear();
	DSTRK1    ->clear();
	DSTRK2    ->clear();

	DSNJ      ->clear();
	DSMAXDR   ->clear();

	DSTRUE    ->clear();
	DSTRUEDR  ->clear();
	DSTRUEIDX ->clear();
	DSFROMB   ->clear();

	LCM       ->clear();
	LCPX      ->clear();
	LCPY      ->clear();
	LCPZ      ->clear();
	LCE       ->clear();
	LCX       ->clear();
	LCY       ->clear();
	LCZ       ->clear();
	LCIP      ->clear();
	LCIPCHI2  ->clear();
	LCFD      ->clear();
	LCFDCHI2  ->clear();
	LCTAU     ->clear();
	LCDIRA    ->clear();
	LCVTXCHI2 ->clear();
	LCVTXNDOF ->clear();
	LCDRJET   ->clear();
	LCDRTAG   ->clear();
	LCPPNN    ->clear();
	LCKPNN    ->clear();
	LCPIPNN   ->clear();

	LCP       ->clear();
	LCPT      ->clear();
	LCPPT     ->clear();
	LCKPT     ->clear();
	LCPIPT    ->clear();

	LCTRK0    ->clear();
	LCTRK1    ->clear();
	LCTRK2    ->clear();

	LCNJ      ->clear();
	LCMAXDR   ->clear();

	K3PIM       ->clear();
	K3PIPX      ->clear();
	K3PIPY      ->clear();
	K3PIPZ      ->clear();
	K3PIE       ->clear();
	K3PIX       ->clear();
	K3PIY       ->clear();
	K3PIZ       ->clear();
	K3PIIP      ->clear();
	K3PIIPCHI2  ->clear();
	K3PIFD      ->clear();
	K3PIFDCHI2  ->clear();
	K3PITAU     ->clear();
	K3PIDIRA    ->clear();
	K3PIVTXCHI2 ->clear();
	K3PIVTXNDOF ->clear();
	K3PIDRJET   ->clear();
	K3PIDRTAG   ->clear();
	K3PIPI1PNN  ->clear();
	K3PIPI2PNN  ->clear();
	K3PIPI3PNN  ->clear();
	K3PIKPNN    ->clear();

	K3PIP       ->clear();
	K3PIPT      ->clear();
	K3PIPI1PT   ->clear();
	K3PIPI2PT   ->clear();
	K3PIPI3PT   ->clear();
	K3PIKPT     ->clear();

	K3PITRK0    ->clear();
	K3PITRK1    ->clear();
	K3PITRK2    ->clear();

	K3PINJ      ->clear();
	K3PIMAXDR   ->clear();
	K3PITRK3    ->clear();

	JPSIM       ->clear();
	JPSIPX      ->clear();
	JPSIPY      ->clear();
	JPSIPZ      ->clear();
	JPSIE       ->clear();
	JPSIX       ->clear();
	JPSIY       ->clear();
	JPSIZ       ->clear();
	JPSIIP      ->clear();
	JPSIIPCHI2  ->clear();
	JPSIFD      ->clear();
	JPSIFDCHI2  ->clear();
	JPSITAU     ->clear();
	JPSIDIRA    ->clear();
	JPSIVTXCHI2 ->clear();
	JPSIVTXNDOF ->clear();
	JPSIDRJET   ->clear();
	JPSIDRTAG   ->clear();

	JPSIP       ->clear();
	JPSIPT      ->clear();
	JPSIPIPT    ->clear();
	JPSIKPT     ->clear();

	JPSITRK0    ->clear();
	JPSITRK1    ->clear();

	JPSINJ      ->clear();
	JPSIMAXDR   ->clear();

	JPSITRUE    ->clear();
	JPSITRUEDR  ->clear();
	JPSITRUEIDX ->clear();
	JPSIFROMB   ->clear();
}
#endif // #ifdef skimTuples_cxx
