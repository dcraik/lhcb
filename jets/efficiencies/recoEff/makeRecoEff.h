//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  4 23:02:48 2017 by ROOT version 5.34/36
// from TTree data/data
// found on file: /eos/lhcb/user/d/dcraik/jets/361/0/output.root
//////////////////////////////////////////////////////////

#ifndef makeRecoEff_h
#define makeRecoEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#include <boost/progress.hpp>

// Fixed size dimensions of array or collections stored in the TTree if any.

class makeRecoEff {
	public :
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

		makeRecoEff();
		virtual ~makeRecoEff();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop(int nmax=-1);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makeRecoEff_cxx
makeRecoEff::makeRecoEff() : fChain(0)
{
	char str[256];
	TChain *t = new TChain("data");
	boost::progress_display show_addfile_progress( 661 );
	for(int i=0; i<30; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/608/%d/output.root",i);//214
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<31; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/609/%d/output.root",i);//214
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<200; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/546/%d/output.root",i);//380
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<200; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/547/%d/output.root",i);//381
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/564/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/568/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/565/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/569/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/566/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/570/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/559/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/562/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/567/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/571/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/596/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/600/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/597/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/601/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/598/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/602/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/560/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/user/d/dcraik/jets/563/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/599/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<10; ++i) {
		++show_addfile_progress;
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/603/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	std::cout << t->GetEntries() << std::endl;
	Init(t);
}

makeRecoEff::~makeRecoEff()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t makeRecoEff::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t makeRecoEff::LoadTree(Long64_t entry)
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

void makeRecoEff::Init(TTree *tree)
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

	Notify();
}

Bool_t makeRecoEff::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void makeRecoEff::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t makeRecoEff::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef makeRecoEff_cxx
