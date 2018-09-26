//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 29 21:32:03 2017 by ROOT version 6.06/00
// from TTree data/data
// found on file: /eos/lhcb/user/d/dcraik/jets/222/1/output.root
//////////////////////////////////////////////////////////

#ifndef charmEfficiencies_h
#define charmEfficiencies_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TSystem.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

#include <boost/progress.hpp>

class charmEfficiencies {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		std::vector<double>  *gen_idx_pvr;
		std::vector<double>  *gen_idx_jet;
		std::vector<double>  *gen_idx_prnt;
		std::vector<double>  *gen_pid;
		std::vector<double>  *gen_q;
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
		std::vector<double>  *svr_px;
		std::vector<double>  *svr_py;
		std::vector<double>  *svr_pz;
		std::vector<double>  *svr_e;
		std::vector<double>  *svr_x;
		std::vector<double>  *svr_y;
		std::vector<double>  *svr_z;
		std::vector<double>  *svr_m;
		std::vector<double>  *svr_m_cor;
		std::vector<double>  *svr_pt;
		std::vector<double>  *svr_fd_min;
		std::vector<double>  *svr_fd_chi2;
		std::vector<double>  *svr_chi2;
		std::vector<double>  *svr_ip_chi2_sum;
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
		std::vector<double>  *jet_px;
		std::vector<double>  *jet_py;
		std::vector<double>  *jet_pz;
		std::vector<double>  *jet_e;
		std::vector<double>  *trk_idx_gen;
		std::vector<double>  *trk_idx_pvr;
		std::vector<double>  *trk_idx_jet;
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
		std::vector<double>  *neu_px;
		std::vector<double>  *neu_py;
		std::vector<double>  *neu_pz;
		std::vector<double>  *neu_e;
		std::vector<double>  *neu_pid;
		std::vector<double>  *d0_idx_pvr;
		std::vector<double>  *d0_idx_jet;
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
		std::vector<double>  *d02k3pi_idx_pvr;
		std::vector<double>  *d02k3pi_idx_jet;
		std::vector<double>  *d02k3pi_px;
		std::vector<double>  *d02k3pi_py;
		std::vector<double>  *d02k3pi_pz;
		std::vector<double>  *d02k3pi_e;
		std::vector<double>  *d02k3pi_x;
		std::vector<double>  *d02k3pi_y;
		std::vector<double>  *d02k3pi_z;
		std::vector<double>  *d02k3pi_m;
		std::vector<double>  *d02k3pi_ip;
		std::vector<double>  *d02k3pi_ip_chi2;
		std::vector<double>  *d02k3pi_vtx_chi2;
		std::vector<double>  *d02k3pi_vtx_ndof;
		std::vector<double>  *d02k3pi_fd;
		std::vector<double>  *d02k3pi_fd_chi2;
		std::vector<double>  *d02k3pi_tau;
		std::vector<double>  *d02k3pi_tau_err;
		std::vector<double>  *d02k3pi_tau_chi2;
		std::vector<double>  *d02k3pi_ntrk_jet;
		std::vector<double>  *d02k3pi_idx_trk0;
		std::vector<double>  *d02k3pi_idx_trk1;
		std::vector<double>  *d02k3pi_idx_trk2;
		std::vector<double>  *d02k3pi_idx_trk3;
		//std::vector<double>  *evt_dec;
		//std::vector<double>  *evt_j1_idx;
		//std::vector<double>  *evt_j1_dR;
		//std::vector<double>  *evt_j1_nsv;
		//std::vector<double>  *evt_j1_nmu;
		//std::vector<double>  *evt_j1_ntrk;
		//std::vector<double>  *evt_j1_nneu;
		//std::vector<double>  *evt_j1_px;
		//std::vector<double>  *evt_j1_py;
		//std::vector<double>  *evt_j1_pz;
		//std::vector<double>  *evt_j1_e;
		//std::vector<double>  *evt_j2_idx;
		//std::vector<double>  *evt_j2_dR;
		//std::vector<double>  *evt_j2_nsv;
		//std::vector<double>  *evt_j2_nmu;
		//std::vector<double>  *evt_j2_ntrk;
		//std::vector<double>  *evt_j2_nneu;
		//std::vector<double>  *evt_j2_px;
		//std::vector<double>  *evt_j2_py;
		//std::vector<double>  *evt_j2_pz;
		//std::vector<double>  *evt_j2_e;
		//Double_t        evt_pvr_n;

		// List of branches
		TBranch        *b_gen_idx_pvr;   //!
		TBranch        *b_gen_idx_jet;   //!
		TBranch        *b_gen_idx_prnt;   //!
		TBranch        *b_gen_pid;   //!
		TBranch        *b_gen_q;   //!
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
		TBranch        *b_svr_px;   //!
		TBranch        *b_svr_py;   //!
		TBranch        *b_svr_pz;   //!
		TBranch        *b_svr_e;   //!
		TBranch        *b_svr_x;   //!
		TBranch        *b_svr_y;   //!
		TBranch        *b_svr_z;   //!
		TBranch        *b_svr_m;   //!
		TBranch        *b_svr_m_cor;   //!
		TBranch        *b_svr_pt;   //!
		TBranch        *b_svr_fd_min;   //!
		TBranch        *b_svr_fd_chi2;   //!
		TBranch        *b_svr_chi2;   //!
		TBranch        *b_svr_ip_chi2_sum;   //!
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
		TBranch        *b_jet_px;   //!
		TBranch        *b_jet_py;   //!
		TBranch        *b_jet_pz;   //!
		TBranch        *b_jet_e;   //!
		TBranch        *b_trk_idx_gen;   //!
		TBranch        *b_trk_idx_pvr;   //!
		TBranch        *b_trk_idx_jet;   //!
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
		TBranch        *b_neu_px;   //!
		TBranch        *b_neu_py;   //!
		TBranch        *b_neu_pz;   //!
		TBranch        *b_neu_e;   //!
		TBranch        *b_neu_pid;   //!
		TBranch        *b_d0_idx_pvr;   //!
		TBranch        *b_d0_idx_jet;   //!
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
		TBranch        *b_d02k3pi_idx_pvr;   //!
		TBranch        *b_d02k3pi_idx_jet;   //!
		TBranch        *b_d02k3pi_px;   //!
		TBranch        *b_d02k3pi_py;   //!
		TBranch        *b_d02k3pi_pz;   //!
		TBranch        *b_d02k3pi_e;   //!
		TBranch        *b_d02k3pi_x;   //!
		TBranch        *b_d02k3pi_y;   //!
		TBranch        *b_d02k3pi_z;   //!
		TBranch        *b_d02k3pi_m;   //!
		TBranch        *b_d02k3pi_ip;   //!
		TBranch        *b_d02k3pi_ip_chi2;   //!
		TBranch        *b_d02k3pi_vtx_chi2;   //!
		TBranch        *b_d02k3pi_vtx_ndof;   //!
		TBranch        *b_d02k3pi_fd;   //!
		TBranch        *b_d02k3pi_fd_chi2;   //!
		TBranch        *b_d02k3pi_tau;   //!
		TBranch        *b_d02k3pi_tau_err;   //!
		TBranch        *b_d02k3pi_tau_chi2;   //!
		TBranch        *b_d02k3pi_ntrk_jet;   //!
		TBranch        *b_d02k3pi_idx_trk0;   //!
		TBranch        *b_d02k3pi_idx_trk1;   //!
		TBranch        *b_d02k3pi_idx_trk2;   //!
		TBranch        *b_d02k3pi_idx_trk3;   //!
		//TBranch        *b_evt_dec;   //!
		//TBranch        *b_evt_j1_idx;   //!
		//TBranch        *b_evt_j1_dR;   //!
		//TBranch        *b_evt_j1_nsv;   //!
		//TBranch        *b_evt_j1_nmu;   //!
		//TBranch        *b_evt_j1_ntrk;   //!
		//TBranch        *b_evt_j1_nneu;   //!
		//TBranch        *b_evt_j1_px;   //!
		//TBranch        *b_evt_j1_py;   //!
		//TBranch        *b_evt_j1_pz;   //!
		//TBranch        *b_evt_j1_e;   //!
		//TBranch        *b_evt_j2_idx;   //!
		//TBranch        *b_evt_j2_dR;   //!
		//TBranch        *b_evt_j2_nsv;   //!
		//TBranch        *b_evt_j2_nmu;   //!
		//TBranch        *b_evt_j2_ntrk;   //!
		//TBranch        *b_evt_j2_nneu;   //!
		//TBranch        *b_evt_j2_px;   //!
		//TBranch        *b_evt_j2_py;   //!
		//TBranch        *b_evt_j2_pz;   //!
		//TBranch        *b_evt_j2_e;   //!
		//TBranch        *b_evt_pvr_n;   //!

		TH2F* kaonPIDHist;
		TH2F* pionPIDHist;
		TH2F* protonPIDHist;

		int type_, origin_, nmax_;

		charmEfficiencies(int type=0, int origin=0, int nmax=-1);
		virtual ~charmEfficiencies();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		bool checkForD0(int idxD, int idx0, int idx1, int idxNxt, int origin, bool ordered=false);
		bool checkForDp(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered=false);
		bool checkForDs(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered=false);
		bool checkForLc(int idxD, int idx0, int idx1, int idx2, int idxNxt, int origin, bool ordered=false);

		bool inLHCb(int idx);

		double getProbPi(int idx);
		double getProbK(int idx);
		double getProbP(int idx);

		void setOrigin(int origin) {origin_=origin;}
};

#endif

#ifdef charmEfficiencies_cxx
charmEfficiencies::charmEfficiencies(int type, int origin, int nmax) : fChain(0), type_(type), origin_(origin), nmax_(nmax) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	TChain *t = new TChain("data");
	char str[256];
	int nFiles(0), nFound(0);

	if(type==0 || type==1) {//ccbar
		nFiles += 364;
		boost::progress_display show_addfile_progress( 364 );
		for(int i=0; i<364; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/235/%d/output.root",i);//214
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==11 || type==30) {//ccbar new
		nFiles += 364;
		boost::progress_display show_addfile_progress( 364 );
		for(int i=0; i<364; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/251/%d/output.root",i);//214
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==0 || type==2) {//bbbar
		nFiles += 40;
		boost::progress_display show_addfile_progress( 40 );
		for(int i=0; i<40; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/234/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==12 || type==30) {//bbbar
		nFiles += 40;
		boost::progress_display show_addfile_progress( 40 );
		for(int i=0; i<40; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/250/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==0 || type==3) {//D0
		nFiles += 44;
		boost::progress_display show_addfile_progress( 44 );
		for(int i=0; i<44; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/237/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==13 || type==30) {//D0
		nFiles += 44;
		boost::progress_display show_addfile_progress( 44 );
		for(int i=0; i<44; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/253/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==0 || type==4) {//D+
		nFiles += 182;
		boost::progress_display show_addfile_progress( 182 );
		for(int i=0; i<182; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/236/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==14 || type==30) {//D+
		nFiles += 182;
		boost::progress_display show_addfile_progress( 182 );
		for(int i=0; i<182; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/252/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==0 || type==5) {//Ds
		nFiles += 73;
		boost::progress_display show_addfile_progress( 73 );
		for(int i=0; i<73; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/238/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==15 || type==30) {//Ds
		nFiles += 73;
		boost::progress_display show_addfile_progress( 73 );
		for(int i=0; i<73; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/254/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==0 || type==6) {//Lc
		nFiles += 803;
		boost::progress_display show_addfile_progress( 803 );
		for(int i=0; i<803; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/241/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==10 || type==16 || type==30) {//Lc
		nFiles += 803;
		boost::progress_display show_addfile_progress( 803 );
		for(int i=0; i<803; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/255/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==20 || type==30) {//jets
		nFiles += 150;
		boost::progress_display show_addfile_progress( 150 );
		for(int j=256; j<266; ++j) {
			for(int i=0; i<15; ++i) {
				++show_addfile_progress;
				sprintf(str,"/eos/lhcb/user/d/dcraik/jets/%d/%d/output.root",j,i);
				if(gSystem->AccessPathName(str)) continue;
				t->Add(str);
				++nFound;
			}
		}
	}
	if(type==40 || type==41) {//D0 special
		nFiles += 170;
		boost::progress_display show_addfile_progress( 170 );
		for(int i=0; i<170; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/293/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==40 || type==42) {//Dp special
		nFiles += 170;
		boost::progress_display show_addfile_progress( 170 );
		for(int i=0; i<170; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/294/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	if(type==40 || type==43) {//Ds special
		nFiles += 170;
		boost::progress_display show_addfile_progress( 170 );
		for(int i=0; i<170; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/295/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			++nFound;
		}
	}
	//if(type<0 || type>6) {
	//	std::cout << "Unknown type " << type << std::endl;
	//}
	//job 250: found 28 of 40
	//job 251: found 200 of 364
	//job 252: found 182 of 182
	//job 253: found 42 of 44
	//job 254: found 43 of 73
	//job 255: found 725 of 803
	std::cout << "Found " << nFound << " of " << nFiles << " input files" << std::endl;
	std::cout << "Found " << t->GetEntries() << " entries" << std::endl;

	Init(t);

	TFile* fKpid  = TFile::Open("../pidcalib/PerfHists_K_Turbo16_MagDown_kaons3_Brunel_P_Brunel_PT.root");
	TFile* fpipid = TFile::Open("../pidcalib/PerfHists_Pi_Turbo16_MagDown_pions_Brunel_P_Brunel_PT.root");
	TFile* fppid = TFile::Open("../pidcalib/PerfHists_P_Turbo16_MagDown_protons6_Brunel_P_Brunel_PT.root");
	kaonPIDHist = ((TH2F*)fKpid->Get("K_Brunel_MC15TuneV1_ProbNNK > 0.3_All"));
	pionPIDHist = ((TH2F*)fpipid->Get("Pi_Brunel_MC15TuneV1_ProbNNpi > 0.2_All"));
	protonPIDHist = ((TH2F*)fppid->Get("P_Brunel_MC15TuneV1_ProbNNp > 0.3_All"));

	kaonPIDHist->Print();
	pionPIDHist->Print();
	protonPIDHist->Print();
}

charmEfficiencies::~charmEfficiencies()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t charmEfficiencies::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t charmEfficiencies::LoadTree(Long64_t entry)
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

void charmEfficiencies::Init(TTree *tree)
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
	svr_px = 0;
	svr_py = 0;
	svr_pz = 0;
	svr_e = 0;
	svr_x = 0;
	svr_y = 0;
	svr_z = 0;
	svr_m = 0;
	svr_m_cor = 0;
	svr_pt = 0;
	svr_fd_min = 0;
	svr_fd_chi2 = 0;
	svr_chi2 = 0;
	svr_ip_chi2_sum = 0;
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
	jet_px = 0;
	jet_py = 0;
	jet_pz = 0;
	jet_e = 0;
	trk_idx_gen = 0;
	trk_idx_pvr = 0;
	trk_idx_jet = 0;
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
	neu_px = 0;
	neu_py = 0;
	neu_pz = 0;
	neu_e = 0;
	neu_pid = 0;
	d0_idx_pvr = 0;
	d0_idx_jet = 0;
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
	d02k3pi_idx_pvr = 0;
	d02k3pi_idx_jet = 0;
	d02k3pi_px = 0;
	d02k3pi_py = 0;
	d02k3pi_pz = 0;
	d02k3pi_e = 0;
	d02k3pi_x = 0;
	d02k3pi_y = 0;
	d02k3pi_z = 0;
	d02k3pi_m = 0;
	d02k3pi_ip = 0;
	d02k3pi_ip_chi2 = 0;
	d02k3pi_vtx_chi2 = 0;
	d02k3pi_vtx_ndof = 0;
	d02k3pi_fd = 0;
	d02k3pi_fd_chi2 = 0;
	d02k3pi_tau = 0;
	d02k3pi_tau_err = 0;
	d02k3pi_tau_chi2 = 0;
	d02k3pi_ntrk_jet = 0;
	d02k3pi_idx_trk0 = 0;
	d02k3pi_idx_trk1 = 0;
	d02k3pi_idx_trk2 = 0;
	d02k3pi_idx_trk3 = 0;
	//evt_dec = 0;
	//evt_j1_idx = 0;
	//evt_j1_dR = 0;
	//evt_j1_nsv = 0;
	//evt_j1_nmu = 0;
	//evt_j1_ntrk = 0;
	//evt_j1_nneu = 0;
	//evt_j1_px = 0;
	//evt_j1_py = 0;
	//evt_j1_pz = 0;
	//evt_j1_e = 0;
	//evt_j2_idx = 0;
	//evt_j2_dR = 0;
	//evt_j2_nsv = 0;
	//evt_j2_nmu = 0;
	//evt_j2_ntrk = 0;
	//evt_j2_nneu = 0;
	//evt_j2_px = 0;
	//evt_j2_py = 0;
	//evt_j2_pz = 0;
	//evt_j2_e = 0;
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
	fChain->SetBranchAddress("svr_px", &svr_px, &b_svr_px);
	fChain->SetBranchAddress("svr_py", &svr_py, &b_svr_py);
	fChain->SetBranchAddress("svr_pz", &svr_pz, &b_svr_pz);
	fChain->SetBranchAddress("svr_e", &svr_e, &b_svr_e);
	fChain->SetBranchAddress("svr_x", &svr_x, &b_svr_x);
	fChain->SetBranchAddress("svr_y", &svr_y, &b_svr_y);
	fChain->SetBranchAddress("svr_z", &svr_z, &b_svr_z);
	fChain->SetBranchAddress("svr_m", &svr_m, &b_svr_m);
	fChain->SetBranchAddress("svr_m_cor", &svr_m_cor, &b_svr_m_cor);
	fChain->SetBranchAddress("svr_pt", &svr_pt, &b_svr_pt);
	fChain->SetBranchAddress("svr_fd_min", &svr_fd_min, &b_svr_fd_min);
	fChain->SetBranchAddress("svr_fd_chi2", &svr_fd_chi2, &b_svr_fd_chi2);
	fChain->SetBranchAddress("svr_chi2", &svr_chi2, &b_svr_chi2);
	fChain->SetBranchAddress("svr_ip_chi2_sum", &svr_ip_chi2_sum, &b_svr_ip_chi2_sum);
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
	fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
	fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
	fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
	fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
	fChain->SetBranchAddress("trk_idx_gen", &trk_idx_gen, &b_trk_idx_gen);
	fChain->SetBranchAddress("trk_idx_pvr", &trk_idx_pvr, &b_trk_idx_pvr);
	fChain->SetBranchAddress("trk_idx_jet", &trk_idx_jet, &b_trk_idx_jet);
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
	fChain->SetBranchAddress("neu_px", &neu_px, &b_neu_px);
	fChain->SetBranchAddress("neu_py", &neu_py, &b_neu_py);
	fChain->SetBranchAddress("neu_pz", &neu_pz, &b_neu_pz);
	fChain->SetBranchAddress("neu_e", &neu_e, &b_neu_e);
	fChain->SetBranchAddress("neu_pid", &neu_pid, &b_neu_pid);
	fChain->SetBranchAddress("d0_idx_pvr", &d0_idx_pvr, &b_d0_idx_pvr);
	fChain->SetBranchAddress("d0_idx_jet", &d0_idx_jet, &b_d0_idx_jet);
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
	fChain->SetBranchAddress("d02k3pi_idx_pvr", &d02k3pi_idx_pvr, &b_d02k3pi_idx_pvr);
	fChain->SetBranchAddress("d02k3pi_idx_jet", &d02k3pi_idx_jet, &b_d02k3pi_idx_jet);
	fChain->SetBranchAddress("d02k3pi_px", &d02k3pi_px, &b_d02k3pi_px);
	fChain->SetBranchAddress("d02k3pi_py", &d02k3pi_py, &b_d02k3pi_py);
	fChain->SetBranchAddress("d02k3pi_pz", &d02k3pi_pz, &b_d02k3pi_pz);
	fChain->SetBranchAddress("d02k3pi_e", &d02k3pi_e, &b_d02k3pi_e);
	fChain->SetBranchAddress("d02k3pi_x", &d02k3pi_x, &b_d02k3pi_x);
	fChain->SetBranchAddress("d02k3pi_y", &d02k3pi_y, &b_d02k3pi_y);
	fChain->SetBranchAddress("d02k3pi_z", &d02k3pi_z, &b_d02k3pi_z);
	fChain->SetBranchAddress("d02k3pi_m", &d02k3pi_m, &b_d02k3pi_m);
	fChain->SetBranchAddress("d02k3pi_ip", &d02k3pi_ip, &b_d02k3pi_ip);
	fChain->SetBranchAddress("d02k3pi_ip_chi2", &d02k3pi_ip_chi2, &b_d02k3pi_ip_chi2);
	fChain->SetBranchAddress("d02k3pi_vtx_chi2", &d02k3pi_vtx_chi2, &b_d02k3pi_vtx_chi2);
	fChain->SetBranchAddress("d02k3pi_vtx_ndof", &d02k3pi_vtx_ndof, &b_d02k3pi_vtx_ndof);
	fChain->SetBranchAddress("d02k3pi_fd", &d02k3pi_fd, &b_d02k3pi_fd);
	fChain->SetBranchAddress("d02k3pi_fd_chi2", &d02k3pi_fd_chi2, &b_d02k3pi_fd_chi2);
	fChain->SetBranchAddress("d02k3pi_tau", &d02k3pi_tau, &b_d02k3pi_tau);
	fChain->SetBranchAddress("d02k3pi_tau_err", &d02k3pi_tau_err, &b_d02k3pi_tau_err);
	fChain->SetBranchAddress("d02k3pi_tau_chi2", &d02k3pi_tau_chi2, &b_d02k3pi_tau_chi2);
	fChain->SetBranchAddress("d02k3pi_ntrk_jet", &d02k3pi_ntrk_jet, &b_d02k3pi_ntrk_jet);
	fChain->SetBranchAddress("d02k3pi_idx_trk0", &d02k3pi_idx_trk0, &b_d02k3pi_idx_trk0);
	fChain->SetBranchAddress("d02k3pi_idx_trk1", &d02k3pi_idx_trk1, &b_d02k3pi_idx_trk1);
	fChain->SetBranchAddress("d02k3pi_idx_trk2", &d02k3pi_idx_trk2, &b_d02k3pi_idx_trk2);
	fChain->SetBranchAddress("d02k3pi_idx_trk3", &d02k3pi_idx_trk3, &b_d02k3pi_idx_trk3);
	//fChain->SetBranchAddress("evt_dec", &evt_dec, &b_evt_dec);
	//fChain->SetBranchAddress("evt_j1_idx", &evt_j1_idx, &b_evt_j1_idx);
	//fChain->SetBranchAddress("evt_j1_dR", &evt_j1_dR, &b_evt_j1_dR);
	//fChain->SetBranchAddress("evt_j1_nsv", &evt_j1_nsv, &b_evt_j1_nsv);
	//fChain->SetBranchAddress("evt_j1_nmu", &evt_j1_nmu, &b_evt_j1_nmu);
	//fChain->SetBranchAddress("evt_j1_ntrk", &evt_j1_ntrk, &b_evt_j1_ntrk);
	//fChain->SetBranchAddress("evt_j1_nneu", &evt_j1_nneu, &b_evt_j1_nneu);
	//fChain->SetBranchAddress("evt_j1_px", &evt_j1_px, &b_evt_j1_px);
	//fChain->SetBranchAddress("evt_j1_py", &evt_j1_py, &b_evt_j1_py);
	//fChain->SetBranchAddress("evt_j1_pz", &evt_j1_pz, &b_evt_j1_pz);
	//fChain->SetBranchAddress("evt_j1_e", &evt_j1_e, &b_evt_j1_e);
	//fChain->SetBranchAddress("evt_j2_idx", &evt_j2_idx, &b_evt_j2_idx);
	//fChain->SetBranchAddress("evt_j2_dR", &evt_j2_dR, &b_evt_j2_dR);
	//fChain->SetBranchAddress("evt_j2_nsv", &evt_j2_nsv, &b_evt_j2_nsv);
	//fChain->SetBranchAddress("evt_j2_nmu", &evt_j2_nmu, &b_evt_j2_nmu);
	//fChain->SetBranchAddress("evt_j2_ntrk", &evt_j2_ntrk, &b_evt_j2_ntrk);
	//fChain->SetBranchAddress("evt_j2_nneu", &evt_j2_nneu, &b_evt_j2_nneu);
	//fChain->SetBranchAddress("evt_j2_px", &evt_j2_px, &b_evt_j2_px);
	//fChain->SetBranchAddress("evt_j2_py", &evt_j2_py, &b_evt_j2_py);
	//fChain->SetBranchAddress("evt_j2_pz", &evt_j2_pz, &b_evt_j2_pz);
	//fChain->SetBranchAddress("evt_j2_e", &evt_j2_e, &b_evt_j2_e);
	//fChain->SetBranchAddress("evt_pvr_n", &evt_pvr_n, &b_evt_pvr_n);
	Notify();
}

Bool_t charmEfficiencies::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void charmEfficiencies::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t charmEfficiencies::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef charmEfficiencies_cxx
