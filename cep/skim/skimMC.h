//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 23 12:03:22 2019 by ROOT version 6.12/06
// from TTree data/data
// found on file: ../davinci/output.root
//////////////////////////////////////////////////////////

#ifndef skimMC_h
#define skimMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class skimMC {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		std::vector<double>  *gen_idx_pvr;
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
		std::vector<double>  *trk_idx_gen;
		std::vector<double>  *trk_idx_pvr;
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
		std::vector<double>  *trk_l0_hadron_dec;
		std::vector<double>  *trk_l0_hadron_tis;
		std::vector<double>  *trk_l0_hadron_tos;
		std::vector<double>  *trk_l0_dihadron_lowmult_dec;
		std::vector<double>  *trk_l0_dihadron_lowmult_tis;
		std::vector<double>  *trk_l0_dihadron_lowmult_tos;
		std::vector<double>  *trk_l0_hrc_dihadron_lowmult_dec;
		std::vector<double>  *trk_l0_hrc_dihadron_lowmult_tis;
		std::vector<double>  *trk_l0_hrc_dihadron_lowmult_tos;
		std::vector<double>  *trk_hlt1_passthru_dec;
		std::vector<double>  *trk_hlt1_passthru_tis;
		std::vector<double>  *trk_hlt1_passthru_tos;
		std::vector<double>  *trk_hlt1_lowmult_dec;
		std::vector<double>  *trk_hlt1_lowmult_tis;
		std::vector<double>  *trk_hlt1_lowmult_tos;
		std::vector<double>  *trk_hlt1_herchel_dec;
		std::vector<double>  *trk_hlt1_herchel_tis;
		std::vector<double>  *trk_hlt1_herchel_tos;
		std::vector<double>  *trk_hlt1_velo_dec;
		std::vector<double>  *trk_hlt1_velo_tis;
		std::vector<double>  *trk_hlt1_velo_tos;
		std::vector<double>  *trk_hlt1_veloherchel_dec;
		std::vector<double>  *trk_hlt1_veloherchel_tis;
		std::vector<double>  *trk_hlt1_veloherchel_tos;
		std::vector<double>  *trk_hlt2_lmr2hh_dec;
		std::vector<double>  *trk_hlt2_lmr2hh_tis;
		std::vector<double>  *trk_hlt2_lmr2hh_tos;
		std::vector<double>  *trk_hlt2_lmr2hhws_dec;
		std::vector<double>  *trk_hlt2_lmr2hhws_tis;
		std::vector<double>  *trk_hlt2_lmr2hhws_tos;
		std::vector<double>  *neu_idx_gen;
		std::vector<double>  *neu_p;
		std::vector<double>  *neu_pt;
		std::vector<double>  *neu_px;
		std::vector<double>  *neu_py;
		std::vector<double>  *neu_pz;
		std::vector<double>  *neu_e;
		std::vector<double>  *neu_pid;
		std::vector<double>  *neu_l0_hadron_dec;
		std::vector<double>  *neu_l0_hadron_tis;
		std::vector<double>  *neu_l0_hadron_tos;
		std::vector<double>  *neu_l0_dihadron_lowmult_dec;
		std::vector<double>  *neu_l0_dihadron_lowmult_tis;
		std::vector<double>  *neu_l0_dihadron_lowmult_tos;
		std::vector<double>  *neu_l0_hrc_dihadron_lowmult_dec;
		std::vector<double>  *neu_l0_hrc_dihadron_lowmult_tis;
		std::vector<double>  *neu_l0_hrc_dihadron_lowmult_tos;
		std::vector<double>  *neu_hlt1_passthru_dec;
		std::vector<double>  *neu_hlt1_passthru_tis;
		std::vector<double>  *neu_hlt1_passthru_tos;
		std::vector<double>  *neu_hlt1_lowmult_dec;
		std::vector<double>  *neu_hlt1_lowmult_tis;
		std::vector<double>  *neu_hlt1_lowmult_tos;
		std::vector<double>  *neu_hlt1_herchel_dec;
		std::vector<double>  *neu_hlt1_herchel_tis;
		std::vector<double>  *neu_hlt1_herchel_tos;
		std::vector<double>  *neu_hlt1_velo_dec;
		std::vector<double>  *neu_hlt1_velo_tis;
		std::vector<double>  *neu_hlt1_velo_tos;
		std::vector<double>  *neu_hlt1_veloherchel_dec;
		std::vector<double>  *neu_hlt1_veloherchel_tis;
		std::vector<double>  *neu_hlt1_veloherchel_tos;
		std::vector<double>  *neu_hlt2_lmr2hh_dec;
		std::vector<double>  *neu_hlt2_lmr2hh_tis;
		std::vector<double>  *neu_hlt2_lmr2hh_tos;
		std::vector<double>  *neu_hlt2_lmr2hhws_dec;
		std::vector<double>  *neu_hlt2_lmr2hhws_tis;
		std::vector<double>  *neu_hlt2_lmr2hhws_tos;
		std::vector<double>  *phi_idx_pvr;
		std::vector<double>  *phi_p;
		std::vector<double>  *phi_pt;
		std::vector<double>  *phi_px;
		std::vector<double>  *phi_py;
		std::vector<double>  *phi_pz;
		std::vector<double>  *phi_e;
		std::vector<double>  *phi_x;
		std::vector<double>  *phi_y;
		std::vector<double>  *phi_z;
		std::vector<double>  *phi_m;
		std::vector<double>  *phi_ip;
		std::vector<double>  *phi_ip_chi2;
		std::vector<double>  *phi_vtx_chi2;
		std::vector<double>  *phi_vtx_ndof;
		std::vector<double>  *phi_fd;
		std::vector<double>  *phi_fd_chi2;
		std::vector<double>  *phi_tau;
		std::vector<double>  *phi_tau_err;
		std::vector<double>  *phi_tau_chi2;
		std::vector<double>  *phi_idx_trk0;
		std::vector<double>  *phi_idx_trk1;
		std::vector<double>  *phi_l0_hadron_dec;
		std::vector<double>  *phi_l0_hadron_tis;
		std::vector<double>  *phi_l0_hadron_tos;
		std::vector<double>  *phi_l0_dihadron_lowmult_dec;
		std::vector<double>  *phi_l0_dihadron_lowmult_tis;
		std::vector<double>  *phi_l0_dihadron_lowmult_tos;
		std::vector<double>  *phi_l0_hrc_dihadron_lowmult_dec;
		std::vector<double>  *phi_l0_hrc_dihadron_lowmult_tis;
		std::vector<double>  *phi_l0_hrc_dihadron_lowmult_tos;
		std::vector<double>  *phi_hlt1_passthru_dec;
		std::vector<double>  *phi_hlt1_passthru_tis;
		std::vector<double>  *phi_hlt1_passthru_tos;
		std::vector<double>  *phi_hlt1_lowmult_dec;
		std::vector<double>  *phi_hlt1_lowmult_tis;
		std::vector<double>  *phi_hlt1_lowmult_tos;
		std::vector<double>  *phi_hlt1_herchel_dec;
		std::vector<double>  *phi_hlt1_herchel_tis;
		std::vector<double>  *phi_hlt1_herchel_tos;
		std::vector<double>  *phi_hlt1_velo_dec;
		std::vector<double>  *phi_hlt1_velo_tis;
		std::vector<double>  *phi_hlt1_velo_tos;
		std::vector<double>  *phi_hlt1_veloherchel_dec;
		std::vector<double>  *phi_hlt1_veloherchel_tis;
		std::vector<double>  *phi_hlt1_veloherchel_tos;
		std::vector<double>  *phi_hlt2_lmr2hh_dec;
		std::vector<double>  *phi_hlt2_lmr2hh_tis;
		std::vector<double>  *phi_hlt2_lmr2hh_tos;
		std::vector<double>  *phi_hlt2_lmr2hhws_dec;
		std::vector<double>  *phi_hlt2_lmr2hhws_tis;
		std::vector<double>  *phi_hlt2_lmr2hhws_tos;
		Double_t        evt_pvr_n;
		Double_t        evt_neu_n;
		Double_t        evt_chg_n;
		Double_t        evt_trk_n;
		Double_t        evt_trk_n_velor;
		Double_t        evt_trk_n_velo;
		Double_t        evt_trk_n_long;
		Double_t        evt_trk_n_up;
		Double_t        evt_trk_n_down;
		Double_t        evt_trk_n_t;
		Double_t        evt_trk_n_mu;
		Double_t        evt_hrc_fom;

		//outputs
		Double_t Kplus_TRUEP;
		Double_t Kplus_TRUEPT;
		Double_t Kplus_TRUEPX;
		Double_t Kplus_TRUEPY;
		Double_t Kplus_TRUEPZ;
		Double_t Kplus_TRUEE;

		Double_t Kminus_TRUEP;
		Double_t Kminus_TRUEPT;
		Double_t Kminus_TRUEPX;
		Double_t Kminus_TRUEPY;
		Double_t Kminus_TRUEPZ;
		Double_t Kminus_TRUEE;

		Double_t phi_TRUEP;
		Double_t phi_TRUEPT;
		Double_t phi_TRUEPX;
		Double_t phi_TRUEPY;
		Double_t phi_TRUEPZ;
		Double_t phi_TRUEE;

		Double_t Kplus_P;
		Double_t Kplus_PT;
		Double_t Kplus_PX;
		Double_t Kplus_PY;
		Double_t Kplus_PZ;
		Double_t Kplus_E;

		Double_t Kminus_P;
		Double_t Kminus_PT;
		Double_t Kminus_PX;
		Double_t Kminus_PY;
		Double_t Kminus_PZ;
		Double_t Kminus_E;

		Double_t phi_P;
		Double_t phi_PT;
		Double_t phi_PX;
		Double_t phi_PY;
		Double_t phi_PZ;
		Double_t phi_E;
		Double_t phi_M;

		// List of branches
		TBranch        *b_gen_idx_pvr;   //!
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
		TBranch        *b_trk_idx_gen;   //!
		TBranch        *b_trk_idx_pvr;   //!
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
		TBranch        *b_trk_l0_hadron_dec;   //!
		TBranch        *b_trk_l0_hadron_tis;   //!
		TBranch        *b_trk_l0_hadron_tos;   //!
		TBranch        *b_trk_l0_dihadron_lowmult_dec;   //!
		TBranch        *b_trk_l0_dihadron_lowmult_tis;   //!
		TBranch        *b_trk_l0_dihadron_lowmult_tos;   //!
		TBranch        *b_trk_l0_hrc_dihadron_lowmult_dec;   //!
		TBranch        *b_trk_l0_hrc_dihadron_lowmult_tis;   //!
		TBranch        *b_trk_l0_hrc_dihadron_lowmult_tos;   //!
		TBranch        *b_trk_hlt1_passthru_dec;   //!
		TBranch        *b_trk_hlt1_passthru_tis;   //!
		TBranch        *b_trk_hlt1_passthru_tos;   //!
		TBranch        *b_trk_hlt1_lowmult_dec;   //!
		TBranch        *b_trk_hlt1_lowmult_tis;   //!
		TBranch        *b_trk_hlt1_lowmult_tos;   //!
		TBranch        *b_trk_hlt1_herchel_dec;   //!
		TBranch        *b_trk_hlt1_herchel_tis;   //!
		TBranch        *b_trk_hlt1_herchel_tos;   //!
		TBranch        *b_trk_hlt1_velo_dec;   //!
		TBranch        *b_trk_hlt1_velo_tis;   //!
		TBranch        *b_trk_hlt1_velo_tos;   //!
		TBranch        *b_trk_hlt1_veloherchel_dec;   //!
		TBranch        *b_trk_hlt1_veloherchel_tis;   //!
		TBranch        *b_trk_hlt1_veloherchel_tos;   //!
		TBranch        *b_trk_hlt2_lmr2hh_dec;   //!
		TBranch        *b_trk_hlt2_lmr2hh_tis;   //!
		TBranch        *b_trk_hlt2_lmr2hh_tos;   //!
		TBranch        *b_trk_hlt2_lmr2hhws_dec;   //!
		TBranch        *b_trk_hlt2_lmr2hhws_tis;   //!
		TBranch        *b_trk_hlt2_lmr2hhws_tos;   //!
		TBranch        *b_neu_idx_gen;   //!
		TBranch        *b_neu_p;   //!
		TBranch        *b_neu_pt;   //!
		TBranch        *b_neu_px;   //!
		TBranch        *b_neu_py;   //!
		TBranch        *b_neu_pz;   //!
		TBranch        *b_neu_e;   //!
		TBranch        *b_neu_pid;   //!
		TBranch        *b_neu_l0_hadron_dec;   //!
		TBranch        *b_neu_l0_hadron_tis;   //!
		TBranch        *b_neu_l0_hadron_tos;   //!
		TBranch        *b_neu_l0_dihadron_lowmult_dec;   //!
		TBranch        *b_neu_l0_dihadron_lowmult_tis;   //!
		TBranch        *b_neu_l0_dihadron_lowmult_tos;   //!
		TBranch        *b_neu_l0_hrc_dihadron_lowmult_dec;   //!
		TBranch        *b_neu_l0_hrc_dihadron_lowmult_tis;   //!
		TBranch        *b_neu_l0_hrc_dihadron_lowmult_tos;   //!
		TBranch        *b_neu_hlt1_passthru_dec;   //!
		TBranch        *b_neu_hlt1_passthru_tis;   //!
		TBranch        *b_neu_hlt1_passthru_tos;   //!
		TBranch        *b_neu_hlt1_lowmult_dec;   //!
		TBranch        *b_neu_hlt1_lowmult_tis;   //!
		TBranch        *b_neu_hlt1_lowmult_tos;   //!
		TBranch        *b_neu_hlt1_herchel_dec;   //!
		TBranch        *b_neu_hlt1_herchel_tis;   //!
		TBranch        *b_neu_hlt1_herchel_tos;   //!
		TBranch        *b_neu_hlt1_velo_dec;   //!
		TBranch        *b_neu_hlt1_velo_tis;   //!
		TBranch        *b_neu_hlt1_velo_tos;   //!
		TBranch        *b_neu_hlt1_veloherchel_dec;   //!
		TBranch        *b_neu_hlt1_veloherchel_tis;   //!
		TBranch        *b_neu_hlt1_veloherchel_tos;   //!
		TBranch        *b_neu_hlt2_lmr2hh_dec;   //!
		TBranch        *b_neu_hlt2_lmr2hh_tis;   //!
		TBranch        *b_neu_hlt2_lmr2hh_tos;   //!
		TBranch        *b_neu_hlt2_lmr2hhws_dec;   //!
		TBranch        *b_neu_hlt2_lmr2hhws_tis;   //!
		TBranch        *b_neu_hlt2_lmr2hhws_tos;   //!
		TBranch        *b_phi_idx_pvr;   //!
		TBranch        *b_phi_p;   //!
		TBranch        *b_phi_pt;   //!
		TBranch        *b_phi_px;   //!
		TBranch        *b_phi_py;   //!
		TBranch        *b_phi_pz;   //!
		TBranch        *b_phi_e;   //!
		TBranch        *b_phi_x;   //!
		TBranch        *b_phi_y;   //!
		TBranch        *b_phi_z;   //!
		TBranch        *b_phi_m;   //!
		TBranch        *b_phi_ip;   //!
		TBranch        *b_phi_ip_chi2;   //!
		TBranch        *b_phi_vtx_chi2;   //!
		TBranch        *b_phi_vtx_ndof;   //!
		TBranch        *b_phi_fd;   //!
		TBranch        *b_phi_fd_chi2;   //!
		TBranch        *b_phi_tau;   //!
		TBranch        *b_phi_tau_err;   //!
		TBranch        *b_phi_tau_chi2;   //!
		TBranch        *b_phi_idx_trk0;   //!
		TBranch        *b_phi_idx_trk1;   //!
		TBranch        *b_phi_l0_hadron_dec;   //!
		TBranch        *b_phi_l0_hadron_tis;   //!
		TBranch        *b_phi_l0_hadron_tos;   //!
		TBranch        *b_phi_l0_dihadron_lowmult_dec;   //!
		TBranch        *b_phi_l0_dihadron_lowmult_tis;   //!
		TBranch        *b_phi_l0_dihadron_lowmult_tos;   //!
		TBranch        *b_phi_l0_hrc_dihadron_lowmult_dec;   //!
		TBranch        *b_phi_l0_hrc_dihadron_lowmult_tis;   //!
		TBranch        *b_phi_l0_hrc_dihadron_lowmult_tos;   //!
		TBranch        *b_phi_hlt1_passthru_dec;   //!
		TBranch        *b_phi_hlt1_passthru_tis;   //!
		TBranch        *b_phi_hlt1_passthru_tos;   //!
		TBranch        *b_phi_hlt1_lowmult_dec;   //!
		TBranch        *b_phi_hlt1_lowmult_tis;   //!
		TBranch        *b_phi_hlt1_lowmult_tos;   //!
		TBranch        *b_phi_hlt1_herchel_dec;   //!
		TBranch        *b_phi_hlt1_herchel_tis;   //!
		TBranch        *b_phi_hlt1_herchel_tos;   //!
		TBranch        *b_phi_hlt1_velo_dec;   //!
		TBranch        *b_phi_hlt1_velo_tis;   //!
		TBranch        *b_phi_hlt1_velo_tos;   //!
		TBranch        *b_phi_hlt1_veloherchel_dec;   //!
		TBranch        *b_phi_hlt1_veloherchel_tis;   //!
		TBranch        *b_phi_hlt1_veloherchel_tos;   //!
		TBranch        *b_phi_hlt2_lmr2hh_dec;   //!
		TBranch        *b_phi_hlt2_lmr2hh_tis;   //!
		TBranch        *b_phi_hlt2_lmr2hh_tos;   //!
		TBranch        *b_phi_hlt2_lmr2hhws_dec;   //!
		TBranch        *b_phi_hlt2_lmr2hhws_tis;   //!
		TBranch        *b_phi_hlt2_lmr2hhws_tos;   //!
		TBranch        *b_evt_pvr_n;   //!
		TBranch        *b_evt_neu_n;   //!
		TBranch        *b_evt_chg_n;   //!
		TBranch        *b_evt_trk_n;   //!
		TBranch        *b_evt_trk_n_velor;   //!
		TBranch        *b_evt_trk_n_velo;   //!
		TBranch        *b_evt_trk_n_long;   //!
		TBranch        *b_evt_trk_n_up;   //!
		TBranch        *b_evt_trk_n_down;   //!
		TBranch        *b_evt_trk_n_t;   //!
		TBranch        *b_evt_trk_n_mu;   //!
		TBranch        *b_evt_hrc_fom;   //!

		skimMC(int job, int sjob, TString dir);
		virtual ~skimMC();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		void fillOutput(int gp, int gkp, int gkm, int rp, int rkp, int rkm);

		TString outName;
		TFile* fout;
		TTree* tout;
};

#endif

#ifdef skimMC_cxx
skimMC::skimMC(int job, int sjob, TString dir) : fChain(0) 
{
	outName=dir;
	outName+="/";
	outName+=job;
	outName+="/";
	outName+=sjob;
	outName+="/skimmed.root";

	fout = TFile::Open(outName,"RECREATE");
	tout = new TTree("T","T");

	char str[256];
	TChain* t(0);
	t = new TChain("data");

	for(int i=0; i<3000; ++i) {
		if(sjob<0 || sjob==i) {
			sprintf(str,"%s/%d/%d/output.root",dir.Data(),job,i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}

	Init(t);
}

skimMC::~skimMC()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t skimMC::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t skimMC::LoadTree(Long64_t entry)
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

void skimMC::Init(TTree *tree)
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
	trk_idx_gen = 0;
	trk_idx_pvr = 0;
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
	trk_l0_hadron_dec = 0;
	trk_l0_hadron_tis = 0;
	trk_l0_hadron_tos = 0;
	trk_l0_dihadron_lowmult_dec = 0;
	trk_l0_dihadron_lowmult_tis = 0;
	trk_l0_dihadron_lowmult_tos = 0;
	trk_l0_hrc_dihadron_lowmult_dec = 0;
	trk_l0_hrc_dihadron_lowmult_tis = 0;
	trk_l0_hrc_dihadron_lowmult_tos = 0;
	trk_hlt1_passthru_dec = 0;
	trk_hlt1_passthru_tis = 0;
	trk_hlt1_passthru_tos = 0;
	trk_hlt1_lowmult_dec = 0;
	trk_hlt1_lowmult_tis = 0;
	trk_hlt1_lowmult_tos = 0;
	trk_hlt1_herchel_dec = 0;
	trk_hlt1_herchel_tis = 0;
	trk_hlt1_herchel_tos = 0;
	trk_hlt1_velo_dec = 0;
	trk_hlt1_velo_tis = 0;
	trk_hlt1_velo_tos = 0;
	trk_hlt1_veloherchel_dec = 0;
	trk_hlt1_veloherchel_tis = 0;
	trk_hlt1_veloherchel_tos = 0;
	trk_hlt2_lmr2hh_dec = 0;
	trk_hlt2_lmr2hh_tis = 0;
	trk_hlt2_lmr2hh_tos = 0;
	trk_hlt2_lmr2hhws_dec = 0;
	trk_hlt2_lmr2hhws_tis = 0;
	trk_hlt2_lmr2hhws_tos = 0;
	neu_idx_gen = 0;
	neu_p = 0;
	neu_pt = 0;
	neu_px = 0;
	neu_py = 0;
	neu_pz = 0;
	neu_e = 0;
	neu_pid = 0;
	neu_l0_hadron_dec = 0;
	neu_l0_hadron_tis = 0;
	neu_l0_hadron_tos = 0;
	neu_l0_dihadron_lowmult_dec = 0;
	neu_l0_dihadron_lowmult_tis = 0;
	neu_l0_dihadron_lowmult_tos = 0;
	neu_l0_hrc_dihadron_lowmult_dec = 0;
	neu_l0_hrc_dihadron_lowmult_tis = 0;
	neu_l0_hrc_dihadron_lowmult_tos = 0;
	neu_hlt1_passthru_dec = 0;
	neu_hlt1_passthru_tis = 0;
	neu_hlt1_passthru_tos = 0;
	neu_hlt1_lowmult_dec = 0;
	neu_hlt1_lowmult_tis = 0;
	neu_hlt1_lowmult_tos = 0;
	neu_hlt1_herchel_dec = 0;
	neu_hlt1_herchel_tis = 0;
	neu_hlt1_herchel_tos = 0;
	neu_hlt1_velo_dec = 0;
	neu_hlt1_velo_tis = 0;
	neu_hlt1_velo_tos = 0;
	neu_hlt1_veloherchel_dec = 0;
	neu_hlt1_veloherchel_tis = 0;
	neu_hlt1_veloherchel_tos = 0;
	neu_hlt2_lmr2hh_dec = 0;
	neu_hlt2_lmr2hh_tis = 0;
	neu_hlt2_lmr2hh_tos = 0;
	neu_hlt2_lmr2hhws_dec = 0;
	neu_hlt2_lmr2hhws_tis = 0;
	neu_hlt2_lmr2hhws_tos = 0;
	phi_idx_pvr = 0;
	phi_p = 0;
	phi_pt = 0;
	phi_px = 0;
	phi_py = 0;
	phi_pz = 0;
	phi_e = 0;
	phi_x = 0;
	phi_y = 0;
	phi_z = 0;
	phi_m = 0;
	phi_ip = 0;
	phi_ip_chi2 = 0;
	phi_vtx_chi2 = 0;
	phi_vtx_ndof = 0;
	phi_fd = 0;
	phi_fd_chi2 = 0;
	phi_tau = 0;
	phi_tau_err = 0;
	phi_tau_chi2 = 0;
	phi_idx_trk0 = 0;
	phi_idx_trk1 = 0;
	phi_l0_hadron_dec = 0;
	phi_l0_hadron_tis = 0;
	phi_l0_hadron_tos = 0;
	phi_l0_dihadron_lowmult_dec = 0;
	phi_l0_dihadron_lowmult_tis = 0;
	phi_l0_dihadron_lowmult_tos = 0;
	phi_l0_hrc_dihadron_lowmult_dec = 0;
	phi_l0_hrc_dihadron_lowmult_tis = 0;
	phi_l0_hrc_dihadron_lowmult_tos = 0;
	phi_hlt1_passthru_dec = 0;
	phi_hlt1_passthru_tis = 0;
	phi_hlt1_passthru_tos = 0;
	phi_hlt1_lowmult_dec = 0;
	phi_hlt1_lowmult_tis = 0;
	phi_hlt1_lowmult_tos = 0;
	phi_hlt1_herchel_dec = 0;
	phi_hlt1_herchel_tis = 0;
	phi_hlt1_herchel_tos = 0;
	phi_hlt1_velo_dec = 0;
	phi_hlt1_velo_tis = 0;
	phi_hlt1_velo_tos = 0;
	phi_hlt1_veloherchel_dec = 0;
	phi_hlt1_veloherchel_tis = 0;
	phi_hlt1_veloherchel_tos = 0;
	phi_hlt2_lmr2hh_dec = 0;
	phi_hlt2_lmr2hh_tis = 0;
	phi_hlt2_lmr2hh_tos = 0;
	phi_hlt2_lmr2hhws_dec = 0;
	phi_hlt2_lmr2hhws_tis = 0;
	phi_hlt2_lmr2hhws_tos = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("gen_idx_pvr", &gen_idx_pvr, &b_gen_idx_pvr);
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
	fChain->SetBranchAddress("trk_idx_gen", &trk_idx_gen, &b_trk_idx_gen);
	fChain->SetBranchAddress("trk_idx_pvr", &trk_idx_pvr, &b_trk_idx_pvr);
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
	fChain->SetBranchAddress("trk_l0_hadron_dec", &trk_l0_hadron_dec, &b_trk_l0_hadron_dec);
	fChain->SetBranchAddress("trk_l0_hadron_tis", &trk_l0_hadron_tis, &b_trk_l0_hadron_tis);
	fChain->SetBranchAddress("trk_l0_hadron_tos", &trk_l0_hadron_tos, &b_trk_l0_hadron_tos);
	fChain->SetBranchAddress("trk_l0_dihadron_lowmult_dec", &trk_l0_dihadron_lowmult_dec, &b_trk_l0_dihadron_lowmult_dec);
	fChain->SetBranchAddress("trk_l0_dihadron_lowmult_tis", &trk_l0_dihadron_lowmult_tis, &b_trk_l0_dihadron_lowmult_tis);
	fChain->SetBranchAddress("trk_l0_dihadron_lowmult_tos", &trk_l0_dihadron_lowmult_tos, &b_trk_l0_dihadron_lowmult_tos);
	fChain->SetBranchAddress("trk_l0_hrc_dihadron_lowmult_dec", &trk_l0_hrc_dihadron_lowmult_dec, &b_trk_l0_hrc_dihadron_lowmult_dec);
	fChain->SetBranchAddress("trk_l0_hrc_dihadron_lowmult_tis", &trk_l0_hrc_dihadron_lowmult_tis, &b_trk_l0_hrc_dihadron_lowmult_tis);
	fChain->SetBranchAddress("trk_l0_hrc_dihadron_lowmult_tos", &trk_l0_hrc_dihadron_lowmult_tos, &b_trk_l0_hrc_dihadron_lowmult_tos);
	fChain->SetBranchAddress("trk_hlt1_passthru_dec", &trk_hlt1_passthru_dec, &b_trk_hlt1_passthru_dec);
	fChain->SetBranchAddress("trk_hlt1_passthru_tis", &trk_hlt1_passthru_tis, &b_trk_hlt1_passthru_tis);
	fChain->SetBranchAddress("trk_hlt1_passthru_tos", &trk_hlt1_passthru_tos, &b_trk_hlt1_passthru_tos);
	fChain->SetBranchAddress("trk_hlt1_lowmult_dec", &trk_hlt1_lowmult_dec, &b_trk_hlt1_lowmult_dec);
	fChain->SetBranchAddress("trk_hlt1_lowmult_tis", &trk_hlt1_lowmult_tis, &b_trk_hlt1_lowmult_tis);
	fChain->SetBranchAddress("trk_hlt1_lowmult_tos", &trk_hlt1_lowmult_tos, &b_trk_hlt1_lowmult_tos);
	fChain->SetBranchAddress("trk_hlt1_herchel_dec", &trk_hlt1_herchel_dec, &b_trk_hlt1_herchel_dec);
	fChain->SetBranchAddress("trk_hlt1_herchel_tis", &trk_hlt1_herchel_tis, &b_trk_hlt1_herchel_tis);
	fChain->SetBranchAddress("trk_hlt1_herchel_tos", &trk_hlt1_herchel_tos, &b_trk_hlt1_herchel_tos);
	fChain->SetBranchAddress("trk_hlt1_velo_dec", &trk_hlt1_velo_dec, &b_trk_hlt1_velo_dec);
	fChain->SetBranchAddress("trk_hlt1_velo_tis", &trk_hlt1_velo_tis, &b_trk_hlt1_velo_tis);
	fChain->SetBranchAddress("trk_hlt1_velo_tos", &trk_hlt1_velo_tos, &b_trk_hlt1_velo_tos);
	fChain->SetBranchAddress("trk_hlt1_veloherchel_dec", &trk_hlt1_veloherchel_dec, &b_trk_hlt1_veloherchel_dec);
	fChain->SetBranchAddress("trk_hlt1_veloherchel_tis", &trk_hlt1_veloherchel_tis, &b_trk_hlt1_veloherchel_tis);
	fChain->SetBranchAddress("trk_hlt1_veloherchel_tos", &trk_hlt1_veloherchel_tos, &b_trk_hlt1_veloherchel_tos);
	fChain->SetBranchAddress("trk_hlt2_lmr2hh_dec", &trk_hlt2_lmr2hh_dec, &b_trk_hlt2_lmr2hh_dec);
	fChain->SetBranchAddress("trk_hlt2_lmr2hh_tis", &trk_hlt2_lmr2hh_tis, &b_trk_hlt2_lmr2hh_tis);
	fChain->SetBranchAddress("trk_hlt2_lmr2hh_tos", &trk_hlt2_lmr2hh_tos, &b_trk_hlt2_lmr2hh_tos);
	fChain->SetBranchAddress("trk_hlt2_lmr2hhws_dec", &trk_hlt2_lmr2hhws_dec, &b_trk_hlt2_lmr2hhws_dec);
	fChain->SetBranchAddress("trk_hlt2_lmr2hhws_tis", &trk_hlt2_lmr2hhws_tis, &b_trk_hlt2_lmr2hhws_tis);
	fChain->SetBranchAddress("trk_hlt2_lmr2hhws_tos", &trk_hlt2_lmr2hhws_tos, &b_trk_hlt2_lmr2hhws_tos);
	fChain->SetBranchAddress("neu_idx_gen", &neu_idx_gen, &b_neu_idx_gen);
	fChain->SetBranchAddress("neu_p", &neu_p, &b_neu_p);
	fChain->SetBranchAddress("neu_pt", &neu_pt, &b_neu_pt);
	fChain->SetBranchAddress("neu_px", &neu_px, &b_neu_px);
	fChain->SetBranchAddress("neu_py", &neu_py, &b_neu_py);
	fChain->SetBranchAddress("neu_pz", &neu_pz, &b_neu_pz);
	fChain->SetBranchAddress("neu_e", &neu_e, &b_neu_e);
	fChain->SetBranchAddress("neu_pid", &neu_pid, &b_neu_pid);
	fChain->SetBranchAddress("neu_l0_hadron_dec", &neu_l0_hadron_dec, &b_neu_l0_hadron_dec);
	fChain->SetBranchAddress("neu_l0_hadron_tis", &neu_l0_hadron_tis, &b_neu_l0_hadron_tis);
	fChain->SetBranchAddress("neu_l0_hadron_tos", &neu_l0_hadron_tos, &b_neu_l0_hadron_tos);
	fChain->SetBranchAddress("neu_l0_dihadron_lowmult_dec", &neu_l0_dihadron_lowmult_dec, &b_neu_l0_dihadron_lowmult_dec);
	fChain->SetBranchAddress("neu_l0_dihadron_lowmult_tis", &neu_l0_dihadron_lowmult_tis, &b_neu_l0_dihadron_lowmult_tis);
	fChain->SetBranchAddress("neu_l0_dihadron_lowmult_tos", &neu_l0_dihadron_lowmult_tos, &b_neu_l0_dihadron_lowmult_tos);
	fChain->SetBranchAddress("neu_l0_hrc_dihadron_lowmult_dec", &neu_l0_hrc_dihadron_lowmult_dec, &b_neu_l0_hrc_dihadron_lowmult_dec);
	fChain->SetBranchAddress("neu_l0_hrc_dihadron_lowmult_tis", &neu_l0_hrc_dihadron_lowmult_tis, &b_neu_l0_hrc_dihadron_lowmult_tis);
	fChain->SetBranchAddress("neu_l0_hrc_dihadron_lowmult_tos", &neu_l0_hrc_dihadron_lowmult_tos, &b_neu_l0_hrc_dihadron_lowmult_tos);
	fChain->SetBranchAddress("neu_hlt1_passthru_dec", &neu_hlt1_passthru_dec, &b_neu_hlt1_passthru_dec);
	fChain->SetBranchAddress("neu_hlt1_passthru_tis", &neu_hlt1_passthru_tis, &b_neu_hlt1_passthru_tis);
	fChain->SetBranchAddress("neu_hlt1_passthru_tos", &neu_hlt1_passthru_tos, &b_neu_hlt1_passthru_tos);
	fChain->SetBranchAddress("neu_hlt1_lowmult_dec", &neu_hlt1_lowmult_dec, &b_neu_hlt1_lowmult_dec);
	fChain->SetBranchAddress("neu_hlt1_lowmult_tis", &neu_hlt1_lowmult_tis, &b_neu_hlt1_lowmult_tis);
	fChain->SetBranchAddress("neu_hlt1_lowmult_tos", &neu_hlt1_lowmult_tos, &b_neu_hlt1_lowmult_tos);
	fChain->SetBranchAddress("neu_hlt1_herchel_dec", &neu_hlt1_herchel_dec, &b_neu_hlt1_herchel_dec);
	fChain->SetBranchAddress("neu_hlt1_herchel_tis", &neu_hlt1_herchel_tis, &b_neu_hlt1_herchel_tis);
	fChain->SetBranchAddress("neu_hlt1_herchel_tos", &neu_hlt1_herchel_tos, &b_neu_hlt1_herchel_tos);
	fChain->SetBranchAddress("neu_hlt1_velo_dec", &neu_hlt1_velo_dec, &b_neu_hlt1_velo_dec);
	fChain->SetBranchAddress("neu_hlt1_velo_tis", &neu_hlt1_velo_tis, &b_neu_hlt1_velo_tis);
	fChain->SetBranchAddress("neu_hlt1_velo_tos", &neu_hlt1_velo_tos, &b_neu_hlt1_velo_tos);
	fChain->SetBranchAddress("neu_hlt1_veloherchel_dec", &neu_hlt1_veloherchel_dec, &b_neu_hlt1_veloherchel_dec);
	fChain->SetBranchAddress("neu_hlt1_veloherchel_tis", &neu_hlt1_veloherchel_tis, &b_neu_hlt1_veloherchel_tis);
	fChain->SetBranchAddress("neu_hlt1_veloherchel_tos", &neu_hlt1_veloherchel_tos, &b_neu_hlt1_veloherchel_tos);
	fChain->SetBranchAddress("neu_hlt2_lmr2hh_dec", &neu_hlt2_lmr2hh_dec, &b_neu_hlt2_lmr2hh_dec);
	fChain->SetBranchAddress("neu_hlt2_lmr2hh_tis", &neu_hlt2_lmr2hh_tis, &b_neu_hlt2_lmr2hh_tis);
	fChain->SetBranchAddress("neu_hlt2_lmr2hh_tos", &neu_hlt2_lmr2hh_tos, &b_neu_hlt2_lmr2hh_tos);
	fChain->SetBranchAddress("neu_hlt2_lmr2hhws_dec", &neu_hlt2_lmr2hhws_dec, &b_neu_hlt2_lmr2hhws_dec);
	fChain->SetBranchAddress("neu_hlt2_lmr2hhws_tis", &neu_hlt2_lmr2hhws_tis, &b_neu_hlt2_lmr2hhws_tis);
	fChain->SetBranchAddress("neu_hlt2_lmr2hhws_tos", &neu_hlt2_lmr2hhws_tos, &b_neu_hlt2_lmr2hhws_tos);
	fChain->SetBranchAddress("phi_idx_pvr", &phi_idx_pvr, &b_phi_idx_pvr);
	fChain->SetBranchAddress("phi_p", &phi_p, &b_phi_p);
	fChain->SetBranchAddress("phi_pt", &phi_pt, &b_phi_pt);
	fChain->SetBranchAddress("phi_px", &phi_px, &b_phi_px);
	fChain->SetBranchAddress("phi_py", &phi_py, &b_phi_py);
	fChain->SetBranchAddress("phi_pz", &phi_pz, &b_phi_pz);
	fChain->SetBranchAddress("phi_e", &phi_e, &b_phi_e);
	fChain->SetBranchAddress("phi_x", &phi_x, &b_phi_x);
	fChain->SetBranchAddress("phi_y", &phi_y, &b_phi_y);
	fChain->SetBranchAddress("phi_z", &phi_z, &b_phi_z);
	fChain->SetBranchAddress("phi_m", &phi_m, &b_phi_m);
	fChain->SetBranchAddress("phi_ip", &phi_ip, &b_phi_ip);
	fChain->SetBranchAddress("phi_ip_chi2", &phi_ip_chi2, &b_phi_ip_chi2);
	fChain->SetBranchAddress("phi_vtx_chi2", &phi_vtx_chi2, &b_phi_vtx_chi2);
	fChain->SetBranchAddress("phi_vtx_ndof", &phi_vtx_ndof, &b_phi_vtx_ndof);
	fChain->SetBranchAddress("phi_fd", &phi_fd, &b_phi_fd);
	fChain->SetBranchAddress("phi_fd_chi2", &phi_fd_chi2, &b_phi_fd_chi2);
	fChain->SetBranchAddress("phi_tau", &phi_tau, &b_phi_tau);
	fChain->SetBranchAddress("phi_tau_err", &phi_tau_err, &b_phi_tau_err);
	fChain->SetBranchAddress("phi_tau_chi2", &phi_tau_chi2, &b_phi_tau_chi2);
	fChain->SetBranchAddress("phi_idx_trk0", &phi_idx_trk0, &b_phi_idx_trk0);
	fChain->SetBranchAddress("phi_idx_trk1", &phi_idx_trk1, &b_phi_idx_trk1);
	fChain->SetBranchAddress("phi_l0_hadron_dec", &phi_l0_hadron_dec, &b_phi_l0_hadron_dec);
	fChain->SetBranchAddress("phi_l0_hadron_tis", &phi_l0_hadron_tis, &b_phi_l0_hadron_tis);
	fChain->SetBranchAddress("phi_l0_hadron_tos", &phi_l0_hadron_tos, &b_phi_l0_hadron_tos);
	fChain->SetBranchAddress("phi_l0_dihadron_lowmult_dec", &phi_l0_dihadron_lowmult_dec, &b_phi_l0_dihadron_lowmult_dec);
	fChain->SetBranchAddress("phi_l0_dihadron_lowmult_tis", &phi_l0_dihadron_lowmult_tis, &b_phi_l0_dihadron_lowmult_tis);
	fChain->SetBranchAddress("phi_l0_dihadron_lowmult_tos", &phi_l0_dihadron_lowmult_tos, &b_phi_l0_dihadron_lowmult_tos);
	fChain->SetBranchAddress("phi_l0_hrc_dihadron_lowmult_dec", &phi_l0_hrc_dihadron_lowmult_dec, &b_phi_l0_hrc_dihadron_lowmult_dec);
	fChain->SetBranchAddress("phi_l0_hrc_dihadron_lowmult_tis", &phi_l0_hrc_dihadron_lowmult_tis, &b_phi_l0_hrc_dihadron_lowmult_tis);
	fChain->SetBranchAddress("phi_l0_hrc_dihadron_lowmult_tos", &phi_l0_hrc_dihadron_lowmult_tos, &b_phi_l0_hrc_dihadron_lowmult_tos);
	fChain->SetBranchAddress("phi_hlt1_passthru_dec", &phi_hlt1_passthru_dec, &b_phi_hlt1_passthru_dec);
	fChain->SetBranchAddress("phi_hlt1_passthru_tis", &phi_hlt1_passthru_tis, &b_phi_hlt1_passthru_tis);
	fChain->SetBranchAddress("phi_hlt1_passthru_tos", &phi_hlt1_passthru_tos, &b_phi_hlt1_passthru_tos);
	fChain->SetBranchAddress("phi_hlt1_lowmult_dec", &phi_hlt1_lowmult_dec, &b_phi_hlt1_lowmult_dec);
	fChain->SetBranchAddress("phi_hlt1_lowmult_tis", &phi_hlt1_lowmult_tis, &b_phi_hlt1_lowmult_tis);
	fChain->SetBranchAddress("phi_hlt1_lowmult_tos", &phi_hlt1_lowmult_tos, &b_phi_hlt1_lowmult_tos);
	fChain->SetBranchAddress("phi_hlt1_herchel_dec", &phi_hlt1_herchel_dec, &b_phi_hlt1_herchel_dec);
	fChain->SetBranchAddress("phi_hlt1_herchel_tis", &phi_hlt1_herchel_tis, &b_phi_hlt1_herchel_tis);
	fChain->SetBranchAddress("phi_hlt1_herchel_tos", &phi_hlt1_herchel_tos, &b_phi_hlt1_herchel_tos);
	fChain->SetBranchAddress("phi_hlt1_velo_dec", &phi_hlt1_velo_dec, &b_phi_hlt1_velo_dec);
	fChain->SetBranchAddress("phi_hlt1_velo_tis", &phi_hlt1_velo_tis, &b_phi_hlt1_velo_tis);
	fChain->SetBranchAddress("phi_hlt1_velo_tos", &phi_hlt1_velo_tos, &b_phi_hlt1_velo_tos);
	fChain->SetBranchAddress("phi_hlt1_veloherchel_dec", &phi_hlt1_veloherchel_dec, &b_phi_hlt1_veloherchel_dec);
	fChain->SetBranchAddress("phi_hlt1_veloherchel_tis", &phi_hlt1_veloherchel_tis, &b_phi_hlt1_veloherchel_tis);
	fChain->SetBranchAddress("phi_hlt1_veloherchel_tos", &phi_hlt1_veloherchel_tos, &b_phi_hlt1_veloherchel_tos);
	fChain->SetBranchAddress("phi_hlt2_lmr2hh_dec", &phi_hlt2_lmr2hh_dec, &b_phi_hlt2_lmr2hh_dec);
	fChain->SetBranchAddress("phi_hlt2_lmr2hh_tis", &phi_hlt2_lmr2hh_tis, &b_phi_hlt2_lmr2hh_tis);
	fChain->SetBranchAddress("phi_hlt2_lmr2hh_tos", &phi_hlt2_lmr2hh_tos, &b_phi_hlt2_lmr2hh_tos);
	fChain->SetBranchAddress("phi_hlt2_lmr2hhws_dec", &phi_hlt2_lmr2hhws_dec, &b_phi_hlt2_lmr2hhws_dec);
	fChain->SetBranchAddress("phi_hlt2_lmr2hhws_tis", &phi_hlt2_lmr2hhws_tis, &b_phi_hlt2_lmr2hhws_tis);
	fChain->SetBranchAddress("phi_hlt2_lmr2hhws_tos", &phi_hlt2_lmr2hhws_tos, &b_phi_hlt2_lmr2hhws_tos);
	fChain->SetBranchAddress("evt_pvr_n", &evt_pvr_n, &b_evt_pvr_n);
	fChain->SetBranchAddress("evt_neu_n", &evt_neu_n, &b_evt_neu_n);
	fChain->SetBranchAddress("evt_chg_n", &evt_chg_n, &b_evt_chg_n);
	fChain->SetBranchAddress("evt_trk_n", &evt_trk_n, &b_evt_trk_n);
	fChain->SetBranchAddress("evt_trk_n_velor", &evt_trk_n_velor, &b_evt_trk_n_velor);
	fChain->SetBranchAddress("evt_trk_n_velo", &evt_trk_n_velo, &b_evt_trk_n_velo);
	fChain->SetBranchAddress("evt_trk_n_long", &evt_trk_n_long, &b_evt_trk_n_long);
	fChain->SetBranchAddress("evt_trk_n_up", &evt_trk_n_up, &b_evt_trk_n_up);
	fChain->SetBranchAddress("evt_trk_n_down", &evt_trk_n_down, &b_evt_trk_n_down);
	fChain->SetBranchAddress("evt_trk_n_t", &evt_trk_n_t, &b_evt_trk_n_t);
	fChain->SetBranchAddress("evt_trk_n_mu", &evt_trk_n_mu, &b_evt_trk_n_mu);
	fChain->SetBranchAddress("evt_hrc_fom", &evt_hrc_fom, &b_evt_hrc_fom);

	tout->Branch("Kplus_TRUEP", 	&Kplus_TRUEP);
	tout->Branch("Kplus_TRUEPT", 	&Kplus_TRUEPT);
	tout->Branch("Kplus_TRUEPX", 	&Kplus_TRUEPX);
	tout->Branch("Kplus_TRUEPY", 	&Kplus_TRUEPY);
	tout->Branch("Kplus_TRUEPZ", 	&Kplus_TRUEPZ);
	tout->Branch("Kplus_TRUEE",  	&Kplus_TRUEE);
	tout->Branch("Kminus_TRUEP",	&Kminus_TRUEP);
	tout->Branch("Kminus_TRUEPT",	&Kminus_TRUEPT);
	tout->Branch("Kminus_TRUEPX",	&Kminus_TRUEPX);
	tout->Branch("Kminus_TRUEPY",	&Kminus_TRUEPY);
	tout->Branch("Kminus_TRUEPZ",	&Kminus_TRUEPZ);
	tout->Branch("Kminus_TRUEE", 	&Kminus_TRUEE);
	tout->Branch("phi_TRUEP",   	&phi_TRUEP);
	tout->Branch("phi_TRUEPT",   	&phi_TRUEPT);
	tout->Branch("phi_TRUEPX",   	&phi_TRUEPX);
	tout->Branch("phi_TRUEPY",   	&phi_TRUEPY);
	tout->Branch("phi_TRUEPZ",   	&phi_TRUEPZ);
	tout->Branch("phi_TRUEE",    	&phi_TRUEE);
	tout->Branch("Kplus_P",     	&Kplus_P);
	tout->Branch("Kplus_PT",     	&Kplus_PT);
	tout->Branch("Kplus_PX",     	&Kplus_PX);
	tout->Branch("Kplus_PY",     	&Kplus_PY);
	tout->Branch("Kplus_PZ",     	&Kplus_PZ);
	tout->Branch("Kplus_E",      	&Kplus_E);
	tout->Branch("Kminus_P",    	&Kminus_P);
	tout->Branch("Kminus_PT",    	&Kminus_PT);
	tout->Branch("Kminus_PX",    	&Kminus_PX);
	tout->Branch("Kminus_PY",    	&Kminus_PY);
	tout->Branch("Kminus_PZ",    	&Kminus_PZ);
	tout->Branch("Kminus_E",     	&Kminus_E);
	tout->Branch("phi_P",       	&phi_P);
	tout->Branch("phi_PT",       	&phi_PT);
	tout->Branch("phi_PX",       	&phi_PX);
	tout->Branch("phi_PY",       	&phi_PY);
	tout->Branch("phi_PZ",       	&phi_PZ);
	tout->Branch("phi_E",        	&phi_E);
	tout->Branch("phi_M",        	&phi_M);

	Notify();
}

Bool_t skimMC::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void skimMC::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t skimMC::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef skimMC_cxx
