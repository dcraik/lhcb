//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  1 12:56:17 2020 by ROOT version 6.18/00
// from TTree data/data
// found on file: /data/alps/4/23/output.root
//////////////////////////////////////////////////////////

#ifndef TreeSkimmer_h
#define TreeSkimmer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

//structures to simplify filling skimmed tree for a single B mode
struct Track {
	double p;
	double pt;
	double px;
	double py;
	double pz;
	double e;
	double pid;
	double q;
	double ip;
	double ip_best;
	double ip_bvtx;
	double ip_chi2;
	double ip_chi2_best;
	double ip_chi2_bvtx;
	double pnn_e;
	double pnn_mu;
	double pnn_pi;
	double pnn_k;
	double pnn_p;
	double pnn_ghost;
	double ecal;
	double hcal;
	double prb_ghost;
	double type;
	double is_mu;
	double vid;
	double x;
	double y;
	double z;
	double l0_hadron_dec;
	double l0_hadron_tis;
	double l0_hadron_tos;
	double l0_photon_dec;
	double l0_photon_tis;
	double l0_photon_tos;
	double l0_electron_dec;
	double l0_electron_tis;
	double l0_electron_tos;
	double l0_muon_dec;
	double l0_muon_tis;
	double l0_muon_tos;
	double l0_dimuon_dec;
	double l0_dimuon_tis;
	double l0_dimuon_tos;
};
struct Neutral {
	double p;
	double pt;
	double px;
	double py;
	double pz;
	double e;
	double pid;
	double cl;
	double l0_hadron_dec;
	double l0_hadron_tis;
	double l0_hadron_tos;
	double l0_photon_dec;
	double l0_photon_tis;
	double l0_photon_tos;
	double l0_electron_dec;
	double l0_electron_tis;
	double l0_electron_tos;
	double l0_muon_dec;
	double l0_muon_tis;
	double l0_muon_tos;
	double l0_dimuon_dec;
	double l0_dimuon_tis;
	double l0_dimuon_tos;
	double m;
};
struct BPart {
	double p;
	double pt;
	double px;
	double py;
	double pz;
	double e;
	double x;
	double y;
	double z;
	double m;
	double vtx_chi2;
	double vtx_ndof;
	double ip;
	double ip_chi2;
	double fd;
	double fd_chi2;
	double tau;
	double tau_err;
	double tau_chi2;
	double l0_hadron_dec;
	double l0_hadron_tis;
	double l0_hadron_tos;
	double l0_photon_dec;
	double l0_photon_tis;
	double l0_photon_tos;
	double l0_electron_dec;
	double l0_electron_tis;
	double l0_electron_tos;
	double l0_muon_dec;
	double l0_muon_tis;
	double l0_muon_tos;
	double l0_dimuon_dec;
	double l0_dimuon_tis;
	double l0_dimuon_tos;
	double dira;
	double mcor;
	double maxdoca;
	double doca_kpi;
	double m_kpi;
};
struct XPart {
	double p;
	double pt;
	double px;
	double py;
	double pz;
	double e;
	double x;
	double y;
	double z;
	double m;
	double vtx_chi2;
	double vtx_ndof;
	double ip;
	double ip_best;
	double ip_bvtx;
	double ip_chi2;
	double ip_chi2_best;
	double ip_chi2_bvtx;
	double fd;
	double fd_best;
	double fd_bvtx;
	double fd_chi2;
	double fd_chi2_best;
	double fd_chi2_bvtx;
	double tau;
	double tau_err;
	double tau_chi2;
	double l0_hadron_dec;
	double l0_hadron_tis;
	double l0_hadron_tos;
	double l0_photon_dec;
	double l0_photon_tis;
	double l0_photon_tos;
	double l0_electron_dec;
	double l0_electron_tis;
	double l0_electron_tos;
	double l0_muon_dec;
	double l0_muon_tis;
	double l0_muon_tos;
	double l0_dimuon_dec;
	double l0_dimuon_tis;
	double l0_dimuon_tos;
	double maxdoca;
};
struct DTFVars {
	std::vector<double> px;
	std::vector<double> py;
	std::vector<double> pz;
	std::vector<double> pe;
	std::vector<double> pm;
	std::vector<double> ctau;
	double chi2;
	double ndof;
	double _ntrk;
	double _nneu;

	DTFVars(int ntrk=8, int nneu=8)
		: px(2+ntrk+nneu,0.), py(2+ntrk+nneu,0.), pz(2+ntrk+nneu,0.),
		  pe(2+ntrk+nneu,0.), pm(2+ntrk+nneu,0.), ctau(2,0.), _ntrk(ntrk), _nneu(nneu) {}
};

class TreeSkimmer {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
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
		std::vector<double>  *trk_ip_best;
		std::vector<double>  *trk_ip_bvtx;
		std::vector<double>  *trk_ip_chi2;
		std::vector<double>  *trk_ip_chi2_best;
		std::vector<double>  *trk_ip_chi2_bvtx;
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
		std::vector<double>  *trk_l0_photon_dec;
		std::vector<double>  *trk_l0_photon_tis;
		std::vector<double>  *trk_l0_photon_tos;
		std::vector<double>  *trk_l0_electron_dec;
		std::vector<double>  *trk_l0_electron_tis;
		std::vector<double>  *trk_l0_electron_tos;
		std::vector<double>  *trk_l0_muon_dec;
		std::vector<double>  *trk_l0_muon_tis;
		std::vector<double>  *trk_l0_muon_tos;
		std::vector<double>  *trk_l0_dimuon_dec;
		std::vector<double>  *trk_l0_dimuon_tis;
		std::vector<double>  *trk_l0_dimuon_tos;
		std::vector<double>  *neu_idx_gen;
		std::vector<double>  *neu_p;
		std::vector<double>  *neu_pt;
		std::vector<double>  *neu_px;
		std::vector<double>  *neu_py;
		std::vector<double>  *neu_pz;
		std::vector<double>  *neu_e;
		std::vector<double>  *neu_pid;
		std::vector<double>  *neu_cl;
		std::vector<double>  *neu_l0_hadron_dec;
		std::vector<double>  *neu_l0_hadron_tis;
		std::vector<double>  *neu_l0_hadron_tos;
		std::vector<double>  *neu_l0_photon_dec;
		std::vector<double>  *neu_l0_photon_tis;
		std::vector<double>  *neu_l0_photon_tos;
		std::vector<double>  *neu_l0_electron_dec;
		std::vector<double>  *neu_l0_electron_tis;
		std::vector<double>  *neu_l0_electron_tos;
		std::vector<double>  *neu_l0_muon_dec;
		std::vector<double>  *neu_l0_muon_tis;
		std::vector<double>  *neu_l0_muon_tos;
		std::vector<double>  *neu_l0_dimuon_dec;
		std::vector<double>  *neu_l0_dimuon_tis;
		std::vector<double>  *neu_l0_dimuon_tos;
		std::vector<double>  *b_idx_pvr;
		std::vector<double>  *b_p;
		std::vector<double>  *b_pt;
		std::vector<double>  *b_px;
		std::vector<double>  *b_py;
		std::vector<double>  *b_pz;
		std::vector<double>  *b_e;
		std::vector<double>  *b_x;
		std::vector<double>  *b_y;
		std::vector<double>  *b_z;
		std::vector<double>  *b_m;
		std::vector<double>  *b_vtx_chi2;
		std::vector<double>  *b_vtx_ndof;
		std::vector<double>  *b_ip;
		std::vector<double>  *b_ip_chi2;
		std::vector<double>  *b_fd;
		std::vector<double>  *b_fd_chi2;
		std::vector<double>  *b_tau;
		std::vector<double>  *b_tau_err;
		std::vector<double>  *b_tau_chi2;
		std::vector<double>  *b_idx_kst;
		std::vector<double>  *b_idx_k;
		std::vector<double>  *b_idx_pi;
		std::vector<double>  *b_idx_x;
		std::vector<double>  *b_l0_hadron_dec;
		std::vector<double>  *b_l0_hadron_tis;
		std::vector<double>  *b_l0_hadron_tos;
		std::vector<double>  *b_l0_photon_dec;
		std::vector<double>  *b_l0_photon_tis;
		std::vector<double>  *b_l0_photon_tos;
		std::vector<double>  *b_l0_electron_dec;
		std::vector<double>  *b_l0_electron_tis;
		std::vector<double>  *b_l0_electron_tos;
		std::vector<double>  *b_l0_muon_dec;
		std::vector<double>  *b_l0_muon_tis;
		std::vector<double>  *b_l0_muon_tos;
		std::vector<double>  *b_l0_dimuon_dec;
		std::vector<double>  *b_l0_dimuon_tis;
		std::vector<double>  *b_l0_dimuon_tos;
		std::vector<double>  *b_mode;
		std::vector<double>  *b_dtf_chi2;
		std::vector<double>  *b_dtf_ndof;
		std::vector<double>  *b_dtf_px;
		std::vector<double>  *b_dtf_py;
		std::vector<double>  *b_dtf_pz;
		std::vector<double>  *b_dtf_e;
		std::vector<double>  *b_dtf_m;
		std::vector<double>  *b_dtf_ctau;
		std::vector<double>  *b_dtf_x_m;
		std::vector<double>  *b_dtf_x_ctau;
		std::vector<double>  *b_dtf_trk0_px;
		std::vector<double>  *b_dtf_trk0_py;
		std::vector<double>  *b_dtf_trk0_pz;
		std::vector<double>  *b_dtf_trk0_e;
		std::vector<double>  *b_dtf_trk0_m;
		std::vector<double>  *b_dtf_trk1_px;
		std::vector<double>  *b_dtf_trk1_py;
		std::vector<double>  *b_dtf_trk1_pz;
		std::vector<double>  *b_dtf_trk1_e;
		std::vector<double>  *b_dtf_trk1_m;
		std::vector<double>  *b_dtf_trk2_px;
		std::vector<double>  *b_dtf_trk2_py;
		std::vector<double>  *b_dtf_trk2_pz;
		std::vector<double>  *b_dtf_trk2_e;
		std::vector<double>  *b_dtf_trk2_m;
		std::vector<double>  *b_dtf_trk3_px;
		std::vector<double>  *b_dtf_trk3_py;
		std::vector<double>  *b_dtf_trk3_pz;
		std::vector<double>  *b_dtf_trk3_e;
		std::vector<double>  *b_dtf_trk3_m;
		std::vector<double>  *b_dtf_trk4_px;
		std::vector<double>  *b_dtf_trk4_py;
		std::vector<double>  *b_dtf_trk4_pz;
		std::vector<double>  *b_dtf_trk4_e;
		std::vector<double>  *b_dtf_trk4_m;
		std::vector<double>  *b_dtf_trk5_px;
		std::vector<double>  *b_dtf_trk5_py;
		std::vector<double>  *b_dtf_trk5_pz;
		std::vector<double>  *b_dtf_trk5_e;
		std::vector<double>  *b_dtf_trk5_m;
		std::vector<double>  *b_dtf_trk6_px;
		std::vector<double>  *b_dtf_trk6_py;
		std::vector<double>  *b_dtf_trk6_pz;
		std::vector<double>  *b_dtf_trk6_e;
		std::vector<double>  *b_dtf_trk6_m;
		std::vector<double>  *b_dtf_trk7_px;
		std::vector<double>  *b_dtf_trk7_py;
		std::vector<double>  *b_dtf_trk7_pz;
		std::vector<double>  *b_dtf_trk7_e;
		std::vector<double>  *b_dtf_trk7_m;
		std::vector<double>  *b_dtf_neu0_px;
		std::vector<double>  *b_dtf_neu0_py;
		std::vector<double>  *b_dtf_neu0_pz;
		std::vector<double>  *b_dtf_neu0_e;
		std::vector<double>  *b_dtf_neu0_m;
		std::vector<double>  *b_dtf_neu1_px;
		std::vector<double>  *b_dtf_neu1_py;
		std::vector<double>  *b_dtf_neu1_pz;
		std::vector<double>  *b_dtf_neu1_e;
		std::vector<double>  *b_dtf_neu1_m;
		std::vector<double>  *b_dtf_neu2_px;
		std::vector<double>  *b_dtf_neu2_py;
		std::vector<double>  *b_dtf_neu2_pz;
		std::vector<double>  *b_dtf_neu2_e;
		std::vector<double>  *b_dtf_neu2_m;
		std::vector<double>  *b_dtf_neu3_px;
		std::vector<double>  *b_dtf_neu3_py;
		std::vector<double>  *b_dtf_neu3_pz;
		std::vector<double>  *b_dtf_neu3_e;
		std::vector<double>  *b_dtf_neu3_m;
		std::vector<double>  *b_dtf_neu4_px;
		std::vector<double>  *b_dtf_neu4_py;
		std::vector<double>  *b_dtf_neu4_pz;
		std::vector<double>  *b_dtf_neu4_e;
		std::vector<double>  *b_dtf_neu4_m;
		std::vector<double>  *b_dtf_neu5_px;
		std::vector<double>  *b_dtf_neu5_py;
		std::vector<double>  *b_dtf_neu5_pz;
		std::vector<double>  *b_dtf_neu5_e;
		std::vector<double>  *b_dtf_neu5_m;
		std::vector<double>  *b_dtf_neu6_px;
		std::vector<double>  *b_dtf_neu6_py;
		std::vector<double>  *b_dtf_neu6_pz;
		std::vector<double>  *b_dtf_neu6_e;
		std::vector<double>  *b_dtf_neu6_m;
		std::vector<double>  *b_dtf_neu7_px;
		std::vector<double>  *b_dtf_neu7_py;
		std::vector<double>  *b_dtf_neu7_pz;
		std::vector<double>  *b_dtf_neu7_e;
		std::vector<double>  *b_dtf_neu7_m;
		std::vector<double>  *kst_idx_pvr;
		std::vector<double>  *kst_p;
		std::vector<double>  *kst_pt;
		std::vector<double>  *kst_px;
		std::vector<double>  *kst_py;
		std::vector<double>  *kst_pz;
		std::vector<double>  *kst_e;
		std::vector<double>  *kst_x;
		std::vector<double>  *kst_y;
		std::vector<double>  *kst_z;
		std::vector<double>  *kst_m;
		std::vector<double>  *kst_vtx_chi2;
		std::vector<double>  *kst_vtx_ndof;
		std::vector<double>  *kst_ip;
		std::vector<double>  *kst_ip_best;
		std::vector<double>  *kst_ip_bvtx;
		std::vector<double>  *kst_ip_chi2;
		std::vector<double>  *kst_ip_chi2_best;
		std::vector<double>  *kst_ip_chi2_bvtx;
		std::vector<double>  *kst_fd;
		std::vector<double>  *kst_fd_best;
		std::vector<double>  *kst_fd_bvtx;
		std::vector<double>  *kst_fd_chi2;
		std::vector<double>  *kst_fd_chi2_best;
		std::vector<double>  *kst_fd_chi2_bvtx;
		std::vector<double>  *kst_tau;
		std::vector<double>  *kst_tau_err;
		std::vector<double>  *kst_tau_chi2;
		std::vector<double>  *kst_idx_k;
		std::vector<double>  *kst_idx_pi;
		std::vector<double>  *kst_l0_hadron_dec;
		std::vector<double>  *kst_l0_hadron_tis;
		std::vector<double>  *kst_l0_hadron_tos;
		std::vector<double>  *kst_l0_photon_dec;
		std::vector<double>  *kst_l0_photon_tis;
		std::vector<double>  *kst_l0_photon_tos;
		std::vector<double>  *kst_l0_electron_dec;
		std::vector<double>  *kst_l0_electron_tis;
		std::vector<double>  *kst_l0_electron_tos;
		std::vector<double>  *kst_l0_muon_dec;
		std::vector<double>  *kst_l0_muon_tis;
		std::vector<double>  *kst_l0_muon_tos;
		std::vector<double>  *kst_l0_dimuon_dec;
		std::vector<double>  *kst_l0_dimuon_tis;
		std::vector<double>  *kst_l0_dimuon_tos;
		std::vector<double>  *x_idx_pvr;
		std::vector<double>  *x_p;
		std::vector<double>  *x_pt;
		std::vector<double>  *x_px;
		std::vector<double>  *x_py;
		std::vector<double>  *x_pz;
		std::vector<double>  *x_e;
		std::vector<double>  *x_x;
		std::vector<double>  *x_y;
		std::vector<double>  *x_z;
		std::vector<double>  *x_m;
		std::vector<double>  *x_vtx_chi2;
		std::vector<double>  *x_vtx_ndof;
		std::vector<double>  *x_ip;
		std::vector<double>  *x_ip_best;
		std::vector<double>  *x_ip_bvtx;
		std::vector<double>  *x_ip_chi2;
		std::vector<double>  *x_ip_chi2_best;
		std::vector<double>  *x_ip_chi2_bvtx;
		std::vector<double>  *x_fd;
		std::vector<double>  *x_fd_best;
		std::vector<double>  *x_fd_bvtx;
		std::vector<double>  *x_fd_chi2;
		std::vector<double>  *x_fd_chi2_best;
		std::vector<double>  *x_fd_chi2_bvtx;
		std::vector<double>  *x_tau;
		std::vector<double>  *x_tau_err;
		std::vector<double>  *x_tau_chi2;
		std::vector<double>  *x_idx_trk0;
		std::vector<double>  *x_idx_trk1;
		std::vector<double>  *x_idx_trk2;
		std::vector<double>  *x_idx_trk3;
		std::vector<double>  *x_idx_trk4;
		std::vector<double>  *x_idx_trk5;
		std::vector<double>  *x_idx_neu0;
		std::vector<double>  *x_idx_neu1;
		std::vector<double>  *x_idx_neu2;
		std::vector<double>  *x_idx_neu3;
		std::vector<double>  *x_idx_neu4;
		std::vector<double>  *x_idx_neu5;
		std::vector<double>  *x_l0_hadron_dec;
		std::vector<double>  *x_l0_hadron_tis;
		std::vector<double>  *x_l0_hadron_tos;
		std::vector<double>  *x_l0_photon_dec;
		std::vector<double>  *x_l0_photon_tis;
		std::vector<double>  *x_l0_photon_tos;
		std::vector<double>  *x_l0_electron_dec;
		std::vector<double>  *x_l0_electron_tis;
		std::vector<double>  *x_l0_electron_tos;
		std::vector<double>  *x_l0_muon_dec;
		std::vector<double>  *x_l0_muon_tis;
		std::vector<double>  *x_l0_muon_tos;
		std::vector<double>  *x_l0_dimuon_dec;
		std::vector<double>  *x_l0_dimuon_tis;
		std::vector<double>  *x_l0_dimuon_tos;
		Double_t        evt_pvr_n;
		Double_t        evt_trk_n;

		// List of branches
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
		TBranch        *b_trk_ip_best;   //!
		TBranch        *b_trk_ip_bvtx;   //!
		TBranch        *b_trk_ip_chi2;   //!
		TBranch        *b_trk_ip_chi2_best;   //!
		TBranch        *b_trk_ip_chi2_bvtx;   //!
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
		TBranch        *b_trk_l0_photon_dec;   //!
		TBranch        *b_trk_l0_photon_tis;   //!
		TBranch        *b_trk_l0_photon_tos;   //!
		TBranch        *b_trk_l0_electron_dec;   //!
		TBranch        *b_trk_l0_electron_tis;   //!
		TBranch        *b_trk_l0_electron_tos;   //!
		TBranch        *b_trk_l0_muon_dec;   //!
		TBranch        *b_trk_l0_muon_tis;   //!
		TBranch        *b_trk_l0_muon_tos;   //!
		TBranch        *b_trk_l0_dimuon_dec;   //!
		TBranch        *b_trk_l0_dimuon_tis;   //!
		TBranch        *b_trk_l0_dimuon_tos;   //!
		TBranch        *b_neu_idx_gen;   //!
		TBranch        *b_neu_p;   //!
		TBranch        *b_neu_pt;   //!
		TBranch        *b_neu_px;   //!
		TBranch        *b_neu_py;   //!
		TBranch        *b_neu_pz;   //!
		TBranch        *b_neu_e;   //!
		TBranch        *b_neu_pid;   //!
		TBranch        *b_neu_cl;   //!
		TBranch        *b_neu_l0_hadron_dec;   //!
		TBranch        *b_neu_l0_hadron_tis;   //!
		TBranch        *b_neu_l0_hadron_tos;   //!
		TBranch        *b_neu_l0_photon_dec;   //!
		TBranch        *b_neu_l0_photon_tis;   //!
		TBranch        *b_neu_l0_photon_tos;   //!
		TBranch        *b_neu_l0_electron_dec;   //!
		TBranch        *b_neu_l0_electron_tis;   //!
		TBranch        *b_neu_l0_electron_tos;   //!
		TBranch        *b_neu_l0_muon_dec;   //!
		TBranch        *b_neu_l0_muon_tis;   //!
		TBranch        *b_neu_l0_muon_tos;   //!
		TBranch        *b_neu_l0_dimuon_dec;   //!
		TBranch        *b_neu_l0_dimuon_tis;   //!
		TBranch        *b_neu_l0_dimuon_tos;   //!
		TBranch        *b_b_idx_pvr;   //!
		TBranch        *b_b_p;   //!
		TBranch        *b_b_pt;   //!
		TBranch        *b_b_px;   //!
		TBranch        *b_b_py;   //!
		TBranch        *b_b_pz;   //!
		TBranch        *b_b_e;   //!
		TBranch        *b_b_x;   //!
		TBranch        *b_b_y;   //!
		TBranch        *b_b_z;   //!
		TBranch        *b_b_m;   //!
		TBranch        *b_b_vtx_chi2;   //!
		TBranch        *b_b_vtx_ndof;   //!
		TBranch        *b_b_ip;   //!
		TBranch        *b_b_ip_chi2;   //!
		TBranch        *b_b_fd;   //!
		TBranch        *b_b_fd_chi2;   //!
		TBranch        *b_b_tau;   //!
		TBranch        *b_b_tau_err;   //!
		TBranch        *b_b_tau_chi2;   //!
		TBranch        *b_b_idx_kst;   //!
		TBranch        *b_b_idx_k;   //!
		TBranch        *b_b_idx_pi;   //!
		TBranch        *b_b_idx_x;   //!
		TBranch        *b_b_l0_hadron_dec;   //!
		TBranch        *b_b_l0_hadron_tis;   //!
		TBranch        *b_b_l0_hadron_tos;   //!
		TBranch        *b_b_l0_photon_dec;   //!
		TBranch        *b_b_l0_photon_tis;   //!
		TBranch        *b_b_l0_photon_tos;   //!
		TBranch        *b_b_l0_electron_dec;   //!
		TBranch        *b_b_l0_electron_tis;   //!
		TBranch        *b_b_l0_electron_tos;   //!
		TBranch        *b_b_l0_muon_dec;   //!
		TBranch        *b_b_l0_muon_tis;   //!
		TBranch        *b_b_l0_muon_tos;   //!
		TBranch        *b_b_l0_dimuon_dec;   //!
		TBranch        *b_b_l0_dimuon_tis;   //!
		TBranch        *b_b_l0_dimuon_tos;   //!
		TBranch        *b_b_mode;   //!
		TBranch        *b_b_dtf_chi2;   //!
		TBranch        *b_b_dtf_ndof;   //!
		TBranch        *b_b_dtf_px;   //!
		TBranch        *b_b_dtf_py;   //!
		TBranch        *b_b_dtf_pz;   //!
		TBranch        *b_b_dtf_e;   //!
		TBranch        *b_b_dtf_m;   //!
		TBranch        *b_b_dtf_ctau;   //!
		TBranch        *b_b_dtf_x_m;   //!
		TBranch        *b_b_dtf_x_ctau;   //!
		TBranch        *b_b_dtf_trk0_px;   //!
		TBranch        *b_b_dtf_trk0_py;   //!
		TBranch        *b_b_dtf_trk0_pz;   //!
		TBranch        *b_b_dtf_trk0_e;   //!
		TBranch        *b_b_dtf_trk0_m;   //!
		TBranch        *b_b_dtf_trk1_px;   //!
		TBranch        *b_b_dtf_trk1_py;   //!
		TBranch        *b_b_dtf_trk1_pz;   //!
		TBranch        *b_b_dtf_trk1_e;   //!
		TBranch        *b_b_dtf_trk1_m;   //!
		TBranch        *b_b_dtf_trk2_px;   //!
		TBranch        *b_b_dtf_trk2_py;   //!
		TBranch        *b_b_dtf_trk2_pz;   //!
		TBranch        *b_b_dtf_trk2_e;   //!
		TBranch        *b_b_dtf_trk2_m;   //!
		TBranch        *b_b_dtf_trk3_px;   //!
		TBranch        *b_b_dtf_trk3_py;   //!
		TBranch        *b_b_dtf_trk3_pz;   //!
		TBranch        *b_b_dtf_trk3_e;   //!
		TBranch        *b_b_dtf_trk3_m;   //!
		TBranch        *b_b_dtf_trk4_px;   //!
		TBranch        *b_b_dtf_trk4_py;   //!
		TBranch        *b_b_dtf_trk4_pz;   //!
		TBranch        *b_b_dtf_trk4_e;   //!
		TBranch        *b_b_dtf_trk4_m;   //!
		TBranch        *b_b_dtf_trk5_px;   //!
		TBranch        *b_b_dtf_trk5_py;   //!
		TBranch        *b_b_dtf_trk5_pz;   //!
		TBranch        *b_b_dtf_trk5_e;   //!
		TBranch        *b_b_dtf_trk5_m;   //!
		TBranch        *b_b_dtf_trk6_px;   //!
		TBranch        *b_b_dtf_trk6_py;   //!
		TBranch        *b_b_dtf_trk6_pz;   //!
		TBranch        *b_b_dtf_trk6_e;   //!
		TBranch        *b_b_dtf_trk6_m;   //!
		TBranch        *b_b_dtf_trk7_px;   //!
		TBranch        *b_b_dtf_trk7_py;   //!
		TBranch        *b_b_dtf_trk7_pz;   //!
		TBranch        *b_b_dtf_trk7_e;   //!
		TBranch        *b_b_dtf_trk7_m;   //!
		TBranch        *b_b_dtf_neu0_px;   //!
		TBranch        *b_b_dtf_neu0_py;   //!
		TBranch        *b_b_dtf_neu0_pz;   //!
		TBranch        *b_b_dtf_neu0_e;   //!
		TBranch        *b_b_dtf_neu0_m;   //!
		TBranch        *b_b_dtf_neu1_px;   //!
		TBranch        *b_b_dtf_neu1_py;   //!
		TBranch        *b_b_dtf_neu1_pz;   //!
		TBranch        *b_b_dtf_neu1_e;   //!
		TBranch        *b_b_dtf_neu1_m;   //!
		TBranch        *b_b_dtf_neu2_px;   //!
		TBranch        *b_b_dtf_neu2_py;   //!
		TBranch        *b_b_dtf_neu2_pz;   //!
		TBranch        *b_b_dtf_neu2_e;   //!
		TBranch        *b_b_dtf_neu2_m;   //!
		TBranch        *b_b_dtf_neu3_px;   //!
		TBranch        *b_b_dtf_neu3_py;   //!
		TBranch        *b_b_dtf_neu3_pz;   //!
		TBranch        *b_b_dtf_neu3_e;   //!
		TBranch        *b_b_dtf_neu3_m;   //!
		TBranch        *b_b_dtf_neu4_px;   //!
		TBranch        *b_b_dtf_neu4_py;   //!
		TBranch        *b_b_dtf_neu4_pz;   //!
		TBranch        *b_b_dtf_neu4_e;   //!
		TBranch        *b_b_dtf_neu4_m;   //!
		TBranch        *b_b_dtf_neu5_px;   //!
		TBranch        *b_b_dtf_neu5_py;   //!
		TBranch        *b_b_dtf_neu5_pz;   //!
		TBranch        *b_b_dtf_neu5_e;   //!
		TBranch        *b_b_dtf_neu5_m;   //!
		TBranch        *b_b_dtf_neu6_px;   //!
		TBranch        *b_b_dtf_neu6_py;   //!
		TBranch        *b_b_dtf_neu6_pz;   //!
		TBranch        *b_b_dtf_neu6_e;   //!
		TBranch        *b_b_dtf_neu6_m;   //!
		TBranch        *b_b_dtf_neu7_px;   //!
		TBranch        *b_b_dtf_neu7_py;   //!
		TBranch        *b_b_dtf_neu7_pz;   //!
		TBranch        *b_b_dtf_neu7_e;   //!
		TBranch        *b_b_dtf_neu7_m;   //!
		TBranch        *b_kst_idx_pvr;   //!
		TBranch        *b_kst_p;   //!
		TBranch        *b_kst_pt;   //!
		TBranch        *b_kst_px;   //!
		TBranch        *b_kst_py;   //!
		TBranch        *b_kst_pz;   //!
		TBranch        *b_kst_e;   //!
		TBranch        *b_kst_x;   //!
		TBranch        *b_kst_y;   //!
		TBranch        *b_kst_z;   //!
		TBranch        *b_kst_m;   //!
		TBranch        *b_kst_vtx_chi2;   //!
		TBranch        *b_kst_vtx_ndof;   //!
		TBranch        *b_kst_ip;   //!
		TBranch        *b_kst_ip_best;   //!
		TBranch        *b_kst_ip_bvtx;   //!
		TBranch        *b_kst_ip_chi2;   //!
		TBranch        *b_kst_ip_chi2_best;   //!
		TBranch        *b_kst_ip_chi2_bvtx;   //!
		TBranch        *b_kst_fd;   //!
		TBranch        *b_kst_fd_best;   //!
		TBranch        *b_kst_fd_bvtx;   //!
		TBranch        *b_kst_fd_chi2;   //!
		TBranch        *b_kst_fd_chi2_best;   //!
		TBranch        *b_kst_fd_chi2_bvtx;   //!
		TBranch        *b_kst_tau;   //!
		TBranch        *b_kst_tau_err;   //!
		TBranch        *b_kst_tau_chi2;   //!
		TBranch        *b_kst_idx_k;   //!
		TBranch        *b_kst_idx_pi;   //!
		TBranch        *b_kst_l0_hadron_dec;   //!
		TBranch        *b_kst_l0_hadron_tis;   //!
		TBranch        *b_kst_l0_hadron_tos;   //!
		TBranch        *b_kst_l0_photon_dec;   //!
		TBranch        *b_kst_l0_photon_tis;   //!
		TBranch        *b_kst_l0_photon_tos;   //!
		TBranch        *b_kst_l0_electron_dec;   //!
		TBranch        *b_kst_l0_electron_tis;   //!
		TBranch        *b_kst_l0_electron_tos;   //!
		TBranch        *b_kst_l0_muon_dec;   //!
		TBranch        *b_kst_l0_muon_tis;   //!
		TBranch        *b_kst_l0_muon_tos;   //!
		TBranch        *b_kst_l0_dimuon_dec;   //!
		TBranch        *b_kst_l0_dimuon_tis;   //!
		TBranch        *b_kst_l0_dimuon_tos;   //!
		TBranch        *b_x_idx_pvr;   //!
		TBranch        *b_x_p;   //!
		TBranch        *b_x_pt;   //!
		TBranch        *b_x_px;   //!
		TBranch        *b_x_py;   //!
		TBranch        *b_x_pz;   //!
		TBranch        *b_x_e;   //!
		TBranch        *b_x_x;   //!
		TBranch        *b_x_y;   //!
		TBranch        *b_x_z;   //!
		TBranch        *b_x_m;   //!
		TBranch        *b_x_vtx_chi2;   //!
		TBranch        *b_x_vtx_ndof;   //!
		TBranch        *b_x_ip;   //!
		TBranch        *b_x_ip_best;   //!
		TBranch        *b_x_ip_bvtx;   //!
		TBranch        *b_x_ip_chi2;   //!
		TBranch        *b_x_ip_chi2_best;   //!
		TBranch        *b_x_ip_chi2_bvtx;   //!
		TBranch        *b_x_fd;   //!
		TBranch        *b_x_fd_best;   //!
		TBranch        *b_x_fd_bvtx;   //!
		TBranch        *b_x_fd_chi2;   //!
		TBranch        *b_x_fd_chi2_best;   //!
		TBranch        *b_x_fd_chi2_bvtx;   //!
		TBranch        *b_x_tau;   //!
		TBranch        *b_x_tau_err;   //!
		TBranch        *b_x_tau_chi2;   //!
		TBranch        *b_x_idx_trk0;   //!
		TBranch        *b_x_idx_trk1;   //!
		TBranch        *b_x_idx_trk2;   //!
		TBranch        *b_x_idx_trk3;   //!
		TBranch        *b_x_idx_trk4;   //!
		TBranch        *b_x_idx_trk5;   //!
		TBranch        *b_x_idx_neu0;   //!
		TBranch        *b_x_idx_neu1;   //!
		TBranch        *b_x_idx_neu2;   //!
		TBranch        *b_x_idx_neu3;   //!
		TBranch        *b_x_idx_neu4;   //!
		TBranch        *b_x_idx_neu5;   //!
		TBranch        *b_x_l0_hadron_dec;   //!
		TBranch        *b_x_l0_hadron_tis;   //!
		TBranch        *b_x_l0_hadron_tos;   //!
		TBranch        *b_x_l0_photon_dec;   //!
		TBranch        *b_x_l0_photon_tis;   //!
		TBranch        *b_x_l0_photon_tos;   //!
		TBranch        *b_x_l0_electron_dec;   //!
		TBranch        *b_x_l0_electron_tis;   //!
		TBranch        *b_x_l0_electron_tos;   //!
		TBranch        *b_x_l0_muon_dec;   //!
		TBranch        *b_x_l0_muon_tis;   //!
		TBranch        *b_x_l0_muon_tos;   //!
		TBranch        *b_x_l0_dimuon_dec;   //!
		TBranch        *b_x_l0_dimuon_tis;   //!
		TBranch        *b_x_l0_dimuon_tos;   //!
		TBranch        *b_evt_pvr_n;   //!
		TBranch        *b_evt_trk_n;   //!

		TreeSkimmer(int dir);
		virtual ~TreeSkimmer();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop_Kpi_pipineu(int mode, int nmax=-1);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da,
				const TVector3 &pb, const TVector3 &db);
		double calcDoca(const TVector3 &pa, const TVector3 &da, const TVector3 &pb);

		BPart* addB(TTree* tree, TString pre);
		XPart* addX(TTree* tree, TString pre);
		Track* addTrk(TTree* tree, TString pre);
		Neutral* addNeu(TTree* tree, TString pre);
		DTFVars* addDTF(TTree* tree, TString pre,
		                std::vector<TString> trks,
			        std::vector<TString> neus);

		void fillB(BPart* part, int idx);
		void fillX(XPart* part, int idx);
		void fillTrk(Track* part, int idx);
		void fillNeu(Neutral* part, int idx);
		void fillDTF(DTFVars* dtf, int idx);

		TString savedir{""};
};

#endif

#ifdef TreeSkimmer_cxx
TreeSkimmer::TreeSkimmer(int dir) : fChain(0) 
{
	TChain* tree = new TChain("data");
	for(int j=47; j<55; ++j) {
		if(dir>=0 && dir!=j) continue;
		for(int i=0; i<700; ++i) {
			TString str= TString::Format("/data/alps/%d/%d/output.root",j,i);
			if(gSystem->AccessPathName(str)) continue;
			tree->Add(str);
		}
	}

	Init(tree);
	savedir+="/data/alps/";
	if(dir>=0) {
		savedir += dir;
		savedir+="/";
	}
}

TreeSkimmer::~TreeSkimmer()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t TreeSkimmer::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t TreeSkimmer::LoadTree(Long64_t entry)
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

void TreeSkimmer::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
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
	trk_ip_best = 0;
	trk_ip_bvtx = 0;
	trk_ip_chi2 = 0;
	trk_ip_chi2_best = 0;
	trk_ip_chi2_bvtx = 0;
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
	trk_l0_photon_dec = 0;
	trk_l0_photon_tis = 0;
	trk_l0_photon_tos = 0;
	trk_l0_electron_dec = 0;
	trk_l0_electron_tis = 0;
	trk_l0_electron_tos = 0;
	trk_l0_muon_dec = 0;
	trk_l0_muon_tis = 0;
	trk_l0_muon_tos = 0;
	trk_l0_dimuon_dec = 0;
	trk_l0_dimuon_tis = 0;
	trk_l0_dimuon_tos = 0;
	neu_idx_gen = 0;
	neu_p = 0;
	neu_pt = 0;
	neu_px = 0;
	neu_py = 0;
	neu_pz = 0;
	neu_e = 0;
	neu_pid = 0;
	neu_cl = 0;
	neu_l0_hadron_dec = 0;
	neu_l0_hadron_tis = 0;
	neu_l0_hadron_tos = 0;
	neu_l0_photon_dec = 0;
	neu_l0_photon_tis = 0;
	neu_l0_photon_tos = 0;
	neu_l0_electron_dec = 0;
	neu_l0_electron_tis = 0;
	neu_l0_electron_tos = 0;
	neu_l0_muon_dec = 0;
	neu_l0_muon_tis = 0;
	neu_l0_muon_tos = 0;
	neu_l0_dimuon_dec = 0;
	neu_l0_dimuon_tis = 0;
	neu_l0_dimuon_tos = 0;
	b_idx_pvr = 0;
	b_p = 0;
	b_pt = 0;
	b_px = 0;
	b_py = 0;
	b_pz = 0;
	b_e = 0;
	b_x = 0;
	b_y = 0;
	b_z = 0;
	b_m = 0;
	b_vtx_chi2 = 0;
	b_vtx_ndof = 0;
	b_ip = 0;
	b_ip_chi2 = 0;
	b_fd = 0;
	b_fd_chi2 = 0;
	b_tau = 0;
	b_tau_err = 0;
	b_tau_chi2 = 0;
	b_idx_kst = 0;
	b_idx_k = 0;
	b_idx_pi = 0;
	b_idx_x = 0;
	b_l0_hadron_dec = 0;
	b_l0_hadron_tis = 0;
	b_l0_hadron_tos = 0;
	b_l0_photon_dec = 0;
	b_l0_photon_tis = 0;
	b_l0_photon_tos = 0;
	b_l0_electron_dec = 0;
	b_l0_electron_tis = 0;
	b_l0_electron_tos = 0;
	b_l0_muon_dec = 0;
	b_l0_muon_tis = 0;
	b_l0_muon_tos = 0;
	b_l0_dimuon_dec = 0;
	b_l0_dimuon_tis = 0;
	b_l0_dimuon_tos = 0;
	b_mode = 0;
	b_dtf_chi2 = 0;
	b_dtf_ndof = 0;
	b_dtf_px = 0;
	b_dtf_py = 0;
	b_dtf_pz = 0;
	b_dtf_e = 0;
	b_dtf_m = 0;
	b_dtf_ctau = 0;
	b_dtf_x_m = 0;
	b_dtf_x_ctau = 0;
	b_dtf_trk0_px = 0;
	b_dtf_trk0_py = 0;
	b_dtf_trk0_pz = 0;
	b_dtf_trk0_e = 0;
	b_dtf_trk0_m = 0;
	b_dtf_trk1_px = 0;
	b_dtf_trk1_py = 0;
	b_dtf_trk1_pz = 0;
	b_dtf_trk1_e = 0;
	b_dtf_trk1_m = 0;
	b_dtf_trk2_px = 0;
	b_dtf_trk2_py = 0;
	b_dtf_trk2_pz = 0;
	b_dtf_trk2_e = 0;
	b_dtf_trk2_m = 0;
	b_dtf_trk3_px = 0;
	b_dtf_trk3_py = 0;
	b_dtf_trk3_pz = 0;
	b_dtf_trk3_e = 0;
	b_dtf_trk3_m = 0;
	b_dtf_trk4_px = 0;
	b_dtf_trk4_py = 0;
	b_dtf_trk4_pz = 0;
	b_dtf_trk4_e = 0;
	b_dtf_trk4_m = 0;
	b_dtf_trk5_px = 0;
	b_dtf_trk5_py = 0;
	b_dtf_trk5_pz = 0;
	b_dtf_trk5_e = 0;
	b_dtf_trk5_m = 0;
	b_dtf_trk6_px = 0;
	b_dtf_trk6_py = 0;
	b_dtf_trk6_pz = 0;
	b_dtf_trk6_e = 0;
	b_dtf_trk6_m = 0;
	b_dtf_trk7_px = 0;
	b_dtf_trk7_py = 0;
	b_dtf_trk7_pz = 0;
	b_dtf_trk7_e = 0;
	b_dtf_trk7_m = 0;
	b_dtf_neu0_px = 0;
	b_dtf_neu0_py = 0;
	b_dtf_neu0_pz = 0;
	b_dtf_neu0_e = 0;
	b_dtf_neu0_m = 0;
	b_dtf_neu1_px = 0;
	b_dtf_neu1_py = 0;
	b_dtf_neu1_pz = 0;
	b_dtf_neu1_e = 0;
	b_dtf_neu1_m = 0;
	b_dtf_neu2_px = 0;
	b_dtf_neu2_py = 0;
	b_dtf_neu2_pz = 0;
	b_dtf_neu2_e = 0;
	b_dtf_neu2_m = 0;
	b_dtf_neu3_px = 0;
	b_dtf_neu3_py = 0;
	b_dtf_neu3_pz = 0;
	b_dtf_neu3_e = 0;
	b_dtf_neu3_m = 0;
	b_dtf_neu4_px = 0;
	b_dtf_neu4_py = 0;
	b_dtf_neu4_pz = 0;
	b_dtf_neu4_e = 0;
	b_dtf_neu4_m = 0;
	b_dtf_neu5_px = 0;
	b_dtf_neu5_py = 0;
	b_dtf_neu5_pz = 0;
	b_dtf_neu5_e = 0;
	b_dtf_neu5_m = 0;
	b_dtf_neu6_px = 0;
	b_dtf_neu6_py = 0;
	b_dtf_neu6_pz = 0;
	b_dtf_neu6_e = 0;
	b_dtf_neu6_m = 0;
	b_dtf_neu7_px = 0;
	b_dtf_neu7_py = 0;
	b_dtf_neu7_pz = 0;
	b_dtf_neu7_e = 0;
	b_dtf_neu7_m = 0;
	kst_idx_pvr = 0;
	kst_p = 0;
	kst_pt = 0;
	kst_px = 0;
	kst_py = 0;
	kst_pz = 0;
	kst_e = 0;
	kst_x = 0;
	kst_y = 0;
	kst_z = 0;
	kst_m = 0;
	kst_vtx_chi2 = 0;
	kst_vtx_ndof = 0;
	kst_ip = 0;
	kst_ip_best = 0;
	kst_ip_bvtx = 0;
	kst_ip_chi2 = 0;
	kst_ip_chi2_best = 0;
	kst_ip_chi2_bvtx = 0;
	kst_fd = 0;
	kst_fd_best = 0;
	kst_fd_bvtx = 0;
	kst_fd_chi2 = 0;
	kst_fd_chi2_best = 0;
	kst_fd_chi2_bvtx = 0;
	kst_tau = 0;
	kst_tau_err = 0;
	kst_tau_chi2 = 0;
	kst_idx_k = 0;
	kst_idx_pi = 0;
	kst_l0_hadron_dec = 0;
	kst_l0_hadron_tis = 0;
	kst_l0_hadron_tos = 0;
	kst_l0_photon_dec = 0;
	kst_l0_photon_tis = 0;
	kst_l0_photon_tos = 0;
	kst_l0_electron_dec = 0;
	kst_l0_electron_tis = 0;
	kst_l0_electron_tos = 0;
	kst_l0_muon_dec = 0;
	kst_l0_muon_tis = 0;
	kst_l0_muon_tos = 0;
	kst_l0_dimuon_dec = 0;
	kst_l0_dimuon_tis = 0;
	kst_l0_dimuon_tos = 0;
	x_idx_pvr = 0;
	x_p = 0;
	x_pt = 0;
	x_px = 0;
	x_py = 0;
	x_pz = 0;
	x_e = 0;
	x_x = 0;
	x_y = 0;
	x_z = 0;
	x_m = 0;
	x_vtx_chi2 = 0;
	x_vtx_ndof = 0;
	x_ip = 0;
	x_ip_best = 0;
	x_ip_bvtx = 0;
	x_ip_chi2 = 0;
	x_ip_chi2_best = 0;
	x_ip_chi2_bvtx = 0;
	x_fd = 0;
	x_fd_best = 0;
	x_fd_bvtx = 0;
	x_fd_chi2 = 0;
	x_fd_chi2_best = 0;
	x_fd_chi2_bvtx = 0;
	x_tau = 0;
	x_tau_err = 0;
	x_tau_chi2 = 0;
	x_idx_trk0 = 0;
	x_idx_trk1 = 0;
	x_idx_trk2 = 0;
	x_idx_trk3 = 0;
	x_idx_trk4 = 0;
	x_idx_trk5 = 0;
	x_idx_neu0 = 0;
	x_idx_neu1 = 0;
	x_idx_neu2 = 0;
	x_idx_neu3 = 0;
	x_idx_neu4 = 0;
	x_idx_neu5 = 0;
	x_l0_hadron_dec = 0;
	x_l0_hadron_tis = 0;
	x_l0_hadron_tos = 0;
	x_l0_photon_dec = 0;
	x_l0_photon_tis = 0;
	x_l0_photon_tos = 0;
	x_l0_electron_dec = 0;
	x_l0_electron_tis = 0;
	x_l0_electron_tos = 0;
	x_l0_muon_dec = 0;
	x_l0_muon_tis = 0;
	x_l0_muon_tos = 0;
	x_l0_dimuon_dec = 0;
	x_l0_dimuon_tis = 0;
	x_l0_dimuon_tos = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

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
	fChain->SetBranchAddress("trk_ip_best", &trk_ip_best, &b_trk_ip_best);
	fChain->SetBranchAddress("trk_ip_bvtx", &trk_ip_bvtx, &b_trk_ip_bvtx);
	fChain->SetBranchAddress("trk_ip_chi2", &trk_ip_chi2, &b_trk_ip_chi2);
	fChain->SetBranchAddress("trk_ip_chi2_best", &trk_ip_chi2_best, &b_trk_ip_chi2_best);
	fChain->SetBranchAddress("trk_ip_chi2_bvtx", &trk_ip_chi2_bvtx, &b_trk_ip_chi2_bvtx);
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
	fChain->SetBranchAddress("trk_l0_photon_dec", &trk_l0_photon_dec, &b_trk_l0_photon_dec);
	fChain->SetBranchAddress("trk_l0_photon_tis", &trk_l0_photon_tis, &b_trk_l0_photon_tis);
	fChain->SetBranchAddress("trk_l0_photon_tos", &trk_l0_photon_tos, &b_trk_l0_photon_tos);
	fChain->SetBranchAddress("trk_l0_electron_dec", &trk_l0_electron_dec, &b_trk_l0_electron_dec);
	fChain->SetBranchAddress("trk_l0_electron_tis", &trk_l0_electron_tis, &b_trk_l0_electron_tis);
	fChain->SetBranchAddress("trk_l0_electron_tos", &trk_l0_electron_tos, &b_trk_l0_electron_tos);
	fChain->SetBranchAddress("trk_l0_muon_dec", &trk_l0_muon_dec, &b_trk_l0_muon_dec);
	fChain->SetBranchAddress("trk_l0_muon_tis", &trk_l0_muon_tis, &b_trk_l0_muon_tis);
	fChain->SetBranchAddress("trk_l0_muon_tos", &trk_l0_muon_tos, &b_trk_l0_muon_tos);
	fChain->SetBranchAddress("trk_l0_dimuon_dec", &trk_l0_dimuon_dec, &b_trk_l0_dimuon_dec);
	fChain->SetBranchAddress("trk_l0_dimuon_tis", &trk_l0_dimuon_tis, &b_trk_l0_dimuon_tis);
	fChain->SetBranchAddress("trk_l0_dimuon_tos", &trk_l0_dimuon_tos, &b_trk_l0_dimuon_tos);
	fChain->SetBranchAddress("neu_idx_gen", &neu_idx_gen, &b_neu_idx_gen);
	fChain->SetBranchAddress("neu_p", &neu_p, &b_neu_p);
	fChain->SetBranchAddress("neu_pt", &neu_pt, &b_neu_pt);
	fChain->SetBranchAddress("neu_px", &neu_px, &b_neu_px);
	fChain->SetBranchAddress("neu_py", &neu_py, &b_neu_py);
	fChain->SetBranchAddress("neu_pz", &neu_pz, &b_neu_pz);
	fChain->SetBranchAddress("neu_e", &neu_e, &b_neu_e);
	fChain->SetBranchAddress("neu_pid", &neu_pid, &b_neu_pid);
	fChain->SetBranchAddress("neu_cl", &neu_cl, &b_neu_cl);
	fChain->SetBranchAddress("neu_l0_hadron_dec", &neu_l0_hadron_dec, &b_neu_l0_hadron_dec);
	fChain->SetBranchAddress("neu_l0_hadron_tis", &neu_l0_hadron_tis, &b_neu_l0_hadron_tis);
	fChain->SetBranchAddress("neu_l0_hadron_tos", &neu_l0_hadron_tos, &b_neu_l0_hadron_tos);
	fChain->SetBranchAddress("neu_l0_photon_dec", &neu_l0_photon_dec, &b_neu_l0_photon_dec);
	fChain->SetBranchAddress("neu_l0_photon_tis", &neu_l0_photon_tis, &b_neu_l0_photon_tis);
	fChain->SetBranchAddress("neu_l0_photon_tos", &neu_l0_photon_tos, &b_neu_l0_photon_tos);
	fChain->SetBranchAddress("neu_l0_electron_dec", &neu_l0_electron_dec, &b_neu_l0_electron_dec);
	fChain->SetBranchAddress("neu_l0_electron_tis", &neu_l0_electron_tis, &b_neu_l0_electron_tis);
	fChain->SetBranchAddress("neu_l0_electron_tos", &neu_l0_electron_tos, &b_neu_l0_electron_tos);
	fChain->SetBranchAddress("neu_l0_muon_dec", &neu_l0_muon_dec, &b_neu_l0_muon_dec);
	fChain->SetBranchAddress("neu_l0_muon_tis", &neu_l0_muon_tis, &b_neu_l0_muon_tis);
	fChain->SetBranchAddress("neu_l0_muon_tos", &neu_l0_muon_tos, &b_neu_l0_muon_tos);
	fChain->SetBranchAddress("neu_l0_dimuon_dec", &neu_l0_dimuon_dec, &b_neu_l0_dimuon_dec);
	fChain->SetBranchAddress("neu_l0_dimuon_tis", &neu_l0_dimuon_tis, &b_neu_l0_dimuon_tis);
	fChain->SetBranchAddress("neu_l0_dimuon_tos", &neu_l0_dimuon_tos, &b_neu_l0_dimuon_tos);
	fChain->SetBranchAddress("b_idx_pvr", &b_idx_pvr, &b_b_idx_pvr);
	fChain->SetBranchAddress("b_p", &b_p, &b_b_p);
	fChain->SetBranchAddress("b_pt", &b_pt, &b_b_pt);
	fChain->SetBranchAddress("b_px", &b_px, &b_b_px);
	fChain->SetBranchAddress("b_py", &b_py, &b_b_py);
	fChain->SetBranchAddress("b_pz", &b_pz, &b_b_pz);
	fChain->SetBranchAddress("b_e", &b_e, &b_b_e);
	fChain->SetBranchAddress("b_x", &b_x, &b_b_x);
	fChain->SetBranchAddress("b_y", &b_y, &b_b_y);
	fChain->SetBranchAddress("b_z", &b_z, &b_b_z);
	fChain->SetBranchAddress("b_m", &b_m, &b_b_m);
	fChain->SetBranchAddress("b_vtx_chi2", &b_vtx_chi2, &b_b_vtx_chi2);
	fChain->SetBranchAddress("b_vtx_ndof", &b_vtx_ndof, &b_b_vtx_ndof);
	fChain->SetBranchAddress("b_ip", &b_ip, &b_b_ip);
	fChain->SetBranchAddress("b_ip_chi2", &b_ip_chi2, &b_b_ip_chi2);
	fChain->SetBranchAddress("b_fd", &b_fd, &b_b_fd);
	fChain->SetBranchAddress("b_fd_chi2", &b_fd_chi2, &b_b_fd_chi2);
	fChain->SetBranchAddress("b_tau", &b_tau, &b_b_tau);
	fChain->SetBranchAddress("b_tau_err", &b_tau_err, &b_b_tau_err);
	fChain->SetBranchAddress("b_tau_chi2", &b_tau_chi2, &b_b_tau_chi2);
	fChain->SetBranchAddress("b_idx_kst", &b_idx_kst, &b_b_idx_kst);
	fChain->SetBranchAddress("b_idx_k", &b_idx_k, &b_b_idx_k);
	fChain->SetBranchAddress("b_idx_pi", &b_idx_pi, &b_b_idx_pi);
	fChain->SetBranchAddress("b_idx_x", &b_idx_x, &b_b_idx_x);
	fChain->SetBranchAddress("b_l0_hadron_dec", &b_l0_hadron_dec, &b_b_l0_hadron_dec);
	fChain->SetBranchAddress("b_l0_hadron_tis", &b_l0_hadron_tis, &b_b_l0_hadron_tis);
	fChain->SetBranchAddress("b_l0_hadron_tos", &b_l0_hadron_tos, &b_b_l0_hadron_tos);
	fChain->SetBranchAddress("b_l0_photon_dec", &b_l0_photon_dec, &b_b_l0_photon_dec);
	fChain->SetBranchAddress("b_l0_photon_tis", &b_l0_photon_tis, &b_b_l0_photon_tis);
	fChain->SetBranchAddress("b_l0_photon_tos", &b_l0_photon_tos, &b_b_l0_photon_tos);
	fChain->SetBranchAddress("b_l0_electron_dec", &b_l0_electron_dec, &b_b_l0_electron_dec);
	fChain->SetBranchAddress("b_l0_electron_tis", &b_l0_electron_tis, &b_b_l0_electron_tis);
	fChain->SetBranchAddress("b_l0_electron_tos", &b_l0_electron_tos, &b_b_l0_electron_tos);
	fChain->SetBranchAddress("b_l0_muon_dec", &b_l0_muon_dec, &b_b_l0_muon_dec);
	fChain->SetBranchAddress("b_l0_muon_tis", &b_l0_muon_tis, &b_b_l0_muon_tis);
	fChain->SetBranchAddress("b_l0_muon_tos", &b_l0_muon_tos, &b_b_l0_muon_tos);
	fChain->SetBranchAddress("b_l0_dimuon_dec", &b_l0_dimuon_dec, &b_b_l0_dimuon_dec);
	fChain->SetBranchAddress("b_l0_dimuon_tis", &b_l0_dimuon_tis, &b_b_l0_dimuon_tis);
	fChain->SetBranchAddress("b_l0_dimuon_tos", &b_l0_dimuon_tos, &b_b_l0_dimuon_tos);
	fChain->SetBranchAddress("b_mode", &b_mode, &b_b_mode);
	fChain->SetBranchAddress("b_dtf_chi2", &b_dtf_chi2, &b_b_dtf_chi2);
	fChain->SetBranchAddress("b_dtf_ndof", &b_dtf_ndof, &b_b_dtf_ndof);
	fChain->SetBranchAddress("b_dtf_px", &b_dtf_px, &b_b_dtf_px);
	fChain->SetBranchAddress("b_dtf_py", &b_dtf_py, &b_b_dtf_py);
	fChain->SetBranchAddress("b_dtf_pz", &b_dtf_pz, &b_b_dtf_pz);
	fChain->SetBranchAddress("b_dtf_e", &b_dtf_e, &b_b_dtf_e);
	fChain->SetBranchAddress("b_dtf_m", &b_dtf_m, &b_b_dtf_m);
	fChain->SetBranchAddress("b_dtf_ctau", &b_dtf_ctau, &b_b_dtf_ctau);
	fChain->SetBranchAddress("b_dtf_x_m", &b_dtf_x_m, &b_b_dtf_x_m);
	fChain->SetBranchAddress("b_dtf_x_ctau", &b_dtf_x_ctau, &b_b_dtf_x_ctau);
	fChain->SetBranchAddress("b_dtf_trk0_px", &b_dtf_trk0_px, &b_b_dtf_trk0_px);
	fChain->SetBranchAddress("b_dtf_trk0_py", &b_dtf_trk0_py, &b_b_dtf_trk0_py);
	fChain->SetBranchAddress("b_dtf_trk0_pz", &b_dtf_trk0_pz, &b_b_dtf_trk0_pz);
	fChain->SetBranchAddress("b_dtf_trk0_e", &b_dtf_trk0_e, &b_b_dtf_trk0_e);
	fChain->SetBranchAddress("b_dtf_trk0_m", &b_dtf_trk0_m, &b_b_dtf_trk0_m);
	fChain->SetBranchAddress("b_dtf_trk1_px", &b_dtf_trk1_px, &b_b_dtf_trk1_px);
	fChain->SetBranchAddress("b_dtf_trk1_py", &b_dtf_trk1_py, &b_b_dtf_trk1_py);
	fChain->SetBranchAddress("b_dtf_trk1_pz", &b_dtf_trk1_pz, &b_b_dtf_trk1_pz);
	fChain->SetBranchAddress("b_dtf_trk1_e", &b_dtf_trk1_e, &b_b_dtf_trk1_e);
	fChain->SetBranchAddress("b_dtf_trk1_m", &b_dtf_trk1_m, &b_b_dtf_trk1_m);
	fChain->SetBranchAddress("b_dtf_trk2_px", &b_dtf_trk2_px, &b_b_dtf_trk2_px);
	fChain->SetBranchAddress("b_dtf_trk2_py", &b_dtf_trk2_py, &b_b_dtf_trk2_py);
	fChain->SetBranchAddress("b_dtf_trk2_pz", &b_dtf_trk2_pz, &b_b_dtf_trk2_pz);
	fChain->SetBranchAddress("b_dtf_trk2_e", &b_dtf_trk2_e, &b_b_dtf_trk2_e);
	fChain->SetBranchAddress("b_dtf_trk2_m", &b_dtf_trk2_m, &b_b_dtf_trk2_m);
	fChain->SetBranchAddress("b_dtf_trk3_px", &b_dtf_trk3_px, &b_b_dtf_trk3_px);
	fChain->SetBranchAddress("b_dtf_trk3_py", &b_dtf_trk3_py, &b_b_dtf_trk3_py);
	fChain->SetBranchAddress("b_dtf_trk3_pz", &b_dtf_trk3_pz, &b_b_dtf_trk3_pz);
	fChain->SetBranchAddress("b_dtf_trk3_e", &b_dtf_trk3_e, &b_b_dtf_trk3_e);
	fChain->SetBranchAddress("b_dtf_trk3_m", &b_dtf_trk3_m, &b_b_dtf_trk3_m);
	fChain->SetBranchAddress("b_dtf_trk4_px", &b_dtf_trk4_px, &b_b_dtf_trk4_px);
	fChain->SetBranchAddress("b_dtf_trk4_py", &b_dtf_trk4_py, &b_b_dtf_trk4_py);
	fChain->SetBranchAddress("b_dtf_trk4_pz", &b_dtf_trk4_pz, &b_b_dtf_trk4_pz);
	fChain->SetBranchAddress("b_dtf_trk4_e", &b_dtf_trk4_e, &b_b_dtf_trk4_e);
	fChain->SetBranchAddress("b_dtf_trk4_m", &b_dtf_trk4_m, &b_b_dtf_trk4_m);
	fChain->SetBranchAddress("b_dtf_trk5_px", &b_dtf_trk5_px, &b_b_dtf_trk5_px);
	fChain->SetBranchAddress("b_dtf_trk5_py", &b_dtf_trk5_py, &b_b_dtf_trk5_py);
	fChain->SetBranchAddress("b_dtf_trk5_pz", &b_dtf_trk5_pz, &b_b_dtf_trk5_pz);
	fChain->SetBranchAddress("b_dtf_trk5_e", &b_dtf_trk5_e, &b_b_dtf_trk5_e);
	fChain->SetBranchAddress("b_dtf_trk5_m", &b_dtf_trk5_m, &b_b_dtf_trk5_m);
	fChain->SetBranchAddress("b_dtf_trk6_px", &b_dtf_trk6_px, &b_b_dtf_trk6_px);
	fChain->SetBranchAddress("b_dtf_trk6_py", &b_dtf_trk6_py, &b_b_dtf_trk6_py);
	fChain->SetBranchAddress("b_dtf_trk6_pz", &b_dtf_trk6_pz, &b_b_dtf_trk6_pz);
	fChain->SetBranchAddress("b_dtf_trk6_e", &b_dtf_trk6_e, &b_b_dtf_trk6_e);
	fChain->SetBranchAddress("b_dtf_trk6_m", &b_dtf_trk6_m, &b_b_dtf_trk6_m);
	fChain->SetBranchAddress("b_dtf_trk7_px", &b_dtf_trk7_px, &b_b_dtf_trk7_px);
	fChain->SetBranchAddress("b_dtf_trk7_py", &b_dtf_trk7_py, &b_b_dtf_trk7_py);
	fChain->SetBranchAddress("b_dtf_trk7_pz", &b_dtf_trk7_pz, &b_b_dtf_trk7_pz);
	fChain->SetBranchAddress("b_dtf_trk7_e", &b_dtf_trk7_e, &b_b_dtf_trk7_e);
	fChain->SetBranchAddress("b_dtf_trk7_m", &b_dtf_trk7_m, &b_b_dtf_trk7_m);
	fChain->SetBranchAddress("b_dtf_neu0_px", &b_dtf_neu0_px, &b_b_dtf_neu0_px);
	fChain->SetBranchAddress("b_dtf_neu0_py", &b_dtf_neu0_py, &b_b_dtf_neu0_py);
	fChain->SetBranchAddress("b_dtf_neu0_pz", &b_dtf_neu0_pz, &b_b_dtf_neu0_pz);
	fChain->SetBranchAddress("b_dtf_neu0_e", &b_dtf_neu0_e, &b_b_dtf_neu0_e);
	fChain->SetBranchAddress("b_dtf_neu0_m", &b_dtf_neu0_m, &b_b_dtf_neu0_m);
	fChain->SetBranchAddress("b_dtf_neu1_px", &b_dtf_neu1_px, &b_b_dtf_neu1_px);
	fChain->SetBranchAddress("b_dtf_neu1_py", &b_dtf_neu1_py, &b_b_dtf_neu1_py);
	fChain->SetBranchAddress("b_dtf_neu1_pz", &b_dtf_neu1_pz, &b_b_dtf_neu1_pz);
	fChain->SetBranchAddress("b_dtf_neu1_e", &b_dtf_neu1_e, &b_b_dtf_neu1_e);
	fChain->SetBranchAddress("b_dtf_neu1_m", &b_dtf_neu1_m, &b_b_dtf_neu1_m);
	fChain->SetBranchAddress("b_dtf_neu2_px", &b_dtf_neu2_px, &b_b_dtf_neu2_px);
	fChain->SetBranchAddress("b_dtf_neu2_py", &b_dtf_neu2_py, &b_b_dtf_neu2_py);
	fChain->SetBranchAddress("b_dtf_neu2_pz", &b_dtf_neu2_pz, &b_b_dtf_neu2_pz);
	fChain->SetBranchAddress("b_dtf_neu2_e", &b_dtf_neu2_e, &b_b_dtf_neu2_e);
	fChain->SetBranchAddress("b_dtf_neu2_m", &b_dtf_neu2_m, &b_b_dtf_neu2_m);
	fChain->SetBranchAddress("b_dtf_neu3_px", &b_dtf_neu3_px, &b_b_dtf_neu3_px);
	fChain->SetBranchAddress("b_dtf_neu3_py", &b_dtf_neu3_py, &b_b_dtf_neu3_py);
	fChain->SetBranchAddress("b_dtf_neu3_pz", &b_dtf_neu3_pz, &b_b_dtf_neu3_pz);
	fChain->SetBranchAddress("b_dtf_neu3_e", &b_dtf_neu3_e, &b_b_dtf_neu3_e);
	fChain->SetBranchAddress("b_dtf_neu3_m", &b_dtf_neu3_m, &b_b_dtf_neu3_m);
	fChain->SetBranchAddress("b_dtf_neu4_px", &b_dtf_neu4_px, &b_b_dtf_neu4_px);
	fChain->SetBranchAddress("b_dtf_neu4_py", &b_dtf_neu4_py, &b_b_dtf_neu4_py);
	fChain->SetBranchAddress("b_dtf_neu4_pz", &b_dtf_neu4_pz, &b_b_dtf_neu4_pz);
	fChain->SetBranchAddress("b_dtf_neu4_e", &b_dtf_neu4_e, &b_b_dtf_neu4_e);
	fChain->SetBranchAddress("b_dtf_neu4_m", &b_dtf_neu4_m, &b_b_dtf_neu4_m);
	fChain->SetBranchAddress("b_dtf_neu5_px", &b_dtf_neu5_px, &b_b_dtf_neu5_px);
	fChain->SetBranchAddress("b_dtf_neu5_py", &b_dtf_neu5_py, &b_b_dtf_neu5_py);
	fChain->SetBranchAddress("b_dtf_neu5_pz", &b_dtf_neu5_pz, &b_b_dtf_neu5_pz);
	fChain->SetBranchAddress("b_dtf_neu5_e", &b_dtf_neu5_e, &b_b_dtf_neu5_e);
	fChain->SetBranchAddress("b_dtf_neu5_m", &b_dtf_neu5_m, &b_b_dtf_neu5_m);
	fChain->SetBranchAddress("b_dtf_neu6_px", &b_dtf_neu6_px, &b_b_dtf_neu6_px);
	fChain->SetBranchAddress("b_dtf_neu6_py", &b_dtf_neu6_py, &b_b_dtf_neu6_py);
	fChain->SetBranchAddress("b_dtf_neu6_pz", &b_dtf_neu6_pz, &b_b_dtf_neu6_pz);
	fChain->SetBranchAddress("b_dtf_neu6_e", &b_dtf_neu6_e, &b_b_dtf_neu6_e);
	fChain->SetBranchAddress("b_dtf_neu6_m", &b_dtf_neu6_m, &b_b_dtf_neu6_m);
	fChain->SetBranchAddress("b_dtf_neu7_px", &b_dtf_neu7_px, &b_b_dtf_neu7_px);
	fChain->SetBranchAddress("b_dtf_neu7_py", &b_dtf_neu7_py, &b_b_dtf_neu7_py);
	fChain->SetBranchAddress("b_dtf_neu7_pz", &b_dtf_neu7_pz, &b_b_dtf_neu7_pz);
	fChain->SetBranchAddress("b_dtf_neu7_e", &b_dtf_neu7_e, &b_b_dtf_neu7_e);
	fChain->SetBranchAddress("b_dtf_neu7_m", &b_dtf_neu7_m, &b_b_dtf_neu7_m);
	fChain->SetBranchAddress("kst_idx_pvr", &kst_idx_pvr, &b_kst_idx_pvr);
	fChain->SetBranchAddress("kst_p", &kst_p, &b_kst_p);
	fChain->SetBranchAddress("kst_pt", &kst_pt, &b_kst_pt);
	fChain->SetBranchAddress("kst_px", &kst_px, &b_kst_px);
	fChain->SetBranchAddress("kst_py", &kst_py, &b_kst_py);
	fChain->SetBranchAddress("kst_pz", &kst_pz, &b_kst_pz);
	fChain->SetBranchAddress("kst_e", &kst_e, &b_kst_e);
	fChain->SetBranchAddress("kst_x", &kst_x, &b_kst_x);
	fChain->SetBranchAddress("kst_y", &kst_y, &b_kst_y);
	fChain->SetBranchAddress("kst_z", &kst_z, &b_kst_z);
	fChain->SetBranchAddress("kst_m", &kst_m, &b_kst_m);
	fChain->SetBranchAddress("kst_vtx_chi2", &kst_vtx_chi2, &b_kst_vtx_chi2);
	fChain->SetBranchAddress("kst_vtx_ndof", &kst_vtx_ndof, &b_kst_vtx_ndof);
	fChain->SetBranchAddress("kst_ip", &kst_ip, &b_kst_ip);
	fChain->SetBranchAddress("kst_ip_best", &kst_ip_best, &b_kst_ip_best);
	fChain->SetBranchAddress("kst_ip_bvtx", &kst_ip_bvtx, &b_kst_ip_bvtx);
	fChain->SetBranchAddress("kst_ip_chi2", &kst_ip_chi2, &b_kst_ip_chi2);
	fChain->SetBranchAddress("kst_ip_chi2_best", &kst_ip_chi2_best, &b_kst_ip_chi2_best);
	fChain->SetBranchAddress("kst_ip_chi2_bvtx", &kst_ip_chi2_bvtx, &b_kst_ip_chi2_bvtx);
	fChain->SetBranchAddress("kst_fd", &kst_fd, &b_kst_fd);
	fChain->SetBranchAddress("kst_fd_best", &kst_fd_best, &b_kst_fd_best);
	fChain->SetBranchAddress("kst_fd_bvtx", &kst_fd_bvtx, &b_kst_fd_bvtx);
	fChain->SetBranchAddress("kst_fd_chi2", &kst_fd_chi2, &b_kst_fd_chi2);
	fChain->SetBranchAddress("kst_fd_chi2_best", &kst_fd_chi2_best, &b_kst_fd_chi2_best);
	fChain->SetBranchAddress("kst_fd_chi2_bvtx", &kst_fd_chi2_bvtx, &b_kst_fd_chi2_bvtx);
	fChain->SetBranchAddress("kst_tau", &kst_tau, &b_kst_tau);
	fChain->SetBranchAddress("kst_tau_err", &kst_tau_err, &b_kst_tau_err);
	fChain->SetBranchAddress("kst_tau_chi2", &kst_tau_chi2, &b_kst_tau_chi2);
	fChain->SetBranchAddress("kst_idx_k", &kst_idx_k, &b_kst_idx_k);
	fChain->SetBranchAddress("kst_idx_pi", &kst_idx_pi, &b_kst_idx_pi);
	fChain->SetBranchAddress("kst_l0_hadron_dec", &kst_l0_hadron_dec, &b_kst_l0_hadron_dec);
	fChain->SetBranchAddress("kst_l0_hadron_tis", &kst_l0_hadron_tis, &b_kst_l0_hadron_tis);
	fChain->SetBranchAddress("kst_l0_hadron_tos", &kst_l0_hadron_tos, &b_kst_l0_hadron_tos);
	fChain->SetBranchAddress("kst_l0_photon_dec", &kst_l0_photon_dec, &b_kst_l0_photon_dec);
	fChain->SetBranchAddress("kst_l0_photon_tis", &kst_l0_photon_tis, &b_kst_l0_photon_tis);
	fChain->SetBranchAddress("kst_l0_photon_tos", &kst_l0_photon_tos, &b_kst_l0_photon_tos);
	fChain->SetBranchAddress("kst_l0_electron_dec", &kst_l0_electron_dec, &b_kst_l0_electron_dec);
	fChain->SetBranchAddress("kst_l0_electron_tis", &kst_l0_electron_tis, &b_kst_l0_electron_tis);
	fChain->SetBranchAddress("kst_l0_electron_tos", &kst_l0_electron_tos, &b_kst_l0_electron_tos);
	fChain->SetBranchAddress("kst_l0_muon_dec", &kst_l0_muon_dec, &b_kst_l0_muon_dec);
	fChain->SetBranchAddress("kst_l0_muon_tis", &kst_l0_muon_tis, &b_kst_l0_muon_tis);
	fChain->SetBranchAddress("kst_l0_muon_tos", &kst_l0_muon_tos, &b_kst_l0_muon_tos);
	fChain->SetBranchAddress("kst_l0_dimuon_dec", &kst_l0_dimuon_dec, &b_kst_l0_dimuon_dec);
	fChain->SetBranchAddress("kst_l0_dimuon_tis", &kst_l0_dimuon_tis, &b_kst_l0_dimuon_tis);
	fChain->SetBranchAddress("kst_l0_dimuon_tos", &kst_l0_dimuon_tos, &b_kst_l0_dimuon_tos);
	fChain->SetBranchAddress("x_idx_pvr", &x_idx_pvr, &b_x_idx_pvr);
	fChain->SetBranchAddress("x_p", &x_p, &b_x_p);
	fChain->SetBranchAddress("x_pt", &x_pt, &b_x_pt);
	fChain->SetBranchAddress("x_px", &x_px, &b_x_px);
	fChain->SetBranchAddress("x_py", &x_py, &b_x_py);
	fChain->SetBranchAddress("x_pz", &x_pz, &b_x_pz);
	fChain->SetBranchAddress("x_e", &x_e, &b_x_e);
	fChain->SetBranchAddress("x_x", &x_x, &b_x_x);
	fChain->SetBranchAddress("x_y", &x_y, &b_x_y);
	fChain->SetBranchAddress("x_z", &x_z, &b_x_z);
	fChain->SetBranchAddress("x_m", &x_m, &b_x_m);
	fChain->SetBranchAddress("x_vtx_chi2", &x_vtx_chi2, &b_x_vtx_chi2);
	fChain->SetBranchAddress("x_vtx_ndof", &x_vtx_ndof, &b_x_vtx_ndof);
	fChain->SetBranchAddress("x_ip", &x_ip, &b_x_ip);
	fChain->SetBranchAddress("x_ip_best", &x_ip_best, &b_x_ip_best);
	fChain->SetBranchAddress("x_ip_bvtx", &x_ip_bvtx, &b_x_ip_bvtx);
	fChain->SetBranchAddress("x_ip_chi2", &x_ip_chi2, &b_x_ip_chi2);
	fChain->SetBranchAddress("x_ip_chi2_best", &x_ip_chi2_best, &b_x_ip_chi2_best);
	fChain->SetBranchAddress("x_ip_chi2_bvtx", &x_ip_chi2_bvtx, &b_x_ip_chi2_bvtx);
	fChain->SetBranchAddress("x_fd", &x_fd, &b_x_fd);
	fChain->SetBranchAddress("x_fd_best", &x_fd_best, &b_x_fd_best);
	fChain->SetBranchAddress("x_fd_bvtx", &x_fd_bvtx, &b_x_fd_bvtx);
	fChain->SetBranchAddress("x_fd_chi2", &x_fd_chi2, &b_x_fd_chi2);
	fChain->SetBranchAddress("x_fd_chi2_best", &x_fd_chi2_best, &b_x_fd_chi2_best);
	fChain->SetBranchAddress("x_fd_chi2_bvtx", &x_fd_chi2_bvtx, &b_x_fd_chi2_bvtx);
	fChain->SetBranchAddress("x_tau", &x_tau, &b_x_tau);
	fChain->SetBranchAddress("x_tau_err", &x_tau_err, &b_x_tau_err);
	fChain->SetBranchAddress("x_tau_chi2", &x_tau_chi2, &b_x_tau_chi2);
	fChain->SetBranchAddress("x_idx_trk0", &x_idx_trk0, &b_x_idx_trk0);
	fChain->SetBranchAddress("x_idx_trk1", &x_idx_trk1, &b_x_idx_trk1);
	fChain->SetBranchAddress("x_idx_trk2", &x_idx_trk2, &b_x_idx_trk2);
	fChain->SetBranchAddress("x_idx_trk3", &x_idx_trk3, &b_x_idx_trk3);
	fChain->SetBranchAddress("x_idx_trk4", &x_idx_trk4, &b_x_idx_trk4);
	fChain->SetBranchAddress("x_idx_trk5", &x_idx_trk5, &b_x_idx_trk5);
	fChain->SetBranchAddress("x_idx_neu0", &x_idx_neu0, &b_x_idx_neu0);
	fChain->SetBranchAddress("x_idx_neu1", &x_idx_neu1, &b_x_idx_neu1);
	fChain->SetBranchAddress("x_idx_neu2", &x_idx_neu2, &b_x_idx_neu2);
	fChain->SetBranchAddress("x_idx_neu3", &x_idx_neu3, &b_x_idx_neu3);
	fChain->SetBranchAddress("x_idx_neu4", &x_idx_neu4, &b_x_idx_neu4);
	fChain->SetBranchAddress("x_idx_neu5", &x_idx_neu5, &b_x_idx_neu5);
	fChain->SetBranchAddress("x_l0_hadron_dec", &x_l0_hadron_dec, &b_x_l0_hadron_dec);
	fChain->SetBranchAddress("x_l0_hadron_tis", &x_l0_hadron_tis, &b_x_l0_hadron_tis);
	fChain->SetBranchAddress("x_l0_hadron_tos", &x_l0_hadron_tos, &b_x_l0_hadron_tos);
	fChain->SetBranchAddress("x_l0_photon_dec", &x_l0_photon_dec, &b_x_l0_photon_dec);
	fChain->SetBranchAddress("x_l0_photon_tis", &x_l0_photon_tis, &b_x_l0_photon_tis);
	fChain->SetBranchAddress("x_l0_photon_tos", &x_l0_photon_tos, &b_x_l0_photon_tos);
	fChain->SetBranchAddress("x_l0_electron_dec", &x_l0_electron_dec, &b_x_l0_electron_dec);
	fChain->SetBranchAddress("x_l0_electron_tis", &x_l0_electron_tis, &b_x_l0_electron_tis);
	fChain->SetBranchAddress("x_l0_electron_tos", &x_l0_electron_tos, &b_x_l0_electron_tos);
	fChain->SetBranchAddress("x_l0_muon_dec", &x_l0_muon_dec, &b_x_l0_muon_dec);
	fChain->SetBranchAddress("x_l0_muon_tis", &x_l0_muon_tis, &b_x_l0_muon_tis);
	fChain->SetBranchAddress("x_l0_muon_tos", &x_l0_muon_tos, &b_x_l0_muon_tos);
	fChain->SetBranchAddress("x_l0_dimuon_dec", &x_l0_dimuon_dec, &b_x_l0_dimuon_dec);
	fChain->SetBranchAddress("x_l0_dimuon_tis", &x_l0_dimuon_tis, &b_x_l0_dimuon_tis);
	fChain->SetBranchAddress("x_l0_dimuon_tos", &x_l0_dimuon_tos, &b_x_l0_dimuon_tos);
	fChain->SetBranchAddress("evt_pvr_n", &evt_pvr_n, &b_evt_pvr_n);
	fChain->SetBranchAddress("evt_trk_n", &evt_trk_n, &b_evt_trk_n);
	Notify();
}

Bool_t TreeSkimmer::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void TreeSkimmer::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t TreeSkimmer::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef TreeSkimmer_cxx
