//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  7 12:11:21 2017 by ROOT version 6.06/00
// from TTree data/data
// found on file: lightjets.part.root
//////////////////////////////////////////////////////////

#ifndef resampleTrackIPs_h
#define resampleTrackIPs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class resampleTrackIPs {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<double>  *gen_idx_pvr;
   std::vector<double>  *gen_idx_jet;
   std::vector<double>  *gen_pid;
   std::vector<double>  *gen_q;
   std::vector<double>  *gen_px;
   std::vector<double>  *gen_py;
   std::vector<double>  *gen_pz;
   std::vector<double>  *gen_e;
   std::vector<double>  *gen_x;
   std::vector<double>  *gen_y;
   std::vector<double>  *gen_z;
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
   std::vector<double>  *jet_px;
   std::vector<double>  *jet_py;
   std::vector<double>  *jet_pz;
   std::vector<double>  *jet_e;
   std::vector<double>  *jet_idx_trk0;
   std::vector<double>  *jet_idx_trk1;
   std::vector<double>  *jet_idx_trk2;
   std::vector<double>  *jet_idx_trk3;
   std::vector<double>  *jet_idx_trk4;
   std::vector<double>  *jet_idx_trk5;
   std::vector<double>  *jet_idx_trk6;
   std::vector<double>  *jet_idx_trk7;
   std::vector<double>  *jet_idx_trk8;
   std::vector<double>  *jet_idx_trk9;
   std::vector<double>  *jet_idx_trk10;
   std::vector<double>  *jet_idx_trk11;
   std::vector<double>  *jet_idx_trk12;
   std::vector<double>  *jet_idx_trk13;
   std::vector<double>  *jet_idx_trk14;
   std::vector<double>  *jet_idx_trk15;
   std::vector<double>  *jet_idx_trk16;
   std::vector<double>  *jet_idx_trk17;
   std::vector<double>  *jet_idx_trk18;
   std::vector<double>  *jet_idx_trk19;
   std::vector<double>  *jet_idx_trk20;
   std::vector<double>  *jet_idx_trk21;
   std::vector<double>  *jet_idx_trk22;
   std::vector<double>  *jet_idx_trk23;
   std::vector<double>  *jet_idx_trk24;
   std::vector<double>  *jet_idx_trk25;
   std::vector<double>  *jet_idx_trk26;
   std::vector<double>  *jet_idx_trk27;
   std::vector<double>  *jet_idx_trk28;
   std::vector<double>  *jet_idx_trk29;
   std::vector<double>  *jet_idx_trk30;
   std::vector<double>  *jet_idx_trk31;
   std::vector<double>  *jet_idx_trk32;
   std::vector<double>  *jet_idx_trk33;
   std::vector<double>  *jet_idx_trk34;
   std::vector<double>  *jet_idx_trk35;
   std::vector<double>  *jet_idx_trk36;
   std::vector<double>  *jet_idx_trk37;
   std::vector<double>  *jet_idx_trk38;
   std::vector<double>  *jet_idx_trk39;
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
   std::vector<double>  *trk_ecal;
   std::vector<double>  *trk_hcal;
   std::vector<double>  *trk_prb_ghost;
   std::vector<double>  *trk_type;
   std::vector<double>  *trk_is_mu;
   std::vector<double>  *trk_vid;
   std::vector<double>  *trk_x;
   std::vector<double>  *trk_y;
   std::vector<double>  *trk_z;
   std::vector<double>  *trk_vhit0;
   std::vector<double>  *trk_vhit1;
   std::vector<double>  *trk_vhit2;
   std::vector<double>  *trk_vhit3;
   std::vector<double>  *trk_vhit4;
   std::vector<double>  *trk_vhit5;
   std::vector<double>  *trk_vhit6;
   std::vector<double>  *trk_vhit7;
   std::vector<double>  *trk_vhit8;
   std::vector<double>  *trk_vhit9;
   std::vector<double>  *trk_vhit10;
   std::vector<double>  *trk_vhit11;
   std::vector<double>  *trk_vhit12;
   std::vector<double>  *trk_vhit13;
   std::vector<double>  *trk_vhit14;
   std::vector<double>  *trk_vhit15;
   std::vector<double>  *trk_vhit16;
   std::vector<double>  *trk_vhit17;
   std::vector<double>  *trk_vhit18;
   std::vector<double>  *trk_vhit19;
   std::vector<double>  *trk_vhit20;
   std::vector<double>  *trk_vhit21;
   std::vector<double>  *trk_vhit22;
   std::vector<double>  *trk_vhit23;
   std::vector<double>  *trk_vhit24;
   std::vector<double>  *trk_vhit25;
   std::vector<double>  *trk_vhit26;
   std::vector<double>  *trk_vhit27;
   std::vector<double>  *trk_vhit28;
   std::vector<double>  *trk_vhit29;
   std::vector<double>  *trk_vhit30;
   std::vector<double>  *trk_vhit31;
   std::vector<double>  *trk_vhit32;
   std::vector<double>  *trk_vhit33;
   std::vector<double>  *trk_vhit34;
   std::vector<double>  *trk_vhit35;
   std::vector<double>  *trk_vhit36;
   std::vector<double>  *trk_vhit37;
   std::vector<double>  *trk_vhit38;
   std::vector<double>  *trk_vhit39;
   std::vector<double>  *trk_vhit40;
   std::vector<double>  *trk_vhit41;
   std::vector<double>  *trk_vhit42;
   std::vector<double>  *trk_vhit43;
   std::vector<double>  *trk_vhit44;
   std::vector<double>  *trk_vhit45;
   std::vector<double>  *trk_vhit46;
   std::vector<double>  *trk_vhit47;
   std::vector<double>  *trk_vhit48;
   std::vector<double>  *trk_vhit49;
   std::vector<double>  *trk_vhit50;
   std::vector<double>  *trk_vhit51;
   std::vector<double>  *trk_vhit52;
   std::vector<double>  *trk_vhit53;
   std::vector<double>  *trk_vhit54;
   std::vector<double>  *trk_vhit55;
   std::vector<double>  *trk_vhit56;
   std::vector<double>  *trk_vhit57;
   std::vector<double>  *trk_vhit58;
   std::vector<double>  *trk_vhit59;
   std::vector<double>  *trk_vhit60;
   std::vector<double>  *neu_idx_gen;
   std::vector<double>  *neu_idx_jet;
   std::vector<double>  *neu_px;
   std::vector<double>  *neu_py;
   std::vector<double>  *neu_pz;
   std::vector<double>  *neu_e;
   std::vector<double>  *neu_pid;
   Double_t        evt_pvr_n;

   // List of branches
   TBranch        *b_gen_idx_pvr;   //!
   TBranch        *b_gen_idx_jet;   //!
   TBranch        *b_gen_pid;   //!
   TBranch        *b_gen_q;   //!
   TBranch        *b_gen_px;   //!
   TBranch        *b_gen_py;   //!
   TBranch        *b_gen_pz;   //!
   TBranch        *b_gen_e;   //!
   TBranch        *b_gen_x;   //!
   TBranch        *b_gen_y;   //!
   TBranch        *b_gen_z;   //!
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
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_idx_trk0;   //!
   TBranch        *b_jet_idx_trk1;   //!
   TBranch        *b_jet_idx_trk2;   //!
   TBranch        *b_jet_idx_trk3;   //!
   TBranch        *b_jet_idx_trk4;   //!
   TBranch        *b_jet_idx_trk5;   //!
   TBranch        *b_jet_idx_trk6;   //!
   TBranch        *b_jet_idx_trk7;   //!
   TBranch        *b_jet_idx_trk8;   //!
   TBranch        *b_jet_idx_trk9;   //!
   TBranch        *b_jet_idx_trk10;   //!
   TBranch        *b_jet_idx_trk11;   //!
   TBranch        *b_jet_idx_trk12;   //!
   TBranch        *b_jet_idx_trk13;   //!
   TBranch        *b_jet_idx_trk14;   //!
   TBranch        *b_jet_idx_trk15;   //!
   TBranch        *b_jet_idx_trk16;   //!
   TBranch        *b_jet_idx_trk17;   //!
   TBranch        *b_jet_idx_trk18;   //!
   TBranch        *b_jet_idx_trk19;   //!
   TBranch        *b_jet_idx_trk20;   //!
   TBranch        *b_jet_idx_trk21;   //!
   TBranch        *b_jet_idx_trk22;   //!
   TBranch        *b_jet_idx_trk23;   //!
   TBranch        *b_jet_idx_trk24;   //!
   TBranch        *b_jet_idx_trk25;   //!
   TBranch        *b_jet_idx_trk26;   //!
   TBranch        *b_jet_idx_trk27;   //!
   TBranch        *b_jet_idx_trk28;   //!
   TBranch        *b_jet_idx_trk29;   //!
   TBranch        *b_jet_idx_trk30;   //!
   TBranch        *b_jet_idx_trk31;   //!
   TBranch        *b_jet_idx_trk32;   //!
   TBranch        *b_jet_idx_trk33;   //!
   TBranch        *b_jet_idx_trk34;   //!
   TBranch        *b_jet_idx_trk35;   //!
   TBranch        *b_jet_idx_trk36;   //!
   TBranch        *b_jet_idx_trk37;   //!
   TBranch        *b_jet_idx_trk38;   //!
   TBranch        *b_jet_idx_trk39;   //!
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
   TBranch        *b_trk_ecal;   //!
   TBranch        *b_trk_hcal;   //!
   TBranch        *b_trk_prb_ghost;   //!
   TBranch        *b_trk_type;   //!
   TBranch        *b_trk_is_mu;   //!
   TBranch        *b_trk_vid;   //!
   TBranch        *b_trk_x;   //!
   TBranch        *b_trk_y;   //!
   TBranch        *b_trk_z;   //!
   TBranch        *b_trk_vhit0;   //!
   TBranch        *b_trk_vhit1;   //!
   TBranch        *b_trk_vhit2;   //!
   TBranch        *b_trk_vhit3;   //!
   TBranch        *b_trk_vhit4;   //!
   TBranch        *b_trk_vhit5;   //!
   TBranch        *b_trk_vhit6;   //!
   TBranch        *b_trk_vhit7;   //!
   TBranch        *b_trk_vhit8;   //!
   TBranch        *b_trk_vhit9;   //!
   TBranch        *b_trk_vhit10;   //!
   TBranch        *b_trk_vhit11;   //!
   TBranch        *b_trk_vhit12;   //!
   TBranch        *b_trk_vhit13;   //!
   TBranch        *b_trk_vhit14;   //!
   TBranch        *b_trk_vhit15;   //!
   TBranch        *b_trk_vhit16;   //!
   TBranch        *b_trk_vhit17;   //!
   TBranch        *b_trk_vhit18;   //!
   TBranch        *b_trk_vhit19;   //!
   TBranch        *b_trk_vhit20;   //!
   TBranch        *b_trk_vhit21;   //!
   TBranch        *b_trk_vhit22;   //!
   TBranch        *b_trk_vhit23;   //!
   TBranch        *b_trk_vhit24;   //!
   TBranch        *b_trk_vhit25;   //!
   TBranch        *b_trk_vhit26;   //!
   TBranch        *b_trk_vhit27;   //!
   TBranch        *b_trk_vhit28;   //!
   TBranch        *b_trk_vhit29;   //!
   TBranch        *b_trk_vhit30;   //!
   TBranch        *b_trk_vhit31;   //!
   TBranch        *b_trk_vhit32;   //!
   TBranch        *b_trk_vhit33;   //!
   TBranch        *b_trk_vhit34;   //!
   TBranch        *b_trk_vhit35;   //!
   TBranch        *b_trk_vhit36;   //!
   TBranch        *b_trk_vhit37;   //!
   TBranch        *b_trk_vhit38;   //!
   TBranch        *b_trk_vhit39;   //!
   TBranch        *b_trk_vhit40;   //!
   TBranch        *b_trk_vhit41;   //!
   TBranch        *b_trk_vhit42;   //!
   TBranch        *b_trk_vhit43;   //!
   TBranch        *b_trk_vhit44;   //!
   TBranch        *b_trk_vhit45;   //!
   TBranch        *b_trk_vhit46;   //!
   TBranch        *b_trk_vhit47;   //!
   TBranch        *b_trk_vhit48;   //!
   TBranch        *b_trk_vhit49;   //!
   TBranch        *b_trk_vhit50;   //!
   TBranch        *b_trk_vhit51;   //!
   TBranch        *b_trk_vhit52;   //!
   TBranch        *b_trk_vhit53;   //!
   TBranch        *b_trk_vhit54;   //!
   TBranch        *b_trk_vhit55;   //!
   TBranch        *b_trk_vhit56;   //!
   TBranch        *b_trk_vhit57;   //!
   TBranch        *b_trk_vhit58;   //!
   TBranch        *b_trk_vhit59;   //!
   TBranch        *b_trk_vhit60;   //!
   TBranch        *b_neu_idx_gen;   //!
   TBranch        *b_neu_idx_jet;   //!
   TBranch        *b_neu_px;   //!
   TBranch        *b_neu_py;   //!
   TBranch        *b_neu_pz;   //!
   TBranch        *b_neu_e;   //!
   TBranch        *b_neu_pid;   //!
   TBranch        *b_evt_pvr_n;   //!

   int _irep;
   TFile* newfile;
   TTree* newtree;

   resampleTrackIPs(int irep=-1);
   virtual ~resampleTrackIPs();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   int getNSharedVeloHits(int idxi, int idxj);
   int getNSharedEarlyVeloHits(int idxi, int idxj);
   int findTrkInJet(int trk_idx, int jet_idx);
   int findBestGenPvrForTrk(int idx, double x=-999., double y=-999.);
   int findBestPvrForSvr(TVector3 x, TVector3 p, bool gen=true);
};

#endif

#ifdef resampleTrackIPs_cxx
resampleTrackIPs::resampleTrackIPs(int irep) : fChain(0), _irep(irep), newfile(0), newtree(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   TTree* tree(0);
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lightjets_filtered_addVars.root");
   if (!f || !f->IsOpen()) {
      f = new TFile("lightjets_filtered_addVars.root");
   }
   f->GetObject("data",tree);

   Init(tree);

   if(_irep>=0) {
   	TString fname("lightjets_filtered_addVars_resampled");
   	fname+=_irep; fname+=".root";
   	newfile = TFile::Open(fname, "RECREATE");
   	newtree = tree->CloneTree(0);
   }
}

resampleTrackIPs::~resampleTrackIPs()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t resampleTrackIPs::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t resampleTrackIPs::LoadTree(Long64_t entry)
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

void resampleTrackIPs::Init(TTree *tree)
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
   gen_pid = 0;
   gen_q = 0;
   gen_px = 0;
   gen_py = 0;
   gen_pz = 0;
   gen_e = 0;
   gen_x = 0;
   gen_y = 0;
   gen_z = 0;
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
   jet_px = 0;
   jet_py = 0;
   jet_pz = 0;
   jet_e = 0;
   jet_idx_trk0 = 0;
   jet_idx_trk1 = 0;
   jet_idx_trk2 = 0;
   jet_idx_trk3 = 0;
   jet_idx_trk4 = 0;
   jet_idx_trk5 = 0;
   jet_idx_trk6 = 0;
   jet_idx_trk7 = 0;
   jet_idx_trk8 = 0;
   jet_idx_trk9 = 0;
   jet_idx_trk10 = 0;
   jet_idx_trk11 = 0;
   jet_idx_trk12 = 0;
   jet_idx_trk13 = 0;
   jet_idx_trk14 = 0;
   jet_idx_trk15 = 0;
   jet_idx_trk16 = 0;
   jet_idx_trk17 = 0;
   jet_idx_trk18 = 0;
   jet_idx_trk19 = 0;
   jet_idx_trk20 = 0;
   jet_idx_trk21 = 0;
   jet_idx_trk22 = 0;
   jet_idx_trk23 = 0;
   jet_idx_trk24 = 0;
   jet_idx_trk25 = 0;
   jet_idx_trk26 = 0;
   jet_idx_trk27 = 0;
   jet_idx_trk28 = 0;
   jet_idx_trk29 = 0;
   jet_idx_trk30 = 0;
   jet_idx_trk31 = 0;
   jet_idx_trk32 = 0;
   jet_idx_trk33 = 0;
   jet_idx_trk34 = 0;
   jet_idx_trk35 = 0;
   jet_idx_trk36 = 0;
   jet_idx_trk37 = 0;
   jet_idx_trk38 = 0;
   jet_idx_trk39 = 0;
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
   trk_ecal = 0;
   trk_hcal = 0;
   trk_prb_ghost = 0;
   trk_type = 0;
   trk_is_mu = 0;
   trk_vid = 0;
   trk_x = 0;
   trk_y = 0;
   trk_z = 0;
   trk_vhit0 = 0;
   trk_vhit1 = 0;
   trk_vhit2 = 0;
   trk_vhit3 = 0;
   trk_vhit4 = 0;
   trk_vhit5 = 0;
   trk_vhit6 = 0;
   trk_vhit7 = 0;
   trk_vhit8 = 0;
   trk_vhit9 = 0;
   trk_vhit10 = 0;
   trk_vhit11 = 0;
   trk_vhit12 = 0;
   trk_vhit13 = 0;
   trk_vhit14 = 0;
   trk_vhit15 = 0;
   trk_vhit16 = 0;
   trk_vhit17 = 0;
   trk_vhit18 = 0;
   trk_vhit19 = 0;
   trk_vhit20 = 0;
   trk_vhit21 = 0;
   trk_vhit22 = 0;
   trk_vhit23 = 0;
   trk_vhit24 = 0;
   trk_vhit25 = 0;
   trk_vhit26 = 0;
   trk_vhit27 = 0;
   trk_vhit28 = 0;
   trk_vhit29 = 0;
   trk_vhit30 = 0;
   trk_vhit31 = 0;
   trk_vhit32 = 0;
   trk_vhit33 = 0;
   trk_vhit34 = 0;
   trk_vhit35 = 0;
   trk_vhit36 = 0;
   trk_vhit37 = 0;
   trk_vhit38 = 0;
   trk_vhit39 = 0;
   trk_vhit40 = 0;
   trk_vhit41 = 0;
   trk_vhit42 = 0;
   trk_vhit43 = 0;
   trk_vhit44 = 0;
   trk_vhit45 = 0;
   trk_vhit46 = 0;
   trk_vhit47 = 0;
   trk_vhit48 = 0;
   trk_vhit49 = 0;
   trk_vhit50 = 0;
   trk_vhit51 = 0;
   trk_vhit52 = 0;
   trk_vhit53 = 0;
   trk_vhit54 = 0;
   trk_vhit55 = 0;
   trk_vhit56 = 0;
   trk_vhit57 = 0;
   trk_vhit58 = 0;
   trk_vhit59 = 0;
   trk_vhit60 = 0;
   neu_idx_gen = 0;
   neu_idx_jet = 0;
   neu_px = 0;
   neu_py = 0;
   neu_pz = 0;
   neu_e = 0;
   neu_pid = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("gen_idx_pvr", &gen_idx_pvr, &b_gen_idx_pvr);
   fChain->SetBranchAddress("gen_idx_jet", &gen_idx_jet, &b_gen_idx_jet);
   fChain->SetBranchAddress("gen_pid", &gen_pid, &b_gen_pid);
   fChain->SetBranchAddress("gen_q", &gen_q, &b_gen_q);
   fChain->SetBranchAddress("gen_px", &gen_px, &b_gen_px);
   fChain->SetBranchAddress("gen_py", &gen_py, &b_gen_py);
   fChain->SetBranchAddress("gen_pz", &gen_pz, &b_gen_pz);
   fChain->SetBranchAddress("gen_e", &gen_e, &b_gen_e);
   fChain->SetBranchAddress("gen_x", &gen_x, &b_gen_x);
   fChain->SetBranchAddress("gen_y", &gen_y, &b_gen_y);
   fChain->SetBranchAddress("gen_z", &gen_z, &b_gen_z);
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
   fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_idx_trk0", &jet_idx_trk0, &b_jet_idx_trk0);
   fChain->SetBranchAddress("jet_idx_trk1", &jet_idx_trk1, &b_jet_idx_trk1);
   fChain->SetBranchAddress("jet_idx_trk2", &jet_idx_trk2, &b_jet_idx_trk2);
   fChain->SetBranchAddress("jet_idx_trk3", &jet_idx_trk3, &b_jet_idx_trk3);
   fChain->SetBranchAddress("jet_idx_trk4", &jet_idx_trk4, &b_jet_idx_trk4);
   fChain->SetBranchAddress("jet_idx_trk5", &jet_idx_trk5, &b_jet_idx_trk5);
   fChain->SetBranchAddress("jet_idx_trk6", &jet_idx_trk6, &b_jet_idx_trk6);
   fChain->SetBranchAddress("jet_idx_trk7", &jet_idx_trk7, &b_jet_idx_trk7);
   fChain->SetBranchAddress("jet_idx_trk8", &jet_idx_trk8, &b_jet_idx_trk8);
   fChain->SetBranchAddress("jet_idx_trk9", &jet_idx_trk9, &b_jet_idx_trk9);
   fChain->SetBranchAddress("jet_idx_trk10", &jet_idx_trk10, &b_jet_idx_trk10);
   fChain->SetBranchAddress("jet_idx_trk11", &jet_idx_trk11, &b_jet_idx_trk11);
   fChain->SetBranchAddress("jet_idx_trk12", &jet_idx_trk12, &b_jet_idx_trk12);
   fChain->SetBranchAddress("jet_idx_trk13", &jet_idx_trk13, &b_jet_idx_trk13);
   fChain->SetBranchAddress("jet_idx_trk14", &jet_idx_trk14, &b_jet_idx_trk14);
   fChain->SetBranchAddress("jet_idx_trk15", &jet_idx_trk15, &b_jet_idx_trk15);
   fChain->SetBranchAddress("jet_idx_trk16", &jet_idx_trk16, &b_jet_idx_trk16);
   fChain->SetBranchAddress("jet_idx_trk17", &jet_idx_trk17, &b_jet_idx_trk17);
   fChain->SetBranchAddress("jet_idx_trk18", &jet_idx_trk18, &b_jet_idx_trk18);
   fChain->SetBranchAddress("jet_idx_trk19", &jet_idx_trk19, &b_jet_idx_trk19);
   fChain->SetBranchAddress("jet_idx_trk20", &jet_idx_trk20, &b_jet_idx_trk20);
   fChain->SetBranchAddress("jet_idx_trk21", &jet_idx_trk21, &b_jet_idx_trk21);
   fChain->SetBranchAddress("jet_idx_trk22", &jet_idx_trk22, &b_jet_idx_trk22);
   fChain->SetBranchAddress("jet_idx_trk23", &jet_idx_trk23, &b_jet_idx_trk23);
   fChain->SetBranchAddress("jet_idx_trk24", &jet_idx_trk24, &b_jet_idx_trk24);
   fChain->SetBranchAddress("jet_idx_trk25", &jet_idx_trk25, &b_jet_idx_trk25);
   fChain->SetBranchAddress("jet_idx_trk26", &jet_idx_trk26, &b_jet_idx_trk26);
   fChain->SetBranchAddress("jet_idx_trk27", &jet_idx_trk27, &b_jet_idx_trk27);
   fChain->SetBranchAddress("jet_idx_trk28", &jet_idx_trk28, &b_jet_idx_trk28);
   fChain->SetBranchAddress("jet_idx_trk29", &jet_idx_trk29, &b_jet_idx_trk29);
   fChain->SetBranchAddress("jet_idx_trk30", &jet_idx_trk30, &b_jet_idx_trk30);
   fChain->SetBranchAddress("jet_idx_trk31", &jet_idx_trk31, &b_jet_idx_trk31);
   fChain->SetBranchAddress("jet_idx_trk32", &jet_idx_trk32, &b_jet_idx_trk32);
   fChain->SetBranchAddress("jet_idx_trk33", &jet_idx_trk33, &b_jet_idx_trk33);
   fChain->SetBranchAddress("jet_idx_trk34", &jet_idx_trk34, &b_jet_idx_trk34);
   fChain->SetBranchAddress("jet_idx_trk35", &jet_idx_trk35, &b_jet_idx_trk35);
   fChain->SetBranchAddress("jet_idx_trk36", &jet_idx_trk36, &b_jet_idx_trk36);
   fChain->SetBranchAddress("jet_idx_trk37", &jet_idx_trk37, &b_jet_idx_trk37);
   fChain->SetBranchAddress("jet_idx_trk38", &jet_idx_trk38, &b_jet_idx_trk38);
   fChain->SetBranchAddress("jet_idx_trk39", &jet_idx_trk39, &b_jet_idx_trk39);
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
   fChain->SetBranchAddress("trk_ecal", &trk_ecal, &b_trk_ecal);
   fChain->SetBranchAddress("trk_hcal", &trk_hcal, &b_trk_hcal);
   fChain->SetBranchAddress("trk_prb_ghost", &trk_prb_ghost, &b_trk_prb_ghost);
   fChain->SetBranchAddress("trk_type", &trk_type, &b_trk_type);
   fChain->SetBranchAddress("trk_is_mu", &trk_is_mu, &b_trk_is_mu);
   fChain->SetBranchAddress("trk_vid", &trk_vid, &b_trk_vid);
   fChain->SetBranchAddress("trk_x", &trk_x, &b_trk_x);
   fChain->SetBranchAddress("trk_y", &trk_y, &b_trk_y);
   fChain->SetBranchAddress("trk_z", &trk_z, &b_trk_z);
   fChain->SetBranchAddress("trk_vhit0", &trk_vhit0, &b_trk_vhit0);
   fChain->SetBranchAddress("trk_vhit1", &trk_vhit1, &b_trk_vhit1);
   fChain->SetBranchAddress("trk_vhit2", &trk_vhit2, &b_trk_vhit2);
   fChain->SetBranchAddress("trk_vhit3", &trk_vhit3, &b_trk_vhit3);
   fChain->SetBranchAddress("trk_vhit4", &trk_vhit4, &b_trk_vhit4);
   fChain->SetBranchAddress("trk_vhit5", &trk_vhit5, &b_trk_vhit5);
   fChain->SetBranchAddress("trk_vhit6", &trk_vhit6, &b_trk_vhit6);
   fChain->SetBranchAddress("trk_vhit7", &trk_vhit7, &b_trk_vhit7);
   fChain->SetBranchAddress("trk_vhit8", &trk_vhit8, &b_trk_vhit8);
   fChain->SetBranchAddress("trk_vhit9", &trk_vhit9, &b_trk_vhit9);
   fChain->SetBranchAddress("trk_vhit10", &trk_vhit10, &b_trk_vhit10);
   fChain->SetBranchAddress("trk_vhit11", &trk_vhit11, &b_trk_vhit11);
   fChain->SetBranchAddress("trk_vhit12", &trk_vhit12, &b_trk_vhit12);
   fChain->SetBranchAddress("trk_vhit13", &trk_vhit13, &b_trk_vhit13);
   fChain->SetBranchAddress("trk_vhit14", &trk_vhit14, &b_trk_vhit14);
   fChain->SetBranchAddress("trk_vhit15", &trk_vhit15, &b_trk_vhit15);
   fChain->SetBranchAddress("trk_vhit16", &trk_vhit16, &b_trk_vhit16);
   fChain->SetBranchAddress("trk_vhit17", &trk_vhit17, &b_trk_vhit17);
   fChain->SetBranchAddress("trk_vhit18", &trk_vhit18, &b_trk_vhit18);
   fChain->SetBranchAddress("trk_vhit19", &trk_vhit19, &b_trk_vhit19);
   fChain->SetBranchAddress("trk_vhit20", &trk_vhit20, &b_trk_vhit20);
   fChain->SetBranchAddress("trk_vhit21", &trk_vhit21, &b_trk_vhit21);
   fChain->SetBranchAddress("trk_vhit22", &trk_vhit22, &b_trk_vhit22);
   fChain->SetBranchAddress("trk_vhit23", &trk_vhit23, &b_trk_vhit23);
   fChain->SetBranchAddress("trk_vhit24", &trk_vhit24, &b_trk_vhit24);
   fChain->SetBranchAddress("trk_vhit25", &trk_vhit25, &b_trk_vhit25);
   fChain->SetBranchAddress("trk_vhit26", &trk_vhit26, &b_trk_vhit26);
   fChain->SetBranchAddress("trk_vhit27", &trk_vhit27, &b_trk_vhit27);
   fChain->SetBranchAddress("trk_vhit28", &trk_vhit28, &b_trk_vhit28);
   fChain->SetBranchAddress("trk_vhit29", &trk_vhit29, &b_trk_vhit29);
   fChain->SetBranchAddress("trk_vhit30", &trk_vhit30, &b_trk_vhit30);
   fChain->SetBranchAddress("trk_vhit31", &trk_vhit31, &b_trk_vhit31);
   fChain->SetBranchAddress("trk_vhit32", &trk_vhit32, &b_trk_vhit32);
   fChain->SetBranchAddress("trk_vhit33", &trk_vhit33, &b_trk_vhit33);
   fChain->SetBranchAddress("trk_vhit34", &trk_vhit34, &b_trk_vhit34);
   fChain->SetBranchAddress("trk_vhit35", &trk_vhit35, &b_trk_vhit35);
   fChain->SetBranchAddress("trk_vhit36", &trk_vhit36, &b_trk_vhit36);
   fChain->SetBranchAddress("trk_vhit37", &trk_vhit37, &b_trk_vhit37);
   fChain->SetBranchAddress("trk_vhit38", &trk_vhit38, &b_trk_vhit38);
   fChain->SetBranchAddress("trk_vhit39", &trk_vhit39, &b_trk_vhit39);
   fChain->SetBranchAddress("trk_vhit40", &trk_vhit40, &b_trk_vhit40);
   fChain->SetBranchAddress("trk_vhit41", &trk_vhit41, &b_trk_vhit41);
   fChain->SetBranchAddress("trk_vhit42", &trk_vhit42, &b_trk_vhit42);
   fChain->SetBranchAddress("trk_vhit43", &trk_vhit43, &b_trk_vhit43);
   fChain->SetBranchAddress("trk_vhit44", &trk_vhit44, &b_trk_vhit44);
   fChain->SetBranchAddress("trk_vhit45", &trk_vhit45, &b_trk_vhit45);
   fChain->SetBranchAddress("trk_vhit46", &trk_vhit46, &b_trk_vhit46);
   fChain->SetBranchAddress("trk_vhit47", &trk_vhit47, &b_trk_vhit47);
   fChain->SetBranchAddress("trk_vhit48", &trk_vhit48, &b_trk_vhit48);
   fChain->SetBranchAddress("trk_vhit49", &trk_vhit49, &b_trk_vhit49);
   fChain->SetBranchAddress("trk_vhit50", &trk_vhit50, &b_trk_vhit50);
   fChain->SetBranchAddress("trk_vhit51", &trk_vhit51, &b_trk_vhit51);
   fChain->SetBranchAddress("trk_vhit52", &trk_vhit52, &b_trk_vhit52);
   fChain->SetBranchAddress("trk_vhit53", &trk_vhit53, &b_trk_vhit53);
   fChain->SetBranchAddress("trk_vhit54", &trk_vhit54, &b_trk_vhit54);
   fChain->SetBranchAddress("trk_vhit55", &trk_vhit55, &b_trk_vhit55);
   fChain->SetBranchAddress("trk_vhit56", &trk_vhit56, &b_trk_vhit56);
   fChain->SetBranchAddress("trk_vhit57", &trk_vhit57, &b_trk_vhit57);
   fChain->SetBranchAddress("trk_vhit58", &trk_vhit58, &b_trk_vhit58);
   fChain->SetBranchAddress("trk_vhit59", &trk_vhit59, &b_trk_vhit59);
   fChain->SetBranchAddress("trk_vhit60", &trk_vhit60, &b_trk_vhit60);
   fChain->SetBranchAddress("neu_idx_gen", &neu_idx_gen, &b_neu_idx_gen);
   fChain->SetBranchAddress("neu_idx_jet", &neu_idx_jet, &b_neu_idx_jet);
   fChain->SetBranchAddress("neu_px", &neu_px, &b_neu_px);
   fChain->SetBranchAddress("neu_py", &neu_py, &b_neu_py);
   fChain->SetBranchAddress("neu_pz", &neu_pz, &b_neu_pz);
   fChain->SetBranchAddress("neu_e", &neu_e, &b_neu_e);
   fChain->SetBranchAddress("neu_pid", &neu_pid, &b_neu_pid);
   fChain->SetBranchAddress("evt_pvr_n", &evt_pvr_n, &b_evt_pvr_n);
   Notify();
}

Bool_t resampleTrackIPs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void resampleTrackIPs::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t resampleTrackIPs::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef resampleTrackIPs_cxx
