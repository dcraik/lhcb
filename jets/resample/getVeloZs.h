//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 16 21:40:33 2017 by ROOT version 5.34/36
// from TTree data/data
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef getVeloZs_h
#define getVeloZs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class getVeloZs {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *gen_idx_pvr;
   vector<double>  *gen_idx_jet;
   vector<double>  *gen_pid;
   vector<double>  *gen_q;
   vector<double>  *gen_px;
   vector<double>  *gen_py;
   vector<double>  *gen_pz;
   vector<double>  *gen_e;
   vector<double>  *gen_x;
   vector<double>  *gen_y;
   vector<double>  *gen_z;
   vector<double>  *pvr_x;
   vector<double>  *pvr_y;
   vector<double>  *pvr_z;
   vector<double>  *pvr_dx;
   vector<double>  *pvr_dy;
   vector<double>  *pvr_dz;
   vector<double>  *pvr_chi2;
   vector<double>  *pvr_ndof;
   vector<double>  *svr_idx_pvr;
   vector<double>  *svr_idx_jet;
   vector<double>  *svr_idx_trk0;
   vector<double>  *svr_idx_trk1;
   vector<double>  *svr_idx_trk2;
   vector<double>  *svr_idx_trk3;
   vector<double>  *svr_idx_trk4;
   vector<double>  *svr_idx_trk5;
   vector<double>  *svr_idx_trk6;
   vector<double>  *svr_idx_trk7;
   vector<double>  *svr_idx_trk8;
   vector<double>  *svr_idx_trk9;
   vector<double>  *svr_px;
   vector<double>  *svr_py;
   vector<double>  *svr_pz;
   vector<double>  *svr_e;
   vector<double>  *svr_x;
   vector<double>  *svr_y;
   vector<double>  *svr_z;
   vector<double>  *svr_fd_min;
   vector<double>  *svr_m_cor;
   vector<double>  *jet_idx_pvr;
   vector<double>  *jet_px;
   vector<double>  *jet_py;
   vector<double>  *jet_pz;
   vector<double>  *jet_e;
   vector<double>  *trk_idx_gen;
   vector<double>  *trk_idx_pvr;
   vector<double>  *trk_idx_jet;
   vector<double>  *trk_px;
   vector<double>  *trk_py;
   vector<double>  *trk_pz;
   vector<double>  *trk_e;
   vector<double>  *trk_pid;
   vector<double>  *trk_q;
   vector<double>  *trk_ip;
   vector<double>  *trk_ip_chi2;
   vector<double>  *trk_pnn_e;
   vector<double>  *trk_pnn_mu;
   vector<double>  *trk_pnn_pi;
   vector<double>  *trk_pnn_k;
   vector<double>  *trk_pnn_p;
   vector<double>  *trk_ecal;
   vector<double>  *trk_hcal;
   vector<double>  *trk_prb_ghost;
   vector<double>  *trk_is_mu;
   vector<double>  *trk_vid;
   vector<double>  *trk_x;
   vector<double>  *trk_y;
   vector<double>  *trk_z;
   vector<double>  *trk_vhit0;
   vector<double>  *trk_vhit1;
   vector<double>  *trk_vhit2;
   vector<double>  *trk_vhit3;
   vector<double>  *trk_vhit4;
   vector<double>  *trk_vhit5;
   vector<double>  *trk_vhit6;
   vector<double>  *trk_vhit7;
   vector<double>  *trk_vhit8;
   vector<double>  *trk_vhit9;
   vector<double>  *trk_vhit10;
   vector<double>  *trk_vhit11;
   vector<double>  *trk_vhit12;
   vector<double>  *trk_vhit13;
   vector<double>  *trk_vhit14;
   vector<double>  *trk_vhit15;
   vector<double>  *trk_vhit16;
   vector<double>  *trk_vhit17;
   vector<double>  *trk_vhit18;
   vector<double>  *trk_vhit19;
   vector<double>  *trk_vhit20;
   vector<double>  *trk_vhit21;
   vector<double>  *trk_vhit22;
   vector<double>  *trk_vhit23;
   vector<double>  *trk_vhit24;
   vector<double>  *trk_vhit25;
   vector<double>  *trk_vhit26;
   vector<double>  *trk_vhit27;
   vector<double>  *trk_vhit28;
   vector<double>  *trk_vhit29;
   vector<double>  *trk_vhit30;
   vector<double>  *trk_vhit31;
   vector<double>  *trk_vhit32;
   vector<double>  *trk_vhit33;
   vector<double>  *trk_vhit34;
   vector<double>  *trk_vhit35;
   vector<double>  *trk_vhit36;
   vector<double>  *trk_vhit37;
   vector<double>  *trk_vhit38;
   vector<double>  *trk_vhit39;
   vector<double>  *trk_vhit40;
   vector<double>  *trk_vhit41;
   vector<double>  *trk_vhit42;
   vector<double>  *trk_vhit43;
   vector<double>  *trk_vhit44;
   vector<double>  *trk_vhit45;
   vector<double>  *trk_vz0;
   vector<double>  *trk_vz1;
   vector<double>  *trk_vz2;
   vector<double>  *trk_vz3;
   vector<double>  *trk_vz4;
   vector<double>  *trk_vz5;
   vector<double>  *trk_vz6;
   vector<double>  *trk_vz7;
   vector<double>  *trk_vz8;
   vector<double>  *trk_vz9;
   vector<double>  *trk_vz10;
   vector<double>  *trk_vz11;
   vector<double>  *trk_vz12;
   vector<double>  *trk_vz13;
   vector<double>  *trk_vz14;
   vector<double>  *trk_vz15;
   vector<double>  *trk_vz16;
   vector<double>  *trk_vz17;
   vector<double>  *trk_vz18;
   vector<double>  *trk_vz19;
   vector<double>  *trk_vz20;
   vector<double>  *trk_vz21;
   vector<double>  *trk_vz22;
   vector<double>  *trk_vz23;
   vector<double>  *trk_vz24;
   vector<double>  *trk_vz25;
   vector<double>  *trk_vz26;
   vector<double>  *trk_vz27;
   vector<double>  *trk_vz28;
   vector<double>  *trk_vz29;
   vector<double>  *trk_vz30;
   vector<double>  *trk_vz31;
   vector<double>  *trk_vz32;
   vector<double>  *trk_vz33;
   vector<double>  *trk_vz34;
   vector<double>  *trk_vz35;
   vector<double>  *trk_vz36;
   vector<double>  *trk_vz37;
   vector<double>  *trk_vz38;
   vector<double>  *trk_vz39;
   vector<double>  *trk_vz40;
   vector<double>  *trk_vz41;
   vector<double>  *trk_vz42;
   vector<double>  *trk_vz43;
   vector<double>  *trk_vz44;
   vector<double>  *trk_vz45;
   vector<double>  *neu_idx_gen;
   vector<double>  *neu_idx_jet;
   vector<double>  *neu_px;
   vector<double>  *neu_py;
   vector<double>  *neu_pz;
   vector<double>  *neu_e;
   vector<double>  *neu_pid;
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
   TBranch        *b_svr_fd_min;   //!
   TBranch        *b_svr_m_cor;   //!
   TBranch        *b_jet_idx_pvr;   //!
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
   TBranch        *b_trk_ecal;   //!
   TBranch        *b_trk_hcal;   //!
   TBranch        *b_trk_prb_ghost;   //!
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
   TBranch        *b_trk_vz0;   //!
   TBranch        *b_trk_vz1;   //!
   TBranch        *b_trk_vz2;   //!
   TBranch        *b_trk_vz3;   //!
   TBranch        *b_trk_vz4;   //!
   TBranch        *b_trk_vz5;   //!
   TBranch        *b_trk_vz6;   //!
   TBranch        *b_trk_vz7;   //!
   TBranch        *b_trk_vz8;   //!
   TBranch        *b_trk_vz9;   //!
   TBranch        *b_trk_vz10;   //!
   TBranch        *b_trk_vz11;   //!
   TBranch        *b_trk_vz12;   //!
   TBranch        *b_trk_vz13;   //!
   TBranch        *b_trk_vz14;   //!
   TBranch        *b_trk_vz15;   //!
   TBranch        *b_trk_vz16;   //!
   TBranch        *b_trk_vz17;   //!
   TBranch        *b_trk_vz18;   //!
   TBranch        *b_trk_vz19;   //!
   TBranch        *b_trk_vz20;   //!
   TBranch        *b_trk_vz21;   //!
   TBranch        *b_trk_vz22;   //!
   TBranch        *b_trk_vz23;   //!
   TBranch        *b_trk_vz24;   //!
   TBranch        *b_trk_vz25;   //!
   TBranch        *b_trk_vz26;   //!
   TBranch        *b_trk_vz27;   //!
   TBranch        *b_trk_vz28;   //!
   TBranch        *b_trk_vz29;   //!
   TBranch        *b_trk_vz30;   //!
   TBranch        *b_trk_vz31;   //!
   TBranch        *b_trk_vz32;   //!
   TBranch        *b_trk_vz33;   //!
   TBranch        *b_trk_vz34;   //!
   TBranch        *b_trk_vz35;   //!
   TBranch        *b_trk_vz36;   //!
   TBranch        *b_trk_vz37;   //!
   TBranch        *b_trk_vz38;   //!
   TBranch        *b_trk_vz39;   //!
   TBranch        *b_trk_vz40;   //!
   TBranch        *b_trk_vz41;   //!
   TBranch        *b_trk_vz42;   //!
   TBranch        *b_trk_vz43;   //!
   TBranch        *b_trk_vz44;   //!
   TBranch        *b_trk_vz45;   //!
   TBranch        *b_neu_idx_gen;   //!
   TBranch        *b_neu_idx_jet;   //!
   TBranch        *b_neu_px;   //!
   TBranch        *b_neu_py;   //!
   TBranch        *b_neu_pz;   //!
   TBranch        *b_neu_e;   //!
   TBranch        *b_neu_pid;   //!
   TBranch        *b_evt_pvr_n;   //!

   getVeloZs(TTree *tree=0);
   virtual ~getVeloZs();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef getVeloZs_cxx
getVeloZs::getVeloZs(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);
}

getVeloZs::~getVeloZs()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t getVeloZs::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t getVeloZs::LoadTree(Long64_t entry)
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

void getVeloZs::Init(TTree *tree)
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
   svr_fd_min = 0;
   svr_m_cor = 0;
   jet_idx_pvr = 0;
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
   trk_ecal = 0;
   trk_hcal = 0;
   trk_prb_ghost = 0;
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
   trk_vz0 = 0;
   trk_vz1 = 0;
   trk_vz2 = 0;
   trk_vz3 = 0;
   trk_vz4 = 0;
   trk_vz5 = 0;
   trk_vz6 = 0;
   trk_vz7 = 0;
   trk_vz8 = 0;
   trk_vz9 = 0;
   trk_vz10 = 0;
   trk_vz11 = 0;
   trk_vz12 = 0;
   trk_vz13 = 0;
   trk_vz14 = 0;
   trk_vz15 = 0;
   trk_vz16 = 0;
   trk_vz17 = 0;
   trk_vz18 = 0;
   trk_vz19 = 0;
   trk_vz20 = 0;
   trk_vz21 = 0;
   trk_vz22 = 0;
   trk_vz23 = 0;
   trk_vz24 = 0;
   trk_vz25 = 0;
   trk_vz26 = 0;
   trk_vz27 = 0;
   trk_vz28 = 0;
   trk_vz29 = 0;
   trk_vz30 = 0;
   trk_vz31 = 0;
   trk_vz32 = 0;
   trk_vz33 = 0;
   trk_vz34 = 0;
   trk_vz35 = 0;
   trk_vz36 = 0;
   trk_vz37 = 0;
   trk_vz38 = 0;
   trk_vz39 = 0;
   trk_vz40 = 0;
   trk_vz41 = 0;
   trk_vz42 = 0;
   trk_vz43 = 0;
   trk_vz44 = 0;
   trk_vz45 = 0;
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
   fChain->SetBranchAddress("svr_fd_min", &svr_fd_min, &b_svr_fd_min);
   fChain->SetBranchAddress("svr_m_cor", &svr_m_cor, &b_svr_m_cor);
   fChain->SetBranchAddress("jet_idx_pvr", &jet_idx_pvr, &b_jet_idx_pvr);
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
   fChain->SetBranchAddress("trk_ecal", &trk_ecal, &b_trk_ecal);
   fChain->SetBranchAddress("trk_hcal", &trk_hcal, &b_trk_hcal);
   fChain->SetBranchAddress("trk_prb_ghost", &trk_prb_ghost, &b_trk_prb_ghost);
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
   fChain->SetBranchAddress("trk_vz0", &trk_vz0, &b_trk_vz0);
   fChain->SetBranchAddress("trk_vz1", &trk_vz1, &b_trk_vz1);
   fChain->SetBranchAddress("trk_vz2", &trk_vz2, &b_trk_vz2);
   fChain->SetBranchAddress("trk_vz3", &trk_vz3, &b_trk_vz3);
   fChain->SetBranchAddress("trk_vz4", &trk_vz4, &b_trk_vz4);
   fChain->SetBranchAddress("trk_vz5", &trk_vz5, &b_trk_vz5);
   fChain->SetBranchAddress("trk_vz6", &trk_vz6, &b_trk_vz6);
   fChain->SetBranchAddress("trk_vz7", &trk_vz7, &b_trk_vz7);
   fChain->SetBranchAddress("trk_vz8", &trk_vz8, &b_trk_vz8);
   fChain->SetBranchAddress("trk_vz9", &trk_vz9, &b_trk_vz9);
   fChain->SetBranchAddress("trk_vz10", &trk_vz10, &b_trk_vz10);
   fChain->SetBranchAddress("trk_vz11", &trk_vz11, &b_trk_vz11);
   fChain->SetBranchAddress("trk_vz12", &trk_vz12, &b_trk_vz12);
   fChain->SetBranchAddress("trk_vz13", &trk_vz13, &b_trk_vz13);
   fChain->SetBranchAddress("trk_vz14", &trk_vz14, &b_trk_vz14);
   fChain->SetBranchAddress("trk_vz15", &trk_vz15, &b_trk_vz15);
   fChain->SetBranchAddress("trk_vz16", &trk_vz16, &b_trk_vz16);
   fChain->SetBranchAddress("trk_vz17", &trk_vz17, &b_trk_vz17);
   fChain->SetBranchAddress("trk_vz18", &trk_vz18, &b_trk_vz18);
   fChain->SetBranchAddress("trk_vz19", &trk_vz19, &b_trk_vz19);
   fChain->SetBranchAddress("trk_vz20", &trk_vz20, &b_trk_vz20);
   fChain->SetBranchAddress("trk_vz21", &trk_vz21, &b_trk_vz21);
   fChain->SetBranchAddress("trk_vz22", &trk_vz22, &b_trk_vz22);
   fChain->SetBranchAddress("trk_vz23", &trk_vz23, &b_trk_vz23);
   fChain->SetBranchAddress("trk_vz24", &trk_vz24, &b_trk_vz24);
   fChain->SetBranchAddress("trk_vz25", &trk_vz25, &b_trk_vz25);
   fChain->SetBranchAddress("trk_vz26", &trk_vz26, &b_trk_vz26);
   fChain->SetBranchAddress("trk_vz27", &trk_vz27, &b_trk_vz27);
   fChain->SetBranchAddress("trk_vz28", &trk_vz28, &b_trk_vz28);
   fChain->SetBranchAddress("trk_vz29", &trk_vz29, &b_trk_vz29);
   fChain->SetBranchAddress("trk_vz30", &trk_vz30, &b_trk_vz30);
   fChain->SetBranchAddress("trk_vz31", &trk_vz31, &b_trk_vz31);
   fChain->SetBranchAddress("trk_vz32", &trk_vz32, &b_trk_vz32);
   fChain->SetBranchAddress("trk_vz33", &trk_vz33, &b_trk_vz33);
   fChain->SetBranchAddress("trk_vz34", &trk_vz34, &b_trk_vz34);
   fChain->SetBranchAddress("trk_vz35", &trk_vz35, &b_trk_vz35);
   fChain->SetBranchAddress("trk_vz36", &trk_vz36, &b_trk_vz36);
   fChain->SetBranchAddress("trk_vz37", &trk_vz37, &b_trk_vz37);
   fChain->SetBranchAddress("trk_vz38", &trk_vz38, &b_trk_vz38);
   fChain->SetBranchAddress("trk_vz39", &trk_vz39, &b_trk_vz39);
   fChain->SetBranchAddress("trk_vz40", &trk_vz40, &b_trk_vz40);
   fChain->SetBranchAddress("trk_vz41", &trk_vz41, &b_trk_vz41);
   fChain->SetBranchAddress("trk_vz42", &trk_vz42, &b_trk_vz42);
   fChain->SetBranchAddress("trk_vz43", &trk_vz43, &b_trk_vz43);
   fChain->SetBranchAddress("trk_vz44", &trk_vz44, &b_trk_vz44);
   fChain->SetBranchAddress("trk_vz45", &trk_vz45, &b_trk_vz45);
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

Bool_t getVeloZs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void getVeloZs::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t getVeloZs::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef getVeloZs_cxx
