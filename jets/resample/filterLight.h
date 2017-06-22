//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  7 12:11:21 2017 by ROOT version 6.06/00
// from TTree data/data
// found on file: lightjets.part.root
//////////////////////////////////////////////////////////

#ifndef filterLight_h
#define filterLight_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class filterLight {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

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
   TBranch        *b_neu_idx_gen;   //!
   TBranch        *b_neu_idx_jet;   //!
   TBranch        *b_neu_px;   //!
   TBranch        *b_neu_py;   //!
   TBranch        *b_neu_pz;   //!
   TBranch        *b_neu_e;   //!
   TBranch        *b_neu_pid;   //!
   TBranch        *b_evt_pvr_n;   //!

   TFile* newfile;
   TTree* newtree;

   filterLight(TTree *tree=0);
   virtual ~filterLight();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef filterLight_cxx
filterLight::filterLight(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lightjets.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("lightjets.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);

   newfile = TFile::Open("lightjets_filtered.root","RECREATE");
   newtree = tree->CloneTree(0);

}

filterLight::~filterLight()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t filterLight::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t filterLight::LoadTree(Long64_t entry)
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

void filterLight::Init(TTree *tree)
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

Bool_t filterLight::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void filterLight::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t filterLight::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef filterLight_cxx
