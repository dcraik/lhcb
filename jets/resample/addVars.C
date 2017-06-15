#define addVars_cxx
#include "addVars.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>

void addVars::Loop()
{
//   In a ROOT session, you can do:
//      root> .L addVars.C
//      root> addVars t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   std::vector<double>* svr_p   = new vector<double>();
   std::vector<double>* svr_fdt = new vector<double>();
   std::vector<double>* svr_fd  = new vector<double>();

   std::vector<double>* trk_p  = new vector<double>();
   std::vector<double>* trk_pt = new vector<double>();

   newtree->Branch("svr_p",   &svr_p);
   newtree->Branch("svr_fdt", &svr_fdt);
   newtree->Branch("svr_fd",  &svr_fd);

   newtree->Branch("trk_p",   &trk_p);
   newtree->Branch("trk_pt",  &trk_pt);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for( int idx=0; idx<svr_idx_pvr->size(); ++idx) {
	      int idx_pvr = svr_idx_pvr->at(idx);

	      svr_p->push_back(TMath::Sqrt(svr_px->at(idx)*svr_px->at(idx) + svr_py->at(idx)*svr_py->at(idx) + svr_pz->at(idx)*svr_pz->at(idx)));
	      double fd_x = svr_x->at(idx) - pvr_x->at(idx_pvr);
	      double fd_y = svr_y->at(idx) - pvr_y->at(idx_pvr);
	      double fd_z = svr_z->at(idx) - pvr_z->at(idx_pvr);
	      svr_fdt->push_back(TMath::Sqrt(fd_x*fd_x + fd_y*fd_y));
	      svr_fd->push_back(TMath::Sqrt(fd_x*fd_x + fd_y*fd_y + fd_z*fd_z));
      }

      for( int idx=0; idx<trk_px->size(); ++idx) {
	      trk_p->push_back(TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx) + trk_pz->at(idx)*trk_pz->at(idx)));
	      trk_pt->push_back(TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx)));
      }

      newtree->Fill();

      svr_p->clear();
      svr_fdt->clear();
      svr_fd->clear();

      trk_p->clear();
      trk_pt->clear();

   }
   newtree->AutoSave();
   newfile->Close();
}
