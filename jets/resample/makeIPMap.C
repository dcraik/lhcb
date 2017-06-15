#define makeIPMap_cxx
#include "makeIPMap.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void makeIPMap::Loop()
{
//   In a ROOT session, you can do:
//      root> .L makeIPMap.C
//      root> makeIPMap t
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
   if (fChain == 0) return;

   double maxX(0.), minX(1.), maxY(0.), minY(1.);
   double binsx[59] = {0.,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,6.,8.,10.,20.,30.,75.,150.};
//   double binsy[5] = {0.000,0.0005,0.001,0.0015,0.002};
   double binsy[4] = {0.000,0.0005,0.001,0.002};
   TH2D h("h","", 58, binsx, 3, binsy);
   TH2D h2("h2","", 58, binsx, 3, binsy);

   int pvr_idx(0);
   double tx(0.), ty(0.), dx(0.), dy(0.), sigma(0.), invpt(0.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for(int idx=0; idx!= trk_z->size(); ++idx) {
		int gen_idx = trk_idx_gen->at(idx);
		if(gen_idx<0) continue;
		int pvr_idx = gen_idx_pvr->at(gen_idx);
		if(pvr_idx<0) continue;
		if(TMath::Abs(gen_x->at(gen_idx) - pvr_x->at(pvr_idx)) > 1e-3 ||
		   TMath::Abs(gen_y->at(gen_idx) - pvr_y->at(pvr_idx)) > 1e-3 ||
		   TMath::Abs(gen_z->at(gen_idx) - pvr_z->at(pvr_idx)) > 1e-3) continue;//only use prompt tracks
	      if(trk_vid->at(idx) == 1) continue; //ignore upstream tracks TODO can we ignore VELO tracks here too?
	      if(trk_pt->at(idx)<500 || trk_prb_ghost->at(idx)>0.3) continue; //check track meets requirements that won't change (Phys/JetTagging/src/LoKiBDTTag.cpp:262)

	      //int pvr_idx = trk_idx_pvr->at(idx);

	      invpt = 1./TMath::Sqrt(trk_px->at(idx)*trk_px->at(idx) + trk_py->at(idx)*trk_py->at(idx));
	      //invpt = 1./TMath::Sqrt(gen_px->at(gen_idx)*gen_px->at(gen_idx) + gen_py->at(gen_idx)*gen_py->at(gen_idx));

	      tx = trk_px->at(idx) / trk_pz->at(idx);
	      ty = trk_py->at(idx) / trk_pz->at(idx);

	      dx = TMath::Abs(trk_x->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*tx - pvr_x->at(pvr_idx));
	      dy = TMath::Abs(trk_y->at(idx) + (pvr_z->at(pvr_idx) - trk_z->at(idx))*ty - pvr_y->at(pvr_idx));

      	      sigma = TMath::Sqrt(trk_ip->at(idx)*trk_ip->at(idx) / trk_ip_chi2->at(idx)); //assume unchanged by shift

	      if(TMath::Sqrt(dx*dx+dy*dy)/sigma>maxX) maxX=TMath::Sqrt(dx*dx+dy*dy)/sigma;
	      if(dx/sigma<minX) minX=dx/sigma;
	      if(dy/sigma<minX) minX=dy/sigma;
	      if(invpt>maxY) maxY=invpt;
	      if(invpt<minY) minY=invpt;

	      h2.Fill(TMath::Sqrt(dx*dx+dy*dy)/sigma, invpt);
	      h.Fill(dx/sigma, invpt);
	      h.Fill(dy/sigma, invpt);
      }

   }

   TFile* f = TFile::Open("ipmap.root","RECREATE");
   h.SetName("ipx");
   h.Write();
   h2.SetName("ipr");
   h2.Write();
   f->Close();

   std::cout << maxX << "\t" << minX << "\t" << maxY << "\t" << minY << std::endl;
   //138.68  5.6146e-06      0.00995349      2.36823e-05
}
