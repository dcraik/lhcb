#define filterLight_cxx
#include "filterLight.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void filterLight::Loop()
{
//   In a ROOT session, you can do:
//      root> .L filterLight.C
//      root> filterLight t
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      std::vector<double>::iterator it = gen_pid->begin();
      bool heavy(false);
      for( ; it!= gen_pid->end(); ++it) {
	      Int_t absid = TMath::Abs(*it);
	      Int_t q1 = absid>99?(absid%100)/10:0;
	      Int_t q2 = (absid%1000)/100;
	      Int_t q3 = (absid%10000)/1000;

	      if(absid==4 || absid==5 || q2==4 || q2==5 || q3==4 || q3==5 || q1==4 || q1==5) {
		      heavy=true;
		      break;
	      }
      }

      if(!heavy) newtree->Fill();

   }
   newtree->AutoSave();
   newfile->Close();
}

int main() {
	filterLight fl;
	fl.Loop();
}
