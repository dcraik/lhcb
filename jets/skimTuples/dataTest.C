#define dataTest_cxx
#include "dataTest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void dataTest::Loop()
{
//   In a ROOT session, you can do:
//      root> .L dataTest.C
//      root> dataTest t
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

   int bothFound(0);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      int nMatch1(0), nMatch2(0), nMatchBoth(0);
      for(int ijet=0; ijet<jet_dR1->size(); ++ijet) {
	      if(jet_dR1->at(ijet) < 0.5) {
		      ++nMatch1;
		      if(jet_dR2->at(ijet) < 0.5) {
			      ++nMatchBoth;
		      }
	      }
	      if(jet_dR2->at(ijet) < 0.5) {
		      ++nMatch2;
	      }
      }
      if(nMatch1!=1 || nMatch2!=1 || nMatchBoth!=0) {
	      std::cout << nMatch1 << "\t" << nMatch2 << "\t" << nMatchBoth << std::endl;
      } else {
	      ++bothFound;
      }
   }

   std::cout << bothFound << "\t" << nentries << std::endl;
}
