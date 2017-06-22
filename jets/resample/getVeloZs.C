#define getVeloZs_cxx
#include "getVeloZs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void getVeloZs::Loop()
{
   if (fChain == 0) return;

   std::set<int> vLocs;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for(uint idx=0; idx<trk_vz0->size(); ++idx) {
	      vLocs.insert(static_cast<int>(1e4*trk_vz0->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz1->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz2->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz3->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz4->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz5->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz6->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz7->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz8->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz9->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz10->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz11->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz12->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz13->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz14->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz15->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz16->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz17->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz18->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz19->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz20->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz21->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz22->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz23->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz24->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz25->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz26->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz27->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz28->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz29->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz30->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz31->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz32->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz33->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz34->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz35->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz36->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz37->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz38->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz39->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz40->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz41->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz42->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz43->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz44->at(idx)));
	      vLocs.insert(static_cast<int>(1e4*trk_vz45->at(idx)));
	      std::cout << vLocs.size() << std::endl;
      }
   }
   std::set<int>::iterator it = vLocs.begin();
   for( ; it!=vLocs.end(); ++it) {
           std::cout << static_cast<double>(*it)/1e4 << std::endl;
   }
}
