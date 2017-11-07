#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"

TChain* loadZChain(int n) {

	char str[256];
	TChain* t(0);
	t = new TChain("data");
	for(int i=0; i<(n+1)/2; ++i) {
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/296/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}
	for(int i=0; i<n/2; ++i) {
		sprintf(str,"/eos/lhcb/user/d/dcraik/jets/297/%d/output.root",i);
		if(gSystem->AccessPathName(str)) continue;
		t->Add(str);
	}

	return t;
}
