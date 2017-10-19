void applySelPart2() {

	Double_t K_ProbNNk(0.), K_ProbNNpi(0.), pi_ProbNNk(0.), pi_ProbNNpi(0);
	Double_t D0K_ProbNNk(0.), D0K_ProbNNpi(0.), D0pi_ProbNNk(0.), D0pi_ProbNNpi(0);

	Double_t D_ENDVERTEX_Z(0.), B_ENDVERTEX_Z(0.);

	Double_t mKDpi(0.), mKpiDpi(0.), mKpi(0.), mKD(0.);

//	Bool_t pi_isMuon(0.), D0pi_isMuon(0.);

	TFile* file = TFile::Open("B2D0Kpi_D02Kpi_selBd_Dsidebands_NND2Kpi_addIMs.root");
	TTree* tree = (TTree*)file->Get("DecayTree");

	tree->SetBranchAddress("K_ProbNNk"     ,  &K_ProbNNk      );
	tree->SetBranchAddress("K_ProbNNpi"    ,  &K_ProbNNpi     );
	tree->SetBranchAddress("pi_ProbNNk"    ,  &pi_ProbNNk     );
	tree->SetBranchAddress("pi_ProbNNpi"   ,  &pi_ProbNNpi    );
	tree->SetBranchAddress("D0K_ProbNNk"   ,  &D0K_ProbNNk    );
	tree->SetBranchAddress("D0K_ProbNNpi"  ,  &D0K_ProbNNpi   );
	tree->SetBranchAddress("D0pi_ProbNNk"  ,  &D0pi_ProbNNk   );
	tree->SetBranchAddress("D0pi_ProbNNpi" ,  &D0pi_ProbNNpi  );

	tree->SetBranchAddress("D_ENDVERTEX_Z" ,   &D_ENDVERTEX_Z );
	tree->SetBranchAddress("B_ENDVERTEX_Z" ,   &B_ENDVERTEX_Z );

	tree->SetBranchAddress("mKDpi" ,           &mKDpi         );
	tree->SetBranchAddress("mKpiDpi" ,         &mKpiDpi       );
	tree->SetBranchAddress("mKpi" ,            &mKpi          );
	tree->SetBranchAddress("mKD" ,             &mKD           );

//	tree->SetBranchAddress("pi_isMuon"  ,      &pi_isMuon     );
//	tree->SetBranchAddress("D0pi_isMuon" ,     &D0pi_isMuon   );

	TFile* newfile = new TFile("B2D0Kpi_D02Kpi_selBd_Dsidebands_vetoes_PID3_NND2Kpi_addIMs.root","recreate");
	TTree* newtree = tree->CloneTree(0);

	Int_t n(tree->GetEntries());
	for(Int_t i=0; i<n; ++i) {
		if(i%10000==0) std::cout << i << " of " << n << "..." << std::endl;
		tree->GetEntry(i);

		if(K_ProbNNk*(1-K_ProbNNpi)>0.3 && pi_ProbNNpi*(1-pi_ProbNNk)>0.2 && D0K_ProbNNk*(1-D0K_ProbNNpi)>0.1 && D0pi_ProbNNpi*(1-D0pi_ProbNNk)>0.1) {//PID
			if((mKDpi<1840||mKDpi>1890)&&(mKpi<1840||mKpi>1890)&&(mKpiDpi<1840||mKpiDpi>1890)) {//wrong D vetoes
				if(D_ENDVERTEX_Z-B_ENDVERTEX_Z>0) {//forwards D
//					if(pi_isMuon + D0pi_isMuon < 2) {//max 1 muon
						if(mKD<5200) {//B->DK veto
							newtree->Fill();
						}
//					}
				}
			}
		}
	}
	newtree->AutoSave();
	newfile->Close();
}
