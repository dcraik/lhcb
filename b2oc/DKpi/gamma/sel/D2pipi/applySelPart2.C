void applySelPart2() {

	Double_t K_ProbNNk(0.), K_ProbNNpi(0.), pi_ProbNNk(0.), pi_ProbNNpi(0);
	Double_t D0p_ProbNNk(0.), D0p_ProbNNpi(0.), D0m_ProbNNk(0.), D0m_ProbNNpi(0);

	Double_t D_ENDVERTEX_Z(0.), B_ENDVERTEX_Z(0.);

	Double_t mDmK(0.), mDmKpi(0.), mDpK(0.), mDpKpi(0.), mKpi(0.), mD0K(0.);

	Bool_t pi_isMuon(0.), D0p_isMuon(0.), D0m_isMuon(0.);

	TFile* file = TFile::Open("B2D0Kpi_D02pipi_selBd_Dsidebands_NND2pipi_addIMs.root");
	TTree* tree = (TTree*)file->Get("DecayTree");

	tree->SetBranchAddress("K_ProbNNk"     ,  &K_ProbNNk      );
	tree->SetBranchAddress("K_ProbNNpi"    ,  &K_ProbNNpi     );
	tree->SetBranchAddress("pi_ProbNNk"    ,  &pi_ProbNNk     );
	tree->SetBranchAddress("pi_ProbNNpi"   ,  &pi_ProbNNpi    );
	tree->SetBranchAddress("D0p_ProbNNk"   ,  &D0p_ProbNNk    );
	tree->SetBranchAddress("D0p_ProbNNpi"  ,  &D0p_ProbNNpi   );
	tree->SetBranchAddress("D0m_ProbNNk"   ,  &D0m_ProbNNk    );
	tree->SetBranchAddress("D0m_ProbNNpi"  ,  &D0m_ProbNNpi   );

	tree->SetBranchAddress("D_ENDVERTEX_Z" ,   &D_ENDVERTEX_Z );
	tree->SetBranchAddress("B_ENDVERTEX_Z" ,   &B_ENDVERTEX_Z );

	tree->SetBranchAddress("mDmK" ,            &mDmK          );
	tree->SetBranchAddress("mDmKpi" ,          &mDmKpi        );
	tree->SetBranchAddress("mDpK" ,            &mDpK          );
	tree->SetBranchAddress("mDpKpi" ,          &mDpKpi        );
	tree->SetBranchAddress("mKpi" ,            &mKpi          );
	tree->SetBranchAddress("mD0K" ,            &mD0K          );

	tree->SetBranchAddress("pi_isMuon"  ,      &pi_isMuon     );
	tree->SetBranchAddress("D0p_isMuon"  ,     &D0p_isMuon    );
	tree->SetBranchAddress("D0m_isMuon"  ,     &D0m_isMuon    );

	TFile* newfile = new TFile("B2D0Kpi_D02pipi_selBd_Dsidebands_vetoes_PID3_NND2pipi_addIMs.root","recreate");
	TTree* newtree = tree->CloneTree(0);

	Int_t n(tree->GetEntries());
	for(Int_t i=0; i<n; ++i) {
		if(i%10000==0) std::cout << i << " of " << n << "..." << std::endl;
		tree->GetEntry(i);

		if(K_ProbNNk*(1-K_ProbNNpi)>0.3 && pi_ProbNNpi*(1-pi_ProbNNk)>0.2 && D0p_ProbNNpi*(1-D0p_ProbNNk)>0.1 && D0m_ProbNNpi*(1-D0m_ProbNNk)>0.1) {//PID
			if((mDpK<1840||mDpK>1890)&&(mDmK<1840||mDmK>1890)&&(mKpi<1840||mKpi>1890)&&(mDpKpi<1840||mDpKpi>1890)&&(mDmKpi<1840||mDmKpi>1890)) {//wrong D vetoes
				if(D_ENDVERTEX_Z-B_ENDVERTEX_Z>0) {//forwards D
					if(pi_isMuon + D0p_isMuon + D0m_isMuon < 2) {//max 1 muon
						//if(mD0K<5200) {//B->DK veto
							newtree->Fill();
						//}
					}
				}
			}
		}
	}
	newtree->AutoSave();
	newfile->Close();
}
