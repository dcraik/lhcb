void applyVetoes() {

	Double_t Bd_CM_m12;

	Double_t B_D0_D_CM_E;
	Double_t B_D0_D_CM_PX;
	Double_t B_D0_D_CM_PY;
	Double_t B_D0_D_CM_PZ;

	Double_t B_D0_p_CM_E;
	Double_t B_D0_p_CM_PX;
	Double_t B_D0_p_CM_PY;
	Double_t B_D0_p_CM_PZ;
	Double_t B_D0_p_CM_ID;

	Double_t B_D0_m_CM_E;
	Double_t B_D0_m_CM_PX;
	Double_t B_D0_m_CM_PY;
	Double_t B_D0_m_CM_PZ;
	Double_t B_D0_m_ID;

	Double_t mDK;

	TFile* file = TFile::Open("B2D0Kpi_D02pipi_selBd_Dst25_PID3_NND2pipi.root");
	TTree* tree = (TTree*)file->Get("DecayTree");

	tree->SetBranchAddress("Bd_CM_m12" ,      &Bd_CM_m12    );
	tree->SetBranchAddress("B_D0_D_CM_E" ,    &B_D0_D_CM_E  );
	tree->SetBranchAddress("B_D0_D_CM_PX" ,   &B_D0_D_CM_PX );
	tree->SetBranchAddress("B_D0_D_CM_PY" ,   &B_D0_D_CM_PY );
	tree->SetBranchAddress("B_D0_D_CM_PZ" ,   &B_D0_D_CM_PZ );
	tree->SetBranchAddress("B_D0_p_CM_E" ,    &B_D0_p_CM_E  );
	tree->SetBranchAddress("B_D0_p_CM_PX" ,   &B_D0_p_CM_PX );
	tree->SetBranchAddress("B_D0_p_CM_PY" ,   &B_D0_p_CM_PY );
	tree->SetBranchAddress("B_D0_p_CM_PZ" ,   &B_D0_p_CM_PZ );
	tree->SetBranchAddress("B_D0_p_CM_ID" ,   &B_D0_p_CM_ID );
	tree->SetBranchAddress("B_D0_m_CM_E" ,    &B_D0_m_CM_E  );
	tree->SetBranchAddress("B_D0_m_CM_PX" ,   &B_D0_m_CM_PX );
	tree->SetBranchAddress("B_D0_m_CM_PY" ,   &B_D0_m_CM_PY );
	tree->SetBranchAddress("B_D0_m_CM_PZ" ,   &B_D0_m_CM_PZ );

	TFile* newfile = new TFile("B2D0Kpi_D02pipi_selBd_Dst25_DD55_DK52_PID3_NND2pipi.root","recreate");
	TTree* newtree = tree->CloneTree(0);

	Int_t n(tree->GetEntries());
	for(Int_t i=0; i<n; ++i) {
		tree->GetEntry(i);

		if(B_D0_p_CM_ID==321) {
			mDK=sqrt((B_D0_D_CM_E+B_D0_p_CM_E)*(B_D0_D_CM_E+B_D0_p_CM_E)-(B_D0_D_CM_PX+B_D0_p_CM_PX)*(B_D0_D_CM_PX+B_D0_p_CM_PX)-(B_D0_D_CM_PY+B_D0_p_CM_PY)*(B_D0_D_CM_PY+B_D0_p_CM_PY)-(B_D0_D_CM_PZ+B_D0_p_CM_PZ)*(B_D0_D_CM_PZ+B_D0_p_CM_PZ));
		} else {
			mDK=sqrt((B_D0_D_CM_E+B_D0_m_CM_E)*(B_D0_D_CM_E+B_D0_m_CM_E)-(B_D0_D_CM_PX+B_D0_m_CM_PX)*(B_D0_D_CM_PX+B_D0_m_CM_PX)-(B_D0_D_CM_PY+B_D0_m_CM_PY)*(B_D0_D_CM_PY+B_D0_m_CM_PY)-(B_D0_D_CM_PZ+B_D0_m_CM_PZ)*(B_D0_D_CM_PZ+B_D0_m_CM_PZ));
		}

		if(Bd_CM_m12<1.835 || Bd_CM_m12>1.890) {//m(Kpi) D0
			if(mDK<5200) {
				newtree->Fill();
			}
		}
	}
	newtree->AutoSave();
	newfile->Close();
}
