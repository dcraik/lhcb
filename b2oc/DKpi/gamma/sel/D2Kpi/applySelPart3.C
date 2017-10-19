void applySelPart3() {

	TLorentzVector Dd1_4mom;
	TLorentzVector Dd2_4mom;
	TLorentzVector D0_4mom;
	TLorentzVector Pi_4mom;
	TLorentzVector D0Pi;

	TLorentzVector K_misID_4mom;
	TLorentzVector D0K_misID;

	Float_t  Bd_K_E[100],   Bd_K_PX[100],   Bd_K_PY[100],   Bd_K_PZ[100];
	Float_t  Bd_pi_E[100],  Bd_pi_PX[100],  Bd_pi_PY[100],  Bd_pi_PZ[100];
	Float_t  Bd_Dd1_E[100], Bd_Dd1_PX[100], Bd_Dd1_PY[100], Bd_Dd1_PZ[100];
	Float_t  Bd_Dd2_E[100], Bd_Dd2_PX[100], Bd_Dd2_PY[100], Bd_Dd2_PZ[100];

	Double_t D_M(0.);

	Float_t NN(0.);

	TFile* file = TFile::Open("B2D0Kpi_D02Kpi_selBd_Dsidebands_vetoes_PID3_NND2Kpi_addIMs.root");
	TTree* tree = (TTree*)file->Get("DecayTree");

	tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PE"        ,    Bd_Dd1_E ); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PX"        ,    Bd_Dd1_PX); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PY"        ,    Bd_Dd1_PY); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PZ"        ,    Bd_Dd1_PZ); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_piplus_PE"      ,    Bd_Dd2_E ); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_piplus_PX"      ,    Bd_Dd2_PX); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_piplus_PY"      ,    Bd_Dd2_PY); 
	tree->SetBranchAddress("B_ConstB0Fit_D0_piplus_PZ"      ,    Bd_Dd2_PZ); 

	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PE"  ,    Bd_K_E   ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PX"  ,    Bd_K_PX  ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PY"  ,    Bd_K_PY  ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PZ"  ,    Bd_K_PZ  ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PE" ,    Bd_pi_E  ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PX" ,    Bd_pi_PX ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PY" ,    Bd_pi_PY ); 
	tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PZ" ,    Bd_pi_PZ );

	tree->SetBranchAddress("D_M"                              ,   &D_M      );

	tree->SetBranchAddress("NN"                               ,   &NN       );

	TFile* newfile = new TFile("B2D0Kpi_D02Kpi_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root","recreate");
	TTree* newtree = tree->CloneTree(0);

	Int_t num(0), denom(0);

	Int_t n(tree->GetEntries());
	for(Int_t i=0; i<n; ++i) {
		if(i%10000==0) std::cout << i << " of " << n << "..." << std::endl;
		tree->GetEntry(i);

		Dd1_4mom.SetPxPyPzE(     Bd_Dd1_PX[0]/1000, Bd_Dd1_PY[0]/1000, Bd_Dd1_PZ[0]/1000, Bd_Dd1_E[0]/1000);
		Dd2_4mom.SetPxPyPzE(     Bd_Dd2_PX[0]/1000, Bd_Dd2_PY[0]/1000, Bd_Dd2_PZ[0]/1000, Bd_Dd2_E[0]/1000);
		Pi_4mom.SetPxPyPzE(      Bd_pi_PX[0]/1000,  Bd_pi_PY[0]/1000,  Bd_pi_PZ[0]/1000,  Bd_pi_E[0]/1000);

		K_misID_4mom.SetPxPyPzE( Bd_K_PX[0]/1000,   Bd_K_PY[0]/1000,   Bd_K_PZ[0]/1000,   sqrt((Bd_K_E[0]/1000.)**2-0.224237195));

		D0_4mom = Dd1_4mom + Dd2_4mom;

		D0Pi = D0_4mom + Pi_4mom;
		D0K_misID = D0_4mom + K_misID_4mom;

		if(D_M>1846. && D_M<1887) {
			if(NN>-0.8) ++denom;
			if(D0Pi.M() - D0_4mom.M() < (145.4-2.5)/1000. || D0Pi.M() - D0_4mom.M() > (145.4+2.5)/1000.) {//m(D pi) D*+
				if(D0K_misID.M() - D0_4mom.M() < (145.4-2.5)/1000. || D0K_misID.M() - D0_4mom.M() > (145.4+2.5)/1000.) {//m(D K) D*+
					if(NN>-0.8) ++num;
					newtree->Fill();
				}
			}
		}
	}
	newtree->AutoSave();
	newfile->Close();

	std::cout << num << " of " << denom << std::endl;
}
