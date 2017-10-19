{
	TFile* f = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2Kpi_Scaled/B2D0Kpi_D02Kpi_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("DecayTree"));

	TFile* fd = new TFile("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2Kpi_Scaled/B2D0Kpi_D02Kpi_MD_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root","RECREATE");
	TTree* td = t->CloneTree(0);
	TFile* fu = new TFile("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2Kpi_Scaled/B2D0Kpi_D02Kpi_MU_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root","RECREATE");
	TTree* tu = t->CloneTree(0);

	Short_t Polarity(0);
	t->SetBranchAddress("Polarity", &Polarity);

	Int_t n = t->GetEntries();

	for(Int_t i=0; i<n; ++i) {
		t->GetEntry(i);

		if(Polarity==-1) td->Fill();
		else if(Polarity==1) tu->Fill();
		else std::cout << "Bad polarity" << std::endl;
	}
	td->AutoSave();
	tu->AutoSave();

	fd->Close();
	fu->Close();
}
