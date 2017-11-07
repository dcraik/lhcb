void plotError(TString fname, TString hname, TString plotname) {
	TFile* f = TFile::Open(fname);
	TH2D* h = dynamic_cast<TH2D*>(f->Get(hname));

	TH2D* herr = new TH2D(*h); //dynamic_cast<TH2*>(h->Clone("error"));

	for(int i=1; i<=h->GetNbinsX(); ++i) {
		for(int j=1; j<=h->GetNbinsY(); ++j) {
			std::cout << i << "\t" << j << "\t" << h->GetBinContent(i,j) << "+/-" << h->GetBinError(i,j) << std::endl;
			herr->SetBinContent(i,j,h->GetBinError(i,j));
		}
	}

	TCanvas c;
	herr->Draw("colz");
	c.Update();
	c.SaveAs(plotname);

	f->Close();
}
