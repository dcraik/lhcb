{
   Int_t nBins(5);

   TH2D* errorHi = new TH2D("errorHi","errorHi",nBins,0.,1.,nBins,0.,1.);
   TH2D* errorLo = new TH2D("errorLo","errorLo",nBins,0.,1.,nBins,0.,1.);

   TFile* f  = TFile::Open("hists/l0h_sq5_eff.root");
   TH2D*  h  = (TH2D*)f->Get("ratio");
   TFile* f0 = TFile::Open("hists/l0h_sq5_errTable.root");
   TH2D*  h0 = (TH2D*)f0->Get("error");
   TFile* f1 = TFile::Open("hists/l0h_sq5_Bd_errStat_asym.root");
   TH2D*  h1 = (TH2D*)f1->Get("errorHi");
   TH2D*  h2 = (TH2D*)f1->Get("errorLo");

   Double_t err0(0.), err1(0.), err2(0.);

   for(Int_t i=1; i<=nBins; ++i) {
      for(Int_t j=1; j<=nBins; ++j) {
         err0 = h0->GetBinContent(i,j);
	 err1 = h1->GetBinContent(i,j);
	 err2 = h2->GetBinContent(i,j);
         errorHi->SetBinContent(i,j,sqrt(err0*err0 + err1*err1));
         errorLo->SetBinContent(i,j,sqrt(err0*err0 + err2*err2));
      }
   }

   outfile = new TFile("hists/l0nothNew_sq5_all.root","RECREATE");
   h->SetName("efficiency");
   h->Write();
   errorHi->Write();
   errorLo->Write();
   outfile->Close();
}
