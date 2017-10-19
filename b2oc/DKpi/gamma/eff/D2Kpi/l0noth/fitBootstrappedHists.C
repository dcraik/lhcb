void fit(TH1D* h, Double_t& sh, Double_t& sl, TString plotName="") {
   RooRealVar eff("eff","",0.8,1.2);
   RooDataHist * data = new RooDataHist("data", "", RooArgSet(eff), h);

   RooRealVar mean("mean","",1.0,0.8,1.2,"");
   //RooRealVar sigmaHi("sigmaHi","",0.01,0.003,0.1,"");
   RooRealVar sigmaLo("sigmaLo","",0.05,0.003,0.15,"");
   RooRealVar ratio("ratio","",1.0,0.9,1.1,"");
   RooProduct sigmaHi("sigmaHi","",RooArgSet(sigmaLo,ratio));

   RooBifurGauss bifurGauss("bifurGauss","",eff,mean,sigmaLo,sigmaHi);
   bifurGauss->fitTo(*data);

   cout << mean.getVal() << "\t" << sigmaHi.getVal() << "\t" << sigmaLo.getVal() << endl;

   sh = sigmaHi.getVal();
   sl = sigmaLo.getVal();

   if(plotName) {
      TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
      plot = eff.frame(100);
      plot->SetTitle("");
      plot->GetYaxis()->SetTitle("Toys");
      plot->GetXaxis()->SetTitle("Efficiency");
      
      data->plotOn(plot);
      bifurGauss.plotOn(plot);
      plot->Draw();
      can->SaveAs(plotName);
   }
}

void fitBootstrappedHists() {
   gSystem->Load("libRooFit");

   Int_t nBins(5);

   TH2D* err  = new TH2D("err", "err", nBins,0.,1.,nBins,0.,1.);
   TH2D* errp = new TH2D("errp","errp",nBins,0.,1.,nBins,0.,1.);
   TH2D* errm = new TH2D("errm","errm",nBins,0.,1.,nBins,0.,1.);

   Double_t sigHi(0.);
   Double_t sigLo(0.);

   TFile* f = TFile::Open("l0h_sq_Bd_bootstrapped.root");
   for(Int_t i=1; i<=nBins; ++i) {
      for(Int_t j=1; j<=nBins; ++j) {
         TString hName = "h_";
	 hName += i;
	 hName += "_";
	 hName += j;
         TH1D* h = f->Get(hName);
	 TString pName = "latexBootstrappedFits/figs/";
	 pName += i;
	 pName += "_";
	 pName += j;
	 pName += ".pdf";
	 fit(h,sigHi,sigLo,pName);
         h->Fit("gaus");
	 cout << i << "\t" << j << "\t" << h->GetFunction("gaus")->GetParameter(1) << "+/-" << h->GetFunction("gaus")->GetParameter(2) << endl;
	 err->SetBinContent(i,j,h->GetFunction("gaus")->GetParameter(2));
	 errp->SetBinContent(i,j,sigHi);
	 errm->SetBinContent(i,j,sigLo);
      }
   }

   TFile* outfile = TFile::Open("hists/l0h_sq5_Bd_errStat_asym.root","recreate");
   err->SetName("errorSym");
   err->Write();
   errp->SetName("errorHi");
   errp->Write();
   errm->SetName("errorLo");
   errm->Write();
   outfile->Close();
}
