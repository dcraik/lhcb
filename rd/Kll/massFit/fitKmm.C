{
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromAlex/BuKmm.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	//B_M 
	RooRealVar B_M("B_M","; m(Kmumu) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",5150,6000);
	RooRealVar Psi_M("Psi_M","; m(mumu) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",500,5000);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,Psi_M));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce("Psi_M> 1050 && Psi_M< 2450"));

// from J/Psi region
//   1  #sigma_{Lo}  1.59171e+01   9.61516e-02   1.80663e-03   6.11760e-02
//   2  M_{B}        5.28397e+03   3.00802e-02   1.66677e-03   3.66768e-01
//   3  a1           1.57752e+00   1.64484e-02   2.65338e-03  -7.53912e-01
//   4  a2          -2.64268e+00   2.11254e-02   2.51938e-03   4.90950e-01
//   5  frac         6.78672e-01   1.29969e-02   7.03329e-03   3.65422e-01
//   6  n1           4.79832e+00   2.84430e-01   2.61785e-02  -4.03463e-02
//   7  n2           1.08224e+00   2.68180e-02   5.47500e-03  -9.00362e-01
//   8  nbkg         5.56890e+03   1.31433e+02   7.62084e-03  -8.36640e-01
//   9  nsig         6.56230e+05   8.17224e+02   4.15943e-03   6.95832e-01
//  10  p0          -6.44379e-02   2.13769e-03   2.57927e-02   4.41139e-01
//  11  ratio        1.60407e+00   9.46569e-03   3.93086e-03  -7.72555e-01


	// B DCB 
	// start, range to from. plus names and titles.
	RooRealVar sigmean("M_{B}","B mass",5281.0,5250.0,5300.0,"MeV/c^{2}");
	RooRealVar sigsigma("#sigma_{Lo}","B sigma",15.9,0.0,30.0,"MeV/c^{2}");
	RooRealVar a1("a1","a1", 1.57752e+00);
	RooRealVar n1("n1","n1", 4.79832e+00);
	RooRealVar a2("a2","a2",-2.64268e+00);
	RooRealVar n2("n2","n2", 1.08224e+00);
	RooRealVar ratio("ratio","Ratio of widths",1.60407e+00);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",6.78672e-01);
	RooCBShape BSig_RF( "Bsig_RF", "Signal CB B RF Mass", B_M, sigmean, sigsigma, a1, n1 );
	RooCBShape BSig_RF2( "Bsig_RF2", "Signal CB B RF Mass", B_M, sigmean, sigsigma2, a2, n2 );
	RooAddPdf B0Sig("B0signal","signal pdf",RooArgList(BSig_RF,BSig_RF2),RooArgList(frac));

	RooRealVar p0("p0","",-6.44379e-02,-0.1,0.1);
	RooExponential comb_bkg("comb_bkg","",B_M,p0);

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",1000,-1000,50000,"Events");
	RooRealVar nbkg("nbkg","#signal events",1000,-1000,50000,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(B0Sig,comb_bkg), RooArgList(nsig,nbkg));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	B_M_RF_Plot = B_M.frame(100);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 8.5 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(K#mu#mu) (MeV/c^{2})");

	data1->plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("comb_bkg"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
        full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("B0signal"), RooFit::LineStyle(kDashed));
	B_M_RF_Plot->Draw();

	can->SaveAs("plots/Kmm_loQ.pdf");

	can->SetLogy();
        B_M_RF_Plot->SetMinimum(1.e-1);
        B_M_RF_Plot->SetMaximum(5.e+2);
	B_M_RF_Plot->Draw();
	can->SaveAs("plots/Kmm_loQ_log.pdf");

	//can->SetLogy(false);
	//B_M_RF_Plot->GetYaxis()->SetTitleOffset(1.5);
	//can->SetLeftMargin(0.15);
	//full_RF_PDF->paramOn(B_M_RF_Plot, data1, "", 0, "NELU" ,0.7,0.99,0.99 );
	//B_M_RF_Plot->getAttText()->SetTextSize(0.03);
	//B_M_RF_Plot->Draw();

	//can->SaveAs("fitMCtrue_DCB_params.pdf");
	//
      //// Try splot stuff
      //// First set all parameters to constant except for yields
      sigmean.setConstant();
      sigsigma.setConstant();
      p0.setConstant();
      
      RooStats::SPlot * sData = new RooStats::SPlot("sData","An SPlot",*data1, &full_RF_PDF, RooArgList(nsig,nbkg));
      sData->GetSDataSet()->write("~/Kll/tuples/Kmm_loQ_sWeights.txt");

}
