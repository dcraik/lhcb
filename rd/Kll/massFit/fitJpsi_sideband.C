{
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromAlex/BuKmm.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	//B_M 
	RooRealVar B_M("B_M","; m(Kmumu) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",5150,5190);
	RooRealVar Psi_M("Psi_M","; m(mumu) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",500,5000);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,Psi_M));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce("Psi_M>3096.916-50&&Psi_M<3096.916+50"));

	// B DCB 
	// start, range to from. plus names and titles.
	RooRealVar sigmean("M_{B}","B mass",5281.0,5250.0,5300.0,"MeV/c^{2}");
	RooRealVar sigsigma("#sigma_{Lo}","B sigma",12.0,0.0,30.0,"MeV/c^{2}");
	RooRealVar a1("a1","a1", 1.85179e+00);
	RooRealVar n1("n1","n1", 1.30452e+00);
	RooRealVar a2("a2","a2",-2.17629e+00);
	RooRealVar n2("n2","n2", 3.45811e+00);
	RooRealVar ratio("ratio","Ratio of widths",1.756,1.,5.);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",0.91,0.0,1.0);
	RooCBShape BSig_RF( "Bsig_RF", "Signal CB B RF Mass", B_M, sigmean, sigsigma, a1, n1 );
	RooCBShape BSig_RF2( "Bsig_RF2", "Signal CB B RF Mass", B_M, sigmean, sigsigma2, a2, n2 );
	RooAddPdf B0Sig("B0signal","signal pdf",RooArgList(BSig_RF,BSig_RF2),RooArgList(frac));

	RooRealVar p0("p0","",-0.01,-0.1,0.1);
	RooExponential comb_bkg("comb_bkg","",B_M,p0);

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",     0,-1000,700000,"Events");
	RooRealVar nbkg("nbkg","#signal events", 10000,-1000,500000,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(B0Sig,comb_bkg), RooArgList(nsig,nbkg));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	B_M_RF_Plot = B_M.frame(100);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 5 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");

	data1->plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("comb_bkg"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
        full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("B0signal"), RooFit::LineStyle(kDashed));
	B_M_RF_Plot->Draw();

	can->SaveAs("plots/jpsiSB.pdf");

	can->SetLogy();
	B_M_RF_Plot->Draw();
	can->SaveAs("plots/jpsiSB_log.pdf");
}
