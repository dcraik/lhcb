void plot(RooRealVar* var, RooDataHist* dh, RooAbsPdf* pdf, TString plotName) {
	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	
	RooPlot* B_M_sig_Plot = var->frame(20);
	B_M_sig_Plot->SetTitle("");
	B_M_sig_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_sig_Plot->GetXaxis()->SetTitle("m(Kee) (MeV/c^{2})");

	dh->plotOn(B_M_sig_Plot);
	pdf->plotOn(B_M_sig_Plot);
	B_M_sig_Plot->Draw();

	can->SaveAs(plotName);
	
}

void fitKeeSig_2Dc() {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);//1111);

	gROOT->ProcessLine(".L RooMyPdf.cxx+");

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto_NN95_hop4k.root");
	TFile * fileMC = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/MC12_Bu2KEE_PreSel_addCorrMass_NNA_TrigVeto_NN95_hop4k.root");

	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));
	TTree * DecayTreeMC = dynamic_cast<TTree*>(fileMC->Get("DecayTree"));

	TH2D* hsig   = new TH2D("hsig",   "",90,4800.,6600.,90,4800.,6600.);

	TCanvas c1;
	//DecayTree->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hsig","B_DTF_PV_JPs_Mmeas>5200&&B_DTF_PV_JPs_Mmeas<5350 && B_M02_Subst0_e2pi > 2000","colz");
	DecayTreeMC->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hsig","B_M02_Subst0_e2pi > 2000","colz");
	c1.SaveAs("sig.pdf");

	//B_M 
	RooRealVar B_M("B_DTF_PV_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar B_Mcorr("B_DTF_PV_Mcorr","; m(Kee)_{corr} (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);


	RooFormulaVar Mdiff("Mdiff","@1-@0",RooArgSet(B_M,B_Mcorr));

//	RooDataSet * dataSig = dynamic_cast<RooDataSet*>(data->reduce("B_DTF_PV_JPs_Mmeas>5200&&B_DTF_PV_JPs_Mmeas<5350"));

	//backgrounds
	RooDataHist * dhsig    = new RooDataHist("dhsig", "", RooArgSet(B_M,B_Mcorr), hsig);

//	RooHistPdf sig( "sig", "", RooArgSet(B_M,B_Mcorr), *dhsig,  2);

	RooRealVar mu0("mu0","",      5300.,5000.,6000.);
	RooRealVar sigma0("sigma0","", 100.,  50.,1200.);
//	RooGaussian sig0("sig0","",B_M,mu0,sigma0);

	RooRealVar mu1("mu1","",       100.,  50.,200.);
	RooRealVar sigma1("sigma1","", 200.,  50.,1200.);
//	RooGaussian sig1("sig1","",Mdiff,mu1,sigma1);
//	RooLandau sig1("sig1","",Mdiff,mu1,sigma1);

	RooRealVar ratio0("ratio0","",1.756,1.0,3.0);
	RooProduct sigma0b("sigma0b","",RooArgSet(sigma0,ratio0));
	RooRealVar frac0("frac0","",0.796,0.0,1.0);

	RooRealVar a0a("a0a","a0a",  1.85179e+00,   0.0, 10.0);
	RooRealVar n0a("n0a","n0a",  1.30452e+00,   0.0, 10.0);
	RooRealVar a0b("a0b","a0b", -2.17629e+00, -10.0,  0.0);
	RooRealVar n0b("n0b","n0b",  3.45811e+00,   0.0, 10.0);
	RooCBShape sig0a("sig0a","",B_M, mu0, sigma0, a0a, n0a);
	RooCBShape sig0b("sig0b","",B_M, mu0, sigma0b, a0b, n0b);
	RooAddPdf sig0("sig0","",RooArgList(sig0a,sig0b),RooArgList(frac0));

	RooRealVar ka("ka","", 2., 0.1, 100.);
	RooRealVar kb("kb","", 2., 0.1, 100.);
	RooRealVar mu1b("mu1b","", 100., 50.,200.);

//	RooGaussian sig1("sig1", "", Mdiff, mu1, sigma0);
	RooLognormal sig1a("sig1a", "", Mdiff, mu1,  ka);
	RooLognormal sig1b("sig1b", "", Mdiff, mu1b, kb);
	RooGenericPdf cutoff("cutoff", "", "(@0 > 0)", RooArgSet(Mdiff));

//	RooRealVar ratio1("ratio1","",1.756,1.0,3.0);
//	RooProduct sigma1b("sigma1b","",RooArgSet(sigma1,ratio1));
	RooRealVar frac1("frac1","",0.796,0.0,1.0);
//
//	RooRealVar a1a("a1a","a1a",  1.85179e+00,   0.0, 10.0);
//	RooRealVar n1a("n1a","n1a",  1.30452e+00,   0.0, 10.0);
////	RooRealVar a1b("a1b","a1b", -2.17629e+00, -10.0,  0.0);
////	RooRealVar n1b("n1b","n1b",  3.45811e+00,   0.0, 10.0);
//	RooCBShape sig1a("sig1a","",Mdiff, mu1, sigma1, a1a, n1a);
//	RooGaussian sig1b("sig1b","",Mdiff,mu1, sigma1b);
////	RooCBShape sig1b("sig1b","",Mdiff, mu1, sigma1b, a1b, n1b);
	RooAddPdf sig1("sig1","",RooArgList(sig1a,sig1b),RooArgList(frac1));

	RooProdPdf sig("sig","",RooArgList(sig0,sig1,cutoff));

	//fit comb parameters to events that fail the NN cut
	RooRealVar nsig("nsig","#signal events",60000,-1000,500000,"Events");
	RooAddPdf Mm_sig_PDF("Mm_sig_PDF","",RooArgList(sig0), RooArgList(nsig));
	RooAddPdf Mc_sig_PDF("Mc_sig_PDF","",RooArgList(sig1), RooArgList(nsig));
	RooAddPdf full_sig_PDF("full_sig_PDF","",RooArgList(sig), RooArgList(nsig));

	Mm_sig_PDF.fitTo(*dhsig,RooFit::Extended());
	//Mc_sig_PDF.fitTo(*dhsig,RooFit::Extended());

	mu0.setConstant();
	sigma0.setConstant();
	ratio0.setConstant();
	frac0.setConstant();
	a0a.setConstant();
	n0a.setConstant();
	a0b.setConstant();
	n0b.setConstant();
//	a1a.setConstant();
//	n1a.setConstant();
////	a1b.setConstant();
////	n1b.setConstant();
	
	full_sig_PDF.fitTo(*dhsig,RooFit::Extended());
	//full_sig_PDF.fitTo(*dhMDiff,RooFit::Extended());
        TH1D* hd = sig.createHistogram("hd",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hd->Draw("colz");
	c1.SaveAs("sigFit.pdf");

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	
	B_M_sig_Plot = B_M.frame(20);
	B_M_sig_Plot->SetTitle("");
	B_M_sig_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_sig_Plot->GetXaxis()->SetTitle("m(Kee) (MeV/c^{2})");

	dhsig->plotOn(B_M_sig_Plot);
	full_sig_PDF.plotOn(B_M_sig_Plot);
	B_M_sig_Plot->Draw();

	can->SaveAs("plots/Kee_sig_M.pdf");

	B_Mc_sig_Plot = B_Mcorr.frame(20);
	B_Mc_sig_Plot->SetTitle("");
	B_Mc_sig_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_Mc_sig_Plot->GetXaxis()->SetTitle("m(Kee)_{corr} (MeV/c^{2})");

	dhsig->plotOn(B_Mc_sig_Plot);
	full_sig_PDF.plotOn(B_Mc_sig_Plot);
	B_Mc_sig_Plot->Draw();

	can->SaveAs("plots/Kee_sig_Mc.pdf");

	RooDataSet* toy = full_sig_PDF.generate(RooArgSet(B_M,B_Mcorr), 1000000);
	toy->write("toy.txt");

//	//Get integrals
//	double mBdm = sigmean.getVal() - 2.5*(sigsigma.getVal());
//	double mBdp = sigmean.getVal() + 2.5*(sigsigma.getVal());
//
//	B_M.setRange("signal",mBdm,mBdp);
//	B_M.setRange("sideband",5400,5970);
//	B_M.setRange("full",5170,5970);
//
//	double fsig1 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal"))->getVal();
//	double fsig2 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband"))->getVal();
//	double fsig0 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full"))->getVal();
//
//	double fbkg1 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal"))->getVal();
//	double fbkg2 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband"))->getVal();
//	double fbkg0 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full"))->getVal();
//
//	std::cout << std::endl;
//	std::cout << sigmean.getVal() << "\t" << sigsigma.getVal() << std::endl << std::endl;
//	std::cout << "\t\tsig\tbkg" << std::endl;
//	std::cout << "window  \t" << nsig.getVal()*fsig1/fsig0 << "\t" << nbkg.getVal()*fbkg1/fbkg0 << std::endl;
//	std::cout << "sideband\t" << nsig.getVal()*fsig2/fsig0 << "\t" << nbkg.getVal()*fbkg2/fbkg0 << std::endl;
//	std::cout << std::endl;
//
//	std::ofstream fout;
//	fout.open("bkgParams/"+binStr+".dat");
//	fout << sigmean.getVal() - 2.5*sigsigma.getVal() << "\t" << sigmean.getVal() + 2.5*sigsigma.getVal() << "\t" << fbkg1/fbkg2 << std::endl;
//	fout.close();
//
//	//// Try splot stuff
//	//// First set all parameters to constant except for yields
//	sigmean.setConstant();
//	sigsigma.setConstant();
//	p0.setConstant();
//
//	RooStats::SPlot * sData = new RooStats::SPlot("sData","An SPlot",*data1, &full_RF_PDF, RooArgList(nsig,nbkg));
//	sData->GetSDataSet()->write("/Home/dcraik/Kll/tuples/fromPatrick/Kmm_Q"+binStr+"_sWeights.txt");

}
