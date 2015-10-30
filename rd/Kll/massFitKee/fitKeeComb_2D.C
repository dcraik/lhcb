void fitKeeComb_2D() {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);//1111);

	gROOT->ProcessLine(".L RooMyPdf.cxx+");
	gROOT->ProcessLine(".L RooMyPdfB.cxx+");
	gROOT->ProcessLine(".L RooMyPdfC.cxx+");
	gROOT->ProcessLine(".L RooMyPdfD.cxx+");
	gROOT->ProcessLine(".L RooMyPdfE.cxx+");

	TFile * fileComb = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto.root");

	TTree * DecayTreeComb = dynamic_cast<TTree*>(fileComb->Get("DecayTree"));

	TH2D* hcomb = new TH2D("hcomb","",90,4800.,6600.,90,4800.,6600.);

	TCanvas c1;
	DecayTreeComb->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hcomb","B_M02_Subst0_e2pi > 2000 && JPs_M > 1000 && JPs_M < 2450 && NN < 0.8","colz");
	c1.SaveAs("comb.pdf");

	//B_M 
	RooRealVar B_M("B_DTF_PV_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar B_Mcorr("B_DTF_PV_Mcorr","; m(Kee)_{corr} (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar Psi_M("JPs_M","; m(ee) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",0,5000);
	RooRealVar B_PsiFit_M("B_DTF_PV_JPs_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",3000,7000);
	RooRealVar Ke_M("B_M02_Subst0_e2pi","; m(Ke) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",0,6000);
	RooRealVar NN("NN","",-1.0,1.0);

//	RooDataSet * dataComb = new RooDataSet("dataComb", "dataset with B_REFITTED_M", DecayTreeComb, RooArgSet(B_M,B_Mcorr,Psi_M,B_PsiFit_M,Ke_M,NN));
//	RooDataSet * dataComb1 = dynamic_cast<RooDataSet*>(data->reduce("B_M02_Subst0_e2pi > 2000 && JPs_M > 1000 && JPs_M < 2450 && NN < 0.8"));

	//backgrounds
	RooDataHist * dhcomb = new RooDataHist("dhcomb", "", RooArgSet(B_M,B_Mcorr), hcomb);

//	RooHistPdf comb("comb","", RooArgSet(B_M,B_Mcorr), *dhcomb, 2);

	RooRealVar p0("p0","", 1.63e-03,-0.1,0.1);
	RooRealVar p1("p1","", 1.63e-03,-0.1,0.1);
	
	RooRealVar a0("a0","", 0.,-10,10);
	RooRealVar a1("a1","", 0.,-10,10);
	RooRealVar a2("a2","", 0.,-10,10);
	RooRealVar b0("b0","", 0.,-10,10);
	RooRealVar b1("b1","", 0.,-10,10);
	RooRealVar b2("b2","", 0.,-10,10);

	RooRealVar mu("mu","",      6000.,5500.,6500.);
	RooRealVar sigma("sigma","", 600., 400.,1200.);
	
	RooRealVar muM("muM","",      5000.,4000.,6000.);
	RooRealVar sigmaM("sigmaM","", 600., 400.,1200.);
	RooRealVar muC("muC","",      6000.,5500.,6500.);
	RooRealVar sigmaC("sigmaC","", 600., 400.,1200.);
	
	RooRealVar muM1("muM1","",      5000.,4000.,6000.);
	RooRealVar sigmaM1("sigmaM1","", 600., 400.,1200.);
	RooRealVar muC1("muC1","",       500.,   0.,1000.);
	RooRealVar sigmaC1("sigmaC1","", 600., 400.,1200.);

	RooMyPdf  combA("combA","",B_M,B_Mcorr,p0,mu,sigma);
	RooMyPdfB combB("combB","",B_M,B_Mcorr,muM,sigmaM,muC,sigmaC);
	//RooMyPdfC combC("combC","",B_M,B_Mcorr,p0,mu,sigma);
	//RooMyPdfD combD("combD","",B_M,B_Mcorr,a0,a1,a2,b0,b1,b2);//,mu,sigma);
	RooMyPdfE combE("combE","",B_M,B_Mcorr,muM1,sigmaM1,muC1,sigmaC1);

	//fit comb parameters to events that fail the NN cut
	RooRealVar nbkg1("nbkg1","#signal events",60000,-1000,500000,"Events");
	RooRealVar nbkg2("nbkg2","#signal events",60000,-1000,500000,"Events");
	RooRealVar nbkg3("nbkg3","#signal events",60000,-1000,500000,"Events");
	RooAddPdf full_combA_PDF("full_combA_PDF","",RooArgList(combA), RooArgList(nbkg1));
	RooAddPdf full_combB_PDF("full_combB_PDF","",RooArgList(combB), RooArgList(nbkg2));
	RooAddPdf full_combC_PDF("full_combC_PDF","",RooArgList(combE), RooArgList(nbkg3));

	full_combA_PDF.fitTo(*dhcomb,RooFit::Extended());
	full_combB_PDF.fitTo(*dhcomb,RooFit::Extended());
	full_combC_PDF.fitTo(*dhcomb,RooFit::Extended());
        TH1D* hA = combA.createHistogram("hA",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hA->Draw("colz");
	c1.SaveAs("combBkgA.pdf");
        TH1D* hB = combB.createHistogram("hB",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hB->Draw("colz");
	c1.SaveAs("combBkgB.pdf");
        TH1D* hD = combE.createHistogram("hD",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hD->Draw("colz");
	c1.SaveAs("combBkgE.pdf");

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	
	B_M_comb_Plot = B_M.frame(20);
	B_M_comb_Plot->SetTitle("");
	B_M_comb_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_comb_Plot->GetXaxis()->SetTitle("m(Kee) (MeV/c^{2})");

	dhcomb->plotOn(B_M_comb_Plot);
	full_combA_PDF.plotOn(B_M_comb_Plot);
	full_combB_PDF.plotOn(B_M_comb_Plot, RooFit::LineColor(kRed));
	full_combC_PDF.plotOn(B_M_comb_Plot, RooFit::LineColor(kMagenta));
	B_M_comb_Plot->Draw();

	can->SaveAs("plots/Kee_comb_M.pdf");

	B_Mc_comb_Plot = B_Mcorr.frame(20);
	B_Mc_comb_Plot->SetTitle("");
	B_Mc_comb_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_Mc_comb_Plot->GetXaxis()->SetTitle("m(Kee)_{corr} (MeV/c^{2})");

	dhcomb->plotOn(B_Mc_comb_Plot);
	full_combA_PDF.plotOn(B_Mc_comb_Plot);
	full_combB_PDF.plotOn(B_Mc_comb_Plot, RooFit::LineColor(kRed));
	full_combC_PDF.plotOn(B_Mc_comb_Plot, RooFit::LineColor(kMagenta));
	B_Mc_comb_Plot->Draw();

	can->SaveAs("plots/Kee_comb_Mc.pdf");


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
