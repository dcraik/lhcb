void fitKee_2D(Int_t bin) {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/BuKee_reduced_addMcorr.root");

	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	TH2D* hpr  = new TH2D("hpr", "",160,4900.,6500.,160,4900.,6500.);
	TH2D* hsig = new TH2D("hsig","",160,4900.,6500.,160,4900.,6500.);

	TCanvas c1;
	DecayTree->Draw("B_NoPsiFit_Mcorr:B_NoPsiFit_M>>hpr", "B_FullFit_M>4800&&B_FullFit_M<5200","colz");
	c1.SaveAs("pr.pdf");
	DecayTree->Draw("B_NoPsiFit_Mcorr:B_NoPsiFit_M>>hsig","B_FullFit_M>5200&&B_FullFit_M<5350","colz");
	c1.SaveAs("sig.pdf");

	TString binStr; binStr+=bin;
	Double_t minQ(0.), maxQ(0.);

	switch(bin) {
	//case 0:
	//	minQ = TMath::Sqrt(0.1e6);
	//	maxQ = TMath::Sqrt(0.98e6);
	//	break;
	//case 1:
	//	minQ = TMath::Sqrt(1.1e6);
	//	maxQ = TMath::Sqrt(2.e6);
	//	break;
	//case 2:
	//	minQ = TMath::Sqrt(2.e6);
	//	maxQ = TMath::Sqrt(3.e6);
	//	break;
	//case 3:
	//	minQ = TMath::Sqrt(3.e6);
	//	maxQ = TMath::Sqrt(4.e6);
	//	break;
	//case 4:
	//	minQ = TMath::Sqrt(4.e6);
	//	maxQ = TMath::Sqrt(5.e6);
	//	break;
	//case 5:
	//	minQ = TMath::Sqrt(5.e6);
	//	maxQ = TMath::Sqrt(6.e6);
	//	break;
	//case 6:
	//	minQ = TMath::Sqrt(6.e6);
	//	maxQ = TMath::Sqrt(7.e6);
	//	break;
	//case 7:
	//	minQ = TMath::Sqrt(7.e6);
	//	maxQ = TMath::Sqrt(8.e6);
	//	break;
	//case 8:
	//	minQ = TMath::Sqrt(11.e6);
	//	maxQ = TMath::Sqrt(11.75e6);
	//	break;
	//case 9:
	//	minQ = TMath::Sqrt(11.75e6);
	//	maxQ = TMath::Sqrt(12.5e6);
	//	break;
	//case 10:
	//	minQ = TMath::Sqrt(15.e6);
	//	maxQ = TMath::Sqrt(16.e6);
	//	break;
	//case 11:
	//	minQ = TMath::Sqrt(16.e6);
	//	maxQ = TMath::Sqrt(17.e6);
	//	break;
	//case 12:
	//	minQ = TMath::Sqrt(17.e6);
	//	maxQ = TMath::Sqrt(18.e6);
	//	break;
	//case 13:
	//	minQ = TMath::Sqrt(18.e6);
	//	maxQ = TMath::Sqrt(19.e6);
	//	break;
	//case 14:
	//	minQ = TMath::Sqrt(19.e6);
	//	maxQ = TMath::Sqrt(20.e6);
	//	break;
	//case 15:
	//	minQ = TMath::Sqrt(20.e6);
	//	maxQ = TMath::Sqrt(21.e6);
	//	break;
	//case 16:
	//	minQ = TMath::Sqrt(21.e6);
	//	maxQ = TMath::Sqrt(22.e6);
	//	break;
	case 17:
		minQ = TMath::Sqrt(1.1e6);
		maxQ = TMath::Sqrt(6.e6);
		break;
	case 18:
		minQ = TMath::Sqrt(15.e6);
		maxQ = TMath::Sqrt(22.e6);
		break;
	case 19:
		minQ = TMath::Sqrt(1.1e6);
		maxQ = TMath::Sqrt(22.e6);
		break;
	default:
		return;
	}
	TString cutStr("B_ConstBNoPsiFit_J_psi_1S_M> "); cutStr += minQ; cutStr += " && B_ConstBNoPsiFit_J_psi_1S_M< "; cutStr += maxQ;

	//B_M 
	RooRealVar B_M("B_NoPsiFit_M","; m(Kmumu) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar B_Mcorr("B_NoPsiFit_Mcorr","; m(Kmumu)_{corr} (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar Psi_M("B_ConstBNoPsiFit_J_psi_1S_M","; m(mumu) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",0,5000);
	RooRealVar B_PsiFit_M("B_FullFit_M","; m(Kmumu) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",3000,7000);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,B_Mcorr,Psi_M,B_PsiFit_M));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));
//	RooDataSet * dataPR  = dynamic_cast<RooDataSet*>(data->reduce("B_FullFit_M>4800&&B_FullFit_M<5200"));
//	RooDataSet * dataSig = dynamic_cast<RooDataSet*>(data->reduce("B_FullFit_M>5200&&B_FullFit_M<5350"));

	//backgrounds
	RooDataHist * dhpr  = new RooDataHist("dhpr",  "", RooArgSet(B_M,B_Mcorr), hpr);
	RooDataHist * dhsig = new RooDataHist("dhsig", "", RooArgSet(B_M,B_Mcorr), hsig);

	RooHistPdf pr("pr","", RooArgSet(B_M,B_Mcorr), *dhpr, 2);
	RooHistPdf sig("sig","", RooArgSet(B_M,B_Mcorr), *dhsig, 2);

	RooRealVar p0("p0","",-4.03708e-03,-0.1,0.1);
	RooRealVar p1("p1","",-4.03708e-03,-0.1,0.1);
	RooExponential comb0("comb0","",B_M,p0);
	
	RooRealVar mu("mu","",      5900.,5500.,6500.);
	RooRealVar sigma("sigma","", 800., 400.,1200.);
	RooGaussian comb1("comb1","",B_Mcorr,mu,sigma);

	RooProdPdf comb("comb","",RooArgList(comb0,comb1));

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",300,-1000,50000,"Events");
	RooRealVar nPR( "nPR" ,"#PR events"    ,300,-1000,50000,"Events");
	RooRealVar nbkg("nbkg","#signal events",300,-1000,50000,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr,comb), RooArgList(nsig,nPR,nbkg));
//	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr), RooArgList(nsig,nPR));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	B_M_RF_Plot = B_M.frame(20);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(K#mu#mu) (MeV/c^{2})");

	data1->plotOn(B_M_RF_Plot);
//	dhsig->plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("pr"),   RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("comb"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("sig"),  RooFit::LineStyle(kDashed));
	B_M_RF_Plot->Draw();

	can->SaveAs("plots/Kee_Q"+binStr+"_M.pdf");

//	can->SetLogy();
//	B_M_RF_Plot->SetMinimum(1.e-1);
//	B_M_RF_Plot->SetMaximum(5.e+2);
//	B_M_RF_Plot->Draw();
//	can->SaveAs("plots/Kee_Q"+binStr+"_M_log.pdf");

	B_Mc_RF_Plot = B_Mcorr.frame(20);
	B_Mc_RF_Plot->SetTitle("");
	B_Mc_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_Mc_RF_Plot->GetXaxis()->SetTitle("m(K#mu#mu)_{corr} (MeV/c^{2})");

	data1->plotOn(B_Mc_RF_Plot);
//	dhsig->plotOn(B_Mc_RF_Plot);
	full_RF_PDF.plotOn(B_Mc_RF_Plot);
	full_RF_PDF.plotOn(B_Mc_RF_Plot, RooFit::Components("pr"),   RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	full_RF_PDF.plotOn(B_Mc_RF_Plot, RooFit::Components("comb"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	full_RF_PDF.plotOn(B_Mc_RF_Plot, RooFit::Components("sig"),  RooFit::LineStyle(kDashed));
	B_Mc_RF_Plot->Draw();

	can->SaveAs("plots/Kee_Q"+binStr+"_Mc.pdf");
//
//	can->SetLogy();
//	B_M_RF_Plot->SetMinimum(1.e-1);
//	B_M_RF_Plot->SetMaximum(5.e+2);
//	B_M_RF_Plot->Draw();
//	can->SaveAs("plots/Kee_Q"+binStr+"_Mc_log.pdf");

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
