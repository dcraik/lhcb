void fitKeeNew_2D(Int_t bin) {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);//1111);

	gROOT->ProcessLine(".L RooMyPdfB.cxx+");

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto_NN80_hop4k.root");
	TFile * fileMC = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/MC12_Bu2KEE_PreSel_addCorrMass_NNA_TrigVeto_NN80_hop4k.root");
	TFile * fileComb = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto.root");

	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));
	TTree * DecayTreeMC = dynamic_cast<TTree*>(fileMC->Get("DecayTree"));
	TTree * DecayTreeComb = dynamic_cast<TTree*>(fileComb->Get("DecayTree"));

	TH2D* hpr   = new TH2D("hpr",  "",90,4800.,6600.,90,4800.,6600.);
	TH2D* hsig  = new TH2D("hsig", "",90,4800.,6600.,90,4800.,6600.);
	TH2D* hcomb = new TH2D("hcomb","",90,4800.,6600.,90,4800.,6600.);

	TH1D* h1pr   = new TH1D("h1pr",  "",90,4800.,6600.);
	TH1D* h1sig  = new TH1D("h1sig", "",90,4800.,6600.);
	TH1D* h1comb = new TH1D("h1comb","",90,4800.,6600.);

	TCanvas c1;
	DecayTree->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hpr", "B_DTF_PV_JPs_Mmeas>4800&&B_DTF_PV_JPs_Mmeas<5200 && B_M02_Subst0_e2pi > 2000","colz");
	c1.SaveAs("pr.pdf");
	//DecayTree->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hsig","B_DTF_PV_JPs_Mmeas>5200&&B_DTF_PV_JPs_Mmeas<5350 && B_M02_Subst0_e2pi > 2000","colz");
	DecayTreeMC->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hsig","B_M02_Subst0_e2pi > 2000 && JPs_M**2. > 1.1e6 && JPs_M**2. < 6.0e6","colz");
	c1.SaveAs("sig.pdf");
	DecayTreeComb->Draw("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas>>hcomb","B_M02_Subst0_e2pi > 2000 && JPs_M > 1000 && JPs_M < 2450 && NN < 0.","colz");
	c1.SaveAs("comb.pdf");

	DecayTree->Draw("B_DTF_PV_Mmeas>>h1pr", "B_DTF_PV_JPs_Mmeas>4800&&B_DTF_PV_JPs_Mmeas<5200 && B_M02_Subst0_e2pi > 2000","colz");
	DecayTreeMC->Draw("B_DTF_PV_Mmeas>>h1sig","B_M02_Subst0_e2pi > 2000","colz");
	DecayTreeComb->Draw("B_DTF_PV_Mmeas>>h1comb","B_M02_Subst0_e2pi > 2000 && JPs_M > 1000 && JPs_M < 2450 && NN < 0.","colz");

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
	TString cutStr("JPs_M> "); cutStr += minQ; cutStr += " && JPs_M< "; cutStr += maxQ;
	cutStr+="&& B_M02_Subst0_e2pi > 2000";

	//B_M 
	RooRealVar B_M("B_DTF_PV_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar B_Mcorr("B_DTF_PV_Mcorr","; m(Kee)_{corr} (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar Psi_M("JPs_M","; m(ee) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",0,5000);
	RooRealVar B_PsiFit_M("B_DTF_PV_JPs_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",3000,7000);
	RooRealVar Ke_M("B_M02_Subst0_e2pi","; m(Ke) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",0,6000);
	RooRealVar NN("NN","",-1.0,1.0);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,B_Mcorr,Psi_M,B_PsiFit_M,Ke_M,NN));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));
//	RooDataSet * dataPR  = dynamic_cast<RooDataSet*>(data->reduce("B_DTF_PV_JPs_Mmeas>4800&&B_DTF_PV_JPs_Mmeas<5200"));
//	RooDataSet * dataSig = dynamic_cast<RooDataSet*>(data->reduce("B_DTF_PV_JPs_Mmeas>5200&&B_DTF_PV_JPs_Mmeas<5350"));
	RooDataSet * dataComb = new RooDataSet("dataComb", "dataset with B_REFITTED_M", DecayTreeComb, RooArgSet(B_M,B_Mcorr,Psi_M,B_PsiFit_M,Ke_M,NN));
	RooDataSet * dataComb1 = dynamic_cast<RooDataSet*>(data->reduce("B_M02_Subst0_e2pi > 2000 && JPs_M > 1000 && JPs_M < 2450 && NN < 0."));

	//backgrounds
	RooDataHist * dhpr   = new RooDataHist("dhpr",   "", RooArgSet(B_M,B_Mcorr), hpr);
	RooDataHist * dhsig  = new RooDataHist("dhsig",  "", RooArgSet(B_M,B_Mcorr), hsig);
	RooDataHist * dhcomb = new RooDataHist("dhcomb", "", RooArgSet(B_M,B_Mcorr), hcomb);

	RooDataHist * dh1pr   = new RooDataHist("dh1pr",   "", RooArgSet(B_M), h1pr);
	RooDataHist * dh1sig  = new RooDataHist("dh1sig",  "", RooArgSet(B_M), h1sig);
	RooDataHist * dh1comb = new RooDataHist("dh1comb", "", RooArgSet(B_M), h1comb);

	RooHistPdf pr(  "pr",  "", RooArgSet(B_M,B_Mcorr), *dhpr,   2);
	RooHistPdf sig( "sig", "", RooArgSet(B_M,B_Mcorr), *dhsig,  2);
	RooHistPdf comb("comb","", RooArgSet(B_M,B_Mcorr), *dhcomb, 2);

//	RooHistPdf pr(  "pr",  "", RooArgSet(B_M), *dh1pr,   2);
//	RooHistPdf sig( "sig", "", RooArgSet(B_M), *dh1sig,  2);
//	RooHistPdf comb("comb","", RooArgSet(B_M), *dh1comb, 2);

//	RooRealVar p0("p0","", 4.03708e-03,-0.1,0.1);
//	RooRealVar p1("p1","",-4.03708e-03,-0.1,0.1);
//	RooExponential comb0("comb0","",B_M,p0);
//	
//	RooRealVar mu("mu","",      5900.,5500.,6500.);
//	RooRealVar sigma("sigma","", 800., 400.,1200.);
//	RooGaussian comb1("comb1","",B_Mcorr,mu,sigma);
//
//	RooProdPdf comb("comb","",RooArgList(comb0,comb1));
//	RooMyPdf comb("comb","",B_M,B_Mcorr,p0,mu,sigma);
	RooRealVar muM("muM","",      5000.,4000.,6000.);
	RooRealVar sigmaM("sigmaM","", 600., 400.,1200.);
	RooRealVar muC("muC","",      6000.,5500.,6500.);
	RooRealVar sigmaC("sigmaC","", 600., 400.,1200.);

//	RooMyPdfB comb("comb","",B_M,B_Mcorr,muM,sigmaM,muC,sigmaC);
//
//	//fit comb parameters to events that fail the NN cut
//	RooRealVar nbkg1("nbkg1","#signal events",60000,-1000,500000,"Events");
//	RooAddPdf full_comb_PDF("full_comb_PDF","",RooArgList(comb), RooArgList(nbkg1));
//
//	full_comb_PDF.fitTo(*dataComb1,RooFit::Extended());
//        TH1D* hd = comb.createHistogram("hd",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
//	hd->Draw("colz");
//	c1.SaveAs("combBkg.pdf");
//
////	return;
//
//	muM.setConstant();
//	sigmaM.setConstant();
//	muC.setConstant();
//	sigmaC.setConstant();
////	p0.setConstant();

	//Now fit to data

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",282,-100,500,"Events");
	RooRealVar nPR( "nPR" ,"#PR events"    ,298,-100,500,"Events");
	RooRealVar nbkg("nbkg","#signal events",168,-100,500,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr,comb), RooArgList(nsig,nPR,nbkg));
//	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr), RooArgList(nsig,nPR));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	//TH1D* ha = comb.createHistogram("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas");
        TH1D* ha = comb.createHistogram("ha",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	ha->Draw("colz");
	c1.SaveAs("comb.pdf");

	//TH1D* hb = pr.createHistogram("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas");
        TH1D* hb = pr.createHistogram("ha",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hb->Draw("colz");
	c1.SaveAs("pr.pdf");

	//TH1D* hc = sig.createHistogram("B_DTF_PV_Mcorr:B_DTF_PV_Mmeas");
        TH1D* hc = sig.createHistogram("ha",B_Mcorr,RooFit::Binning(50,5000,6500), RooFit::YVar(B_M,RooFit::Binning(50,4800,5700)));
	hc->Draw("colz");
	c1.SaveAs("sig.pdf");

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	
//	B_M_comb_Plot = B_M.frame(20);
//	B_M_comb_Plot->SetTitle("");
//	B_M_comb_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
//	B_M_comb_Plot->GetXaxis()->SetTitle("m(Kee) (MeV/c^{2})");
//
//	dataComb1->plotOn(B_M_comb_Plot);
////	full_comb_PDF.plotOn(B_M_comb_Plot);
//	B_M_comb_Plot->Draw();
//
//	can->SaveAs("plots/Kee_Q"+binStr+"_M_comb.pdf");
//
//	B_Mc_comb_Plot = B_Mcorr.frame(20);
//	B_Mc_comb_Plot->SetTitle("");
//	B_Mc_comb_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
//	B_Mc_comb_Plot->GetXaxis()->SetTitle("m(Kee)_{corr} (MeV/c^{2})");
//
//	dataComb1->plotOn(B_Mc_comb_Plot);
////	full_comb_PDF.plotOn(B_Mc_comb_Plot);
//	B_Mc_comb_Plot->Draw();
//
//	can->SaveAs("plots/Kee_Q"+binStr+"_Mc_comb.pdf");

	B_M_RF_Plot = B_M.frame(30);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(Kee) (MeV/c^{2})");

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

	B_Mc_RF_Plot = B_Mcorr.frame(30);
	B_Mc_RF_Plot->SetTitle("");
	B_Mc_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_Mc_RF_Plot->GetXaxis()->SetTitle("m(Kee)_{corr} (MeV/c^{2})");

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
