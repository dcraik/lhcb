void fitKmm_allQ_P_floatRatios(Int_t bin) {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);

	TFile * file = TFile::Open("/Home/dcraik/lhcb/rd/Kll/tuples/fromPatrick/Kmm.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("finalTree_KMuMu"));

	TString binStr; binStr+=bin;
	Double_t minQ(0.), maxQ(0.);

	switch(bin) {
	case 0:
		minQ = TMath::Sqrt(0.1e6);
		maxQ = TMath::Sqrt(0.98e6);
		break;
	case 1:
		minQ = TMath::Sqrt(1.1e6);
		maxQ = TMath::Sqrt(2.e6);
		break;
	case 2:
		minQ = TMath::Sqrt(2.e6);
		maxQ = TMath::Sqrt(3.e6);
		break;
	case 3:
		minQ = TMath::Sqrt(3.e6);
		maxQ = TMath::Sqrt(4.e6);
		break;
	case 4:
		minQ = TMath::Sqrt(4.e6);
		maxQ = TMath::Sqrt(5.e6);
		break;
	case 5:
		minQ = TMath::Sqrt(5.e6);
		maxQ = TMath::Sqrt(6.e6);
		break;
	case 6:
		minQ = TMath::Sqrt(6.e6);
		maxQ = TMath::Sqrt(7.e6);
		break;
	case 7:
		minQ = TMath::Sqrt(7.e6);
		maxQ = TMath::Sqrt(8.e6);
		break;
	case 8:
		minQ = TMath::Sqrt(11.e6);
		maxQ = TMath::Sqrt(11.75e6);
		break;
	case 9:
		minQ = TMath::Sqrt(11.75e6);
		maxQ = TMath::Sqrt(12.5e6);
		break;
	case 10:
		minQ = TMath::Sqrt(15.e6);
		maxQ = TMath::Sqrt(16.e6);
		break;
	case 11:
		minQ = TMath::Sqrt(16.e6);
		maxQ = TMath::Sqrt(17.e6);
		break;
	case 12:
		minQ = TMath::Sqrt(17.e6);
		maxQ = TMath::Sqrt(18.e6);
		break;
	case 13:
		minQ = TMath::Sqrt(18.e6);
		maxQ = TMath::Sqrt(19.e6);
		break;
	case 14:
		minQ = TMath::Sqrt(19.e6);
		maxQ = TMath::Sqrt(20.e6);
		break;
	case 15:
		minQ = TMath::Sqrt(20.e6);
		maxQ = TMath::Sqrt(21.e6);
		break;
	case 16:
		minQ = TMath::Sqrt(21.e6);
		maxQ = TMath::Sqrt(22.e6);
		break;
	case 17:
		minQ = TMath::Sqrt(1.1e6);
		maxQ = TMath::Sqrt(6.e6);
		break;
	case 18:
		minQ = TMath::Sqrt(15.e6);
		maxQ = TMath::Sqrt(22.e6);
		break;
	case 19:
		minQ = TMath::Sqrt(8.e6);
		maxQ = TMath::Sqrt(11.e6);
		break;
	case 20:
		minQ = TMath::Sqrt(12.5e6);
		maxQ = TMath::Sqrt(15.e6);
		break;
	default:
		return;
	}
	TString cutStr("Jpsi_M> "); cutStr += minQ; cutStr += " && Jpsi_M< "; cutStr += maxQ;

	//B_M 
	RooRealVar B_M("Bplus_M","; m(Kmumu) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",5170,5970);
	RooRealVar Psi_M("Jpsi_M","; m(mumu) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",0,5000);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,Psi_M));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));

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
	RooRealVar ratio("ratio","Ratio of widths",1.60407e+00, 1.0, 3.0);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",6.78672e-01, 0.0, 1.0);
	RooCBShape BSig_RF( "Bsig_RF", "Signal CB B RF Mass", B_M, sigmean, sigsigma, a1, n1 );
	RooCBShape BSig_RF2( "Bsig_RF2", "Signal CB B RF Mass", B_M, sigmean, sigsigma2, a2, n2 );
	RooAddPdf B0Sig("B0signal","signal pdf",RooArgList(BSig_RF,BSig_RF2),RooArgList(frac));

	RooRealVar p0("p0","",-6.44379e-02,-0.1,0.1);
	RooExponential comb_bkg("comb_bkg","",B_M,p0);

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",300,-1000,500000,"Events");
	RooRealVar nbkg("nbkg","#signal events",300,-1000,500000,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(B0Sig,comb_bkg), RooArgList(nsig,nbkg));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	B_M_RF_Plot = B_M.frame(50);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(K#mu#mu) (MeV/c^{2})");

	data1->plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("comb_bkg"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
        full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("B0signal"), RooFit::LineStyle(kDashed));
	B_M_RF_Plot->Draw();

	can->SaveAs("plots/fromPatrick/floatRatios/Kmm_Q"+binStr+".pdf");

	can->SetLogy();
        B_M_RF_Plot->SetMinimum(1.e-1);
        B_M_RF_Plot->SetMaximum(5.e+2);
	B_M_RF_Plot->Draw();
	can->SaveAs("plots/fromPatrick/floatRatios/Kmm_Q"+binStr+"_log.pdf");

	//Get integrals
	double mBdm = sigmean.getVal() - 2.5*(sigsigma.getVal());
	double mBdp = sigmean.getVal() + 2.5*(sigsigma.getVal());

	B_M.setRange("signal",mBdm,mBdp);
	B_M.setRange("sideband",5400,5970);
	B_M.setRange("full",5170,5970);

	double fsig1 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal"))->getVal();
	double fsig2 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband"))->getVal();
	double fsig0 = B0Sig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full"))->getVal();

	double fbkg1 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal"))->getVal();
	double fbkg2 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband"))->getVal();
	double fbkg0 = comb_bkg.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full"))->getVal();

	std::cout << std::endl;
	std::cout << sigmean.getVal() << "\t" << sigsigma.getVal() << std::endl << std::endl;
	std::cout << "\t\tsig\tbkg" << std::endl;
	std::cout << "window  \t" << nsig.getVal()*fsig1/fsig0 << "\t" << nbkg.getVal()*fbkg1/fbkg0 << std::endl;
	std::cout << "sideband\t" << nsig.getVal()*fsig2/fsig0 << "\t" << nbkg.getVal()*fbkg2/fbkg0 << std::endl;
	std::cout << std::endl;

	std::ofstream fout;
	fout.open("bkgParams/floatRatios/"+binStr+".dat");
	fout << sigmean.getVal() - 2.5*sigsigma.getVal() << "\t" << sigmean.getVal() + 2.5*sigsigma.getVal() << "\t" << fbkg1/fbkg2 << std::endl;
	fout.close();

	fout.open("yields/floatRatios/"+binStr+".dat");
	fout << nsig.getVal()*fsig1/fsig0 << "\t" << nsig.getError()*fsig1/fsig0 << "\t" << nbkg.getVal()*fbkg1/fbkg0 << "\t" << nbkg.getError()*fbkg1/fbkg0 << "\t"
	     << nsig.getVal()*fsig2/fsig0 << "\t" << nsig.getError()*fsig2/fsig0 << "\t" << nbkg.getVal()*fbkg2/fbkg0 << "\t" << nbkg.getError()*fbkg2/fbkg0 << std::endl;
	fout.close();

	//// Try splot stuff
	//// First set all parameters to constant except for yields
	sigmean.setConstant();
	sigsigma.setConstant();
	p0.setConstant();

	//RooStats::SPlot * sData = new RooStats::SPlot("sData","An SPlot",*data1, &full_RF_PDF, RooArgList(nsig,nbkg));
	//sData->GetSDataSet()->write("/Home/dcraik/Kll/tuples/fromPatrick/Kmm_Q"+binStr+"_sWeights.txt");

}
