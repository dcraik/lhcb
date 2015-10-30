{
	gSystem.Load("libRooFit");
	gROOT.SetStyle("Plain");
	//gStyle.SetOptStat(1111);

	gStyle->SetOptStat(0000);
	gROOT->ProcessLine(".L ~/lhcb/lhcbStyle.C");
	lhcbStyle();

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/Kee/B2Kee_Strip21_data_folded_presel_JPsi_noPID_addMass.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));
	DecayTree->SetBranchStatus("*",0);
	DecayTree->SetBranchStatus("mKJpsi",1);
	DecayTree->SetBranchStatus("BDTKeeBig",1);
	DecayTree->SetBranchStatus("K_Kst_ProbNNk",1);
	DecayTree->SetBranchStatus("e_plus_ProbNNe",1);
	DecayTree->SetBranchStatus("e_minus_ProbNNe",1);
	DecayTree->SetBranchStatus("category",1);

	//B_M 
	RooRealVar B_M("mKJpsi","; m(K J/psi) (MeV); Candidates / 10 MeV", 5200,6200);
	RooRealVar B_M("mKJpsi","; m(K J/psi) (MeV); Candidates / 10 MeV", 5200,6200);
	RooRealVar BDT_KeeBig("BDTKeeBig","",0.,1.);
	RooRealVar K_Kst_ProbNNk("K_Kst_ProbNNk","",0.,1.);
	RooRealVar e_plus_ProbNNe("e_plus_ProbNNe","",0.,1.);
	RooRealVar e_minus_ProbNNe("e_minus_ProbNNe","",0.,1.);

	RooDataSet * data = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,BDT_KeeBig,K_Kst_ProbNNk,e_plus_ProbNNe,e_minus_ProbNNe));

	//1  #sigma_{B}   1.30451e+01   8.17659e-02   1.05141e-04  -1.30699e-01
	//2  M_{B}        5.28265e+03   5.32422e-02   2.66635e-04  -1.29827e-01
	//3  frac         5.14008e-01   6.55765e-03   2.78474e-04  -3.16961e+00
	//4  nbkg         1.48484e+05   5.17449e+02   1.29617e-04   1.92290e-01
	//5  nsig         1.68127e+05   5.36094e+02   1.42446e-04   3.55038e-01
	//6  p0          -2.42757e-03   1.33109e-05   8.39673e-08  -4.85513e-04
	//7  ratio        2.00000e+00   1.77712e-04   4.06499e-03   1.57086e+00
	//
	//
	//1  #sigma_{B}   9.74446e+00   1.24021e-01   2.66080e-03  -3.57965e-01
	//2  M_{B}        5.28192e+03   5.43310e-02   1.23448e-03  -1.49288e-01
	//3  a1           2.18708e+00   2.26963e-02   1.52773e-03  -5.97508e-01
	//4  a2          -7.89258e-01   1.34786e-02   1.48876e-03   1.00125e+00
	//5  frac         3.45165e-01   6.80455e-03   4.86091e-03  -3.14845e-01
	//6  n1           4.07486e-09   5.17795e-03   2.29921e-02  -1.57076e+00
	//7  n2           3.14713e+00   1.82263e-01   8.14130e-03  -3.79627e-01
	//8  nbkg         1.02246e+05   1.59357e+03   2.91632e-03  -1.78264e-01
	//9  nsig         2.14366e+05   1.62597e+03   5.35264e-03   7.98151e-01
	//10  p0         -1.36164e-03   3.41550e-05   2.54453e-06  -2.72329e-04
	//11  ratio       2.17495e+00   3.27107e-02   7.34430e-03   3.99386e-01
	// B DCB 
	// start, range to from. plus names and titles.
	RooRealVar sigmean("M_{B}","B mass",5282.7,5250.0,5325.0,"MeV");
	RooRealVar sigsigma("#sigma_{B}","B sigma",9.7,0.0,30.0,"MeV");
	RooRealVar a1("a1","a1",  2.2,  0.0, 10.0);
	RooRealVar n1("n1","n1",  0.0,  0.0, 10.0);
	RooRealVar a2("a2","a2", -0.8,-10.0,  0.0);
	RooRealVar n2("n2","n2",  3.1,  0.0, 10.0);
	RooRealVar ratio("ratio","Ratio of widths",2.17,0.3,3.);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",0.315,0.0,1.0);
	RooCBShape BSig_RF( "Bsig_RF", "Signal CB B RF Mass", B_M, sigmean, sigsigma, a1, n1 );
	RooCBShape BSig_RF2( "Bsig_RF2", "Signal CB B RF Mass", B_M, sigmean, sigsigma2, a2, n2 );
	//RooGaussian BSig_RF( "Bsig_RF", "Signal B RF Mass", B_M, sigmean, sigsigma);
	//RooGaussian BSig_RF2( "Bsig_RF2", "Signal B RF Mass", B_M, sigmean, sigsigma2);
	RooAddPdf BdSig("Bdsignal","signal pdf",RooArgList(BSig_RF,BSig_RF2),RooArgList(frac));
	
	//# Flat Background - Combinatorial
	//RooRealVar p0("p0","p0 of background",-0.0001,-0.001,-0.0);
	//RooPolynomial flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,RooArgList(p0));
	RooRealVar p0("p0","p0 of background",-0.00243,-5.,5.);
	RooExponential flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,p0);
	
	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events",200000,-1000,250000,"Events");
	RooRealVar nbkg("nbkg","#C background events",100000,-1000,250000,"Events");
	RooRealVar nbkgSB("nbkgSB","#C background events",10000,-1000,250000,"Events");
	
	RooAddPdf bkg_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(flatbkg_RF
				),
			RooArgList(nbkgSB
				));

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(BdSig
				,flatbkg_RF
				),
			RooArgList(nsig
				,nbkg
				));
	
	B_M.setRange("sideband",5600,6200);
	B_M.setRange("full",5200,6200);

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);

	bkg_RF_PDF->fitTo(*data,RooFit::Extended(),RooFit::Range("sideband"));
	p0.setConstant();
	nbkg.setVal(nbkgSB.getVal()*flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband")).getVal()/(Double_t)flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full")).getVal());
	nbkg.setConstant();
	full_RF_PDF->fitTo(*data,RooFit::Extended());
//	a1.setConstant();
//	a2.setConstant();
//	n1.setConstant();
//	n2.setConstant();
//	frac.setConstant();
//	ratio.setConstant();

	B_M_RF_Plot = B_M->frame(100);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 10 MeV");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(K J/#psi) (MeV)");
	//B_M_RF_Plot->GetYaxis()->SetTitleOffset(1.1);
	
	data->plotOn(B_M_RF_Plot);
	//full_RF_PDF->paramOn(B_M_RF_Plot, data, "", 0, "NELU" ,0.7,0.99,0.99 );
	//B_M_RF_Plot->getAttText()->SetTextSize(0.03);
	full_RF_PDF->plotOn(B_M_RF_Plot);
	cout << "Fit chi2:\t" << B_M_RF_Plot.chiSquare(14) << endl;
	full_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("flatbkg_RF"), RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	full_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("Bdsignal"), RooFit::LineStyle(2));
	B_M_RF_Plot->Draw();
	//printLHCb("R");
	
	can->SaveAs("fit.pdf");
	
	B_M_RF_Plot->SetMinimum(1.e-0);
	//B_M_RF_Plot->SetMaximum(1.e+3);
	can->SetLogy();
	B_M_RF_Plot->Draw();
	//printLHCb("R");
	can->SaveAs("fit_log.pdf");
	can->SetLogy(0);
			

	//RooProdPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(full_RF_PDF1,constrain_ratio,constrain_frac));
			
	Double_t bdtMin(-0.1), bdtMax(0.3);
	Double_t kMin(0.), kMax(0.3);
	Double_t eMin(0.), eMax(0.3);

	Double_t cutbdt(0.), cutk(0.), cute(0.);

	Int_t ni(5), nj(4), nk(4);

	Double_t S[ni][nj][nk];
	Double_t B[ni][nj][nk];

	for(Int_t i=0; i<ni; ++i) {
		for(Int_t j=0; j<nj; ++j) {
			for(Int_t k=0; k<nk; ++k) {

				TString cutStr="";
				
				cutbdt = bdtMin + (bdtMax - bdtMin) * static_cast<Double_t>(i)/(ni-1);
				cutk   = kMin   + (  kMax - kMin  ) * static_cast<Double_t>(j)/(nj-1);
				cute   = eMin   + (  eMax - eMin  ) * static_cast<Double_t>(k)/(nk-1);

				std::cout << cutbdt << "\t" << cutk << "\t" << cute << std::endl;

				cutStr = "BDTKeeBig>";          cutStr+= cutbdt;
				cutStr+= " && K_Kst_ProbNNk>";   cutStr+= cutk;
				cutStr+= " && e_plus_ProbNNe>";  cutStr+= cute;
				cutStr+= " && e_minus_ProbNNe>"; cutStr+= cute;

				RooDataSet* data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));

				p0.setConstant(kFALSE);
				bkg_RF_PDF->fitTo(*data,RooFit::Extended(),RooFit::Range("sideband"));
				p0.setConstant();
				nbkg.setVal(nbkgSB.getVal()*flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband")).getVal()/(Double_t)flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full")).getVal());
				nbkg.setConstant();
			//	full_RF_PDF->fitTo(*data,RooFit::Constrain(RooArgSet(ratio,frac)),RooFit::Extended());
				full_RF_PDF->fitTo(*data1,RooFit::Extended());
			
				B_M_RF_Plot = B_M->frame(100);
				B_M_RF_Plot->SetTitle("");
				B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 10 MeV");
				B_M_RF_Plot->GetXaxis()->SetTitle("m(K J/#psi) (MeV)");
				//B_M_RF_Plot->GetYaxis()->SetTitleOffset(1.1);
			
				data->plotOn(B_M_RF_Plot);
				//full_RF_PDF->paramOn(B_M_RF_Plot, data, "", 0, "NELU" ,0.7,0.99,0.99 );
				//B_M_RF_Plot->getAttText()->SetTextSize(0.03);
				full_RF_PDF->plotOn(B_M_RF_Plot);
				cout << "Fit chi2:\t" << B_M_RF_Plot.chiSquare(14) << endl;
				full_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("flatbkg_RF"), RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
				full_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("Bdsignal"), RooFit::LineStyle(2));
				B_M_RF_Plot->Draw();
				//printLHCb("R");
			
				TString fitName("fits/");
				fitName+=i; fitName+="_";
				fitName+=j; fitName+="_";
				fitName+=k;

				can->SaveAs(fitName+".pdf");
			
				B_M_RF_Plot->SetMinimum(1.e-0);
			//	B_M_RF_Plot->SetMaximum(1.e+3);
				can->SetLogy();
				B_M_RF_Plot->Draw();
			//	printLHCb("R");
				can->SaveAs(fitName+"_log.pdf");
				can->SetLogy(0);
			
				Double_t v25m = sigmean.getVal() - 2.5*sigsigma.getVal();
				Double_t v25p = sigmean.getVal() + 2.5*sigsigma.getVal();
			
				B_M.setRange("signal",v25m,v25p);
			
				double fBd2 = BdSig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal")).getVal();
			//	double fBd5 = BdSig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband")).getVal();
				double fBd0 = BdSig.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full")).getVal();
			
				double fComb2 = flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal")).getVal();
			//	double fComb5 = flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband")).getVal();
				double fComb0 = flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("full")).getVal();
			
				cout << "Mean: " << sigmean.getVal() << ", sigma: " << sigsigma.getVal() << endl;
				cout << "nominal"      << endl;
				cout << "Bd\t"         << nsig.getVal()*fBd2/fBd0        << endl;
				cout << "Comb\t"       << nbkg.getVal()*fComb2/fComb0    << endl;

				S[i][j][k] = nsig.getVal()*fBd2/fBd0;
				B[i][j][k] = nbkg.getVal()*fComb2/fComb0;

				delete data1;
			}
		}
	}
	for(Int_t i=0; i<ni; ++i) {
		for(Int_t j=0; j<nj; ++j) {
			for(Int_t k=0; k<nk; ++k) {
				std::cout << bdtMin + (bdtMax - bdtMin) * static_cast<Double_t>(i)/(ni-1)
				  << "\t" << kMin   + (  kMax - kMin  ) * static_cast<Double_t>(j)/(nj-1)
				  << "\t" << eMin   + (  eMax - eMin  ) * static_cast<Double_t>(k)/(nk-1)
				  << "\t" << S[i][j][k] << "\t" << B[i][j][k] << std::endl;
				//(4351./7933.)/(5.5e-7/(0.05971*1.026e-3))/(19627./60970.)
			}
		}
	}
}
