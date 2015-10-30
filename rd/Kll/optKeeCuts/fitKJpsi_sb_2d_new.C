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
	DecayTree->SetBranchStatus("B_plus_M",1);
	DecayTree->SetBranchStatus("mKJpsi",1);
	DecayTree->SetBranchStatus("BDTKeeBig",1);
	DecayTree->SetBranchStatus("K_Kst_ProbNNk",1);
	DecayTree->SetBranchStatus("K_Kst_PIDe",1);
	DecayTree->SetBranchStatus("e_plus_ProbNNe",1);
	DecayTree->SetBranchStatus("e_minus_ProbNNe",1);
	DecayTree->SetBranchStatus("category",1);

	//B_M 
	RooRealVar B_M("B_plus_M","; m(K J/psi) (MeV); Candidates / 10 MeV", 5100,6200);
	RooRealVar BDT_KeeBig("BDTKeeBig","",0.,1.);
	RooRealVar K_Kst_ProbNNk("K_Kst_ProbNNk","",0.,1.);
	RooRealVar K_Kst_PIDe("K_Kst_PIDe","",-100.,100.);
	RooRealVar e_plus_ProbNNe("e_plus_ProbNNe","",0.,1.);
	RooRealVar e_minus_ProbNNe("e_minus_ProbNNe","",0.,1.);

	RooDataSet * data = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,BDT_KeeBig,K_Kst_ProbNNk,e_plus_ProbNNe,e_minus_ProbNNe,K_Kst_PIDe));

	////# Flat Background - Combinatorial
	////RooRealVar p0("p0","p0 of background",-0.0001,-0.001,-0.0);
	////RooPolynomial flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,RooArgList(p0));
	//RooRealVar p0("p0","p0 of background",-0.00243,-5.,5.);
	//RooExponential flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,p0);
	//
	//// Number of signal & background events
	//RooRealVar nbkgSB("nbkgSB","#C background events",10000,-1000,250000,"Events");
	//
	//RooAddPdf bkg_RF_PDF("bkg_RF_PDF","RF PDF of everything",RooArgList(flatbkg_RF
	//			),
	//		RooArgList(nbkgSB
	//			));

	//B_M.setRange("signal",5200,5400);
	B_M.setRange("signal",5100,5400);
	B_M.setRange("sideband",5400,6200);

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);

	//bkg_RF_PDF->fitTo(*data,RooFit::Extended(),RooFit::Range("sideband"));

	//B_M_RF_Plot = B_M->frame(100);
	//B_M_RF_Plot->SetTitle("");
	//B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 10 MeV");
	//B_M_RF_Plot->GetXaxis()->SetTitle("m(K J/#psi) (MeV)");
	////B_M_RF_Plot->GetYaxis()->SetTitleOffset(1.1);
	//
	//data->plotOn(B_M_RF_Plot);
	////bkg_RF_PDF->paramOn(B_M_RF_Plot, data, "", 0, "NELU" ,0.7,0.99,0.99 );
	////B_M_RF_Plot->getAttText()->SetTextSize(0.03);
	//bkg_RF_PDF->plotOn(B_M_RF_Plot);
	//cout << "Fit chi2:\t" << B_M_RF_Plot.chiSquare(14) << endl;
	//bkg_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("flatbkg_RF"), RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	//bkg_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("Bdsignal"), RooFit::LineStyle(2));
	//B_M_RF_Plot->Draw();
	////printLHCb("R");
	//
	//can->SaveAs("fit.pdf");
	//
	//B_M_RF_Plot->SetMinimum(1.e-0);
	////B_M_RF_Plot->SetMaximum(1.e+3);
	//can->SetLogy();
	//B_M_RF_Plot->Draw();
	////printLHCb("R");
	//can->SaveAs("fit_log.pdf");
	//can->SetLogy(0);
			
	//1
	//Double_t bdtMin(-0.1), bdtMax(0.3);
	//Double_t kMin(0.), kMax(0.3);
	//Double_t eMin(0.), eMax(0.3);
	//Int_t ni(5), nj(4), nk(4);
	//2
	Double_t bdtMin(0.0), bdtMax(0.4);
	Double_t eMin(0.0), eMax(0.4);
	Int_t ni(21), nj(21);

	Double_t cutbdt(0.), cute(0.);

	Double_t S[ni][nj];
	Double_t B[ni][nj];

	for(Int_t i=0; i<ni; ++i) {
		for(Int_t j=0; j<nj; ++j) {

			TString cutStr="";
			
			cutbdt = bdtMin + (bdtMax - bdtMin) * static_cast<Double_t>(i)/(ni-1);
			cute   = eMin   + (  eMax - eMin  ) * static_cast<Double_t>(j)/(nj-1);

			std::cout << cutbdt << "\t" << cute << std::endl;

			cutStr = "K_Kst_ProbNNk>0.2 && K_Kst_PIDe<0 && BDTKeeBig>"; cutStr+= cutbdt;
			cutStr+= " && e_plus_ProbNNe>";  cutStr+= cute;
			cutStr+= " && e_minus_ProbNNe>"; cutStr+= cute;

			RooDataSet* data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));
	
			//# Flat Background - Combinatorial
			//RooRealVar p0("p0","p0 of background",-0.0001,-0.001,-0.0);
			//RooPolynomial flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,RooArgList(p0));
			RooRealVar p0("p0","p0 of background",-0.00243,-5.,5.);
			RooExponential flatbkg_RF("flatbkg_RF","flat background B RF Mass p.d.f",B_M,p0);
			
			// Number of signal & background events
			RooRealVar nbkgSB("nbkgSB","#C background events",10000,-1000,250000,"Events");
			
			RooAddPdf bkg_RF_PDF("bkg_RF_PDF","RF PDF of everything",RooArgList(flatbkg_RF
						),
					RooArgList(nbkgSB
						));


			bkg_RF_PDF->fitTo(*data1,RooFit::Extended(),RooFit::Range("sideband"));
		
			B_M_RF_Plot = B_M->frame(100);
			B_M_RF_Plot->SetTitle("");
			B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 10 MeV");
			B_M_RF_Plot->GetXaxis()->SetTitle("m(K J/#psi) (MeV)");
			//B_M_RF_Plot->GetYaxis()->SetTitleOffset(1.1);
		
			data1->plotOn(B_M_RF_Plot);
			//bkg_RF_PDF->paramOn(B_M_RF_Plot, data, "", 0, "NELU" ,0.7,0.99,0.99 );
			//B_M_RF_Plot->getAttText()->SetTextSize(0.03);
			bkg_RF_PDF->plotOn(B_M_RF_Plot);
			cout << "Fit chi2:\t" << B_M_RF_Plot.chiSquare(14) << endl;
			bkg_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("flatbkg_RF"), RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
			bkg_RF_PDF->plotOn(B_M_RF_Plot, RooFit::Components("Bdsignal"), RooFit::LineStyle(2));
			B_M_RF_Plot->Draw();
			//printLHCb("R");
		
			TString fitName("fits2DNew/");
			fitName+=i; fitName+="_";
			fitName+=j;

			can->SaveAs(fitName+".pdf");
		
			B_M_RF_Plot->SetMinimum(1.e-0);
		//	B_M_RF_Plot->SetMaximum(1.e+3);
			can->SetLogy();
			B_M_RF_Plot->Draw();
		//	printLHCb("R");
			can->SaveAs(fitName+"_log.pdf");
			can->SetLogy(0);

//			Double_t nsig = DecayTree->GetEntries(cutStr+"&&B_plus_M>5200&&B_plus_M<5400");
			Double_t nsig = DecayTree->GetEntries(cutStr+"&&B_plus_M>5100&&B_plus_M<5400&&mKJpsi>5200&&mKJpsi<5400");
		
			double fComb2 = flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("signal")).getVal();
			double fComb5 = flatbkg_RF.createIntegral(RooArgSet(B_M),RooFit::NormSet(B_M),RooFit::Range("sideband")).getVal();

			B[i][j] = nbkgSB.getVal()*fComb2/fComb5;
			S[i][j] = nsig;

			std::cout << S[i][j] << "\t" << B[i][j] << std::endl;

			delete data1;
		}
	}
	for(Int_t i=0; i<ni; ++i) {
		for(Int_t j=0; j<nj; ++j) {
				std::cout << bdtMin + (bdtMax - bdtMin) * static_cast<Double_t>(i)/(ni-1)
				  << "\t" << eMin   + (  eMax - eMin  ) * static_cast<Double_t>(j)/(nj-1)
				  << "\t" << S[i][j] << "\t" << B[i][j] << std::endl;
				//(4351./7933.)/(5.5e-7/(0.05971*1.026e-3))/(19627./60970.)
		}
	}
}
