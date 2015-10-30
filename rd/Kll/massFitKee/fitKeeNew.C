void fitKeeNew(Int_t bin, Int_t NNmin) {
	gSystem->Load("libRooFit");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);

	TFile * file = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto_NN00_hop4k.root");

	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	TString binStr; binStr+=bin;
	TString NNStr;  NNStr+=NNmin;
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
	cutStr += "&& NN > 0."; cutStr+= NNmin;
	cutStr+="&& B_M02_Subst0_e2pi > 2000";

	//B_M 
	RooRealVar B_M("B_DTF_PV_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",4900,6500);
	RooRealVar Psi_M("JPs_M","; m(ee) (MeV/c^{2}); Candidates / 45 MeV/c^{2}",0,5000);
	RooRealVar B_PsiFit_M("B_DTF_PV_JPs_Mmeas","; m(Kee) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",3000,7000);
	RooRealVar Ke_M("B_M02_Subst0_e2pi","; m(Ke) (MeV/c^{2}); Candidates / 12 MeV/c^{2}",0,6000);

	RooRealVar NN("NN","",-1.0,1.0);

	RooDataSet * data  = new RooDataSet("data", "dataset with B_REFITTED_M", DecayTree, RooArgSet(B_M,Psi_M,B_PsiFit_M,Ke_M,NN));
	RooDataSet * data1 = dynamic_cast<RooDataSet*>(data->reduce(cutStr));
	RooDataSet * dataPR  = dynamic_cast<RooDataSet*>(data->reduce("B_DTF_PV_JPs_Mmeas>4800&&B_DTF_PV_JPs_Mmeas<5200"));
	RooDataSet * dataSig = dynamic_cast<RooDataSet*>(data->reduce("B_DTF_PV_JPs_Mmeas>5200&&B_DTF_PV_JPs_Mmeas<5350"));
	RooDataSet * dataSL  = dynamic_cast<RooDataSet*>(data->reduce(cutStr+"&&B_M02_Subst0_e2pi<2000"));

    	RooKeysPdf pr( "pr" ,"pr" ,B_M,*dataPR);
    	RooKeysPdf sig("sig","sig",B_M,*dataSig);
    	RooKeysPdf sl( "sl" ,"sl" ,B_M,*dataSL);

	RooRealVar p0("p0","",-4.03708e-03,-0.1,0.1);
	RooExponential comb("comb","",B_M,p0);

//	RooRealVar m0("m0","",5200.);//,5000.,5400.);//-4.03708e-03,-0.1,0.1);
//	RooRealVar c("c","",0.1);//,0.01,1.);//-4.03708e-03,-0.1,0.1);
//	RooArgusBG sl("sl","",B_M,m0,c);

	// Number of signal & background events
	RooRealVar nsig("nsig","#signal events"       ,300,-1000,50000,"Events");
	RooRealVar nPR( "nPR" ,"#PR events"           ,300,-1000,50000,"Events");
	RooRealVar nbkg("nbkg","#comb events"         ,300,-1000,50000,"Events");
	RooRealVar nSL( "nSL" ,"#semi-leptonic events",300,-1000,50000,"Events");

	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr,comb,sl), RooArgList(nsig,nPR,nbkg,nSL));
	RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig,pr,comb), RooArgList(nsig,nPR,nbkg));

	//# Do the fit on REFITTED Mass
	full_RF_PDF.fitTo(*data1,RooFit::Extended());

	TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
	B_M_RF_Plot = B_M.frame(20);
	B_M_RF_Plot->SetTitle("");
	B_M_RF_Plot->GetYaxis()->SetTitle("Candidates / 16 MeV/c^{2}");
	B_M_RF_Plot->GetXaxis()->SetTitle("m(K#mu#mu) (MeV/c^{2})");

	//dataSL->plotOn(B_M_RF_Plot);
	//B_M_RF_Plot->Draw();
	//can->SaveAs("tmp.pdf");

	data1->plotOn(B_M_RF_Plot);
	full_RF_PDF.plotOn(B_M_RF_Plot);
//	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("sl"),   RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+2));
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("pr"),   RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("comb"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
        full_RF_PDF.plotOn(B_M_RF_Plot, RooFit::Components("sig"),  RooFit::LineStyle(kDashed));
	B_M_RF_Plot->Draw();

	can->SaveAs("plots/Kee_Q"+binStr+"_NN"+NNStr+".pdf");

	can->SetLogy();
        B_M_RF_Plot->SetMinimum(1.e-1);
        B_M_RF_Plot->SetMaximum(5.e+2);
	B_M_RF_Plot->Draw();
	can->SaveAs("plots/Kee_Q"+binStr+"_NN"+NNStr+"_log.pdf");

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
