{
	gSystem.Load("libRooFit");
	gROOT.SetStyle("Plain");
//	gStyle.SetOptStat(1111);

	gROOT->ProcessLine(".L RooMyPdf.cxx+");

	gStyle->SetOptStat(0000);
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();

	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2KK_Scaled/B2D0Kpi_D02KK_selBd_Dsignal_Dst25_vetoes_PID3_NND2KK_addIMs.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	//B_D0_B_CM_M 
	RooRealVar B_D0_B_CM_M("B_D0_B_CM_M","",5100,5900);
	RooRealVar NN("NN","",-1.0,1.0);
	RooDataSet * data1 = new RooDataSet("data1", "", DecayTree, RooArgSet(B_D0_B_CM_M,NN));

	RooDataSet * adata = data1->reduce("NN>-0.80 && NN<=-0.08");
	RooDataSet * bdata = data1->reduce("NN>-0.08 && NN<= 0.43");
	RooDataSet * cdata = data1->reduce("NN> 0.43 && NN<= 0.59");
	RooDataSet * ddata = data1->reduce("NN> 0.59 && NN<= 0.67");
	RooDataSet * edata = data1->reduce("NN> 0.67 && NN<= 1.00");

	//Bd -> D* K pi
	//TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/Bd2DKpi/Bd2Dst0Kpi/Bd2DKpi_Bd2Dst0Kpi_P8_selBd_Dst25_NNf1_PID_misIDIMs.root");
	//TTree  * treeK1  = dynamic_cast<TTree*>(fileK1->Get("DecayTree"));
	//RooDataSet* dataK1 = new RooDataSet("dataK1","",treeK1,RooArgSet(B_D0_B_CM_M));
	//RooRealVar shift("shift","",-5.,-20.,-0.1);
	//RooMyPdf bd2dstkpi("bd2dstkpi","",B_D0_B_CM_M,shift,*dataK1,RooMyPdf::MirrorBoth,2);
	TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Bd2Dst0Kpi_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK1  = dynamic_cast<TH2D*>(fileK1->Get("Bd2Dst0Kpi"));
	TH1D  * histK1a = histK1->ProjectionX("Bd2D0stkpi_A",11,100);
	histK1a->Smooth(1);
	RooDataHist * dataK1 = new RooDataHist("dataK1", "", RooArgSet(B_D0_B_CM_M), histK1a);
	RooHistPdf bd2dstkpi("bd2dstkpi","", RooArgSet(B_D0_B_CM_M), *dataK1, 2);

	//// B0 -> D(*) pi pi
	TFile * fileK2 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Bd2D0pipi_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK2  = dynamic_cast<TH2D*>(fileK2->Get("Bd2D0pipi"));
	TH1D  * histK2a = histK2->ProjectionX("Bd2D0pipi_A",11,100);
	histK2a->Smooth(1);
	RooDataHist * dataK2 = new RooDataHist("dataK2", "", RooArgSet(B_D0_B_CM_M), histK2a);
	RooHistPdf bd2dpipi("bd2dpipi","", RooArgSet(B_D0_B_CM_M), *dataK2, 2);
	
	//// Lambda b -> D p pi
	TFile * fileK3 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Lb2D0ppi_smooth3_bins50_SDPweight_hist.root");
	TH2D  * histK3  = dynamic_cast<TH2D*>(fileK3->Get("Lb2D0ppi"));
	TH1D  * histK3a = histK3->ProjectionX("Lb2D0ppi_A",11,100);
	histK3a->Smooth(3);
	RooDataHist * dataK3 = new RooDataHist("dataK3", "", RooArgSet(B_D0_B_CM_M), histK3a);
	RooHistPdf lb2dppi("lb2dppi","", RooArgSet(B_D0_B_CM_M), *dataK3, 2);
	
	//// Lambda b -> D p K
	TFile * fileK4 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Lb2D0pK_smooth3_bins50_SDPweight_hist.root");
	TH2D  * histK4  = dynamic_cast<TH2D*>(fileK4->Get("Lb2D0pK"));
	TH1D  * histK4a = histK4->ProjectionX("Lb2D0pK_A",11,100);
	histK4a->Smooth(3);
	RooDataHist * dataK4 = new RooDataHist("dataK4", "", RooArgSet(B_D0_B_CM_M), histK4a);
	RooHistPdf lb2dpk("lb2dpk","", RooArgSet(B_D0_B_CM_M), *dataK4, 2);

	////Bd -> D K K
	TFile * fileK5 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Bd2D0KK_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK5  = dynamic_cast<TH2D*>(fileK5->Get("Bd2D0KK"));
	TH1D  * histK5a = histK5->ProjectionX("Bd2D0KK_A",11,100);
	histK5a->Smooth(1);
	RooDataHist * dataK5 = new RooDataHist("dataK5", "", RooArgSet(B_D0_B_CM_M), histK5a);
	RooHistPdf bd2dkk("bd2dkk","", RooArgSet(B_D0_B_CM_M), *dataK5, 2);

	//Bs -> D K K
	TFile * fileK6 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bs2D0KK_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK6  = dynamic_cast<TH2D*>(fileK6->Get("Bs2D0KK"));
	TH1D  * histK6a = histK6->ProjectionX("Bs2D0KK_A",11,100);
	histK6a->Smooth(1);
	RooDataHist * dataK6 = new RooDataHist("dataK6", "", RooArgSet(B_D0_B_CM_M), histK6a);
	RooHistPdf bs2dkk("bs2dkk","", RooArgSet(B_D0_B_CM_M), *dataK6, 2);

	////Bu -> D* K
	//TFile * fileK7 = TFile::Open("../../../reweightBackgrounds/DpipiBackgrounds/Bu2D0pi/B2Dpipi_Bu2Dstpi_selDpipi_hist_kill.root");
	//TH1F  * histK7  = dynamic_cast<TH1F*>(fileK7->Get("all"));
	//RooDataHist * dataK7 = new RooDataHist("dataK7", "", RooArgSet(B_D0_B_CM_M), histK7);
	//RooHistPdf bu2dstk("bu2dstk","", RooArgSet(B_D0_B_CM_M), *dataK7, 2);

	//Bs -> D* K pi
	//TFile * fileK8 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/Bd2DKpi/Bd2Dst0Kpi/Bd2DKpi_Bs2Dst0Kpi_P8_selBd_Dst25_NNf1_PID_misIDIMs.root");
	//TTree  * treeK8  = dynamic_cast<TTree*>(fileK8->Get("DecayTree"));
	//RooDataSet* dataK8 = new RooDataSet("dataK8","",treeK8,RooArgSet(B_D0_B_CM_M));
	//RooMyPdf bs2dstkpi("bs2dstkpi","",B_D0_B_CM_M,shift,*dataK8,RooMyPdf::MirrorBoth,2);
//	TFile * fileK7 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/hists/D2KK_Bs2Dst0Kpi_smooth1_bins50_SDPweight_hist.root");
//	TFile * fileK7 = TFile::Open("/data/lhcb/phrkbf/Bd2DKpiAnalysis/gammaDPBkgs/D2KK_Bs2DstKpi_smooth0_bins50_SDPweight_hist.root");
	TFile * fileK7 = TFile::Open("/data/lhcb/phrkbf/Bd2DKpiAnalysis/gammaDPBkgs/D2KK_Bs2DstKpi_smooth0_bins50_SDPweight_truth_hist.root");
	TH2D  * histK7  = dynamic_cast<TH2D*>(fileK7->Get("Bs2Dst0Kpi"));
	TH1D  * histK7a = histK7->ProjectionX("Bs2D0stkpi_A",11,100);
	histK7a->Smooth(1);
	RooDataHist * dataK7 = new RooDataHist("dataK7", "", RooArgSet(B_D0_B_CM_M), histK7a);
	RooHistPdf bs2dstkpi("bs2dstkpi","", RooArgSet(B_D0_B_CM_M), *dataK7, 2);

	// B Gaussian 
	// start, range to from. plus names and titles.
	RooRealVar sigmean("M_{B}","B mass",5281.0,5270.0,5290.0,"MeV/c^{2}");
	RooRealVar sigsigma("#sigma_{B}","B sigma",12.0,0.0,50.0,"MeV/c^{2}");
	RooRealVar ratio("ratio","Ratio of widths",1.756,1.0,3.0);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",0.796,0.0,1.0);

	RooRealVar a1("a1","",1.8);
	RooRealVar n1("n1","",1.2);
	RooRealVar a2("a2","",-2.1);
	RooRealVar n2("n2","",3.5);
	RooCBShape BSig_RF("Bsig_RF","",B_D0_B_CM_M, sigmean, sigsigma, a1, n1);
	RooCBShape BSig_RF2("Bsig_RF2","",B_D0_B_CM_M, sigmean, sigsigma2, a2, n2);

	RooAddPdf BSig("BSig","",RooArgList(BSig_RF,BSig_RF2),RooArgList(frac));

	RooGaussian constrain_ratio("cr","cr",ratio,RooFit::RooConst(1.756),RooFit::RooConst(0.048));
	RooGaussian constrain_frac("cf","cf",frac,RooFit::RooConst(0.796),RooFit::RooConst(0.018));

	RooRealVar BsMassDiff("M_{Bs}-M_{B}","Bs mass difference",87.19);
	RooAddition sigmeanBs("M_{Bs}","Bs mass",RooArgSet(sigmean,BsMassDiff));
	RooCBShape BsSig_RF( "Bssig_RF", "Signal CB Bs RF Mass", B_D0_B_CM_M, sigmeanBs, sigsigma, a1, n1 );
	RooCBShape BsSig_RF2( "Bssig_RF2", "Signal CB Bs RF Mass", B_D0_B_CM_M, sigmeanBs, sigsigma2, a2, n2 );
	RooAddPdf BsSig("BsSig","signal pdf",RooArgList(BsSig_RF,BsSig_RF2),RooArgList(frac));

	//# Flat Background - Combinatorial
	//RooRealVar p0("p0","p0 of background",-0.0001,-0.001,0.0);
	RooRealVar p0("p0","",-0.001,-0.1,0.1);
	//RooPolynomial comb_bkg1("comb_bkg1","flat background B RF Mass p.d.f",B_D0_B_CM_M,RooArgList(p0));
	RooExponential comb_bkg("comb_bkg","",B_D0_B_CM_M,p0);

//	RooRealVar combfrac("combfrac","",0.5,0.0,1.0);
//
//	RooAddPdf comb_bkg("comb_bkg","",RooArgList(comb_bkg1,bu2dstk),RooArgList(combfrac));

	// Number of signal & background events
	RooRealVar nbkg("nbkg","",7500,-1000,10000);
	RooRealVar nsig("nsig","",2000,-1000,6500);
	RooRealVar nsig2("nsig2","",5000,-1000,10000);
	RooRealVar nbkg2("nbkg2","",4000,-1000,8000);
	RooRealVar nbkg3("nbkg3","",1000,-1000,8000);
	RooRealVar nbkg4("nbkg4","",2000*0.370,-1000,1500);
	RooRealVar nbkg5("nbkg5","",2000*0.630,-1000,1500);
	RooRealVar nbkg6("nbkg6","",2000*0.089,-1000,1500);
	RooRealVar nbkg7("nbkg7","",2000*0.170,-1000,1500);
	RooRealVar nbkg8("nbkg8","",2000*0.067,-1000,1500);

	//// Constrain Dpipi background
	RooFormulaVar bd2dpipiratio("bd2dpipiratio","@0/@1",RooArgList(nbkg4,nsig));
	RooGaussian constrain_bd2dpipiratio("cbd2dpipi","cbd2dpipi",bd2dpipiratio,RooFit::RooConst(0.370),RooFit::RooConst(0.040));

	//// Constrain Dppi background
	RooFormulaVar lb2dppiratio("lb2dppiratio","@0/@1",RooArgList(nbkg5,nsig));
	RooGaussian constrain_lb2dppiratio("clb2dppi","clb2dppi",lb2dppiratio,RooFit::RooConst(0.630),RooFit::RooConst(0.140));

	//// Constrain DpK background
	RooFormulaVar lb2dpkratio("lb2dpkratio","@0/@1",RooArgList(nbkg6,nsig));
	RooGaussian constrain_lb2dpkratio("clb2dpk","clb2dpk",lb2dpkratio,RooFit::RooConst(0.089),RooFit::RooConst(0.023));
	
	//// Constrain DKK background
	RooFormulaVar bd2dkkratio("bd2dkkratio","@0/@1",RooArgList(nbkg7,nsig));
	RooGaussian constrain_bd2dkkratio("cbd2dkk","cbd2dkk",bd2dkkratio,RooFit::RooConst(0.170),RooFit::RooConst(0.040));
	
	//// Constrain BsDKK background
	RooFormulaVar bs2dkkratio("bs2dkkratio","@0/@1",RooArgList(nbkg8,nsig));
	RooGaussian constrain_bs2dkkratio("cbs2dkk","cbs2dkk",bs2dkkratio,RooFit::RooConst(0.067),RooFit::RooConst(0.029));


	//RooRealVar nbkg4("nbkg4","",2000*0.215,-1000,1500);
	//RooRealVar nbkg5("nbkg5","",2000*0.970,-1000,1500);
	//RooRealVar nbkg6("nbkg6","",2000*0.079,-1000,1500);
	//RooRealVar nbkg7("nbkg7","",2000*0.103,-1000,1500);
	//RooRealVar nbkg8("nbkg8","",2000*0.036,-1000,1500);

	////// Constrain Dpipi background
	//RooFormulaVar bd2dpipiratio("bd2dpipiratio","@0/@1",RooArgList(nbkg4,nsig));
	//RooGaussian constrain_bd2dpipiratio("cbd2dpipi","cbd2dpipi",bd2dpipiratio,RooFit::RooConst(0.215),RooFit::RooConst(0.021));

	////// Constrain Dppi background
	//RooFormulaVar lb2dppiratio("lb2dppiratio","@0/@1",RooArgList(nbkg5,nsig));
	//RooGaussian constrain_lb2dppiratio("clb2dppi","clb2dppi",lb2dppiratio,RooFit::RooConst(0.97),RooFit::RooConst(0.21));

	////// Constrain DpK background
	//RooFormulaVar lb2dpkratio("lb2dpkratio","@0/@1",RooArgList(nbkg6,nsig));
	//RooGaussian constrain_lb2dpkratio("clb2dpk","clb2dpk",lb2dpkratio,RooFit::RooConst(0.079),RooFit::RooConst(0.020));
	//
	////// Constrain DKK background
	//RooFormulaVar bd2dkkratio("bd2dkkratio","@0/@1",RooArgList(nbkg7,nsig));
	//RooGaussian constrain_bd2dkkratio("cbd2dkk","cbd2dkk",bd2dkkratio,RooFit::RooConst(0.103),RooFit::RooConst(0.025));
	//
	////// Constrain BsDKK background
	//RooFormulaVar bs2dkkratio("bs2dkkratio","@0/@1",RooArgList(nbkg8,nsig));
	//RooGaussian constrain_bs2dkkratio("cbs2dkk","cbs2dkk",bs2dkkratio,RooFit::RooConst(0.036),RooFit::RooConst(0.016));


	RooRealVar afrac("afrac","",0.2,0.0,1.0);
	RooRealVar bfrac("bfrac","",0.2,0.0,1.0);
	RooRealVar cfrac("cfrac","",0.2,0.0,1.0);
	RooRealVar dfrac("dfrac","",0.2,0.0,1.0);
	RooFormulaVar efrac("efrac","","1-afrac-bfrac-cfrac-dfrac",RooArgList(afrac,bfrac,cfrac,dfrac));

	RooRealVar aprfrac("aprfrac","",0.2,0.0,1.0);
	RooRealVar bprfrac("bprfrac","",0.2,0.0,1.0);
	RooRealVar cprfrac("cprfrac","",0.2,0.0,1.0);
	RooRealVar dprfrac("dprfrac","",0.2,0.0,1.0);
	RooFormulaVar eprfrac("eprfrac","","1-aprfrac-bprfrac-cprfrac-dprfrac",RooArgList(aprfrac,bprfrac,cprfrac,dprfrac));

	RooRealVar abkgfrac("abkgfrac","",0.7,0.0,1.0);
	RooRealVar bbkgfrac("bbkgfrac","",0.1,0.0,1.0);
	RooRealVar cbkgfrac("cbkgfrac","",0.05,0.0,1.0);
	RooFormulaVar dbkgfrac("dbkgfrac","","1-abkgfrac-bbkgfrac-cbkgfrac",RooArgList(abkgfrac,bbkgfrac,cbkgfrac));
	RooRealVar ebkgfrac("ebkgfrac","",0.0);
	//RooRealVar dbkgfrac("dbkgfrac","",0.01,0.0,1.0);
	//RooFormulaVar ebkgfrac("ebkgfrac","","1-abkgfrac-bbkgfrac-cbkgfrac-dbkgfrac",RooArgList(abkgfrac,bbkgfrac,cbkgfrac,dbkgfrac));

	RooFormulaVar anbkg("anbkg","abkgfrac*nbkg",RooArgList(abkgfrac,nbkg));
	RooFormulaVar bnbkg("bnbkg","bbkgfrac*nbkg",RooArgList(bbkgfrac,nbkg));
	RooFormulaVar cnbkg("cnbkg","cbkgfrac*nbkg",RooArgList(cbkgfrac,nbkg));
	RooFormulaVar dnbkg("dnbkg","dbkgfrac*nbkg",RooArgList(dbkgfrac,nbkg));
	RooFormulaVar enbkg("enbkg","ebkgfrac*nbkg",RooArgList(ebkgfrac,nbkg));

	RooFormulaVar ansig("ansig","","afrac*nsig",RooArgList(afrac,nsig));
	RooFormulaVar bnsig("bnsig","","bfrac*nsig",RooArgList(bfrac,nsig));
	RooFormulaVar cnsig("cnsig","","cfrac*nsig",RooArgList(cfrac,nsig));
	RooFormulaVar dnsig("dnsig","","dfrac*nsig",RooArgList(dfrac,nsig));
	RooFormulaVar ensig("ensig","","efrac*nsig",RooArgList(efrac,nsig));

	RooFormulaVar ansig2("ansig2","","afrac*nsig2",RooArgList(afrac,nsig2));
	RooFormulaVar bnsig2("bnsig2","","bfrac*nsig2",RooArgList(bfrac,nsig2));
	RooFormulaVar cnsig2("cnsig2","","cfrac*nsig2",RooArgList(cfrac,nsig2));
	RooFormulaVar dnsig2("dnsig2","","dfrac*nsig2",RooArgList(dfrac,nsig2));
	RooFormulaVar ensig2("ensig2","","efrac*nsig2",RooArgList(efrac,nsig2));

	RooFormulaVar anbkg2("anbkg2","","aprfrac*nbkg2",RooArgList(aprfrac,nbkg2));
	RooFormulaVar bnbkg2("bnbkg2","","bprfrac*nbkg2",RooArgList(bprfrac,nbkg2));
	RooFormulaVar cnbkg2("cnbkg2","","cprfrac*nbkg2",RooArgList(cprfrac,nbkg2));
	RooFormulaVar dnbkg2("dnbkg2","","dprfrac*nbkg2",RooArgList(dprfrac,nbkg2));
	RooFormulaVar enbkg2("enbkg2","","eprfrac*nbkg2",RooArgList(eprfrac,nbkg2));

	RooFormulaVar anbkg3("anbkg3","","aprfrac*nbkg3",RooArgList(aprfrac,nbkg3));
	RooFormulaVar bnbkg3("bnbkg3","","bprfrac*nbkg3",RooArgList(bprfrac,nbkg3));
	RooFormulaVar cnbkg3("cnbkg3","","cprfrac*nbkg3",RooArgList(cprfrac,nbkg3));
	RooFormulaVar dnbkg3("dnbkg3","","dprfrac*nbkg3",RooArgList(dprfrac,nbkg3));
	RooFormulaVar enbkg3("enbkg3","","eprfrac*nbkg3",RooArgList(eprfrac,nbkg3));

	RooFormulaVar anbkg4("anbkg4","","afrac*nbkg4",RooArgList(afrac,nbkg4));
	RooFormulaVar bnbkg4("bnbkg4","","bfrac*nbkg4",RooArgList(bfrac,nbkg4));
	RooFormulaVar cnbkg4("cnbkg4","","cfrac*nbkg4",RooArgList(cfrac,nbkg4));
	RooFormulaVar dnbkg4("dnbkg4","","dfrac*nbkg4",RooArgList(dfrac,nbkg4));
	RooFormulaVar enbkg4("enbkg4","","efrac*nbkg4",RooArgList(efrac,nbkg4));

	RooFormulaVar anbkg5("anbkg5","","afrac*nbkg5",RooArgList(afrac,nbkg5));
	RooFormulaVar bnbkg5("bnbkg5","","bfrac*nbkg5",RooArgList(bfrac,nbkg5));
	RooFormulaVar cnbkg5("cnbkg5","","cfrac*nbkg5",RooArgList(cfrac,nbkg5));
	RooFormulaVar dnbkg5("dnbkg5","","dfrac*nbkg5",RooArgList(dfrac,nbkg5));
	RooFormulaVar enbkg5("enbkg5","","efrac*nbkg5",RooArgList(efrac,nbkg5));

	RooFormulaVar anbkg6("anbkg6","","afrac*nbkg6",RooArgList(afrac,nbkg6));
	RooFormulaVar bnbkg6("bnbkg6","","bfrac*nbkg6",RooArgList(bfrac,nbkg6));
	RooFormulaVar cnbkg6("cnbkg6","","cfrac*nbkg6",RooArgList(cfrac,nbkg6));
	RooFormulaVar dnbkg6("dnbkg6","","dfrac*nbkg6",RooArgList(dfrac,nbkg6));
	RooFormulaVar enbkg6("enbkg6","","efrac*nbkg6",RooArgList(efrac,nbkg6));

	RooFormulaVar anbkg7("anbkg7","","afrac*nbkg7",RooArgList(afrac,nbkg7));
	RooFormulaVar bnbkg7("bnbkg7","","bfrac*nbkg7",RooArgList(bfrac,nbkg7));
	RooFormulaVar cnbkg7("cnbkg7","","cfrac*nbkg7",RooArgList(cfrac,nbkg7));
	RooFormulaVar dnbkg7("dnbkg7","","dfrac*nbkg7",RooArgList(dfrac,nbkg7));
	RooFormulaVar enbkg7("enbkg7","","efrac*nbkg7",RooArgList(efrac,nbkg7));

	RooFormulaVar anbkg8("anbkg8","","afrac*nbkg8",RooArgList(afrac,nbkg8));
	RooFormulaVar bnbkg8("bnbkg8","","bfrac*nbkg8",RooArgList(bfrac,nbkg8));
	RooFormulaVar cnbkg8("cnbkg8","","cfrac*nbkg8",RooArgList(cfrac,nbkg8));
	RooFormulaVar dnbkg8("dnbkg8","","dfrac*nbkg8",RooArgList(dfrac,nbkg8));
	RooFormulaVar enbkg8("enbkg8","","efrac*nbkg8",RooArgList(efrac,nbkg8));

	RooArgList pdfs(BSig,BsSig,comb_bkg,bd2dstkpi,bs2dstkpi,bd2dpipi,lb2dppi,lb2dpk,bd2dkk);
	pdfs.add(bs2dkk);
	RooArgList coeffsA(ansig,ansig2,anbkg,anbkg2,anbkg3,anbkg4,anbkg5,anbkg6,anbkg7);
	coeffsA.add(anbkg8);
	RooArgList coeffsB(bnsig,bnsig2,bnbkg,bnbkg2,bnbkg3,bnbkg4,bnbkg5,bnbkg6,bnbkg7);
	coeffsB.add(bnbkg8);
	RooArgList coeffsC(cnsig,cnsig2,cnbkg,cnbkg2,cnbkg3,cnbkg4,cnbkg5,cnbkg6,cnbkg7);
	coeffsC.add(cnbkg8);
	RooArgList coeffsD(dnsig,dnsig2,dnbkg,dnbkg2,dnbkg3,dnbkg4,dnbkg5,dnbkg6,dnbkg7);
	coeffsD.add(dnbkg8);
	RooArgList coeffsE(ensig,ensig2,enbkg,enbkg2,enbkg3,enbkg4,enbkg5,enbkg6,enbkg7);
	coeffsE.add(enbkg8);

	RooAddPdf afull_PDF("afull_PDF","",pdfs,coeffsA);
	RooAddPdf bfull_PDF("bfull_PDF","",pdfs,coeffsB);
	RooAddPdf cfull_PDF("cfull_PDF","",pdfs,coeffsC);
	RooAddPdf dfull_PDF("dfull_PDF","",pdfs,coeffsD);
	RooAddPdf efull_PDF("efull_PDF","",pdfs,coeffsE);

	// Simultaneous fit
	RooCategory cat("cat","cat");
	cat->defineType("cata");
	cat->defineType("catb");
	cat->defineType("catc");
	cat->defineType("catd");
	cat->defineType("cate");

	RooDataSet combo("combo","",B_D0_B_CM_M,RooFit::Index(cat),RooFit::Import("cata",*adata),RooFit::Import("catb",*bdata),RooFit::Import("catc",*cdata),RooFit::Import("catd",*ddata),RooFit::Import("cate",*edata));

	RooSimultaneous Simpdf("Simpdf","",cat);
	Simpdf->addPdf(afull_PDF,"cata");
	Simpdf->addPdf(bfull_PDF,"catb");
	Simpdf->addPdf(cfull_PDF,"catc");
	Simpdf->addPdf(dfull_PDF,"catd");
	Simpdf->addPdf(efull_PDF,"cate");

	RooProdPdf fullpdf("fullpdf","",RooArgList(Simpdf,constrain_ratio,constrain_frac,constrain_bd2dpipiratio,constrain_lb2dppiratio,constrain_lb2dpkratio,constrain_bd2dkkratio,constrain_bs2dkkratio));

	//# Do the fit on REFITTED Mass
	RooFitResult * r = fullpdf->fitTo(combo,RooFit::Extended(),RooFit::Save());

	TCanvas * can = new TCanvas("can","",2000,1500);
	can->Divide(3,2);
	can.cd(1);
	aPlot = B_D0_B_CM_M->frame(50);
	aPlot->SetTitle("");
	aPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	aPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	aPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(aPlot,RooFit::Cut("cat==cat::cata"));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 A:\t" << aPlot.chiSquare(0) << endl;
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	aPlot->Draw();
	can.cd(2);
	bPlot = B_D0_B_CM_M->frame(50);
	bPlot->SetTitle("");
	bPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	bPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	bPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(bPlot,RooFit::Cut("cat==cat::catb"));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 B:\t" << bPlot.chiSquare(0) << endl;
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	bPlot->Draw();
	can.cd(3);
	cPlot = B_D0_B_CM_M->frame(50);
	cPlot->SetTitle("");
	cPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	cPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	cPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(cPlot,RooFit::Cut("cat==cat::catc"));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 C:\t" << cPlot.chiSquare(0) << endl;
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	cPlot->Draw();
	can.cd(4);
	dPlot = B_D0_B_CM_M->frame(50);
	dPlot->SetTitle("");
	dPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	dPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	dPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(dPlot,RooFit::Cut("cat==cat::catd"));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 D:\t" << dPlot.chiSquare(0) << endl;
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	dPlot->Draw();
	can.cd(5);
	ePlot = B_D0_B_CM_M->frame(50);
	ePlot->SetTitle("");
	ePlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	ePlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	ePlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(ePlot,RooFit::Cut("cat==cat::cate"));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 E:\t" << ePlot.chiSquare(0) << endl;
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	ePlot->Draw();
	can.cd(6);
	fPlot = B_D0_B_CM_M->frame(50);
	fPlot->SetTitle("");
	fPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	fPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	fPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn( fPlot);
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),                RooFit::LineStyle(kDashed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),               RooFit::LineStyle(kDashed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),            RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi,bs2dstkpi)), RooFit::LineStyle(9), RooFit::LineColor(kRed));
//	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dstkpi)),           RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),            RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dppi)),             RooFit::LineStyle(3), RooFit::LineColor(kBlack));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),              RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),       RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bs2dkk)),              RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	fPlot->Draw();
	
	// Add the 3 sigma lines to the plot
	double xMin = 5100.;
	double xMax = 5900.;

	TLine* midLine = 0;
	TLine* uppLine = 0;
	TLine* lowLine = 0;

	midLine = new TLine( xMin,  0., xMax,  0. );
	uppLine = new TLine( xMin,  3., xMax,  3. );
	lowLine = new TLine( xMin, -3., xMax, -3. );

	uppLine->SetLineColor( kRed );
	lowLine->SetLineColor( kRed );
	midLine->SetLineColor( kBlue );
	midLine->SetLineStyle( kDashed );

	TCanvas * can2 = new TCanvas("can2","",2000,1500);
	can2->Divide(3,2);
	can2.cd(1);
	Plot2 = B_D0_B_CM_M->frame(50);
	Plot2->SetTitle("");
	Plot2->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(Plot2,RooFit::Cut("cat==cat::cata"));
	Simpdf->plotOn(Plot2, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo));
	Plot2->Draw();
	aPlot3 = B_D0_B_CM_M->frame(50);
	aPlot3->SetTitle("");
	aPlot3->SetMaximum( 5.);
	aPlot3->SetMinimum(-5.);
	aPlot3->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	aPlot3->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	aPlot3->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	aPlot3->addPlotable(Plot2->pullHist());
	aPlot3->addObject( lowLine );
	aPlot3->addObject( midLine );
	aPlot3->addObject( uppLine );
	aPlot3->Draw();

	can2.cd(2);
	Plot2 = B_D0_B_CM_M->frame(50);
	Plot2->SetTitle("");
	Plot2->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(Plot2,RooFit::Cut("cat==cat::catb"));
	Simpdf->plotOn(Plot2, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo));
	Plot2->Draw();
	bPlot3 = B_D0_B_CM_M->frame(50);
	bPlot3->SetTitle("");
	bPlot3->SetMaximum( 5.);
	bPlot3->SetMinimum(-5.);
	bPlot3->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	bPlot3->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	bPlot3->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	bPlot3->addPlotable(Plot2->pullHist());
	bPlot3->addObject( lowLine );
	bPlot3->addObject( midLine );
	bPlot3->addObject( uppLine );
	bPlot3->Draw();

	can2.cd(3);
	Plot2 = B_D0_B_CM_M->frame(50);
	Plot2->SetTitle("");
	Plot2->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(Plot2,RooFit::Cut("cat==cat::catc"));
	Simpdf->plotOn(Plot2, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo));
	Plot2->Draw();
	cPlot3 = B_D0_B_CM_M->frame(50);
	cPlot3->SetTitle("");
	cPlot3->SetMaximum( 5.);
	cPlot3->SetMinimum(-5.);
	cPlot3->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	cPlot3->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	cPlot3->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	cPlot3->addPlotable(Plot2->pullHist());
	cPlot3->addObject( lowLine );
	cPlot3->addObject( midLine );
	cPlot3->addObject( uppLine );
	cPlot3->Draw();

	can2.cd(4);
	Plot2 = B_D0_B_CM_M->frame(50);
	Plot2->SetTitle("");
	Plot2->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(Plot2,RooFit::Cut("cat==cat::catd"));
	Simpdf->plotOn(Plot2, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo));
	Plot2->Draw();
	dPlot3 = B_D0_B_CM_M->frame(50);
	dPlot3->SetTitle("");
	dPlot3->SetMaximum( 5.);
	dPlot3->SetMinimum(-5.);
	dPlot3->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	dPlot3->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	dPlot3->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	dPlot3->addPlotable(Plot2->pullHist());
	dPlot3->addObject( lowLine );
	dPlot3->addObject( midLine );
	dPlot3->addObject( uppLine );
	dPlot3->Draw();

	can2.cd(5);
	Plot2 = B_D0_B_CM_M->frame(50);
	Plot2->SetTitle("");
	Plot2->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	Plot2->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(Plot2,RooFit::Cut("cat==cat::cate"));
	Simpdf->plotOn(Plot2, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo));
	Plot2->Draw();
	ePlot3 = B_D0_B_CM_M->frame(50);
	ePlot3->SetTitle("");
	ePlot3->SetMaximum( 5.);
	ePlot3->SetMinimum(-5.);
	ePlot3->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	ePlot3->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	ePlot3->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	ePlot3->addPlotable(Plot2->pullHist());
	ePlot3->addObject( lowLine );
	ePlot3->addObject( midLine );
	ePlot3->addObject( uppLine );
	ePlot3->Draw();

	can->SaveAs("D2KKfit_fit.C");
	can2->SaveAs("D2KKfit_pull.C");

	uppLineA = new TLine( sigmean.getVal() + 2.5*sigsigma.getVal(), aPlot->GetMinimum(), 
	        	      sigmean.getVal() + 2.5*sigsigma.getVal(), aPlot->GetMaximum() );
	lowLineA = new TLine( sigmean.getVal() - 2.5*sigsigma.getVal(), aPlot->GetMinimum(), 
			      sigmean.getVal() - 2.5*sigsigma.getVal(), aPlot->GetMaximum() );

	uppLineB = new TLine( sigmean.getVal() + 2.5*sigsigma.getVal(), bPlot->GetMinimum(), 
	        	      sigmean.getVal() + 2.5*sigsigma.getVal(), bPlot->GetMaximum() );
	lowLineB = new TLine( sigmean.getVal() - 2.5*sigsigma.getVal(), bPlot->GetMinimum(), 
			      sigmean.getVal() - 2.5*sigsigma.getVal(), bPlot->GetMaximum() );

	uppLineC = new TLine( sigmean.getVal() + 2.5*sigsigma.getVal(), cPlot->GetMinimum(), 
	        	      sigmean.getVal() + 2.5*sigsigma.getVal(), cPlot->GetMaximum() );
	lowLineC = new TLine( sigmean.getVal() - 2.5*sigsigma.getVal(), cPlot->GetMinimum(), 
			      sigmean.getVal() - 2.5*sigsigma.getVal(), cPlot->GetMaximum() );

	uppLineD = new TLine( sigmean.getVal() + 2.5*sigsigma.getVal(), dPlot->GetMinimum(), 
	        	      sigmean.getVal() + 2.5*sigsigma.getVal(), dPlot->GetMaximum() );
	lowLineD = new TLine( sigmean.getVal() - 2.5*sigsigma.getVal(), dPlot->GetMinimum(), 
			      sigmean.getVal() - 2.5*sigsigma.getVal(), dPlot->GetMaximum() );

	uppLineE = new TLine( sigmean.getVal() + 2.5*sigsigma.getVal(), ePlot->GetMinimum(), 
	        	      sigmean.getVal() + 2.5*sigsigma.getVal(), ePlot->GetMaximum() );
	lowLineE = new TLine( sigmean.getVal() - 2.5*sigsigma.getVal(), ePlot->GetMinimum(), 
			      sigmean.getVal() - 2.5*sigsigma.getVal(), ePlot->GetMaximum() );

	uppLineA->SetLineStyle( kDotted );
	lowLineA->SetLineStyle( kDotted );
	uppLineB->SetLineStyle( kDotted );
	lowLineB->SetLineStyle( kDotted );
	uppLineC->SetLineStyle( kDotted );
	lowLineC->SetLineStyle( kDotted );
	uppLineD->SetLineStyle( kDotted );
	lowLineD->SetLineStyle( kDotted );
	uppLineE->SetLineStyle( kDotted );
	lowLineE->SetLineStyle( kDotted );

	aPlot4 = dynamic_cast<RooPlot*>(aPlot->Clone("aPlot4"));
	bPlot4 = dynamic_cast<RooPlot*>(bPlot->Clone("bPlot4"));
	cPlot4 = dynamic_cast<RooPlot*>(cPlot->Clone("cPlot4"));
	dPlot4 = dynamic_cast<RooPlot*>(dPlot->Clone("dPlot4"));
	ePlot4 = dynamic_cast<RooPlot*>(ePlot->Clone("ePlot4"));

	aPlot4->addObject( uppLineA );
	aPlot4->addObject( lowLineA );
	bPlot4->addObject( uppLineB );
	bPlot4->addObject( lowLineB );
	cPlot4->addObject( uppLineC );
	cPlot4->addObject( lowLineC );
	dPlot4->addObject( uppLineD );
	dPlot4->addObject( lowLineD );
	ePlot4->addObject( uppLineE );
	ePlot4->addObject( lowLineE );

	TCanvas * temp = new TCanvas("temp","");
	aPlot->Draw(); temp->SaveAs("figs/D2Kpifit_1.pdf"); temp->SaveAs("figs/D2Kpifit_1.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_1_log.pdf"); temp->SetLogy(0);
	TCanvas * temp = new TCanvas("temp","");
	aPlot->Draw(); temp->SaveAs("figs/D2KKfit_1.pdf"); temp->SaveAs("figs/D2KKfit_1.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_1_log.pdf"); temp->SetLogy(0);
	bPlot->Draw(); printLHCb("R","Other","LHCb (b)"); temp->SaveAs("figs/D2KKfit_2.pdf"); temp->SaveAs("figs/D2KKfit_2.C");
	temp->SaveAs("figs/fig7b.pdf");
	temp->SaveAs("figs/fig7b.png");
	temp->SaveAs("figs/fig7b.eps");
	temp->SaveAs("figs/fig7b.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_2_log.pdf"); temp->SetLogy(0);
	cPlot->Draw(); printLHCb("R","Other","LHCb (c)"); temp->SaveAs("figs/D2KKfit_3.pdf"); temp->SaveAs("figs/D2KKfit_3.C");
	temp->SaveAs("figs/fig7c.pdf");
	temp->SaveAs("figs/fig7c.png");
	temp->SaveAs("figs/fig7c.eps");
	temp->SaveAs("figs/fig7c.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_3_log.pdf"); temp->SetLogy(0);
	dPlot->Draw(); printLHCb("R","Other","LHCb (d)"); temp->SaveAs("figs/D2KKfit_4.pdf"); temp->SaveAs("figs/D2KKfit_4.C");
	temp->SaveAs("figs/fig7d.pdf");
	temp->SaveAs("figs/fig7d.png");
	temp->SaveAs("figs/fig7d.eps");
	temp->SaveAs("figs/fig7d.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_4_log.pdf"); temp->SetLogy(0);
	ePlot->Draw(); printLHCb("R","Other","LHCb (e)"); temp->SaveAs("figs/D2KKfit_5.pdf"); temp->SaveAs("figs/D2KKfit_5.C");
	temp->SaveAs("figs/fig7e.pdf");
	temp->SaveAs("figs/fig7e.png");
	temp->SaveAs("figs/fig7e.eps");
	temp->SaveAs("figs/fig7e.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_5_log.pdf"); temp->SetLogy(0);
	fPlot->Draw(); temp->SaveAs("figs/D2KKfit_1-5.pdf"); temp->SaveAs("figs/D2KKfit_1-5.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2KKfit_1-5_log.pdf"); temp->SetLogy(0);
	aPlot3->Draw(); temp->SaveAs("figs/D2KKfit_1_pull.pdf"); temp->SaveAs("figs/D2KKfit_1_pull.C");
	bPlot3->Draw(); temp->SaveAs("figs/D2KKfit_2_pull.pdf"); temp->SaveAs("figs/D2KKfit_2_pull.C");
	cPlot3->Draw(); temp->SaveAs("figs/D2KKfit_3_pull.pdf"); temp->SaveAs("figs/D2KKfit_3_pull.C");
	dPlot3->Draw(); temp->SaveAs("figs/D2KKfit_4_pull.pdf"); temp->SaveAs("figs/D2KKfit_4_pull.C");
	ePlot3->Draw(); temp->SaveAs("figs/D2KKfit_5_pull.pdf"); temp->SaveAs("figs/D2KKfit_5_pull.C");
	aPlot4->Draw(); printLHCb("R","Other","LHCb (a)"); temp->SaveAs("figs/D2KKfit_1_window.pdf"); temp->SaveAs("figs/D2KKfit_1_window.C");
	temp->SaveAs("figs/fig7a.pdf");
	temp->SaveAs("figs/fig7a.png");
	temp->SaveAs("figs/fig7a.eps");
	temp->SaveAs("figs/fig7a.C");
	bPlot4->Draw(); temp->SaveAs("figs/D2KKfit_2_window.pdf"); temp->SaveAs("figs/D2KKfit_2_window.C");
	cPlot4->Draw(); temp->SaveAs("figs/D2KKfit_3_window.pdf"); temp->SaveAs("figs/D2KKfit_3_window.C");
	dPlot4->Draw(); temp->SaveAs("figs/D2KKfit_4_window.pdf"); temp->SaveAs("figs/D2KKfit_4_window.C");
	ePlot4->Draw(); temp->SaveAs("figs/D2KKfit_5_window.pdf"); temp->SaveAs("figs/D2KKfit_5_window.C");

	double mBdm = sigmean.getVal() - 2.5*(sigsigma.getVal());
	double mBdp = sigmean.getVal() + 2.5*(sigsigma.getVal());
	double mBsm = sigmeanBs.getVal() - 2.5*(sigsigma.getVal());
	double mBsp = sigmeanBs.getVal() + 2.5*(sigsigma.getVal());

	B_D0_B_CM_M.setRange("Bd",mBdm,mBdp);
	B_D0_B_CM_M.setRange("Bs",mBsm,mBsp);
	B_D0_B_CM_M.setRange("full",5100,5900);
	B_D0_B_CM_M.setRange("sideband",5450,5900);

	double fBd1 = BSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBd2 = BSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBd3 = BSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBd0 = BSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBs1 = BsSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBs2 = BsSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBs3 = BsSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBs0 = BsSig.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fComb1 = comb_bkg.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fComb2 = comb_bkg.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fComb3 = comb_bkg.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fComb0 = comb_bkg.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBd2DstKpi1 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBd2DstKpi2 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBd2DstKpi3 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBd2DstKpi0 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBs2DstKpi1 = bs2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBs2DstKpi2 = bs2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBs2DstKpi3 = bs2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBs2DstKpi0 = bs2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBd2Dpipi1 = bd2dpipi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBd2Dpipi2 = bd2dpipi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBd2Dpipi3 = bd2dpipi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBd2Dpipi0 = bd2dpipi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBd2DKK1 = bd2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBd2DKK2 = bd2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBd2DKK3 = bd2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBd2DKK0 = bd2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBs2DKK1 = bs2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBs2DKK2 = bs2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBs2DKK3 = bs2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBs2DKK0 = bs2dkk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fLb2Dppi1 = lb2dppi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fLb2Dppi2 = lb2dppi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fLb2Dppi3 = lb2dppi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fLb2Dppi0 = lb2dppi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fLb2DpK1 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fLb2DpK2 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fLb2DpK3 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fLb2DpK0 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	cout << "Mean: " << sigmean.getVal() << ", sigma: " << sigsigma.getVal() << endl;
	cout << endl;

	Double_t nBd1        = nsig.getVal()*fBd1/fBd0;               
        Double_t nBs1        = nsig2.getVal()*fBs1/fBs0;              
        Double_t nComb1      = nbkg.getVal()*fComb1/fComb0;           
        Double_t nBd2DstKpi1 = nbkg2.getVal()*fBd2DstKpi1/fBd2DstKpi0;
        Double_t nBs2DstKpi1 = nbkg3.getVal()*fBs2DstKpi1/fBs2DstKpi0;
        Double_t nBd2Dpipi1  = nbkg4.getVal()*fBd2Dpipi1/fBd2Dpipi0;  
        Double_t nLb2Dppi1   = nbkg5.getVal()*fLb2Dppi1/fLb2Dppi0;      
        Double_t nLb2DpK1    = nbkg6.getVal()*fLb2DpK1/fLb2DpK0;      
        Double_t nBd2DKK1    = nbkg7.getVal()*fBd2DKK1/fBd2DKK0;      
        Double_t nBs2DKK1    = nbkg8.getVal()*fBs2DKK1/fBs2DKK0;      

	Double_t nBd2        = nsig.getVal()*fBd2/fBd0;               
        Double_t nBs2        = nsig2.getVal()*fBs2/fBs0;              
        Double_t nComb2      = nbkg.getVal()*fComb2/fComb0;           
        Double_t nBd2DstKpi2 = nbkg2.getVal()*fBd2DstKpi2/fBd2DstKpi0;
        Double_t nBs2DstKpi2 = nbkg3.getVal()*fBs2DstKpi2/fBs2DstKpi0;
        Double_t nBd2Dpipi2  = nbkg4.getVal()*fBd2Dpipi2/fBd2Dpipi0;  
        Double_t nLb2Dppi2   = nbkg5.getVal()*fLb2Dppi2/fLb2Dppi0;      
        Double_t nLb2DpK2    = nbkg6.getVal()*fLb2DpK2/fLb2DpK0;      
        Double_t nBd2DKK2    = nbkg7.getVal()*fBd2DKK2/fBd2DKK0;      
        Double_t nBs2DKK2    = nbkg8.getVal()*fBs2DKK2/fBs2DKK0;      

	Double_t nBd3        = nsig.getVal()*fBd3/fBd0;               
        Double_t nBs3        = nsig2.getVal()*fBs3/fBs0;              
        Double_t nComb3      = nbkg.getVal()*fComb3/fComb0;           
        Double_t nBd2DstKpi3 = nbkg2.getVal()*fBd2DstKpi3/fBd2DstKpi0;
        Double_t nBs2DstKpi3 = nbkg3.getVal()*fBs2DstKpi3/fBs2DstKpi0;
        Double_t nBd2Dpipi3  = nbkg4.getVal()*fBd2Dpipi3/fBd2Dpipi0;  
        Double_t nLb2Dppi3   = nbkg5.getVal()*fLb2Dppi3/fLb2Dppi0;      
        Double_t nLb2DpK3    = nbkg6.getVal()*fLb2DpK3/fLb2DpK0;      
        Double_t nBd2DKK3    = nbkg7.getVal()*fBd2DKK3/fBd2DKK0;      
        Double_t nBs2DKK3    = nbkg8.getVal()*fBs2DKK3/fBs2DKK0;      

	cout << "Bd window"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd1*afrac->getValV()          , nBd1*bfrac->getValV()          , nBd1*cfrac->getValV()          , nBd1*dfrac->getValV()          , nBd1*efrac->getValV());
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBs1*afrac->getValV()          , nBs1*bfrac->getValV()          , nBs1*cfrac->getValV()          , nBs1*dfrac->getValV()          , nBs1*efrac->getValV());
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nComb1*abkgfrac->getValV()     , nComb1*bbkgfrac->getValV()     , nComb1*cbkgfrac->getValV()     , nComb1*dbkgfrac->getValV()     , nComb1*ebkgfrac->getValV());
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2DstKpi1*aprfrac->getValV() , nBd2DstKpi1*bprfrac->getValV() , nBd2DstKpi1*cprfrac->getValV() , nBd2DstKpi1*dprfrac->getValV() , nBd2DstKpi1*eprfrac->getValV());
	printf("Bs2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBs2DstKpi1*aprfrac->getValV() , nBs2DstKpi1*bprfrac->getValV() , nBs2DstKpi1*cprfrac->getValV() , nBs2DstKpi1*dprfrac->getValV() , nBs2DstKpi1*eprfrac->getValV());
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2Dpipi1*afrac->getValV()    , nBd2Dpipi1*bfrac->getValV()    , nBd2Dpipi1*cfrac->getValV()    , nBd2Dpipi1*dfrac->getValV()    , nBd2Dpipi1*efrac->getValV());
	printf("Lb2Dppi    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nLb2Dppi1*afrac->getValV()     , nLb2Dppi1*bfrac->getValV()     , nLb2Dppi1*cfrac->getValV()     , nLb2Dppi1*dfrac->getValV()     , nLb2Dppi1*efrac->getValV());
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nLb2DpK1*afrac->getValV()      , nLb2DpK1*bfrac->getValV()      , nLb2DpK1*cfrac->getValV()      , nLb2DpK1*dfrac->getValV()      , nLb2DpK1*efrac->getValV());
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2DKK1*afrac->getValV()      , nBd2DKK1*bfrac->getValV()      , nBd2DKK1*cfrac->getValV()      , nBd2DKK1*dfrac->getValV()      , nBd2DKK1*efrac->getValV());
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBs2DKK1*afrac->getValV()      , nBs2DKK1*bfrac->getValV()      , nBs2DKK1*cfrac->getValV()      , nBs2DKK1*dfrac->getValV()      , nBs2DKK1*efrac->getValV());
	cout << endl;

	cout << "Bs window"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2           , bfrac->getValV()*nBd2           , cfrac->getValV()*nBd2           , dfrac->getValV()*nBd2           , efrac->getVal()*nBd2          );
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2           , bfrac->getValV()*nBs2           , cfrac->getValV()*nBs2           , dfrac->getValV()*nBs2           , efrac->getVal()*nBs2          );
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", abkgfrac->getValV()*nComb2      , bbkgfrac->getValV()*nComb2      , cbkgfrac->getValV()*nComb2      , dbkgfrac->getValV()*nComb2      , ebkgfrac->getVal()*nComb2      );
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBd2DstKpi2  , bprfrac->getValV()*nBd2DstKpi2  , cprfrac->getValV()*nBd2DstKpi2  , dprfrac->getValV()*nBd2DstKpi2  , eprfrac->getVal()*nBd2DstKpi2);
	printf("Bs2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBs2DstKpi2  , bprfrac->getValV()*nBs2DstKpi2  , cprfrac->getValV()*nBs2DstKpi2  , dprfrac->getValV()*nBs2DstKpi2  , eprfrac->getVal()*nBs2DstKpi2);
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2Dpipi2     , bfrac->getValV()*nBd2Dpipi2     , cfrac->getValV()*nBd2Dpipi2     , dfrac->getValV()*nBd2Dpipi2     , efrac->getVal()*nBd2Dpipi2    );
	printf("Lb2Dppi    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2Dppi2      , bfrac->getValV()*nLb2Dppi2      , cfrac->getValV()*nLb2Dppi2      , dfrac->getValV()*nLb2Dppi2      , efrac->getVal()*nLb2Dppi2      );
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2DpK2       , bfrac->getValV()*nLb2DpK2       , cfrac->getValV()*nLb2DpK2       , dfrac->getValV()*nLb2DpK2       , efrac->getVal()*nLb2DpK2      );
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2DKK2       , bfrac->getValV()*nBd2DKK2       , cfrac->getValV()*nBd2DKK2       , dfrac->getValV()*nBd2DKK2       , efrac->getVal()*nBd2DKK2      );
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2DKK2       , bfrac->getValV()*nBs2DKK2       , cfrac->getValV()*nBs2DKK2       , dfrac->getValV()*nBs2DKK2       , efrac->getVal()*nBs2DKK2      );
	cout << endl;

	cout << "sideband"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd3           , bfrac->getValV()*nBd3           , cfrac->getValV()*nBd3           , dfrac->getValV()*nBd3           , efrac->getVal()*nBd3          );
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs3           , bfrac->getValV()*nBs3           , cfrac->getValV()*nBs3           , dfrac->getValV()*nBs3           , efrac->getVal()*nBs3          );
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", abkgfrac->getValV()*nComb3      , bbkgfrac->getValV()*nComb3      , cbkgfrac->getValV()*nComb3      , dbkgfrac->getValV()*nComb3      , ebkgfrac->getVal()*nComb3      );
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBd2DstKpi3  , bprfrac->getValV()*nBd2DstKpi3  , cprfrac->getValV()*nBd2DstKpi3  , dprfrac->getValV()*nBd2DstKpi3  , eprfrac->getVal()*nBd2DstKpi3);
	printf("Bs2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBs2DstKpi3  , bprfrac->getValV()*nBs2DstKpi3  , cprfrac->getValV()*nBs2DstKpi3  , dprfrac->getValV()*nBs2DstKpi3  , eprfrac->getVal()*nBs2DstKpi3);
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2Dpipi3     , bfrac->getValV()*nBd2Dpipi3     , cfrac->getValV()*nBd2Dpipi3     , dfrac->getValV()*nBd2Dpipi3     , efrac->getVal()*nBd2Dpipi3    );
	printf("Lb2Dppi    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2Dppi3      , bfrac->getValV()*nLb2Dppi3      , cfrac->getValV()*nLb2Dppi3      , dfrac->getValV()*nLb2Dppi3      , efrac->getVal()*nLb2Dppi3      );
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2DpK3       , bfrac->getValV()*nLb2DpK3       , cfrac->getValV()*nLb2DpK3       , dfrac->getValV()*nLb2DpK3       , efrac->getVal()*nLb2DpK3      );
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2DKK3       , bfrac->getValV()*nBd2DKK3       , cfrac->getValV()*nBd2DKK3       , dfrac->getValV()*nBd2DKK3       , efrac->getVal()*nBd2DKK3      );
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2DKK3       , bfrac->getValV()*nBs2DKK3       , cfrac->getValV()*nBs2DKK3       , dfrac->getValV()*nBs2DKK3       , efrac->getVal()*nBs2DKK3      );
	cout << endl;

	printf("$B^0$ in window       & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ \\\\\n", 
		afrac->getValV()*nBd1, (afrac->getValV()*nBd1/ansig->getVal())*ansig->getPropagatedError(*r),
		bfrac->getValV()*nBd1, (bfrac->getValV()*nBd1/bnsig->getVal())*bnsig->getPropagatedError(*r),
		cfrac->getValV()*nBd1, (cfrac->getValV()*nBd1/cnsig->getVal())*cnsig->getPropagatedError(*r),
		dfrac->getValV()*nBd1, (dfrac->getValV()*nBd1/dnsig->getVal())*dnsig->getPropagatedError(*r),
		efrac->getValV()*nBd1, (efrac->getValV()*nBd1/ensig->getVal())*ensig->getPropagatedError(*r));
	printf("$B^0$ window purity   & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%%\\\\\n",
		100*afrac->getValV()*nBd1/(afrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2Dppi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+abkgfrac->getValV()*nComb1+aprfrac->getValV()*(nBd2DstKpi1+nBs2DstKpi1)),
		100*bfrac->getValV()*nBd1/(bfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2Dppi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+bbkgfrac->getValV()*nComb1+bprfrac->getValV()*(nBd2DstKpi1+nBs2DstKpi1)),
		100*cfrac->getValV()*nBd1/(cfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2Dppi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+cbkgfrac->getValV()*nComb1+cprfrac->getValV()*(nBd2DstKpi1+nBs2DstKpi1)),
		100*dfrac->getValV()*nBd1/(dfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2Dppi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+dbkgfrac->getValV()*nComb1+dprfrac->getValV()*(nBd2DstKpi1+nBs2DstKpi1)),
		100*efrac->getValV()*nBd1/(efrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2Dppi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+ebkgfrac->getValV()*nComb1+eprfrac->getValV()*(nBd2DstKpi1+nBs2DstKpi1)));
	printf("$B_s^0$ in window     & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ \\\\\n",
		afrac->getValV()*nBs2, (afrac->getValV()*nBs2/ansig2->getVal())*ansig2->getPropagatedError(*r),
		bfrac->getValV()*nBs2, (bfrac->getValV()*nBs2/bnsig2->getVal())*bnsig2->getPropagatedError(*r),
		cfrac->getValV()*nBs2, (cfrac->getValV()*nBs2/cnsig2->getVal())*cnsig2->getPropagatedError(*r),
		dfrac->getValV()*nBs2, (dfrac->getValV()*nBs2/dnsig2->getVal())*dnsig2->getPropagatedError(*r),
		efrac->getValV()*nBs2, (efrac->getValV()*nBs2/ensig2->getVal())*ensig2->getPropagatedError(*r));
	printf("$B_s^0$ window purity & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%%\\\\\n",
		100*afrac->getValV()*nBs2/(afrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2Dppi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+abkgfrac->getValV()*nComb2+aprfrac->getValV()*(nBd2DstKpi2+nBs2DstKpi2)),
		100*bfrac->getValV()*nBs2/(bfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2Dppi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+bbkgfrac->getValV()*nComb2+bprfrac->getValV()*(nBd2DstKpi2+nBs2DstKpi2)),
		100*cfrac->getValV()*nBs2/(cfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2Dppi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+cbkgfrac->getValV()*nComb2+cprfrac->getValV()*(nBd2DstKpi2+nBs2DstKpi2)),
		100*dfrac->getValV()*nBs2/(dfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2Dppi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+dbkgfrac->getValV()*nComb2+dprfrac->getValV()*(nBd2DstKpi2+nBs2DstKpi2)),
		100*efrac->getValV()*nBs2/(efrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2Dppi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+ebkgfrac->getValV()*nComb2+eprfrac->getValV()*(nBd2DstKpi2+nBs2DstKpi2)));
	cout << endl;

	std::cout << afrac->getValV() << " +/- " << afrac->getError() << std::endl;
	std::cout << bfrac->getValV() << " +/- " << bfrac->getError() << std::endl;
	std::cout << cfrac->getValV() << " +/- " << cfrac->getError() << std::endl;
	std::cout << dfrac->getValV() << " +/- " << dfrac->getError() << std::endl;
	std::cout << efrac->getVal()  << " +/- " << efrac->getPropagatedError(*r) << std::endl;
	std::cout << std::endl;
	std::cout << ansig->getVal()  << " +/- " << ansig->getPropagatedError(*r) << std::endl;
	std::cout << bnsig->getVal()  << " +/- " << bnsig->getPropagatedError(*r) << std::endl;
	std::cout << cnsig->getVal()  << " +/- " << cnsig->getPropagatedError(*r) << std::endl;
	std::cout << dnsig->getVal()  << " +/- " << dnsig->getPropagatedError(*r) << std::endl;
	std::cout << ensig->getVal()  << " +/- " << ensig->getPropagatedError(*r) << std::endl;
	std::cout << std::endl;
	std::cout << ansig->getVal()*ansig->getVal()/ansig->getPropagatedError(*r)/ansig->getPropagatedError(*r) << std::endl;
	std::cout << bnsig->getVal()*bnsig->getVal()/bnsig->getPropagatedError(*r)/bnsig->getPropagatedError(*r) << std::endl;
	std::cout << cnsig->getVal()*cnsig->getVal()/cnsig->getPropagatedError(*r)/cnsig->getPropagatedError(*r) << std::endl;
	std::cout << dnsig->getVal()*dnsig->getVal()/dnsig->getPropagatedError(*r)/dnsig->getPropagatedError(*r) << std::endl;
	std::cout << ensig->getVal()*ensig->getVal()/ensig->getPropagatedError(*r)/ensig->getPropagatedError(*r) << std::endl;
	std::cout << std::endl;
	std::cout << nsig->getValV() << " +/- " << nsig->getError() << std::endl;
	std::cout << std::endl;
	std::cout << nsig->getValV()*nsig->getValV()/nsig->getError()/nsig->getError() << std::endl;

	printf("& $%5.3f \\pm %5.3f$ & $%5.3f \\pm %5.3f$ & $%5.3f \\pm %5.3f$ & $%5.3f \\pm %5.3f$ & $%5.3f \\pm %5.3f$ & \\\\\n",afrac->getValV(),afrac->getError(),bfrac->getValV(),bfrac->getError(),cfrac->getValV(),cfrac->getError(),dfrac->getValV(),dfrac->getError(),efrac->getVal(),efrac->getPropagatedError(*r));
	printf("& $%5.0f \\pm %5.0f$ & $%5.0f \\pm %5.0f$ & $%5.0f \\pm %5.0f$ & $%5.0f \\pm %5.0f$ & $%5.0f \\pm %5.0f$ & $%5.0f \\pm %5.0f$ \\\\\n",ansig->getVal(),ansig->getPropagatedError(*r),bnsig->getVal(),bnsig->getPropagatedError(*r),cnsig->getVal(),cnsig->getPropagatedError(*r),dnsig->getVal(),dnsig->getPropagatedError(*r),ensig->getVal(),ensig->getPropagatedError(*r),nsig->getValV(),nsig->getError());
	printf("& %5.0f           & %5.0f           & %5.0f           & %5.0f           & %5.0f           & %5.0f \\\\\n",ansig->getVal()*ansig->getVal()/ansig->getPropagatedError(*r)/ansig->getPropagatedError(*r),bnsig->getVal()*bnsig->getVal()/bnsig->getPropagatedError(*r)/bnsig->getPropagatedError(*r),cnsig->getVal()*cnsig->getVal()/cnsig->getPropagatedError(*r)/cnsig->getPropagatedError(*r),dnsig->getVal()*dnsig->getVal()/dnsig->getPropagatedError(*r)/dnsig->getPropagatedError(*r),ensig->getVal()*ensig->getVal()/ensig->getPropagatedError(*r)/ensig->getPropagatedError(*r),nsig->getValV()*nsig->getValV()/nsig->getError()/nsig->getError());

	TH1D* h1 = new TH1D("h1","",400,0.,200.);
	TH1D* h2 = new TH1D("h2","",400,0.,200.);
	TH1D* h3 = new TH1D("h3","",400,0.,200.);
	TH1D* h4 = new TH1D("h4","",400,0.,200.);
	TH1D* h5 = new TH1D("h5","",400,0.,200.);

	RooArgSet* vars = ansig->getVariables(); 

	for(Int_t i=0; i<10000; ++i) {
		*vars = r->randomizePars();
		h1->Fill(ansig->getVal());
		h2->Fill(bnsig->getVal());
		h3->Fill(cnsig->getVal());
		h4->Fill(dnsig->getVal());
		h5->Fill(ensig->getVal());
	}

	std::cout << h1->GetMean() << " +/- " << h1->GetRMS() << std::endl;
	std::cout << h2->GetMean() << " +/- " << h2->GetRMS() << std::endl;
	std::cout << h3->GetMean() << " +/- " << h3->GetRMS() << std::endl;
	std::cout << h4->GetMean() << " +/- " << h4->GetRMS() << std::endl;
	std::cout << h5->GetMean() << " +/- " << h5->GetRMS() << std::endl;

	TF1 gaus("gaus","gaus(0)",0,200);
	gaus.SetParameter(0,10000);

	gaus.SetParameter(1,h1->GetMean());
	gaus.SetParameter(2,h1->GetRMS());
	h1->Fit(&gaus,"S");
	h1->GetFunction("gaus")->SetLineColor(kRed);
	h1->GetXaxis()->SetTitle("N_{sig}");
	h1->GetYaxis()->SetTitle("N_{toys}");
	h1->Draw(); temp->SaveAs("yieldA.pdf");
	std::cout << h1->GetFunction("gaus")->GetParameter(1) << " +/- " << h1->GetFunction("gaus")->GetParameter(2) << std::endl;

	gaus.SetParameter(1,h2->GetMean());
	gaus.SetParameter(2,h2->GetRMS());
	h2->Fit(&gaus,"S");
	h2->GetFunction("gaus")->SetLineColor(kRed);
	h2->GetXaxis()->SetTitle("N_{sig}");
	h2->GetYaxis()->SetTitle("N_{toys}");
	h2->Draw(); temp->SaveAs("yieldB.pdf");
	std::cout << h2->GetFunction("gaus")->GetParameter(1) << " +/- " << h2->GetFunction("gaus")->GetParameter(2) << std::endl;

	gaus.SetParameter(1,h3->GetMean());
	gaus.SetParameter(2,h3->GetRMS());
	h3->Fit(&gaus,"S");
	h3->GetFunction("gaus")->SetLineColor(kRed);
	h3->GetXaxis()->SetTitle("N_{sig}");
	h3->GetYaxis()->SetTitle("N_{toys}");
	h3->Draw(); temp->SaveAs("yieldC.pdf");
	std::cout << h3->GetFunction("gaus")->GetParameter(1) << " +/- " << h3->GetFunction("gaus")->GetParameter(2) << std::endl;

	gaus.SetParameter(1,h4->GetMean());
	gaus.SetParameter(2,h4->GetRMS());
	h4->Fit(&gaus,"S");
	h4->GetFunction("gaus")->SetLineColor(kRed);
	h4->GetXaxis()->SetTitle("N_{sig}");
	h4->GetYaxis()->SetTitle("N_{toys}");
	h4->Draw(); temp->SaveAs("yieldD.pdf");
	std::cout << h4->GetFunction("gaus")->GetParameter(1) << " +/- " << h4->GetFunction("gaus")->GetParameter(2) << std::endl;

	gaus.SetParameter(1,h5->GetMean());
	gaus.SetParameter(2,h5->GetRMS());
	h5->Fit(&gaus,"S");
	h5->GetFunction("gaus")->SetLineColor(kRed);
	h5->GetXaxis()->SetTitle("N_{sig}");
	h5->GetYaxis()->SetTitle("N_{toys}");
	h5->Draw(); temp->SaveAs("yieldE.pdf");
	std::cout << h5->GetFunction("gaus")->GetParameter(1) << " +/- " << h5->GetFunction("gaus")->GetParameter(2) << std::endl;

	printf("& $% 7.1f\\pm% -6.1f$ \n", sigmean->getVal(),  sigmean->getError());
	printf("& $% 7.1f\\pm% -6.1f$ \n", sigsigma->getVal(), sigsigma->getError());
	printf("& $% 7.3f\\pm% -6.3f$ \n", frac->getVal(),     frac->getError());
	printf("& $% 7.2f\\pm% -6.2f$ \n", ratio->getVal(),    ratio->getError());
	printf("& $% 7.2f\\pm% -6.2f$ \n", 1e3*p0->getVal(),   1e3*p0->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nsig->getVal(),     nsig->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nsig2->getVal(),    nsig2->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg->getVal(),     nbkg->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg2->getVal(),    nbkg2->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg3->getVal(),    nbkg3->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg4->getVal(),    nbkg4->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg5->getVal(),    nbkg5->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg6->getVal(),    nbkg6->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg7->getVal(),    nbkg7->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg8->getVal(),    nbkg8->getError());
	printf("\n");
	printf("& $% 7.3f\\pm% -6.3f$ \n", afrac->getVal(),    afrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", bfrac->getVal(),    bfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", cfrac->getVal(),    cfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", dfrac->getVal(),    dfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", efrac->getVal(),    efrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", aprfrac->getVal(),  aprfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", bprfrac->getVal(),  bprfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", cprfrac->getVal(),  cprfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", dprfrac->getVal(),  dprfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", eprfrac->getVal(),  eprfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", abkgfrac->getVal(), abkgfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", bbkgfrac->getVal(), bbkgfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", cbkgfrac->getVal(), cbkgfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", dbkgfrac->getVal(), dbkgfrac->getPropagatedError(*r));
	printf("& $% 7.3f\\pm% -6.3f$ \n", ebkgfrac->getVal(), ebkgfrac->getPropagatedError(*r));

}
