{
	gSystem.Load("libRooFit");
	gROOT.SetStyle("Plain");
//	gStyle.SetOptStat(1111);

	gROOT->ProcessLine(".L RooMyPdf.cxx+");

	gStyle->SetOptStat(0000);
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();

	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2Kpi_Scaled/B2D0Kpi_D02Kpi_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root");
	//TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_DKpi_Scaled/TSS_baseline_DKpi_selBd_allVetoes_Dst25_DD55_DK52_PID3p_NNfix.root");
	//TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_DKpi_Scaled/TSS_baseline_DKpi_selBd_allVetoes_Dst25_DD55_PID3p_NNfix.root");
	TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

	//B_D0_B_CM_M 
	RooRealVar B_D0_B_CM_M("B_D0_B_CM_M","",5100,5900);
	RooRealVar NN("NN","",-1.0,1.0);
	RooDataSet * data1 = new RooDataSet("data1", "", DecayTree, RooArgSet(B_D0_B_CM_M,NN));

	RooDataSet * adata = data1->reduce("NN>-0.80 && NN<= 0.00");
	RooDataSet * bdata = data1->reduce("NN> 0.00 && NN<= 0.50");
	RooDataSet * cdata = data1->reduce("NN> 0.50 && NN<= 0.72");
	RooDataSet * ddata = data1->reduce("NN> 0.72 && NN<= 0.81");
	RooDataSet * edata = data1->reduce("NN> 0.81 && NN<= 1.00");

	//Bd -> D* K pi
//	TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/Bd2DKpi/Bd2Dst0Kpi/Bd2DKpi_Bd2Dst0Kpi_P8_selBd_Dst25_NNf1_PID_misIDIMs.root");
	//TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/Bd2Dst0Kpi/D2Kpi_Bd2Dst0Kpi_selBd_Dsignal_vetoes_PID_NND2Kpi_addIMs_addMisIDIMs.root");
	//TTree  * treeK1  = dynamic_cast<TTree*>(fileK1->Get("DecayTree"));
	//RooDataSet* dataK1 = new RooDataSet("dataK1","",treeK1,RooArgSet(B_D0_B_CM_M,NN));
	//RooDataSet* dataK1a = dataK1->reduce("NN>-0.8");
	//RooRealVar shift("shift","",-5.,-20.,-0.1);
	//RooMyPdf bd2dstkpi("bd2dstkpi","",B_D0_B_CM_M,shift,*dataK1a,RooMyPdf::MirrorBoth,2);
	TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bd2Dst0Kpi_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK1  = dynamic_cast<TH2D*>(fileK1->Get("Bd2Dst0Kpi"));
	TH1D  * histK1a = histK1->ProjectionX("Bd2D0stkpi_A",11,100);
	histK1a->Smooth(1);
	RooDataHist * dataK1 = new RooDataHist("dataK1", "", RooArgSet(B_D0_B_CM_M), histK1a);
	RooHistPdf bd2dstkpi("bd2dstkpi","", RooArgSet(B_D0_B_CM_M), *dataK1, 2);

	//// B0 -> D(*) pi pi
//	TFile * fileK2 = TFile::Open("../../../reweightBackgrounds/BdBackgroundsNewDst25_NNf1/Dpipi/B2DKpi_Bd2D_st_0pipi_hist_Dst25_smooth_new.root");
	TFile * fileK2 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bd2D0pipi_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK2  = dynamic_cast<TH2D*>(fileK2->Get("Bd2D0pipi"));
	TH1D  * histK2a = histK2->ProjectionX("Bd2D0pipi_A",11,100);
	histK2a->Smooth(1);
	RooDataHist * dataK2 = new RooDataHist("dataK2", "", RooArgSet(B_D0_B_CM_M), histK2a);
	RooHistPdf bd2dpipi("bd2dpipi","", RooArgSet(B_D0_B_CM_M), *dataK2, 2);
	
	//// Lambda b -> D p K
//	TFile * fileK3 = TFile::Open("../../../reweightBackgrounds/BdBackgroundsNewDst25_NNf1/DpK/B2DKpi_Lb2D_st_0pK_hist_smooth_kill.root");
	TFile * fileK3 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Lb2D0pK_smooth3_bins50_SDPweight_hist.root");
	TH2D  * histK3  = dynamic_cast<TH2D*>(fileK3->Get("Lb2D0pK"));
	TH1D  * histK3a = histK3->ProjectionX("Lb2D0pK_A",11,100);
	histK3a->Smooth(3);
	RooDataHist * dataK3 = new RooDataHist("dataK3", "", RooArgSet(B_D0_B_CM_M), histK3a);
	RooHistPdf lb2dpk("lb2dpk","", RooArgSet(B_D0_B_CM_M), *dataK3, 2);

	////Bd -> D K K
//	TFile * fileK4 = TFile::Open("../../../reweightBackgrounds/BdBackgroundsNewDst25_NNf1/DKK/B2DKpi_Bd2D_st_0KK_hist_Dst25_newSDP_smooth_kill.root");
	TFile * fileK4 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bd2D0KK_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK4  = dynamic_cast<TH2D*>(fileK4->Get("Bd2D0KK"));
	TH1D  * histK4a = histK4->ProjectionX("Bd2D0KK_A",11,100);
	histK4a->Smooth(1);
	RooDataHist * dataK4 = new RooDataHist("dataK4", "", RooArgSet(B_D0_B_CM_M), histK4a);
	RooHistPdf bd2dkk("bd2dkk","", RooArgSet(B_D0_B_CM_M), *dataK4, 2);

	//Bs -> D K K
//	TFile * fileK5 = TFile::Open("../../../reweightBackgrounds/BdBackgroundsNewDst25_NNf1/Bs2DKK/B2DKpi_Bs2D_st_0KK_hist_Dst25_newSDP_smooth_kill.root");
	TFile * fileK5 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bs2D0KK_smooth1_bins50_SDPweight_hist.root");
	TH2D  * histK5  = dynamic_cast<TH2D*>(fileK5->Get("Bs2D0KK"));
	TH1D  * histK5a = histK5->ProjectionX("Bs2D0KK_A",11,100);
	histK5a->Smooth(1);
	RooDataHist * dataK5 = new RooDataHist("dataK5", "", RooArgSet(B_D0_B_CM_M), histK5a);
	RooHistPdf bs2dkk("bs2dkk","", RooArgSet(B_D0_B_CM_M), *dataK5, 2);

	//Bu -> D* K
//	TFile * fileK6 = TFile::Open("../../../reweightBackgrounds/DpipiBackgrounds/Bu2D0pi/B2Dpipi_Bu2Dstpi_selDpipi_hist_kill.root");
	TFile * fileK6 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/hists/D2Kpi_Bu2Dst0K_smooth10_bins40_hist.root");
	TH2D  * histK6  = dynamic_cast<TH2D*>(fileK6->Get("Bu2Dst0K"));
	TH1D  * histK6a = histK6->ProjectionX("Bu2Dst0K_A",11,100);//11,50);
	TH1D  * histK6b = histK6->ProjectionX("Bu2Dst0K_B",11,100);//51,86);
	TH1D  * histK6c = histK6->ProjectionX("Bu2Dst0K_C",11,100);//87,100);
	histK6a->Smooth(10);
	histK6b->Smooth(10);
	histK6c->Smooth(10);
	RooDataHist * dataK6a = new RooDataHist("dataK6a", "", RooArgSet(B_D0_B_CM_M), histK6a);
	RooHistPdf bu2dstkA("bu2dstkA","", RooArgSet(B_D0_B_CM_M), *dataK6a, 2);
	RooDataHist * dataK6b = new RooDataHist("dataK6b", "", RooArgSet(B_D0_B_CM_M), histK6b);
	RooHistPdf bu2dstkB("bu2dstkB","", RooArgSet(B_D0_B_CM_M), *dataK6b, 2);
	RooDataHist * dataK6c = new RooDataHist("dataK6c", "", RooArgSet(B_D0_B_CM_M), histK6c);
	RooHistPdf bu2dstkC("bu2dstkC","", RooArgSet(B_D0_B_CM_M), *dataK6c, 2);

	//Bs -> D* K pi
//	TFile * fileK7 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/Bd2DKpi/Bd2Dst0Kpi/Bd2DKpi_Bs2Dst0Kpi_P8_selBd_Dst25_NNf1_PID_misIDIMs.root");
//	TFile * fileK1 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2Kpi/Bd2Dst0Kpi/D2Kpi_Bd2Dst0Kpi_selBd_Dsignal_vetoes_PID_NND2Kpi_addIMs_addMisIDIMs.root");
//	TTree  * treeK7  = dynamic_cast<TTree*>(fileK7->Get("DecayTree"));
//	RooDataSet* dataK7 = new RooDataSet("dataK7","",treeK7,RooArgSet(B_D0_B_CM_M));
//	RooMyPdf bs2dstkpi("bs2dstkpi","",B_D0_B_CM_M,shift,*dataK7,RooMyPdf::MirrorBoth,2);

	// B Gaussian 
	// start, range to from. plus names and titles.
	RooRealVar sigmean("M_{B}","B mass",5281.0,5270.0,5290.0,"MeV/c^{2}");
	RooRealVar sigsigma("#sigma_{B}","B sigma",12.0,0.0,50.0,"MeV/c^{2}");
	RooRealVar ratio("ratio","Ratio of widths",1.756,1.0,3.0);
	RooProduct sigsigma2("#sigma_{B}2","B sigma2",RooArgSet(sigsigma,ratio));
	RooRealVar frac("frac","fraction of events in each gaussian",0.796,0.0,1.0);

	RooRealVar a1("a1","a1",  1.85179e+00);
	RooRealVar n1("n1","n1",  1.30452e+00);
	RooRealVar a2("a2","a2", -2.17629e+00);
	RooRealVar n2("n2","n2",  3.45811e+00);
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

	//RooRealVar combfracA("combfracA","",0.5,0.0,1.0);
	//RooRealVar combfracB("combfracB","",0.5,0.0,1.0);
	//RooRealVar combfracC("combfracC","",0.5,0.0,1.0);
	//RooRealVar combfracD("combfracD","",0.5,0.0,1.0);
	//RooRealVar combfracE("combfracE","",0.5,0.0,1.0);

	//RooAddPdf comb_bkgA("comb_bkgA","",RooArgList(comb_bkg1,bu2dstkA),RooArgList(combfracA));
	//RooAddPdf comb_bkgB("comb_bkgB","",RooArgList(comb_bkg1,bu2dstkB),RooArgList(combfracB));
	//RooAddPdf comb_bkgC("comb_bkgC","",RooArgList(comb_bkg1,bu2dstkB),RooArgList(combfracC));
	//RooAddPdf comb_bkgD("comb_bkgD","",RooArgList(comb_bkg1,bu2dstkC),RooArgList(combfracD));
	//RooAddPdf comb_bkgE("comb_bkgE","",RooArgList(comb_bkg1,bu2dstkC),RooArgList(combfracE));

	// Number of signal & background events
	RooRealVar nsig("nsig","",4500,-1000,10000);
	RooRealVar nsig2("nsig2","",200,-1000,1000);
	RooRealVar nbkg("nbkg","",7500,-1000,50000);
	RooRealVar nbkg2("nbkg2","",7500,-1000,50000);
	RooRealVar nbkg3("nbkg3","",4500,-1000,10000);
	RooRealVar nbkg4("nbkg4","",4500*0.192,-1000,5000);
	RooRealVar nbkg5("nbkg5","",4500*0.088,-1000,5000);
	RooRealVar nbkg6("nbkg6","",4500*0.089,-1000,5000);
	RooRealVar nbkg7("nbkg7","",4500*0.036,-1000,5000);
//	RooRealVar nbkg8("nbkg8","",1200,-1000,5000);

	//// Constrain Dpipi background
	RooFormulaVar bd2dpipiratio("bd2dpipiratio","@0/@1",RooArgList(nbkg4,nsig));
	RooGaussian constrain_bd2dpipiratio("cbd2dpipi","cbd2dpipi",bd2dpipiratio,RooFit::RooConst(0.192),RooFit::RooConst(0.019));

	//// Constrain DpK background
	RooFormulaVar lb2dpkratio("lb2dpkratio","@0/@1",RooArgList(nbkg5,nsig));
	RooGaussian constrain_lb2dpkratio("clb2dpk","clb2dpk",lb2dpkratio,RooFit::RooConst(0.088),RooFit::RooConst(0.022));
	
	//// Constrain DKK background
	RooFormulaVar bd2dkkratio("bd2dkkratio","@0/@1",RooArgList(nbkg6,nsig));
	RooGaussian constrain_bd2dkkratio("cbd2dkk","cbd2dkk",bd2dkkratio,RooFit::RooConst(0.089),RooFit::RooConst(0.022));
	
	//// Constrain BsDKK background
	RooFormulaVar bs2dkkratio("bs2dkkratio","@0/@1",RooArgList(nbkg7,nsig));
	RooGaussian constrain_bs2dkkratio("cbs2dkk","cbs2dkk",bs2dkkratio,RooFit::RooConst(0.036),RooFit::RooConst(0.016));

//	RooRealVar nbkg4("nbkg4","",4500*0.176,-1000,5000);
//	RooRealVar nbkg5("nbkg5","",4500*0.037,-1000,5000);
//	RooRealVar nbkg6("nbkg6","",4500*0.045,-1000,5000);
//	RooRealVar nbkg7("nbkg7","",4500*0.017,-1000,5000);
//	RooRealVar nbkg8("nbkg8","",1200,-1000,5000);
//
//	//// Constrain Dpipi background
//	RooFormulaVar bd2dpipiratio("bd2dpipiratio","@0/@1",RooArgList(nbkg4,nsig));
//	RooGaussian constrain_bd2dpipiratio("cbd2dpipi","cbd2dpipi",bd2dpipiratio,RooFit::RooConst(0.176),RooFit::RooConst(0.018));
//
//	//// Constrain DpK background
//	RooFormulaVar lb2dpkratio("lb2dpkratio","@0/@1",RooArgList(nbkg5,nsig));
//	RooGaussian constrain_lb2dpkratio("clb2dpk","clb2dpk",lb2dpkratio,RooFit::RooConst(0.037),RooFit::RooConst(0.010));
//	
//	//// Constrain DKK background
//	RooFormulaVar bd2dkkratio("bd2dkkratio","@0/@1",RooArgList(nbkg6,nsig));
//	RooGaussian constrain_bd2dkkratio("cbd2dkk","cbd2dkk",bd2dkkratio,RooFit::RooConst(0.045),RooFit::RooConst(0.011));
//	
//	//// Constrain BsDKK background
//	RooFormulaVar bs2dkkratio("bs2dkkratio","@0/@1",RooArgList(nbkg7,nsig));
//	RooGaussian constrain_bs2dkkratio("cbs2dkk","cbs2dkk",bs2dkkratio,RooFit::RooConst(0.017),RooFit::RooConst(0.007));
//
	//RooRealVar nbkg4("nbkg4","",4500*0.160,-1000,5000);
	//RooRealVar nbkg5("nbkg5","",4500*0.085,-1000,5000);
	//RooRealVar nbkg6("nbkg6","",4500*0.067,-1000,5000);
	//RooRealVar nbkg7("nbkg7","",4500*0.020,-1000,5000);
	//RooRealVar nbkg8("nbkg8","",1200,-1000,5000);

	////// Constrain Dpipi background
	//RooFormulaVar bd2dpipiratio("bd2dpipiratio","@0/@1",RooArgList(nbkg4,nsig));
	//RooGaussian constrain_bd2dpipiratio("cbd2dpipi","cbd2dpipi",bd2dpipiratio,RooFit::RooConst(0.160),RooFit::RooConst(0.016));

	////// Constrain DpK background
	//RooFormulaVar lb2dpkratio("lb2dpkratio","@0/@1",RooArgList(nbkg5,nsig));
	//RooGaussian constrain_lb2dpkratio("clb2dpk","clb2dpk",lb2dpkratio,RooFit::RooConst(0.085),RooFit::RooConst(0.022));
	//
	////// Constrain DKK background
	//RooFormulaVar bd2dkkratio("bd2dkkratio","@0/@1",RooArgList(nbkg6,nsig));
	//RooGaussian constrain_bd2dkkratio("cbd2dkk","cbd2dkk",bd2dkkratio,RooFit::RooConst(0.067),RooFit::RooConst(0.016));
	//
	////// Constrain BsDKK background
	//RooFormulaVar bs2dkkratio("bs2dkkratio","@0/@1",RooArgList(nbkg7,nsig));
	//RooGaussian constrain_bs2dkkratio("cbs2dkk","cbs2dkk",bs2dkkratio,RooFit::RooConst(0.020),RooFit::RooConst(0.009));

	RooRealVar afrac("afrac","",0.2,0.0,1.0);
	RooRealVar bfrac("bfrac","",0.2,0.0,1.0);
	RooRealVar cfrac("cfrac","",0.2,0.0,1.0);
	RooRealVar dfrac("dfrac","",0.2,0.0,1.0);
	RooFormulaVar efrac("efrac","","1-afrac-bfrac-cfrac-dfrac",RooArgList(afrac,bfrac,cfrac,dfrac));
	//RooFormulaVar efrac("efrac","","1-dfrac",RooArgList(dfrac));

	RooRealVar aprfrac("aprfrac","",0.2,0.0,1.0);
	RooRealVar bprfrac("bprfrac","",0.2,0.0,1.0);
	RooRealVar cprfrac("cprfrac","",0.2,0.0,1.0);
	RooRealVar dprfrac("dprfrac","",0.2,0.0,1.0);
	RooFormulaVar eprfrac("eprfrac","","1-aprfrac-bprfrac-cprfrac-dprfrac",RooArgList(aprfrac,bprfrac,cprfrac,dprfrac));
	//RooFormulaVar eprfrac("eprfrac","","1-dprfrac",RooArgList(dprfrac));

	RooRealVar abkgfrac("abkgfrac","",0.70,0.0,1.0);
	RooRealVar bbkgfrac("bbkgfrac","",0.20,0.0,1.0);
	RooRealVar cbkgfrac("cbkgfrac","",0.08,0.0,1.0);
	RooRealVar dbkgfrac("dbkgfrac","",0.01,0.0,1.0);
	RooFormulaVar ebkgfrac("ebkgfrac","","1-abkgfrac-bbkgfrac-cbkgfrac-dbkgfrac",RooArgList(abkgfrac,bbkgfrac,cbkgfrac,dbkgfrac));
	//RooFormulaVar ebkgfrac("ebkgfrac","","1-dbkgfrac",RooArgList(dbkgfrac));

	RooRealVar adstkfrac("adstkfrac","",0.70,0.0,1.0);
	RooRealVar bdstkfrac("bdstkfrac","",0.20,0.0,1.0);
	RooRealVar cdstkfrac("cdstkfrac","",0.08,0.0,1.0);
	RooRealVar ddstkfrac("ddstkfrac","",0.01,0.0,1.0);
	RooFormulaVar edstkfrac("edstkfrac","","1-adstkfrac-bdstkfrac-cdstkfrac-ddstkfrac",RooArgList(adstkfrac,bdstkfrac,cdstkfrac,ddstkfrac));
	//RooFormulaVar edstkfrac("edstkfrac","","1-ddstkfrac",RooArgList(ddstkfrac));

	RooFormulaVar anbkg("anbkg","abkgfrac*nbkg",RooArgList(abkgfrac,nbkg));
	RooFormulaVar bnbkg("bnbkg","bbkgfrac*nbkg",RooArgList(bbkgfrac,nbkg));
	RooFormulaVar cnbkg("cnbkg","cbkgfrac*nbkg",RooArgList(cbkgfrac,nbkg));
	RooFormulaVar dnbkg("dnbkg","dbkgfrac*nbkg",RooArgList(dbkgfrac,nbkg));
	RooFormulaVar enbkg("enbkg","ebkgfrac*nbkg",RooArgList(ebkgfrac,nbkg));

	RooFormulaVar anbkg2("anbkg2","abkgfrac*nbkg2",RooArgList(abkgfrac,nbkg2));
	RooFormulaVar bnbkg2("bnbkg2","bbkgfrac*nbkg2",RooArgList(bbkgfrac,nbkg2));
	RooFormulaVar cnbkg2("cnbkg2","cbkgfrac*nbkg2",RooArgList(cbkgfrac,nbkg2));
	RooFormulaVar dnbkg2("dnbkg2","dbkgfrac*nbkg2",RooArgList(dbkgfrac,nbkg2));
	RooFormulaVar enbkg2("enbkg2","ebkgfrac*nbkg2",RooArgList(ebkgfrac,nbkg2));

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

	RooAddPdf afull_PDF("afull_PDF","",RooArgList(BSig,BsSig,comb_bkg,bu2dstkA,bd2dstkpi,bd2dpipi,lb2dpk,bd2dkk,bs2dkk),//,bs2dstkpi),
			                   RooArgList(ansig,ansig2,anbkg,anbkg2,anbkg3,anbkg4,anbkg5,anbkg6,anbkg7));
	RooAddPdf bfull_PDF("bfull_PDF","",RooArgList(BSig,BsSig,comb_bkg,bu2dstkB,bd2dstkpi,bd2dpipi,lb2dpk,bd2dkk,bs2dkk),//,bs2dstkpi),
			                   RooArgList(bnsig,bnsig2,bnbkg,bnbkg2,bnbkg3,bnbkg4,bnbkg5,bnbkg6,bnbkg7));
	RooAddPdf cfull_PDF("cfull_PDF","",RooArgList(BSig,BsSig,comb_bkg,bu2dstkB,bd2dstkpi,bd2dpipi,lb2dpk,bd2dkk,bs2dkk),//,bs2dstkpi),
			                   RooArgList(cnsig,cnsig2,cnbkg,cnbkg2,cnbkg3,cnbkg4,cnbkg5,cnbkg6,cnbkg7));
	RooAddPdf dfull_PDF("dfull_PDF","",RooArgList(BSig,BsSig,comb_bkg,bu2dstkC,bd2dstkpi,bd2dpipi,lb2dpk,bd2dkk,bs2dkk),//,bs2dstkpi),
			                   RooArgList(dnsig,dnsig2,dnbkg,dnbkg2,dnbkg3,dnbkg4,dnbkg5,dnbkg6,dnbkg7));
	RooAddPdf efull_PDF("efull_PDF","",RooArgList(BSig,BsSig,comb_bkg,bu2dstkC,bd2dstkpi,bd2dpipi,lb2dpk,bd2dkk,bs2dkk),//,bs2dstkpi),
			                   RooArgList(ensig,ensig2,enbkg,enbkg2,enbkg3,enbkg4,enbkg5,enbkg6,enbkg7));

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

	RooProdPdf fullpdf("fullpdf","",RooArgList(Simpdf,constrain_ratio,constrain_frac,constrain_bd2dpipiratio,constrain_lb2dpkratio,constrain_bd2dkkratio,constrain_bs2dkkratio));

	//# Do the fit on REFITTED Mass
	RooFitResult * r = fullpdf->fitTo(combo,RooFit::Extended(),RooFit::Save());
	//Simpdf->fitTo(combo,RooFit::Extended());

	TCanvas * can = new TCanvas("can","",2000,1500);
	can->Divide(3,2);
	can.cd(1);
	aPlot = B_D0_B_CM_M->frame(50);
	aPlot->SetTitle("");
	aPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//aPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	aPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(aPlot,RooFit::Cut("cat==cat::cata"));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 A:\t" << aPlot.chiSquare(0) << endl;
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkA)),         RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(aPlot, RooFit::Slice(cat,"cata"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	aPlot->Draw();
	can.cd(2);
	bPlot = B_D0_B_CM_M->frame(50);
	bPlot->SetTitle("");
	bPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//bPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	bPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(bPlot,RooFit::Cut("cat==cat::catb"));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 B:\t" << bPlot.chiSquare(0) << endl;
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkB)),         RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(bPlot, RooFit::Slice(cat,"catb"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	bPlot->Draw();
	can.cd(3);
	cPlot = B_D0_B_CM_M->frame(50);
	cPlot->SetTitle("");
	cPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//cPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	cPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(cPlot,RooFit::Cut("cat==cat::catc"));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 C:\t" << cPlot.chiSquare(0) << endl;
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkB)),         RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(cPlot, RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	cPlot->Draw();
	can.cd(4);
	dPlot = B_D0_B_CM_M->frame(50);
	dPlot->SetTitle("");
	dPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//dPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	dPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(dPlot,RooFit::Cut("cat==cat::catd"));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 D:\t" << dPlot.chiSquare(0) << endl;
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkC)),         RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(dPlot, RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	dPlot->Draw();
	can.cd(5);
	ePlot = B_D0_B_CM_M->frame(50);
	ePlot->SetTitle("");
	ePlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//ePlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	ePlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(ePlot,RooFit::Cut("cat==cat::cate"));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo));
	cout << "Fit chi2 E:\t" << ePlot.chiSquare(0) << endl;
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkC)),         RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	ePlot->Draw();
	can.cd(6);
	fPlot = B_D0_B_CM_M->frame(50);
	fPlot->SetTitle("");
	fPlot->GetYaxis()->SetTitle("Candidates / (16 MeV/c^{2})");
	//fPlot->GetXaxis()->SetTitle("m(DK#pi) (MeV/c^{2})");
	fPlot->GetXaxis()->SetTitle("#it{m}(#it{D}#it{K}^{#font[122]{+}}#pi^{#font[122]{-}}) [MeV/#it{c}^{2}]");
	combo->plotOn(fPlot);
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BSig)),             RooFit::LineStyle(kDashed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(BsSig)),            RooFit::LineStyle(kDashed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(comb_bkg)),         RooFit::LineStyle(10),RooFit::LineColor(kMagenta));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bu2dstkA,bu2dstkB,bu2dstkC)), RooFit::LineStyle(8), RooFit::LineColor(kTeal+2));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dstkpi)),        RooFit::LineStyle(9), RooFit::LineColor(kRed));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dpipi)),         RooFit::LineStyle(5), RooFit::LineColor(kGreen+2));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(lb2dpk)),           RooFit::LineStyle(7), RooFit::LineColor(kGray+2));
	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components(RooArgSet(bd2dkk,bs2dkk)),    RooFit::LineStyle(6), RooFit::LineColor(kOrange-3));
//	Simpdf->plotOn(fPlot, RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//combo->plotOn(fPlot,RooFit::Cut("cat==cat::catc||cat==cat::catd||cat==cat::cate"));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_all"),      RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_BSig"),     RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("BSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_BsSig"),    RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("BsSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_comb_bkg"), RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("comb_bkg"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_bu2dstk"),  RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bu2dstkB"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal+2));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_bd2dstkpi"),RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_bd2dpipi"), RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dpipi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+2));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_lb2dpk"),   RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("lb2dpk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray+2));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_bd2dkk"),   RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(fPlot, RooFit::Name("curve_bs2dkk"),   RooFit::Invisible(),RooFit::Slice(cat,"catc"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));

	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_all"),      RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_BSig"),     RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("BSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_BsSig"),    RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("BsSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_comb_bkg"), RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("comb_bkg"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bu2dstk"),  RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bu2dstkC"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dstkpi"),RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dpipi"), RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dpipi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_lb2dpk"),   RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("lb2dpk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dkk"),   RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bs2dkk"),   RooFit::Invisible(),RooFit::Slice(cat,"catd"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));

	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_all"),      RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_BSig"),     RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("BSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_BsSig"),    RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("BsSig"),RooFit::LineStyle(kDashed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_comb_bkg"), RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("comb_bkg"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bu2dstk"),  RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bu2dstkC"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dstkpi"),RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dpipi"), RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dpipi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_lb2dpk"),   RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("lb2dpk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray+2));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bd2dkk"),   RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bd2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(fPlot, RooFit::AddTo("curve_bs2dkk"),   RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dkk"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3));
	//Simpdf->plotOn(ePlot, RooFit::Slice(cat,"cate"),RooFit::ProjWData(cat,combo), RooFit::Components("bs2dstkpi"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
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
	aPlot3->addPlotable(Plot2->pullHist(),"PE1L");
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
	bPlot3->addPlotable(Plot2->pullHist(),"PE1L");
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
	cPlot3->addPlotable(Plot2->pullHist(),"PE1L");
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
	dPlot3->addPlotable(Plot2->pullHist(),"PE1L");
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
	ePlot3->addPlotable(Plot2->pullHist(),"PE1L");
	ePlot3->addObject( lowLine );
	ePlot3->addObject( midLine );
	ePlot3->addObject( uppLine );
	ePlot3->Draw();

	can->SaveAs("D2Kpifit_fit.C");
	can2->SaveAs("D2Kpifit_pull.C");

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
	bPlot->Draw(); printLHCb("R","Other","LHCb (b)"); temp->SaveAs("figs/D2Kpifit_2.pdf"); temp->SaveAs("figs/D2Kpifit_2.C");
	temp->SaveAs("figs/fig6b.pdf");
	temp->SaveAs("figs/fig6b.png");
	temp->SaveAs("figs/fig6b.eps");
	temp->SaveAs("figs/fig6b.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_2_log.pdf"); temp->SetLogy(0);
	cPlot->Draw(); printLHCb("R","Other","LHCb (c)"); temp->SaveAs("figs/D2Kpifit_3.pdf"); temp->SaveAs("figs/D2Kpifit_3.C");
	temp->SaveAs("figs/fig6c.pdf");
	temp->SaveAs("figs/fig6c.png");
	temp->SaveAs("figs/fig6c.eps");
	temp->SaveAs("figs/fig6c.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_3_log.pdf"); temp->SetLogy(0);
	dPlot->Draw(); printLHCb("R","Other","LHCb (d)"); temp->SaveAs("figs/D2Kpifit_4.pdf"); temp->SaveAs("figs/D2Kpifit_4.C");
	temp->SaveAs("figs/fig6d.pdf");
	temp->SaveAs("figs/fig6d.png");
	temp->SaveAs("figs/fig6d.eps");
	temp->SaveAs("figs/fig6d.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_4_log.pdf"); temp->SetLogy(0);
	ePlot->Draw(); printLHCb("R","Other","LHCb (e)"); temp->SaveAs("figs/D2Kpifit_5.pdf"); temp->SaveAs("figs/D2Kpifit_5.C");
	temp->SaveAs("figs/fig6e.pdf");
	temp->SaveAs("figs/fig6e.png");
	temp->SaveAs("figs/fig6e.eps");
	temp->SaveAs("figs/fig6e.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_5_log.pdf"); temp->SetLogy(0);
	fPlot->Draw(); temp->SaveAs("figs/D2Kpifit_1-5.pdf"); temp->SaveAs("figs/D2Kpifit_1-5.C");
	temp->SetLogy(1); temp->SaveAs("figs/D2Kpifit_1-5_log.pdf"); temp->SetLogy(0);
	aPlot3->Draw(); temp->SaveAs("figs/D2Kpifit_1_pull.pdf"); temp->SaveAs("figs/D2Kpifit_1_pull.C");
	bPlot3->Draw(); temp->SaveAs("figs/D2Kpifit_2_pull.pdf"); temp->SaveAs("figs/D2Kpifit_2_pull.C");
	cPlot3->Draw(); temp->SaveAs("figs/D2Kpifit_3_pull.pdf"); temp->SaveAs("figs/D2Kpifit_3_pull.C");
	dPlot3->Draw(); temp->SaveAs("figs/D2Kpifit_4_pull.pdf"); temp->SaveAs("figs/D2Kpifit_4_pull.C");
	ePlot3->Draw(); temp->SaveAs("figs/D2Kpifit_5_pull.pdf"); temp->SaveAs("figs/D2Kpifit_5_pull.C");
	aPlot4->Draw(); printLHCb("R","Other","LHCb (a)"); temp->SaveAs("figs/D2Kpifit_1_window.pdf"); temp->SaveAs("figs/D2Kpifit_1_window.C");
	temp->SaveAs("figs/fig6a.pdf");
	temp->SaveAs("figs/fig6a.png");
	temp->SaveAs("figs/fig6a.eps");
	temp->SaveAs("figs/fig6a.C");
	bPlot4->Draw(); temp->SaveAs("figs/D2Kpifit_2_window.pdf"); temp->SaveAs("figs/D2Kpifit_2_window.C");
	cPlot4->Draw(); temp->SaveAs("figs/D2Kpifit_3_window.pdf"); temp->SaveAs("figs/D2Kpifit_3_window.C");
	dPlot4->Draw(); temp->SaveAs("figs/D2Kpifit_4_window.pdf"); temp->SaveAs("figs/D2Kpifit_4_window.C");
	ePlot4->Draw(); temp->SaveAs("figs/D2Kpifit_5_window.pdf"); temp->SaveAs("figs/D2Kpifit_5_window.C");

	double mBdm = sigmean.getVal() - 2.5*(sigsigma.getVal());
	double mBdp = sigmean.getVal() + 2.5*(sigsigma.getVal());
	double mBsm = sigmeanBs.getVal() - 2.5*(sigsigma.getVal());
	double mBsp = sigmeanBs.getVal() + 2.5*(sigsigma.getVal());

	B_D0_B_CM_M.setRange("Bd",mBdm,mBdp);
	B_D0_B_CM_M.setRange("Bs",mBsm,mBsp);
	B_D0_B_CM_M.setRange("full",5100,5900);
	B_D0_B_CM_M.setRange("sideband",5400,5900);

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

	double fBu2DstKA1 = bu2dstkA.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBu2DstKA2 = bu2dstkA.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBu2DstKA3 = bu2dstkA.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBu2DstKA0 = bu2dstkA.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBu2DstKB1 = bu2dstkB.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBu2DstKB2 = bu2dstkB.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBu2DstKB3 = bu2dstkB.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBu2DstKB0 = bu2dstkB.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBu2DstKC1 = bu2dstkC.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBu2DstKC2 = bu2dstkC.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBu2DstKC3 = bu2dstkC.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBu2DstKC0 = bu2dstkC.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	double fBd2DstKpi1 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fBd2DstKpi2 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fBd2DstKpi3 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fBd2DstKpi0 = bd2dstkpi.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

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

	double fLb2DpK1 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bd")).getVal();
	double fLb2DpK2 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("Bs")).getVal();
	double fLb2DpK3 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("sideband")).getVal();
	double fLb2DpK0 = lb2dpk.createIntegral(RooArgSet(B_D0_B_CM_M),RooFit::NormSet(B_D0_B_CM_M),RooFit::Range("full")).getVal();

	cout << "Mean: " << sigmean.getVal() << ", sigma: " << sigsigma.getVal() << endl;
	cout << endl;

	Double_t nBd1        = nsig.getVal()*fBd1/fBd0;               
        Double_t nBs1        = nsig2.getVal()*fBs1/fBs0;              
        Double_t nComb1      = nbkg.getVal()*fComb1/fComb0;           
        Double_t nBu2DstKA1  = anbkg2.getVal()*fBu2DstKA1/fBu2DstKA0;           
        Double_t nBu2DstKB1  = bnbkg2.getVal()*fBu2DstKB1/fBu2DstKB0;           
        Double_t nBu2DstKC1  = cnbkg2.getVal()*fBu2DstKB1/fBu2DstKB0;           
        Double_t nBu2DstKD1  = dnbkg2.getVal()*fBu2DstKC1/fBu2DstKC0;           
        Double_t nBu2DstKE1  = enbkg2.getVal()*fBu2DstKC1/fBu2DstKC0;           
        Double_t nBd2DstKpi1 = nbkg3.getVal()*fBd2DstKpi1/fBd2DstKpi0;
        Double_t nBd2Dpipi1  = nbkg4.getVal()*fBd2Dpipi1/fBd2Dpipi0;  
        Double_t nLb2DpK1    = nbkg5.getVal()*fLb2DpK1/fLb2DpK0;      
        Double_t nBd2DKK1    = nbkg6.getVal()*fBd2DKK1/fBd2DKK0;      
        Double_t nBs2DKK1    = nbkg7.getVal()*fBs2DKK1/fBs2DKK0;      

	Double_t nBd2        = nsig.getVal()*fBd2/fBd0;               
        Double_t nBs2        = nsig2.getVal()*fBs2/fBs0;              
        Double_t nComb2     =  nbkg.getVal()*fComb2/fComb0;           
        Double_t nBu2DstKA2  = anbkg2.getVal()*fBu2DstKA2/fBu2DstKA0;           
        Double_t nBu2DstKB2  = bnbkg2.getVal()*fBu2DstKB2/fBu2DstKB0;           
        Double_t nBu2DstKC2  = cnbkg2.getVal()*fBu2DstKB2/fBu2DstKB0;           
        Double_t nBu2DstKD2  = dnbkg2.getVal()*fBu2DstKC2/fBu2DstKC0;           
        Double_t nBu2DstKE2  = enbkg2.getVal()*fBu2DstKC2/fBu2DstKC0;           
        Double_t nBd2DstKpi2 = nbkg3.getVal()*fBd2DstKpi2/fBd2DstKpi0;
        Double_t nBd2Dpipi2  = nbkg4.getVal()*fBd2Dpipi2/fBd2Dpipi0;  
        Double_t nLb2DpK2    = nbkg5.getVal()*fLb2DpK2/fLb2DpK0;      
        Double_t nBd2DKK2    = nbkg6.getVal()*fBd2DKK2/fBd2DKK0;      
        Double_t nBs2DKK2    = nbkg7.getVal()*fBs2DKK2/fBs2DKK0;      

	Double_t nBd3        = nsig.getVal()*fBd3/fBd0;               
        Double_t nBs3        = nsig2.getVal()*fBs3/fBs0;              
        Double_t nComb3     =  nbkg.getVal()*fComb3/fComb0;           
        Double_t nBu2DstKA3  = anbkg2.getVal()*fBu2DstKA3/fBu2DstKA0;           
        Double_t nBu2DstKB3  = bnbkg2.getVal()*fBu2DstKB3/fBu2DstKB0;           
        Double_t nBu2DstKC3  = cnbkg2.getVal()*fBu2DstKB3/fBu2DstKB0;           
        Double_t nBu2DstKD3  = dnbkg2.getVal()*fBu2DstKC3/fBu2DstKC0;           
        Double_t nBu2DstKE3  = enbkg2.getVal()*fBu2DstKC3/fBu2DstKC0;           
        Double_t nBd2DstKpi3 = nbkg3.getVal()*fBd2DstKpi3/fBd2DstKpi0;
        Double_t nBd2Dpipi3  = nbkg4.getVal()*fBd2Dpipi3/fBd2Dpipi0;  
        Double_t nLb2DpK3    = nbkg5.getVal()*fLb2DpK3/fLb2DpK0;      
        Double_t nBd2DKK3    = nbkg6.getVal()*fBd2DKK3/fBd2DKK0;      
        Double_t nBs2DKK3    = nbkg7.getVal()*fBs2DKK3/fBs2DKK0;      

	cout << "Bd window"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd1*afrac->getValV()          , nBd1*bfrac->getValV()          , nBd1*cfrac->getValV()          , nBd1*dfrac->getValV()          , nBd1*efrac->getValV());
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBs1*afrac->getValV()          , nBs1*bfrac->getValV()          , nBs1*cfrac->getValV()          , nBs1*dfrac->getValV()          , nBs1*efrac->getValV());
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nComb1*abkgfrac->getValV()     , nComb1*bbkgfrac->getValV()     , nComb1*cbkgfrac->getValV()     , nComb1*dbkgfrac->getValV()     , nComb1*ebkgfrac->getValV());
	printf("Bu2DstK    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBu2DstKA1                     , nBu2DstKB1                     , nBu2DstKC1                     , nBu2DstKD1                     , nBu2DstKE1);
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2DstKpi1*aprfrac->getValV() , nBd2DstKpi1*bprfrac->getValV() , nBd2DstKpi1*cprfrac->getValV() , nBd2DstKpi1*dprfrac->getValV() , nBd2DstKpi1*eprfrac->getValV());
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2Dpipi1*afrac->getValV()    , nBd2Dpipi1*bfrac->getValV()    , nBd2Dpipi1*cfrac->getValV()    , nBd2Dpipi1*dfrac->getValV()    , nBd2Dpipi1*efrac->getValV());
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nLb2DpK1*afrac->getValV()      , nLb2DpK1*bfrac->getValV()      , nLb2DpK1*cfrac->getValV()      , nLb2DpK1*dfrac->getValV()      , nLb2DpK1*efrac->getValV());
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBd2DKK1*afrac->getValV()      , nBd2DKK1*bfrac->getValV()      , nBd2DKK1*cfrac->getValV()      , nBd2DKK1*dfrac->getValV()      , nBd2DKK1*efrac->getValV());
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBs2DKK1*afrac->getValV()      , nBs2DKK1*bfrac->getValV()      , nBs2DKK1*cfrac->getValV()      , nBs2DKK1*dfrac->getValV()      , nBs2DKK1*efrac->getValV());
	cout << endl;

	cout << "Bs window"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2           , bfrac->getValV()*nBd2           , cfrac->getValV()*nBd2           , dfrac->getValV()*nBd2           , efrac->getVal()*nBd2          );
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2           , bfrac->getValV()*nBs2           , cfrac->getValV()*nBs2           , dfrac->getValV()*nBs2           , efrac->getVal()*nBs2          );
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", abkgfrac->getValV()*nComb2      , bbkgfrac->getValV()*nComb2      , cbkgfrac->getValV()*nComb2      , dbkgfrac->getValV()*nComb2      , ebkgfrac->getVal()*nComb2     );
	printf("Bu2DstK    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBu2DstKA2                     , nBu2DstKB2                     , nBu2DstKC2                     , nBu2DstKD2                     , nBu2DstKE2);
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBd2DstKpi2  , bprfrac->getValV()*nBd2DstKpi2  , cprfrac->getValV()*nBd2DstKpi2  , dprfrac->getValV()*nBd2DstKpi2  , eprfrac->getVal()*nBd2DstKpi2);
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2Dpipi2     , bfrac->getValV()*nBd2Dpipi2     , cfrac->getValV()*nBd2Dpipi2     , dfrac->getValV()*nBd2Dpipi2     , efrac->getVal()*nBd2Dpipi2    );
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2DpK2       , bfrac->getValV()*nLb2DpK2       , cfrac->getValV()*nLb2DpK2       , dfrac->getValV()*nLb2DpK2       , efrac->getVal()*nLb2DpK2      );
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2DKK2       , bfrac->getValV()*nBd2DKK2       , cfrac->getValV()*nBd2DKK2       , dfrac->getValV()*nBd2DKK2       , efrac->getVal()*nBd2DKK2      );
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2DKK2       , bfrac->getValV()*nBs2DKK2       , cfrac->getValV()*nBs2DKK2       , dfrac->getValV()*nBs2DKK2       , efrac->getVal()*nBs2DKK2      );
	cout << endl;
	cout << endl;

	cout << "sideband"    << endl;
	printf("Bd         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd3           , bfrac->getValV()*nBd3           , cfrac->getValV()*nBd3           , dfrac->getValV()*nBd3           , efrac->getVal()*nBd3          );
	printf("Bs         %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs3           , bfrac->getValV()*nBs3           , cfrac->getValV()*nBs3           , dfrac->getValV()*nBs3           , efrac->getVal()*nBs3          );
	printf("Comb       %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", abkgfrac->getValV()*nComb3      , bbkgfrac->getValV()*nComb3      , cbkgfrac->getValV()*nComb3      , dbkgfrac->getValV()*nComb3      , ebkgfrac->getVal()*nComb3     );
	printf("Bu2DstK    %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", nBu2DstKA3                     , nBu2DstKB3                     , nBu2DstKC3                     , nBu2DstKD3                     , nBu2DstKE3);
	printf("Bd2DstKpi  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", aprfrac->getValV()*nBd2DstKpi3  , bprfrac->getValV()*nBd2DstKpi3  , cprfrac->getValV()*nBd2DstKpi3  , dprfrac->getValV()*nBd2DstKpi3  , eprfrac->getVal()*nBd2DstKpi3);
	printf("Bd2Dpipi   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2Dpipi3     , bfrac->getValV()*nBd2Dpipi3     , cfrac->getValV()*nBd2Dpipi3     , dfrac->getValV()*nBd2Dpipi3     , efrac->getVal()*nBd2Dpipi3    );
	printf("Lb2DpK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nLb2DpK3       , bfrac->getValV()*nLb2DpK3       , cfrac->getValV()*nLb2DpK3       , dfrac->getValV()*nLb2DpK3       , efrac->getVal()*nLb2DpK3      );
	printf("Bd2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBd2DKK3       , bfrac->getValV()*nBd2DKK3       , cfrac->getValV()*nBd2DKK3       , dfrac->getValV()*nBd2DKK3       , efrac->getVal()*nBd2DKK3      );
	printf("Bs2DKK     %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", afrac->getValV()*nBs2DKK3       , bfrac->getValV()*nBs2DKK3       , cfrac->getValV()*nBs2DKK3       , dfrac->getValV()*nBs2DKK3       , efrac->getVal()*nBs2DKK3      );
	cout << endl;
	cout << endl;

	printf("$B^0$ in window       & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ \\\\\n", 
		afrac->getValV()*nBd1, (afrac->getValV()*nBd1/ansig->getVal())*ansig->getPropagatedError(*r),
		bfrac->getValV()*nBd1, (bfrac->getValV()*nBd1/bnsig->getVal())*bnsig->getPropagatedError(*r),
		cfrac->getValV()*nBd1, (cfrac->getValV()*nBd1/cnsig->getVal())*cnsig->getPropagatedError(*r),
		dfrac->getValV()*nBd1, (dfrac->getValV()*nBd1/dnsig->getVal())*dnsig->getPropagatedError(*r),
		efrac->getValV()*nBd1, (efrac->getValV()*nBd1/ensig->getVal())*ensig->getPropagatedError(*r));
	printf("$B^0$ window purity   & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%%\\\\\n",
		100*afrac->getValV()*nBd1/(afrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+abkgfrac->getValV()*nComb1+nBu2DstKA1+aprfrac->getValV()*nBd2DstKpi1),
		100*bfrac->getValV()*nBd1/(bfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+bbkgfrac->getValV()*nComb1+nBu2DstKB1+bprfrac->getValV()*nBd2DstKpi1),
		100*cfrac->getValV()*nBd1/(cfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+cbkgfrac->getValV()*nComb1+nBu2DstKC1+cprfrac->getValV()*nBd2DstKpi1),
		100*dfrac->getValV()*nBd1/(dfrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+dbkgfrac->getValV()*nComb1+nBu2DstKD1+dprfrac->getValV()*nBd2DstKpi1),
		100*efrac->getValV()*nBd1/(efrac->getValV()*(nBd1+nBs1+nBd2Dpipi1+nLb2DpK1+nBd2DKK1+nBs2DKK1)+ebkgfrac->getValV()*nComb1+nBu2DstKE1+eprfrac->getValV()*nBd2DstKpi1));
	printf("$B_s^0$ in window     & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ & $%5.1f \\pm %5.1f$ \\\\\n",
		afrac->getValV()*nBs2, (afrac->getValV()*nBs2/ansig2->getVal())*ansig2->getPropagatedError(*r),
		bfrac->getValV()*nBs2, (bfrac->getValV()*nBs2/bnsig2->getVal())*bnsig2->getPropagatedError(*r),
		cfrac->getValV()*nBs2, (cfrac->getValV()*nBs2/cnsig2->getVal())*cnsig2->getPropagatedError(*r),
		dfrac->getValV()*nBs2, (dfrac->getValV()*nBs2/dnsig2->getVal())*dnsig2->getPropagatedError(*r),
		efrac->getValV()*nBs2, (efrac->getValV()*nBs2/ensig2->getVal())*ensig2->getPropagatedError(*r));
	printf("$B_s^0$ window purity & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%% & %5.1f\\,\\%%\\\\\n",
		100*afrac->getValV()*nBs2/(afrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+abkgfrac->getValV()*nComb2+nBu2DstKA2+aprfrac->getValV()*nBd2DstKpi2),
		100*bfrac->getValV()*nBs2/(bfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+bbkgfrac->getValV()*nComb2+nBu2DstKB2+bprfrac->getValV()*nBd2DstKpi2),
		100*cfrac->getValV()*nBs2/(cfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+cbkgfrac->getValV()*nComb2+nBu2DstKC2+cprfrac->getValV()*nBd2DstKpi2),
		100*dfrac->getValV()*nBs2/(dfrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+dbkgfrac->getValV()*nComb2+nBu2DstKD2+dprfrac->getValV()*nBd2DstKpi2),
		100*efrac->getValV()*nBs2/(efrac->getValV()*(nBd2+nBs2+nBd2Dpipi2+nLb2DpK2+nBd2DKK2+nBs2DKK2)+ebkgfrac->getValV()*nComb2+nBu2DstKE2+eprfrac->getValV()*nBd2DstKpi2));
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

	printf("& $% 7.1f\\pm% -6.1f$ \n", sigmean->getVal(),  sigmean->getError());
	printf("& $% 7.1f\\pm% -6.1f$ \n", sigsigma->getVal(), sigsigma->getError());
	printf("& $% 7.3f\\pm% -6.3f$ \n", frac->getVal(),     frac->getError());
	printf("& $% 7.2f\\pm% -6.2f$ \n", ratio->getVal(),    ratio->getError());
	printf("& $% 7.2f\\pm% -6.2f$ \n", 1e3*p0->getVal(),   1e3*p0->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nsig->getVal(),     nsig->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nsig2->getVal(),    nsig2->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg->getVal(),     nbkg->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg3->getVal(),    nbkg3->getError());
	printf("\n");
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg4->getVal(),    nbkg4->getError());
	printf("\n");
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg5->getVal(),    nbkg5->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg6->getVal(),    nbkg6->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg7->getVal(),    nbkg7->getError());
	printf("& $% 7.0f\\pm% -6.0f$ \n", nbkg2->getVal(),    nbkg2->getError());
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
