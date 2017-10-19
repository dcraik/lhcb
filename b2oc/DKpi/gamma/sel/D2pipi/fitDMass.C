{
gSystem.Load("libRooFit");
gROOT.SetStyle("Plain");
//gStyle.SetOptStat(1111);

gStyle->SetOptStat(0000);

TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2pipi_Scaled/B2D0Kpi_D02pipi_selBd_Dsidebands_D50_DK52_PID3_NND2pipi_addIMs.root");
TTree * DecayTree = dynamic_cast<TTree*>(file->Get("DecayTree"));

RooRealVar D_M("D_M","; m(#pi#pi) (MeV/c^{2}); Candidates / 0.9 MeV/c^{2}",  1770,1960);

RooDataSet * data1 = new RooDataSet("data1", "", DecayTree, RooArgSet(D_M));

// D Gaussian 
RooRealVar sigmean("sigmean","",1864.,1844.,1884.,"MeV/c^{2}");
RooRealVar sigsigma("#sigsigma","",5.,0.0,30.0,"MeV/c^{2}");
RooRealVar ratio("ratio","",2.,1.0,3.0);
RooProduct sigsigma2("#sigsigma2","",RooArgSet(sigsigma,ratio));
RooRealVar frac("frac","",0.7,0.0,1.0);

RooGaussian sig( "sig", "", D_M, sigmean, sigsigma );
RooGaussian sig2( "sig2", "", D_M, sigmean, sigsigma2 );
//RooAddPdf sig("sig","",RooArgList(sig1,sig2),RooArgList(frac));

//# Flat Background - Combinatorial
RooRealVar p0("p0","",-0.0001,-0.001,-0.0);
RooPolynomial comb("comb","",D_M,RooArgList(p0));

// D->Kpi Gaussian
RooRealVar Kpimean("Kpimean","",1764.,1714.,1814.,"MeV/c^{2}");
RooRealVar Kpisigma("#Kpisigma","",20.,0.0,50.0,"MeV/c^{2}");
RooRealVar ratio("ratio","",2.,1.0,3.0);
RooProduct Kpisigma2("#Kpisigma2","",RooArgSet(Kpisigma,ratio));
RooRealVar frac("frac","",0.7,0.0,1.0);

RooGaussian Kpi( "Kpi", "", D_M, Kpimean, Kpisigma );
RooGaussian Kpi2( "Kpi2", "", D_M, Kpimean, Kpisigma2 );
//RooAddPdf Kpi("Kpi","",RooArgList(Kpi1,Kpi2),RooArgList(frac));

// Number of signal & background events
RooRealVar nsig("nsig","", 5000.,-1000,20000,"Events");
RooRealVar nbkg("nbkg","",15000.,-1000,30000,"Events");
RooRealVar nKpi("nKpi","",20000.,-1000,40000,"Events");

RooAddPdf full_RF_PDF("full_RF_PDF","RF PDF of everything",RooArgList(sig
			                                             ,comb
								     ,Kpi
								     ),
		                                           RooArgList(nsig
								     ,nbkg
								     ,nKpi
								     ));

//# Do the fit on REFITTED Mass
full_RF_PDF->fitTo(*data1,RooFit::Extended());

TCanvas * can = new TCanvas("can","Mass fits Data",800,600);
can.SetLeftMargin(0.16);
D_Plot = D_M->frame(100);
D_Plot->SetTitle("");
D_Plot->GetYaxis()->SetTitle("Candidates / 0.9 MeV/c^{2}");
D_Plot->GetXaxis()->SetTitle("m(#pi#pi) (MeV/c^{2})");
D_Plot->GetYaxis()->SetTitleOffset(1.6);

data1->plotOn(D_Plot);
full_RF_PDF->plotOn(D_Plot);
full_RF_PDF->plotOn(D_Plot, RooFit::Components("comb"), RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
full_RF_PDF->plotOn(D_Plot, RooFit::Components("Kpi"),  RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
full_RF_PDF->plotOn(D_Plot, RooFit::Components("sig"),  RooFit::LineStyle(kDashed));
D_Plot->Draw();

can->SaveAs("DMass.pdf");
}
