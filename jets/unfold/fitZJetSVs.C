#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "RooAbsDataStore.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooStats/SPlot.h"

#include "RooPromptShape.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"

#include <boost/progress.hpp>

enum svFitType{
	fitBoth,
	fitMCorr,
	fitNTrk
};

//globals to save passing these around
svFitType whichSVFit(fitBoth);
TString savedir("output-fitZj-statOnly-24bins");

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
TString charmSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_14X.root";
TString beautySimFile = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_15X.root";
TString lightSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_101.root";
TString dataFile      = "/tmp/dcraik/skimmed.root";

// function to make plot of 1D projection of the fit
void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
	TCanvas c1;
	c1.SetBottomMargin(0.19);
	RooPlot* plot = var.frame(min, max, nbins);
	dh->plotOn( plot );
	int iCol(0);
	int colours[6] = {kBlue, kRed, kGreen+2, kMagenta, kOrange, kGray};
	int fills[6] = {3245, 3454, 3644, 3205, 3495, 3690};

	//RooCmdArg normCmd = RooFit::NormRange("FIT");
	RooCmdArg normCmd = RooFit::Normalization(1.);

	for (std::vector<std::string>::iterator it = bkg_pdfs.begin(); it != bkg_pdfs.end(); ++it){
		if(iCol>=6) iCol=0;
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(9) );
		++iCol;
	}

	for (std::vector<std::string>::iterator it = sig_pdfs.begin(); it != sig_pdfs.end(); ++it){
		if(iCol>=6) iCol=0;
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(kDashed) );
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::VLines(), RooFit::FillStyle(fills[iCol]), RooFit::FillColor( colours[iCol] ), RooFit::DrawOption("LF") );
		++iCol;
	}

	pdf.plotOn( plot, normCmd, RooFit::LineColor(kBlack) );
	dh->plotOn( plot );
	plot->SetXTitle(title);
	plot->SetTitle("");
	plot->GetXaxis()->SetLabelOffset(0.02);
	plot->GetXaxis()->SetTitleOffset(1.18);
	plot->Draw();

	c1.SaveAs(savedir+"/"+name+"_fit.pdf");
	c1.SaveAs(savedir+"/"+name+"_fit.png");
	c1.SetLogy();
	c1.SaveAs(savedir+"/"+name+"_fit_log.pdf");

	delete plot;
}

// function to print parameters to a file
void printParams(TString file, const RooArgList& params) {
	FILE * pFile = fopen((savedir+"/"+file).Data(), "w");

	int nPar = params.getSize();

	for ( int i=0; i<nPar; ++i) {
		RooRealVar* par = dynamic_cast<RooRealVar*>(params.at(i));

		TString title = par->getTitle();
		float value = par->getValV();
		float error = par->getError();

		int exponent(0);
		int exponentErr(0);
		int precision(0);

		while(TMath::Abs(value) < TMath::Power(10,exponent)) exponent-=1;
		while(TMath::Abs(value) > TMath::Power(10,exponent+1)) exponent+=1;
		while(TMath::Abs(error) < TMath::Power(10,exponentErr)) exponentErr-=1;
		while(TMath::Abs(error) > TMath::Power(10,exponentErr+1)) exponentErr+=1;

		if(error < 3.5*TMath::Power(10,exponentErr)) precision = 1-exponentErr;
		else precision = -exponentErr;
		if(precision<0) precision=0;

		if(exponent<-2) {
			title+=" ($10^{";
			title+=exponent;
			title+="}$)";

			precision+=exponent;

			value/=TMath::Power(10,exponent);
			error/=TMath::Power(10,exponent);
		}

		TString format = "%-30s";
		format+="\t& $% 6.";
		format+=precision;
		format+="f \\pm % 6.";
		format+=precision;
		format+="f$ \\\\\n";

		fprintf(pFile, format, title.Data(), value, error);
	}
	fclose(pFile);
}


TString fitSV(double& NC, double& eC, double& NJ, double& eJ, double minPT=20000, double maxPT=30000, double minY=3.0, double maxY=3.25){
	std::cout << "INFO : fitting SV - 1D corrected mass and Ntrk fit, pT in (" << minPT << "," << maxPT << "); y(Z) in (" << minY << "," << maxY << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT; ptStr+="_";
	ptStr+=minY; ptStr+="-"; ptStr+=maxY;

	TFile* f0 = TFile::Open(lightSimFile);
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f4 = TFile::Open(charmSimFile);
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));

	TFile* f5 = TFile::Open(beautySimFile);
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	TFile* fd = TFile::Open(dataFile);
	TTree* td = dynamic_cast<TTree*>(fd->Get("T"));


	if(!t0 || !t4 || !t5 || !td) return "FAIL";

	double mmin(500.), mmax(5000.);
	//int nmbins(46.);
	int nmbins(20.);
	int nbin(nmbins+4);
	double scaleNtrk=4./nmbins;//scale down Ntrks part for nicer plots
	double min(0.), max(nbin);
	double scale(1./(1.+scaleNtrk));//scale to correct for ntrks part when extracting yields
	if(whichSVFit==fitMCorr) {
		nbin=nmbins;
		max=nmbins;
		scale=1.;
	} else if(whichSVFit==fitNTrk) {
		nbin=4;
		min=nmbins;
		scale=1.;
		scaleNtrk=1.;
	}

	TH1D svmcorsvn_0("svmcorsvn_0","",nbin,min,max);
	TH1D svmcorsvn_4("svmcorsvn_4","",nbin,min,max);
	TH1D svmcorsvn_5("svmcorsvn_5","",nbin,min,max);
	TH1D svmcorsvn_d("svmcorsvn_d","",nbin,min,max);

	svmcorsvn_4.Sumw2();
	svmcorsvn_5.Sumw2();
	svmcorsvn_d.Sumw2();

	std::vector<double>* SVN = new std::vector<double>();
	std::vector<double>* SVMCor = new std::vector<double>();
	double JetPT, JetTruePT, ZE, ZPZ;

	t0->SetBranchAddress("SVMCor", &SVMCor);
	t4->SetBranchAddress("SVMCor", &SVMCor);
	t5->SetBranchAddress("SVMCor", &SVMCor);
	td->SetBranchAddress("SVMCor", &SVMCor);
	t0->SetBranchAddress("SVN", &SVN);
	t4->SetBranchAddress("SVN", &SVN);
	t5->SetBranchAddress("SVN", &SVN);
	td->SetBranchAddress("SVN", &SVN);
	t0->SetBranchAddress("JetPT", &JetPT);
	t4->SetBranchAddress("JetPT", &JetPT);
	t5->SetBranchAddress("JetPT", &JetPT);
	td->SetBranchAddress("JetPT", &JetPT);
	td->SetBranchAddress("ZE", &ZE);
	td->SetBranchAddress("ZPZ", &ZPZ);
	t4->SetBranchAddress("JetTruePT", &JetTruePT);
	t5->SetBranchAddress("JetTruePT", &JetTruePT);

	//get MC/backwards efficiencies for estimating total yields
	int eff0denom(0.), eff4denom(0.), eff5denom(0.);
	int eff0num(0.), eff4num(0.), eff5num(0.);
	double weight(1.);

	unsigned int npt=4;
	double* ptInputWeightBins  = new double[npt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* mcWeights4 = new TH1D("mcWeights4","",npt,ptInputWeightBins);
	TH1D* mcWeights5 = new TH1D("mcWeights5","",npt,ptInputWeightBins);
	TH1D* mcWeightsD = new TH1D("mcWeightsD","",npt,ptInputWeightBins);//average to use for mixed sample
	mcWeights4->SetBinContent(1, 1.);
	mcWeights4->SetBinContent(2, 0.3);
	mcWeights4->SetBinContent(3, 0.12);
	mcWeights4->SetBinContent(4, 0.008);
	mcWeights5->SetBinContent(1, 1.);
	mcWeights5->SetBinContent(2, 0.2);
	mcWeights5->SetBinContent(3, 0.08);
	mcWeights5->SetBinContent(4, 0.006);
	mcWeightsD->SetBinContent(1, 1.);
	mcWeightsD->SetBinContent(2, 0.25);
	mcWeightsD->SetBinContent(3, 0.10);
	mcWeightsD->SetBinContent(4, 0.007);

	for(int i=0; i<t0->GetEntries(); ++i) {
		t0->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff0denom;
		if(SVN->size()<1) continue;
		++eff0num;
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_0.Fill(SVN->at(0)+nmbins-2,scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_0.Fill(nmbins+3,scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_0.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin));
		else if(SVMCor->at(0)>mmax) svmcorsvn_0.Fill(nmbins-1);
	}
	std::cout << svmcorsvn_0.Integral(1,nmbins) << "\t" << svmcorsvn_0.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	for(int i=0; i<t4->GetEntries(); ++i) {
		t4->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff4denom;
		if(SVN->size()<1) continue;
		++eff4num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights4->GetBinContent(mcWeights4->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_4.Fill(SVN->at(0)+nmbins-2,weight*scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_4.Fill(nmbins+3,weight*scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_4.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin),weight);
		else if(SVMCor->at(0)>mmax) svmcorsvn_4.Fill(nmbins-1,weight);
	}
	std::cout << svmcorsvn_4.Integral(1,nmbins) << "\t" << svmcorsvn_4.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	for(int i=0; i<t5->GetEntries(); ++i) {
		t5->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff5denom;
		if(SVN->size()<1) continue;
		++eff5num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights5->GetBinContent(mcWeights5->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_5.Fill(SVN->at(0)+nmbins-2,weight*scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_5.Fill(nmbins+3,weight*scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_5.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin),weight);
		else if(SVMCor->at(0)>mmax) svmcorsvn_5.Fill(nmbins-1,weight);
	}
	std::cout << svmcorsvn_5.Integral(1,nmbins) << "\t" << svmcorsvn_5.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	NJ=0;
	for(int i=0; i<td->GetEntries(); ++i) {
		td->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		if(0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ))<minY || 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ))>maxY) continue;
		++NJ;
		if(SVN->size()<1) continue;
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_d.Fill(SVN->at(0)+nmbins-2,scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_d.Fill(nmbins+3,scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_d.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin));
		else if(SVMCor->at(0)>mmax) svmcorsvn_d.Fill(nmbins-1);
	}
	std::cout << svmcorsvn_d.Integral(1,nmbins) << "\t" << svmcorsvn_d.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;
	eJ = TMath::Sqrt(NJ);

	// -- variables from datasets
	RooRealVar SVComb(  "SVComb",  "SVComb",  min, max,  ""); 
	SVComb.setBins(nbin);

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.4*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.3*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());

	RooDataHist histSVMSVNB("histCombB", "histCombB", RooArgList(SVComb), &svmcorsvn_5);
	RooDataHist histSVMSVNC("histCombC", "histCombC", RooArgList(SVComb), &svmcorsvn_4);
	RooDataHist histSVMSVNQ("histCombQ", "histCombQ", RooArgList(SVComb), &svmcorsvn_0);

	// -- simulation PDFs for each category
	RooHistPdf pdfB( "pdfB", "pdfB", RooArgSet(SVComb), histSVMSVNB ); 
	RooHistPdf pdfC( "pdfC", "pdfC", RooArgSet(SVComb), histSVMSVNC ); 
	RooHistPdf pdfQ( "pdfQ", "pdfQ", RooArgSet(SVComb), histSVMSVNQ ); 

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(pdfB, pdfC, pdfQ), RooArgList(yieldB, yieldC, yieldQ) );

	// -- add all feature observables to dataset
	RooArgSet obs;
	obs.add(SVComb);

	RooDataHist dh("dh", "dh", RooArgList(SVComb), &svmcorsvn_d);

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(savedir+"/SVComb_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
	gSystem->RedirectOutput(0);
	///*RooFitResult * result =*/ data_pdf.chi2FitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	TString plotName = "SVComb_"; plotName+=ptStr;
	plotFit(SVComb, min, max, nbin, &dh, data_pdf, sig_pdfs, bkg_pdfs, plotName, "M_{cor}, N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_comb_"; paramsName+=minPT; paramsName+="-"; paramsName+=maxPT; paramsName+=".dat";
	printParams(paramsName,params);

	double NQ, NB, Ntot;
	double eQ, eB;

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV()*scale;
	NC = yieldC.getValV()*scale;
	NQ = yieldQ.getValV()*scale;
	eB = yieldB.getError()*scale;
	eC = yieldC.getError()*scale;
	eQ = yieldQ.getError()*scale;
	Ntot*=scale;
	double eff5 = static_cast<double>(eff5num)/eff5denom;
	double eff4 = static_cast<double>(eff4num)/eff4denom;
	double eff0 = static_cast<double>(eff0num)/eff0denom;
	double erreff5 = TMath::Sqrt(static_cast<double>(eff5num))/eff5denom;
	double erreff4 = TMath::Sqrt(static_cast<double>(eff4num))/eff4denom;
	double erreff0 = TMath::Sqrt(static_cast<double>(eff0num))/eff0denom;
	double corr5 = NB/eff5;
	double corr4 = NC/eff4;
	double corr0 = NQ/eff0;
	double errcorr5 = corr5 * TMath::Sqrt(eB*eB/NB/NB + erreff5*erreff5/eff5/eff5);
	double errcorr4 = corr4 * TMath::Sqrt(eC*eC/NC/NC + erreff4*erreff4/eff4/eff4);
	double errcorr0 = corr0 * TMath::Sqrt(eQ*eQ/NQ/NQ + erreff0*erreff0/eff0/eff0);
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ, corr5, errcorr5, corr4, errcorr4, corr0, errcorr0);

	return true;
}

int main() {
	gSystem->Exec("mkdir -p "+savedir);

	double nJ(0.), eJ(0.), nC(0.), eC(0.);

	const int npt(3);
	double ptBounds[npt+1] = {15000.,20000.,30000.,100000.};
	//const int nzy(2);
	//double zyBounds[nzy+1] = {2.0,3.0,4.5};
	const int nzy(5);
	double zyBounds[nzy+1] = {2.0,2.5,3.0,3.5,4.0,4.5};

	TH2D hC("hC","",npt,ptBounds,nzy,zyBounds);
	TH2D hJ("hJ","",npt,ptBounds,nzy,zyBounds);
	TH2D hE("hE","",npt,ptBounds,nzy,zyBounds);
	TH2D hR("hR","",npt,ptBounds,nzy,zyBounds);

	hC.Sumw2();
	hJ.Sumw2();
	hE.Sumw2();
	hR.Sumw2();

	for(int j=1; j<=nzy; ++j) {
		hE.SetBinContent(1,j,0.191);
		hE.SetBinError(1,j,0.);//0.021);
		hE.SetBinContent(2,j,0.237);
		hE.SetBinError(2,j,0.);//0.016);
		hE.SetBinContent(3,j,0.241);
		hE.SetBinError(3,j,0.);//0.022);
	}

	for(int i=1; i<=npt; ++i) {
		for(int j=1; j<=nzy; ++j) {
			if(!fitSV(nC,eC,nJ,eJ,ptBounds[i-1],ptBounds[i],zyBounds[j-1],zyBounds[j])) return 1;
			hC.SetBinContent(i,j,nC);
			hC.SetBinError(  i,j,eC);
			hJ.SetBinContent(i,j,nJ);
			hJ.SetBinError(  i,j,eJ);
		}
	}

	hR.Divide(&hC,&hJ);
	hR.Divide(&hE);

	TCanvas c;

	TH1D* h1520 = hR.ProjectionY("h1520",1,1);
	TH1D* h2030 = hR.ProjectionY("h2030",2,2);
	TH1D* h30up = hR.ProjectionY("h30up",3,3);

	h30up->SetLineColor(kRed);
	h30up->SetMarkerColor(kRed);
	h1520->SetLineColor(kGreen+2);
	h1520->SetMarkerColor(kGreen+2);
	h2030->SetLineColor(kBlue);
	h2030->SetMarkerColor(kBlue);
	h30up->SetMinimum(0.);
	h30up->SetMaximum(0.15);
	gStyle->SetOptStat(0);
	h30up->GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	h30up->GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{c}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	h30up->Draw("E1 P");
	h2030->Draw("E1 P same");
	h1520->Draw("E1 P same");
	c.SaveAs(savedir+"/fracCharmZj.pdf");

	gStyle->SetOptStat(0);
	hC.Draw("colz text45");
	c.SaveAs(savedir+"/charmYields.pdf");
	hJ.Draw("colz text45");
	c.SaveAs(savedir+"/totalYields.pdf");
	hR.Draw("colz text45");
	c.SaveAs(savedir+"/charmFractions.pdf");

	for(int i=1; i<=npt; ++i) {
		for(int j=1; j<=nzy; ++j) {
			std::cout << ptBounds[i-1] << "-" << ptBounds[i] << "," << zyBounds[j-1] << "-" << zyBounds[j] << ":\t" << hJ.GetBinContent(i,j) << "+/-" << hJ.GetBinError(i,j) << "\t" << hC.GetBinContent(i,j) << "+/-" << hC.GetBinError(i,j) << "\t" << hE.GetBinContent(i,j) << "+/-" << hE.GetBinError(i,j) << "\t" << hR.GetBinContent(i,j) << "+/-" << hR.GetBinError(i,j) << std::endl;
		}
	}

	TFile* fout = TFile::Open(savedir+"/hists.root","RECREATE");
	hC.Write();
	hJ.Write();
	hE.Write();
	hR.Write();
	fout->Close();

	return 0;
}
