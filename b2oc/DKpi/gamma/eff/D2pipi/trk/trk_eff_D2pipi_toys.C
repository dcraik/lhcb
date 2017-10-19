#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <TColor.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>

#include <Laura++/LauKinematics.hh>
#include <Laura++/LauDaughters.hh>

void trk_Eff(Int_t seed){

	TString xTitle="m'";
	TString yTitle="#theta'";
	Int_t nbins(0);

	nbins = 10;

	//double mB(5.27958), mBs(5.3662998046875), mD(1.86484), mK(0.493677001953125000), mPi(0.139570175170830879);
	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);

	LauKinematics* kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * hpiCorr     = new TH2D("hpiCorr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiNoCorr   = new TH2D("hpiNoCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiEff      = new TH2D("hpiEff",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiErr      = new TH2D("hpiErr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKCorr      = new TH2D("hKCorr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKNoCorr    = new TH2D("hKNoCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKEff       = new TH2D("hKEff",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKErr       = new TH2D("hKErr",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pCorr   = new TH2D("hD0pCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pNoCorr = new TH2D("hD0pNoCorr","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pEff    = new TH2D("hD0pEff",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pErr    = new TH2D("hD0pErr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mCorr    = new TH2D("hD0mCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mNoCorr  = new TH2D("hD0mNoCorr", "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mEff     = new TH2D("hD0mEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mErr     = new TH2D("hD0mErr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hAllEff     = new TH2D("hAllEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);

	TH2D * petaPi  = new TH2D("petaPi", "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaK   = new TH2D("petaK",  "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDpi = new TH2D("petaDpi","",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDK  = new TH2D("petaDK", "",50,0.0,100.0,20,1.0,5.0);
	
	hpiNoCorr->Sumw2();
	hpiCorr->Sumw2();
	hKNoCorr->Sumw2();
	hKCorr->Sumw2();
	hD0pNoCorr->Sumw2();
	hD0pCorr->Sumw2();
	hD0mNoCorr->Sumw2();
	hD0mCorr->Sumw2();

	TFile * file = TFile::Open("D2pipi_Bd2D0Kpi_trackWeight.root");

	TTree * tree = dynamic_cast<TTree*>( file->Get("DecayTree") );
	Int_t nEntries = tree->GetEntries();

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0m_TRUEID,D0p_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0m_MC_GD_MOTHER_KEY,D0p_MC_GD_MOTHER_KEY;
	Float_t NN;
	Double_t D0m_ProbNNk, D0m_ProbNNpi, D0p_ProbNNk, D0p_ProbNNpi;
	Double_t piWeight(0.), KWeight(0.), D0pWeight(0.), D0mWeight(0.);
	Double_t piWeightErr(0.), KWeightErr(0.), D0pWeightErr(0.), D0mWeightErr(0.);
	
	tree->SetBranchAddress("Bd_CM_m13Sq",             &Bd_CM_m13Sq);
	tree->SetBranchAddress("Bd_CM_m23Sq",             &Bd_CM_m23Sq);
	tree->SetBranchAddress("Bd_CM_m12Sq",             &Bd_CM_m12Sq);
	tree->SetBranchAddress("B_TRUEID",                &B_TRUEID);
	tree->SetBranchAddress("D_TRUEID",                &D_TRUEID);
	tree->SetBranchAddress("D0m_TRUEID",              &D0m_TRUEID);
	tree->SetBranchAddress("D0p_TRUEID",             &D0p_TRUEID);
	tree->SetBranchAddress("K_TRUEID",                &K_TRUEID);
	tree->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);

	tree->SetBranchAddress("D0m_MC_GD_MOTHER_KEY",    &D0m_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("D0p_MC_GD_MOTHER_KEY",   &D0p_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("K_MC_MOTHER_KEY",         &K_MC_MOTHER_KEY);
	tree->SetBranchAddress("pi_MC_MOTHER_KEY",        &pi_MC_MOTHER_KEY);

	tree->SetBranchAddress("NN",                      &NN);
	tree->SetBranchAddress("D0m_ProbNNk",             &D0m_ProbNNk);
	tree->SetBranchAddress("D0m_ProbNNpi",            &D0m_ProbNNpi);
	tree->SetBranchAddress("D0p_ProbNNk",            &D0p_ProbNNk);
	tree->SetBranchAddress("D0p_ProbNNpi",           &D0p_ProbNNpi);

	tree->SetBranchAddress("piWeight",                &piWeight);
	tree->SetBranchAddress("KWeight",                 &KWeight);
	tree->SetBranchAddress("D0pWeight",              &D0pWeight);
	tree->SetBranchAddress("D0mWeight",               &D0mWeight);

	tree->SetBranchAddress("piWeightErr",             &piWeightErr);
	tree->SetBranchAddress("KWeightErr",              &KWeightErr);
	tree->SetBranchAddress("D0pWeightErr",           &D0pWeightErr);
	tree->SetBranchAddress("D0mWeightErr",            &D0mWeightErr);

	Double_t pi_PX(0);
	Double_t pi_PY(0);
	Double_t pi_PZ(0);
	Double_t K_PX(0);
	Double_t K_PY(0);
	Double_t K_PZ(0);
	Double_t D0p_PX(0);
	Double_t D0p_PY(0);
	Double_t D0p_PZ(0);
	Double_t D0m_PX(0);
	Double_t D0m_PY(0);
	Double_t D0m_PZ(0);

	Double_t pi_P(0);
	Double_t pi_eta(0);
	Double_t K_P(0);
	Double_t K_eta(0);
	Double_t D0p_P(0);
	Double_t D0p_eta(0);
	Double_t D0m_P(0);
	Double_t D0m_eta(0);

	Double_t piCorr(0), KCorr(0.), D0pCorr(0.), D0mCorr(0.);
	Double_t piCorrErr(0), KCorrErr(0.), D0pCorrErr(0.), D0mCorrErr(0.);

	// PID for systematics
	tree->SetBranchAddress("pi_TRUEP_X",        &pi_PX);
	tree->SetBranchAddress("pi_TRUEP_Y",        &pi_PY);
	tree->SetBranchAddress("pi_TRUEP_Z",        &pi_PZ);
	tree->SetBranchAddress("K_TRUEP_X",         &K_PX);
	tree->SetBranchAddress("K_TRUEP_Y",         &K_PY);
	tree->SetBranchAddress("K_TRUEP_Z",         &K_PZ);
	tree->SetBranchAddress("D0p_TRUEP_X",      &D0p_PX);
	tree->SetBranchAddress("D0p_TRUEP_Y",      &D0p_PY);
	tree->SetBranchAddress("D0p_TRUEP_Z",      &D0p_PZ);
	tree->SetBranchAddress("D0m_TRUEP_X",       &D0m_PX);
	tree->SetBranchAddress("D0m_TRUEP_Y",       &D0m_PY);
	tree->SetBranchAddress("D0m_TRUEP_Z",       &D0m_PZ);

	// Get Histograms
	TFile * histfile = TFile::Open("ratio2012S20.root");
	TH2F  * hist     = (TH2F*)histfile->Get("Ratio")->Clone();
	Int_t nHistBinsX = hist->GetXaxis()->GetNbins();
	Int_t nHistBinsY = hist->GetYaxis()->GetNbins();

	TRandom3 r(60000+seed);
	for(int j=1; j<=nHistBinsY; ++j) {
		for(int i=1; i<=nHistBinsX; ++i) {
			Float_t mean = hist->GetBinContent(i,j);
			Float_t sigma = hist->GetBinError(i,j);
			hist->SetBinContent(i,j,r.Gaus(mean,sigma));
		}
	}

	// Variables required for bin matching
	Int_t binpi(-1), binK(-1), binD0p(-1), binD0m(-1);

	Double_t piWeight2(0.), KWeight2(0.), D0pWeight2(0.), D0mWeight2(0.);
	Double_t piCorrErr2(0), KCorrErr2(0.), D0pCorrErr2(0.), D0mCorrErr2(0.);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
		  continue;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0m_TRUEID)==211&&abs(D0p_TRUEID)==211) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0m_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0p_MC_GD_MOTHER_KEY) {
				if(NN>-0.8) {
					//if(D0m_ProbNNk*(1-D0m_ProbNNpi)>0.1 && D0p_ProbNNpi*(1-D0p_ProbNNk)>0.1) {
						piCorr   = 0.;
						KCorr    = 0.;
						D0pCorr = 0.;
						D0mCorr  = 0.;
						pi_P       = sqrt(pi_PX * pi_PX + pi_PY * pi_PY + pi_PZ * pi_PZ);
						pi_eta     = 0.5*log( (pi_P + pi_PZ) / (pi_P - pi_PZ) );
						K_P        = sqrt(K_PX  * K_PX  + K_PY  * K_PY  + K_PZ  * K_PZ);
						K_eta      = 0.5*log( (K_P  + K_PZ)  / (K_P  - K_PZ)  );

						D0p_P     = sqrt(D0p_PX * D0p_PX + D0p_PY * D0p_PY + D0p_PZ * D0p_PZ);
						D0p_eta   = 0.5*log( (D0p_P + D0p_PZ) / (D0p_P - D0p_PZ) );
						D0m_P      = sqrt(D0m_PX  * D0m_PX  + D0m_PY  * D0m_PY  + D0m_PZ  * D0m_PZ);
						D0m_eta    = 0.5*log( (D0m_P  + D0m_PZ)  / (D0m_P  - D0m_PZ)  );

						binpi   = hist->FindBin(pi_P/1000.,   pi_eta  );
						binK    = hist->FindBin(K_P/1000.,    K_eta   );
						binD0p = hist->FindBin(D0p_P/1000., D0p_eta);
						binD0m  = hist->FindBin(D0m_P/1000.,  D0m_eta );

						if(binpi   %  (nHistBinsX+2) == 0)            binpi+=1;
						if(binpi   %  (nHistBinsX+2) == nHistBinsX+1) binpi-=1;
						if(binpi   <  (nHistBinsX+2))                 binpi+=(nHistBinsX+2);
						if(binpi   >= (nHistBinsX+2)*(nHistBinsY+1))  binpi-=(nHistBinsX+2);
						if(binK    %  (nHistBinsX+2) == 0)            binK+=1;
						if(binK    %  (nHistBinsX+2) == nHistBinsX+1) binK-=1;
						if(binK    <  (nHistBinsX+2))                 binK+=(nHistBinsX+2);
						if(binK    >= (nHistBinsX+2)*(nHistBinsY+1))  binK-=(nHistBinsX+2);
						if(binD0p %  (nHistBinsX+2) == 0)            binD0p+=1;
						if(binD0p %  (nHistBinsX+2) == nHistBinsX+1) binD0p-=1;
						if(binD0p <  (nHistBinsX+2))                 binD0p+=(nHistBinsX+2);
						if(binD0p >= (nHistBinsX+2)*(nHistBinsY+1))  binD0p-=(nHistBinsX+2);
						if(binD0m  %  (nHistBinsX+2) == 0)            binD0m+=1;
						if(binD0m  %  (nHistBinsX+2) == nHistBinsX+1) binD0m-=1;
						if(binD0m  <  (nHistBinsX+2))                 binD0m+=(nHistBinsX+2);
						if(binD0m  >= (nHistBinsX+2)*(nHistBinsY+1))  binD0m-=(nHistBinsX+2);


						piCorr   = hist->GetBinContent(binpi);
						KCorr    = hist->GetBinContent(binK);
						D0pCorr = hist->GetBinContent(binD0p);
						D0mCorr  = hist->GetBinContent(binD0m);

						piCorrErr   = hist->GetBinError(binpi);
						KCorrErr    = hist->GetBinError(binK);
						D0pCorrErr = hist->GetBinError(binD0p);
						D0mCorrErr  = hist->GetBinError(binD0m);

	                                        piWeight2      = piWeight*piWeight;
						KWeight2       = KWeight*KWeight;
						D0pWeight2    = D0pWeight*D0pWeight;
						D0mWeight2     = D0mWeight*D0mWeight;
	                                        piCorrErr2     = piCorrErr*piCorrErr;
						KCorrErr2      = KCorrErr*KCorrErr;
						D0pCorrErr2   = D0pCorrErr*D0pCorrErr;
						D0mCorrErr2    = D0mCorrErr*D0mCorrErr;

	//					if(weight<0.9) cout << pi_P << "\t" << pi_eta << "\t" << K_P << "\t" << K_eta << "\t" << D0p_P << "\t" << D0p_eta << "\t" << D0m_P << "\t" << D0m_eta << endl;   
	//					if(weight<0.9) cout << hist->FindBin(pi_P/1000.,pi_eta) << "\t" << hist->FindBin(K_P/1000.,K_eta) << "\t" 
	//						            << hist->FindBin(D0p_P/1000.,D0p_eta) << "\t" << hist->FindBin(D0m_P/1000.,D0m_eta ) << endl;   
	//					if(weight<0.9) cout << hist->GetBinContent(hist->FindBin(pi_P/1000.,pi_eta)) << "\t" << hist->GetBinContent(hist->FindBin(K_P/1000.,K_eta)) << "\t" 
	//						            << hist->GetBinContent(hist->FindBin(D0p_P/1000.,D0p_eta)) << "\t" << hist->GetBinContent(hist->FindBin(D0m_P/1000.,D0m_eta )) << endl;   

						hpiNoCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight);
						hpiCorr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight   * piCorr);
						hKNoCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight);
						hKCorr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight    * KCorr);
						hD0pNoCorr->Fill( kinematics->getmPrime(), kinematics->getThetaPrime(), D0pWeight);
						hD0pCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), D0pWeight * D0pCorr);
						hD0mNoCorr->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), D0mWeight);
						hD0mCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), D0mWeight  * D0mCorr);

						hpiErr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), /*piWeightErr2   * piCorr2   +*/ piWeight2   * piCorrErr2);
						hKErr->Fill(       kinematics->getmPrime(), kinematics->getThetaPrime(), /*KWeightErr2    * KCorr2    +*/ KWeight2    * KCorrErr2);
						hD0pErr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), /*D0pWeightErr2 * D0pCorr2 +*/ D0pWeight2 * D0pCorrErr2);
						hD0mErr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), /*D0mWeightErr2  * D0mCorr2  +*/ D0mWeight2  * D0mCorrErr2);

						petaPi->Fill(pi_P/1000.,pi_eta);
						petaK->Fill(K_P/1000.,K_eta);
						petaDpi->Fill(D0p_P/1000.,D0p_eta);
						petaDK->Fill(D0m_P/1000.,D0m_eta);
					//}
				}
			}
		}	
	}

	file->Close();

	hpiEff->Divide(  hpiCorr,  hpiNoCorr);
	hKEff->Divide(   hKCorr,   hKNoCorr);
	hD0pEff->Divide(hD0pCorr,hD0pNoCorr);
	hD0mEff->Divide( hD0mCorr, hD0mNoCorr);

	for (Int_t k=0; k<nbins; ++k) {
	   for (Int_t l=0;l<nbins;++l) {
	      Float_t errPi2 = hpiErr->GetBinContent(k+1,l+1);
	      hpiEff->SetBinError(k+1,l+1,sqrt(errPi2)/hpiNoCorr->GetBinContent(k+1,l+1));

	      Float_t errK2 = hKErr->GetBinContent(k+1,l+1);
	      hKEff->SetBinError(k+1,l+1,sqrt(errK2)/hKNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0p2 = hD0pErr->GetBinContent(k+1,l+1);
	      hD0pEff->SetBinError(k+1,l+1,sqrt(errD0p2)/hD0pNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0m2 = hD0mErr->GetBinContent(k+1,l+1);
	      hD0mEff->SetBinError(k+1,l+1,sqrt(errD0m2)/hD0mNoCorr->GetBinContent(k+1,l+1));
	   }
	}

	hAllEff->Add(hpiEff);
	hAllEff->Multiply(hKEff);
	hAllEff->Multiply(hD0pEff);
	hAllEff->Multiply(hD0mEff);

	TString outFileName("toy/trk_sq_Bd_effs");
	outFileName+=seed;
	outFileName+=".root";
	TFile* outfile = new TFile(outFileName,"RECREATE");
	hAllEff->SetName("efficiency");
	hAllEff->Write();
	outfile->Close();

	delete kinematics;
}

//void trk_Eff_Bs_square() {
//  gROOT->ProcessLine(".L lhcbStyle.C");
//  gSystem->Load( "libEG" );
//  gSystem->Load( "libFitEff.so" );
//  lhcbStyle();
//  gStyle->SetOptStat(0000);
//  gStyle->SetPalette(1,0);
//  const Int_t NRGBs = 5;
//  const Int_t NCont = 255;
//  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00};
//  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51};
//  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00};
//  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00};
//  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
//  gStyle->SetNumberContours(NCont);
//
//  trk_Eff();
//}

int main(int argc, char** argv) {
  Int_t first(0), number(100);

  if(argc>1) {
    number = atoi(argv[1]);
  }
  if(argc>2) {
    first  = atoi(argv[2]);
  }

  gStyle->SetOptStat(0000);
  gStyle->SetPalette(1,0);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  for(Int_t seed=first; seed<(first+number); ++seed) {
	  std::cout << seed << std::endl;
	  trk_Eff(seed);
  }
  return 0;
}
