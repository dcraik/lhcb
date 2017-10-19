#include <TLorentzVector.h>

class LauKinematics;

void trk_Eff(){

	TString xTitle="m'";
	TString yTitle="#theta'";
	Int_t nbins(0);

	nbins = 10;

	//double mB(5.27958), mBs(5.3662998046875), mD(1.86484), mK(0.493677001953125000), mPi(0.139570175170830879);
	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * hpiCorr     = new TH2D("hpiCorr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiNoCorr   = new TH2D("hpiNoCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiEff      = new TH2D("hpiEff",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiErr      = new TH2D("hpiErr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKCorr      = new TH2D("hKCorr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKNoCorr    = new TH2D("hKNoCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKEff       = new TH2D("hKEff",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKErr       = new TH2D("hKErr",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mCorr   = new TH2D("hD0mCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mNoCorr = new TH2D("hD0mNoCorr","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mEff    = new TH2D("hD0mEff",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0mErr    = new TH2D("hD0mErr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pCorr    = new TH2D("hD0pCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pNoCorr  = new TH2D("hD0pNoCorr", "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pEff     = new TH2D("hD0pEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0pErr     = new TH2D("hD0pErr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hAllEff     = new TH2D("hAllEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);

	TH2D * petaPi  = new TH2D("petaPi", "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaK   = new TH2D("petaK",  "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDpi = new TH2D("petaDpi","",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDK  = new TH2D("petaDK", "",50,0.0,100.0,20,1.0,5.0);
	
	hpiNoCorr->Sumw2();
	hpiCorr->Sumw2();
	hKNoCorr->Sumw2();
	hKCorr->Sumw2();
	hD0mNoCorr->Sumw2();
	hD0mCorr->Sumw2();
	hD0pNoCorr->Sumw2();
	hD0pCorr->Sumw2();

	TFile * file = TFile::Open("D2pipi_Bd2D0Kpi_trackWeight.root");

	TTree * tree = dynamic_cast<TTree*>( file->Get("DecayTree") );
	Int_t nEntries = tree->GetEntries();

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0p_TRUEID,D0m_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0p_MC_GD_MOTHER_KEY,D0m_MC_GD_MOTHER_KEY;
	Double_t piWeight(0.), KWeight(0.), D0mWeight(0.), D0pWeight(0.);
	Double_t piWeightErr(0.), KWeightErr(0.), D0mWeightErr(0.), D0pWeightErr(0.);
	
	tree->SetBranchAddress("Bd_CM_m13Sq",             &Bd_CM_m13Sq);
	tree->SetBranchAddress("Bd_CM_m23Sq",             &Bd_CM_m23Sq);
	tree->SetBranchAddress("Bd_CM_m12Sq",             &Bd_CM_m12Sq);
	tree->SetBranchAddress("B_TRUEID",                &B_TRUEID);
	tree->SetBranchAddress("D_TRUEID",                &D_TRUEID);
	tree->SetBranchAddress("D0p_TRUEID",              &D0p_TRUEID);
	tree->SetBranchAddress("D0m_TRUEID",             &D0m_TRUEID);
	tree->SetBranchAddress("K_TRUEID",                &K_TRUEID);
	tree->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);

	tree->SetBranchAddress("D0p_MC_GD_MOTHER_KEY",    &D0p_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("D0m_MC_GD_MOTHER_KEY",   &D0m_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("K_MC_MOTHER_KEY",         &K_MC_MOTHER_KEY);
	tree->SetBranchAddress("pi_MC_MOTHER_KEY",        &pi_MC_MOTHER_KEY);

	tree->SetBranchAddress("piWeight",                &piWeight);
	tree->SetBranchAddress("KWeight",                 &KWeight);
	tree->SetBranchAddress("D0mWeight",              &D0mWeight);
	tree->SetBranchAddress("D0pWeight",               &D0pWeight);

	tree->SetBranchAddress("piWeightErr",             &piWeightErr);
	tree->SetBranchAddress("KWeightErr",              &KWeightErr);
	tree->SetBranchAddress("D0mWeightErr",           &D0mWeightErr);
	tree->SetBranchAddress("D0pWeightErr",            &D0pWeightErr);

	Double_t pi_PX(0);
	Double_t pi_PY(0);
	Double_t pi_PZ(0);
	Double_t K_PX(0);
	Double_t K_PY(0);
	Double_t K_PZ(0);
	Double_t D0m_PX(0);
	Double_t D0m_PY(0);
	Double_t D0m_PZ(0);
	Double_t D0p_PX(0);
	Double_t D0p_PY(0);
	Double_t D0p_PZ(0);

	Double_t pi_P(0);
	Double_t pi_eta(0);
	Double_t K_P(0);
	Double_t K_eta(0);
	Double_t D0m_P(0);
	Double_t D0m_eta(0);
	Double_t D0p_P(0);
	Double_t D0p_eta(0);

	Double_t piCorr(0), KCorr(0.), D0mCorr(0.), D0pCorr(0.);
	Double_t piCorrErr(0), KCorrErr(0.), D0mCorrErr(0.), D0pCorrErr(0.);

	// PID for systematics
	tree->SetBranchAddress("pi_TRUEP_X",        &pi_PX);
	tree->SetBranchAddress("pi_TRUEP_Y",        &pi_PY);
	tree->SetBranchAddress("pi_TRUEP_Z",        &pi_PZ);
	tree->SetBranchAddress("K_TRUEP_X",         &K_PX);
	tree->SetBranchAddress("K_TRUEP_Y",         &K_PY);
	tree->SetBranchAddress("K_TRUEP_Z",         &K_PZ);
	tree->SetBranchAddress("D0m_TRUEP_X",      &D0m_PX);
	tree->SetBranchAddress("D0m_TRUEP_Y",      &D0m_PY);
	tree->SetBranchAddress("D0m_TRUEP_Z",      &D0m_PZ);
	tree->SetBranchAddress("D0p_TRUEP_X",       &D0p_PX);
	tree->SetBranchAddress("D0p_TRUEP_Y",       &D0p_PY);
	tree->SetBranchAddress("D0p_TRUEP_Z",       &D0p_PZ);

	// Get Histograms
	TFile * histfile = TFile::Open("ratio2012S20.root");
	TH2F  * hist     = (TH2F*)histfile->Get("Ratio");
	Int_t nHistBinsX = hist->GetXaxis()->GetNbins();
	Int_t nHistBinsY = hist->GetYaxis()->GetNbins();

	// Variables required for bin matching
	Int_t binpi(-1), binK(-1), binD0m(-1), binD0p(-1);

	Double_t piWeight2(0.), KWeight2(0.), D0mWeight2(0.), D0pWeight2(0.);
	Double_t piWeightErr2(0.), KWeightErr2(0.), D0mWeightErr2(0.), D0pWeightErr2(0.);
	Double_t piCorr2(0), KCorr2(0.), D0mCorr2(0.), D0pCorr2(0.);
	Double_t piCorrErr2(0), KCorrErr2(0.), D0mCorrErr2(0.), D0pCorrErr2(0.);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
		  continue;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0p_TRUEID)==211&&abs(D0m_TRUEID)==211) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0p_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0m_MC_GD_MOTHER_KEY) {
						piCorr   = 0.;
						KCorr    = 0.;
						D0mCorr = 0.;
						D0pCorr  = 0.;
						pi_P       = sqrt(pi_PX * pi_PX + pi_PY * pi_PY + pi_PZ * pi_PZ);
						pi_eta     = 0.5*log( (pi_P + pi_PZ) / (pi_P - pi_PZ) );
						K_P        = sqrt(K_PX  * K_PX  + K_PY  * K_PY  + K_PZ  * K_PZ);
						K_eta      = 0.5*log( (K_P  + K_PZ)  / (K_P  - K_PZ)  );

						D0m_P     = sqrt(D0m_PX * D0m_PX + D0m_PY * D0m_PY + D0m_PZ * D0m_PZ);
						D0m_eta   = 0.5*log( (D0m_P + D0m_PZ) / (D0m_P - D0m_PZ) );
						D0p_P      = sqrt(D0p_PX  * D0p_PX  + D0p_PY  * D0p_PY  + D0p_PZ  * D0p_PZ);
						D0p_eta    = 0.5*log( (D0p_P  + D0p_PZ)  / (D0p_P  - D0p_PZ)  );

						binpi   = hist->FindBin(pi_P/1000.,   pi_eta  );
						binK    = hist->FindBin(K_P/1000.,    K_eta   );
						binD0m = hist->FindBin(D0m_P/1000., D0m_eta);
						binD0p  = hist->FindBin(D0p_P/1000.,  D0p_eta );

						if(binpi   %  (nHistBinsX+2) == 0)            binpi+=1;
						if(binpi   %  (nHistBinsX+2) == nHistBinsX+1) binpi-=1;
						if(binpi   <  (nHistBinsX+2))                 binpi+=(nHistBinsX+2);
						if(binpi   >= (nHistBinsX+2)*(nHistBinsY+1))  binpi-=(nHistBinsX+2);
						if(binK    %  (nHistBinsX+2) == 0)            binK+=1;
						if(binK    %  (nHistBinsX+2) == nHistBinsX+1) binK-=1;
						if(binK    <  (nHistBinsX+2))                 binK+=(nHistBinsX+2);
						if(binK    >= (nHistBinsX+2)*(nHistBinsY+1))  binK-=(nHistBinsX+2);
						if(binD0m %  (nHistBinsX+2) == 0)            binD0m+=1;
						if(binD0m %  (nHistBinsX+2) == nHistBinsX+1) binD0m-=1;
						if(binD0m <  (nHistBinsX+2))                 binD0m+=(nHistBinsX+2);
						if(binD0m >= (nHistBinsX+2)*(nHistBinsY+1))  binD0m-=(nHistBinsX+2);
						if(binD0p  %  (nHistBinsX+2) == 0)            binD0p+=1;
						if(binD0p  %  (nHistBinsX+2) == nHistBinsX+1) binD0p-=1;
						if(binD0p  <  (nHistBinsX+2))                 binD0p+=(nHistBinsX+2);
						if(binD0p  >= (nHistBinsX+2)*(nHistBinsY+1))  binD0p-=(nHistBinsX+2);


						piCorr   = hist->GetBinContent(binpi);
						KCorr    = hist->GetBinContent(binK);
						D0mCorr = hist->GetBinContent(binD0m);
						D0pCorr  = hist->GetBinContent(binD0p);

						piCorrErr   = hist->GetBinError(binpi);
						KCorrErr    = hist->GetBinError(binK);
						D0mCorrErr = hist->GetBinError(binD0m);
						D0pCorrErr  = hist->GetBinError(binD0p);

	                                        piWeight2      = piWeight*piWeight;
						KWeight2       = KWeight*KWeight;
						D0mWeight2    = D0mWeight*D0mWeight;
						D0pWeight2     = D0pWeight*D0pWeight;
	                                        piWeightErr2   = piWeightErr*piWeightErr;
						KWeightErr2    = KWeightErr*KWeightErr;
						D0mWeightErr2 = D0mWeightErr*D0mWeightErr;
						D0pWeightErr2  = D0pWeightErr*D0pWeightErr;
	                                        piCorr2        = piCorr*piCorr;
						KCorr2         = KCorr*KCorr;
						D0mCorr2      = D0mCorr*D0mCorr;
						D0pCorr2       = D0pCorr*D0pCorr;
	                                        piCorrErr2     = piCorrErr*piCorrErr;
						KCorrErr2      = KCorrErr*KCorrErr;
						D0mCorrErr2   = D0mCorrErr*D0mCorrErr;
						D0pCorrErr2    = D0pCorrErr*D0pCorrErr;

	//					if(weight<0.9) cout << pi_P << "\t" << pi_eta << "\t" << K_P << "\t" << K_eta << "\t" << D0m_P << "\t" << D0m_eta << "\t" << D0p_P << "\t" << D0p_eta << endl;   
	//					if(weight<0.9) cout << hist->FindBin(pi_P/1000.,pi_eta) << "\t" << hist->FindBin(K_P/1000.,K_eta) << "\t" 
	//						            << hist->FindBin(D0m_P/1000.,D0m_eta) << "\t" << hist->FindBin(D0p_P/1000.,D0p_eta ) << endl;   
	//					if(weight<0.9) cout << hist->GetBinContent(hist->FindBin(pi_P/1000.,pi_eta)) << "\t" << hist->GetBinContent(hist->FindBin(K_P/1000.,K_eta)) << "\t" 
	//						            << hist->GetBinContent(hist->FindBin(D0m_P/1000.,D0m_eta)) << "\t" << hist->GetBinContent(hist->FindBin(D0p_P/1000.,D0p_eta )) << endl;   

						hpiNoCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight);
						hpiCorr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight   * piCorr);
						hKNoCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight);
						hKCorr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight    * KCorr);
						hD0mNoCorr->Fill( kinematics->getmPrime(), kinematics->getThetaPrime(), D0mWeight);
						hD0mCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), D0mWeight * D0mCorr);
						hD0pNoCorr->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), D0pWeight);
						hD0pCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), D0pWeight  * D0pCorr);

						hpiErr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), /*piWeightErr2   * piCorr2   +*/ piWeight2   * piCorrErr2);
						hKErr->Fill(       kinematics->getmPrime(), kinematics->getThetaPrime(), /*KWeightErr2    * KCorr2    +*/ KWeight2    * KCorrErr2);
						hD0mErr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), /*D0mWeightErr2 * D0mCorr2 +*/ D0mWeight2 * D0mCorrErr2);
						hD0pErr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), /*D0pWeightErr2  * D0pCorr2  +*/ D0pWeight2  * D0pCorrErr2);

						//if(kinematics->getmPrime()>0.92 && kinematics->getThetaPrime()>0.92) {
						//	std::cout << Bd_CM_m23Sq << "\t" << Bd_CM_m12Sq << "\t" << kinematics->getmPrime() << "\t" << kinematics->getThetaPrime() << std::endl;
						//}

						petaPi->Fill(pi_P/1000.,pi_eta);
						petaK->Fill(K_P/1000.,K_eta);
						petaDpi->Fill(D0m_P/1000.,D0m_eta);
						petaDK->Fill(D0p_P/1000.,D0p_eta);
			}
		}	
	}

	file->Close();

//	cout << sdpnocuts->GetEntries() << endl;
//	cout << sdpcut->GetEntries() << endl;
	hpiEff->Divide(  hpiCorr,  hpiNoCorr);
	hKEff->Divide(   hKCorr,   hKNoCorr);
	hD0mEff->Divide(hD0mCorr,hD0mNoCorr);
	hD0pEff->Divide( hD0pCorr, hD0pNoCorr);

	for (Int_t k=0; k<nbins; ++k) {
	   for (Int_t l=0;l<nbins;++l) {
	      Float_t errPi2 = hpiErr->GetBinContent(k+1,l+1);
	      hpiEff->SetBinError(k+1,l+1,sqrt(errPi2)/hpiNoCorr->GetBinContent(k+1,l+1));

	      Float_t errK2 = hKErr->GetBinContent(k+1,l+1);
	      hKEff->SetBinError(k+1,l+1,sqrt(errK2)/hKNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0m2 = hD0mErr->GetBinContent(k+1,l+1);
	      hD0mEff->SetBinError(k+1,l+1,sqrt(errD0m2)/hD0mNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0p2 = hD0pErr->GetBinContent(k+1,l+1);
	      hD0pEff->SetBinError(k+1,l+1,sqrt(errD0p2)/hD0pNoCorr->GetBinContent(k+1,l+1));
	   }
	}

	hAllEff->Add(hpiEff);
	hAllEff->Multiply(hKEff);
	hAllEff->Multiply(hD0mEff);
	hAllEff->Multiply(hD0pEff);

	for (Int_t k=0; k<nbins; ++k) {
	   for (Int_t l=0;l<nbins;++l) {
	      cout << k+1 << "\t" << l+1 << "\t" << hAllEff->GetBinContent(k+1,l+1) << "\t" << hAllEff->GetBinError(k+1,l+1) << endl;
	   }
	}

//	Float_t counter(0);
//	Float_t sum(0);
//
//	for (Int_t k=0;k<200;k++){
//		for (Int_t l=0;l<200;l++) {
//			Float_t temp2(-1);
//			temp2 = sdpeff->GetBinContent(k+1,l+1);
//			if (temp2 >0) {
//				sum = sum + temp2;
//				counter = counter + 1;
//			}
//		}
//	}
//
//	cout << counter << " " << sum << " " << 100*sum/counter << "%" << endl;

	// Add a file to save the plots
	TString fileName("hists/trk_sq_Bd_effs.root");
	outfile = new TFile(fileName,"RECREATE");
	hpiEff->SetName("pi");
	hpiEff->Write();
	hKEff->SetName("K");
	hKEff->Write();
	hD0mEff->SetName("D0m");
	hD0mEff->Write();
	hD0pEff->SetName("D0p");
	hD0pEff->Write();
	hAllEff->SetName("all");
	hAllEff->Write();
	outfile->Close();

	lhcbStyle();
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

	TCanvas * canvas = new TCanvas("c1","c1");
	petaPi->Draw("colz");
	canvas->SaveAs("plots/pi_p_eta.pdf");

	TCanvas * canvas = new TCanvas("c1","c1");
	petaK->Draw("colz");
	canvas->SaveAs("plots/K_p_eta.pdf");

	TCanvas * canvas = new TCanvas("c1","c1");
	petaDpi->Draw("colz");
	canvas->SaveAs("plots/Dpi_p_eta.pdf");

	TCanvas * canvas = new TCanvas("c1","c1");
	petaDK->Draw("colz");
	canvas->SaveAs("plots/DK_p_eta.pdf");

	TCanvas * canvas = new TCanvas("c2","c2");
	canvas->SetRightMargin(0.15);
	hAllEff->SetLabelFont(62,"x");
	hAllEff->SetLabelFont(62,"y");
	hAllEff->SetTitleFont(62,"x");
	hAllEff->SetTitleFont(62,"y");
	hAllEff->SetTitleSize(0.06,"x");
	hAllEff->SetTitleSize(0.06,"y");
	hAllEff->SetLabelSize(0.05,"x");
	hAllEff->SetLabelSize(0.05,"y");
	hAllEff->SetXTitle(xTitle);
	hAllEff->SetYTitle(yTitle);
	hAllEff->GetZaxis()->SetRangeUser(0.95,1.1);
	hAllEff->Draw("colz");
//	printLHCb("R","S");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_eff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hpiEff->SetLabelFont(62,"x");
	hpiEff->SetLabelFont(62,"y");
	hpiEff->SetTitleFont(62,"x");
	hpiEff->SetTitleFont(62,"y");
	hpiEff->SetTitleSize(0.06,"x");
	hpiEff->SetTitleSize(0.06,"y");
	hpiEff->SetLabelSize(0.05,"x");
	hpiEff->SetLabelSize(0.05,"y");
	hpiEff->SetXTitle(xTitle);
	hpiEff->SetYTitle(yTitle);
	hpiEff->GetZaxis()->SetRangeUser(0.98,1.06);
	hpiEff->Draw("colz");
//	printLHCb("R","S");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_pieff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hKEff->SetLabelFont(62,"x");
	hKEff->SetLabelFont(62,"y");
	hKEff->SetTitleFont(62,"x");
	hKEff->SetTitleFont(62,"y");
	hKEff->SetTitleSize(0.06,"x");
	hKEff->SetTitleSize(0.06,"y");
	hKEff->SetLabelSize(0.05,"x");
	hKEff->SetLabelSize(0.05,"y");
	hKEff->SetXTitle(xTitle);
	hKEff->SetYTitle(yTitle);
	hKEff->GetZaxis()->SetRangeUser(0.99,1.03);
	hKEff->Draw("colz");
//	printLHCb("R","S");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_Keff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hD0mEff->SetLabelFont(62,"x");
	hD0mEff->SetLabelFont(62,"y");
	hD0mEff->SetTitleFont(62,"x");
	hD0mEff->SetTitleFont(62,"y");
	hD0mEff->SetTitleSize(0.06,"x");
	hD0mEff->SetTitleSize(0.06,"y");
	hD0mEff->SetLabelSize(0.05,"x");
	hD0mEff->SetLabelSize(0.05,"y");
	hD0mEff->SetXTitle(xTitle);
	hD0mEff->SetYTitle(yTitle);
	hD0mEff->GetZaxis()->SetRangeUser(0.99,1.03);
	hD0mEff->Draw("colz");
//	printLHCb("R","S");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_D0meff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hD0pEff->SetLabelFont(62,"x");
	hD0pEff->SetLabelFont(62,"y");
	hD0pEff->SetTitleFont(62,"x");
	hD0pEff->SetTitleFont(62,"y");
	hD0pEff->SetTitleSize(0.06,"x");
	hD0pEff->SetTitleSize(0.06,"y");
	hD0pEff->SetLabelSize(0.05,"x");
	hD0pEff->SetLabelSize(0.05,"y");
	hD0pEff->SetXTitle(xTitle);
	hD0pEff->SetYTitle(yTitle);
	hD0pEff->GetZaxis()->SetRangeUser(0.99,1.03);
	hD0pEff->Draw("colz");
//	printLHCb("R","S");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_D0peff.pdf";
	canvas->SaveAs(sdpPlotName);
//
//	dpeff->SetLabelFont(62,"x");
//	dpeff->SetLabelFont(62,"y");
//	dpeff->SetTitleFont(62,"x");
//	dpeff->SetTitleFont(62,"y");
//	dpeff->SetTitleSize(0.06,"x");
//	dpeff->SetTitleSize(0.06,"y");
//	dpeff->SetLabelSize(0.05,"x");
//	dpeff->SetLabelSize(0.05,"y");
//	dpeff->SetXTitle("m^{2}(DK) (GeV/c^{2})^{2}");
//	dpeff->SetYTitle("m^{2}(D#pi) (GeV/c^{2})^{2}");
//	dpeff->GetZaxis()->SetRangeUser(0.5,0.9);
//	gStyle->SetPalette(1);
//	canvas->SetRightMargin(0.15);
//	dpeff->Draw("colz");
//	kinematics->drawDPContour();
//	printLHCb("R","S");
//	TString dpPlotName("plots/trk_sq");
//	dpPlotName+=nbins;
//	dpPlotName+="_dp_eff.pdf";
//	canvas->SaveAs(dpPlotName);

	delete kinematics;
}

void trk_eff_D2pipi() {
  gROOT->ProcessLine(".L lhcbStyle.C");
  gSystem->Load( "libEG" );
  gSystem->Load( "libFitEff.so" );
//  gROOT->ProcessLine(".L DPKinematics/LauKinematics.cc");
  lhcbStyle();
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
  trk_Eff();
}
