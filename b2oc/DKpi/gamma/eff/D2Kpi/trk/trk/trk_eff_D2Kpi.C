#include <TLorentzVector.h>

class LauKinematics;

void trk_Eff(){

	TString xTitle="m'";
	TString yTitle="#theta'";
	Int_t nbins(0);

	nbins = 10;

	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	//histograms for numerator, denominator, efficiency and error for each track
	TH2D * hpiCorr     = new TH2D("hpiCorr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiNoCorr   = new TH2D("hpiNoCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiEff      = new TH2D("hpiEff",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hpiErr      = new TH2D("hpiErr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKCorr      = new TH2D("hKCorr",     "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKNoCorr    = new TH2D("hKNoCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKEff       = new TH2D("hKEff",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hKErr       = new TH2D("hKErr",      "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0piCorr   = new TH2D("hD0piCorr",  "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0piNoCorr = new TH2D("hD0piNoCorr","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0piEff    = new TH2D("hD0piEff",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0piErr    = new TH2D("hD0piErr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0KCorr    = new TH2D("hD0KCorr",   "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0KNoCorr  = new TH2D("hD0KNoCorr", "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0KEff     = new TH2D("hD0KEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hD0KErr     = new TH2D("hD0KErr",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	//total efficiency histogram
	TH2D * hAllEff     = new TH2D("hAllEff",    "",nbins,0.0,1.0,nbins,0.0,1.0);

	//p vs eta histograms for each track
	TH2D * petaPi  = new TH2D("petaPi", "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaK   = new TH2D("petaK",  "",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDpi = new TH2D("petaDpi","",50,0.0,100.0,20,1.0,5.0);
	TH2D * petaDK  = new TH2D("petaDK", "",50,0.0,100.0,20,1.0,5.0);
	
	hpiNoCorr->Sumw2();
	hpiCorr->Sumw2();
	hKNoCorr->Sumw2();
	hKCorr->Sumw2();
	hD0piNoCorr->Sumw2();
	hD0piCorr->Sumw2();
	hD0KNoCorr->Sumw2();
	hD0KCorr->Sumw2();

	//open the MC with weights to give same (p, eta, nTracks) distribution as data
	TFile * file = TFile::Open("D2Kpi_Bd2D0Kpi_trackWeight.root");

	TTree * tree = dynamic_cast<TTree*>( file->Get("DecayTree") );
	Int_t nEntries = tree->GetEntries();

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0K_TRUEID,D0pi_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0K_MC_GD_MOTHER_KEY,D0pi_MC_GD_MOTHER_KEY;
	Double_t piWeight(0.), KWeight(0.), D0piWeight(0.), D0KWeight(0.);
	Double_t piWeightErr(0.), KWeightErr(0.), D0piWeightErr(0.), D0KWeightErr(0.);
	
	tree->SetBranchAddress("Bd_CM_m13Sq",             &Bd_CM_m13Sq);
	tree->SetBranchAddress("Bd_CM_m23Sq",             &Bd_CM_m23Sq);
	tree->SetBranchAddress("Bd_CM_m12Sq",             &Bd_CM_m12Sq);
	tree->SetBranchAddress("B_TRUEID",                &B_TRUEID);
	tree->SetBranchAddress("D_TRUEID",                &D_TRUEID);
	tree->SetBranchAddress("D0K_TRUEID",              &D0K_TRUEID);
	tree->SetBranchAddress("D0pi_TRUEID",             &D0pi_TRUEID);
	tree->SetBranchAddress("K_TRUEID",                &K_TRUEID);
	tree->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);

	tree->SetBranchAddress("D0K_MC_GD_MOTHER_KEY",    &D0K_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("D0pi_MC_GD_MOTHER_KEY",   &D0pi_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("K_MC_MOTHER_KEY",         &K_MC_MOTHER_KEY);
	tree->SetBranchAddress("pi_MC_MOTHER_KEY",        &pi_MC_MOTHER_KEY);

	tree->SetBranchAddress("piWeight",                &piWeight);
	tree->SetBranchAddress("KWeight",                 &KWeight);
	tree->SetBranchAddress("D0piWeight",              &D0piWeight);
	tree->SetBranchAddress("D0KWeight",               &D0KWeight);

	tree->SetBranchAddress("piWeightErr",             &piWeightErr);
	tree->SetBranchAddress("KWeightErr",              &KWeightErr);
	tree->SetBranchAddress("D0piWeightErr",           &D0piWeightErr);
	tree->SetBranchAddress("D0KWeightErr",            &D0KWeightErr);

	Double_t pi_PX(0);
	Double_t pi_PY(0);
	Double_t pi_PZ(0);
	Double_t K_PX(0);
	Double_t K_PY(0);
	Double_t K_PZ(0);
	Double_t D0pi_PX(0);
	Double_t D0pi_PY(0);
	Double_t D0pi_PZ(0);
	Double_t D0K_PX(0);
	Double_t D0K_PY(0);
	Double_t D0K_PZ(0);

	Double_t pi_P(0);
	Double_t pi_eta(0);
	Double_t K_P(0);
	Double_t K_eta(0);
	Double_t D0pi_P(0);
	Double_t D0pi_eta(0);
	Double_t D0K_P(0);
	Double_t D0K_eta(0);

	Double_t piCorr(0), KCorr(0.), D0piCorr(0.), D0KCorr(0.);
	Double_t piCorrErr(0), KCorrErr(0.), D0piCorrErr(0.), D0KCorrErr(0.);

	//true 3-momenta of tracks to determine p and eta
	tree->SetBranchAddress("pi_TRUEP_X",        &pi_PX);
	tree->SetBranchAddress("pi_TRUEP_Y",        &pi_PY);
	tree->SetBranchAddress("pi_TRUEP_Z",        &pi_PZ);
	tree->SetBranchAddress("K_TRUEP_X",         &K_PX);
	tree->SetBranchAddress("K_TRUEP_Y",         &K_PY);
	tree->SetBranchAddress("K_TRUEP_Z",         &K_PZ);
	tree->SetBranchAddress("D0pi_TRUEP_X",      &D0pi_PX);
	tree->SetBranchAddress("D0pi_TRUEP_Y",      &D0pi_PY);
	tree->SetBranchAddress("D0pi_TRUEP_Z",      &D0pi_PZ);
	tree->SetBranchAddress("D0K_TRUEP_X",       &D0K_PX);
	tree->SetBranchAddress("D0K_TRUEP_Y",       &D0K_PY);
	tree->SetBranchAddress("D0K_TRUEP_Z",       &D0K_PZ);

	// Get Histogram containing data/MC correction factors in p and eta
	TFile * histfile = TFile::Open("ratio2012S20.root");
	TH2F  * hist     = (TH2F*)histfile->Get("Ratio");
	Int_t nHistBinsX = hist->GetXaxis()->GetNbins();
	Int_t nHistBinsY = hist->GetYaxis()->GetNbins();

	// Variables required for bin matching
	Int_t binpi(-1), binK(-1), binD0pi(-1), binD0K(-1);

	Double_t piWeight2(0.), KWeight2(0.), D0piWeight2(0.), D0KWeight2(0.);
	Double_t piWeightErr2(0.), KWeightErr2(0.), D0piWeightErr2(0.), D0KWeightErr2(0.);
	Double_t piCorr2(0), KCorr2(0.), D0piCorr2(0.), D0KCorr2(0.);
	Double_t piCorrErr2(0), KCorrErr2(0.), D0piCorrErr2(0.), D0KCorrErr2(0.);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		//only use candidates within the DP
		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
		  continue;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		//only use truth-matched candidates
		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0K_TRUEID)==321&&abs(D0pi_TRUEID)==211) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0K_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0pi_MC_GD_MOTHER_KEY) {
						piCorr   = 0.;
						KCorr    = 0.;
						D0piCorr = 0.;
						D0KCorr  = 0.;
						pi_P       = sqrt(pi_PX * pi_PX + pi_PY * pi_PY + pi_PZ * pi_PZ);
						pi_eta     = 0.5*log( (pi_P + pi_PZ) / (pi_P - pi_PZ) );
						K_P        = sqrt(K_PX  * K_PX  + K_PY  * K_PY  + K_PZ  * K_PZ);
						K_eta      = 0.5*log( (K_P  + K_PZ)  / (K_P  - K_PZ)  );

						D0pi_P     = sqrt(D0pi_PX * D0pi_PX + D0pi_PY * D0pi_PY + D0pi_PZ * D0pi_PZ);
						D0pi_eta   = 0.5*log( (D0pi_P + D0pi_PZ) / (D0pi_P - D0pi_PZ) );
						D0K_P      = sqrt(D0K_PX  * D0K_PX  + D0K_PY  * D0K_PY  + D0K_PZ  * D0K_PZ);
						D0K_eta    = 0.5*log( (D0K_P  + D0K_PZ)  / (D0K_P  - D0K_PZ)  );

						binpi   = hist->FindBin(pi_P/1000.,   pi_eta  );
						binK    = hist->FindBin(K_P/1000.,    K_eta   );
						binD0pi = hist->FindBin(D0pi_P/1000., D0pi_eta);
						binD0K  = hist->FindBin(D0K_P/1000.,  D0K_eta );

						//ratios not defined over all values of (p,eta)
						//move over/underflow candidates into the nearest bin
						if(binpi   %  (nHistBinsX+2) == 0)            binpi+=1;
						if(binpi   %  (nHistBinsX+2) == nHistBinsX+1) binpi-=1;
						if(binpi   <  (nHistBinsX+2))                 binpi+=(nHistBinsX+2);
						if(binpi   >= (nHistBinsX+2)*(nHistBinsY+1))  binpi-=(nHistBinsX+2);
						if(binK    %  (nHistBinsX+2) == 0)            binK+=1;
						if(binK    %  (nHistBinsX+2) == nHistBinsX+1) binK-=1;
						if(binK    <  (nHistBinsX+2))                 binK+=(nHistBinsX+2);
						if(binK    >= (nHistBinsX+2)*(nHistBinsY+1))  binK-=(nHistBinsX+2);
						if(binD0pi %  (nHistBinsX+2) == 0)            binD0pi+=1;
						if(binD0pi %  (nHistBinsX+2) == nHistBinsX+1) binD0pi-=1;
						if(binD0pi <  (nHistBinsX+2))                 binD0pi+=(nHistBinsX+2);
						if(binD0pi >= (nHistBinsX+2)*(nHistBinsY+1))  binD0pi-=(nHistBinsX+2);
						if(binD0K  %  (nHistBinsX+2) == 0)            binD0K+=1;
						if(binD0K  %  (nHistBinsX+2) == nHistBinsX+1) binD0K-=1;
						if(binD0K  <  (nHistBinsX+2))                 binD0K+=(nHistBinsX+2);
						if(binD0K  >= (nHistBinsX+2)*(nHistBinsY+1))  binD0K-=(nHistBinsX+2);


						piCorr   = hist->GetBinContent(binpi);
						KCorr    = hist->GetBinContent(binK);
						D0piCorr = hist->GetBinContent(binD0pi);
						D0KCorr  = hist->GetBinContent(binD0K);

						piCorrErr   = hist->GetBinError(binpi);
						KCorrErr    = hist->GetBinError(binK);
						D0piCorrErr = hist->GetBinError(binD0pi);
						D0KCorrErr  = hist->GetBinError(binD0K);

	                                        piWeight2      = piWeight*piWeight;
						KWeight2       = KWeight*KWeight;
						D0piWeight2    = D0piWeight*D0piWeight;
						D0KWeight2     = D0KWeight*D0KWeight;
	                                        piWeightErr2   = piWeightErr*piWeightErr;
						KWeightErr2    = KWeightErr*KWeightErr;
						D0piWeightErr2 = D0piWeightErr*D0piWeightErr;
						D0KWeightErr2  = D0KWeightErr*D0KWeightErr;
	                                        piCorr2        = piCorr*piCorr;
						KCorr2         = KCorr*KCorr;
						D0piCorr2      = D0piCorr*D0piCorr;
						D0KCorr2       = D0KCorr*D0KCorr;
	                                        piCorrErr2     = piCorrErr*piCorrErr;
						KCorrErr2      = KCorrErr*KCorrErr;
						D0piCorrErr2   = D0piCorrErr*D0piCorrErr;
						D0KCorrErr2    = D0KCorrErr*D0KCorrErr;


						//fill numerators and denominators for correction due to each track weighted according to the (p,eta,nTracks) weight of the relevant track
						hpiNoCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight);
						hpiNoCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight);
						hpiCorr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight   * piCorr);
						hKNoCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight);
						hKCorr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight    * KCorr);
						hD0piNoCorr->Fill( kinematics->getmPrime(), kinematics->getThetaPrime(), D0piWeight);
						hD0piCorr->Fill(   kinematics->getmPrime(), kinematics->getThetaPrime(), D0piWeight * D0piCorr);
						hD0KNoCorr->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), D0KWeight);
						hD0KCorr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), D0KWeight  * D0KCorr);

						hpiErr->Fill(      kinematics->getmPrime(), kinematics->getThetaPrime(), piWeight2   * piCorrErr2);
						hKErr->Fill(       kinematics->getmPrime(), kinematics->getThetaPrime(), KWeight2    * KCorrErr2);
						hD0piErr->Fill(    kinematics->getmPrime(), kinematics->getThetaPrime(), D0piWeight2 * D0piCorrErr2);
						hD0KErr->Fill(     kinematics->getmPrime(), kinematics->getThetaPrime(), D0KWeight2  * D0KCorrErr2);

						petaPi->Fill(pi_P/1000.,pi_eta);
						petaK->Fill(K_P/1000.,K_eta);
						petaDpi->Fill(D0pi_P/1000.,D0pi_eta);
						petaDK->Fill(D0K_P/1000.,D0K_eta);
			}
		}	
	}

	file->Close();

	//produce efficiencies
	hpiEff->Divide(  hpiCorr,  hpiNoCorr);
	hKEff->Divide(   hKCorr,   hKNoCorr);
	hD0piEff->Divide(hD0piCorr,hD0piNoCorr);
	hD0KEff->Divide( hD0KCorr, hD0KNoCorr);

	for (Int_t k=0; k<nbins; ++k) {
	   for (Int_t l=0;l<nbins;++l) {
	      Float_t errPi2 = hpiErr->GetBinContent(k+1,l+1);
	      hpiEff->SetBinError(k+1,l+1,sqrt(errPi2)/hpiNoCorr->GetBinContent(k+1,l+1));

	      Float_t errK2 = hKErr->GetBinContent(k+1,l+1);
	      hKEff->SetBinError(k+1,l+1,sqrt(errK2)/hKNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0pi2 = hD0piErr->GetBinContent(k+1,l+1);
	      hD0piEff->SetBinError(k+1,l+1,sqrt(errD0pi2)/hD0piNoCorr->GetBinContent(k+1,l+1));

	      Float_t errD0K2 = hD0KErr->GetBinContent(k+1,l+1);
	      hD0KEff->SetBinError(k+1,l+1,sqrt(errD0K2)/hD0KNoCorr->GetBinContent(k+1,l+1));
	   }
	}

	hAllEff->Add(hpiEff);
	hAllEff->Multiply(hKEff);
	hAllEff->Multiply(hD0piEff);
	hAllEff->Multiply(hD0KEff);

	for (Int_t k=0; k<nbins; ++k) {
	   for (Int_t l=0;l<nbins;++l) {
	      cout << k+1 << "\t" << l+1 << "\t" << hAllEff->GetBinContent(k+1,l+1) << "\t" << hAllEff->GetBinError(k+1,l+1) << endl;
	   }
	}

	// Add a file to save the histograms
	TString fileName("hists/trk_sq_Bd_effs.root");
	outfile = new TFile(fileName,"RECREATE");
	hpiEff->SetName("pi");
	hpiEff->Write();
	hKEff->SetName("K");
	hKEff->Write();
	hD0piEff->SetName("D0pi");
	hD0piEff->Write();
	hD0KEff->SetName("D0K");
	hD0KEff->Write();
	hAllEff->SetName("all");
	hAllEff->Write();
	outfile->Close();

	//make plots
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
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_Keff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hD0piEff->SetLabelFont(62,"x");
	hD0piEff->SetLabelFont(62,"y");
	hD0piEff->SetTitleFont(62,"x");
	hD0piEff->SetTitleFont(62,"y");
	hD0piEff->SetTitleSize(0.06,"x");
	hD0piEff->SetTitleSize(0.06,"y");
	hD0piEff->SetLabelSize(0.05,"x");
	hD0piEff->SetLabelSize(0.05,"y");
	hD0piEff->SetXTitle(xTitle);
	hD0piEff->SetYTitle(yTitle);
	hD0piEff->GetZaxis()->SetRangeUser(0.99,1.03);
	hD0piEff->Draw("colz");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_D0pieff.pdf";
	canvas->SaveAs(sdpPlotName);

	canvas->SetRightMargin(0.15);
	hD0KEff->SetLabelFont(62,"x");
	hD0KEff->SetLabelFont(62,"y");
	hD0KEff->SetTitleFont(62,"x");
	hD0KEff->SetTitleFont(62,"y");
	hD0KEff->SetTitleSize(0.06,"x");
	hD0KEff->SetTitleSize(0.06,"y");
	hD0KEff->SetLabelSize(0.05,"x");
	hD0KEff->SetLabelSize(0.05,"y");
	hD0KEff->SetXTitle(xTitle);
	hD0KEff->SetYTitle(yTitle);
	hD0KEff->GetZaxis()->SetRangeUser(0.99,1.03);
	hD0KEff->Draw("colz");
	TString sdpPlotName("plots/trk_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_D0Keff.pdf";
	canvas->SaveAs(sdpPlotName);

	delete kinematics;
}

void trk_eff_D2Kpi() {
  gROOT->ProcessLine(".L lhcbStyle.C");
  gSystem->Load( "libEG" );
  gSystem->Load( "libFitEff.so" );
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
