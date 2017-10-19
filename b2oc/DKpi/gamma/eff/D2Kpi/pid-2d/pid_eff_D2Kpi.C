#include <TLorentzVector.h>

class LauKinematics;

void pid_Eff(TString jobName, TString fileName, TString treeName, Double_t NNmin=-0.8, Double_t NNmax=1., TString binLabel="",
             TString xTitle="m'", TString yTitle="#theta'"){

	Int_t nbins(0);

	nbins = 10;

	Int_t NNmin1k = NNmin*1000;
	Int_t NNmax1k = NNmax*1000;

	//double mB(5.27958), mBs(5.3662998046875), mD(1.86484), mK(0.493677001953125000), mPi(0.139570175170830879);
	double mB(5.27958), Bs(5.3663), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * sdpcut = new TH2D("Dpipi CUT SDP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdpnocuts = new TH2D("Dpipi NO CUTS SDP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdpeff = new TH2D("sdpeff","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdperr2 = new TH2D("sdperr2","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdperr = new TH2D("sdperr","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdperrpc = new TH2D("sdperrpc","",nbins,0.0,1.0,nbins,0.0,1.0);
//	TH2D * dpcut = new TH2D("Dpipi CUT DP","",nbins,5.0,28.0,nbins,3.0,25.0);
//	TH2D * dpnocuts = new TH2D("Dpipi NO CUTS DP","",nbins,5.0,28.0,nbins,3.0,25.0);
//	TH2D * dpeff = new TH2D("dpeff","",nbins,5.0,28.0,nbins,3.0,25.0);

	sdpnocuts->Sumw2();
	sdpcut->Sumw2();
//	dpnocuts->Sumw2();
//	dpcut->Sumw2();

	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job"+jobName+"/"+fileName+".root");

	TTree * tree = dynamic_cast<TTree*>( file->Get(treeName) );
	Int_t nEntries = tree->GetEntries();

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0K_TRUEID,D0pi_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0K_MC_GD_MOTHER_KEY,D0pi_MC_GD_MOTHER_KEY;
	Float_t NN;
	Double_t pythia;
	
	tree->SetBranchAddress("Bd_CM_m13Sq",           &Bd_CM_m13Sq);
	tree->SetBranchAddress("Bd_CM_m23Sq",           &Bd_CM_m23Sq);
	tree->SetBranchAddress("Bd_CM_m12Sq",           &Bd_CM_m12Sq);
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

	tree->SetBranchAddress("NN",                      &NN);

	Double_t D0_p_CM_PX(0);
	Double_t D0_p_CM_PY(0);
	Double_t D0_p_CM_PZ(0);
	Double_t D0_m_CM_PX(0);
	Double_t D0_m_CM_PY(0);
	Double_t D0_m_CM_PZ(0);

	Double_t D0_p_CM_P(0);
	Double_t D0_p_CM_PT(0);
	Double_t D0_m_CM_P(0);
	Double_t D0_m_CM_PT(0);

	Double_t D0_D0p_CM_PX(0);
	Double_t D0_D0p_CM_PY(0);
	Double_t D0_D0p_CM_PZ(0);
	Double_t D0_D0m_CM_PX(0);
	Double_t D0_D0m_CM_PY(0);
	Double_t D0_D0m_CM_PZ(0);

	Double_t D0_D0p_CM_P(0);
	Double_t D0_D0p_CM_PT(0);
	Double_t D0_D0m_CM_P(0);
	Double_t D0_D0m_CM_PT(0);

	Int_t NOCM_K_ID(0);

	Double_t B_D0_B_CM_PX(0.), B_D0_B_CM_PY(0.), B_D0_B_CM_PT(0.);
	
	Bool_t B_L0Global_TIS;

	Double_t weight(0);
	Double_t err2(0);

	Double_t piEff(0);
	Double_t piErr(0);
	Double_t KEff(0);
	Double_t KErr(0);
	Double_t DpiEff(0);
	Double_t DpiErr(0);
	Double_t DKEff(0);
	Double_t DKErr(0);

	// PID for systematics
	tree->SetBranchAddress("B_D0_p_CM_PX",      &D0_p_CM_PX);
	tree->SetBranchAddress("B_D0_p_CM_PY",      &D0_p_CM_PY);
	tree->SetBranchAddress("B_D0_p_CM_PZ",      &D0_p_CM_PZ);
	tree->SetBranchAddress("B_D0_m_CM_PX",      &D0_m_CM_PX);
	tree->SetBranchAddress("B_D0_m_CM_PY",      &D0_m_CM_PY);
	tree->SetBranchAddress("B_D0_m_CM_PZ",      &D0_m_CM_PZ);

	tree->SetBranchAddress("B_D0_D0p_CM_PX",    &D0_D0p_CM_PX);
	tree->SetBranchAddress("B_D0_D0p_CM_PY",    &D0_D0p_CM_PY);
	tree->SetBranchAddress("B_D0_D0p_CM_PZ",    &D0_D0p_CM_PZ);
	tree->SetBranchAddress("B_D0_D0m_CM_PX",    &D0_D0m_CM_PX);
	tree->SetBranchAddress("B_D0_D0m_CM_PY",    &D0_D0m_CM_PY);
	tree->SetBranchAddress("B_D0_D0m_CM_PZ",    &D0_D0m_CM_PZ);

	tree->SetBranchAddress("K_ID",              &NOCM_K_ID);

	tree->SetBranchAddress("B_D0_B_CM_PX",      &B_D0_B_CM_PX);
	tree->SetBranchAddress("B_D0_B_CM_PY",      &B_D0_B_CM_PY);

	tree->SetBranchAddress("B_L0Global_TIS",  &B_L0Global_TIS);

	// Get Histograms
	TFile * histfile    = TFile::Open("/home/phrkbf/PIDCalibHists/new/PerfHists_K_Strip20rX_MagBoth_DKpi_PPT4a_P_PT.root");
	TH2F* K_hist        = (TH2F*)histfile->Get("K_ProbNNK * (1 - ProbNNpi ) > 0.3_All");
	TH2F* DK_hist       = (TH2F*)histfile->Get("K_ProbNNK * (1 - ProbNNpi ) > 0.1_All");

	TFile * histfile2   = TFile::Open("/home/phrkbf/PIDCalibHists/new/PerfHists_Pi_Strip20rX_MagBoth_DKpi_PPT4a_P_PT.root");
	TH2F* pi_hist       = (TH2F*)histfile2->Get("Pi_ProbNNpi * (1 - ProbNNK ) > 0.2_All");
	TH2F* Dpi_hist      = (TH2F*)histfile2->Get("Pi_ProbNNpi * (1 - ProbNNK ) > 0.1_All");

	// Variables required for bin matching
	Int_t binPi(-1), binK(-1), binDpi(-1), binDK(-1);

	for ( Int_t i(0); i < nEntries; ++i ) {
		if(i%10000==0) std::cout << i << " of " << nEntries << std::endl;

		tree->GetEntry( i );

		if(NN<NNmin || NN>NNmax) continue;

		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
		  continue;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0K_TRUEID)==321&&abs(D0pi_TRUEID)==211) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0K_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0pi_MC_GD_MOTHER_KEY) {
				weight=0;
				err2=0;
				piEff=0;
				piErr=0;
				KEff=0;
				KErr=0;
	                        DpiEff=0;
	                        DpiErr=0;
	                        DKEff=0;
	                        DKErr=0;


				D0_p_CM_P  = sqrt(D0_p_CM_PX*D0_p_CM_PX + D0_p_CM_PY*D0_p_CM_PY + D0_p_CM_PZ*D0_p_CM_PZ);
				D0_p_CM_PT = sqrt(D0_p_CM_PX*D0_p_CM_PX + D0_p_CM_PY*D0_p_CM_PY);
				D0_m_CM_P  = sqrt(D0_m_CM_PX*D0_m_CM_PX + D0_m_CM_PY*D0_m_CM_PY + D0_m_CM_PZ*D0_m_CM_PZ);
				D0_m_CM_PT = sqrt(D0_m_CM_PX*D0_m_CM_PX + D0_m_CM_PY*D0_m_CM_PY);

				D0_D0p_CM_P  = sqrt(D0_D0p_CM_PX*D0_D0p_CM_PX + D0_D0p_CM_PY*D0_D0p_CM_PY + D0_D0p_CM_PZ*D0_D0p_CM_PZ);
				D0_D0p_CM_PT = sqrt(D0_D0p_CM_PX*D0_D0p_CM_PX + D0_D0p_CM_PY*D0_D0p_CM_PY);
				D0_D0m_CM_P  = sqrt(D0_D0m_CM_PX*D0_D0m_CM_PX + D0_D0m_CM_PY*D0_D0m_CM_PY + D0_D0m_CM_PZ*D0_D0m_CM_PZ);
				D0_D0m_CM_PT = sqrt(D0_D0m_CM_PX*D0_D0m_CM_PX + D0_D0m_CM_PY*D0_D0m_CM_PY);

					// Get eff for K track
				if(NOCM_K_ID>0) {
				   binPi  = pi_hist->FindBin(D0_m_CM_P,D0_m_CM_PT);
				   binK   = K_hist->FindBin(D0_p_CM_P,D0_p_CM_PT);
				   binDpi = Dpi_hist->FindBin(D0_D0m_CM_P,D0_D0m_CM_PT);
				   binDK  = DK_hist->FindBin(D0_D0p_CM_P,D0_D0p_CM_PT);

				} else {
				   binPi  = pi_hist->FindBin(D0_p_CM_P,D0_p_CM_PT);
				   binK   = K_hist->FindBin(D0_m_CM_P,D0_m_CM_PT);
				   binDpi = Dpi_hist->FindBin(D0_D0p_CM_P,D0_D0p_CM_PT);
				   binDK  = DK_hist->FindBin(D0_D0m_CM_P,D0_D0m_CM_PT);
				}

				piEff  = pi_hist->GetBinContent(binPi);
				piErr  = pi_hist->GetBinError(binPi);

				KEff   = K_hist->GetBinContent(binK);
				KErr   = K_hist->GetBinError(binK);

				DpiEff = Dpi_hist->GetBinContent(binDpi);
				DpiErr = Dpi_hist->GetBinError(binDpi);

				DKEff  = DK_hist->GetBinContent(binDK);
				DKErr  = DK_hist->GetBinError(binDK);


				weight = piEff * KEff * DpiEff * DKEff;
				err2   = piErr*piErr*KEff*KEff*DpiEff*DpiEff*DKEff*DKEff +
					 piEff*piEff*KErr*KErr*DpiEff*DpiEff*DKEff*DKEff +
					 piEff*piEff*KEff*KEff*DpiErr*DpiErr*DKEff*DKEff +
					 piEff*piEff*KEff*KEff*DpiEff*DpiEff*DKErr*DKErr;
				//std::cout << weight << "\t" << err2 << piEff << "\t" << KEff << "\t" << DpiEff << "\t" << DKEff << std::endl;

				sdpnocuts->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
				sdpcut->Fill(kinematics->getmPrime(),kinematics->getThetaPrime(),weight);
				sdperr2->Fill(kinematics->getmPrime(),kinematics->getThetaPrime(),err2);
//				dpnocuts->Fill(Bd_CM_m23Sq,Bd_CM_m12Sq);
//				dpcut->Fill(Bd_CM_m23Sq,Bd_CM_m12Sq,weight);

			}
		}	
	}

	file->Close();

	cout << sdpnocuts->GetEntries() << endl;
	cout << sdpcut->GetEntries() << endl;
	sdpeff->Divide(sdpcut,sdpnocuts);
	sdperr2->Divide(sdpnocuts);
//	dpeff->Divide(dpcut,dpnocuts);

	for(Int_t j=1; j<=nbins; ++j) {
		for(Int_t i=1; i<=nbins; ++i) {
			sdpeff->SetBinError(i,j,sqrt(sdperr2->GetBinContent(i,j)));
			sdperr->SetBinContent(i,j,sqrt(sdperr2->GetBinContent(i,j)));
			sdperrpc->SetBinContent(i,j,100.*sqrt(sdperr2->GetBinContent(i,j))/sdpeff->GetBinContent(i,j));
		}
	}

	Float_t counter(0);
	Float_t sum(0);

	for (Int_t k=0;k<200;k++){
		for (Int_t l=0;l<200;l++) {
			Float_t temp2(-1);
			temp2 = sdpeff->GetBinContent(k+1,l+1);
			if (temp2 >0) {
				sum = sum + temp2;
				counter = counter + 1;
			}
		}
	}

	cout << "Ignored " << ignored << " of " << total << " events." << endl;
	cout << counter << " " << sum << " " << 100*sum/counter << "%" << endl;

	// Add a file to save the plots
	TString fileName("hists/pid");
	//fileName+=nbins;
	//fileName+="_"; fileName+=NNmin1k; fileName+="_"; fileName+=NNmax1k;
	fileName += "_"; fileName += binLabel;
	fileName+=".root";
	outfile = new TFile(fileName,"RECREATE");
	sdpeff->SetName("efficiency");
	sdpeff->Write();
	outfile->Close();

	TString fileName("hists/pid");
	//fileName+=nbins;
	//fileName+="_"; fileName+=NNmin1k; fileName+="_"; fileName+=NNmax1k;
	fileName += "_"; fileName += binLabel;
	fileName+="_errPid.root";
	outfile = new TFile(fileName,"RECREATE");
	sdperr->SetName("error");
	sdperr->Write();
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

//	TCanvas * canvas = new TCanvas("c0","c0");
//	dpnocuts->Draw("colz");
//	TCanvas * canvas = new TCanvas("c1","c1");
//      dpcut->Draw("colz");
	TCanvas * canvas = new TCanvas("c2","c2");
	sdpeff->SetLabelFont(62,"x");
	sdpeff->SetLabelFont(62,"y");
	sdpeff->SetTitleFont(62,"x");
	sdpeff->SetTitleFont(62,"y");
	sdpeff->SetTitleSize(0.06,"x");
	sdpeff->SetTitleSize(0.06,"y");
	sdpeff->SetLabelSize(0.05,"x");
	sdpeff->SetLabelSize(0.05,"y");
	sdpeff->SetXTitle(xTitle);
	sdpeff->SetYTitle(yTitle);
	sdpeff->GetZaxis()->SetRangeUser(0.4,0.8);
	canvas->SetRightMargin(0.15);
	sdpeff->Draw("colztext");
	//kinematics->drawDPContour();
	//printLHCb("R","S");
	TString sdpPlotName("plots/pid_sq");
	sdpPlotName+=nbins;
//	sdpPlotName+="_"; sdpPlotName+=NNmin1k; sdpPlotName+="_"; sdpPlotName+=NNmax1k;
	sdpPlotName += "_"; sdpPlotName += binLabel;
	sdpPlotName+="_eff.pdf";
	canvas->SaveAs(sdpPlotName);

	gStyle->SetPaintTextFormat("0.2f %%");

	sdperrpc->SetLabelFont(62,"x");
	sdperrpc->SetLabelFont(62,"y");
	sdperrpc->SetTitleFont(62,"x");
	sdperrpc->SetTitleFont(62,"y");
	sdperrpc->SetTitleSize(0.06,"x");
	sdperrpc->SetTitleSize(0.06,"y");
	sdperrpc->SetLabelSize(0.05,"x");
	sdperrpc->SetLabelSize(0.05,"y");
	sdperrpc->SetXTitle(xTitle);
	sdperrpc->SetYTitle(yTitle);
//	sdperrpc->GetZaxis()->SetRangeUser(0.5,0.9);
//	canvas->SetRightMargin(0.15);
	sdperrpc->Draw("colztext");
	TString sdpPlotName("plots/pid_sq");
	sdpPlotName+=nbins;
//	sdpPlotName+="_"; sdpPlotName+=NNmin1k; sdpPlotName+="_"; sdpPlotName+=NNmax1k;
	sdpPlotName += "_"; sdpPlotName += binLabel;
	sdpPlotName+="_err.pdf";
	canvas->SaveAs(sdpPlotName);

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
////	if(fileName2=="DKpi_Bs2DKpi_selBs_noNN_noPID_noTrig") eff3->GetZaxis()->SetRangeUser(0,0.02*1.33);
//	gStyle->SetPalette(1);
//	canvas->SetRightMargin(0.15);
//	dpeff->Draw("colz");
//	kinematics->drawDPContour();
//	printLHCb("R","S");
//	TString dpPlotName("plots/pid_sq");
//	dpPlotName+=nbins;
//	dpPlotName+="_dp_eff.pdf";
//	canvas->SaveAs(dpPlotName);

	delete kinematics;
}

void pid_eff_D2Kpi() {
  gROOT->ProcessLine(".L lhcbStyle.C");
  gSystem->Load( "libEG" );
  gSystem->Load( "libFitEff.so" );
  lhcbStyle();
  gStyle->SetOptStat(0000);
//  pid_Eff("Sim08aMC","Bd2DKpi/Bd2D0Kpi/Bd2DKpi_Bd2D0Kpi_P8_selBd_noDPVetoes_NN","DecayTree");
  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",-0.80, 1.00,"all");
//  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",-0.80, 0.00,"a");
//  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", 0.00, 0.50,"b");
//  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", 0.50, 0.72,"c");
//  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", 0.72, 0.81,"d");
//  pid_Eff("Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", 0.81, 1.00,"e");

}
