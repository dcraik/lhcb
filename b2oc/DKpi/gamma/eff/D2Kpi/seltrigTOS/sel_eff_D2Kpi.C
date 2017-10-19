#include <TLorentzVector.h>

class LauKinematics;

//TH2D* correction(const TH2D* numerator, const TH2D* denomenator, Double_t tolerance);

void sel_Eff(TString jobName, TString fileName, TString treeName, TString jobName2, TString fileName2, TString treeName2, Double_t NNmin=-0.8, Double_t NNmax=1.0, TString binLabel="",
             TString xTitle="m'", TString yTitle="#theta'"){

	Int_t nbins(0);

	nbins = 10;

	Int_t NNmin1k = NNmin * 1000;
	Int_t NNmax1k = NNmax * 1000;

	//double mB(5.27958), mBs(5.3662998046875), mD(1.86484), mK(0.493677001953125000), mPi(0.139570175170830879);
	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * dpcut = new TH2D("Dpipi CUT DP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * dpnocuts = new TH2D("Dpipi NO CUTS DP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * eff = new TH2D("eff","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * errp = new TH2D("errp","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * errm = new TH2D("errm","",nbins,0.0,1.0,nbins,0.0,1.0);
//	TH2D * eff2;
//	TH2D * eff3;
//	TH2D * eff4;

	//TFile * file = TFile::Open("Dpipi_MC11_Truth.root");
	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job"+jobName+"/"+fileName+".root");

	//TTree * tree = dynamic_cast<TTree*>( file->Get("Bd2pipiD0/tupleMCTruth") );
	TTree * tree = dynamic_cast<TTree*>( file->Get(treeName) );
	Int_t nEntries = tree->GetEntries();

//	Double_t m13Sq_MC(-1.0);
//	Double_t m12Sq_MC(-1.0);

//	tree->SetBranchAddress("m23Sq_MC",&m23Sq_MC);
//	tree->SetBranchAddress("m12Sq_MC",&m12Sq_MC);

	Double_t D0K_TRUETHETA(-1.0);
	Double_t D0pi_TRUETHETA(-1.0);
	Double_t h1_TRUETHETA(-1.0);
	Double_t h2_TRUETHETA(-1.0);

	Double_t h1_CORRPE(-1.0);
	Double_t h1_CORRPX(-1.0);
	Double_t h1_CORRPY(-1.0);
	Double_t h1_CORRPZ(-1.0);

	Double_t h2_CORRPE(-1.0);
	Double_t h2_CORRPX(-1.0);
	Double_t h2_CORRPY(-1.0);
	Double_t h2_CORRPZ(-1.0);

	Double_t D0_CORRPE(-1.0);
	Double_t D0_CORRPX(-1.0);
	Double_t D0_CORRPY(-1.0);
	Double_t D0_CORRPZ(-1.0);

//	Double_t D_TRUEID(0);

	TLorentzVector D0;
	TLorentzVector h1;
	TLorentzVector h2;
	TLorentzVector D0h1;
	TLorentzVector D0h2;

	tree->SetBranchAddress("D0K_TRUETHETA",&D0K_TRUETHETA);
	tree->SetBranchAddress("D0pi_TRUETHETA",&D0pi_TRUETHETA);
	tree->SetBranchAddress("h1_TRUETHETA",&h1_TRUETHETA);
	tree->SetBranchAddress("h2_TRUETHETA",&h2_TRUETHETA);

	tree->SetBranchAddress("D0_CORRPX",&D0_CORRPX);
	tree->SetBranchAddress("D0_CORRPY",&D0_CORRPY);
	tree->SetBranchAddress("D0_CORRPZ",&D0_CORRPZ);
	tree->SetBranchAddress("D0_CORRPE",&D0_CORRPE);

	tree->SetBranchAddress("h1_CORRPX",&h1_CORRPX);
	tree->SetBranchAddress("h1_CORRPY",&h1_CORRPY);
	tree->SetBranchAddress("h1_CORRPZ",&h1_CORRPZ);
	tree->SetBranchAddress("h1_CORRPE",&h1_CORRPE);

	tree->SetBranchAddress("h2_CORRPX",&h2_CORRPX);
	tree->SetBranchAddress("h2_CORRPY",&h2_CORRPY);
	tree->SetBranchAddress("h2_CORRPZ",&h2_CORRPZ);
	tree->SetBranchAddress("h2_CORRPE",&h2_CORRPE);

//	tree->SetBranchAddress("D_TRUEID",&D_TRUEID);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		h1.SetPxPyPzE(h1_CORRPX, h1_CORRPY, h1_CORRPZ, h1_CORRPE);
		h2.SetPxPyPzE(h2_CORRPX, h2_CORRPY, h2_CORRPZ, h2_CORRPE);
		D0.SetPxPyPzE(D0_CORRPX, D0_CORRPY, D0_CORRPZ, D0_CORRPE);

		D0h1 = D0 + h1;
		D0h2 = D0 + h2;
		h1h2 = h1 + h2;
		//cout << (D0 + h1 + h2).M() << "\t" << D0.M() << "\t" << h1.M() << "\t" << h2.M() << endl;

		++total;
		if(!kinematics->withinDPLimits(D0h1.M2()/1000000,h1h2.M2()/1000000)) {
//                  cout << "event ignored at " << D0h1.M2() << ", " << D0h2.M2() << endl;
                  ++ignored;
		  continue;
//		} else {
//                  cout << "event accepted at " << D0h1.M2() << ", " << D0h2.M2() << endl;
		}
                kinematics->updateKinematics(D0h1.M2()/1000000,h1h2.M2()/1000000);

//		if(jobName=="2011_Bd2Dpipi_eff" && D_TRUEID<0 ) {
//		  dpnocuts->Fill(D0h2.M2()/1000000,D0h1.M2()/1000000);
//		} else {
		dpnocuts->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
//		}
	}

	file->Close();


	// Now create and fill the other Histogram!
	//TFile * file2 = TFile::Open("Dpipi_MC11_SelEff_VETO.root");
	TFile * file2 = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job"+jobName2+"/"+fileName2+".root");

	TTree * tree2 = dynamic_cast<TTree*>( file2->Get(treeName2) );
	Int_t nEntries2 = tree2->GetEntries();
//	cout<<nEntries2<<endl;

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
//	Int_t Matched;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0K_TRUEID,D0pi_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0K_MC_GD_MOTHER_KEY,D0pi_MC_GD_MOTHER_KEY;
	Double_t D0K_ProbNNk, D0K_ProbNNpi, D0pi_ProbNNk, D0pi_ProbNNpi;
	Bool_t B_L0HadronDecision_TOS;
	Float_t NN;
//	Int_t hlt2global(-1);
//	Int_t l0decision(-1);
//	Int_t l0hadron(-1);
//	Int_t hlt2topo2TOS(-1);
//	Int_t hlt2topo3TOS(-1);
//	Int_t hlt2topo4TOS(-1);
	//here 123 are KpiD
//	if(bs) {
//          tree2->SetBranchAddress("Bs_CM_m13Sq",&Bd_CM_m13Sq);
//          tree2->SetBranchAddress("Bs_CM_m23Sq",&Bd_CM_m23Sq);
//	} else {
	  tree2->SetBranchAddress("Bd_CM_m13Sq",           &Bd_CM_m13Sq); //DK
	  tree2->SetBranchAddress("Bd_CM_m23Sq",           &Bd_CM_m23Sq); //Dpi
	  tree2->SetBranchAddress("Bd_CM_m12Sq",           &Bd_CM_m12Sq); //Kpi
//	}
//	tree2->SetBranchAddress("Matched",&Matched);
	tree2->SetBranchAddress("B_TRUEID",                &B_TRUEID);
	tree2->SetBranchAddress("D_TRUEID",                &D_TRUEID);
	tree2->SetBranchAddress("D0K_TRUEID",              &D0K_TRUEID);
	tree2->SetBranchAddress("D0pi_TRUEID",             &D0pi_TRUEID);
	tree2->SetBranchAddress("K_TRUEID",                &K_TRUEID);
	tree2->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);

	tree2->SetBranchAddress("D0K_MC_GD_MOTHER_KEY",    &D0K_MC_GD_MOTHER_KEY);
	tree2->SetBranchAddress("D0pi_MC_GD_MOTHER_KEY",   &D0pi_MC_GD_MOTHER_KEY);
	tree2->SetBranchAddress("K_MC_MOTHER_KEY",         &K_MC_MOTHER_KEY);
	tree2->SetBranchAddress("pi_MC_MOTHER_KEY",        &pi_MC_MOTHER_KEY);

	tree2->SetBranchAddress("D0K_ProbNNk",             &D0K_ProbNNk);
	tree2->SetBranchAddress("D0K_ProbNNpi",            &D0K_ProbNNpi);
	tree2->SetBranchAddress("D0pi_ProbNNk",            &D0pi_ProbNNk);
	tree2->SetBranchAddress("D0pi_ProbNNpi",           &D0pi_ProbNNpi);

	tree2->SetBranchAddress("B_L0HadronDecision_TOS",  &B_L0HadronDecision_TOS);
//	tree2->SetBranchAddress("hlt2global",&hlt2global);
//	tree2->SetBranchAddress("l0hadron",&l0hadron);
//	tree2->SetBranchAddress("l0decision",&l0decision);
//	tree2->SetBranchAddress("hlt2topo2TOS",&hlt2topo2TOS);
//	tree2->SetBranchAddress("hlt2topo3TOS",&hlt2topo3TOS);
//	tree2->SetBranchAddress("hlt2topo4TOS",&hlt2topo4TOS);
//
	tree2->SetBranchAddress("NN", &NN);

	for ( Int_t i(0); i < nEntries2; ++i ) {

		tree2->GetEntry( i );

		if(NN<NNmin || NN>NNmax) continue;

		if(B_L0HadronDecision_TOS!=1) continue;

//		++total;
		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
//                  cout << "event ignored at " << Bd_CM_m13Sq << ", " << Bd_CM_m23Sq << endl;
//                  ++ignored;
		  continue;
//		} else {
//                  cout << "event accepted at " << D0h1.M2() << ", " << D0h2.M2() << endl;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

//		if (Matched>0){
//		std::cout<<"A"<<std::endl;
		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0K_TRUEID)==321&&abs(D0pi_TRUEID)==211) {
//			std::cout<<"B"<<std::endl;
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0K_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0pi_MC_GD_MOTHER_KEY) {
//				std::cout<<"C"<<std::endl;
				dpcut->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
			}
		}	
	}

	file2->Close();

	cout << dpnocuts->GetEntries() << endl;
	cout << dpcut->GetEntries() << endl;
	dpnocuts->Sumw2();
	dpcut->Sumw2();
	eff->Divide(dpcut,dpnocuts);

	TEfficiency eff2(*dpcut,*dpnocuts);
	//eff2->SetStatisticOption(kFCP); //Default
	
	for(Int_t i=1; i<=nbins; ++i) {
		for(Int_t j=1; j<=nbins; ++j) {
			std::cout << i << "\t" << j << ":\t" << eff->GetBinContent(i,j) << "="
			<< eff2.GetEfficiency(eff2.GetGlobalBin(i,j)) << "\t+"
			<< eff2.GetEfficiencyErrorUp(eff2.GetGlobalBin(i,j)) << "\t-"
			<< eff2.GetEfficiencyErrorLow(eff2.GetGlobalBin(i,j)) << "\t" << std::endl;
		}
	}

	for(Int_t i=1; i<=nbins; ++i) {
		for(Int_t j=1; j<=nbins; ++j) {
			eff->SetBinContent(i,j,eff2.GetEfficiency(eff2.GetGlobalBin(i,j)));
			errp->SetBinContent(i,j,eff2.GetEfficiencyErrorUp(eff2.GetGlobalBin(i,j)));
			errm->SetBinContent(i,j,eff2.GetEfficiencyErrorLow(eff2.GetGlobalBin(i,j)));
		}
	}

	Float_t counter(0);
	Float_t sum(0);

	for (Int_t k=0;k<200;k++){
		for (Int_t l=0;l<200;l++) {
			Float_t temp2(-1);
			temp2 = eff->GetBinContent(k+1,l+1);
			if (temp2 >0) {
				sum = sum + temp2;
				counter = counter + 1;
			}
		}
	}

	cout << "Ignored " << ignored << " of " << total << " events." << endl;
	cout << counter << " " << sum << " " << 100*sum/counter << "%" << endl;
//	cout << counter2<< " " << sum2<< " " << 100*sum2/counter2 << "%" << endl;

	// Add a file to save the plots
	TString fname = "hists/seltrigTOS_";
	//fname += nbins;
	//fname += "_"; fname += NNmin1k; fname += "_"; fname += NNmax1k;
	fname += binLabel;
	fname += ".root";
	outfile = new TFile(fname,"RECREATE");
	eff->SetName("efficiency");
	eff->Write();
	errp->SetName("errorHi");
	errp->Write();
	errm->SetName("errorLo");
	errm->Write();
	outfile->Close();

	lhcbStyle();
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00};
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51};
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00};
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
	gStyle->SetOptStat(0000);

	TCanvas * canvas = new TCanvas("c0","c0");
	dpnocuts->Draw("colz");
	canvas->SetRightMargin(0.15);
	canvas->SaveAs("plots/seltrig_denom.pdf");
	TCanvas * canvas = new TCanvas("c1","c1");
        dpcut->Draw("colz");
	canvas->SetRightMargin(0.15);
	canvas->SaveAs("plots/seltrig_num.pdf");
	TCanvas * canvas = new TCanvas("c2","c2");
	eff->SetLabelFont(62,"x");
	eff->SetLabelFont(62,"y");
	eff->SetTitleFont(62,"x");
	eff->SetTitleFont(62,"y");
	eff->SetTitleSize(0.06,"x");
	eff->SetTitleSize(0.06,"y");
	eff->SetLabelSize(0.05,"x");
	eff->SetLabelSize(0.05,"y");
	//eff->SetXTitle("m^{2}(D^{0}#pi^{-}) (GeV/c^{2})^{2}");
	//eff->SetYTitle("m^{2}(D^{0}#pi^{+}) (GeV/c^{2})^{2}");
	eff->SetXTitle(xTitle);
	eff->SetYTitle(yTitle);
	eff->GetZaxis()->SetRangeUser(0,0.01);
//	if(fileName2=="DKpi_Bs2DKpi_selBs_noNN_noPID_noTrig") eff3->GetZaxis()->SetRangeUser(0,0.02*1.33);
	//gStyle->SetPalette(1);
	canvas->SetRightMargin(0.15);
	eff->Draw("colztext30");
//	kinematics->drawDPContour();
	//printLHCb("R","S");
	TString fname2 = "plots/seltrig_sq";
	fname2 += nbins;
	//fname2 += "_"; fname2 += NNmin1k; fname2 += "_"; fname2 += NNmax1k;
	fname2 += "_NN"; fname2 += binLabel;
	fname2 += "_eff.pdf";
	canvas->SaveAs(fname2);

	gStyle->SetPaintTextFormat("0.2f %%");

	errp->Draw("colztext30");
	TString fname2 = "plots/seltrig_sq";
	fname2 += nbins;
	//fname2 += "_"; fname2 += NNmin1k; fname2 += "_"; fname2 += NNmax1k;
	fname2 += "_NN"; fname2 += binLabel;
	fname2 += "_errp.pdf";
	canvas->SaveAs(fname2);

	errm->Draw("colztext30");
	TString fname2 = "plots/seltrig_sq";
	fname2 += nbins;
	//fname2 += "_"; fname2 += NNmin1k; fname2 += "_"; fname2 += NNmax1k;
	fname2 += "_NN"; fname2 += binLabel;
	fname2 += "_errm.pdf";
	canvas->SaveAs(fname2);

	delete kinematics;
}

void sel_eff_D2Kpi() {
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
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", -0.80,  1.00,"all");
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree", -0.80,  0.00,"a");
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",  0.00,  0.50,"b");
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",  0.50,  0.72,"c");
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",  0.72,  0.81,"d");
  sel_Eff("Sim08aMCEff","Bd2DKpi/Bd2DKpi_truth_P8","Bd2KpiD0/tupleMCTruth","Sim08aMC","D2Kpi/Bd2D0Kpi/D2Kpi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2Kpi_addIMs_addMisIDIMs","DecayTree",  0.81,  1.00,"e");
}
