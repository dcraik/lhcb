#include <TLorentzVector.h>

class LauKinematics;

void DecProdCut(TString jobName, TString fileName, TString treeName,
		TString xTitle="m^{2}(D^{0}K^{-}) (GeV/c^{2})^{2}", TString yTitle="m^{2}(D^{0}#pi^{+}) (GeV/c^{2})^{2}"){

	//double mB(5.27958), mBs(5.3662998046875), mD(1.86484), mK(0.493677001953125000), mPi(0.139570175170830879);
	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);
	int before(0), after(0);

	//kinematics = new LauKinematics(mK,mPi,mD,mBs,true);
	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job"+jobName+"/"+fileName+".root");

	TTree * tree = dynamic_cast<TTree*>( file->Get(treeName) );
	Int_t nEntries = tree->GetEntries();

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

	Double_t D_TRUEID(0);

	TLorentzVector D0;
	TLorentzVector h1;
	TLorentzVector h2;
	TLorentzVector D0h1;
	TLorentzVector D0h2;
	TLorentzVector h1h2;

	tree->SetBranchAddress("D0m_TRUETHETA",&D0K_TRUETHETA);
	tree->SetBranchAddress("D0p_TRUETHETA",&D0pi_TRUETHETA);
	tree->SetBranchAddress("h1_TRUETHETA",&h1_TRUETHETA);
	tree->SetBranchAddress("h2_TRUETHETA",&h2_TRUETHETA);

	tree->SetBranchAddress("D_CORRPX",&D0_CORRPX);
	tree->SetBranchAddress("D_CORRPY",&D0_CORRPY);
	tree->SetBranchAddress("D_CORRPZ",&D0_CORRPZ);
	tree->SetBranchAddress("D_CORRPE",&D0_CORRPE);

	tree->SetBranchAddress("h1_CORRPX",&h1_CORRPX);
	tree->SetBranchAddress("h1_CORRPY",&h1_CORRPY);
	tree->SetBranchAddress("h1_CORRPZ",&h1_CORRPZ);
	tree->SetBranchAddress("h1_CORRPE",&h1_CORRPE);

	tree->SetBranchAddress("h2_CORRPX",&h2_CORRPX);
	tree->SetBranchAddress("h2_CORRPY",&h2_CORRPY);
	tree->SetBranchAddress("h2_CORRPZ",&h2_CORRPZ);
	tree->SetBranchAddress("h2_CORRPE",&h2_CORRPE);

	tree->SetBranchAddress("D_TRUEID",&D_TRUEID);

	Int_t nbins(0);
//	Int_t maxX(0);
//	Int_t minX(0);
//	Int_t maxY(0);
//	Int_t minY(0);

	nbins = 15;
//	maxX = 28;
//	minX = 5;
//	maxY = 24;
//	minY = 3;

	TH2D * dpcut = new TH2D("Dpipi CUT SDP","",nbins,0.,1.,nbins,0.,1.);
	TH2D * dpnocuts = new TH2D("Dpipi NO CUTS SDP","",nbins,0.,1.,nbins,0.,1.);
	TH2D * eff = new TH2D("eff","",nbins,0.,1.,nbins,0.,1.);
	TH2D * errp = new TH2D("errp","",nbins,0.,1.,nbins,0.,1.);
	TH2D * errm = new TH2D("errm","",nbins,0.,1.,nbins,0.,1.);
	//TH2D * eff2 = new TH2D("eff2","",12,3,30,12,3,30);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		D0.SetPxPyPzE(D0_CORRPX, D0_CORRPY, D0_CORRPZ, D0_CORRPE);
		h1.SetPxPyPzE(h1_CORRPX, h1_CORRPY, h1_CORRPZ, h1_CORRPE);
		h2.SetPxPyPzE(h2_CORRPX, h2_CORRPY, h2_CORRPZ, h2_CORRPE);

		D0h1 = D0 + h1;
		D0h2 = D0 + h2;
		h1h2 = h1 + h2;

		//for DKpi SDP need Dpi and Kpi

//		if(jobName=="2011_Bs2DKpi_eff") {
//		  if(dpnocuts->FindBin(D0h1.M2()/1000000,D0h2.M2()/1000000)==587) {
//                     cout << "590\t" << D0h1.M2()/1000000 << "\t" << D0h2.M2()/1000000 << "\t" << kinematics->withinDPLimits(D0h1.M2()/1000000,D0h2.M2()/1000000) << endl;
//		  }
//		}

		++total;
		if(!kinematics->withinDPLimits(D0h1.M2()/1000000,h1h2.M2()/1000000)) {
//                  cout << "event ignored at " << D0h1.M2() << ", " << D0h2.M2() << endl;
                  ++ignored;
		  continue;
//		} else {
//                  cout << "event accepted at " << D0h1.M2() << ", " << D0h2.M2() << endl;
		}

//	        if(D0h1.M2()/1000000<12 && D0h2.M2()/1000000<8) cout << "foobar" <<endl;
                kinematics->updateKinematics(D0h1.M2()/1000000,h1h2.M2()/1000000);
		dpnocuts->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
		before++;

		//Make DecProdCut
		if (TMath::Abs( TMath::Sin( D0K_TRUETHETA ) ) <= TMath::Abs( TMath::Sin( 0.4 ))){
			if (TMath::Abs( TMath::Sin( D0K_TRUETHETA ) ) >= TMath::Abs( TMath::Sin( 0.01 ))){
				if (TMath::Abs( TMath::Sin( D0pi_TRUETHETA ) ) <= TMath::Abs( TMath::Sin( 0.4 ))){
					if (TMath::Abs( TMath::Sin( D0pi_TRUETHETA ) ) >= TMath::Abs( TMath::Sin( 0.01 ))){
						if (TMath::Abs( TMath::Sin( h1_TRUETHETA ) ) <= TMath::Abs( TMath::Sin( 0.4))){ 
							if (TMath::Abs( TMath::Sin( h1_TRUETHETA ) ) >= TMath::Abs( TMath::Sin( 0.01))){ 
								if (TMath::Abs( TMath::Sin( h2_TRUETHETA ) ) <= TMath::Abs( TMath::Sin( 0.4 ))){
									if (TMath::Abs( TMath::Sin( h2_TRUETHETA ) ) >= TMath::Abs( TMath::Sin( 0.01 ))){
										dpcut->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
										after++;
									}
								}
							}
						}
					}
				}
			}
		}
	}

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


//	for (Int_t i=0;i<200;i++){
//		for (Int_t j=0;j<200;j++) {
//			if(dpnocuts->GetBinContent(i+1,j+1) > 0 || dpcut->GetBinContent(i+1,j+1) > 0) {
//				if(dpcut->GetBinContent(i+1,j+1)/dpnocuts->GetBinContent(i+1,j+1) > 0.49) {
//				std::cout << i << "\t" << j << "\t" << dpcut->GetBinContent(i+1,j+1) << "/" << dpnocuts->GetBinContent(i+1,j+1) << "="
//				<< dpcut->GetBinContent(i+1,j+1)/dpnocuts->GetBinContent(i+1,j+1) << std::endl;
//				}
//			}
//		}
//	}

/*	for (Int_t i=0;i<200;i++){
		for (Int_t j=0;j<200;j++) {
			Int_t temp(-1);
			temp = dpcut->GetBinContent(i+1,j+1);
			if (temp<minInBin) {
				if(temp>0) cout << i+1 << "\t" << j+1 << "\t" << dpcut->GetBinContent(i+1,j+1) << "\t" << dpnocuts->GetBinContent(i+1,j+1) << endl;
				dpcut->SetBinContent(i+1,j+1,0);
			}
		}
	}
*/
	Float_t counter(0);
	Float_t sum(0);

//	for (Int_t i=0;i<200;i++){
//		for (Int_t j=0;j<200;j++) {
//		  if(eff3->GetBinContent(i+1,j+1)<1. && eff3->GetBinContent(i+1,j+1)>0.49) std::cout << i << "\t" << j << "\t" << eff3->GetBinContent(i+1,j+1) << std::endl;
//		}
//	}

//	cout << eff3->GetBinContent(587) <<"\t"<< dpnocuts->GetBinContent(587) <<"\t"<< dpcut->GetBinContent(587) << "\t" << eff3->GetBinContent(15,11) << endl;
//	Int_t *x(new Int_t), *y(new Int_t), *z(new Int_t);
//	eff3->GetBinXYZ(587,*x,*y,*z);
//	cout << "bin 587 is " << *x << "," << *y << "," << *z << endl;

//	eff3->SetBinContent(eff3->FindBin(11.75.,7.5),1.);
//        eff3->SetBinContent(591,1.);
//	cout << eff3->FindBin(11.75,7.5) << dpnocuts->FindBin(11.75,7.5) << dpcut->FindBin(11.75,7.5) << endl;

	lhcbStyle();
	gStyle->SetOptStat(0000);

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


	cout << before << "\t" << after << endl;

	cout << "Ignored " << ignored << " of " << total << " events." << endl;

//	cout << counter << " " << sum << " " << 100*sum/counter   << "%" << endl;

	TCanvas * canvas = new TCanvas("c1","c1");
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
	canvas->SetRightMargin(0.15);
	dpnocuts->Draw("colz");
	canvas->SaveAs("plots/sq15_gen.pdf");

	TCanvas * canvas = new TCanvas("c2","c2");
	eff->SetLabelFont(62,"x");
	eff->SetLabelFont(62,"y");
	eff->SetTitleFont(62,"x");
	eff->SetTitleFont(62,"y");
	eff->SetTitleSize(0.06,"x");
	eff->SetTitleSize(0.06,"y");
	eff->SetLabelSize(0.05,"x");
	eff->SetLabelSize(0.05,"y");
	eff->SetXTitle(xTitle);
	eff->SetYTitle(yTitle);
	eff->GetZaxis()->SetRangeUser(0.25,0.45);
	//eff->SetXTitle("m^{2}(D^{0}K^{-}) (GeV/c^{2})^{2}");
	//eff->SetYTitle("m^{2}(D^{0}K^{+}) (GeV/c^{2})^{2}");
	gStyle->SetPalette(1);
	canvas->SetRightMargin(0.15);
	//canvas->SetRightMargin(1);
	eff->Draw("colztext30");
	kinematics->drawDPContour();
	printLHCb("R","S");
	canvas->SaveAs("plots/geom_eff_sq15_eff.pdf");
//	TCanvas * canvas = new TCanvas("c3","c3");
//	eff->Draw("colz");

	gStyle->SetPaintTextFormat("0.2f %%");

	errp->Draw("colztext30");
	TString fname2 = "plots/geom_sq";
	fname2 += nbins;
	fname2 += "_errp.pdf";
	canvas->SaveAs(fname2);

	errm->Draw("colztext30");
	TString fname2 = "plots/geom_sq";
	fname2 += nbins;
	fname2 += "_errm.pdf";
	canvas->SaveAs(fname2);

	outfile = new TFile("hists/geom_sq15_all.root","RECREATE");
	eff->SetName("efficiency");
	eff->Write();
	errp->SetName("errorHi");
	errp->Write();
	errm->SetName("errorLo");
	errm->Write();
	outfile->Close();

	delete kinematics;
}

void geo_Eff_Bd_square() {
  gROOT->ProcessLine(".L lhcbStyle.C");
  gSystem->Load( "libEG" );
  gSystem->Load( "libFitEff.so" );
  lhcbStyle();
  gStyle->SetOptStat(0000);
  DecProdCut("Sim08aMCEff","D2KK/Bd2DKpi_truth_gen","Bd2KpiD0/tupleMCTruth");
}
