#include <TLorentzVector.h>

class LauKinematics;

void l0h_Eff(Double_t NNmin=-0.80, Double_t NNmax=1.00){

	TString xTitle="m'";
	TString yTitle="#theta'";
	Int_t nbins(0);

	nbins = 5;

	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);
	int total(0), ignored(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * hSel          = new TH2D("hSel",          "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hCorr         = new TH2D("hCorr",         "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hNoOverlap    = new TH2D("hNoOverlap",    "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hNoCorr       = new TH2D("hNoCorr",       "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hEff          = new TH2D("hEff",          "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hEffSel       = new TH2D("hEffSel",       "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hRatio        = new TH2D("hRatio",        "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hErr2         = new TH2D("hErr2",         "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hErr          = new TH2D("hErr",          "",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * hOverlapRatio = new TH2D("hOverlapRatio", "",nbins,0.0,1.0,nbins,0.0,1.0);

	hNoCorr->Sumw2();
	hCorr->Sumw2();

	TFile * file = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/Bd2D0Kpi/D2KK_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2KK_addIMs_addMisIDIMs.root");

	TTree * tree = dynamic_cast<TTree*>( file->Get("DecayTree") );
	Int_t nEntries = tree->GetEntries();

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0p_TRUEID,D0m_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0p_MC_GD_MOTHER_KEY,D0m_MC_GD_MOTHER_KEY;
	Double_t pi_et, K_et, Dm_et, Dp_et;
	Double_t pi_etAdd, K_etAdd, Dm_etAdd, Dp_etAdd;
	Int_t pi_HCAL_region, K_HCAL_region, Dm_HCAL_region, Dp_HCAL_region;
	Double_t pi_x, pi_y, K_x, K_y, Dm_x, Dm_y, Dp_x, Dp_y;
	Bool_t l0_GlobalTIS, l0_HadronTOS;
	Int_t pi_ID, D0m_ID;
	Short_t Polarity;
	Float_t NN;

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

	tree->SetBranchAddress("pi_L0Calo_HCAL_realET",        &pi_et);
	tree->SetBranchAddress("pi_L0Calo_HCAL_region",        &pi_HCAL_region);
	tree->SetBranchAddress("pi_L0Calo_HCAL_xProjection",   &pi_x);
	tree->SetBranchAddress("pi_L0Calo_HCAL_yProjection",   &pi_y);
	tree->SetBranchAddress("K_L0Calo_HCAL_realET",         &K_et);
	tree->SetBranchAddress("K_L0Calo_HCAL_region",         &K_HCAL_region);
	tree->SetBranchAddress("K_L0Calo_HCAL_xProjection",    &K_x);
	tree->SetBranchAddress("K_L0Calo_HCAL_yProjection",    &K_y);
	tree->SetBranchAddress("D0m_L0Calo_HCAL_realET",      &Dm_et);
	tree->SetBranchAddress("D0m_L0Calo_HCAL_region",      &Dm_HCAL_region);
	tree->SetBranchAddress("D0m_L0Calo_HCAL_xProjection", &Dm_x);
	tree->SetBranchAddress("D0m_L0Calo_HCAL_yProjection", &Dm_y);
	tree->SetBranchAddress("D0p_L0Calo_HCAL_realET",       &Dp_et);
	tree->SetBranchAddress("D0p_L0Calo_HCAL_region",       &Dp_HCAL_region);
	tree->SetBranchAddress("D0p_L0Calo_HCAL_xProjection",  &Dp_x);
	tree->SetBranchAddress("D0p_L0Calo_HCAL_yProjection",  &Dp_y);
	tree->SetBranchAddress("B_L0Global_TIS",               &l0_GlobalTIS);
	tree->SetBranchAddress("B_L0HadronDecision_TOS",       &l0_HadronTOS);

	tree->SetBranchAddress("pi_ID",                   &pi_ID);
	tree->SetBranchAddress("D0m_ID",                 &D0m_ID);

	tree->SetBranchAddress("Polarity",                &Polarity);

	tree->SetBranchAddress("NN", &NN);

	// Get Histograms
	TFile * piFile1 = TFile::Open("tables/effs_data2012S20_MagnetDown_PiP_HadronThrs.root");
	TH1F  * pip_md_in  = (TH1F*)piFile1->Get("inner");
	TH1F  * pip_md_out = (TH1F*)piFile1->Get("outer");
	TFile * piFile2 = TFile::Open("tables/effs_data2012S20_MagnetDown_PiM_HadronThrs.root");
	TH1F  * pim_md_in  = (TH1F*)piFile2->Get("inner");
	TH1F  * pim_md_out = (TH1F*)piFile2->Get("outer");
	TFile * piFile3 = TFile::Open("tables/effs_data2012S20_MagnetUp_PiP_HadronThrs.root");
	TH1F  * pip_mu_in  = (TH1F*)piFile3->Get("inner");
	TH1F  * pip_mu_out = (TH1F*)piFile3->Get("outer");
	TFile * piFile4 = TFile::Open("tables/effs_data2012S20_MagnetUp_PiM_HadronThrs.root");
	TH1F  * pim_mu_in  = (TH1F*)piFile4->Get("inner");
	TH1F  * pim_mu_out = (TH1F*)piFile4->Get("outer");
	TFile * kFile1 = TFile::Open("tables/effs_data2012S20_MagnetDown_KaonP_HadronThrs.root");
	TH1F  * kp_md_in  = (TH1F*)kFile1->Get("inner");
	TH1F  * kp_md_out = (TH1F*)kFile1->Get("outer");
	TFile * kFile2 = TFile::Open("tables/effs_data2012S20_MagnetDown_KaonM_HadronThrs.root");
	TH1F  * km_md_in  = (TH1F*)kFile2->Get("inner");
	TH1F  * km_md_out = (TH1F*)kFile2->Get("outer");
	TFile * kFile3 = TFile::Open("tables/effs_data2012S20_MagnetUp_KaonP_HadronThrs.root");
	TH1F  * kp_mu_in  = (TH1F*)kFile3->Get("inner");
	TH1F  * kp_mu_out = (TH1F*)kFile3->Get("outer");
	TFile * kFile4 = TFile::Open("tables/effs_data2012S20_MagnetUp_KaonM_HadronThrs.root");
	TH1F  * km_mu_in  = (TH1F*)kFile4->Get("inner");
	TH1F  * km_mu_out = (TH1F*)kFile4->Get("outer");

	//new variables
	Double_t eff(0.), err2(0.);
	TH1F *h1(0), *h2(0), *h3(0), *h4(0);
	Double_t et1(0.), et2(0.), et3(0.), et4(0.);
	Double_t p1(0.), p2(0.), p1plus(0.), p2plus(0.), p3(0.), p4(0.), p12(0.);
	Double_t e1(0.), e2(0.), e1plus(0.), e2plus(0.), e3(0.), e4(0.), e12sq(0.);
	Double_t pFIRE1(0.), pFIRE2(0.), pFIRE1plus(0.), pFIRE2plus(0.), pFIRE3(0.), pFIRE4(0.), pHCAL(0.8);
	Int_t bin1(0.), bin2(0.), bin1plus(0.), bin2plus(0.), bin3(0.), bin4(0.);

	for ( Int_t i(0); i < nEntries; ++i ) {

		tree->GetEntry( i );

		if(NN<NNmin || NN>NNmax) continue;

		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
			continue;
		}
		kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0p_TRUEID)==321&&abs(D0m_TRUEID)==321) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0p_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0m_MC_GD_MOTHER_KEY) {
				if(l0_GlobalTIS==1) {
					if(l0_HadronTOS==1) {
						hSel->Fill(kinematics->getmPrime(), kinematics->getThetaPrime(), 1.0);
					}
					eff     = 0.;

					Int_t closestPair=0; //0 = piK, 1 = piDm, 2 = piDp, 3 = KDm, 4 = KDp, 5 = DmDp
					Double_t shortestDistX=abs(pi_x-K_x);
					Double_t shortestDistY=abs(pi_y-K_y);
					Double_t shortestDistSq=(pi_x-K_x)*(pi_x-K_x) + (pi_y-K_y)*(pi_y-K_y);
					Int_t regionA=pi_HCAL_region;
					Int_t regionB=K_HCAL_region;

					if((pi_x-Dm_x)*(pi_x-Dm_x) + (pi_y-Dm_y)*(pi_y-Dm_y) < shortestDistSq) {
						closestPair=1;
						regionA=pi_HCAL_region;
						regionB=Dm_HCAL_region;
						shortestDistX=abs(pi_x-Dm_x);
						shortestDistY=abs(pi_y-Dm_y);
						shortestDistSq=(pi_x-Dm_x)*(pi_x-Dm_x) + (pi_y-Dm_y)*(pi_y-Dm_y);
					}
					if((pi_x-Dp_x)*(pi_x-Dp_x) + (pi_y-Dp_y)*(pi_y-Dp_y) < shortestDistSq) {
						closestPair=2;
						regionA=pi_HCAL_region;
						regionB=Dp_HCAL_region;
						shortestDistX=abs(pi_x-Dp_x);
						shortestDistY=abs(pi_y-Dp_y);
						shortestDistSq=(pi_x-Dp_x)*(pi_x-Dp_x) + (pi_y-Dp_y)*(pi_y-Dp_y);
					}
					if((K_x-Dm_x)*(K_x-Dm_x) + (K_y-Dm_y)*(K_y-Dm_y) < shortestDistSq) {
						closestPair=3;
						regionA=K_HCAL_region;
						regionB=Dm_HCAL_region;
						shortestDistX=abs(K_x-Dm_x);
						shortestDistY=abs(K_y-Dm_y);
						shortestDistSq=(K_x-Dm_x)*(K_x-Dm_x) + (K_y-Dm_y)*(K_y-Dm_y);
					}
					if((K_x-Dp_x)*(K_x-Dp_x) + (K_y-Dp_y)*(K_y-Dp_y) < shortestDistSq) {
						closestPair=4;
						regionA=K_HCAL_region;
						regionB=Dp_HCAL_region;
						shortestDistX=abs(K_x-Dp_x);
						shortestDistY=abs(K_y-Dp_y);
						shortestDistSq=(K_x-Dp_x)*(K_x-Dp_x) + (K_y-Dp_y)*(K_y-Dp_y);
					}
					if((Dm_x-Dp_x)*(Dm_x-Dp_x) + (Dm_y-Dp_y)*(Dm_y-Dp_y) < shortestDistSq) {
						closestPair=5;
						regionA=Dm_HCAL_region;
						regionB=Dp_HCAL_region;
						shortestDistX=abs(Dm_x-Dp_x);
						shortestDistY=abs(Dm_y-Dp_y);
						shortestDistSq=(Dm_x-Dp_x)*(Dm_x-Dp_x) + (Dm_y-Dp_y)*(Dm_y-Dp_y);
					}

					Double_t cellSize = 0.;
					Bool_t combine=false;
					Double_t alpha=0.;
					if(regionA==1 && regionB==1) cellSize = 131.3; //inner
					else if(regionA==0 && regionB==0) cellSize = 262.6; //outer

					if(shortestDistX<0.5*cellSize && shortestDistY<0.5*cellSize) { //combine
						combine=true;
					} else if((shortestDistX<0.5*cellSize && shortestDistY<1.0*cellSize) ||
							(shortestDistX<1.0*cellSize && shortestDistY<0.5*cellSize)) { //overlap 50%
						alpha=0.5;
					} else if(shortestDistX<1.0*cellSize && shortestDistY<1.0*cellSize) { //overlap 25%
						alpha=0.25;
					}

					switch(closestPair) {
						case 0:
							et1 = pi_et;
							et2 = K_et;
							et3 = Dm_et;
							et4 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_md_out;
									else h1 = pip_md_in;
									if(K_HCAL_region==0) h2 = km_md_out;
									else h2 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_md_out;
									else h1 = pim_md_in;
									if(K_HCAL_region==0) h2 = kp_md_out;
									else h2 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h3 = pim_md_out;
								else h3 = pim_md_in;
								if(Dp_HCAL_region==0) h4 = pip_md_out;
								else h4 = pip_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_mu_out;
									else h1 = pip_mu_in;
									if(K_HCAL_region==0) h2 = km_mu_out;
									else h2 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_mu_out;
									else h1 = pim_mu_in;
									if(K_HCAL_region==0) h2 = kp_mu_out;
									else h2 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h3 = pim_mu_out;
								else h3 = pim_mu_in;
								if(Dp_HCAL_region==0) h4 = pip_mu_out;
								else h4 = pip_mu_in;
							}
							break;
						case 1:
							et1 = pi_et;
							et3 = K_et;
							et2 = Dm_et;
							et4 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_md_out;
									else h1 = pip_md_in;
									if(K_HCAL_region==0) h3 = km_md_out;
									else h3 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_md_out;
									else h1 = pim_md_in;
									if(K_HCAL_region==0) h3 = kp_md_out;
									else h3 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h2 = pim_md_out;
								else h2 = pim_md_in;
								if(Dp_HCAL_region==0) h4 = pip_md_out;
								else h4 = pip_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_mu_out;
									else h1 = pip_mu_in;
									if(K_HCAL_region==0) h3 = km_mu_out;
									else h3 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_mu_out;
									else h1 = pim_mu_in;
									if(K_HCAL_region==0) h3 = kp_mu_out;
									else h3 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h2 = km_mu_out;
								else h2 = km_mu_in;
								if(Dp_HCAL_region==0) h4 = kp_mu_out;
								else h4 = kp_mu_in;
							}
							break;
						case 2:
							et1 = pi_et;
							et4 = K_et;
							et3 = Dm_et;
							et2 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_md_out;
									else h1 = pip_md_in;
									if(K_HCAL_region==0) h4 = km_md_out;
									else h4 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_md_out;
									else h1 = pim_md_in;
									if(K_HCAL_region==0) h4 = kp_md_out;
									else h4 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h3 = km_md_out;
								else h3 = km_md_in;
								if(Dp_HCAL_region==0) h2 = kp_md_out;
								else h2 = kp_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h1 = pip_mu_out;
									else h1 = pip_mu_in;
									if(K_HCAL_region==0) h4 = km_mu_out;
									else h4 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h1 = pim_mu_out;
									else h1 = pim_mu_in;
									if(K_HCAL_region==0) h4 = kp_mu_out;
									else h4 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h3 = km_mu_out;
								else h3 = km_mu_in;
								if(Dp_HCAL_region==0) h2 = kp_mu_out;
								else h2 = kp_mu_in;
							}
							break;
						case 3:
							et3 = pi_et;
							et2 = K_et;
							et1 = Dm_et;
							et4 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h3 = pip_md_out;
									else h3 = pip_md_in;
									if(K_HCAL_region==0) h2 = km_md_out;
									else h2 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h3 = pim_md_out;
									else h3 = pim_md_in;
									if(K_HCAL_region==0) h2 = kp_md_out;
									else h2 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h1 = km_md_out;
								else h1 = km_md_in;
								if(Dp_HCAL_region==0) h4 = kp_md_out;
								else h4 = kp_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h3 = pip_mu_out;
									else h3 = pip_mu_in;
									if(K_HCAL_region==0) h2 = km_mu_out;
									else h2 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h3 = pim_mu_out;
									else h3 = pim_mu_in;
									if(K_HCAL_region==0) h2 = kp_mu_out;
									else h2 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h1 = km_mu_out;
								else h1 = km_mu_in;
								if(Dp_HCAL_region==0) h4 = kp_mu_out;
								else h4 = kp_mu_in;
							}
							break;
						case 4:
							et4 = pi_et;
							et2 = K_et;
							et3 = Dm_et;
							et1 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h4 = pip_md_out;
									else h4 = pip_md_in;
									if(K_HCAL_region==0) h2 = km_md_out;
									else h2 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h4 = pim_md_out;
									else h4 = pim_md_in;
									if(K_HCAL_region==0) h2 = kp_md_out;
									else h2 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h3 = km_md_out;
								else h3 = km_md_in;
								if(Dp_HCAL_region==0) h1 = kp_md_out;
								else h1 = kp_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h4 = pip_mu_out;
									else h4 = pip_mu_in;
									if(K_HCAL_region==0) h2 = km_mu_out;
									else h2 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h4 = pim_mu_out;
									else h4 = pim_mu_in;
									if(K_HCAL_region==0) h2 = kp_mu_out;
									else h2 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h3 = km_mu_out;
								else h3 = km_mu_in;
								if(Dp_HCAL_region==0) h1 = kp_mu_out;
								else h1 = kp_mu_in;
							}
							break;
						case 5:
							et3 = pi_et;
							et4 = K_et;
							et1 = Dm_et;
							et2 = Dp_et;
							if(Polarity==-1) {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h3 = pip_md_out;
									else h3 = pip_md_in;
									if(K_HCAL_region==0) h4 = km_md_out;
									else h4 = km_md_in;
								} else {
									if(pi_HCAL_region==0) h3 = pim_md_out;
									else h3 = pim_md_in;
									if(K_HCAL_region==0) h4 = kp_md_out;
									else h4 = kp_md_in;
								}
								if(Dm_HCAL_region==0) h1 = km_md_out;
								else h1 = km_md_in;
								if(Dp_HCAL_region==0) h2 = kp_md_out;
								else h2 = kp_md_in;
							} else {
								if(pi_ID>0) {
									if(pi_HCAL_region==0) h3 = pip_mu_out;
									else h3 = pip_mu_in;
									if(K_HCAL_region==0) h4 = km_mu_out;
									else h4 = km_mu_in;
								} else {
									if(pi_HCAL_region==0) h3 = pim_mu_out;
									else h3 = pim_mu_in;
									if(K_HCAL_region==0) h4 = kp_mu_out;
									else h4 = kp_mu_in;
								}
								if(Dm_HCAL_region==0) h1 = km_mu_out;
								else h1 = km_mu_in;
								if(Dp_HCAL_region==0) h2 = kp_mu_out;
								else h2 = kp_mu_in;
							}
							break;
					}

					bin1 = h1->FindBin(et1);
					p1   = h1->GetBinContent(bin1);
					e1   = h1->GetBinError(bin1);

					pFIRE1 = p1/pHCAL;
					if(pFIRE1>1.0) pFIRE1=1.0;

					bin2 = h2->FindBin(et2);
					p2   = h2->GetBinContent(bin2);
					e2   = h2->GetBinError(bin2);

					pFIRE2 = p2/pHCAL;
					if(pFIRE2>1.0) pFIRE2=1.0;

					bin3 = h3->FindBin(et3);
					p3   = h3->GetBinContent(bin3);
					e3   = h3->GetBinError(bin3);

					bin4 = h4->FindBin(et4);
					p4   = h4->GetBinContent(bin4);
					e4   = h4->GetBinError(bin4);

					if(combine) {
						bin1plus = h1->FindBin(et1+et2);
						p1plus   = h1->GetBinContent(bin1plus);
						e1plus   = h1->GetBinError(bin1plus);

						pFIRE1plus = p1plus/pHCAL;
						if(pFIRE1plus>1.0) pFIRE1plus=1.0;

						bin2plus = h2->FindBin(et1+et2);
						p2plus   = h2->GetBinContent(bin2plus);
						e2plus   = h2->GetBinError(bin2plus);

						pFIRE2plus = p2plus/pHCAL;
						if(pFIRE2plus>1.0) pFIRE2plus=1.0;

						//    P(fire on 1, 2 misses)   P(fire on 2, 1 misses)   P(fire on combined)
						p12   = pHCAL*(1-pHCAL)*pFIRE1 + pHCAL*(1-pHCAL)*pFIRE2 + pHCAL*pHCAL*(pFIRE1plus+pFIRE2plus)/2.0;
						e12sq = 0.;
						// note that pX = pFIREX*pHCAL
						e12sq = (1-pHCAL)*(1-pHCAL)*e1*e1 + (1-pHCAL)*(1-pHCAL)*e2*e2 
							+ pHCAL*pHCAL*e1plus*e1plus/(2.0*2.0)
							+ pHCAL*pHCAL*e2plus*e2plus/(2.0*2.0);

					} else {
						bin1plus = h1->FindBin(et1 + alpha*et2);
						p1plus   = h1->GetBinContent(bin1plus);
						e1plus   = h1->GetBinError(bin1plus);

						pFIRE1plus = p1plus/pHCAL;
						if(pFIRE1plus>1.0) pFIRE1plus=1.0;

						bin2plus = h2->FindBin(et2 + alpha*et1);
						p2plus   = h2->GetBinContent(bin2plus);
						e2plus   = h2->GetBinError(bin2plus);

						pFIRE2plus = p2plus/pHCAL;
						if(pFIRE2plus>1.0) pFIRE2plus=1.0;

						//    P(fire on 1, 2 misses)   P(fire on 2, 1 misses)   P(fire on 1 or 2)
						p12   = pHCAL*(1-pHCAL)*pFIRE1 + pHCAL*(1-pHCAL)*pFIRE2 + pHCAL*pHCAL*(1.-(1.-pFIRE1plus)*(1.-pFIRE2plus));
						// note that pX = pFIREX*pHCAL
						e12sq = (1-pHCAL)*(1-pHCAL)*e1*e1 + (1-pHCAL)*(1-pHCAL)*e2*e2 
							+ pHCAL*pHCAL*(1.-pFIRE2plus)*(1.-pFIRE2plus)*e1plus*e1plus
							+ pHCAL*pHCAL*(1.-pFIRE1plus)*(1.-pFIRE1plus)*e2plus*e2plus;

					}

					eff = 1.-(1.-p12)*(1.-p3)*(1.-p4);
					err2 = e12sq               * (1.-p3) * (1.-p3) * (1.-p4) * (1.-p4)
						+ (1.-p12) * (1.-p12) * e3      * e3      * (1.-p4) * (1.-p4)
						+ (1.-p12) * (1.-p12) * (1.-p3) * (1.-p3) * e4      * e4;

					hNoCorr->Fill(kinematics->getmPrime(), kinematics->getThetaPrime(), 1.0);
					hCorr->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), eff);
					hNoOverlap->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), 1.-(1.-p1)*(1.-p2)*(1.-p3)*(1.-p4));
					hErr2->Fill(  kinematics->getmPrime(), kinematics->getThetaPrime(), err2);

				}

			}
		}	
	}

	file->Close();

	hEff->Divide(hCorr,hNoCorr);
	hEffSel->Divide(hSel,hNoCorr);
	hRatio->Divide(hCorr,hSel);
	hOverlapRatio->Divide(hCorr,hNoOverlap);

	for (Int_t k=0; k<nbins; ++k) {
		for (Int_t l=0;l<nbins;++l) {
			Float_t temp2 = sqrt(hErr2->GetBinContent(k+1,l+1));
			hCorr->SetBinError(k+1,l+1,temp2);
			hEff->SetBinError(k+1,l+1,temp2/hNoCorr->GetBinContent(k+1,l+1));
			hRatio->SetBinError(k+1,l+1,temp2/hSel->GetBinContent(k+1,l+1));
			hErr->SetBinContent(k+1,l+1,temp2/hSel->GetBinContent(k+1,l+1));
		}
	}

	for (Int_t k=0; k<nbins; ++k) {
		for (Int_t l=0;l<nbins;++l) {
			cout << k+1 << "\t" << l+1 << "\t" << hCorr->GetBinContent(k+1,l+1) << "\t" << hCorr->GetBinError(k+1,l+1) << "\t" << hEff->GetBinContent(k+1,l+1) << "\t" << hEff->GetBinError(k+1,l+1) << "\t" << hRatio->GetBinContent(k+1,l+1) << "\t" << hRatio->GetBinError(k+1,l+1) << endl;
		}
	}

	// Add a file to save the plots
	TString fileName("hists/l0h_sq");
	fileName+=nbins;
	fileName+="_eff.root";
	outfile = new TFile(fileName,"RECREATE");
	hRatio->SetName("ratio");
	hRatio->Write();
	outfile->Close();

	TString fileName("hists/l0h_sq");
	fileName+=nbins;
	fileName+="_errTable.root";
	outfile = new TFile(fileName,"RECREATE");
	hErr->SetName("error");
	hErr->Write();
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

	TCanvas * canvas = new TCanvas("c2","c2");
	canvas->SetRightMargin(0.15);
	hEff->SetLabelFont(62,"x");
	hEff->SetLabelFont(62,"y");
	hEff->SetTitleFont(62,"x");
	hEff->SetTitleFont(62,"y");
	hEff->SetTitleSize(0.06,"x");
	hEff->SetTitleSize(0.06,"y");
	hEff->SetLabelSize(0.05,"x");
	hEff->SetLabelSize(0.05,"y");
	hEff->SetXTitle(xTitle);
	hEff->SetYTitle(yTitle);
	hEff->Draw("colz");
	TString sdpPlotName("plots/l0h_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_effData.pdf";
	canvas->SaveAs(sdpPlotName);

	TCanvas * canvas = new TCanvas("c2","c2");
	canvas->SetRightMargin(0.15);
	hEffSel->SetLabelFont(62,"x");
	hEffSel->SetLabelFont(62,"y");
	hEffSel->SetTitleFont(62,"x");
	hEffSel->SetTitleFont(62,"y");
	hEffSel->SetTitleSize(0.06,"x");
	hEffSel->SetTitleSize(0.06,"y");
	hEffSel->SetLabelSize(0.05,"x");
	hEffSel->SetLabelSize(0.05,"y");
	hEffSel->SetXTitle(xTitle);
	hEffSel->SetYTitle(yTitle);
	hEffSel->Draw("colz");
	TString sdpPlotName("plots/l0h_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_effMC.pdf";
	canvas->SaveAs(sdpPlotName);

	TCanvas * canvas = new TCanvas("c2","c2");
	canvas->SetRightMargin(0.15);
	hRatio->SetLabelFont(62,"x");
	hRatio->SetLabelFont(62,"y");
	hRatio->SetTitleFont(62,"x");
	hRatio->SetTitleFont(62,"y");
	hRatio->SetTitleSize(0.06,"x");
	hRatio->SetTitleSize(0.06,"y");
	hRatio->SetLabelSize(0.05,"x");
	hRatio->SetLabelSize(0.05,"y");
	hRatio->SetXTitle(xTitle);
	hRatio->SetYTitle(yTitle);
	hRatio->Draw("colz");
	TString sdpPlotName("plots/l0h_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_effRatio.pdf";
	canvas->SaveAs(sdpPlotName);

	TCanvas * canvas = new TCanvas("c2","c2");
	canvas->SetRightMargin(0.15);
	hOverlapRatio->SetLabelFont(62,"x");
	hOverlapRatio->SetLabelFont(62,"y");
	hOverlapRatio->SetTitleFont(62,"x");
	hOverlapRatio->SetTitleFont(62,"y");
	hOverlapRatio->SetTitleSize(0.06,"x");
	hOverlapRatio->SetTitleSize(0.06,"y");
	hOverlapRatio->SetLabelSize(0.05,"x");
	hOverlapRatio->SetLabelSize(0.05,"y");
	hOverlapRatio->SetXTitle(xTitle);
	hOverlapRatio->SetYTitle(yTitle);
	hOverlapRatio->Draw("colz");
	TString sdpPlotName("plots/l0h_sq");
	sdpPlotName+=nbins;
	sdpPlotName+="_effOverlapCorrection.pdf";
	canvas->SaveAs(sdpPlotName);

	delete kinematics;
}

void l0h_eff_D2KK() {
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
	l0h_Eff();
}
