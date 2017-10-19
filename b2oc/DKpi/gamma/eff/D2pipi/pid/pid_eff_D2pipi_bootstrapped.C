#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <TColor.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>

#include <Laura++/LauKinematics.hh>
#include <Laura++/LauDaughters.hh>

//class LauKinematics;

void pid_Eff(Int_t seed) {
	TString jobName="Sim08aMC";
	TString fileName="D2pipi/Bd2D0Kpi/D2pipi_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2pipi_addIMs_addMisIDIMs";
	
	TString treeName="DecayTree";

	TString xTitle="m'";
	TString yTitle="#theta'";
	
	Double_t NNmin=-0.8;
	Double_t NNmax=1.;
	TString binLabel="all";
	
	Int_t nbins(0);

	nbins = 10;

	double mB(5.27958), mD(1.86486), mK(0.49368), mPi(0.13957);
	LauKinematics* kinematics(0);

	kinematics = new LauKinematics(mD,mPi,mK,mB,true);

	TH2D * sdpcut = new TH2D("Dpipi CUT SDP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdpnocuts = new TH2D("Dpipi NO CUTS SDP","",nbins,0.0,1.0,nbins,0.0,1.0);
	TH2D * sdpeff = new TH2D("sdpeff","",nbins,0.0,1.0,nbins,0.0,1.0);
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

	std::vector<Int_t> entries;
	TRandom3 r(30000+seed);
	for(Int_t i=0; i<nEntries; ++i) {
           entries.push_back((Int_t)(r.Rndm()*nEntries));
	}

	std::sort(entries.begin(),entries.end());

	Double_t Bd_CM_m13Sq;
	Double_t Bd_CM_m23Sq;
	Double_t Bd_CM_m12Sq;
	Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0p_TRUEID,D0m_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0p_MC_GD_MOTHER_KEY,D0m_MC_GD_MOTHER_KEY;
	Float_t NN;
	
	tree->SetBranchAddress("Bd_CM_m13Sq",           &Bd_CM_m13Sq);
	tree->SetBranchAddress("Bd_CM_m23Sq",           &Bd_CM_m23Sq);
	tree->SetBranchAddress("Bd_CM_m12Sq",           &Bd_CM_m12Sq);
	tree->SetBranchAddress("B_TRUEID",                &B_TRUEID);
	tree->SetBranchAddress("D_TRUEID",                &D_TRUEID);
	tree->SetBranchAddress("D0p_TRUEID",              &D0p_TRUEID);
	tree->SetBranchAddress("D0m_TRUEID",              &D0m_TRUEID);
	tree->SetBranchAddress("K_TRUEID",                &K_TRUEID);
	tree->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);

	tree->SetBranchAddress("D0p_MC_GD_MOTHER_KEY",    &D0p_MC_GD_MOTHER_KEY);
	tree->SetBranchAddress("D0m_MC_GD_MOTHER_KEY",    &D0m_MC_GD_MOTHER_KEY);
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

	Int_t nTracks(0);

	Int_t NOCM_K_ID(0);
	
	Bool_t B_L0Global_TIS;

	Double_t weight(0);

	Double_t piEff(0);
	Double_t KEff(0);
	Double_t DpEff(0);
	Double_t DmEff(0);

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

	tree->SetBranchAddress("nTracks",           &nTracks);

	tree->SetBranchAddress("K_ID",              &NOCM_K_ID);

	tree->SetBranchAddress("B_L0Global_TIS",  &B_L0Global_TIS);

	// Get Histograms
	TFile * histfile    = TFile::Open("/home/phrkbf/PIDCalibHists/new/PerfHists_K_Strip20rX_MagBoth_DKpi_PPTN3_P_PT_nTracks.root");
	TH3F* K_hist        = (TH3F*)histfile->Get("K_ProbNNK * (1 - ProbNNpi ) > 0.3_All");
//	TH3F* DK_hist       = (TH3F*)histfile->Get("K_ProbNNK * (1 - ProbNNpi ) > 0.1_All");

	TFile * histfile2   = TFile::Open("/home/phrkbf/PIDCalibHists/new/PerfHists_Pi_Strip20rX_MagBoth_DKpi_PPTN3_P_PT_nTracks.root");
//	TH3F* pi_hist       = (TH3F*)histfile2->Get("Pi_ProbNNpi * (1 - ProbNNK ) * (1 - ProbNNp ) > 0.2_All");
	TH3F* pi_hist       = (TH3F*)histfile2->Get("Pi_ProbNNpi * (1 - ProbNNK ) > 0.2_All");
	TH3F* Dpi_hist      = (TH3F*)histfile2->Get("Pi_ProbNNpi * (1 - ProbNNK ) > 0.1_All");

	// Variables required for bin matching
	Int_t binPi(-1), binK(-1), binDp(-1), binDm(-1);

	for ( Int_t i(0); i < nEntries; ++i ) {
		if(i%10000==0) std::cout << i << " of " << nEntries << std::endl;

		tree->GetEntry( entries[i] );

		if(NN<NNmin || NN>NNmax) continue;

		if(!kinematics->withinDPLimits(Bd_CM_m13Sq,Bd_CM_m12Sq)) {
		  continue;
		}
	        kinematics->updateKinematics(Bd_CM_m13Sq,Bd_CM_m12Sq);

		if(abs(B_TRUEID)==511&&abs(D_TRUEID)==421&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(D0p_TRUEID)==211&&abs(D0m_TRUEID)==211) {
			if(K_MC_MOTHER_KEY==pi_MC_MOTHER_KEY&&K_MC_MOTHER_KEY==D0p_MC_GD_MOTHER_KEY&&K_MC_MOTHER_KEY==D0m_MC_GD_MOTHER_KEY) {
				weight=0;
				piEff=0;
				KEff=0;
	                        DpEff=0;
	                        DmEff=0;


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
				   binPi  = pi_hist->FindBin(D0_m_CM_P,D0_m_CM_PT,1.25*nTracks);
				   binK   = K_hist->FindBin(D0_p_CM_P,D0_p_CM_PT,1.25*nTracks);

				} else {
				   binPi  = pi_hist->FindBin(D0_p_CM_P,D0_p_CM_PT,1.25*nTracks);
				   binK   = K_hist->FindBin(D0_m_CM_P,D0_m_CM_PT,1.25*nTracks);
				}
				binDp  = Dpi_hist->FindBin(D0_D0p_CM_P,D0_D0p_CM_PT,1.25*nTracks);
				binDm  = Dpi_hist->FindBin(D0_D0m_CM_P,D0_D0m_CM_PT,1.25*nTracks);

				piEff  = pi_hist->GetBinContent(binPi);

				KEff   = K_hist->GetBinContent(binK);

				DpEff = Dpi_hist->GetBinContent(binDp);
				
				DmEff  = Dpi_hist->GetBinContent(binDm);

				weight = piEff * KEff * DpEff * DmEff;

				sdpnocuts->Fill(kinematics->getmPrime(),kinematics->getThetaPrime());
				sdpcut->Fill(kinematics->getmPrime(),kinematics->getThetaPrime(),weight);
			}
		}	
	}

	file->Close();

	sdpeff->Divide(sdpcut,sdpnocuts);

	// Add a file to save the plots
	TString outFileName("bootstrapped"); outFileName+="/pid_sq_Bd_effs";
	outFileName+=seed;
	outFileName+=".root";
	TFile* outfile = new TFile(outFileName,"RECREATE");
	sdpeff->SetName("efficiency");
	sdpeff->Write();
	outfile->Close();

	delete kinematics;
}

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
	  pid_Eff(seed);
  }
  return 0;
}
