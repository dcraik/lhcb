//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 12 19:16:36 2018 by ROOT version 5.34/38
// from TTree DecayTree/DecayTree
// found on file: /tmp/dcraik/Tuples.root
//////////////////////////////////////////////////////////

#ifndef skim_h
#define skim_h

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>

//#include <boost/progress.hpp>

// Header file for the classes stored in the TTree if any.

class skim {
	public :
		// Fixed size dimensions of array or collections stored in the TTree if any.
		static const Int_t kMaxphi_ENDVERTEX_COV = 1;
		static const Int_t kMaxphi_OWNPV_COV = 1;
		static const Int_t kMaxKplus_OWNPV_COV = 1;
		static const Int_t kMaxKplus_ORIVX_COV = 1;
		static const Int_t kMaxKminus_OWNPV_COV = 1;
		static const Int_t kMaxKminus_ORIVX_COV = 1;

		static constexpr Int_t kMaxZ_ENDVERTEX_COV = 1;
		static constexpr Int_t kMaxZ_OWNPV_COV = 1;
		static constexpr Int_t kMaxmuplus_OWNPV_COV = 1;
		static constexpr Int_t kMaxmuplus_ORIVX_COV = 1;
		static constexpr Int_t kMaxmuminus_OWNPV_COV = 1;
		static constexpr Int_t kMaxmuminus_ORIVX_COV = 1;

		// Declaration of leaf types
		Double_t        phi_ENDVERTEX_X;
		Double_t        phi_ENDVERTEX_Y;
		Double_t        phi_ENDVERTEX_Z;
		Double_t        phi_ENDVERTEX_XERR;
		Double_t        phi_ENDVERTEX_YERR;
		Double_t        phi_ENDVERTEX_ZERR;
		Double_t        phi_ENDVERTEX_CHI2;
		Int_t           phi_ENDVERTEX_NDOF;
		Float_t         phi_ENDVERTEX_COV_[3][3];
		Double_t        phi_OWNPV_X;
		Double_t        phi_OWNPV_Y;
		Double_t        phi_OWNPV_Z;
		Double_t        phi_OWNPV_XERR;
		Double_t        phi_OWNPV_YERR;
		Double_t        phi_OWNPV_ZERR;
		Double_t        phi_OWNPV_CHI2;
		Int_t           phi_OWNPV_NDOF;
		Float_t         phi_OWNPV_COV_[3][3];
		Double_t        phi_IP_OWNPV;
		Double_t        phi_IPCHI2_OWNPV;
		Double_t        phi_FD_OWNPV;
		Double_t        phi_FDCHI2_OWNPV;
		Double_t        phi_DIRA_OWNPV;
		Double_t        phi_P;
		Double_t        phi_PT;
		Double_t        phi_PE;
		Double_t        phi_PX;
		Double_t        phi_PY;
		Double_t        phi_PZ;
		Double_t        phi_MM;
		Double_t        phi_MMERR;
		Double_t        phi_M;
		Int_t           phi_ID;
		Bool_t          phi_L0Global_Dec;
		Bool_t          phi_L0Global_TIS;
		Bool_t          phi_L0Global_TOS;
		Bool_t          phi_Hlt1Global_Dec;
		Bool_t          phi_Hlt1Global_TIS;
		Bool_t          phi_Hlt1Global_TOS;
		Bool_t          phi_Hlt1Phys_Dec;
		Bool_t          phi_Hlt1Phys_TIS;
		Bool_t          phi_Hlt1Phys_TOS;
		Bool_t          phi_Hlt2Global_Dec;
		Bool_t          phi_Hlt2Global_TIS;
		Bool_t          phi_Hlt2Global_TOS;
		Bool_t          phi_Hlt2Phys_Dec;
		Bool_t          phi_Hlt2Phys_TIS;
		Bool_t          phi_Hlt2Phys_TOS;
		//Bool_t          phi_L0DiHadron_lowMultDecision_Dec;
		//Bool_t          phi_L0DiHadron_lowMultDecision_TIS;
		//Bool_t          phi_L0DiHadron_lowMultDecision_TOS;
		//Bool_t          phi_L0HRCDiHadron_lowMultDecision_Dec;
		//Bool_t          phi_L0HRCDiHadron_lowMultDecision_TIS;
		//Bool_t          phi_L0HRCDiHadron_lowMultDecision_TOS;
		//Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		//Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		//Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		//Bool_t          phi_Hlt1LowMultHerschelDecision_Dec;
		//Bool_t          phi_Hlt1LowMultHerschelDecision_TIS;
		//Bool_t          phi_Hlt1LowMultHerschelDecision_TOS;
		//Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		//Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		//Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		//Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		//Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		//Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		//Bool_t          phi_Hlt2LowMultLMR2HHDecision_Dec;
		//Bool_t          phi_Hlt2LowMultLMR2HHDecision_TIS;
		//Bool_t          phi_Hlt2LowMultLMR2HHDecision_TOS;
		//Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_Dec;
		//Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_TIS;
		//Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_TOS;
		Bool_t		phi_L0DiHadron_lowMultDecision_Dec;
		Bool_t		phi_L0Muon_lowMultDecision_Dec;
		Bool_t		phi_L0DiMuon_lowMultDecision_Dec;
		Bool_t		phi_L0Electron_lowMultDecision_Dec;
		Bool_t		phi_L0Photon_lowMultDecision_Dec;
		Bool_t		phi_L0DiEM_lowMultDecision_Dec;
		Bool_t		phi_L0SPDDecision_Dec;
		Bool_t		phi_L0CALODecision_Dec;
		Bool_t		phi_L0PUDecision_Dec;
		Bool_t		phi_L0MuonDecision_Dec;
		Bool_t		phi_Hlt1LowMultPassThroughDecision_Dec;
		Bool_t		phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t		phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t		phi_Hlt1BBMicroBiasVeloDecision_Dec;
		Bool_t		phi_Hlt1BBHighMultDecision_Dec;
		Bool_t		phi_Hlt1MBNoBiasDecision_Dec;
		Bool_t		phi_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t		phi_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t		phi_Hlt2MBMicroBiasVeloDecision_Dec;
		Bool_t		phi_Hlt2MBHighMultDecision_Dec;
		Bool_t		phi_Hlt2MBNoBiasDecision_Dec;
		Bool_t		phi_Hlt2PassThroughDecision_Dec;
		Bool_t		phi_L0HadronDecision_Dec;
//		Bool_t		phi_L0MuonDecision_Dec;
//		Bool_t		phi_L0SPDDecision_Dec;
		Bool_t		phi_L0HadronLowMultDecision_Dec;
		Bool_t		phi_L0MuonLowMultDecision_Dec;
		Bool_t		phi_L0ElectronLowMultDecision_Dec;
		Bool_t		phi_L0PhotonLowMultDecision_Dec;
		Bool_t		phi_L0DiEMLowMultDecision_Dec;
		Bool_t		phi_L0SPDLowMultDecision_Dec;
		Bool_t		phi_L0SoftCEPDecision_Dec;
		Bool_t		phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec;
		Bool_t		phi_Hlt1BBMicroBiasSoftCEPDecision_Dec;
		Bool_t		phi_Hlt1BBHasTrackDecision_Dec;
		Bool_t		phi_Hlt2BBPassThroughDecision_Dec;
		Bool_t		phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec;
		Bool_t		phi_Hlt2BBLongTrackDecision_Dec;
		Bool_t		phi_Hlt2SingleTrackDecision_Dec;
//		Bool_t		phi_Hlt2MBMicroBiasVeloDecision_Dec;
		Double_t        Kplus_MC12TuneV2_ProbNNe;
		Double_t        Kplus_MC12TuneV2_ProbNNmu;
		Double_t        Kplus_MC12TuneV2_ProbNNpi;
		Double_t        Kplus_MC12TuneV2_ProbNNk;
		Double_t        Kplus_MC12TuneV2_ProbNNp;
		Double_t        Kplus_MC12TuneV2_ProbNNghost;
		Double_t        Kplus_MC12TuneV3_ProbNNe;
		Double_t        Kplus_MC12TuneV3_ProbNNmu;
		Double_t        Kplus_MC12TuneV3_ProbNNpi;
		Double_t        Kplus_MC12TuneV3_ProbNNk;
		Double_t        Kplus_MC12TuneV3_ProbNNp;
		Double_t        Kplus_MC12TuneV3_ProbNNghost;
		Double_t        Kplus_MC12TuneV4_ProbNNe;
		Double_t        Kplus_MC12TuneV4_ProbNNmu;
		Double_t        Kplus_MC12TuneV4_ProbNNpi;
		Double_t        Kplus_MC12TuneV4_ProbNNk;
		Double_t        Kplus_MC12TuneV4_ProbNNp;
		Double_t        Kplus_MC12TuneV4_ProbNNghost;
		Double_t        Kplus_MC15TuneV1_ProbNNe;
		Double_t        Kplus_MC15TuneV1_ProbNNmu;
		Double_t        Kplus_MC15TuneV1_ProbNNpi;
		Double_t        Kplus_MC15TuneV1_ProbNNk;
		Double_t        Kplus_MC15TuneV1_ProbNNp;
		Double_t        Kplus_MC15TuneV1_ProbNNghost;
		Double_t        Kplus_OWNPV_X;
		Double_t        Kplus_OWNPV_Y;
		Double_t        Kplus_OWNPV_Z;
		Double_t        Kplus_OWNPV_XERR;
		Double_t        Kplus_OWNPV_YERR;
		Double_t        Kplus_OWNPV_ZERR;
		Double_t        Kplus_OWNPV_CHI2;
		Int_t           Kplus_OWNPV_NDOF;
		Float_t         Kplus_OWNPV_COV_[3][3];
		Double_t        Kplus_IP_OWNPV;
		Double_t        Kplus_IPCHI2_OWNPV;
		Double_t        Kplus_ORIVX_X;
		Double_t        Kplus_ORIVX_Y;
		Double_t        Kplus_ORIVX_Z;
		Double_t        Kplus_ORIVX_XERR;
		Double_t        Kplus_ORIVX_YERR;
		Double_t        Kplus_ORIVX_ZERR;
		Double_t        Kplus_ORIVX_CHI2;
		Int_t           Kplus_ORIVX_NDOF;
		Float_t         Kplus_ORIVX_COV_[3][3];
		Double_t        Kplus_P;
		Double_t        Kplus_PT;
		Double_t        Kplus_PE;
		Double_t        Kplus_PX;
		Double_t        Kplus_PY;
		Double_t        Kplus_PZ;
		Double_t        Kplus_M;
		Int_t           Kplus_ID;
		Double_t        Kplus_PIDe;
		Double_t        Kplus_PIDmu;
		Double_t        Kplus_PIDK;
		Double_t        Kplus_PIDp;
		Double_t        Kplus_ProbNNe;
		Double_t        Kplus_ProbNNk;
		Double_t        Kplus_ProbNNp;
		Double_t        Kplus_ProbNNpi;
		Double_t        Kplus_ProbNNmu;
		Double_t        Kplus_ProbNNghost;
		Bool_t          Kplus_hasMuon;
		Bool_t          Kplus_isMuon;
		Bool_t          Kplus_hasRich;
		Bool_t          Kplus_UsedRichAerogel;
		Bool_t          Kplus_UsedRich1Gas;
		Bool_t          Kplus_UsedRich2Gas;
		Bool_t          Kplus_RichAboveElThres;
		Bool_t          Kplus_RichAboveMuThres;
		Bool_t          Kplus_RichAbovePiThres;
		Bool_t          Kplus_RichAboveKaThres;
		Bool_t          Kplus_RichAbovePrThres;
		Bool_t          Kplus_hasCalo;
		Double_t        Kplus_PP_CombDLLe;
		Double_t        Kplus_PP_CombDLLmu;
		Double_t        Kplus_PP_CombDLLpi;
		Double_t        Kplus_PP_CombDLLk;
		Double_t        Kplus_PP_CombDLLp;
		Double_t        Kplus_PP_CombDLLd;
		Double_t        Kplus_PP_ProbNNe;
		Double_t        Kplus_PP_ProbNNmu;
		Double_t        Kplus_PP_ProbNNpi;
		Double_t        Kplus_PP_ProbNNk;
		Double_t        Kplus_PP_ProbNNp;
		Double_t        Kplus_PP_ProbNNghost;
		Double_t        Kplus_PP_ProbNNd;
		Bool_t          Kplus_L0Global_Dec;
		Bool_t          Kplus_L0Global_TIS;
		Bool_t          Kplus_L0Global_TOS;
		Bool_t          Kplus_Hlt1Global_Dec;
		Bool_t          Kplus_Hlt1Global_TIS;
		Bool_t          Kplus_Hlt1Global_TOS;
		Bool_t          Kplus_Hlt1Phys_Dec;
		Bool_t          Kplus_Hlt1Phys_TIS;
		Bool_t          Kplus_Hlt1Phys_TOS;
		Bool_t          Kplus_Hlt2Global_Dec;
		Bool_t          Kplus_Hlt2Global_TIS;
		Bool_t          Kplus_Hlt2Global_TOS;
		Bool_t          Kplus_Hlt2Phys_Dec;
		Bool_t          Kplus_Hlt2Phys_TIS;
		Bool_t          Kplus_Hlt2Phys_TOS;
		Bool_t          Kplus_L0DiHadron_lowMultDecision_Dec;
		Bool_t          Kplus_L0DiHadron_lowMultDecision_TIS;
		Bool_t          Kplus_L0DiHadron_lowMultDecision_TOS;
		Bool_t          Kplus_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          Kplus_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          Kplus_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          Kplus_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          Kplus_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          Kplus_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          Kplus_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          Kplus_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          Kplus_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          Kplus_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          Kplus_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          Kplus_Hlt2LowMultLMR2HHWSDecision_TOS;
		Int_t           Kplus_TRACK_Type;
		Int_t           Kplus_TRACK_Key;
		Double_t        Kplus_TRACK_CHI2NDOF;
		Double_t        Kplus_TRACK_PCHI2;
		Double_t        Kplus_TRACK_MatchCHI2;
		Double_t        Kplus_TRACK_GhostProb;
		Double_t        Kplus_TRACK_CloneDist;
		Double_t        Kplus_TRACK_Likelihood;
		Double_t        Kminus_MC12TuneV2_ProbNNe;
		Double_t        Kminus_MC12TuneV2_ProbNNmu;
		Double_t        Kminus_MC12TuneV2_ProbNNpi;
		Double_t        Kminus_MC12TuneV2_ProbNNk;
		Double_t        Kminus_MC12TuneV2_ProbNNp;
		Double_t        Kminus_MC12TuneV2_ProbNNghost;
		Double_t        Kminus_MC12TuneV3_ProbNNe;
		Double_t        Kminus_MC12TuneV3_ProbNNmu;
		Double_t        Kminus_MC12TuneV3_ProbNNpi;
		Double_t        Kminus_MC12TuneV3_ProbNNk;
		Double_t        Kminus_MC12TuneV3_ProbNNp;
		Double_t        Kminus_MC12TuneV3_ProbNNghost;
		Double_t        Kminus_MC12TuneV4_ProbNNe;
		Double_t        Kminus_MC12TuneV4_ProbNNmu;
		Double_t        Kminus_MC12TuneV4_ProbNNpi;
		Double_t        Kminus_MC12TuneV4_ProbNNk;
		Double_t        Kminus_MC12TuneV4_ProbNNp;
		Double_t        Kminus_MC12TuneV4_ProbNNghost;
		Double_t        Kminus_MC15TuneV1_ProbNNe;
		Double_t        Kminus_MC15TuneV1_ProbNNmu;
		Double_t        Kminus_MC15TuneV1_ProbNNpi;
		Double_t        Kminus_MC15TuneV1_ProbNNk;
		Double_t        Kminus_MC15TuneV1_ProbNNp;
		Double_t        Kminus_MC15TuneV1_ProbNNghost;
		Double_t        Kminus_OWNPV_X;
		Double_t        Kminus_OWNPV_Y;
		Double_t        Kminus_OWNPV_Z;
		Double_t        Kminus_OWNPV_XERR;
		Double_t        Kminus_OWNPV_YERR;
		Double_t        Kminus_OWNPV_ZERR;
		Double_t        Kminus_OWNPV_CHI2;
		Int_t           Kminus_OWNPV_NDOF;
		Float_t         Kminus_OWNPV_COV_[3][3];
		Double_t        Kminus_IP_OWNPV;
		Double_t        Kminus_IPCHI2_OWNPV;
		Double_t        Kminus_ORIVX_X;
		Double_t        Kminus_ORIVX_Y;
		Double_t        Kminus_ORIVX_Z;
		Double_t        Kminus_ORIVX_XERR;
		Double_t        Kminus_ORIVX_YERR;
		Double_t        Kminus_ORIVX_ZERR;
		Double_t        Kminus_ORIVX_CHI2;
		Int_t           Kminus_ORIVX_NDOF;
		Float_t         Kminus_ORIVX_COV_[3][3];
		Double_t        Kminus_P;
		Double_t        Kminus_PT;
		Double_t        Kminus_PE;
		Double_t        Kminus_PX;
		Double_t        Kminus_PY;
		Double_t        Kminus_PZ;
		Double_t        Kminus_M;
		Int_t           Kminus_ID;
		Double_t        Kminus_PIDe;
		Double_t        Kminus_PIDmu;
		Double_t        Kminus_PIDK;
		Double_t        Kminus_PIDp;
		Double_t        Kminus_ProbNNe;
		Double_t        Kminus_ProbNNk;
		Double_t        Kminus_ProbNNp;
		Double_t        Kminus_ProbNNpi;
		Double_t        Kminus_ProbNNmu;
		Double_t        Kminus_ProbNNghost;
		Bool_t          Kminus_hasMuon;
		Bool_t          Kminus_isMuon;
		Bool_t          Kminus_hasRich;
		Bool_t          Kminus_UsedRichAerogel;
		Bool_t          Kminus_UsedRich1Gas;
		Bool_t          Kminus_UsedRich2Gas;
		Bool_t          Kminus_RichAboveElThres;
		Bool_t          Kminus_RichAboveMuThres;
		Bool_t          Kminus_RichAbovePiThres;
		Bool_t          Kminus_RichAboveKaThres;
		Bool_t          Kminus_RichAbovePrThres;
		Bool_t          Kminus_hasCalo;
		Double_t        Kminus_PP_CombDLLe;
		Double_t        Kminus_PP_CombDLLmu;
		Double_t        Kminus_PP_CombDLLpi;
		Double_t        Kminus_PP_CombDLLk;
		Double_t        Kminus_PP_CombDLLp;
		Double_t        Kminus_PP_CombDLLd;
		Double_t        Kminus_PP_ProbNNe;
		Double_t        Kminus_PP_ProbNNmu;
		Double_t        Kminus_PP_ProbNNpi;
		Double_t        Kminus_PP_ProbNNk;
		Double_t        Kminus_PP_ProbNNp;
		Double_t        Kminus_PP_ProbNNghost;
		Double_t        Kminus_PP_ProbNNd;
		Bool_t          Kminus_L0Global_Dec;
		Bool_t          Kminus_L0Global_TIS;
		Bool_t          Kminus_L0Global_TOS;
		Bool_t          Kminus_Hlt1Global_Dec;
		Bool_t          Kminus_Hlt1Global_TIS;
		Bool_t          Kminus_Hlt1Global_TOS;
		Bool_t          Kminus_Hlt1Phys_Dec;
		Bool_t          Kminus_Hlt1Phys_TIS;
		Bool_t          Kminus_Hlt1Phys_TOS;
		Bool_t          Kminus_Hlt2Global_Dec;
		Bool_t          Kminus_Hlt2Global_TIS;
		Bool_t          Kminus_Hlt2Global_TOS;
		Bool_t          Kminus_Hlt2Phys_Dec;
		Bool_t          Kminus_Hlt2Phys_TIS;
		Bool_t          Kminus_Hlt2Phys_TOS;
		Bool_t          Kminus_L0DiHadron_lowMultDecision_Dec;
		Bool_t          Kminus_L0DiHadron_lowMultDecision_TIS;
		Bool_t          Kminus_L0DiHadron_lowMultDecision_TOS;
		Bool_t          Kminus_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          Kminus_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          Kminus_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          Kminus_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          Kminus_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          Kminus_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          Kminus_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          Kminus_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          Kminus_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          Kminus_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          Kminus_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          Kminus_Hlt2LowMultLMR2HHWSDecision_TOS;
		Int_t           Kminus_TRACK_Type;
		Int_t           Kminus_TRACK_Key;
		Double_t        Kminus_TRACK_CHI2NDOF;
		Double_t        Kminus_TRACK_PCHI2;
		Double_t        Kminus_TRACK_MatchCHI2;
		Double_t        Kminus_TRACK_GhostProb;
		Double_t        Kminus_TRACK_CloneDist;
		Double_t        Kminus_TRACK_Likelihood;
		Double_t        Z_ENDVERTEX_X;
		Double_t        Z_ENDVERTEX_Y;
		Double_t        Z_ENDVERTEX_Z;
		Double_t        Z_ENDVERTEX_XERR;
		Double_t        Z_ENDVERTEX_YERR;
		Double_t        Z_ENDVERTEX_ZERR;
		Double_t        Z_ENDVERTEX_CHI2;
		Int_t           Z_ENDVERTEX_NDOF;
		Float_t         Z_ENDVERTEX_COV_[3][3];
		Double_t        Z_OWNPV_X;
		Double_t        Z_OWNPV_Y;
		Double_t        Z_OWNPV_Z;
		Double_t        Z_OWNPV_XERR;
		Double_t        Z_OWNPV_YERR;
		Double_t        Z_OWNPV_ZERR;
		Double_t        Z_OWNPV_CHI2;
		Int_t           Z_OWNPV_NDOF;
		Float_t         Z_OWNPV_COV_[3][3];
		Double_t        Z_IP_OWNPV;
		Double_t        Z_IPCHI2_OWNPV;
		Double_t        Z_FD_OWNPV;
		Double_t        Z_FDCHI2_OWNPV;
		Double_t        Z_DIRA_OWNPV;
		Double_t        Z_P;
		Double_t        Z_PT;
		Double_t        Z_PE;
		Double_t        Z_PX;
		Double_t        Z_PY;
		Double_t        Z_PZ;
		Double_t        Z_MM;
		Double_t        Z_MMERR;
		Double_t        Z_M;
		Int_t           Z_ID;
		Bool_t          Z_L0Global_Dec;
		Bool_t          Z_L0Global_TIS;
		Bool_t          Z_L0Global_TOS;
		Bool_t          Z_Hlt1Global_Dec;
		Bool_t          Z_Hlt1Global_TIS;
		Bool_t          Z_Hlt1Global_TOS;
		Bool_t          Z_Hlt1Phys_Dec;
		Bool_t          Z_Hlt1Phys_TIS;
		Bool_t          Z_Hlt1Phys_TOS;
		Bool_t          Z_Hlt2Global_Dec;
		Bool_t          Z_Hlt2Global_TIS;
		Bool_t          Z_Hlt2Global_TOS;
		Bool_t          Z_Hlt2Phys_Dec;
		Bool_t          Z_Hlt2Phys_TIS;
		Bool_t          Z_Hlt2Phys_TOS;
		Bool_t          Z_L0DiHadron_lowMultDecision_Dec;
		Bool_t          Z_L0DiHadron_lowMultDecision_TIS;
		Bool_t          Z_L0DiHadron_lowMultDecision_TOS;
		Bool_t          Z_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          Z_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          Z_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          Z_L0MuonDecision_Dec;
		Bool_t          Z_L0MuonDecision_TIS;
		Bool_t          Z_L0MuonDecision_TOS;
		Bool_t          Z_L0MuonEWDecision_Dec;
		Bool_t          Z_L0MuonEWDecision_TIS;
		Bool_t          Z_L0MuonEWDecision_TOS;
		Bool_t          Z_L0DiMuonDecision_Dec;
		Bool_t          Z_L0DiMuonDecision_TIS;
		Bool_t          Z_L0DiMuonDecision_TOS;
		Bool_t          Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          Z_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          Z_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          Z_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          Z_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          Z_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          Z_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          Z_Hlt1TrackMVADecision_Dec;
		Bool_t          Z_Hlt1TrackMVADecision_TIS;
		Bool_t          Z_Hlt1TrackMVADecision_TOS;
		Bool_t          Z_Hlt1TwoTrackMVADecision_Dec;
		Bool_t          Z_Hlt1TwoTrackMVADecision_TIS;
		Bool_t          Z_Hlt1TwoTrackMVADecision_TOS;
		Bool_t          Z_Hlt1SingleMuonHighPTDecision_Dec;
		Bool_t          Z_Hlt1SingleMuonHighPTDecision_TIS;
		Bool_t          Z_Hlt1SingleMuonHighPTDecision_TOS;
		Bool_t          Z_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          Z_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          Z_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          Z_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          Z_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          Z_Hlt2LowMultLMR2HHWSDecision_TOS;
		Bool_t          Z_Hlt2SingleMuonHighPTDecision_Dec;
		Bool_t          Z_Hlt2SingleMuonHighPTDecision_TIS;
		Bool_t          Z_Hlt2SingleMuonHighPTDecision_TOS;
		Double_t        muplus_MC12TuneV2_ProbNNe;
		Double_t        muplus_MC12TuneV2_ProbNNmu;
		Double_t        muplus_MC12TuneV2_ProbNNpi;
		Double_t        muplus_MC12TuneV2_ProbNNk;
		Double_t        muplus_MC12TuneV2_ProbNNp;
		Double_t        muplus_MC12TuneV2_ProbNNghost;
		Double_t        muplus_MC12TuneV3_ProbNNe;
		Double_t        muplus_MC12TuneV3_ProbNNmu;
		Double_t        muplus_MC12TuneV3_ProbNNpi;
		Double_t        muplus_MC12TuneV3_ProbNNk;
		Double_t        muplus_MC12TuneV3_ProbNNp;
		Double_t        muplus_MC12TuneV3_ProbNNghost;
		Double_t        muplus_MC12TuneV4_ProbNNe;
		Double_t        muplus_MC12TuneV4_ProbNNmu;
		Double_t        muplus_MC12TuneV4_ProbNNpi;
		Double_t        muplus_MC12TuneV4_ProbNNk;
		Double_t        muplus_MC12TuneV4_ProbNNp;
		Double_t        muplus_MC12TuneV4_ProbNNghost;
		Double_t        muplus_MC15TuneV1_ProbNNe;
		Double_t        muplus_MC15TuneV1_ProbNNmu;
		Double_t        muplus_MC15TuneV1_ProbNNpi;
		Double_t        muplus_MC15TuneV1_ProbNNk;
		Double_t        muplus_MC15TuneV1_ProbNNp;
		Double_t        muplus_MC15TuneV1_ProbNNghost;
		Double_t        muplus_OWNPV_X;
		Double_t        muplus_OWNPV_Y;
		Double_t        muplus_OWNPV_Z;
		Double_t        muplus_OWNPV_XERR;
		Double_t        muplus_OWNPV_YERR;
		Double_t        muplus_OWNPV_ZERR;
		Double_t        muplus_OWNPV_CHI2;
		Int_t           muplus_OWNPV_NDOF;
		Float_t         muplus_OWNPV_COV_[3][3];
		Double_t        muplus_IP_OWNPV;
		Double_t        muplus_IPCHI2_OWNPV;
		Double_t        muplus_ORIVX_X;
		Double_t        muplus_ORIVX_Y;
		Double_t        muplus_ORIVX_Z;
		Double_t        muplus_ORIVX_XERR;
		Double_t        muplus_ORIVX_YERR;
		Double_t        muplus_ORIVX_ZERR;
		Double_t        muplus_ORIVX_CHI2;
		Int_t           muplus_ORIVX_NDOF;
		Float_t         muplus_ORIVX_COV_[3][3];
		Double_t        muplus_P;
		Double_t        muplus_PT;
		Double_t        muplus_PE;
		Double_t        muplus_PX;
		Double_t        muplus_PY;
		Double_t        muplus_PZ;
		Double_t        muplus_M;
		Int_t           muplus_ID;
		Double_t        muplus_PIDe;
		Double_t        muplus_PIDmu;
		Double_t        muplus_PIDK;
		Double_t        muplus_PIDp;
		Double_t        muplus_PIDd;
		Double_t        muplus_ProbNNe;
		Double_t        muplus_ProbNNk;
		Double_t        muplus_ProbNNp;
		Double_t        muplus_ProbNNpi;
		Double_t        muplus_ProbNNmu;
		Double_t        muplus_ProbNNd;
		Double_t        muplus_ProbNNghost;
		Bool_t          muplus_hasMuon;
		Bool_t          muplus_isMuon;
		Bool_t          muplus_hasRich;
		Bool_t          muplus_UsedRichAerogel;
		Bool_t          muplus_UsedRich1Gas;
		Bool_t          muplus_UsedRich2Gas;
		Bool_t          muplus_RichAboveElThres;
		Bool_t          muplus_RichAboveMuThres;
		Bool_t          muplus_RichAbovePiThres;
		Bool_t          muplus_RichAboveKaThres;
		Bool_t          muplus_RichAbovePrThres;
		Bool_t          muplus_hasCalo;
		Double_t        muplus_PP_CombDLLe;
		Double_t        muplus_PP_CombDLLmu;
		Double_t        muplus_PP_CombDLLpi;
		Double_t        muplus_PP_CombDLLk;
		Double_t        muplus_PP_CombDLLp;
		Double_t        muplus_PP_CombDLLd;
		Double_t        muplus_PP_ProbNNe;
		Double_t        muplus_PP_ProbNNmu;
		Double_t        muplus_PP_ProbNNpi;
		Double_t        muplus_PP_ProbNNk;
		Double_t        muplus_PP_ProbNNp;
		Double_t        muplus_PP_ProbNNghost;
		Double_t        muplus_PP_ProbNNd;
		Bool_t          muplus_L0Global_Dec;
		Bool_t          muplus_L0Global_TIS;
		Bool_t          muplus_L0Global_TOS;
		Bool_t          muplus_Hlt1Global_Dec;
		Bool_t          muplus_Hlt1Global_TIS;
		Bool_t          muplus_Hlt1Global_TOS;
		Bool_t          muplus_Hlt1Phys_Dec;
		Bool_t          muplus_Hlt1Phys_TIS;
		Bool_t          muplus_Hlt1Phys_TOS;
		Bool_t          muplus_Hlt2Global_Dec;
		Bool_t          muplus_Hlt2Global_TIS;
		Bool_t          muplus_Hlt2Global_TOS;
		Bool_t          muplus_Hlt2Phys_Dec;
		Bool_t          muplus_Hlt2Phys_TIS;
		Bool_t          muplus_Hlt2Phys_TOS;
		Bool_t          muplus_L0DiHadron_lowMultDecision_Dec;
		Bool_t          muplus_L0DiHadron_lowMultDecision_TIS;
		Bool_t          muplus_L0DiHadron_lowMultDecision_TOS;
		Bool_t          muplus_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          muplus_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          muplus_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          muplus_L0MuonDecision_Dec;
		Bool_t          muplus_L0MuonDecision_TIS;
		Bool_t          muplus_L0MuonDecision_TOS;
		Bool_t          muplus_L0MuonEWDecision_Dec;
		Bool_t          muplus_L0MuonEWDecision_TIS;
		Bool_t          muplus_L0MuonEWDecision_TOS;
		Bool_t          muplus_L0DiMuonDecision_Dec;
		Bool_t          muplus_L0DiMuonDecision_TIS;
		Bool_t          muplus_L0DiMuonDecision_TOS;
		Bool_t          muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          muplus_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          muplus_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          muplus_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          muplus_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          muplus_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          muplus_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          muplus_Hlt1TrackMVADecision_Dec;
		Bool_t          muplus_Hlt1TrackMVADecision_TIS;
		Bool_t          muplus_Hlt1TrackMVADecision_TOS;
		Bool_t          muplus_Hlt1TwoTrackMVADecision_Dec;
		Bool_t          muplus_Hlt1TwoTrackMVADecision_TIS;
		Bool_t          muplus_Hlt1TwoTrackMVADecision_TOS;
		Bool_t          muplus_Hlt1SingleMuonHighPTDecision_Dec;
		Bool_t          muplus_Hlt1SingleMuonHighPTDecision_TIS;
		Bool_t          muplus_Hlt1SingleMuonHighPTDecision_TOS;
		Bool_t          muplus_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          muplus_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          muplus_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          muplus_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          muplus_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          muplus_Hlt2LowMultLMR2HHWSDecision_TOS;
		Bool_t          muplus_Hlt2SingleMuonHighPTDecision_Dec;
		Bool_t          muplus_Hlt2SingleMuonHighPTDecision_TIS;
		Bool_t          muplus_Hlt2SingleMuonHighPTDecision_TOS;
		Int_t           muplus_TRACK_Type;
		Int_t           muplus_TRACK_Key;
		Double_t        muplus_TRACK_CHI2NDOF;
		Double_t        muplus_TRACK_PCHI2;
		Double_t        muplus_TRACK_MatchCHI2;
		Double_t        muplus_TRACK_GhostProb;
		Double_t        muplus_TRACK_CloneDist;
		Double_t        muplus_TRACK_Likelihood;
		Double_t        muminus_MC12TuneV2_ProbNNe;
		Double_t        muminus_MC12TuneV2_ProbNNmu;
		Double_t        muminus_MC12TuneV2_ProbNNpi;
		Double_t        muminus_MC12TuneV2_ProbNNk;
		Double_t        muminus_MC12TuneV2_ProbNNp;
		Double_t        muminus_MC12TuneV2_ProbNNghost;
		Double_t        muminus_MC12TuneV3_ProbNNe;
		Double_t        muminus_MC12TuneV3_ProbNNmu;
		Double_t        muminus_MC12TuneV3_ProbNNpi;
		Double_t        muminus_MC12TuneV3_ProbNNk;
		Double_t        muminus_MC12TuneV3_ProbNNp;
		Double_t        muminus_MC12TuneV3_ProbNNghost;
		Double_t        muminus_MC12TuneV4_ProbNNe;
		Double_t        muminus_MC12TuneV4_ProbNNmu;
		Double_t        muminus_MC12TuneV4_ProbNNpi;
		Double_t        muminus_MC12TuneV4_ProbNNk;
		Double_t        muminus_MC12TuneV4_ProbNNp;
		Double_t        muminus_MC12TuneV4_ProbNNghost;
		Double_t        muminus_MC15TuneV1_ProbNNe;
		Double_t        muminus_MC15TuneV1_ProbNNmu;
		Double_t        muminus_MC15TuneV1_ProbNNpi;
		Double_t        muminus_MC15TuneV1_ProbNNk;
		Double_t        muminus_MC15TuneV1_ProbNNp;
		Double_t        muminus_MC15TuneV1_ProbNNghost;
		Double_t        muminus_OWNPV_X;
		Double_t        muminus_OWNPV_Y;
		Double_t        muminus_OWNPV_Z;
		Double_t        muminus_OWNPV_XERR;
		Double_t        muminus_OWNPV_YERR;
		Double_t        muminus_OWNPV_ZERR;
		Double_t        muminus_OWNPV_CHI2;
		Int_t           muminus_OWNPV_NDOF;
		Float_t         muminus_OWNPV_COV_[3][3];
		Double_t        muminus_IP_OWNPV;
		Double_t        muminus_IPCHI2_OWNPV;
		Double_t        muminus_ORIVX_X;
		Double_t        muminus_ORIVX_Y;
		Double_t        muminus_ORIVX_Z;
		Double_t        muminus_ORIVX_XERR;
		Double_t        muminus_ORIVX_YERR;
		Double_t        muminus_ORIVX_ZERR;
		Double_t        muminus_ORIVX_CHI2;
		Int_t           muminus_ORIVX_NDOF;
		Float_t         muminus_ORIVX_COV_[3][3];
		Double_t        muminus_P;
		Double_t        muminus_PT;
		Double_t        muminus_PE;
		Double_t        muminus_PX;
		Double_t        muminus_PY;
		Double_t        muminus_PZ;
		Double_t        muminus_M;
		Int_t           muminus_ID;
		Double_t        muminus_PIDe;
		Double_t        muminus_PIDmu;
		Double_t        muminus_PIDK;
		Double_t        muminus_PIDp;
		Double_t        muminus_PIDd;
		Double_t        muminus_ProbNNe;
		Double_t        muminus_ProbNNk;
		Double_t        muminus_ProbNNp;
		Double_t        muminus_ProbNNpi;
		Double_t        muminus_ProbNNmu;
		Double_t        muminus_ProbNNd;
		Double_t        muminus_ProbNNghost;
		Bool_t          muminus_hasMuon;
		Bool_t          muminus_isMuon;
		Bool_t          muminus_hasRich;
		Bool_t          muminus_UsedRichAerogel;
		Bool_t          muminus_UsedRich1Gas;
		Bool_t          muminus_UsedRich2Gas;
		Bool_t          muminus_RichAboveElThres;
		Bool_t          muminus_RichAboveMuThres;
		Bool_t          muminus_RichAbovePiThres;
		Bool_t          muminus_RichAboveKaThres;
		Bool_t          muminus_RichAbovePrThres;
		Bool_t          muminus_hasCalo;
		Double_t        muminus_PP_CombDLLe;
		Double_t        muminus_PP_CombDLLmu;
		Double_t        muminus_PP_CombDLLpi;
		Double_t        muminus_PP_CombDLLk;
		Double_t        muminus_PP_CombDLLp;
		Double_t        muminus_PP_CombDLLd;
		Double_t        muminus_PP_ProbNNe;
		Double_t        muminus_PP_ProbNNmu;
		Double_t        muminus_PP_ProbNNpi;
		Double_t        muminus_PP_ProbNNk;
		Double_t        muminus_PP_ProbNNp;
		Double_t        muminus_PP_ProbNNghost;
		Double_t        muminus_PP_ProbNNd;
		Bool_t          muminus_L0Global_Dec;
		Bool_t          muminus_L0Global_TIS;
		Bool_t          muminus_L0Global_TOS;
		Bool_t          muminus_Hlt1Global_Dec;
		Bool_t          muminus_Hlt1Global_TIS;
		Bool_t          muminus_Hlt1Global_TOS;
		Bool_t          muminus_Hlt1Phys_Dec;
		Bool_t          muminus_Hlt1Phys_TIS;
		Bool_t          muminus_Hlt1Phys_TOS;
		Bool_t          muminus_Hlt2Global_Dec;
		Bool_t          muminus_Hlt2Global_TIS;
		Bool_t          muminus_Hlt2Global_TOS;
		Bool_t          muminus_Hlt2Phys_Dec;
		Bool_t          muminus_Hlt2Phys_TIS;
		Bool_t          muminus_Hlt2Phys_TOS;
		Bool_t          muminus_L0DiHadron_lowMultDecision_Dec;
		Bool_t          muminus_L0DiHadron_lowMultDecision_TIS;
		Bool_t          muminus_L0DiHadron_lowMultDecision_TOS;
		Bool_t          muminus_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          muminus_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          muminus_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          muminus_L0MuonDecision_Dec;
		Bool_t          muminus_L0MuonDecision_TIS;
		Bool_t          muminus_L0MuonDecision_TOS;
		Bool_t          muminus_L0MuonEWDecision_Dec;
		Bool_t          muminus_L0MuonEWDecision_TIS;
		Bool_t          muminus_L0MuonEWDecision_TOS;
		Bool_t          muminus_L0DiMuonDecision_Dec;
		Bool_t          muminus_L0DiMuonDecision_TIS;
		Bool_t          muminus_L0DiMuonDecision_TOS;
		Bool_t          muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          muminus_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          muminus_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          muminus_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          muminus_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          muminus_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          muminus_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          muminus_Hlt1TrackMVADecision_Dec;
		Bool_t          muminus_Hlt1TrackMVADecision_TIS;
		Bool_t          muminus_Hlt1TrackMVADecision_TOS;
		Bool_t          muminus_Hlt1TwoTrackMVADecision_Dec;
		Bool_t          muminus_Hlt1TwoTrackMVADecision_TIS;
		Bool_t          muminus_Hlt1TwoTrackMVADecision_TOS;
		Bool_t          muminus_Hlt1SingleMuonHighPTDecision_Dec;
		Bool_t          muminus_Hlt1SingleMuonHighPTDecision_TIS;
		Bool_t          muminus_Hlt1SingleMuonHighPTDecision_TOS;
		Bool_t          muminus_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          muminus_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          muminus_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          muminus_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          muminus_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          muminus_Hlt2LowMultLMR2HHWSDecision_TOS;
		Bool_t          muminus_Hlt2SingleMuonHighPTDecision_Dec;
		Bool_t          muminus_Hlt2SingleMuonHighPTDecision_TIS;
		Bool_t          muminus_Hlt2SingleMuonHighPTDecision_TOS;
		Int_t           muminus_TRACK_Type;
		Int_t           muminus_TRACK_Key;
		Double_t        muminus_TRACK_CHI2NDOF;
		Double_t        muminus_TRACK_PCHI2;
		Double_t        muminus_TRACK_MatchCHI2;
		Double_t        muminus_TRACK_GhostProb;
		Double_t        muminus_TRACK_CloneDist;
		Double_t        muminus_TRACK_Likelihood;
		UInt_t          nCandidate;
		ULong64_t       totCandidates;
		ULong64_t       EventInSequence;
		Double_t        LoKi_nCharged;
		Double_t        LoKi_nITClusters;
		Double_t        LoKi_nNeutrals;
		Double_t        LoKi_nOThits;
		Double_t        LoKi_nPVs;
		Double_t        LoKi_nSpdMult;
		Double_t        LoKi_nTTClusters;
		Double_t        LoKi_nVeloClusters;
		Double_t        LoKi_nVeloLiteClusters;
		UInt_t          runNumber;
		ULong64_t       eventNumber;
		UInt_t          BCID;
		Int_t           BCType;
		UInt_t          OdinTCK;
		UInt_t          L0DUTCK;
		UInt_t          HLT1TCK;
		UInt_t          HLT2TCK;
		ULong64_t       GpsTime;
		Short_t         Polarity;
		Int_t           B00;
		Int_t           B01;
		Int_t           B02;
		Int_t           B03;
		Int_t           B10;
		Int_t           B11;
		Int_t           B12;
		Int_t           B13;
		Int_t           B20;
		Int_t           B21;
		Int_t           B22;
		Int_t           B23;
		Int_t           F10;
		Int_t           F11;
		Int_t           F12;
		Int_t           F13;
		Int_t           F20;
		Int_t           F21;
		Int_t           F22;
		Int_t           F23;
		Double_t        log_hrc_fom_v3;
		Double_t        log_hrc_fom_B_v3;
		Double_t        log_hrc_fom_F_v3;
		Int_t           nchB;
		Float_t         adc_B[1000];   //[nchB]
		Int_t           nchF;
		Float_t         adc_F[1000];   //[nchF]
		Int_t           nPV;
		Float_t         PVX[100];   //[nPV]
		Float_t         PVY[100];   //[nPV]
		Float_t         PVZ[100];   //[nPV]
		Float_t         PVXERR[100];   //[nPV]
		Float_t         PVYERR[100];   //[nPV]
		Float_t         PVZERR[100];   //[nPV]
		Float_t         PVCHI2[100];   //[nPV]
		Float_t         PVNDOF[100];   //[nPV]
		Float_t         PVNTRACKS[100];   //[nPV]
		Int_t           nPVs;
		Int_t           nTracks;
		Int_t           nLongTracks;
		Int_t           nDownstreamTracks;
		Int_t           nUpstreamTracks;
		Int_t           nVeloTracks;
		Int_t           nTTracks;
		Int_t           nBackTracks;
		Int_t           nRich1Hits;
		Int_t           nRich2Hits;
		Int_t           nVeloClusters;
		Int_t           nITClusters;
		Int_t           nTTClusters;
		Int_t           nOTClusters;
		Int_t           nSPDHits;
		Int_t           nMuonCoordsS0;
		Int_t           nMuonCoordsS1;
		Int_t           nMuonCoordsS2;
		Int_t           nMuonCoordsS3;
		Int_t           nMuonCoordsS4;
		Int_t           nMuonTracks;

		// List of branches
		TBranch        *b_phi_ENDVERTEX_X;   //!
		TBranch        *b_phi_ENDVERTEX_Y;   //!
		TBranch        *b_phi_ENDVERTEX_Z;   //!
		TBranch        *b_phi_ENDVERTEX_XERR;   //!
		TBranch        *b_phi_ENDVERTEX_YERR;   //!
		TBranch        *b_phi_ENDVERTEX_ZERR;   //!
		TBranch        *b_phi_ENDVERTEX_CHI2;   //!
		TBranch        *b_phi_ENDVERTEX_NDOF;   //!
		TBranch        *b_phi_ENDVERTEX_COV_;   //!
		TBranch        *b_phi_OWNPV_X;   //!
		TBranch        *b_phi_OWNPV_Y;   //!
		TBranch        *b_phi_OWNPV_Z;   //!
		TBranch        *b_phi_OWNPV_XERR;   //!
		TBranch        *b_phi_OWNPV_YERR;   //!
		TBranch        *b_phi_OWNPV_ZERR;   //!
		TBranch        *b_phi_OWNPV_CHI2;   //!
		TBranch        *b_phi_OWNPV_NDOF;   //!
		TBranch        *b_phi_OWNPV_COV_;   //!
		TBranch        *b_phi_IP_OWNPV;   //!
		TBranch        *b_phi_IPCHI2_OWNPV;   //!
		TBranch        *b_phi_FD_OWNPV;   //!
		TBranch        *b_phi_FDCHI2_OWNPV;   //!
		TBranch        *b_phi_DIRA_OWNPV;   //!
		TBranch        *b_phi_P;   //!
		TBranch        *b_phi_PT;   //!
		TBranch        *b_phi_PE;   //!
		TBranch        *b_phi_PX;   //!
		TBranch        *b_phi_PY;   //!
		TBranch        *b_phi_PZ;   //!
		TBranch        *b_phi_MM;   //!
		TBranch        *b_phi_MMERR;   //!
		TBranch        *b_phi_M;   //!
		TBranch        *b_phi_ID;   //!
		TBranch        *b_phi_L0Global_Dec;   //!
		TBranch        *b_phi_L0Global_TIS;   //!
		TBranch        *b_phi_L0Global_TOS;   //!
		TBranch        *b_phi_Hlt1Global_Dec;   //!
		TBranch        *b_phi_Hlt1Global_TIS;   //!
		TBranch        *b_phi_Hlt1Global_TOS;   //!
		TBranch        *b_phi_Hlt1Phys_Dec;   //!
		TBranch        *b_phi_Hlt1Phys_TIS;   //!
		TBranch        *b_phi_Hlt1Phys_TOS;   //!
		TBranch        *b_phi_Hlt2Global_Dec;   //!
		TBranch        *b_phi_Hlt2Global_TIS;   //!
		TBranch        *b_phi_Hlt2Global_TOS;   //!
		TBranch        *b_phi_Hlt2Phys_Dec;   //!
		TBranch        *b_phi_Hlt2Phys_TIS;   //!
		TBranch        *b_phi_Hlt2Phys_TOS;   //!
		//TBranch        *b_phi_L0DiHadron_lowMultDecision_Dec;   //!
		//TBranch        *b_phi_L0DiHadron_lowMultDecision_TIS;   //!
		//TBranch        *b_phi_L0DiHadron_lowMultDecision_TOS;   //!
		//TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_Dec;   //!
		//TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_TIS;   //!
		//TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_TOS;   //!
		//TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		//TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		//TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		//TBranch        *b_phi_Hlt1LowMultHerschelDecision_Dec;   //!
		//TBranch        *b_phi_Hlt1LowMultHerschelDecision_TIS;   //!
		//TBranch        *b_phi_Hlt1LowMultHerschelDecision_TOS;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		//TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_Dec;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_TIS;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_TOS;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		//TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_phi_L0DiHadron_lowMultDecision_Dec;
		TBranch        *b_phi_L0Muon_lowMultDecision_Dec;
		TBranch        *b_phi_L0DiMuon_lowMultDecision_Dec;
		TBranch        *b_phi_L0Electron_lowMultDecision_Dec;
		TBranch        *b_phi_L0Photon_lowMultDecision_Dec;
		TBranch        *b_phi_L0DiEM_lowMultDecision_Dec;
		TBranch        *b_phi_L0SPDDecision_Dec;
		TBranch        *b_phi_L0CALODecision_Dec;
		TBranch        *b_phi_L0PUDecision_Dec;
		TBranch        *b_phi_L0MuonDecision_Dec;
		TBranch        *b_phi_Hlt1LowMultPassThroughDecision_Dec;
		TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		TBranch        *b_phi_Hlt1BBMicroBiasVeloDecision_Dec;
		TBranch        *b_phi_Hlt1BBHighMultDecision_Dec;
		TBranch        *b_phi_Hlt1MBNoBiasDecision_Dec;
		TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_Dec;
		TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_Dec;
		TBranch        *b_phi_Hlt2MBMicroBiasVeloDecision_Dec;
		TBranch        *b_phi_Hlt2MBHighMultDecision_Dec;
		TBranch        *b_phi_Hlt2MBNoBiasDecision_Dec;
		TBranch        *b_phi_Hlt2PassThroughDecision_Dec;
		TBranch        *b_phi_L0HadronDecision_Dec;
//		TBranch        *b_phi_L0MuonDecision_Dec;
//		TBranch        *b_phi_L0SPDDecision_Dec;
		TBranch        *b_phi_L0HadronLowMultDecision_Dec;
		TBranch        *b_phi_L0MuonLowMultDecision_Dec;
		TBranch        *b_phi_L0ElectronLowMultDecision_Dec;
		TBranch        *b_phi_L0PhotonLowMultDecision_Dec;
		TBranch        *b_phi_L0DiEMLowMultDecision_Dec;
		TBranch        *b_phi_L0SPDLowMultDecision_Dec;
		TBranch        *b_phi_L0SoftCEPDecision_Dec;
		TBranch        *b_phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec;
		TBranch        *b_phi_Hlt1BBMicroBiasSoftCEPDecision_Dec;
		TBranch        *b_phi_Hlt1BBHasTrackDecision_Dec;
		TBranch        *b_phi_Hlt2BBPassThroughDecision_Dec;
		TBranch        *b_phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec;
		TBranch        *b_phi_Hlt2BBLongTrackDecision_Dec;
		TBranch        *b_phi_Hlt2SingleTrackDecision_Dec;
//		TBranch        *b_phi_Hlt2MBMicroBiasVeloDecision_Dec;
		TBranch        *b_Kplus_MC12TuneV2_ProbNNe;   //!
		TBranch        *b_Kplus_MC12TuneV2_ProbNNmu;   //!
		TBranch        *b_Kplus_MC12TuneV2_ProbNNpi;   //!
		TBranch        *b_Kplus_MC12TuneV2_ProbNNk;   //!
		TBranch        *b_Kplus_MC12TuneV2_ProbNNp;   //!
		TBranch        *b_Kplus_MC12TuneV2_ProbNNghost;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNe;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNmu;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNpi;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNk;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNp;   //!
		TBranch        *b_Kplus_MC12TuneV3_ProbNNghost;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNe;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNmu;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNpi;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNk;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNp;   //!
		TBranch        *b_Kplus_MC12TuneV4_ProbNNghost;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNe;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNmu;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNpi;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNk;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNp;   //!
		TBranch        *b_Kplus_MC15TuneV1_ProbNNghost;   //!
		TBranch        *b_Kplus_OWNPV_X;   //!
		TBranch        *b_Kplus_OWNPV_Y;   //!
		TBranch        *b_Kplus_OWNPV_Z;   //!
		TBranch        *b_Kplus_OWNPV_XERR;   //!
		TBranch        *b_Kplus_OWNPV_YERR;   //!
		TBranch        *b_Kplus_OWNPV_ZERR;   //!
		TBranch        *b_Kplus_OWNPV_CHI2;   //!
		TBranch        *b_Kplus_OWNPV_NDOF;   //!
		TBranch        *b_Kplus_OWNPV_COV_;   //!
		TBranch        *b_Kplus_IP_OWNPV;   //!
		TBranch        *b_Kplus_IPCHI2_OWNPV;   //!
		TBranch        *b_Kplus_ORIVX_X;   //!
		TBranch        *b_Kplus_ORIVX_Y;   //!
		TBranch        *b_Kplus_ORIVX_Z;   //!
		TBranch        *b_Kplus_ORIVX_XERR;   //!
		TBranch        *b_Kplus_ORIVX_YERR;   //!
		TBranch        *b_Kplus_ORIVX_ZERR;   //!
		TBranch        *b_Kplus_ORIVX_CHI2;   //!
		TBranch        *b_Kplus_ORIVX_NDOF;   //!
		TBranch        *b_Kplus_ORIVX_COV_;   //!
		TBranch        *b_Kplus_P;   //!
		TBranch        *b_Kplus_PT;   //!
		TBranch        *b_Kplus_PE;   //!
		TBranch        *b_Kplus_PX;   //!
		TBranch        *b_Kplus_PY;   //!
		TBranch        *b_Kplus_PZ;   //!
		TBranch        *b_Kplus_M;   //!
		TBranch        *b_Kplus_ID;   //!
		TBranch        *b_Kplus_PIDe;   //!
		TBranch        *b_Kplus_PIDmu;   //!
		TBranch        *b_Kplus_PIDK;   //!
		TBranch        *b_Kplus_PIDp;   //!
		TBranch        *b_Kplus_ProbNNe;   //!
		TBranch        *b_Kplus_ProbNNk;   //!
		TBranch        *b_Kplus_ProbNNp;   //!
		TBranch        *b_Kplus_ProbNNpi;   //!
		TBranch        *b_Kplus_ProbNNmu;   //!
		TBranch        *b_Kplus_ProbNNghost;   //!
		TBranch        *b_Kplus_hasMuon;   //!
		TBranch        *b_Kplus_isMuon;   //!
		TBranch        *b_Kplus_hasRich;   //!
		TBranch        *b_Kplus_UsedRichAerogel;   //!
		TBranch        *b_Kplus_UsedRich1Gas;   //!
		TBranch        *b_Kplus_UsedRich2Gas;   //!
		TBranch        *b_Kplus_RichAboveElThres;   //!
		TBranch        *b_Kplus_RichAboveMuThres;   //!
		TBranch        *b_Kplus_RichAbovePiThres;   //!
		TBranch        *b_Kplus_RichAboveKaThres;   //!
		TBranch        *b_Kplus_RichAbovePrThres;   //!
		TBranch        *b_Kplus_hasCalo;   //!
		TBranch        *b_Kplus_PP_CombDLLe;   //!
		TBranch        *b_Kplus_PP_CombDLLmu;   //!
		TBranch        *b_Kplus_PP_CombDLLpi;   //!
		TBranch        *b_Kplus_PP_CombDLLk;   //!
		TBranch        *b_Kplus_PP_CombDLLp;   //!
		TBranch        *b_Kplus_PP_CombDLLd;   //!
		TBranch        *b_Kplus_PP_ProbNNe;   //!
		TBranch        *b_Kplus_PP_ProbNNmu;   //!
		TBranch        *b_Kplus_PP_ProbNNpi;   //!
		TBranch        *b_Kplus_PP_ProbNNk;   //!
		TBranch        *b_Kplus_PP_ProbNNp;   //!
		TBranch        *b_Kplus_PP_ProbNNghost;   //!
		TBranch        *b_Kplus_PP_ProbNNd;   //!
		TBranch        *b_Kplus_L0Global_Dec;   //!
		TBranch        *b_Kplus_L0Global_TIS;   //!
		TBranch        *b_Kplus_L0Global_TOS;   //!
		TBranch        *b_Kplus_Hlt1Global_Dec;   //!
		TBranch        *b_Kplus_Hlt1Global_TIS;   //!
		TBranch        *b_Kplus_Hlt1Global_TOS;   //!
		TBranch        *b_Kplus_Hlt1Phys_Dec;   //!
		TBranch        *b_Kplus_Hlt1Phys_TIS;   //!
		TBranch        *b_Kplus_Hlt1Phys_TOS;   //!
		TBranch        *b_Kplus_Hlt2Global_Dec;   //!
		TBranch        *b_Kplus_Hlt2Global_TIS;   //!
		TBranch        *b_Kplus_Hlt2Global_TOS;   //!
		TBranch        *b_Kplus_Hlt2Phys_Dec;   //!
		TBranch        *b_Kplus_Hlt2Phys_TIS;   //!
		TBranch        *b_Kplus_Hlt2Phys_TOS;   //!
		TBranch        *b_Kplus_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Kplus_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Kplus_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Kplus_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Kplus_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Kplus_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_Kplus_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_Kplus_TRACK_Type;   //!
		TBranch        *b_Kplus_TRACK_Key;   //!
		TBranch        *b_Kplus_TRACK_CHI2NDOF;   //!
		TBranch        *b_Kplus_TRACK_PCHI2;   //!
		TBranch        *b_Kplus_TRACK_MatchCHI2;   //!
		TBranch        *b_Kplus_TRACK_GhostProb;   //!
		TBranch        *b_Kplus_TRACK_CloneDist;   //!
		TBranch        *b_Kplus_TRACK_Likelihood;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNe;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNmu;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNpi;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNk;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNp;   //!
		TBranch        *b_Kminus_MC12TuneV2_ProbNNghost;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNe;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNmu;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNpi;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNk;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNp;   //!
		TBranch        *b_Kminus_MC12TuneV3_ProbNNghost;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNe;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNmu;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNpi;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNk;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNp;   //!
		TBranch        *b_Kminus_MC12TuneV4_ProbNNghost;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNe;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNmu;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNpi;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNk;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNp;   //!
		TBranch        *b_Kminus_MC15TuneV1_ProbNNghost;   //!
		TBranch        *b_Kminus_OWNPV_X;   //!
		TBranch        *b_Kminus_OWNPV_Y;   //!
		TBranch        *b_Kminus_OWNPV_Z;   //!
		TBranch        *b_Kminus_OWNPV_XERR;   //!
		TBranch        *b_Kminus_OWNPV_YERR;   //!
		TBranch        *b_Kminus_OWNPV_ZERR;   //!
		TBranch        *b_Kminus_OWNPV_CHI2;   //!
		TBranch        *b_Kminus_OWNPV_NDOF;   //!
		TBranch        *b_Kminus_OWNPV_COV_;   //!
		TBranch        *b_Kminus_IP_OWNPV;   //!
		TBranch        *b_Kminus_IPCHI2_OWNPV;   //!
		TBranch        *b_Kminus_ORIVX_X;   //!
		TBranch        *b_Kminus_ORIVX_Y;   //!
		TBranch        *b_Kminus_ORIVX_Z;   //!
		TBranch        *b_Kminus_ORIVX_XERR;   //!
		TBranch        *b_Kminus_ORIVX_YERR;   //!
		TBranch        *b_Kminus_ORIVX_ZERR;   //!
		TBranch        *b_Kminus_ORIVX_CHI2;   //!
		TBranch        *b_Kminus_ORIVX_NDOF;   //!
		TBranch        *b_Kminus_ORIVX_COV_;   //!
		TBranch        *b_Kminus_P;   //!
		TBranch        *b_Kminus_PT;   //!
		TBranch        *b_Kminus_PE;   //!
		TBranch        *b_Kminus_PX;   //!
		TBranch        *b_Kminus_PY;   //!
		TBranch        *b_Kminus_PZ;   //!
		TBranch        *b_Kminus_M;   //!
		TBranch        *b_Kminus_ID;   //!
		TBranch        *b_Kminus_PIDe;   //!
		TBranch        *b_Kminus_PIDmu;   //!
		TBranch        *b_Kminus_PIDK;   //!
		TBranch        *b_Kminus_PIDp;   //!
		TBranch        *b_Kminus_ProbNNe;   //!
		TBranch        *b_Kminus_ProbNNk;   //!
		TBranch        *b_Kminus_ProbNNp;   //!
		TBranch        *b_Kminus_ProbNNpi;   //!
		TBranch        *b_Kminus_ProbNNmu;   //!
		TBranch        *b_Kminus_ProbNNghost;   //!
		TBranch        *b_Kminus_hasMuon;   //!
		TBranch        *b_Kminus_isMuon;   //!
		TBranch        *b_Kminus_hasRich;   //!
		TBranch        *b_Kminus_UsedRichAerogel;   //!
		TBranch        *b_Kminus_UsedRich1Gas;   //!
		TBranch        *b_Kminus_UsedRich2Gas;   //!
		TBranch        *b_Kminus_RichAboveElThres;   //!
		TBranch        *b_Kminus_RichAboveMuThres;   //!
		TBranch        *b_Kminus_RichAbovePiThres;   //!
		TBranch        *b_Kminus_RichAboveKaThres;   //!
		TBranch        *b_Kminus_RichAbovePrThres;   //!
		TBranch        *b_Kminus_hasCalo;   //!
		TBranch        *b_Kminus_PP_CombDLLe;   //!
		TBranch        *b_Kminus_PP_CombDLLmu;   //!
		TBranch        *b_Kminus_PP_CombDLLpi;   //!
		TBranch        *b_Kminus_PP_CombDLLk;   //!
		TBranch        *b_Kminus_PP_CombDLLp;   //!
		TBranch        *b_Kminus_PP_CombDLLd;   //!
		TBranch        *b_Kminus_PP_ProbNNe;   //!
		TBranch        *b_Kminus_PP_ProbNNmu;   //!
		TBranch        *b_Kminus_PP_ProbNNpi;   //!
		TBranch        *b_Kminus_PP_ProbNNk;   //!
		TBranch        *b_Kminus_PP_ProbNNp;   //!
		TBranch        *b_Kminus_PP_ProbNNghost;   //!
		TBranch        *b_Kminus_PP_ProbNNd;   //!
		TBranch        *b_Kminus_L0Global_Dec;   //!
		TBranch        *b_Kminus_L0Global_TIS;   //!
		TBranch        *b_Kminus_L0Global_TOS;   //!
		TBranch        *b_Kminus_Hlt1Global_Dec;   //!
		TBranch        *b_Kminus_Hlt1Global_TIS;   //!
		TBranch        *b_Kminus_Hlt1Global_TOS;   //!
		TBranch        *b_Kminus_Hlt1Phys_Dec;   //!
		TBranch        *b_Kminus_Hlt1Phys_TIS;   //!
		TBranch        *b_Kminus_Hlt1Phys_TOS;   //!
		TBranch        *b_Kminus_Hlt2Global_Dec;   //!
		TBranch        *b_Kminus_Hlt2Global_TIS;   //!
		TBranch        *b_Kminus_Hlt2Global_TOS;   //!
		TBranch        *b_Kminus_Hlt2Phys_Dec;   //!
		TBranch        *b_Kminus_Hlt2Phys_TIS;   //!
		TBranch        *b_Kminus_Hlt2Phys_TOS;   //!
		TBranch        *b_Kminus_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Kminus_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Kminus_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Kminus_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Kminus_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Kminus_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_Kminus_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_Kminus_TRACK_Type;   //!
		TBranch        *b_Kminus_TRACK_Key;   //!
		TBranch        *b_Kminus_TRACK_CHI2NDOF;   //!
		TBranch        *b_Kminus_TRACK_PCHI2;   //!
		TBranch        *b_Kminus_TRACK_MatchCHI2;   //!
		TBranch        *b_Kminus_TRACK_GhostProb;   //!
		TBranch        *b_Kminus_TRACK_CloneDist;   //!
		TBranch        *b_Kminus_TRACK_Likelihood;   //!
		TBranch        *b_Z_ENDVERTEX_X;   //!
		TBranch        *b_Z_ENDVERTEX_Y;   //!
		TBranch        *b_Z_ENDVERTEX_Z;   //!
		TBranch        *b_Z_ENDVERTEX_XERR;   //!
		TBranch        *b_Z_ENDVERTEX_YERR;   //!
		TBranch        *b_Z_ENDVERTEX_ZERR;   //!
		TBranch        *b_Z_ENDVERTEX_CHI2;   //!
		TBranch        *b_Z_ENDVERTEX_NDOF;   //!
		TBranch        *b_Z_ENDVERTEX_COV_;   //!
		TBranch        *b_Z_OWNPV_X;   //!
		TBranch        *b_Z_OWNPV_Y;   //!
		TBranch        *b_Z_OWNPV_Z;   //!
		TBranch        *b_Z_OWNPV_XERR;   //!
		TBranch        *b_Z_OWNPV_YERR;   //!
		TBranch        *b_Z_OWNPV_ZERR;   //!
		TBranch        *b_Z_OWNPV_CHI2;   //!
		TBranch        *b_Z_OWNPV_NDOF;   //!
		TBranch        *b_Z_OWNPV_COV_;   //!
		TBranch        *b_Z_IP_OWNPV;   //!
		TBranch        *b_Z_IPCHI2_OWNPV;   //!
		TBranch        *b_Z_FD_OWNPV;   //!
		TBranch        *b_Z_FDCHI2_OWNPV;   //!
		TBranch        *b_Z_DIRA_OWNPV;   //!
		TBranch        *b_Z_P;   //!
		TBranch        *b_Z_PT;   //!
		TBranch        *b_Z_PE;   //!
		TBranch        *b_Z_PX;   //!
		TBranch        *b_Z_PY;   //!
		TBranch        *b_Z_PZ;   //!
		TBranch        *b_Z_MM;   //!
		TBranch        *b_Z_MMERR;   //!
		TBranch        *b_Z_M;   //!
		TBranch        *b_Z_ID;   //!
		TBranch        *b_Z_L0Global_Dec;   //!
		TBranch        *b_Z_L0Global_TIS;   //!
		TBranch        *b_Z_L0Global_TOS;   //!
		TBranch        *b_Z_Hlt1Global_Dec;   //!
		TBranch        *b_Z_Hlt1Global_TIS;   //!
		TBranch        *b_Z_Hlt1Global_TOS;   //!
		TBranch        *b_Z_Hlt1Phys_Dec;   //!
		TBranch        *b_Z_Hlt1Phys_TIS;   //!
		TBranch        *b_Z_Hlt1Phys_TOS;   //!
		TBranch        *b_Z_Hlt2Global_Dec;   //!
		TBranch        *b_Z_Hlt2Global_TIS;   //!
		TBranch        *b_Z_Hlt2Global_TOS;   //!
		TBranch        *b_Z_Hlt2Phys_Dec;   //!
		TBranch        *b_Z_Hlt2Phys_TIS;   //!
		TBranch        *b_Z_Hlt2Phys_TOS;   //!
		TBranch        *b_Z_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Z_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Z_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Z_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_Z_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_Z_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_Z_L0MuonDecision_Dec;   //!
		TBranch        *b_Z_L0MuonDecision_TIS;   //!
		TBranch        *b_Z_L0MuonDecision_TOS;   //!
		TBranch        *b_Z_L0MuonEWDecision_Dec;   //!
		TBranch        *b_Z_L0MuonEWDecision_TIS;   //!
		TBranch        *b_Z_L0MuonEWDecision_TOS;   //!
		TBranch        *b_Z_L0DiMuonDecision_Dec;   //!
		TBranch        *b_Z_L0DiMuonDecision_TIS;   //!
		TBranch        *b_Z_L0DiMuonDecision_TOS;   //!
		TBranch        *b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_Z_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_Z_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_Z_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_Z_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_Z_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_Z_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_Z_Hlt1TrackMVADecision_Dec;   //!
		TBranch        *b_Z_Hlt1TrackMVADecision_TIS;   //!
		TBranch        *b_Z_Hlt1TrackMVADecision_TOS;   //!
		TBranch        *b_Z_Hlt1TwoTrackMVADecision_Dec;   //!
		TBranch        *b_Z_Hlt1TwoTrackMVADecision_TIS;   //!
		TBranch        *b_Z_Hlt1TwoTrackMVADecision_TOS;   //!
		TBranch        *b_Z_Hlt1SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_Z_Hlt1SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_Z_Hlt1SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_Z_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_Z_Hlt2SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_Z_Hlt2SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_Z_Hlt2SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNe;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNmu;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNpi;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNk;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNp;   //!
		TBranch        *b_muplus_MC12TuneV2_ProbNNghost;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNe;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNmu;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNpi;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNk;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNp;   //!
		TBranch        *b_muplus_MC12TuneV3_ProbNNghost;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNe;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNmu;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNpi;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNk;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNp;   //!
		TBranch        *b_muplus_MC12TuneV4_ProbNNghost;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNe;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNmu;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNpi;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNk;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNp;   //!
		TBranch        *b_muplus_MC15TuneV1_ProbNNghost;   //!
		TBranch        *b_muplus_OWNPV_X;   //!
		TBranch        *b_muplus_OWNPV_Y;   //!
		TBranch        *b_muplus_OWNPV_Z;   //!
		TBranch        *b_muplus_OWNPV_XERR;   //!
		TBranch        *b_muplus_OWNPV_YERR;   //!
		TBranch        *b_muplus_OWNPV_ZERR;   //!
		TBranch        *b_muplus_OWNPV_CHI2;   //!
		TBranch        *b_muplus_OWNPV_NDOF;   //!
		TBranch        *b_muplus_OWNPV_COV_;   //!
		TBranch        *b_muplus_IP_OWNPV;   //!
		TBranch        *b_muplus_IPCHI2_OWNPV;   //!
		TBranch        *b_muplus_ORIVX_X;   //!
		TBranch        *b_muplus_ORIVX_Y;   //!
		TBranch        *b_muplus_ORIVX_Z;   //!
		TBranch        *b_muplus_ORIVX_XERR;   //!
		TBranch        *b_muplus_ORIVX_YERR;   //!
		TBranch        *b_muplus_ORIVX_ZERR;   //!
		TBranch        *b_muplus_ORIVX_CHI2;   //!
		TBranch        *b_muplus_ORIVX_NDOF;   //!
		TBranch        *b_muplus_ORIVX_COV_;   //!
		TBranch        *b_muplus_P;   //!
		TBranch        *b_muplus_PT;   //!
		TBranch        *b_muplus_PE;   //!
		TBranch        *b_muplus_PX;   //!
		TBranch        *b_muplus_PY;   //!
		TBranch        *b_muplus_PZ;   //!
		TBranch        *b_muplus_M;   //!
		TBranch        *b_muplus_ID;   //!
		TBranch        *b_muplus_PIDe;   //!
		TBranch        *b_muplus_PIDmu;   //!
		TBranch        *b_muplus_PIDK;   //!
		TBranch        *b_muplus_PIDp;   //!
		TBranch        *b_muplus_PIDd;   //!
		TBranch        *b_muplus_ProbNNe;   //!
		TBranch        *b_muplus_ProbNNk;   //!
		TBranch        *b_muplus_ProbNNp;   //!
		TBranch        *b_muplus_ProbNNpi;   //!
		TBranch        *b_muplus_ProbNNmu;   //!
		TBranch        *b_muplus_ProbNNd;   //!
		TBranch        *b_muplus_ProbNNghost;   //!
		TBranch        *b_muplus_hasMuon;   //!
		TBranch        *b_muplus_isMuon;   //!
		TBranch        *b_muplus_hasRich;   //!
		TBranch        *b_muplus_UsedRichAerogel;   //!
		TBranch        *b_muplus_UsedRich1Gas;   //!
		TBranch        *b_muplus_UsedRich2Gas;   //!
		TBranch        *b_muplus_RichAboveElThres;   //!
		TBranch        *b_muplus_RichAboveMuThres;   //!
		TBranch        *b_muplus_RichAbovePiThres;   //!
		TBranch        *b_muplus_RichAboveKaThres;   //!
		TBranch        *b_muplus_RichAbovePrThres;   //!
		TBranch        *b_muplus_hasCalo;   //!
		TBranch        *b_muplus_PP_CombDLLe;   //!
		TBranch        *b_muplus_PP_CombDLLmu;   //!
		TBranch        *b_muplus_PP_CombDLLpi;   //!
		TBranch        *b_muplus_PP_CombDLLk;   //!
		TBranch        *b_muplus_PP_CombDLLp;   //!
		TBranch        *b_muplus_PP_CombDLLd;   //!
		TBranch        *b_muplus_PP_ProbNNe;   //!
		TBranch        *b_muplus_PP_ProbNNmu;   //!
		TBranch        *b_muplus_PP_ProbNNpi;   //!
		TBranch        *b_muplus_PP_ProbNNk;   //!
		TBranch        *b_muplus_PP_ProbNNp;   //!
		TBranch        *b_muplus_PP_ProbNNghost;   //!
		TBranch        *b_muplus_PP_ProbNNd;   //!
		TBranch        *b_muplus_L0Global_Dec;   //!
		TBranch        *b_muplus_L0Global_TIS;   //!
		TBranch        *b_muplus_L0Global_TOS;   //!
		TBranch        *b_muplus_Hlt1Global_Dec;   //!
		TBranch        *b_muplus_Hlt1Global_TIS;   //!
		TBranch        *b_muplus_Hlt1Global_TOS;   //!
		TBranch        *b_muplus_Hlt1Phys_Dec;   //!
		TBranch        *b_muplus_Hlt1Phys_TIS;   //!
		TBranch        *b_muplus_Hlt1Phys_TOS;   //!
		TBranch        *b_muplus_Hlt2Global_Dec;   //!
		TBranch        *b_muplus_Hlt2Global_TIS;   //!
		TBranch        *b_muplus_Hlt2Global_TOS;   //!
		TBranch        *b_muplus_Hlt2Phys_Dec;   //!
		TBranch        *b_muplus_Hlt2Phys_TIS;   //!
		TBranch        *b_muplus_Hlt2Phys_TOS;   //!
		TBranch        *b_muplus_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_muplus_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_muplus_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_muplus_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_muplus_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_muplus_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_muplus_L0MuonDecision_Dec;   //!
		TBranch        *b_muplus_L0MuonDecision_TIS;   //!
		TBranch        *b_muplus_L0MuonDecision_TOS;   //!
		TBranch        *b_muplus_L0MuonEWDecision_Dec;   //!
		TBranch        *b_muplus_L0MuonEWDecision_TIS;   //!
		TBranch        *b_muplus_L0MuonEWDecision_TOS;   //!
		TBranch        *b_muplus_L0DiMuonDecision_Dec;   //!
		TBranch        *b_muplus_L0DiMuonDecision_TIS;   //!
		TBranch        *b_muplus_L0DiMuonDecision_TOS;   //!
		TBranch        *b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_muplus_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_muplus_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_muplus_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_muplus_Hlt1TrackMVADecision_Dec;   //!
		TBranch        *b_muplus_Hlt1TrackMVADecision_TIS;   //!
		TBranch        *b_muplus_Hlt1TrackMVADecision_TOS;   //!
		TBranch        *b_muplus_Hlt1TwoTrackMVADecision_Dec;   //!
		TBranch        *b_muplus_Hlt1TwoTrackMVADecision_TIS;   //!
		TBranch        *b_muplus_Hlt1TwoTrackMVADecision_TOS;   //!
		TBranch        *b_muplus_Hlt1SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_muplus_Hlt1SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_muplus_Hlt1SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_muplus_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_muplus_Hlt2SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_muplus_Hlt2SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_muplus_Hlt2SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_muplus_TRACK_Type;   //!
		TBranch        *b_muplus_TRACK_Key;   //!
		TBranch        *b_muplus_TRACK_CHI2NDOF;   //!
		TBranch        *b_muplus_TRACK_PCHI2;   //!
		TBranch        *b_muplus_TRACK_MatchCHI2;   //!
		TBranch        *b_muplus_TRACK_GhostProb;   //!
		TBranch        *b_muplus_TRACK_CloneDist;   //!
		TBranch        *b_muplus_TRACK_Likelihood;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNe;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNmu;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNpi;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNk;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNp;   //!
		TBranch        *b_muminus_MC12TuneV2_ProbNNghost;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNe;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNmu;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNpi;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNk;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNp;   //!
		TBranch        *b_muminus_MC12TuneV3_ProbNNghost;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNe;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNmu;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNpi;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNk;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNp;   //!
		TBranch        *b_muminus_MC12TuneV4_ProbNNghost;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNe;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNmu;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNpi;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNk;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNp;   //!
		TBranch        *b_muminus_MC15TuneV1_ProbNNghost;   //!
		TBranch        *b_muminus_OWNPV_X;   //!
		TBranch        *b_muminus_OWNPV_Y;   //!
		TBranch        *b_muminus_OWNPV_Z;   //!
		TBranch        *b_muminus_OWNPV_XERR;   //!
		TBranch        *b_muminus_OWNPV_YERR;   //!
		TBranch        *b_muminus_OWNPV_ZERR;   //!
		TBranch        *b_muminus_OWNPV_CHI2;   //!
		TBranch        *b_muminus_OWNPV_NDOF;   //!
		TBranch        *b_muminus_OWNPV_COV_;   //!
		TBranch        *b_muminus_IP_OWNPV;   //!
		TBranch        *b_muminus_IPCHI2_OWNPV;   //!
		TBranch        *b_muminus_ORIVX_X;   //!
		TBranch        *b_muminus_ORIVX_Y;   //!
		TBranch        *b_muminus_ORIVX_Z;   //!
		TBranch        *b_muminus_ORIVX_XERR;   //!
		TBranch        *b_muminus_ORIVX_YERR;   //!
		TBranch        *b_muminus_ORIVX_ZERR;   //!
		TBranch        *b_muminus_ORIVX_CHI2;   //!
		TBranch        *b_muminus_ORIVX_NDOF;   //!
		TBranch        *b_muminus_ORIVX_COV_;   //!
		TBranch        *b_muminus_P;   //!
		TBranch        *b_muminus_PT;   //!
		TBranch        *b_muminus_PE;   //!
		TBranch        *b_muminus_PX;   //!
		TBranch        *b_muminus_PY;   //!
		TBranch        *b_muminus_PZ;   //!
		TBranch        *b_muminus_M;   //!
		TBranch        *b_muminus_ID;   //!
		TBranch        *b_muminus_PIDe;   //!
		TBranch        *b_muminus_PIDmu;   //!
		TBranch        *b_muminus_PIDK;   //!
		TBranch        *b_muminus_PIDp;   //!
		TBranch        *b_muminus_PIDd;   //!
		TBranch        *b_muminus_ProbNNe;   //!
		TBranch        *b_muminus_ProbNNk;   //!
		TBranch        *b_muminus_ProbNNp;   //!
		TBranch        *b_muminus_ProbNNpi;   //!
		TBranch        *b_muminus_ProbNNmu;   //!
		TBranch        *b_muminus_ProbNNd;   //!
		TBranch        *b_muminus_ProbNNghost;   //!
		TBranch        *b_muminus_hasMuon;   //!
		TBranch        *b_muminus_isMuon;   //!
		TBranch        *b_muminus_hasRich;   //!
		TBranch        *b_muminus_UsedRichAerogel;   //!
		TBranch        *b_muminus_UsedRich1Gas;   //!
		TBranch        *b_muminus_UsedRich2Gas;   //!
		TBranch        *b_muminus_RichAboveElThres;   //!
		TBranch        *b_muminus_RichAboveMuThres;   //!
		TBranch        *b_muminus_RichAbovePiThres;   //!
		TBranch        *b_muminus_RichAboveKaThres;   //!
		TBranch        *b_muminus_RichAbovePrThres;   //!
		TBranch        *b_muminus_hasCalo;   //!
		TBranch        *b_muminus_PP_CombDLLe;   //!
		TBranch        *b_muminus_PP_CombDLLmu;   //!
		TBranch        *b_muminus_PP_CombDLLpi;   //!
		TBranch        *b_muminus_PP_CombDLLk;   //!
		TBranch        *b_muminus_PP_CombDLLp;   //!
		TBranch        *b_muminus_PP_CombDLLd;   //!
		TBranch        *b_muminus_PP_ProbNNe;   //!
		TBranch        *b_muminus_PP_ProbNNmu;   //!
		TBranch        *b_muminus_PP_ProbNNpi;   //!
		TBranch        *b_muminus_PP_ProbNNk;   //!
		TBranch        *b_muminus_PP_ProbNNp;   //!
		TBranch        *b_muminus_PP_ProbNNghost;   //!
		TBranch        *b_muminus_PP_ProbNNd;   //!
		TBranch        *b_muminus_L0Global_Dec;   //!
		TBranch        *b_muminus_L0Global_TIS;   //!
		TBranch        *b_muminus_L0Global_TOS;   //!
		TBranch        *b_muminus_Hlt1Global_Dec;   //!
		TBranch        *b_muminus_Hlt1Global_TIS;   //!
		TBranch        *b_muminus_Hlt1Global_TOS;   //!
		TBranch        *b_muminus_Hlt1Phys_Dec;   //!
		TBranch        *b_muminus_Hlt1Phys_TIS;   //!
		TBranch        *b_muminus_Hlt1Phys_TOS;   //!
		TBranch        *b_muminus_Hlt2Global_Dec;   //!
		TBranch        *b_muminus_Hlt2Global_TIS;   //!
		TBranch        *b_muminus_Hlt2Global_TOS;   //!
		TBranch        *b_muminus_Hlt2Phys_Dec;   //!
		TBranch        *b_muminus_Hlt2Phys_TIS;   //!
		TBranch        *b_muminus_Hlt2Phys_TOS;   //!
		TBranch        *b_muminus_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_muminus_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_muminus_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_muminus_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_muminus_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_muminus_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_muminus_L0MuonDecision_Dec;   //!
		TBranch        *b_muminus_L0MuonDecision_TIS;   //!
		TBranch        *b_muminus_L0MuonDecision_TOS;   //!
		TBranch        *b_muminus_L0MuonEWDecision_Dec;   //!
		TBranch        *b_muminus_L0MuonEWDecision_TIS;   //!
		TBranch        *b_muminus_L0MuonEWDecision_TOS;   //!
		TBranch        *b_muminus_L0DiMuonDecision_Dec;   //!
		TBranch        *b_muminus_L0DiMuonDecision_TIS;   //!
		TBranch        *b_muminus_L0DiMuonDecision_TOS;   //!
		TBranch        *b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_muminus_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_muminus_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_muminus_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_muminus_Hlt1TrackMVADecision_Dec;   //!
		TBranch        *b_muminus_Hlt1TrackMVADecision_TIS;   //!
		TBranch        *b_muminus_Hlt1TrackMVADecision_TOS;   //!
		TBranch        *b_muminus_Hlt1TwoTrackMVADecision_Dec;   //!
		TBranch        *b_muminus_Hlt1TwoTrackMVADecision_TIS;   //!
		TBranch        *b_muminus_Hlt1TwoTrackMVADecision_TOS;   //!
		TBranch        *b_muminus_Hlt1SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_muminus_Hlt1SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_muminus_Hlt1SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_muminus_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
		TBranch        *b_muminus_Hlt2SingleMuonHighPTDecision_Dec;   //!
		TBranch        *b_muminus_Hlt2SingleMuonHighPTDecision_TIS;   //!
		TBranch        *b_muminus_Hlt2SingleMuonHighPTDecision_TOS;   //!
		TBranch        *b_muminus_TRACK_Type;   //!
		TBranch        *b_muminus_TRACK_Key;   //!
		TBranch        *b_muminus_TRACK_CHI2NDOF;   //!
		TBranch        *b_muminus_TRACK_PCHI2;   //!
		TBranch        *b_muminus_TRACK_MatchCHI2;   //!
		TBranch        *b_muminus_TRACK_GhostProb;   //!
		TBranch        *b_muminus_TRACK_CloneDist;   //!
		TBranch        *b_muminus_TRACK_Likelihood;   //!
		TBranch        *b_nCandidate;   //!
		TBranch        *b_totCandidates;   //!
		TBranch        *b_EventInSequence;   //!
		TBranch        *b_LoKi_nCharged;   //!
		TBranch        *b_LoKi_nITClusters;   //!
		TBranch        *b_LoKi_nNeutrals;   //!
		TBranch        *b_LoKi_nOThits;   //!
		TBranch        *b_LoKi_nPVs;   //!
		TBranch        *b_LoKi_nSpdMult;   //!
		TBranch        *b_LoKi_nTTClusters;   //!
		TBranch        *b_LoKi_nVeloClusters;   //!
		TBranch        *b_LoKi_nVeloLiteClusters;   //!
		TBranch        *b_runNumber;   //!
		TBranch        *b_eventNumber;   //!
		TBranch        *b_BCID;   //!
		TBranch        *b_BCType;   //!
		TBranch        *b_OdinTCK;   //!
		TBranch        *b_L0DUTCK;   //!
		TBranch        *b_HLT1TCK;   //!
		TBranch        *b_HLT2TCK;   //!
		TBranch        *b_GpsTime;   //!
		TBranch        *b_Polarity;   //!
		TBranch        *b_B00;   //!
		TBranch        *b_B01;   //!
		TBranch        *b_B02;   //!
		TBranch        *b_B03;   //!
		TBranch        *b_B10;   //!
		TBranch        *b_B11;   //!
		TBranch        *b_B12;   //!
		TBranch        *b_B13;   //!
		TBranch        *b_B20;   //!
		TBranch        *b_B21;   //!
		TBranch        *b_B22;   //!
		TBranch        *b_B23;   //!
		TBranch        *b_F10;   //!
		TBranch        *b_F11;   //!
		TBranch        *b_F12;   //!
		TBranch        *b_F13;   //!
		TBranch        *b_F20;   //!
		TBranch        *b_F21;   //!
		TBranch        *b_F22;   //!
		TBranch        *b_F23;   //!
		TBranch        *b_log_hrc_fom_v3;   //!
		TBranch        *b_log_hrc_fom_B_v3;   //!
		TBranch        *b_log_hrc_fom_F_v3;   //!
		TBranch        *b_nchB;   //!
		TBranch        *b_adc_B;   //!
		TBranch        *b_nchF;   //!
		TBranch        *b_adc_F;   //!
		TBranch        *b_nPV;   //!
		TBranch        *b_PVX;   //!
		TBranch        *b_PVY;   //!
		TBranch        *b_PVZ;   //!
		TBranch        *b_PVXERR;   //!
		TBranch        *b_PVYERR;   //!
		TBranch        *b_PVZERR;   //!
		TBranch        *b_PVCHI2;   //!
		TBranch        *b_PVNDOF;   //!
		TBranch        *b_PVNTRACKS;   //!
		TBranch        *b_nPVs;   //!
		TBranch        *b_nTracks;   //!
		TBranch        *b_nLongTracks;   //!
		TBranch        *b_nDownstreamTracks;   //!
		TBranch        *b_nUpstreamTracks;   //!
		TBranch        *b_nVeloTracks;   //!
		TBranch        *b_nTTracks;   //!
		TBranch        *b_nBackTracks;   //!
		TBranch        *b_nRich1Hits;   //!
		TBranch        *b_nRich2Hits;   //!
		TBranch        *b_nVeloClusters;   //!
		TBranch        *b_nITClusters;   //!
		TBranch        *b_nTTClusters;   //!
		TBranch        *b_nOTClusters;   //!
		TBranch        *b_nSPDHits;   //!
		TBranch        *b_nMuonCoordsS0;   //!
		TBranch        *b_nMuonCoordsS1;   //!
		TBranch        *b_nMuonCoordsS2;   //!
		TBranch        *b_nMuonCoordsS3;   //!
		TBranch        *b_nMuonCoordsS4;   //!
		TBranch        *b_nMuonTracks;   //!

		skim(int job, int sjob=-1, TString dir="/tmp/dcraik", TString which="");
		virtual ~skim();
		virtual void     Init(TTree *tree);
		virtual void     InitZ(TTree *tree);
		virtual void     Loop();

		TString outPath;

		std::vector<TString> phiNames;
		std::vector<TChain*> phiChains;

		//TChain *zChain;
		TChain *lumi;
};

#endif

#ifdef skim_cxx
skim::skim(int job, int sjob, TString dir, TString which)
{
	outPath=dir;
	outPath+="/";
	outPath+=job;
	outPath+="/";
	outPath+=sjob;

	char str[256];

	phiNames.push_back("");
	//phiNames.push_back("WS");
	phiNames.push_back("PiPi");
	//phiNames.push_back("PiPiWS");

	phiChains.push_back(new TChain("Tuple/DecayTree"));
	//phiChains.push_back(new TChain("WSTuple/DecayTree"));
	//phiChains.push_back(new TChain("PiPiTuple/DecayTree"));
	phiChains.push_back(new TChain("Tuple/DecayTree"));//no PID cut at DV now
	//phiChains.push_back(new TChain("PiPiWSTuple/DecayTree"));

	//zChain = new TChain("ZTuple/DecayTree");
	lumi = new TChain("GetIntegratedLuminosity/LumiTuple");
	//	boost::progress_display show_addfile_progress( 700 );
	for(int i=0; i<3000; ++i) {
		//		++show_addfile_progress;
		if(sjob<0 || sjob==i) {
			sprintf(str,"%s/%d/%d/Tuples.root",dir.Data(),job,i);
			if(gSystem->AccessPathName(str)) continue;
			for(auto it=phiChains.begin(); it!=phiChains.end(); ++it) {
				(*it)->Add(str);
			}
			//zChain->Add(str);
			lumi->Add(str);
		}
	}

	for(auto it=phiChains.begin(); it!=phiChains.end(); ++it) {
		Init(*it);
		std::cout << (*it)->GetEntries() << std::endl;
	}
	//InitZ(zChain);
	//std::cout << zChain->GetEntries() << std::endl;
}

skim::~skim()
{
}

void skim::Init(TTree *tree)
{
	if (!tree) return;
	tree->SetMakeClass(1);

	tree->SetBranchAddress("phi_ENDVERTEX_X", &phi_ENDVERTEX_X, &b_phi_ENDVERTEX_X);
	tree->SetBranchAddress("phi_ENDVERTEX_Y", &phi_ENDVERTEX_Y, &b_phi_ENDVERTEX_Y);
	tree->SetBranchAddress("phi_ENDVERTEX_Z", &phi_ENDVERTEX_Z, &b_phi_ENDVERTEX_Z);
	tree->SetBranchAddress("phi_ENDVERTEX_XERR", &phi_ENDVERTEX_XERR, &b_phi_ENDVERTEX_XERR);
	tree->SetBranchAddress("phi_ENDVERTEX_YERR", &phi_ENDVERTEX_YERR, &b_phi_ENDVERTEX_YERR);
	tree->SetBranchAddress("phi_ENDVERTEX_ZERR", &phi_ENDVERTEX_ZERR, &b_phi_ENDVERTEX_ZERR);
	tree->SetBranchAddress("phi_ENDVERTEX_CHI2", &phi_ENDVERTEX_CHI2, &b_phi_ENDVERTEX_CHI2);
	tree->SetBranchAddress("phi_ENDVERTEX_NDOF", &phi_ENDVERTEX_NDOF, &b_phi_ENDVERTEX_NDOF);
	tree->SetBranchAddress("phi_ENDVERTEX_COV_", phi_ENDVERTEX_COV_, &b_phi_ENDVERTEX_COV_);
	tree->SetBranchAddress("phi_OWNPV_X", &phi_OWNPV_X, &b_phi_OWNPV_X);
	tree->SetBranchAddress("phi_OWNPV_Y", &phi_OWNPV_Y, &b_phi_OWNPV_Y);
	tree->SetBranchAddress("phi_OWNPV_Z", &phi_OWNPV_Z, &b_phi_OWNPV_Z);
	tree->SetBranchAddress("phi_OWNPV_XERR", &phi_OWNPV_XERR, &b_phi_OWNPV_XERR);
	tree->SetBranchAddress("phi_OWNPV_YERR", &phi_OWNPV_YERR, &b_phi_OWNPV_YERR);
	tree->SetBranchAddress("phi_OWNPV_ZERR", &phi_OWNPV_ZERR, &b_phi_OWNPV_ZERR);
	tree->SetBranchAddress("phi_OWNPV_CHI2", &phi_OWNPV_CHI2, &b_phi_OWNPV_CHI2);
	tree->SetBranchAddress("phi_OWNPV_NDOF", &phi_OWNPV_NDOF, &b_phi_OWNPV_NDOF);
	tree->SetBranchAddress("phi_OWNPV_COV_", phi_OWNPV_COV_, &b_phi_OWNPV_COV_);
	tree->SetBranchAddress("phi_IP_OWNPV", &phi_IP_OWNPV, &b_phi_IP_OWNPV);
	tree->SetBranchAddress("phi_IPCHI2_OWNPV", &phi_IPCHI2_OWNPV, &b_phi_IPCHI2_OWNPV);
	tree->SetBranchAddress("phi_FD_OWNPV", &phi_FD_OWNPV, &b_phi_FD_OWNPV);
	tree->SetBranchAddress("phi_FDCHI2_OWNPV", &phi_FDCHI2_OWNPV, &b_phi_FDCHI2_OWNPV);
	tree->SetBranchAddress("phi_DIRA_OWNPV", &phi_DIRA_OWNPV, &b_phi_DIRA_OWNPV);
	tree->SetBranchAddress("phi_P", &phi_P, &b_phi_P);
	tree->SetBranchAddress("phi_PT", &phi_PT, &b_phi_PT);
	tree->SetBranchAddress("phi_PE", &phi_PE, &b_phi_PE);
	tree->SetBranchAddress("phi_PX", &phi_PX, &b_phi_PX);
	tree->SetBranchAddress("phi_PY", &phi_PY, &b_phi_PY);
	tree->SetBranchAddress("phi_PZ", &phi_PZ, &b_phi_PZ);
	tree->SetBranchAddress("phi_MM", &phi_MM, &b_phi_MM);
	tree->SetBranchAddress("phi_MMERR", &phi_MMERR, &b_phi_MMERR);
	tree->SetBranchAddress("phi_M", &phi_M, &b_phi_M);
	tree->SetBranchAddress("phi_ID", &phi_ID, &b_phi_ID);
	tree->SetBranchAddress("phi_L0Global_Dec", &phi_L0Global_Dec, &b_phi_L0Global_Dec);
	tree->SetBranchAddress("phi_L0Global_TIS", &phi_L0Global_TIS, &b_phi_L0Global_TIS);
	tree->SetBranchAddress("phi_L0Global_TOS", &phi_L0Global_TOS, &b_phi_L0Global_TOS);
	tree->SetBranchAddress("phi_Hlt1Global_Dec", &phi_Hlt1Global_Dec, &b_phi_Hlt1Global_Dec);
	tree->SetBranchAddress("phi_Hlt1Global_TIS", &phi_Hlt1Global_TIS, &b_phi_Hlt1Global_TIS);
	tree->SetBranchAddress("phi_Hlt1Global_TOS", &phi_Hlt1Global_TOS, &b_phi_Hlt1Global_TOS);
	tree->SetBranchAddress("phi_Hlt1Phys_Dec", &phi_Hlt1Phys_Dec, &b_phi_Hlt1Phys_Dec);
	tree->SetBranchAddress("phi_Hlt1Phys_TIS", &phi_Hlt1Phys_TIS, &b_phi_Hlt1Phys_TIS);
	tree->SetBranchAddress("phi_Hlt1Phys_TOS", &phi_Hlt1Phys_TOS, &b_phi_Hlt1Phys_TOS);
	tree->SetBranchAddress("phi_Hlt2Global_Dec", &phi_Hlt2Global_Dec, &b_phi_Hlt2Global_Dec);
	tree->SetBranchAddress("phi_Hlt2Global_TIS", &phi_Hlt2Global_TIS, &b_phi_Hlt2Global_TIS);
	tree->SetBranchAddress("phi_Hlt2Global_TOS", &phi_Hlt2Global_TOS, &b_phi_Hlt2Global_TOS);
	tree->SetBranchAddress("phi_Hlt2Phys_Dec", &phi_Hlt2Phys_Dec, &b_phi_Hlt2Phys_Dec);
	tree->SetBranchAddress("phi_Hlt2Phys_TIS", &phi_Hlt2Phys_TIS, &b_phi_Hlt2Phys_TIS);
	tree->SetBranchAddress("phi_Hlt2Phys_TOS", &phi_Hlt2Phys_TOS, &b_phi_Hlt2Phys_TOS);
	//tree->SetBranchAddress("phi_L0DiHadron,lowMultDecision_Dec", &phi_L0DiHadron_lowMultDecision_Dec, &b_phi_L0DiHadron_lowMultDecision_Dec);
	//tree->SetBranchAddress("phi_L0DiHadron,lowMultDecision_TIS", &phi_L0DiHadron_lowMultDecision_TIS, &b_phi_L0DiHadron_lowMultDecision_TIS);
	//tree->SetBranchAddress("phi_L0DiHadron,lowMultDecision_TOS", &phi_L0DiHadron_lowMultDecision_TOS, &b_phi_L0DiHadron_lowMultDecision_TOS);
	//tree->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_Dec", &phi_L0HRCDiHadron_lowMultDecision_Dec, &b_phi_L0HRCDiHadron_lowMultDecision_Dec);
	//tree->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_TIS", &phi_L0HRCDiHadron_lowMultDecision_TIS, &b_phi_L0HRCDiHadron_lowMultDecision_TIS);
	//tree->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_TOS", &phi_L0HRCDiHadron_lowMultDecision_TOS, &b_phi_L0HRCDiHadron_lowMultDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_Dec", &phi_Hlt1LowMultHerschelDecision_Dec, &b_phi_Hlt1LowMultHerschelDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_TIS", &phi_Hlt1LowMultHerschelDecision_TIS, &b_phi_Hlt1LowMultHerschelDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_TOS", &phi_Hlt1LowMultHerschelDecision_TOS, &b_phi_Hlt1LowMultHerschelDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_Dec", &phi_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_TIS", &phi_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_TOS", &phi_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_Dec", &phi_Hlt2LowMultLMR2HHDecision_Dec, &b_phi_Hlt2LowMultLMR2HHDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_TIS", &phi_Hlt2LowMultLMR2HHDecision_TIS, &b_phi_Hlt2LowMultLMR2HHDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_TOS", &phi_Hlt2LowMultLMR2HHDecision_TOS, &b_phi_Hlt2LowMultLMR2HHDecision_TOS);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_Dec", &phi_Hlt2LowMultLMR2HHWSDecision_Dec, &b_phi_Hlt2LowMultLMR2HHWSDecision_Dec);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_TIS", &phi_Hlt2LowMultLMR2HHWSDecision_TIS, &b_phi_Hlt2LowMultLMR2HHWSDecision_TIS);
	//tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_TOS", &phi_Hlt2LowMultLMR2HHWSDecision_TOS, &b_phi_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("phi_L0DiHadron,lowMultDecision_Dec",                &phi_L0DiHadron_lowMultDecision_Dec,                &b_phi_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0Muon,lowMultDecision_Dec",                    &phi_L0Muon_lowMultDecision_Dec,                    &b_phi_L0Muon_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0DiMuon,lowMultDecision_Dec",                  &phi_L0DiMuon_lowMultDecision_Dec,                  &b_phi_L0DiMuon_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0Electron,lowMultDecision_Dec",                &phi_L0Electron_lowMultDecision_Dec,                &b_phi_L0Electron_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0Photon,lowMultDecision_Dec",                  &phi_L0Photon_lowMultDecision_Dec,                  &b_phi_L0Photon_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0DiEM,lowMultDecision_Dec",                    &phi_L0DiEM_lowMultDecision_Dec,                    &b_phi_L0DiEM_lowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0SPDDecision_Dec",                             &phi_L0SPDDecision_Dec,                             &b_phi_L0SPDDecision_Dec);
	tree->SetBranchAddress("phi_L0CALODecision_Dec",                            &phi_L0CALODecision_Dec,                            &b_phi_L0CALODecision_Dec);
	tree->SetBranchAddress("phi_L0PUDecision_Dec",                              &phi_L0PUDecision_Dec,                              &b_phi_L0PUDecision_Dec);
	tree->SetBranchAddress("phi_L0MuonDecision_Dec",                            &phi_L0MuonDecision_Dec,                            &b_phi_L0MuonDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1LowMultPassThroughDecision_Dec",            &phi_Hlt1LowMultPassThroughDecision_Dec,            &b_phi_Hlt1LowMultPassThroughDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_Dec",        &phi_Hlt1LowMultVeloCut_HadronsDecision_Dec,        &b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec",&phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec,&b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1BBMicroBiasVeloDecision_Dec",               &phi_Hlt1BBMicroBiasVeloDecision_Dec,               &b_phi_Hlt1BBMicroBiasVeloDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1BBHighMultDecision_Dec",                    &phi_Hlt1BBHighMultDecision_Dec,                    &b_phi_Hlt1BBHighMultDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1MBNoBiasDecision_Dec",                      &phi_Hlt1MBNoBiasDecision_Dec,                      &b_phi_Hlt1MBNoBiasDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_Dec",                 &phi_Hlt2LowMultLMR2HHDecision_Dec,                 &b_phi_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_Dec",               &phi_Hlt2LowMultLMR2HHWSDecision_Dec,               &b_phi_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2MBMicroBiasVeloDecision_Dec",               &phi_Hlt2MBMicroBiasVeloDecision_Dec,               &b_phi_Hlt2MBMicroBiasVeloDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2MBHighMultDecision_Dec",                    &phi_Hlt2MBHighMultDecision_Dec,                    &b_phi_Hlt2MBHighMultDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2MBNoBiasDecision_Dec",                      &phi_Hlt2MBNoBiasDecision_Dec,                      &b_phi_Hlt2MBNoBiasDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2PassThroughDecision_Dec",                   &phi_Hlt2PassThroughDecision_Dec,                   &b_phi_Hlt2PassThroughDecision_Dec);
	tree->SetBranchAddress("phi_L0HadronDecision_Dec",                          &phi_L0HadronDecision_Dec,                          &b_phi_L0HadronDecision_Dec);
//	tree->SetBranchAddress("phi_L0MuonDecision_Dec",                            &phi_L0MuonDecision_Dec,                            &b_phi_L0MuonDecision_Dec);
//	tree->SetBranchAddress("phi_L0SPDDecision_Dec",                             &phi_L0SPDDecision_Dec,                             &b_phi_L0SPDDecision_Dec);
	tree->SetBranchAddress("phi_L0HadronLowMultDecision_Dec",                   &phi_L0HadronLowMultDecision_Dec,                   &b_phi_L0HadronLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0MuonLowMultDecision_Dec",                     &phi_L0MuonLowMultDecision_Dec,                     &b_phi_L0MuonLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0ElectronLowMultDecision_Dec",                 &phi_L0ElectronLowMultDecision_Dec,                 &b_phi_L0ElectronLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0PhotonLowMultDecision_Dec",                   &phi_L0PhotonLowMultDecision_Dec,                   &b_phi_L0PhotonLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0DiEMLowMultDecision_Dec",                     &phi_L0DiEMLowMultDecision_Dec,                     &b_phi_L0DiEMLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0SPDLowMultDecision_Dec",                      &phi_L0SPDLowMultDecision_Dec,                      &b_phi_L0SPDLowMultDecision_Dec);
	tree->SetBranchAddress("phi_L0SoftCEPDecision_Dec",                         &phi_L0SoftCEPDecision_Dec,                         &b_phi_L0SoftCEPDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec",        &phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec,        &b_phi_Hlt1BBMicroBiasLowMultVeloDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1BBMicroBiasSoftCEPDecision_Dec",            &phi_Hlt1BBMicroBiasSoftCEPDecision_Dec,            &b_phi_Hlt1BBMicroBiasSoftCEPDecision_Dec);
	tree->SetBranchAddress("phi_Hlt1BBHasTrackDecision_Dec",                    &phi_Hlt1BBHasTrackDecision_Dec,                    &b_phi_Hlt1BBHasTrackDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2BBPassThroughDecision_Dec",                 &phi_Hlt2BBPassThroughDecision_Dec,                 &b_phi_Hlt2BBPassThroughDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec",        &phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec,        &b_phi_Hlt2MBMicroBiasLowMultVeloDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2BBLongTrackDecision_Dec",                   &phi_Hlt2BBLongTrackDecision_Dec,                   &b_phi_Hlt2BBLongTrackDecision_Dec);
	tree->SetBranchAddress("phi_Hlt2SingleTrackDecision_Dec",                   &phi_Hlt2SingleTrackDecision_Dec,                   &b_phi_Hlt2SingleTrackDecision_Dec);
//	tree->SetBranchAddress("phi_Hlt2MBMicroBiasVeloDecision_Dec",               &phi_Hlt2MBMicroBiasVeloDecision_Dec,               &b_phi_Hlt2MBMicroBiasVeloDecision_Dec);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNe", &Kplus_MC12TuneV2_ProbNNe, &b_Kplus_MC12TuneV2_ProbNNe);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNmu", &Kplus_MC12TuneV2_ProbNNmu, &b_Kplus_MC12TuneV2_ProbNNmu);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNpi", &Kplus_MC12TuneV2_ProbNNpi, &b_Kplus_MC12TuneV2_ProbNNpi);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNk", &Kplus_MC12TuneV2_ProbNNk, &b_Kplus_MC12TuneV2_ProbNNk);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNp", &Kplus_MC12TuneV2_ProbNNp, &b_Kplus_MC12TuneV2_ProbNNp);
	tree->SetBranchAddress("Kplus_MC12TuneV2_ProbNNghost", &Kplus_MC12TuneV2_ProbNNghost, &b_Kplus_MC12TuneV2_ProbNNghost);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNe", &Kplus_MC12TuneV3_ProbNNe, &b_Kplus_MC12TuneV3_ProbNNe);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNmu", &Kplus_MC12TuneV3_ProbNNmu, &b_Kplus_MC12TuneV3_ProbNNmu);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNpi", &Kplus_MC12TuneV3_ProbNNpi, &b_Kplus_MC12TuneV3_ProbNNpi);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNk", &Kplus_MC12TuneV3_ProbNNk, &b_Kplus_MC12TuneV3_ProbNNk);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNp", &Kplus_MC12TuneV3_ProbNNp, &b_Kplus_MC12TuneV3_ProbNNp);
	tree->SetBranchAddress("Kplus_MC12TuneV3_ProbNNghost", &Kplus_MC12TuneV3_ProbNNghost, &b_Kplus_MC12TuneV3_ProbNNghost);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNe", &Kplus_MC12TuneV4_ProbNNe, &b_Kplus_MC12TuneV4_ProbNNe);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNmu", &Kplus_MC12TuneV4_ProbNNmu, &b_Kplus_MC12TuneV4_ProbNNmu);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNpi", &Kplus_MC12TuneV4_ProbNNpi, &b_Kplus_MC12TuneV4_ProbNNpi);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNk", &Kplus_MC12TuneV4_ProbNNk, &b_Kplus_MC12TuneV4_ProbNNk);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNp", &Kplus_MC12TuneV4_ProbNNp, &b_Kplus_MC12TuneV4_ProbNNp);
	tree->SetBranchAddress("Kplus_MC12TuneV4_ProbNNghost", &Kplus_MC12TuneV4_ProbNNghost, &b_Kplus_MC12TuneV4_ProbNNghost);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNe", &Kplus_MC15TuneV1_ProbNNe, &b_Kplus_MC15TuneV1_ProbNNe);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNmu", &Kplus_MC15TuneV1_ProbNNmu, &b_Kplus_MC15TuneV1_ProbNNmu);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNpi", &Kplus_MC15TuneV1_ProbNNpi, &b_Kplus_MC15TuneV1_ProbNNpi);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNk", &Kplus_MC15TuneV1_ProbNNk, &b_Kplus_MC15TuneV1_ProbNNk);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNp", &Kplus_MC15TuneV1_ProbNNp, &b_Kplus_MC15TuneV1_ProbNNp);
	tree->SetBranchAddress("Kplus_MC15TuneV1_ProbNNghost", &Kplus_MC15TuneV1_ProbNNghost, &b_Kplus_MC15TuneV1_ProbNNghost);
	tree->SetBranchAddress("Kplus_OWNPV_X", &Kplus_OWNPV_X, &b_Kplus_OWNPV_X);
	tree->SetBranchAddress("Kplus_OWNPV_Y", &Kplus_OWNPV_Y, &b_Kplus_OWNPV_Y);
	tree->SetBranchAddress("Kplus_OWNPV_Z", &Kplus_OWNPV_Z, &b_Kplus_OWNPV_Z);
	tree->SetBranchAddress("Kplus_OWNPV_XERR", &Kplus_OWNPV_XERR, &b_Kplus_OWNPV_XERR);
	tree->SetBranchAddress("Kplus_OWNPV_YERR", &Kplus_OWNPV_YERR, &b_Kplus_OWNPV_YERR);
	tree->SetBranchAddress("Kplus_OWNPV_ZERR", &Kplus_OWNPV_ZERR, &b_Kplus_OWNPV_ZERR);
	tree->SetBranchAddress("Kplus_OWNPV_CHI2", &Kplus_OWNPV_CHI2, &b_Kplus_OWNPV_CHI2);
	tree->SetBranchAddress("Kplus_OWNPV_NDOF", &Kplus_OWNPV_NDOF, &b_Kplus_OWNPV_NDOF);
	tree->SetBranchAddress("Kplus_OWNPV_COV_", Kplus_OWNPV_COV_, &b_Kplus_OWNPV_COV_);
	tree->SetBranchAddress("Kplus_IP_OWNPV", &Kplus_IP_OWNPV, &b_Kplus_IP_OWNPV);
	tree->SetBranchAddress("Kplus_IPCHI2_OWNPV", &Kplus_IPCHI2_OWNPV, &b_Kplus_IPCHI2_OWNPV);
	tree->SetBranchAddress("Kplus_ORIVX_X", &Kplus_ORIVX_X, &b_Kplus_ORIVX_X);
	tree->SetBranchAddress("Kplus_ORIVX_Y", &Kplus_ORIVX_Y, &b_Kplus_ORIVX_Y);
	tree->SetBranchAddress("Kplus_ORIVX_Z", &Kplus_ORIVX_Z, &b_Kplus_ORIVX_Z);
	tree->SetBranchAddress("Kplus_ORIVX_XERR", &Kplus_ORIVX_XERR, &b_Kplus_ORIVX_XERR);
	tree->SetBranchAddress("Kplus_ORIVX_YERR", &Kplus_ORIVX_YERR, &b_Kplus_ORIVX_YERR);
	tree->SetBranchAddress("Kplus_ORIVX_ZERR", &Kplus_ORIVX_ZERR, &b_Kplus_ORIVX_ZERR);
	tree->SetBranchAddress("Kplus_ORIVX_CHI2", &Kplus_ORIVX_CHI2, &b_Kplus_ORIVX_CHI2);
	tree->SetBranchAddress("Kplus_ORIVX_NDOF", &Kplus_ORIVX_NDOF, &b_Kplus_ORIVX_NDOF);
	tree->SetBranchAddress("Kplus_ORIVX_COV_", Kplus_ORIVX_COV_, &b_Kplus_ORIVX_COV_);
	tree->SetBranchAddress("Kplus_P", &Kplus_P, &b_Kplus_P);
	tree->SetBranchAddress("Kplus_PT", &Kplus_PT, &b_Kplus_PT);
	tree->SetBranchAddress("Kplus_PE", &Kplus_PE, &b_Kplus_PE);
	tree->SetBranchAddress("Kplus_PX", &Kplus_PX, &b_Kplus_PX);
	tree->SetBranchAddress("Kplus_PY", &Kplus_PY, &b_Kplus_PY);
	tree->SetBranchAddress("Kplus_PZ", &Kplus_PZ, &b_Kplus_PZ);
	tree->SetBranchAddress("Kplus_M", &Kplus_M, &b_Kplus_M);
	tree->SetBranchAddress("Kplus_ID", &Kplus_ID, &b_Kplus_ID);
	tree->SetBranchAddress("Kplus_PIDe", &Kplus_PIDe, &b_Kplus_PIDe);
	tree->SetBranchAddress("Kplus_PIDmu", &Kplus_PIDmu, &b_Kplus_PIDmu);
	tree->SetBranchAddress("Kplus_PIDK", &Kplus_PIDK, &b_Kplus_PIDK);
	tree->SetBranchAddress("Kplus_PIDp", &Kplus_PIDp, &b_Kplus_PIDp);
	tree->SetBranchAddress("Kplus_ProbNNe", &Kplus_ProbNNe, &b_Kplus_ProbNNe);
	tree->SetBranchAddress("Kplus_ProbNNk", &Kplus_ProbNNk, &b_Kplus_ProbNNk);
	tree->SetBranchAddress("Kplus_ProbNNp", &Kplus_ProbNNp, &b_Kplus_ProbNNp);
	tree->SetBranchAddress("Kplus_ProbNNpi", &Kplus_ProbNNpi, &b_Kplus_ProbNNpi);
	tree->SetBranchAddress("Kplus_ProbNNmu", &Kplus_ProbNNmu, &b_Kplus_ProbNNmu);
	tree->SetBranchAddress("Kplus_ProbNNghost", &Kplus_ProbNNghost, &b_Kplus_ProbNNghost);
	tree->SetBranchAddress("Kplus_hasMuon", &Kplus_hasMuon, &b_Kplus_hasMuon);
	tree->SetBranchAddress("Kplus_isMuon", &Kplus_isMuon, &b_Kplus_isMuon);
	tree->SetBranchAddress("Kplus_hasRich", &Kplus_hasRich, &b_Kplus_hasRich);
	tree->SetBranchAddress("Kplus_UsedRichAerogel", &Kplus_UsedRichAerogel, &b_Kplus_UsedRichAerogel);
	tree->SetBranchAddress("Kplus_UsedRich1Gas", &Kplus_UsedRich1Gas, &b_Kplus_UsedRich1Gas);
	tree->SetBranchAddress("Kplus_UsedRich2Gas", &Kplus_UsedRich2Gas, &b_Kplus_UsedRich2Gas);
	tree->SetBranchAddress("Kplus_RichAboveElThres", &Kplus_RichAboveElThres, &b_Kplus_RichAboveElThres);
	tree->SetBranchAddress("Kplus_RichAboveMuThres", &Kplus_RichAboveMuThres, &b_Kplus_RichAboveMuThres);
	tree->SetBranchAddress("Kplus_RichAbovePiThres", &Kplus_RichAbovePiThres, &b_Kplus_RichAbovePiThres);
	tree->SetBranchAddress("Kplus_RichAboveKaThres", &Kplus_RichAboveKaThres, &b_Kplus_RichAboveKaThres);
	tree->SetBranchAddress("Kplus_RichAbovePrThres", &Kplus_RichAbovePrThres, &b_Kplus_RichAbovePrThres);
	tree->SetBranchAddress("Kplus_hasCalo", &Kplus_hasCalo, &b_Kplus_hasCalo);
	tree->SetBranchAddress("Kplus_PP_CombDLLe", &Kplus_PP_CombDLLe, &b_Kplus_PP_CombDLLe);
	tree->SetBranchAddress("Kplus_PP_CombDLLmu", &Kplus_PP_CombDLLmu, &b_Kplus_PP_CombDLLmu);
	tree->SetBranchAddress("Kplus_PP_CombDLLpi", &Kplus_PP_CombDLLpi, &b_Kplus_PP_CombDLLpi);
	tree->SetBranchAddress("Kplus_PP_CombDLLk", &Kplus_PP_CombDLLk, &b_Kplus_PP_CombDLLk);
	tree->SetBranchAddress("Kplus_PP_CombDLLp", &Kplus_PP_CombDLLp, &b_Kplus_PP_CombDLLp);
	tree->SetBranchAddress("Kplus_PP_CombDLLd", &Kplus_PP_CombDLLd, &b_Kplus_PP_CombDLLd);
	tree->SetBranchAddress("Kplus_PP_ProbNNe", &Kplus_PP_ProbNNe, &b_Kplus_PP_ProbNNe);
	tree->SetBranchAddress("Kplus_PP_ProbNNmu", &Kplus_PP_ProbNNmu, &b_Kplus_PP_ProbNNmu);
	tree->SetBranchAddress("Kplus_PP_ProbNNpi", &Kplus_PP_ProbNNpi, &b_Kplus_PP_ProbNNpi);
	tree->SetBranchAddress("Kplus_PP_ProbNNk", &Kplus_PP_ProbNNk, &b_Kplus_PP_ProbNNk);
	tree->SetBranchAddress("Kplus_PP_ProbNNp", &Kplus_PP_ProbNNp, &b_Kplus_PP_ProbNNp);
	tree->SetBranchAddress("Kplus_PP_ProbNNghost", &Kplus_PP_ProbNNghost, &b_Kplus_PP_ProbNNghost);
	tree->SetBranchAddress("Kplus_PP_ProbNNd", &Kplus_PP_ProbNNd, &b_Kplus_PP_ProbNNd);
	tree->SetBranchAddress("Kplus_L0Global_Dec", &Kplus_L0Global_Dec, &b_Kplus_L0Global_Dec);
	tree->SetBranchAddress("Kplus_L0Global_TIS", &Kplus_L0Global_TIS, &b_Kplus_L0Global_TIS);
	tree->SetBranchAddress("Kplus_L0Global_TOS", &Kplus_L0Global_TOS, &b_Kplus_L0Global_TOS);
	tree->SetBranchAddress("Kplus_Hlt1Global_Dec", &Kplus_Hlt1Global_Dec, &b_Kplus_Hlt1Global_Dec);
	tree->SetBranchAddress("Kplus_Hlt1Global_TIS", &Kplus_Hlt1Global_TIS, &b_Kplus_Hlt1Global_TIS);
	tree->SetBranchAddress("Kplus_Hlt1Global_TOS", &Kplus_Hlt1Global_TOS, &b_Kplus_Hlt1Global_TOS);
	tree->SetBranchAddress("Kplus_Hlt1Phys_Dec", &Kplus_Hlt1Phys_Dec, &b_Kplus_Hlt1Phys_Dec);
	tree->SetBranchAddress("Kplus_Hlt1Phys_TIS", &Kplus_Hlt1Phys_TIS, &b_Kplus_Hlt1Phys_TIS);
	tree->SetBranchAddress("Kplus_Hlt1Phys_TOS", &Kplus_Hlt1Phys_TOS, &b_Kplus_Hlt1Phys_TOS);
	tree->SetBranchAddress("Kplus_Hlt2Global_Dec", &Kplus_Hlt2Global_Dec, &b_Kplus_Hlt2Global_Dec);
	tree->SetBranchAddress("Kplus_Hlt2Global_TIS", &Kplus_Hlt2Global_TIS, &b_Kplus_Hlt2Global_TIS);
	tree->SetBranchAddress("Kplus_Hlt2Global_TOS", &Kplus_Hlt2Global_TOS, &b_Kplus_Hlt2Global_TOS);
	tree->SetBranchAddress("Kplus_Hlt2Phys_Dec", &Kplus_Hlt2Phys_Dec, &b_Kplus_Hlt2Phys_Dec);
	tree->SetBranchAddress("Kplus_Hlt2Phys_TIS", &Kplus_Hlt2Phys_TIS, &b_Kplus_Hlt2Phys_TIS);
	tree->SetBranchAddress("Kplus_Hlt2Phys_TOS", &Kplus_Hlt2Phys_TOS, &b_Kplus_Hlt2Phys_TOS);
	tree->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_Dec", &Kplus_L0DiHadron_lowMultDecision_Dec, &b_Kplus_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_TIS", &Kplus_L0DiHadron_lowMultDecision_TIS, &b_Kplus_L0DiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_TOS", &Kplus_L0DiHadron_lowMultDecision_TOS, &b_Kplus_L0DiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_Dec", &Kplus_L0HRCDiHadron_lowMultDecision_Dec, &b_Kplus_L0HRCDiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_TIS", &Kplus_L0HRCDiHadron_lowMultDecision_TIS, &b_Kplus_L0HRCDiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_TOS", &Kplus_L0HRCDiHadron_lowMultDecision_TOS, &b_Kplus_L0HRCDiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_Dec", &Kplus_Hlt1LowMultHerschelDecision_Dec, &b_Kplus_Hlt1LowMultHerschelDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_TIS", &Kplus_Hlt1LowMultHerschelDecision_TIS, &b_Kplus_Hlt1LowMultHerschelDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_TOS", &Kplus_Hlt1LowMultHerschelDecision_TOS, &b_Kplus_Hlt1LowMultHerschelDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_Dec", &Kplus_Hlt2LowMultLMR2HHDecision_Dec, &b_Kplus_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_TIS", &Kplus_Hlt2LowMultLMR2HHDecision_TIS, &b_Kplus_Hlt2LowMultLMR2HHDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_TOS", &Kplus_Hlt2LowMultLMR2HHDecision_TOS, &b_Kplus_Hlt2LowMultLMR2HHDecision_TOS);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_Dec", &Kplus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_TIS", &Kplus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_TIS);
	tree->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_TOS", &Kplus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("Kplus_TRACK_Type", &Kplus_TRACK_Type, &b_Kplus_TRACK_Type);
	tree->SetBranchAddress("Kplus_TRACK_Key", &Kplus_TRACK_Key, &b_Kplus_TRACK_Key);
	tree->SetBranchAddress("Kplus_TRACK_CHI2NDOF", &Kplus_TRACK_CHI2NDOF, &b_Kplus_TRACK_CHI2NDOF);
	tree->SetBranchAddress("Kplus_TRACK_PCHI2", &Kplus_TRACK_PCHI2, &b_Kplus_TRACK_PCHI2);
	tree->SetBranchAddress("Kplus_TRACK_MatchCHI2", &Kplus_TRACK_MatchCHI2, &b_Kplus_TRACK_MatchCHI2);
	tree->SetBranchAddress("Kplus_TRACK_GhostProb", &Kplus_TRACK_GhostProb, &b_Kplus_TRACK_GhostProb);
	tree->SetBranchAddress("Kplus_TRACK_CloneDist", &Kplus_TRACK_CloneDist, &b_Kplus_TRACK_CloneDist);
	tree->SetBranchAddress("Kplus_TRACK_Likelihood", &Kplus_TRACK_Likelihood, &b_Kplus_TRACK_Likelihood);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNe", &Kminus_MC12TuneV2_ProbNNe, &b_Kminus_MC12TuneV2_ProbNNe);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNmu", &Kminus_MC12TuneV2_ProbNNmu, &b_Kminus_MC12TuneV2_ProbNNmu);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNpi", &Kminus_MC12TuneV2_ProbNNpi, &b_Kminus_MC12TuneV2_ProbNNpi);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNk", &Kminus_MC12TuneV2_ProbNNk, &b_Kminus_MC12TuneV2_ProbNNk);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNp", &Kminus_MC12TuneV2_ProbNNp, &b_Kminus_MC12TuneV2_ProbNNp);
	tree->SetBranchAddress("Kminus_MC12TuneV2_ProbNNghost", &Kminus_MC12TuneV2_ProbNNghost, &b_Kminus_MC12TuneV2_ProbNNghost);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNe", &Kminus_MC12TuneV3_ProbNNe, &b_Kminus_MC12TuneV3_ProbNNe);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNmu", &Kminus_MC12TuneV3_ProbNNmu, &b_Kminus_MC12TuneV3_ProbNNmu);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNpi", &Kminus_MC12TuneV3_ProbNNpi, &b_Kminus_MC12TuneV3_ProbNNpi);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNk", &Kminus_MC12TuneV3_ProbNNk, &b_Kminus_MC12TuneV3_ProbNNk);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNp", &Kminus_MC12TuneV3_ProbNNp, &b_Kminus_MC12TuneV3_ProbNNp);
	tree->SetBranchAddress("Kminus_MC12TuneV3_ProbNNghost", &Kminus_MC12TuneV3_ProbNNghost, &b_Kminus_MC12TuneV3_ProbNNghost);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNe", &Kminus_MC12TuneV4_ProbNNe, &b_Kminus_MC12TuneV4_ProbNNe);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNmu", &Kminus_MC12TuneV4_ProbNNmu, &b_Kminus_MC12TuneV4_ProbNNmu);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNpi", &Kminus_MC12TuneV4_ProbNNpi, &b_Kminus_MC12TuneV4_ProbNNpi);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNk", &Kminus_MC12TuneV4_ProbNNk, &b_Kminus_MC12TuneV4_ProbNNk);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNp", &Kminus_MC12TuneV4_ProbNNp, &b_Kminus_MC12TuneV4_ProbNNp);
	tree->SetBranchAddress("Kminus_MC12TuneV4_ProbNNghost", &Kminus_MC12TuneV4_ProbNNghost, &b_Kminus_MC12TuneV4_ProbNNghost);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNe", &Kminus_MC15TuneV1_ProbNNe, &b_Kminus_MC15TuneV1_ProbNNe);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNmu", &Kminus_MC15TuneV1_ProbNNmu, &b_Kminus_MC15TuneV1_ProbNNmu);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNpi", &Kminus_MC15TuneV1_ProbNNpi, &b_Kminus_MC15TuneV1_ProbNNpi);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNk", &Kminus_MC15TuneV1_ProbNNk, &b_Kminus_MC15TuneV1_ProbNNk);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNp", &Kminus_MC15TuneV1_ProbNNp, &b_Kminus_MC15TuneV1_ProbNNp);
	tree->SetBranchAddress("Kminus_MC15TuneV1_ProbNNghost", &Kminus_MC15TuneV1_ProbNNghost, &b_Kminus_MC15TuneV1_ProbNNghost);
	tree->SetBranchAddress("Kminus_OWNPV_X", &Kminus_OWNPV_X, &b_Kminus_OWNPV_X);
	tree->SetBranchAddress("Kminus_OWNPV_Y", &Kminus_OWNPV_Y, &b_Kminus_OWNPV_Y);
	tree->SetBranchAddress("Kminus_OWNPV_Z", &Kminus_OWNPV_Z, &b_Kminus_OWNPV_Z);
	tree->SetBranchAddress("Kminus_OWNPV_XERR", &Kminus_OWNPV_XERR, &b_Kminus_OWNPV_XERR);
	tree->SetBranchAddress("Kminus_OWNPV_YERR", &Kminus_OWNPV_YERR, &b_Kminus_OWNPV_YERR);
	tree->SetBranchAddress("Kminus_OWNPV_ZERR", &Kminus_OWNPV_ZERR, &b_Kminus_OWNPV_ZERR);
	tree->SetBranchAddress("Kminus_OWNPV_CHI2", &Kminus_OWNPV_CHI2, &b_Kminus_OWNPV_CHI2);
	tree->SetBranchAddress("Kminus_OWNPV_NDOF", &Kminus_OWNPV_NDOF, &b_Kminus_OWNPV_NDOF);
	tree->SetBranchAddress("Kminus_OWNPV_COV_", Kminus_OWNPV_COV_, &b_Kminus_OWNPV_COV_);
	tree->SetBranchAddress("Kminus_IP_OWNPV", &Kminus_IP_OWNPV, &b_Kminus_IP_OWNPV);
	tree->SetBranchAddress("Kminus_IPCHI2_OWNPV", &Kminus_IPCHI2_OWNPV, &b_Kminus_IPCHI2_OWNPV);
	tree->SetBranchAddress("Kminus_ORIVX_X", &Kminus_ORIVX_X, &b_Kminus_ORIVX_X);
	tree->SetBranchAddress("Kminus_ORIVX_Y", &Kminus_ORIVX_Y, &b_Kminus_ORIVX_Y);
	tree->SetBranchAddress("Kminus_ORIVX_Z", &Kminus_ORIVX_Z, &b_Kminus_ORIVX_Z);
	tree->SetBranchAddress("Kminus_ORIVX_XERR", &Kminus_ORIVX_XERR, &b_Kminus_ORIVX_XERR);
	tree->SetBranchAddress("Kminus_ORIVX_YERR", &Kminus_ORIVX_YERR, &b_Kminus_ORIVX_YERR);
	tree->SetBranchAddress("Kminus_ORIVX_ZERR", &Kminus_ORIVX_ZERR, &b_Kminus_ORIVX_ZERR);
	tree->SetBranchAddress("Kminus_ORIVX_CHI2", &Kminus_ORIVX_CHI2, &b_Kminus_ORIVX_CHI2);
	tree->SetBranchAddress("Kminus_ORIVX_NDOF", &Kminus_ORIVX_NDOF, &b_Kminus_ORIVX_NDOF);
	tree->SetBranchAddress("Kminus_ORIVX_COV_", Kminus_ORIVX_COV_, &b_Kminus_ORIVX_COV_);
	tree->SetBranchAddress("Kminus_P", &Kminus_P, &b_Kminus_P);
	tree->SetBranchAddress("Kminus_PT", &Kminus_PT, &b_Kminus_PT);
	tree->SetBranchAddress("Kminus_PE", &Kminus_PE, &b_Kminus_PE);
	tree->SetBranchAddress("Kminus_PX", &Kminus_PX, &b_Kminus_PX);
	tree->SetBranchAddress("Kminus_PY", &Kminus_PY, &b_Kminus_PY);
	tree->SetBranchAddress("Kminus_PZ", &Kminus_PZ, &b_Kminus_PZ);
	tree->SetBranchAddress("Kminus_M", &Kminus_M, &b_Kminus_M);
	tree->SetBranchAddress("Kminus_ID", &Kminus_ID, &b_Kminus_ID);
	tree->SetBranchAddress("Kminus_PIDe", &Kminus_PIDe, &b_Kminus_PIDe);
	tree->SetBranchAddress("Kminus_PIDmu", &Kminus_PIDmu, &b_Kminus_PIDmu);
	tree->SetBranchAddress("Kminus_PIDK", &Kminus_PIDK, &b_Kminus_PIDK);
	tree->SetBranchAddress("Kminus_PIDp", &Kminus_PIDp, &b_Kminus_PIDp);
	tree->SetBranchAddress("Kminus_ProbNNe", &Kminus_ProbNNe, &b_Kminus_ProbNNe);
	tree->SetBranchAddress("Kminus_ProbNNk", &Kminus_ProbNNk, &b_Kminus_ProbNNk);
	tree->SetBranchAddress("Kminus_ProbNNp", &Kminus_ProbNNp, &b_Kminus_ProbNNp);
	tree->SetBranchAddress("Kminus_ProbNNpi", &Kminus_ProbNNpi, &b_Kminus_ProbNNpi);
	tree->SetBranchAddress("Kminus_ProbNNmu", &Kminus_ProbNNmu, &b_Kminus_ProbNNmu);
	tree->SetBranchAddress("Kminus_ProbNNghost", &Kminus_ProbNNghost, &b_Kminus_ProbNNghost);
	tree->SetBranchAddress("Kminus_hasMuon", &Kminus_hasMuon, &b_Kminus_hasMuon);
	tree->SetBranchAddress("Kminus_isMuon", &Kminus_isMuon, &b_Kminus_isMuon);
	tree->SetBranchAddress("Kminus_hasRich", &Kminus_hasRich, &b_Kminus_hasRich);
	tree->SetBranchAddress("Kminus_UsedRichAerogel", &Kminus_UsedRichAerogel, &b_Kminus_UsedRichAerogel);
	tree->SetBranchAddress("Kminus_UsedRich1Gas", &Kminus_UsedRich1Gas, &b_Kminus_UsedRich1Gas);
	tree->SetBranchAddress("Kminus_UsedRich2Gas", &Kminus_UsedRich2Gas, &b_Kminus_UsedRich2Gas);
	tree->SetBranchAddress("Kminus_RichAboveElThres", &Kminus_RichAboveElThres, &b_Kminus_RichAboveElThres);
	tree->SetBranchAddress("Kminus_RichAboveMuThres", &Kminus_RichAboveMuThres, &b_Kminus_RichAboveMuThres);
	tree->SetBranchAddress("Kminus_RichAbovePiThres", &Kminus_RichAbovePiThres, &b_Kminus_RichAbovePiThres);
	tree->SetBranchAddress("Kminus_RichAboveKaThres", &Kminus_RichAboveKaThres, &b_Kminus_RichAboveKaThres);
	tree->SetBranchAddress("Kminus_RichAbovePrThres", &Kminus_RichAbovePrThres, &b_Kminus_RichAbovePrThres);
	tree->SetBranchAddress("Kminus_hasCalo", &Kminus_hasCalo, &b_Kminus_hasCalo);
	tree->SetBranchAddress("Kminus_PP_CombDLLe", &Kminus_PP_CombDLLe, &b_Kminus_PP_CombDLLe);
	tree->SetBranchAddress("Kminus_PP_CombDLLmu", &Kminus_PP_CombDLLmu, &b_Kminus_PP_CombDLLmu);
	tree->SetBranchAddress("Kminus_PP_CombDLLpi", &Kminus_PP_CombDLLpi, &b_Kminus_PP_CombDLLpi);
	tree->SetBranchAddress("Kminus_PP_CombDLLk", &Kminus_PP_CombDLLk, &b_Kminus_PP_CombDLLk);
	tree->SetBranchAddress("Kminus_PP_CombDLLp", &Kminus_PP_CombDLLp, &b_Kminus_PP_CombDLLp);
	tree->SetBranchAddress("Kminus_PP_CombDLLd", &Kminus_PP_CombDLLd, &b_Kminus_PP_CombDLLd);
	tree->SetBranchAddress("Kminus_PP_ProbNNe", &Kminus_PP_ProbNNe, &b_Kminus_PP_ProbNNe);
	tree->SetBranchAddress("Kminus_PP_ProbNNmu", &Kminus_PP_ProbNNmu, &b_Kminus_PP_ProbNNmu);
	tree->SetBranchAddress("Kminus_PP_ProbNNpi", &Kminus_PP_ProbNNpi, &b_Kminus_PP_ProbNNpi);
	tree->SetBranchAddress("Kminus_PP_ProbNNk", &Kminus_PP_ProbNNk, &b_Kminus_PP_ProbNNk);
	tree->SetBranchAddress("Kminus_PP_ProbNNp", &Kminus_PP_ProbNNp, &b_Kminus_PP_ProbNNp);
	tree->SetBranchAddress("Kminus_PP_ProbNNghost", &Kminus_PP_ProbNNghost, &b_Kminus_PP_ProbNNghost);
	tree->SetBranchAddress("Kminus_PP_ProbNNd", &Kminus_PP_ProbNNd, &b_Kminus_PP_ProbNNd);
	tree->SetBranchAddress("Kminus_L0Global_Dec", &Kminus_L0Global_Dec, &b_Kminus_L0Global_Dec);
	tree->SetBranchAddress("Kminus_L0Global_TIS", &Kminus_L0Global_TIS, &b_Kminus_L0Global_TIS);
	tree->SetBranchAddress("Kminus_L0Global_TOS", &Kminus_L0Global_TOS, &b_Kminus_L0Global_TOS);
	tree->SetBranchAddress("Kminus_Hlt1Global_Dec", &Kminus_Hlt1Global_Dec, &b_Kminus_Hlt1Global_Dec);
	tree->SetBranchAddress("Kminus_Hlt1Global_TIS", &Kminus_Hlt1Global_TIS, &b_Kminus_Hlt1Global_TIS);
	tree->SetBranchAddress("Kminus_Hlt1Global_TOS", &Kminus_Hlt1Global_TOS, &b_Kminus_Hlt1Global_TOS);
	tree->SetBranchAddress("Kminus_Hlt1Phys_Dec", &Kminus_Hlt1Phys_Dec, &b_Kminus_Hlt1Phys_Dec);
	tree->SetBranchAddress("Kminus_Hlt1Phys_TIS", &Kminus_Hlt1Phys_TIS, &b_Kminus_Hlt1Phys_TIS);
	tree->SetBranchAddress("Kminus_Hlt1Phys_TOS", &Kminus_Hlt1Phys_TOS, &b_Kminus_Hlt1Phys_TOS);
	tree->SetBranchAddress("Kminus_Hlt2Global_Dec", &Kminus_Hlt2Global_Dec, &b_Kminus_Hlt2Global_Dec);
	tree->SetBranchAddress("Kminus_Hlt2Global_TIS", &Kminus_Hlt2Global_TIS, &b_Kminus_Hlt2Global_TIS);
	tree->SetBranchAddress("Kminus_Hlt2Global_TOS", &Kminus_Hlt2Global_TOS, &b_Kminus_Hlt2Global_TOS);
	tree->SetBranchAddress("Kminus_Hlt2Phys_Dec", &Kminus_Hlt2Phys_Dec, &b_Kminus_Hlt2Phys_Dec);
	tree->SetBranchAddress("Kminus_Hlt2Phys_TIS", &Kminus_Hlt2Phys_TIS, &b_Kminus_Hlt2Phys_TIS);
	tree->SetBranchAddress("Kminus_Hlt2Phys_TOS", &Kminus_Hlt2Phys_TOS, &b_Kminus_Hlt2Phys_TOS);
	tree->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_Dec", &Kminus_L0DiHadron_lowMultDecision_Dec, &b_Kminus_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_TIS", &Kminus_L0DiHadron_lowMultDecision_TIS, &b_Kminus_L0DiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_TOS", &Kminus_L0DiHadron_lowMultDecision_TOS, &b_Kminus_L0DiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_Dec", &Kminus_L0HRCDiHadron_lowMultDecision_Dec, &b_Kminus_L0HRCDiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_TIS", &Kminus_L0HRCDiHadron_lowMultDecision_TIS, &b_Kminus_L0HRCDiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_TOS", &Kminus_L0HRCDiHadron_lowMultDecision_TOS, &b_Kminus_L0HRCDiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_Dec", &Kminus_Hlt1LowMultHerschelDecision_Dec, &b_Kminus_Hlt1LowMultHerschelDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_TIS", &Kminus_Hlt1LowMultHerschelDecision_TIS, &b_Kminus_Hlt1LowMultHerschelDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_TOS", &Kminus_Hlt1LowMultHerschelDecision_TOS, &b_Kminus_Hlt1LowMultHerschelDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_Dec", &Kminus_Hlt2LowMultLMR2HHDecision_Dec, &b_Kminus_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_TIS", &Kminus_Hlt2LowMultLMR2HHDecision_TIS, &b_Kminus_Hlt2LowMultLMR2HHDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_TOS", &Kminus_Hlt2LowMultLMR2HHDecision_TOS, &b_Kminus_Hlt2LowMultLMR2HHDecision_TOS);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_Dec", &Kminus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_TIS", &Kminus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_TIS);
	tree->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_TOS", &Kminus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("Kminus_TRACK_Type", &Kminus_TRACK_Type, &b_Kminus_TRACK_Type);
	tree->SetBranchAddress("Kminus_TRACK_Key", &Kminus_TRACK_Key, &b_Kminus_TRACK_Key);
	tree->SetBranchAddress("Kminus_TRACK_CHI2NDOF", &Kminus_TRACK_CHI2NDOF, &b_Kminus_TRACK_CHI2NDOF);
	tree->SetBranchAddress("Kminus_TRACK_PCHI2", &Kminus_TRACK_PCHI2, &b_Kminus_TRACK_PCHI2);
	tree->SetBranchAddress("Kminus_TRACK_MatchCHI2", &Kminus_TRACK_MatchCHI2, &b_Kminus_TRACK_MatchCHI2);
	tree->SetBranchAddress("Kminus_TRACK_GhostProb", &Kminus_TRACK_GhostProb, &b_Kminus_TRACK_GhostProb);
	tree->SetBranchAddress("Kminus_TRACK_CloneDist", &Kminus_TRACK_CloneDist, &b_Kminus_TRACK_CloneDist);
	tree->SetBranchAddress("Kminus_TRACK_Likelihood", &Kminus_TRACK_Likelihood, &b_Kminus_TRACK_Likelihood);
	tree->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
	tree->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
	tree->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
	tree->SetBranchAddress("LoKi_nCharged", &LoKi_nCharged, &b_LoKi_nCharged);
	tree->SetBranchAddress("LoKi_nITClusters", &LoKi_nITClusters, &b_LoKi_nITClusters);
	tree->SetBranchAddress("LoKi_nNeutrals", &LoKi_nNeutrals, &b_LoKi_nNeutrals);
	tree->SetBranchAddress("LoKi_nOThits", &LoKi_nOThits, &b_LoKi_nOThits);
	tree->SetBranchAddress("LoKi_nPVs", &LoKi_nPVs, &b_LoKi_nPVs);
	tree->SetBranchAddress("LoKi_nSpdMult", &LoKi_nSpdMult, &b_LoKi_nSpdMult);
	tree->SetBranchAddress("LoKi_nTTClusters", &LoKi_nTTClusters, &b_LoKi_nTTClusters);
	tree->SetBranchAddress("LoKi_nVeloClusters", &LoKi_nVeloClusters, &b_LoKi_nVeloClusters);
	tree->SetBranchAddress("LoKi_nVeloLiteClusters", &LoKi_nVeloLiteClusters, &b_LoKi_nVeloLiteClusters);
	tree->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	tree->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
	tree->SetBranchAddress("BCID", &BCID, &b_BCID);
	tree->SetBranchAddress("BCType", &BCType, &b_BCType);
	tree->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
	tree->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
	tree->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
	tree->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
	tree->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
	tree->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
	tree->SetBranchAddress("B00", &B00, &b_B00);
	tree->SetBranchAddress("B01", &B01, &b_B01);
	tree->SetBranchAddress("B02", &B02, &b_B02);
	tree->SetBranchAddress("B03", &B03, &b_B03);
	tree->SetBranchAddress("B10", &B10, &b_B10);
	tree->SetBranchAddress("B11", &B11, &b_B11);
	tree->SetBranchAddress("B12", &B12, &b_B12);
	tree->SetBranchAddress("B13", &B13, &b_B13);
	tree->SetBranchAddress("B20", &B20, &b_B20);
	tree->SetBranchAddress("B21", &B21, &b_B21);
	tree->SetBranchAddress("B22", &B22, &b_B22);
	tree->SetBranchAddress("B23", &B23, &b_B23);
	tree->SetBranchAddress("F10", &F10, &b_F10);
	tree->SetBranchAddress("F11", &F11, &b_F11);
	tree->SetBranchAddress("F12", &F12, &b_F12);
	tree->SetBranchAddress("F13", &F13, &b_F13);
	tree->SetBranchAddress("F20", &F20, &b_F20);
	tree->SetBranchAddress("F21", &F21, &b_F21);
	tree->SetBranchAddress("F22", &F22, &b_F22);
	tree->SetBranchAddress("F23", &F23, &b_F23);
	tree->SetBranchAddress("log_hrc_fom_v3", &log_hrc_fom_v3, &b_log_hrc_fom_v3);
	tree->SetBranchAddress("log_hrc_fom_B_v3", &log_hrc_fom_B_v3, &b_log_hrc_fom_B_v3);
	tree->SetBranchAddress("log_hrc_fom_F_v3", &log_hrc_fom_F_v3, &b_log_hrc_fom_F_v3);
	tree->SetBranchAddress("nchB", &nchB, &b_nchB);
	tree->SetBranchAddress("adc_B", adc_B, &b_adc_B);
	tree->SetBranchAddress("nchF", &nchF, &b_nchF);
	tree->SetBranchAddress("adc_F", adc_F, &b_adc_F);
	tree->SetBranchAddress("nPV", &nPV, &b_nPV);
	tree->SetBranchAddress("PVX", PVX, &b_PVX);
	tree->SetBranchAddress("PVY", PVY, &b_PVY);
	tree->SetBranchAddress("PVZ", PVZ, &b_PVZ);
	tree->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
	tree->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
	tree->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
	tree->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
	tree->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
	tree->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
	tree->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
	tree->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	tree->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
	tree->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
	tree->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
	tree->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
	tree->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
	tree->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
	tree->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
	tree->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
	tree->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
	tree->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
	tree->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
	tree->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
	tree->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
	tree->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
	tree->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
	tree->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
	tree->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
	tree->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
	tree->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
}

void skim::InitZ(TTree *tree)
{
	if (!tree) return;
	tree->SetMakeClass(1);

	tree->SetBranchAddress("Z_ENDVERTEX_X", &Z_ENDVERTEX_X, &b_Z_ENDVERTEX_X);
	tree->SetBranchAddress("Z_ENDVERTEX_Y", &Z_ENDVERTEX_Y, &b_Z_ENDVERTEX_Y);
	tree->SetBranchAddress("Z_ENDVERTEX_Z", &Z_ENDVERTEX_Z, &b_Z_ENDVERTEX_Z);
	tree->SetBranchAddress("Z_ENDVERTEX_XERR", &Z_ENDVERTEX_XERR, &b_Z_ENDVERTEX_XERR);
	tree->SetBranchAddress("Z_ENDVERTEX_YERR", &Z_ENDVERTEX_YERR, &b_Z_ENDVERTEX_YERR);
	tree->SetBranchAddress("Z_ENDVERTEX_ZERR", &Z_ENDVERTEX_ZERR, &b_Z_ENDVERTEX_ZERR);
	tree->SetBranchAddress("Z_ENDVERTEX_CHI2", &Z_ENDVERTEX_CHI2, &b_Z_ENDVERTEX_CHI2);
	tree->SetBranchAddress("Z_ENDVERTEX_NDOF", &Z_ENDVERTEX_NDOF, &b_Z_ENDVERTEX_NDOF);
	tree->SetBranchAddress("Z_ENDVERTEX_COV_", Z_ENDVERTEX_COV_, &b_Z_ENDVERTEX_COV_);
	tree->SetBranchAddress("Z_OWNPV_X", &Z_OWNPV_X, &b_Z_OWNPV_X);
	tree->SetBranchAddress("Z_OWNPV_Y", &Z_OWNPV_Y, &b_Z_OWNPV_Y);
	tree->SetBranchAddress("Z_OWNPV_Z", &Z_OWNPV_Z, &b_Z_OWNPV_Z);
	tree->SetBranchAddress("Z_OWNPV_XERR", &Z_OWNPV_XERR, &b_Z_OWNPV_XERR);
	tree->SetBranchAddress("Z_OWNPV_YERR", &Z_OWNPV_YERR, &b_Z_OWNPV_YERR);
	tree->SetBranchAddress("Z_OWNPV_ZERR", &Z_OWNPV_ZERR, &b_Z_OWNPV_ZERR);
	tree->SetBranchAddress("Z_OWNPV_CHI2", &Z_OWNPV_CHI2, &b_Z_OWNPV_CHI2);
	tree->SetBranchAddress("Z_OWNPV_NDOF", &Z_OWNPV_NDOF, &b_Z_OWNPV_NDOF);
	tree->SetBranchAddress("Z_OWNPV_COV_", Z_OWNPV_COV_, &b_Z_OWNPV_COV_);
	tree->SetBranchAddress("Z_IP_OWNPV", &Z_IP_OWNPV, &b_Z_IP_OWNPV);
	tree->SetBranchAddress("Z_IPCHI2_OWNPV", &Z_IPCHI2_OWNPV, &b_Z_IPCHI2_OWNPV);
	tree->SetBranchAddress("Z_FD_OWNPV", &Z_FD_OWNPV, &b_Z_FD_OWNPV);
	tree->SetBranchAddress("Z_FDCHI2_OWNPV", &Z_FDCHI2_OWNPV, &b_Z_FDCHI2_OWNPV);
	tree->SetBranchAddress("Z_DIRA_OWNPV", &Z_DIRA_OWNPV, &b_Z_DIRA_OWNPV);
	tree->SetBranchAddress("Z_P", &Z_P, &b_Z_P);
	tree->SetBranchAddress("Z_PT", &Z_PT, &b_Z_PT);
	tree->SetBranchAddress("Z_PE", &Z_PE, &b_Z_PE);
	tree->SetBranchAddress("Z_PX", &Z_PX, &b_Z_PX);
	tree->SetBranchAddress("Z_PY", &Z_PY, &b_Z_PY);
	tree->SetBranchAddress("Z_PZ", &Z_PZ, &b_Z_PZ);
	tree->SetBranchAddress("Z_MM", &Z_MM, &b_Z_MM);
	tree->SetBranchAddress("Z_MMERR", &Z_MMERR, &b_Z_MMERR);
	tree->SetBranchAddress("Z_M", &Z_M, &b_Z_M);
	tree->SetBranchAddress("Z_ID", &Z_ID, &b_Z_ID);
	tree->SetBranchAddress("Z_L0Global_Dec", &Z_L0Global_Dec, &b_Z_L0Global_Dec);
	tree->SetBranchAddress("Z_L0Global_TIS", &Z_L0Global_TIS, &b_Z_L0Global_TIS);
	tree->SetBranchAddress("Z_L0Global_TOS", &Z_L0Global_TOS, &b_Z_L0Global_TOS);
	tree->SetBranchAddress("Z_Hlt1Global_Dec", &Z_Hlt1Global_Dec, &b_Z_Hlt1Global_Dec);
	tree->SetBranchAddress("Z_Hlt1Global_TIS", &Z_Hlt1Global_TIS, &b_Z_Hlt1Global_TIS);
	tree->SetBranchAddress("Z_Hlt1Global_TOS", &Z_Hlt1Global_TOS, &b_Z_Hlt1Global_TOS);
	tree->SetBranchAddress("Z_Hlt1Phys_Dec", &Z_Hlt1Phys_Dec, &b_Z_Hlt1Phys_Dec);
	tree->SetBranchAddress("Z_Hlt1Phys_TIS", &Z_Hlt1Phys_TIS, &b_Z_Hlt1Phys_TIS);
	tree->SetBranchAddress("Z_Hlt1Phys_TOS", &Z_Hlt1Phys_TOS, &b_Z_Hlt1Phys_TOS);
	tree->SetBranchAddress("Z_Hlt2Global_Dec", &Z_Hlt2Global_Dec, &b_Z_Hlt2Global_Dec);
	tree->SetBranchAddress("Z_Hlt2Global_TIS", &Z_Hlt2Global_TIS, &b_Z_Hlt2Global_TIS);
	tree->SetBranchAddress("Z_Hlt2Global_TOS", &Z_Hlt2Global_TOS, &b_Z_Hlt2Global_TOS);
	tree->SetBranchAddress("Z_Hlt2Phys_Dec", &Z_Hlt2Phys_Dec, &b_Z_Hlt2Phys_Dec);
	tree->SetBranchAddress("Z_Hlt2Phys_TIS", &Z_Hlt2Phys_TIS, &b_Z_Hlt2Phys_TIS);
	tree->SetBranchAddress("Z_Hlt2Phys_TOS", &Z_Hlt2Phys_TOS, &b_Z_Hlt2Phys_TOS);
	tree->SetBranchAddress("Z_L0DiHadron,lowMultDecision_Dec", &Z_L0DiHadron_lowMultDecision_Dec, &b_Z_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Z_L0DiHadron,lowMultDecision_TIS", &Z_L0DiHadron_lowMultDecision_TIS, &b_Z_L0DiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Z_L0DiHadron,lowMultDecision_TOS", &Z_L0DiHadron_lowMultDecision_TOS, &b_Z_L0DiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Z_L0HRCDiHadron,lowMultDecision_Dec", &Z_L0HRCDiHadron_lowMultDecision_Dec, &b_Z_L0HRCDiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("Z_L0HRCDiHadron,lowMultDecision_TIS", &Z_L0HRCDiHadron_lowMultDecision_TIS, &b_Z_L0HRCDiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("Z_L0HRCDiHadron,lowMultDecision_TOS", &Z_L0HRCDiHadron_lowMultDecision_TOS, &b_Z_L0HRCDiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("Z_L0MuonDecision_Dec", &Z_L0MuonDecision_Dec, &b_Z_L0MuonDecision_Dec);
	tree->SetBranchAddress("Z_L0MuonDecision_TIS", &Z_L0MuonDecision_TIS, &b_Z_L0MuonDecision_TIS);
	tree->SetBranchAddress("Z_L0MuonDecision_TOS", &Z_L0MuonDecision_TOS, &b_Z_L0MuonDecision_TOS);
	tree->SetBranchAddress("Z_L0MuonEWDecision_Dec", &Z_L0MuonEWDecision_Dec, &b_Z_L0MuonEWDecision_Dec);
	tree->SetBranchAddress("Z_L0MuonEWDecision_TIS", &Z_L0MuonEWDecision_TIS, &b_Z_L0MuonEWDecision_TIS);
	tree->SetBranchAddress("Z_L0MuonEWDecision_TOS", &Z_L0MuonEWDecision_TOS, &b_Z_L0MuonEWDecision_TOS);
	tree->SetBranchAddress("Z_L0DiMuonDecision_Dec", &Z_L0DiMuonDecision_Dec, &b_Z_L0DiMuonDecision_Dec);
	tree->SetBranchAddress("Z_L0DiMuonDecision_TIS", &Z_L0DiMuonDecision_TIS, &b_Z_L0DiMuonDecision_TIS);
	tree->SetBranchAddress("Z_L0DiMuonDecision_TOS", &Z_L0DiMuonDecision_TOS, &b_Z_L0DiMuonDecision_TOS);
	tree->SetBranchAddress("Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	tree->SetBranchAddress("Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	tree->SetBranchAddress("Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_Z_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	tree->SetBranchAddress("Z_Hlt1LowMultHerschelDecision_Dec", &Z_Hlt1LowMultHerschelDecision_Dec, &b_Z_Hlt1LowMultHerschelDecision_Dec);
	tree->SetBranchAddress("Z_Hlt1LowMultHerschelDecision_TIS", &Z_Hlt1LowMultHerschelDecision_TIS, &b_Z_Hlt1LowMultHerschelDecision_TIS);
	tree->SetBranchAddress("Z_Hlt1LowMultHerschelDecision_TOS", &Z_Hlt1LowMultHerschelDecision_TOS, &b_Z_Hlt1LowMultHerschelDecision_TOS);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloCut_HadronsDecision_Dec", &Z_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_Z_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloCut_HadronsDecision_TIS", &Z_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_Z_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloCut_HadronsDecision_TOS", &Z_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_Z_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	tree->SetBranchAddress("Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_Z_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	tree->SetBranchAddress("Z_Hlt1TrackMVADecision_Dec", &Z_Hlt1TrackMVADecision_Dec, &b_Z_Hlt1TrackMVADecision_Dec);
	tree->SetBranchAddress("Z_Hlt1TrackMVADecision_TIS", &Z_Hlt1TrackMVADecision_TIS, &b_Z_Hlt1TrackMVADecision_TIS);
	tree->SetBranchAddress("Z_Hlt1TrackMVADecision_TOS", &Z_Hlt1TrackMVADecision_TOS, &b_Z_Hlt1TrackMVADecision_TOS);
	tree->SetBranchAddress("Z_Hlt1TwoTrackMVADecision_Dec", &Z_Hlt1TwoTrackMVADecision_Dec, &b_Z_Hlt1TwoTrackMVADecision_Dec);
	tree->SetBranchAddress("Z_Hlt1TwoTrackMVADecision_TIS", &Z_Hlt1TwoTrackMVADecision_TIS, &b_Z_Hlt1TwoTrackMVADecision_TIS);
	tree->SetBranchAddress("Z_Hlt1TwoTrackMVADecision_TOS", &Z_Hlt1TwoTrackMVADecision_TOS, &b_Z_Hlt1TwoTrackMVADecision_TOS);
	tree->SetBranchAddress("Z_Hlt1SingleMuonHighPTDecision_Dec", &Z_Hlt1SingleMuonHighPTDecision_Dec, &b_Z_Hlt1SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("Z_Hlt1SingleMuonHighPTDecision_TIS", &Z_Hlt1SingleMuonHighPTDecision_TIS, &b_Z_Hlt1SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("Z_Hlt1SingleMuonHighPTDecision_TOS", &Z_Hlt1SingleMuonHighPTDecision_TOS, &b_Z_Hlt1SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHDecision_Dec", &Z_Hlt2LowMultLMR2HHDecision_Dec, &b_Z_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHDecision_TIS", &Z_Hlt2LowMultLMR2HHDecision_TIS, &b_Z_Hlt2LowMultLMR2HHDecision_TIS);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHDecision_TOS", &Z_Hlt2LowMultLMR2HHDecision_TOS, &b_Z_Hlt2LowMultLMR2HHDecision_TOS);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHWSDecision_Dec", &Z_Hlt2LowMultLMR2HHWSDecision_Dec, &b_Z_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHWSDecision_TIS", &Z_Hlt2LowMultLMR2HHWSDecision_TIS, &b_Z_Hlt2LowMultLMR2HHWSDecision_TIS);
	tree->SetBranchAddress("Z_Hlt2LowMultLMR2HHWSDecision_TOS", &Z_Hlt2LowMultLMR2HHWSDecision_TOS, &b_Z_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("Z_Hlt2SingleMuonHighPTDecision_Dec", &Z_Hlt2SingleMuonHighPTDecision_Dec, &b_Z_Hlt2SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("Z_Hlt2SingleMuonHighPTDecision_TIS", &Z_Hlt2SingleMuonHighPTDecision_TIS, &b_Z_Hlt2SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("Z_Hlt2SingleMuonHighPTDecision_TOS", &Z_Hlt2SingleMuonHighPTDecision_TOS, &b_Z_Hlt2SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNe", &muplus_MC12TuneV2_ProbNNe, &b_muplus_MC12TuneV2_ProbNNe);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNmu", &muplus_MC12TuneV2_ProbNNmu, &b_muplus_MC12TuneV2_ProbNNmu);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNpi", &muplus_MC12TuneV2_ProbNNpi, &b_muplus_MC12TuneV2_ProbNNpi);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNk", &muplus_MC12TuneV2_ProbNNk, &b_muplus_MC12TuneV2_ProbNNk);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNp", &muplus_MC12TuneV2_ProbNNp, &b_muplus_MC12TuneV2_ProbNNp);
	tree->SetBranchAddress("muplus_MC12TuneV2_ProbNNghost", &muplus_MC12TuneV2_ProbNNghost, &b_muplus_MC12TuneV2_ProbNNghost);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNe", &muplus_MC12TuneV3_ProbNNe, &b_muplus_MC12TuneV3_ProbNNe);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNmu", &muplus_MC12TuneV3_ProbNNmu, &b_muplus_MC12TuneV3_ProbNNmu);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNpi", &muplus_MC12TuneV3_ProbNNpi, &b_muplus_MC12TuneV3_ProbNNpi);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNk", &muplus_MC12TuneV3_ProbNNk, &b_muplus_MC12TuneV3_ProbNNk);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNp", &muplus_MC12TuneV3_ProbNNp, &b_muplus_MC12TuneV3_ProbNNp);
	tree->SetBranchAddress("muplus_MC12TuneV3_ProbNNghost", &muplus_MC12TuneV3_ProbNNghost, &b_muplus_MC12TuneV3_ProbNNghost);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNe", &muplus_MC12TuneV4_ProbNNe, &b_muplus_MC12TuneV4_ProbNNe);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNmu", &muplus_MC12TuneV4_ProbNNmu, &b_muplus_MC12TuneV4_ProbNNmu);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNpi", &muplus_MC12TuneV4_ProbNNpi, &b_muplus_MC12TuneV4_ProbNNpi);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNk", &muplus_MC12TuneV4_ProbNNk, &b_muplus_MC12TuneV4_ProbNNk);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNp", &muplus_MC12TuneV4_ProbNNp, &b_muplus_MC12TuneV4_ProbNNp);
	tree->SetBranchAddress("muplus_MC12TuneV4_ProbNNghost", &muplus_MC12TuneV4_ProbNNghost, &b_muplus_MC12TuneV4_ProbNNghost);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNe", &muplus_MC15TuneV1_ProbNNe, &b_muplus_MC15TuneV1_ProbNNe);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNmu", &muplus_MC15TuneV1_ProbNNmu, &b_muplus_MC15TuneV1_ProbNNmu);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNpi", &muplus_MC15TuneV1_ProbNNpi, &b_muplus_MC15TuneV1_ProbNNpi);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNk", &muplus_MC15TuneV1_ProbNNk, &b_muplus_MC15TuneV1_ProbNNk);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNp", &muplus_MC15TuneV1_ProbNNp, &b_muplus_MC15TuneV1_ProbNNp);
	tree->SetBranchAddress("muplus_MC15TuneV1_ProbNNghost", &muplus_MC15TuneV1_ProbNNghost, &b_muplus_MC15TuneV1_ProbNNghost);
	tree->SetBranchAddress("muplus_OWNPV_X", &muplus_OWNPV_X, &b_muplus_OWNPV_X);
	tree->SetBranchAddress("muplus_OWNPV_Y", &muplus_OWNPV_Y, &b_muplus_OWNPV_Y);
	tree->SetBranchAddress("muplus_OWNPV_Z", &muplus_OWNPV_Z, &b_muplus_OWNPV_Z);
	tree->SetBranchAddress("muplus_OWNPV_XERR", &muplus_OWNPV_XERR, &b_muplus_OWNPV_XERR);
	tree->SetBranchAddress("muplus_OWNPV_YERR", &muplus_OWNPV_YERR, &b_muplus_OWNPV_YERR);
	tree->SetBranchAddress("muplus_OWNPV_ZERR", &muplus_OWNPV_ZERR, &b_muplus_OWNPV_ZERR);
	tree->SetBranchAddress("muplus_OWNPV_CHI2", &muplus_OWNPV_CHI2, &b_muplus_OWNPV_CHI2);
	tree->SetBranchAddress("muplus_OWNPV_NDOF", &muplus_OWNPV_NDOF, &b_muplus_OWNPV_NDOF);
	tree->SetBranchAddress("muplus_OWNPV_COV_", muplus_OWNPV_COV_, &b_muplus_OWNPV_COV_);
	tree->SetBranchAddress("muplus_IP_OWNPV", &muplus_IP_OWNPV, &b_muplus_IP_OWNPV);
	tree->SetBranchAddress("muplus_IPCHI2_OWNPV", &muplus_IPCHI2_OWNPV, &b_muplus_IPCHI2_OWNPV);
	tree->SetBranchAddress("muplus_ORIVX_X", &muplus_ORIVX_X, &b_muplus_ORIVX_X);
	tree->SetBranchAddress("muplus_ORIVX_Y", &muplus_ORIVX_Y, &b_muplus_ORIVX_Y);
	tree->SetBranchAddress("muplus_ORIVX_Z", &muplus_ORIVX_Z, &b_muplus_ORIVX_Z);
	tree->SetBranchAddress("muplus_ORIVX_XERR", &muplus_ORIVX_XERR, &b_muplus_ORIVX_XERR);
	tree->SetBranchAddress("muplus_ORIVX_YERR", &muplus_ORIVX_YERR, &b_muplus_ORIVX_YERR);
	tree->SetBranchAddress("muplus_ORIVX_ZERR", &muplus_ORIVX_ZERR, &b_muplus_ORIVX_ZERR);
	tree->SetBranchAddress("muplus_ORIVX_CHI2", &muplus_ORIVX_CHI2, &b_muplus_ORIVX_CHI2);
	tree->SetBranchAddress("muplus_ORIVX_NDOF", &muplus_ORIVX_NDOF, &b_muplus_ORIVX_NDOF);
	tree->SetBranchAddress("muplus_ORIVX_COV_", muplus_ORIVX_COV_, &b_muplus_ORIVX_COV_);
	tree->SetBranchAddress("muplus_P", &muplus_P, &b_muplus_P);
	tree->SetBranchAddress("muplus_PT", &muplus_PT, &b_muplus_PT);
	tree->SetBranchAddress("muplus_PE", &muplus_PE, &b_muplus_PE);
	tree->SetBranchAddress("muplus_PX", &muplus_PX, &b_muplus_PX);
	tree->SetBranchAddress("muplus_PY", &muplus_PY, &b_muplus_PY);
	tree->SetBranchAddress("muplus_PZ", &muplus_PZ, &b_muplus_PZ);
	tree->SetBranchAddress("muplus_M", &muplus_M, &b_muplus_M);
	tree->SetBranchAddress("muplus_ID", &muplus_ID, &b_muplus_ID);
	tree->SetBranchAddress("muplus_PIDe", &muplus_PIDe, &b_muplus_PIDe);
	tree->SetBranchAddress("muplus_PIDmu", &muplus_PIDmu, &b_muplus_PIDmu);
	tree->SetBranchAddress("muplus_PIDK", &muplus_PIDK, &b_muplus_PIDK);
	tree->SetBranchAddress("muplus_PIDp", &muplus_PIDp, &b_muplus_PIDp);
	tree->SetBranchAddress("muplus_PIDd", &muplus_PIDd, &b_muplus_PIDd);
	tree->SetBranchAddress("muplus_ProbNNe", &muplus_ProbNNe, &b_muplus_ProbNNe);
	tree->SetBranchAddress("muplus_ProbNNk", &muplus_ProbNNk, &b_muplus_ProbNNk);
	tree->SetBranchAddress("muplus_ProbNNp", &muplus_ProbNNp, &b_muplus_ProbNNp);
	tree->SetBranchAddress("muplus_ProbNNpi", &muplus_ProbNNpi, &b_muplus_ProbNNpi);
	tree->SetBranchAddress("muplus_ProbNNmu", &muplus_ProbNNmu, &b_muplus_ProbNNmu);
	tree->SetBranchAddress("muplus_ProbNNd", &muplus_ProbNNd, &b_muplus_ProbNNd);
	tree->SetBranchAddress("muplus_ProbNNghost", &muplus_ProbNNghost, &b_muplus_ProbNNghost);
	tree->SetBranchAddress("muplus_hasMuon", &muplus_hasMuon, &b_muplus_hasMuon);
	tree->SetBranchAddress("muplus_isMuon", &muplus_isMuon, &b_muplus_isMuon);
	tree->SetBranchAddress("muplus_hasRich", &muplus_hasRich, &b_muplus_hasRich);
	tree->SetBranchAddress("muplus_UsedRichAerogel", &muplus_UsedRichAerogel, &b_muplus_UsedRichAerogel);
	tree->SetBranchAddress("muplus_UsedRich1Gas", &muplus_UsedRich1Gas, &b_muplus_UsedRich1Gas);
	tree->SetBranchAddress("muplus_UsedRich2Gas", &muplus_UsedRich2Gas, &b_muplus_UsedRich2Gas);
	tree->SetBranchAddress("muplus_RichAboveElThres", &muplus_RichAboveElThres, &b_muplus_RichAboveElThres);
	tree->SetBranchAddress("muplus_RichAboveMuThres", &muplus_RichAboveMuThres, &b_muplus_RichAboveMuThres);
	tree->SetBranchAddress("muplus_RichAbovePiThres", &muplus_RichAbovePiThres, &b_muplus_RichAbovePiThres);
	tree->SetBranchAddress("muplus_RichAboveKaThres", &muplus_RichAboveKaThres, &b_muplus_RichAboveKaThres);
	tree->SetBranchAddress("muplus_RichAbovePrThres", &muplus_RichAbovePrThres, &b_muplus_RichAbovePrThres);
	tree->SetBranchAddress("muplus_hasCalo", &muplus_hasCalo, &b_muplus_hasCalo);
	tree->SetBranchAddress("muplus_PP_CombDLLe", &muplus_PP_CombDLLe, &b_muplus_PP_CombDLLe);
	tree->SetBranchAddress("muplus_PP_CombDLLmu", &muplus_PP_CombDLLmu, &b_muplus_PP_CombDLLmu);
	tree->SetBranchAddress("muplus_PP_CombDLLpi", &muplus_PP_CombDLLpi, &b_muplus_PP_CombDLLpi);
	tree->SetBranchAddress("muplus_PP_CombDLLk", &muplus_PP_CombDLLk, &b_muplus_PP_CombDLLk);
	tree->SetBranchAddress("muplus_PP_CombDLLp", &muplus_PP_CombDLLp, &b_muplus_PP_CombDLLp);
	tree->SetBranchAddress("muplus_PP_CombDLLd", &muplus_PP_CombDLLd, &b_muplus_PP_CombDLLd);
	tree->SetBranchAddress("muplus_PP_ProbNNe", &muplus_PP_ProbNNe, &b_muplus_PP_ProbNNe);
	tree->SetBranchAddress("muplus_PP_ProbNNmu", &muplus_PP_ProbNNmu, &b_muplus_PP_ProbNNmu);
	tree->SetBranchAddress("muplus_PP_ProbNNpi", &muplus_PP_ProbNNpi, &b_muplus_PP_ProbNNpi);
	tree->SetBranchAddress("muplus_PP_ProbNNk", &muplus_PP_ProbNNk, &b_muplus_PP_ProbNNk);
	tree->SetBranchAddress("muplus_PP_ProbNNp", &muplus_PP_ProbNNp, &b_muplus_PP_ProbNNp);
	tree->SetBranchAddress("muplus_PP_ProbNNghost", &muplus_PP_ProbNNghost, &b_muplus_PP_ProbNNghost);
	tree->SetBranchAddress("muplus_PP_ProbNNd", &muplus_PP_ProbNNd, &b_muplus_PP_ProbNNd);
	tree->SetBranchAddress("muplus_L0Global_Dec", &muplus_L0Global_Dec, &b_muplus_L0Global_Dec);
	tree->SetBranchAddress("muplus_L0Global_TIS", &muplus_L0Global_TIS, &b_muplus_L0Global_TIS);
	tree->SetBranchAddress("muplus_L0Global_TOS", &muplus_L0Global_TOS, &b_muplus_L0Global_TOS);
	tree->SetBranchAddress("muplus_Hlt1Global_Dec", &muplus_Hlt1Global_Dec, &b_muplus_Hlt1Global_Dec);
	tree->SetBranchAddress("muplus_Hlt1Global_TIS", &muplus_Hlt1Global_TIS, &b_muplus_Hlt1Global_TIS);
	tree->SetBranchAddress("muplus_Hlt1Global_TOS", &muplus_Hlt1Global_TOS, &b_muplus_Hlt1Global_TOS);
	tree->SetBranchAddress("muplus_Hlt1Phys_Dec", &muplus_Hlt1Phys_Dec, &b_muplus_Hlt1Phys_Dec);
	tree->SetBranchAddress("muplus_Hlt1Phys_TIS", &muplus_Hlt1Phys_TIS, &b_muplus_Hlt1Phys_TIS);
	tree->SetBranchAddress("muplus_Hlt1Phys_TOS", &muplus_Hlt1Phys_TOS, &b_muplus_Hlt1Phys_TOS);
	tree->SetBranchAddress("muplus_Hlt2Global_Dec", &muplus_Hlt2Global_Dec, &b_muplus_Hlt2Global_Dec);
	tree->SetBranchAddress("muplus_Hlt2Global_TIS", &muplus_Hlt2Global_TIS, &b_muplus_Hlt2Global_TIS);
	tree->SetBranchAddress("muplus_Hlt2Global_TOS", &muplus_Hlt2Global_TOS, &b_muplus_Hlt2Global_TOS);
	tree->SetBranchAddress("muplus_Hlt2Phys_Dec", &muplus_Hlt2Phys_Dec, &b_muplus_Hlt2Phys_Dec);
	tree->SetBranchAddress("muplus_Hlt2Phys_TIS", &muplus_Hlt2Phys_TIS, &b_muplus_Hlt2Phys_TIS);
	tree->SetBranchAddress("muplus_Hlt2Phys_TOS", &muplus_Hlt2Phys_TOS, &b_muplus_Hlt2Phys_TOS);
	tree->SetBranchAddress("muplus_L0DiHadron,lowMultDecision_Dec", &muplus_L0DiHadron_lowMultDecision_Dec, &b_muplus_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("muplus_L0DiHadron,lowMultDecision_TIS", &muplus_L0DiHadron_lowMultDecision_TIS, &b_muplus_L0DiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("muplus_L0DiHadron,lowMultDecision_TOS", &muplus_L0DiHadron_lowMultDecision_TOS, &b_muplus_L0DiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("muplus_L0HRCDiHadron,lowMultDecision_Dec", &muplus_L0HRCDiHadron_lowMultDecision_Dec, &b_muplus_L0HRCDiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("muplus_L0HRCDiHadron,lowMultDecision_TIS", &muplus_L0HRCDiHadron_lowMultDecision_TIS, &b_muplus_L0HRCDiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("muplus_L0HRCDiHadron,lowMultDecision_TOS", &muplus_L0HRCDiHadron_lowMultDecision_TOS, &b_muplus_L0HRCDiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("muplus_L0MuonDecision_Dec", &muplus_L0MuonDecision_Dec, &b_muplus_L0MuonDecision_Dec);
	tree->SetBranchAddress("muplus_L0MuonDecision_TIS", &muplus_L0MuonDecision_TIS, &b_muplus_L0MuonDecision_TIS);
	tree->SetBranchAddress("muplus_L0MuonDecision_TOS", &muplus_L0MuonDecision_TOS, &b_muplus_L0MuonDecision_TOS);
	tree->SetBranchAddress("muplus_L0MuonEWDecision_Dec", &muplus_L0MuonEWDecision_Dec, &b_muplus_L0MuonEWDecision_Dec);
	tree->SetBranchAddress("muplus_L0MuonEWDecision_TIS", &muplus_L0MuonEWDecision_TIS, &b_muplus_L0MuonEWDecision_TIS);
	tree->SetBranchAddress("muplus_L0MuonEWDecision_TOS", &muplus_L0MuonEWDecision_TOS, &b_muplus_L0MuonEWDecision_TOS);
	tree->SetBranchAddress("muplus_L0DiMuonDecision_Dec", &muplus_L0DiMuonDecision_Dec, &b_muplus_L0DiMuonDecision_Dec);
	tree->SetBranchAddress("muplus_L0DiMuonDecision_TIS", &muplus_L0DiMuonDecision_TIS, &b_muplus_L0DiMuonDecision_TIS);
	tree->SetBranchAddress("muplus_L0DiMuonDecision_TOS", &muplus_L0DiMuonDecision_TOS, &b_muplus_L0DiMuonDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_muplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1LowMultHerschelDecision_Dec", &muplus_Hlt1LowMultHerschelDecision_Dec, &b_muplus_Hlt1LowMultHerschelDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1LowMultHerschelDecision_TIS", &muplus_Hlt1LowMultHerschelDecision_TIS, &b_muplus_Hlt1LowMultHerschelDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1LowMultHerschelDecision_TOS", &muplus_Hlt1LowMultHerschelDecision_TOS, &b_muplus_Hlt1LowMultHerschelDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &muplus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_muplus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &muplus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_muplus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &muplus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_muplus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_muplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1TrackMVADecision_Dec", &muplus_Hlt1TrackMVADecision_Dec, &b_muplus_Hlt1TrackMVADecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1TrackMVADecision_TIS", &muplus_Hlt1TrackMVADecision_TIS, &b_muplus_Hlt1TrackMVADecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1TrackMVADecision_TOS", &muplus_Hlt1TrackMVADecision_TOS, &b_muplus_Hlt1TrackMVADecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1TwoTrackMVADecision_Dec", &muplus_Hlt1TwoTrackMVADecision_Dec, &b_muplus_Hlt1TwoTrackMVADecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1TwoTrackMVADecision_TIS", &muplus_Hlt1TwoTrackMVADecision_TIS, &b_muplus_Hlt1TwoTrackMVADecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1TwoTrackMVADecision_TOS", &muplus_Hlt1TwoTrackMVADecision_TOS, &b_muplus_Hlt1TwoTrackMVADecision_TOS);
	tree->SetBranchAddress("muplus_Hlt1SingleMuonHighPTDecision_Dec", &muplus_Hlt1SingleMuonHighPTDecision_Dec, &b_muplus_Hlt1SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt1SingleMuonHighPTDecision_TIS", &muplus_Hlt1SingleMuonHighPTDecision_TIS, &b_muplus_Hlt1SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt1SingleMuonHighPTDecision_TOS", &muplus_Hlt1SingleMuonHighPTDecision_TOS, &b_muplus_Hlt1SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHDecision_Dec", &muplus_Hlt2LowMultLMR2HHDecision_Dec, &b_muplus_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHDecision_TIS", &muplus_Hlt2LowMultLMR2HHDecision_TIS, &b_muplus_Hlt2LowMultLMR2HHDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHDecision_TOS", &muplus_Hlt2LowMultLMR2HHDecision_TOS, &b_muplus_Hlt2LowMultLMR2HHDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHWSDecision_Dec", &muplus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_muplus_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHWSDecision_TIS", &muplus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_muplus_Hlt2LowMultLMR2HHWSDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt2LowMultLMR2HHWSDecision_TOS", &muplus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_muplus_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("muplus_Hlt2SingleMuonHighPTDecision_Dec", &muplus_Hlt2SingleMuonHighPTDecision_Dec, &b_muplus_Hlt2SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("muplus_Hlt2SingleMuonHighPTDecision_TIS", &muplus_Hlt2SingleMuonHighPTDecision_TIS, &b_muplus_Hlt2SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("muplus_Hlt2SingleMuonHighPTDecision_TOS", &muplus_Hlt2SingleMuonHighPTDecision_TOS, &b_muplus_Hlt2SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("muplus_TRACK_Type", &muplus_TRACK_Type, &b_muplus_TRACK_Type);
	tree->SetBranchAddress("muplus_TRACK_Key", &muplus_TRACK_Key, &b_muplus_TRACK_Key);
	tree->SetBranchAddress("muplus_TRACK_CHI2NDOF", &muplus_TRACK_CHI2NDOF, &b_muplus_TRACK_CHI2NDOF);
	tree->SetBranchAddress("muplus_TRACK_PCHI2", &muplus_TRACK_PCHI2, &b_muplus_TRACK_PCHI2);
	tree->SetBranchAddress("muplus_TRACK_MatchCHI2", &muplus_TRACK_MatchCHI2, &b_muplus_TRACK_MatchCHI2);
	tree->SetBranchAddress("muplus_TRACK_GhostProb", &muplus_TRACK_GhostProb, &b_muplus_TRACK_GhostProb);
	tree->SetBranchAddress("muplus_TRACK_CloneDist", &muplus_TRACK_CloneDist, &b_muplus_TRACK_CloneDist);
	tree->SetBranchAddress("muplus_TRACK_Likelihood", &muplus_TRACK_Likelihood, &b_muplus_TRACK_Likelihood);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNe", &muminus_MC12TuneV2_ProbNNe, &b_muminus_MC12TuneV2_ProbNNe);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNmu", &muminus_MC12TuneV2_ProbNNmu, &b_muminus_MC12TuneV2_ProbNNmu);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNpi", &muminus_MC12TuneV2_ProbNNpi, &b_muminus_MC12TuneV2_ProbNNpi);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNk", &muminus_MC12TuneV2_ProbNNk, &b_muminus_MC12TuneV2_ProbNNk);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNp", &muminus_MC12TuneV2_ProbNNp, &b_muminus_MC12TuneV2_ProbNNp);
	tree->SetBranchAddress("muminus_MC12TuneV2_ProbNNghost", &muminus_MC12TuneV2_ProbNNghost, &b_muminus_MC12TuneV2_ProbNNghost);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNe", &muminus_MC12TuneV3_ProbNNe, &b_muminus_MC12TuneV3_ProbNNe);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNmu", &muminus_MC12TuneV3_ProbNNmu, &b_muminus_MC12TuneV3_ProbNNmu);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNpi", &muminus_MC12TuneV3_ProbNNpi, &b_muminus_MC12TuneV3_ProbNNpi);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNk", &muminus_MC12TuneV3_ProbNNk, &b_muminus_MC12TuneV3_ProbNNk);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNp", &muminus_MC12TuneV3_ProbNNp, &b_muminus_MC12TuneV3_ProbNNp);
	tree->SetBranchAddress("muminus_MC12TuneV3_ProbNNghost", &muminus_MC12TuneV3_ProbNNghost, &b_muminus_MC12TuneV3_ProbNNghost);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNe", &muminus_MC12TuneV4_ProbNNe, &b_muminus_MC12TuneV4_ProbNNe);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNmu", &muminus_MC12TuneV4_ProbNNmu, &b_muminus_MC12TuneV4_ProbNNmu);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNpi", &muminus_MC12TuneV4_ProbNNpi, &b_muminus_MC12TuneV4_ProbNNpi);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNk", &muminus_MC12TuneV4_ProbNNk, &b_muminus_MC12TuneV4_ProbNNk);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNp", &muminus_MC12TuneV4_ProbNNp, &b_muminus_MC12TuneV4_ProbNNp);
	tree->SetBranchAddress("muminus_MC12TuneV4_ProbNNghost", &muminus_MC12TuneV4_ProbNNghost, &b_muminus_MC12TuneV4_ProbNNghost);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNe", &muminus_MC15TuneV1_ProbNNe, &b_muminus_MC15TuneV1_ProbNNe);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNmu", &muminus_MC15TuneV1_ProbNNmu, &b_muminus_MC15TuneV1_ProbNNmu);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNpi", &muminus_MC15TuneV1_ProbNNpi, &b_muminus_MC15TuneV1_ProbNNpi);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNk", &muminus_MC15TuneV1_ProbNNk, &b_muminus_MC15TuneV1_ProbNNk);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNp", &muminus_MC15TuneV1_ProbNNp, &b_muminus_MC15TuneV1_ProbNNp);
	tree->SetBranchAddress("muminus_MC15TuneV1_ProbNNghost", &muminus_MC15TuneV1_ProbNNghost, &b_muminus_MC15TuneV1_ProbNNghost);
	tree->SetBranchAddress("muminus_OWNPV_X", &muminus_OWNPV_X, &b_muminus_OWNPV_X);
	tree->SetBranchAddress("muminus_OWNPV_Y", &muminus_OWNPV_Y, &b_muminus_OWNPV_Y);
	tree->SetBranchAddress("muminus_OWNPV_Z", &muminus_OWNPV_Z, &b_muminus_OWNPV_Z);
	tree->SetBranchAddress("muminus_OWNPV_XERR", &muminus_OWNPV_XERR, &b_muminus_OWNPV_XERR);
	tree->SetBranchAddress("muminus_OWNPV_YERR", &muminus_OWNPV_YERR, &b_muminus_OWNPV_YERR);
	tree->SetBranchAddress("muminus_OWNPV_ZERR", &muminus_OWNPV_ZERR, &b_muminus_OWNPV_ZERR);
	tree->SetBranchAddress("muminus_OWNPV_CHI2", &muminus_OWNPV_CHI2, &b_muminus_OWNPV_CHI2);
	tree->SetBranchAddress("muminus_OWNPV_NDOF", &muminus_OWNPV_NDOF, &b_muminus_OWNPV_NDOF);
	tree->SetBranchAddress("muminus_OWNPV_COV_", muminus_OWNPV_COV_, &b_muminus_OWNPV_COV_);
	tree->SetBranchAddress("muminus_IP_OWNPV", &muminus_IP_OWNPV, &b_muminus_IP_OWNPV);
	tree->SetBranchAddress("muminus_IPCHI2_OWNPV", &muminus_IPCHI2_OWNPV, &b_muminus_IPCHI2_OWNPV);
	tree->SetBranchAddress("muminus_ORIVX_X", &muminus_ORIVX_X, &b_muminus_ORIVX_X);
	tree->SetBranchAddress("muminus_ORIVX_Y", &muminus_ORIVX_Y, &b_muminus_ORIVX_Y);
	tree->SetBranchAddress("muminus_ORIVX_Z", &muminus_ORIVX_Z, &b_muminus_ORIVX_Z);
	tree->SetBranchAddress("muminus_ORIVX_XERR", &muminus_ORIVX_XERR, &b_muminus_ORIVX_XERR);
	tree->SetBranchAddress("muminus_ORIVX_YERR", &muminus_ORIVX_YERR, &b_muminus_ORIVX_YERR);
	tree->SetBranchAddress("muminus_ORIVX_ZERR", &muminus_ORIVX_ZERR, &b_muminus_ORIVX_ZERR);
	tree->SetBranchAddress("muminus_ORIVX_CHI2", &muminus_ORIVX_CHI2, &b_muminus_ORIVX_CHI2);
	tree->SetBranchAddress("muminus_ORIVX_NDOF", &muminus_ORIVX_NDOF, &b_muminus_ORIVX_NDOF);
	tree->SetBranchAddress("muminus_ORIVX_COV_", muminus_ORIVX_COV_, &b_muminus_ORIVX_COV_);
	tree->SetBranchAddress("muminus_P", &muminus_P, &b_muminus_P);
	tree->SetBranchAddress("muminus_PT", &muminus_PT, &b_muminus_PT);
	tree->SetBranchAddress("muminus_PE", &muminus_PE, &b_muminus_PE);
	tree->SetBranchAddress("muminus_PX", &muminus_PX, &b_muminus_PX);
	tree->SetBranchAddress("muminus_PY", &muminus_PY, &b_muminus_PY);
	tree->SetBranchAddress("muminus_PZ", &muminus_PZ, &b_muminus_PZ);
	tree->SetBranchAddress("muminus_M", &muminus_M, &b_muminus_M);
	tree->SetBranchAddress("muminus_ID", &muminus_ID, &b_muminus_ID);
	tree->SetBranchAddress("muminus_PIDe", &muminus_PIDe, &b_muminus_PIDe);
	tree->SetBranchAddress("muminus_PIDmu", &muminus_PIDmu, &b_muminus_PIDmu);
	tree->SetBranchAddress("muminus_PIDK", &muminus_PIDK, &b_muminus_PIDK);
	tree->SetBranchAddress("muminus_PIDp", &muminus_PIDp, &b_muminus_PIDp);
	tree->SetBranchAddress("muminus_PIDd", &muminus_PIDd, &b_muminus_PIDd);
	tree->SetBranchAddress("muminus_ProbNNe", &muminus_ProbNNe, &b_muminus_ProbNNe);
	tree->SetBranchAddress("muminus_ProbNNk", &muminus_ProbNNk, &b_muminus_ProbNNk);
	tree->SetBranchAddress("muminus_ProbNNp", &muminus_ProbNNp, &b_muminus_ProbNNp);
	tree->SetBranchAddress("muminus_ProbNNpi", &muminus_ProbNNpi, &b_muminus_ProbNNpi);
	tree->SetBranchAddress("muminus_ProbNNmu", &muminus_ProbNNmu, &b_muminus_ProbNNmu);
	tree->SetBranchAddress("muminus_ProbNNd", &muminus_ProbNNd, &b_muminus_ProbNNd);
	tree->SetBranchAddress("muminus_ProbNNghost", &muminus_ProbNNghost, &b_muminus_ProbNNghost);
	tree->SetBranchAddress("muminus_hasMuon", &muminus_hasMuon, &b_muminus_hasMuon);
	tree->SetBranchAddress("muminus_isMuon", &muminus_isMuon, &b_muminus_isMuon);
	tree->SetBranchAddress("muminus_hasRich", &muminus_hasRich, &b_muminus_hasRich);
	tree->SetBranchAddress("muminus_UsedRichAerogel", &muminus_UsedRichAerogel, &b_muminus_UsedRichAerogel);
	tree->SetBranchAddress("muminus_UsedRich1Gas", &muminus_UsedRich1Gas, &b_muminus_UsedRich1Gas);
	tree->SetBranchAddress("muminus_UsedRich2Gas", &muminus_UsedRich2Gas, &b_muminus_UsedRich2Gas);
	tree->SetBranchAddress("muminus_RichAboveElThres", &muminus_RichAboveElThres, &b_muminus_RichAboveElThres);
	tree->SetBranchAddress("muminus_RichAboveMuThres", &muminus_RichAboveMuThres, &b_muminus_RichAboveMuThres);
	tree->SetBranchAddress("muminus_RichAbovePiThres", &muminus_RichAbovePiThres, &b_muminus_RichAbovePiThres);
	tree->SetBranchAddress("muminus_RichAboveKaThres", &muminus_RichAboveKaThres, &b_muminus_RichAboveKaThres);
	tree->SetBranchAddress("muminus_RichAbovePrThres", &muminus_RichAbovePrThres, &b_muminus_RichAbovePrThres);
	tree->SetBranchAddress("muminus_hasCalo", &muminus_hasCalo, &b_muminus_hasCalo);
	tree->SetBranchAddress("muminus_PP_CombDLLe", &muminus_PP_CombDLLe, &b_muminus_PP_CombDLLe);
	tree->SetBranchAddress("muminus_PP_CombDLLmu", &muminus_PP_CombDLLmu, &b_muminus_PP_CombDLLmu);
	tree->SetBranchAddress("muminus_PP_CombDLLpi", &muminus_PP_CombDLLpi, &b_muminus_PP_CombDLLpi);
	tree->SetBranchAddress("muminus_PP_CombDLLk", &muminus_PP_CombDLLk, &b_muminus_PP_CombDLLk);
	tree->SetBranchAddress("muminus_PP_CombDLLp", &muminus_PP_CombDLLp, &b_muminus_PP_CombDLLp);
	tree->SetBranchAddress("muminus_PP_CombDLLd", &muminus_PP_CombDLLd, &b_muminus_PP_CombDLLd);
	tree->SetBranchAddress("muminus_PP_ProbNNe", &muminus_PP_ProbNNe, &b_muminus_PP_ProbNNe);
	tree->SetBranchAddress("muminus_PP_ProbNNmu", &muminus_PP_ProbNNmu, &b_muminus_PP_ProbNNmu);
	tree->SetBranchAddress("muminus_PP_ProbNNpi", &muminus_PP_ProbNNpi, &b_muminus_PP_ProbNNpi);
	tree->SetBranchAddress("muminus_PP_ProbNNk", &muminus_PP_ProbNNk, &b_muminus_PP_ProbNNk);
	tree->SetBranchAddress("muminus_PP_ProbNNp", &muminus_PP_ProbNNp, &b_muminus_PP_ProbNNp);
	tree->SetBranchAddress("muminus_PP_ProbNNghost", &muminus_PP_ProbNNghost, &b_muminus_PP_ProbNNghost);
	tree->SetBranchAddress("muminus_PP_ProbNNd", &muminus_PP_ProbNNd, &b_muminus_PP_ProbNNd);
	tree->SetBranchAddress("muminus_L0Global_Dec", &muminus_L0Global_Dec, &b_muminus_L0Global_Dec);
	tree->SetBranchAddress("muminus_L0Global_TIS", &muminus_L0Global_TIS, &b_muminus_L0Global_TIS);
	tree->SetBranchAddress("muminus_L0Global_TOS", &muminus_L0Global_TOS, &b_muminus_L0Global_TOS);
	tree->SetBranchAddress("muminus_Hlt1Global_Dec", &muminus_Hlt1Global_Dec, &b_muminus_Hlt1Global_Dec);
	tree->SetBranchAddress("muminus_Hlt1Global_TIS", &muminus_Hlt1Global_TIS, &b_muminus_Hlt1Global_TIS);
	tree->SetBranchAddress("muminus_Hlt1Global_TOS", &muminus_Hlt1Global_TOS, &b_muminus_Hlt1Global_TOS);
	tree->SetBranchAddress("muminus_Hlt1Phys_Dec", &muminus_Hlt1Phys_Dec, &b_muminus_Hlt1Phys_Dec);
	tree->SetBranchAddress("muminus_Hlt1Phys_TIS", &muminus_Hlt1Phys_TIS, &b_muminus_Hlt1Phys_TIS);
	tree->SetBranchAddress("muminus_Hlt1Phys_TOS", &muminus_Hlt1Phys_TOS, &b_muminus_Hlt1Phys_TOS);
	tree->SetBranchAddress("muminus_Hlt2Global_Dec", &muminus_Hlt2Global_Dec, &b_muminus_Hlt2Global_Dec);
	tree->SetBranchAddress("muminus_Hlt2Global_TIS", &muminus_Hlt2Global_TIS, &b_muminus_Hlt2Global_TIS);
	tree->SetBranchAddress("muminus_Hlt2Global_TOS", &muminus_Hlt2Global_TOS, &b_muminus_Hlt2Global_TOS);
	tree->SetBranchAddress("muminus_Hlt2Phys_Dec", &muminus_Hlt2Phys_Dec, &b_muminus_Hlt2Phys_Dec);
	tree->SetBranchAddress("muminus_Hlt2Phys_TIS", &muminus_Hlt2Phys_TIS, &b_muminus_Hlt2Phys_TIS);
	tree->SetBranchAddress("muminus_Hlt2Phys_TOS", &muminus_Hlt2Phys_TOS, &b_muminus_Hlt2Phys_TOS);
	tree->SetBranchAddress("muminus_L0DiHadron,lowMultDecision_Dec", &muminus_L0DiHadron_lowMultDecision_Dec, &b_muminus_L0DiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("muminus_L0DiHadron,lowMultDecision_TIS", &muminus_L0DiHadron_lowMultDecision_TIS, &b_muminus_L0DiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("muminus_L0DiHadron,lowMultDecision_TOS", &muminus_L0DiHadron_lowMultDecision_TOS, &b_muminus_L0DiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("muminus_L0HRCDiHadron,lowMultDecision_Dec", &muminus_L0HRCDiHadron_lowMultDecision_Dec, &b_muminus_L0HRCDiHadron_lowMultDecision_Dec);
	tree->SetBranchAddress("muminus_L0HRCDiHadron,lowMultDecision_TIS", &muminus_L0HRCDiHadron_lowMultDecision_TIS, &b_muminus_L0HRCDiHadron_lowMultDecision_TIS);
	tree->SetBranchAddress("muminus_L0HRCDiHadron,lowMultDecision_TOS", &muminus_L0HRCDiHadron_lowMultDecision_TOS, &b_muminus_L0HRCDiHadron_lowMultDecision_TOS);
	tree->SetBranchAddress("muminus_L0MuonDecision_Dec", &muminus_L0MuonDecision_Dec, &b_muminus_L0MuonDecision_Dec);
	tree->SetBranchAddress("muminus_L0MuonDecision_TIS", &muminus_L0MuonDecision_TIS, &b_muminus_L0MuonDecision_TIS);
	tree->SetBranchAddress("muminus_L0MuonDecision_TOS", &muminus_L0MuonDecision_TOS, &b_muminus_L0MuonDecision_TOS);
	tree->SetBranchAddress("muminus_L0MuonEWDecision_Dec", &muminus_L0MuonEWDecision_Dec, &b_muminus_L0MuonEWDecision_Dec);
	tree->SetBranchAddress("muminus_L0MuonEWDecision_TIS", &muminus_L0MuonEWDecision_TIS, &b_muminus_L0MuonEWDecision_TIS);
	tree->SetBranchAddress("muminus_L0MuonEWDecision_TOS", &muminus_L0MuonEWDecision_TOS, &b_muminus_L0MuonEWDecision_TOS);
	tree->SetBranchAddress("muminus_L0DiMuonDecision_Dec", &muminus_L0DiMuonDecision_Dec, &b_muminus_L0DiMuonDecision_Dec);
	tree->SetBranchAddress("muminus_L0DiMuonDecision_TIS", &muminus_L0DiMuonDecision_TIS, &b_muminus_L0DiMuonDecision_TIS);
	tree->SetBranchAddress("muminus_L0DiMuonDecision_TOS", &muminus_L0DiMuonDecision_TOS, &b_muminus_L0DiMuonDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_muminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1LowMultHerschelDecision_Dec", &muminus_Hlt1LowMultHerschelDecision_Dec, &b_muminus_Hlt1LowMultHerschelDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1LowMultHerschelDecision_TIS", &muminus_Hlt1LowMultHerschelDecision_TIS, &b_muminus_Hlt1LowMultHerschelDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1LowMultHerschelDecision_TOS", &muminus_Hlt1LowMultHerschelDecision_TOS, &b_muminus_Hlt1LowMultHerschelDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &muminus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_muminus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &muminus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_muminus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &muminus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_muminus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_muminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1TrackMVADecision_Dec", &muminus_Hlt1TrackMVADecision_Dec, &b_muminus_Hlt1TrackMVADecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1TrackMVADecision_TIS", &muminus_Hlt1TrackMVADecision_TIS, &b_muminus_Hlt1TrackMVADecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1TrackMVADecision_TOS", &muminus_Hlt1TrackMVADecision_TOS, &b_muminus_Hlt1TrackMVADecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1TwoTrackMVADecision_Dec", &muminus_Hlt1TwoTrackMVADecision_Dec, &b_muminus_Hlt1TwoTrackMVADecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1TwoTrackMVADecision_TIS", &muminus_Hlt1TwoTrackMVADecision_TIS, &b_muminus_Hlt1TwoTrackMVADecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1TwoTrackMVADecision_TOS", &muminus_Hlt1TwoTrackMVADecision_TOS, &b_muminus_Hlt1TwoTrackMVADecision_TOS);
	tree->SetBranchAddress("muminus_Hlt1SingleMuonHighPTDecision_Dec", &muminus_Hlt1SingleMuonHighPTDecision_Dec, &b_muminus_Hlt1SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt1SingleMuonHighPTDecision_TIS", &muminus_Hlt1SingleMuonHighPTDecision_TIS, &b_muminus_Hlt1SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt1SingleMuonHighPTDecision_TOS", &muminus_Hlt1SingleMuonHighPTDecision_TOS, &b_muminus_Hlt1SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHDecision_Dec", &muminus_Hlt2LowMultLMR2HHDecision_Dec, &b_muminus_Hlt2LowMultLMR2HHDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHDecision_TIS", &muminus_Hlt2LowMultLMR2HHDecision_TIS, &b_muminus_Hlt2LowMultLMR2HHDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHDecision_TOS", &muminus_Hlt2LowMultLMR2HHDecision_TOS, &b_muminus_Hlt2LowMultLMR2HHDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHWSDecision_Dec", &muminus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_muminus_Hlt2LowMultLMR2HHWSDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHWSDecision_TIS", &muminus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_muminus_Hlt2LowMultLMR2HHWSDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt2LowMultLMR2HHWSDecision_TOS", &muminus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_muminus_Hlt2LowMultLMR2HHWSDecision_TOS);
	tree->SetBranchAddress("muminus_Hlt2SingleMuonHighPTDecision_Dec", &muminus_Hlt2SingleMuonHighPTDecision_Dec, &b_muminus_Hlt2SingleMuonHighPTDecision_Dec);
	tree->SetBranchAddress("muminus_Hlt2SingleMuonHighPTDecision_TIS", &muminus_Hlt2SingleMuonHighPTDecision_TIS, &b_muminus_Hlt2SingleMuonHighPTDecision_TIS);
	tree->SetBranchAddress("muminus_Hlt2SingleMuonHighPTDecision_TOS", &muminus_Hlt2SingleMuonHighPTDecision_TOS, &b_muminus_Hlt2SingleMuonHighPTDecision_TOS);
	tree->SetBranchAddress("muminus_TRACK_Type", &muminus_TRACK_Type, &b_muminus_TRACK_Type);
	tree->SetBranchAddress("muminus_TRACK_Key", &muminus_TRACK_Key, &b_muminus_TRACK_Key);
	tree->SetBranchAddress("muminus_TRACK_CHI2NDOF", &muminus_TRACK_CHI2NDOF, &b_muminus_TRACK_CHI2NDOF);
	tree->SetBranchAddress("muminus_TRACK_PCHI2", &muminus_TRACK_PCHI2, &b_muminus_TRACK_PCHI2);
	tree->SetBranchAddress("muminus_TRACK_MatchCHI2", &muminus_TRACK_MatchCHI2, &b_muminus_TRACK_MatchCHI2);
	tree->SetBranchAddress("muminus_TRACK_GhostProb", &muminus_TRACK_GhostProb, &b_muminus_TRACK_GhostProb);
	tree->SetBranchAddress("muminus_TRACK_CloneDist", &muminus_TRACK_CloneDist, &b_muminus_TRACK_CloneDist);
	tree->SetBranchAddress("muminus_TRACK_Likelihood", &muminus_TRACK_Likelihood, &b_muminus_TRACK_Likelihood);
	tree->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
	tree->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
	tree->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
	tree->SetBranchAddress("LoKi_nCharged", &LoKi_nCharged, &b_LoKi_nCharged);
	tree->SetBranchAddress("LoKi_nITClusters", &LoKi_nITClusters, &b_LoKi_nITClusters);
	tree->SetBranchAddress("LoKi_nNeutrals", &LoKi_nNeutrals, &b_LoKi_nNeutrals);
	tree->SetBranchAddress("LoKi_nOThits", &LoKi_nOThits, &b_LoKi_nOThits);
	tree->SetBranchAddress("LoKi_nPVs", &LoKi_nPVs, &b_LoKi_nPVs);
	tree->SetBranchAddress("LoKi_nSpdMult", &LoKi_nSpdMult, &b_LoKi_nSpdMult);
	tree->SetBranchAddress("LoKi_nTTClusters", &LoKi_nTTClusters, &b_LoKi_nTTClusters);
	tree->SetBranchAddress("LoKi_nVeloClusters", &LoKi_nVeloClusters, &b_LoKi_nVeloClusters);
	tree->SetBranchAddress("LoKi_nVeloLiteClusters", &LoKi_nVeloLiteClusters, &b_LoKi_nVeloLiteClusters);
	tree->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	tree->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
	tree->SetBranchAddress("BCID", &BCID, &b_BCID);
	tree->SetBranchAddress("BCType", &BCType, &b_BCType);
	tree->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
	tree->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
	tree->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
	tree->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
	tree->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
	tree->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
	tree->SetBranchAddress("B00", &B00, &b_B00);
	tree->SetBranchAddress("B01", &B01, &b_B01);
	tree->SetBranchAddress("B02", &B02, &b_B02);
	tree->SetBranchAddress("B03", &B03, &b_B03);
	tree->SetBranchAddress("B10", &B10, &b_B10);
	tree->SetBranchAddress("B11", &B11, &b_B11);
	tree->SetBranchAddress("B12", &B12, &b_B12);
	tree->SetBranchAddress("B13", &B13, &b_B13);
	tree->SetBranchAddress("B20", &B20, &b_B20);
	tree->SetBranchAddress("B21", &B21, &b_B21);
	tree->SetBranchAddress("B22", &B22, &b_B22);
	tree->SetBranchAddress("B23", &B23, &b_B23);
	tree->SetBranchAddress("F10", &F10, &b_F10);
	tree->SetBranchAddress("F11", &F11, &b_F11);
	tree->SetBranchAddress("F12", &F12, &b_F12);
	tree->SetBranchAddress("F13", &F13, &b_F13);
	tree->SetBranchAddress("F20", &F20, &b_F20);
	tree->SetBranchAddress("F21", &F21, &b_F21);
	tree->SetBranchAddress("F22", &F22, &b_F22);
	tree->SetBranchAddress("F23", &F23, &b_F23);
	tree->SetBranchAddress("log_hrc_fom_v3", &log_hrc_fom_v3, &b_log_hrc_fom_v3);
	tree->SetBranchAddress("log_hrc_fom_B_v3", &log_hrc_fom_B_v3, &b_log_hrc_fom_B_v3);
	tree->SetBranchAddress("log_hrc_fom_F_v3", &log_hrc_fom_F_v3, &b_log_hrc_fom_F_v3);
	tree->SetBranchAddress("nchB", &nchB, &b_nchB);
	tree->SetBranchAddress("adc_B", adc_B, &b_adc_B);
	tree->SetBranchAddress("nchF", &nchF, &b_nchF);
	tree->SetBranchAddress("adc_F", adc_F, &b_adc_F);
	tree->SetBranchAddress("nPV", &nPV, &b_nPV);
	tree->SetBranchAddress("PVX", PVX, &b_PVX);
	tree->SetBranchAddress("PVY", PVY, &b_PVY);
	tree->SetBranchAddress("PVZ", PVZ, &b_PVZ);
	tree->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
	tree->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
	tree->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
	tree->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
	tree->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
	tree->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
	tree->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
	tree->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	tree->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
	tree->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
	tree->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
	tree->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
	tree->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
	tree->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
	tree->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
	tree->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
	tree->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
	tree->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
	tree->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
	tree->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
	tree->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
	tree->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
	tree->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
	tree->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
	tree->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
	tree->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
	tree->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
}

#endif // #ifdef skim_cxx
