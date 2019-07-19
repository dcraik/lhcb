//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 12 19:16:36 2018 by ROOT version 5.34/38
// from TTree DecayTree/DecayTree
// found on file: /tmp/dcraik/Tuples.root
//////////////////////////////////////////////////////////

#ifndef skim_h
#define skim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>

#include <boost/progress.hpp>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxphi_ENDVERTEX_COV = 1;
const Int_t kMaxphi_OWNPV_COV = 1;
const Int_t kMaxKplus_OWNPV_COV = 1;
const Int_t kMaxKplus_ORIVX_COV = 1;
const Int_t kMaxKminus_OWNPV_COV = 1;
const Int_t kMaxKminus_ORIVX_COV = 1;

class skim {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

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
		Bool_t          phi_L0DiHadron_lowMultDecision_Dec;
		Bool_t          phi_L0DiHadron_lowMultDecision_TIS;
		Bool_t          phi_L0DiHadron_lowMultDecision_TOS;
		Bool_t          phi_L0HRCDiHadron_lowMultDecision_Dec;
		Bool_t          phi_L0HRCDiHadron_lowMultDecision_TIS;
		Bool_t          phi_L0HRCDiHadron_lowMultDecision_TOS;
		Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;
		Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;
		Bool_t          phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;
		Bool_t          phi_Hlt1LowMultHerschelDecision_Dec;
		Bool_t          phi_Hlt1LowMultHerschelDecision_TIS;
		Bool_t          phi_Hlt1LowMultHerschelDecision_TOS;
		Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;
		Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_TIS;
		Bool_t          phi_Hlt1LowMultVeloCut_HadronsDecision_TOS;
		Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;
		Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;
		Bool_t          phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;
		Bool_t          phi_Hlt2LowMultLMR2HHDecision_Dec;
		Bool_t          phi_Hlt2LowMultLMR2HHDecision_TIS;
		Bool_t          phi_Hlt2LowMultLMR2HHDecision_TOS;
		Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_Dec;
		Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_TIS;
		Bool_t          phi_Hlt2LowMultLMR2HHWSDecision_TOS;
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
		TBranch        *b_phi_L0DiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_phi_L0DiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_phi_L0DiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_Dec;   //!
		TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_TIS;   //!
		TBranch        *b_phi_L0HRCDiHadron_lowMultDecision_TOS;   //!
		TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec;   //!
		TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS;   //!
		TBranch        *b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS;   //!
		TBranch        *b_phi_Hlt1LowMultHerschelDecision_Dec;   //!
		TBranch        *b_phi_Hlt1LowMultHerschelDecision_TIS;   //!
		TBranch        *b_phi_Hlt1LowMultHerschelDecision_TOS;   //!
		TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec;   //!
		TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_TIS;   //!
		TBranch        *b_phi_Hlt1LowMultVeloCut_HadronsDecision_TOS;   //!
		TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec;   //!
		TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS;   //!
		TBranch        *b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_Dec;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_TIS;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHDecision_TOS;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_Dec;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_TIS;   //!
		TBranch        *b_phi_Hlt2LowMultLMR2HHWSDecision_TOS;   //!
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

		skim(int job, int sjob=-1, TString dir="/tmp/dcraik");
		virtual ~skim();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		TString outName;

		TChain* lumi;
};

#endif

#ifdef skim_cxx
skim::skim(int job, int sjob, TString dir) : fChain(0) 
{
	outName=dir;
	outName+="/";
	outName+=job;
	outName+="/";
	outName+=sjob;
	outName+="/skimmed.root";

	//outName="/tmp/dcraik/skimmed.root";

	char str[256];
	TChain* t(0);
	t = new TChain("Tuple/DecayTree");
	lumi = new TChain("GetIntegratedLuminosity/LumiTuple");
	//t->Add("/tmp/dcraik/Tuples.root");
	//lumi->Add("/tmp/dcraik/Tuples.root");
	boost::progress_display show_addfile_progress( 700 );
	for(int i=0; i<700; ++i) {
		++show_addfile_progress;
		if(sjob<0 || sjob==i) {
			sprintf(str,"%s/%d/%d/Tuples.root",dir.Data(),job,i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
			lumi->Add(str);
		}
	}

	Init(t);

	std::cout << t->GetEntries() << std::endl;
}

skim::~skim()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t skim::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t skim::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void skim::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("phi_ENDVERTEX_X", &phi_ENDVERTEX_X, &b_phi_ENDVERTEX_X);
	fChain->SetBranchAddress("phi_ENDVERTEX_Y", &phi_ENDVERTEX_Y, &b_phi_ENDVERTEX_Y);
	fChain->SetBranchAddress("phi_ENDVERTEX_Z", &phi_ENDVERTEX_Z, &b_phi_ENDVERTEX_Z);
	fChain->SetBranchAddress("phi_ENDVERTEX_XERR", &phi_ENDVERTEX_XERR, &b_phi_ENDVERTEX_XERR);
	fChain->SetBranchAddress("phi_ENDVERTEX_YERR", &phi_ENDVERTEX_YERR, &b_phi_ENDVERTEX_YERR);
	fChain->SetBranchAddress("phi_ENDVERTEX_ZERR", &phi_ENDVERTEX_ZERR, &b_phi_ENDVERTEX_ZERR);
	fChain->SetBranchAddress("phi_ENDVERTEX_CHI2", &phi_ENDVERTEX_CHI2, &b_phi_ENDVERTEX_CHI2);
	fChain->SetBranchAddress("phi_ENDVERTEX_NDOF", &phi_ENDVERTEX_NDOF, &b_phi_ENDVERTEX_NDOF);
	fChain->SetBranchAddress("phi_ENDVERTEX_COV_", phi_ENDVERTEX_COV_, &b_phi_ENDVERTEX_COV_);
	fChain->SetBranchAddress("phi_OWNPV_X", &phi_OWNPV_X, &b_phi_OWNPV_X);
	fChain->SetBranchAddress("phi_OWNPV_Y", &phi_OWNPV_Y, &b_phi_OWNPV_Y);
	fChain->SetBranchAddress("phi_OWNPV_Z", &phi_OWNPV_Z, &b_phi_OWNPV_Z);
	fChain->SetBranchAddress("phi_OWNPV_XERR", &phi_OWNPV_XERR, &b_phi_OWNPV_XERR);
	fChain->SetBranchAddress("phi_OWNPV_YERR", &phi_OWNPV_YERR, &b_phi_OWNPV_YERR);
	fChain->SetBranchAddress("phi_OWNPV_ZERR", &phi_OWNPV_ZERR, &b_phi_OWNPV_ZERR);
	fChain->SetBranchAddress("phi_OWNPV_CHI2", &phi_OWNPV_CHI2, &b_phi_OWNPV_CHI2);
	fChain->SetBranchAddress("phi_OWNPV_NDOF", &phi_OWNPV_NDOF, &b_phi_OWNPV_NDOF);
	fChain->SetBranchAddress("phi_OWNPV_COV_", phi_OWNPV_COV_, &b_phi_OWNPV_COV_);
	fChain->SetBranchAddress("phi_IP_OWNPV", &phi_IP_OWNPV, &b_phi_IP_OWNPV);
	fChain->SetBranchAddress("phi_IPCHI2_OWNPV", &phi_IPCHI2_OWNPV, &b_phi_IPCHI2_OWNPV);
	fChain->SetBranchAddress("phi_FD_OWNPV", &phi_FD_OWNPV, &b_phi_FD_OWNPV);
	fChain->SetBranchAddress("phi_FDCHI2_OWNPV", &phi_FDCHI2_OWNPV, &b_phi_FDCHI2_OWNPV);
	fChain->SetBranchAddress("phi_DIRA_OWNPV", &phi_DIRA_OWNPV, &b_phi_DIRA_OWNPV);
	fChain->SetBranchAddress("phi_P", &phi_P, &b_phi_P);
	fChain->SetBranchAddress("phi_PT", &phi_PT, &b_phi_PT);
	fChain->SetBranchAddress("phi_PE", &phi_PE, &b_phi_PE);
	fChain->SetBranchAddress("phi_PX", &phi_PX, &b_phi_PX);
	fChain->SetBranchAddress("phi_PY", &phi_PY, &b_phi_PY);
	fChain->SetBranchAddress("phi_PZ", &phi_PZ, &b_phi_PZ);
	fChain->SetBranchAddress("phi_MM", &phi_MM, &b_phi_MM);
	fChain->SetBranchAddress("phi_MMERR", &phi_MMERR, &b_phi_MMERR);
	fChain->SetBranchAddress("phi_M", &phi_M, &b_phi_M);
	fChain->SetBranchAddress("phi_ID", &phi_ID, &b_phi_ID);
	fChain->SetBranchAddress("phi_L0Global_Dec", &phi_L0Global_Dec, &b_phi_L0Global_Dec);
	fChain->SetBranchAddress("phi_L0Global_TIS", &phi_L0Global_TIS, &b_phi_L0Global_TIS);
	fChain->SetBranchAddress("phi_L0Global_TOS", &phi_L0Global_TOS, &b_phi_L0Global_TOS);
	fChain->SetBranchAddress("phi_Hlt1Global_Dec", &phi_Hlt1Global_Dec, &b_phi_Hlt1Global_Dec);
	fChain->SetBranchAddress("phi_Hlt1Global_TIS", &phi_Hlt1Global_TIS, &b_phi_Hlt1Global_TIS);
	fChain->SetBranchAddress("phi_Hlt1Global_TOS", &phi_Hlt1Global_TOS, &b_phi_Hlt1Global_TOS);
	fChain->SetBranchAddress("phi_Hlt1Phys_Dec", &phi_Hlt1Phys_Dec, &b_phi_Hlt1Phys_Dec);
	fChain->SetBranchAddress("phi_Hlt1Phys_TIS", &phi_Hlt1Phys_TIS, &b_phi_Hlt1Phys_TIS);
	fChain->SetBranchAddress("phi_Hlt1Phys_TOS", &phi_Hlt1Phys_TOS, &b_phi_Hlt1Phys_TOS);
	fChain->SetBranchAddress("phi_Hlt2Global_Dec", &phi_Hlt2Global_Dec, &b_phi_Hlt2Global_Dec);
	fChain->SetBranchAddress("phi_Hlt2Global_TIS", &phi_Hlt2Global_TIS, &b_phi_Hlt2Global_TIS);
	fChain->SetBranchAddress("phi_Hlt2Global_TOS", &phi_Hlt2Global_TOS, &b_phi_Hlt2Global_TOS);
	fChain->SetBranchAddress("phi_Hlt2Phys_Dec", &phi_Hlt2Phys_Dec, &b_phi_Hlt2Phys_Dec);
	fChain->SetBranchAddress("phi_Hlt2Phys_TIS", &phi_Hlt2Phys_TIS, &b_phi_Hlt2Phys_TIS);
	fChain->SetBranchAddress("phi_Hlt2Phys_TOS", &phi_Hlt2Phys_TOS, &b_phi_Hlt2Phys_TOS);
	fChain->SetBranchAddress("phi_L0DiHadron,lowMultDecision_Dec", &phi_L0DiHadron_lowMultDecision_Dec, &b_phi_L0DiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("phi_L0DiHadron,lowMultDecision_TIS", &phi_L0DiHadron_lowMultDecision_TIS, &b_phi_L0DiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("phi_L0DiHadron,lowMultDecision_TOS", &phi_L0DiHadron_lowMultDecision_TOS, &b_phi_L0DiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_Dec", &phi_L0HRCDiHadron_lowMultDecision_Dec, &b_phi_L0HRCDiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_TIS", &phi_L0HRCDiHadron_lowMultDecision_TIS, &b_phi_L0HRCDiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("phi_L0HRCDiHadron,lowMultDecision_TOS", &phi_L0HRCDiHadron_lowMultDecision_TOS, &b_phi_L0HRCDiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_phi_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_Dec", &phi_Hlt1LowMultHerschelDecision_Dec, &b_phi_Hlt1LowMultHerschelDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_TIS", &phi_Hlt1LowMultHerschelDecision_TIS, &b_phi_Hlt1LowMultHerschelDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt1LowMultHerschelDecision_TOS", &phi_Hlt1LowMultHerschelDecision_TOS, &b_phi_Hlt1LowMultHerschelDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_Dec", &phi_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_TIS", &phi_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloCut_HadronsDecision_TOS", &phi_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_phi_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_phi_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_Dec", &phi_Hlt2LowMultLMR2HHDecision_Dec, &b_phi_Hlt2LowMultLMR2HHDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_TIS", &phi_Hlt2LowMultLMR2HHDecision_TIS, &b_phi_Hlt2LowMultLMR2HHDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHDecision_TOS", &phi_Hlt2LowMultLMR2HHDecision_TOS, &b_phi_Hlt2LowMultLMR2HHDecision_TOS);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_Dec", &phi_Hlt2LowMultLMR2HHWSDecision_Dec, &b_phi_Hlt2LowMultLMR2HHWSDecision_Dec);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_TIS", &phi_Hlt2LowMultLMR2HHWSDecision_TIS, &b_phi_Hlt2LowMultLMR2HHWSDecision_TIS);
	fChain->SetBranchAddress("phi_Hlt2LowMultLMR2HHWSDecision_TOS", &phi_Hlt2LowMultLMR2HHWSDecision_TOS, &b_phi_Hlt2LowMultLMR2HHWSDecision_TOS);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNe", &Kplus_MC12TuneV2_ProbNNe, &b_Kplus_MC12TuneV2_ProbNNe);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNmu", &Kplus_MC12TuneV2_ProbNNmu, &b_Kplus_MC12TuneV2_ProbNNmu);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNpi", &Kplus_MC12TuneV2_ProbNNpi, &b_Kplus_MC12TuneV2_ProbNNpi);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNk", &Kplus_MC12TuneV2_ProbNNk, &b_Kplus_MC12TuneV2_ProbNNk);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNp", &Kplus_MC12TuneV2_ProbNNp, &b_Kplus_MC12TuneV2_ProbNNp);
	fChain->SetBranchAddress("Kplus_MC12TuneV2_ProbNNghost", &Kplus_MC12TuneV2_ProbNNghost, &b_Kplus_MC12TuneV2_ProbNNghost);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNe", &Kplus_MC12TuneV3_ProbNNe, &b_Kplus_MC12TuneV3_ProbNNe);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNmu", &Kplus_MC12TuneV3_ProbNNmu, &b_Kplus_MC12TuneV3_ProbNNmu);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNpi", &Kplus_MC12TuneV3_ProbNNpi, &b_Kplus_MC12TuneV3_ProbNNpi);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNk", &Kplus_MC12TuneV3_ProbNNk, &b_Kplus_MC12TuneV3_ProbNNk);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNp", &Kplus_MC12TuneV3_ProbNNp, &b_Kplus_MC12TuneV3_ProbNNp);
	fChain->SetBranchAddress("Kplus_MC12TuneV3_ProbNNghost", &Kplus_MC12TuneV3_ProbNNghost, &b_Kplus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNe", &Kplus_MC12TuneV4_ProbNNe, &b_Kplus_MC12TuneV4_ProbNNe);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNmu", &Kplus_MC12TuneV4_ProbNNmu, &b_Kplus_MC12TuneV4_ProbNNmu);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNpi", &Kplus_MC12TuneV4_ProbNNpi, &b_Kplus_MC12TuneV4_ProbNNpi);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNk", &Kplus_MC12TuneV4_ProbNNk, &b_Kplus_MC12TuneV4_ProbNNk);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNp", &Kplus_MC12TuneV4_ProbNNp, &b_Kplus_MC12TuneV4_ProbNNp);
	fChain->SetBranchAddress("Kplus_MC12TuneV4_ProbNNghost", &Kplus_MC12TuneV4_ProbNNghost, &b_Kplus_MC12TuneV4_ProbNNghost);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNe", &Kplus_MC15TuneV1_ProbNNe, &b_Kplus_MC15TuneV1_ProbNNe);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNmu", &Kplus_MC15TuneV1_ProbNNmu, &b_Kplus_MC15TuneV1_ProbNNmu);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNpi", &Kplus_MC15TuneV1_ProbNNpi, &b_Kplus_MC15TuneV1_ProbNNpi);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNk", &Kplus_MC15TuneV1_ProbNNk, &b_Kplus_MC15TuneV1_ProbNNk);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNp", &Kplus_MC15TuneV1_ProbNNp, &b_Kplus_MC15TuneV1_ProbNNp);
	fChain->SetBranchAddress("Kplus_MC15TuneV1_ProbNNghost", &Kplus_MC15TuneV1_ProbNNghost, &b_Kplus_MC15TuneV1_ProbNNghost);
	fChain->SetBranchAddress("Kplus_OWNPV_X", &Kplus_OWNPV_X, &b_Kplus_OWNPV_X);
	fChain->SetBranchAddress("Kplus_OWNPV_Y", &Kplus_OWNPV_Y, &b_Kplus_OWNPV_Y);
	fChain->SetBranchAddress("Kplus_OWNPV_Z", &Kplus_OWNPV_Z, &b_Kplus_OWNPV_Z);
	fChain->SetBranchAddress("Kplus_OWNPV_XERR", &Kplus_OWNPV_XERR, &b_Kplus_OWNPV_XERR);
	fChain->SetBranchAddress("Kplus_OWNPV_YERR", &Kplus_OWNPV_YERR, &b_Kplus_OWNPV_YERR);
	fChain->SetBranchAddress("Kplus_OWNPV_ZERR", &Kplus_OWNPV_ZERR, &b_Kplus_OWNPV_ZERR);
	fChain->SetBranchAddress("Kplus_OWNPV_CHI2", &Kplus_OWNPV_CHI2, &b_Kplus_OWNPV_CHI2);
	fChain->SetBranchAddress("Kplus_OWNPV_NDOF", &Kplus_OWNPV_NDOF, &b_Kplus_OWNPV_NDOF);
	fChain->SetBranchAddress("Kplus_OWNPV_COV_", Kplus_OWNPV_COV_, &b_Kplus_OWNPV_COV_);
	fChain->SetBranchAddress("Kplus_IP_OWNPV", &Kplus_IP_OWNPV, &b_Kplus_IP_OWNPV);
	fChain->SetBranchAddress("Kplus_IPCHI2_OWNPV", &Kplus_IPCHI2_OWNPV, &b_Kplus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("Kplus_ORIVX_X", &Kplus_ORIVX_X, &b_Kplus_ORIVX_X);
	fChain->SetBranchAddress("Kplus_ORIVX_Y", &Kplus_ORIVX_Y, &b_Kplus_ORIVX_Y);
	fChain->SetBranchAddress("Kplus_ORIVX_Z", &Kplus_ORIVX_Z, &b_Kplus_ORIVX_Z);
	fChain->SetBranchAddress("Kplus_ORIVX_XERR", &Kplus_ORIVX_XERR, &b_Kplus_ORIVX_XERR);
	fChain->SetBranchAddress("Kplus_ORIVX_YERR", &Kplus_ORIVX_YERR, &b_Kplus_ORIVX_YERR);
	fChain->SetBranchAddress("Kplus_ORIVX_ZERR", &Kplus_ORIVX_ZERR, &b_Kplus_ORIVX_ZERR);
	fChain->SetBranchAddress("Kplus_ORIVX_CHI2", &Kplus_ORIVX_CHI2, &b_Kplus_ORIVX_CHI2);
	fChain->SetBranchAddress("Kplus_ORIVX_NDOF", &Kplus_ORIVX_NDOF, &b_Kplus_ORIVX_NDOF);
	fChain->SetBranchAddress("Kplus_ORIVX_COV_", Kplus_ORIVX_COV_, &b_Kplus_ORIVX_COV_);
	fChain->SetBranchAddress("Kplus_P", &Kplus_P, &b_Kplus_P);
	fChain->SetBranchAddress("Kplus_PT", &Kplus_PT, &b_Kplus_PT);
	fChain->SetBranchAddress("Kplus_PE", &Kplus_PE, &b_Kplus_PE);
	fChain->SetBranchAddress("Kplus_PX", &Kplus_PX, &b_Kplus_PX);
	fChain->SetBranchAddress("Kplus_PY", &Kplus_PY, &b_Kplus_PY);
	fChain->SetBranchAddress("Kplus_PZ", &Kplus_PZ, &b_Kplus_PZ);
	fChain->SetBranchAddress("Kplus_M", &Kplus_M, &b_Kplus_M);
	fChain->SetBranchAddress("Kplus_ID", &Kplus_ID, &b_Kplus_ID);
	fChain->SetBranchAddress("Kplus_PIDe", &Kplus_PIDe, &b_Kplus_PIDe);
	fChain->SetBranchAddress("Kplus_PIDmu", &Kplus_PIDmu, &b_Kplus_PIDmu);
	fChain->SetBranchAddress("Kplus_PIDK", &Kplus_PIDK, &b_Kplus_PIDK);
	fChain->SetBranchAddress("Kplus_PIDp", &Kplus_PIDp, &b_Kplus_PIDp);
	fChain->SetBranchAddress("Kplus_ProbNNe", &Kplus_ProbNNe, &b_Kplus_ProbNNe);
	fChain->SetBranchAddress("Kplus_ProbNNk", &Kplus_ProbNNk, &b_Kplus_ProbNNk);
	fChain->SetBranchAddress("Kplus_ProbNNp", &Kplus_ProbNNp, &b_Kplus_ProbNNp);
	fChain->SetBranchAddress("Kplus_ProbNNpi", &Kplus_ProbNNpi, &b_Kplus_ProbNNpi);
	fChain->SetBranchAddress("Kplus_ProbNNmu", &Kplus_ProbNNmu, &b_Kplus_ProbNNmu);
	fChain->SetBranchAddress("Kplus_ProbNNghost", &Kplus_ProbNNghost, &b_Kplus_ProbNNghost);
	fChain->SetBranchAddress("Kplus_hasMuon", &Kplus_hasMuon, &b_Kplus_hasMuon);
	fChain->SetBranchAddress("Kplus_isMuon", &Kplus_isMuon, &b_Kplus_isMuon);
	fChain->SetBranchAddress("Kplus_hasRich", &Kplus_hasRich, &b_Kplus_hasRich);
	fChain->SetBranchAddress("Kplus_UsedRichAerogel", &Kplus_UsedRichAerogel, &b_Kplus_UsedRichAerogel);
	fChain->SetBranchAddress("Kplus_UsedRich1Gas", &Kplus_UsedRich1Gas, &b_Kplus_UsedRich1Gas);
	fChain->SetBranchAddress("Kplus_UsedRich2Gas", &Kplus_UsedRich2Gas, &b_Kplus_UsedRich2Gas);
	fChain->SetBranchAddress("Kplus_RichAboveElThres", &Kplus_RichAboveElThres, &b_Kplus_RichAboveElThres);
	fChain->SetBranchAddress("Kplus_RichAboveMuThres", &Kplus_RichAboveMuThres, &b_Kplus_RichAboveMuThres);
	fChain->SetBranchAddress("Kplus_RichAbovePiThres", &Kplus_RichAbovePiThres, &b_Kplus_RichAbovePiThres);
	fChain->SetBranchAddress("Kplus_RichAboveKaThres", &Kplus_RichAboveKaThres, &b_Kplus_RichAboveKaThres);
	fChain->SetBranchAddress("Kplus_RichAbovePrThres", &Kplus_RichAbovePrThres, &b_Kplus_RichAbovePrThres);
	fChain->SetBranchAddress("Kplus_hasCalo", &Kplus_hasCalo, &b_Kplus_hasCalo);
	fChain->SetBranchAddress("Kplus_PP_CombDLLe", &Kplus_PP_CombDLLe, &b_Kplus_PP_CombDLLe);
	fChain->SetBranchAddress("Kplus_PP_CombDLLmu", &Kplus_PP_CombDLLmu, &b_Kplus_PP_CombDLLmu);
	fChain->SetBranchAddress("Kplus_PP_CombDLLpi", &Kplus_PP_CombDLLpi, &b_Kplus_PP_CombDLLpi);
	fChain->SetBranchAddress("Kplus_PP_CombDLLk", &Kplus_PP_CombDLLk, &b_Kplus_PP_CombDLLk);
	fChain->SetBranchAddress("Kplus_PP_CombDLLp", &Kplus_PP_CombDLLp, &b_Kplus_PP_CombDLLp);
	fChain->SetBranchAddress("Kplus_PP_CombDLLd", &Kplus_PP_CombDLLd, &b_Kplus_PP_CombDLLd);
	fChain->SetBranchAddress("Kplus_PP_ProbNNe", &Kplus_PP_ProbNNe, &b_Kplus_PP_ProbNNe);
	fChain->SetBranchAddress("Kplus_PP_ProbNNmu", &Kplus_PP_ProbNNmu, &b_Kplus_PP_ProbNNmu);
	fChain->SetBranchAddress("Kplus_PP_ProbNNpi", &Kplus_PP_ProbNNpi, &b_Kplus_PP_ProbNNpi);
	fChain->SetBranchAddress("Kplus_PP_ProbNNk", &Kplus_PP_ProbNNk, &b_Kplus_PP_ProbNNk);
	fChain->SetBranchAddress("Kplus_PP_ProbNNp", &Kplus_PP_ProbNNp, &b_Kplus_PP_ProbNNp);
	fChain->SetBranchAddress("Kplus_PP_ProbNNghost", &Kplus_PP_ProbNNghost, &b_Kplus_PP_ProbNNghost);
	fChain->SetBranchAddress("Kplus_PP_ProbNNd", &Kplus_PP_ProbNNd, &b_Kplus_PP_ProbNNd);
	fChain->SetBranchAddress("Kplus_L0Global_Dec", &Kplus_L0Global_Dec, &b_Kplus_L0Global_Dec);
	fChain->SetBranchAddress("Kplus_L0Global_TIS", &Kplus_L0Global_TIS, &b_Kplus_L0Global_TIS);
	fChain->SetBranchAddress("Kplus_L0Global_TOS", &Kplus_L0Global_TOS, &b_Kplus_L0Global_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1Global_Dec", &Kplus_Hlt1Global_Dec, &b_Kplus_Hlt1Global_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1Global_TIS", &Kplus_Hlt1Global_TIS, &b_Kplus_Hlt1Global_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1Global_TOS", &Kplus_Hlt1Global_TOS, &b_Kplus_Hlt1Global_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1Phys_Dec", &Kplus_Hlt1Phys_Dec, &b_Kplus_Hlt1Phys_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1Phys_TIS", &Kplus_Hlt1Phys_TIS, &b_Kplus_Hlt1Phys_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1Phys_TOS", &Kplus_Hlt1Phys_TOS, &b_Kplus_Hlt1Phys_TOS);
	fChain->SetBranchAddress("Kplus_Hlt2Global_Dec", &Kplus_Hlt2Global_Dec, &b_Kplus_Hlt2Global_Dec);
	fChain->SetBranchAddress("Kplus_Hlt2Global_TIS", &Kplus_Hlt2Global_TIS, &b_Kplus_Hlt2Global_TIS);
	fChain->SetBranchAddress("Kplus_Hlt2Global_TOS", &Kplus_Hlt2Global_TOS, &b_Kplus_Hlt2Global_TOS);
	fChain->SetBranchAddress("Kplus_Hlt2Phys_Dec", &Kplus_Hlt2Phys_Dec, &b_Kplus_Hlt2Phys_Dec);
	fChain->SetBranchAddress("Kplus_Hlt2Phys_TIS", &Kplus_Hlt2Phys_TIS, &b_Kplus_Hlt2Phys_TIS);
	fChain->SetBranchAddress("Kplus_Hlt2Phys_TOS", &Kplus_Hlt2Phys_TOS, &b_Kplus_Hlt2Phys_TOS);
	fChain->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_Dec", &Kplus_L0DiHadron_lowMultDecision_Dec, &b_Kplus_L0DiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_TIS", &Kplus_L0DiHadron_lowMultDecision_TIS, &b_Kplus_L0DiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("Kplus_L0DiHadron,lowMultDecision_TOS", &Kplus_L0DiHadron_lowMultDecision_TOS, &b_Kplus_L0DiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_Dec", &Kplus_L0HRCDiHadron_lowMultDecision_Dec, &b_Kplus_L0HRCDiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_TIS", &Kplus_L0HRCDiHadron_lowMultDecision_TIS, &b_Kplus_L0HRCDiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("Kplus_L0HRCDiHadron,lowMultDecision_TOS", &Kplus_L0HRCDiHadron_lowMultDecision_TOS, &b_Kplus_L0HRCDiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_Kplus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_Dec", &Kplus_Hlt1LowMultHerschelDecision_Dec, &b_Kplus_Hlt1LowMultHerschelDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_TIS", &Kplus_Hlt1LowMultHerschelDecision_TIS, &b_Kplus_Hlt1LowMultHerschelDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultHerschelDecision_TOS", &Kplus_Hlt1LowMultHerschelDecision_TOS, &b_Kplus_Hlt1LowMultHerschelDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_Kplus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_Kplus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_Dec", &Kplus_Hlt2LowMultLMR2HHDecision_Dec, &b_Kplus_Hlt2LowMultLMR2HHDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_TIS", &Kplus_Hlt2LowMultLMR2HHDecision_TIS, &b_Kplus_Hlt2LowMultLMR2HHDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHDecision_TOS", &Kplus_Hlt2LowMultLMR2HHDecision_TOS, &b_Kplus_Hlt2LowMultLMR2HHDecision_TOS);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_Dec", &Kplus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_Dec);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_TIS", &Kplus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_TIS);
	fChain->SetBranchAddress("Kplus_Hlt2LowMultLMR2HHWSDecision_TOS", &Kplus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_Kplus_Hlt2LowMultLMR2HHWSDecision_TOS);
	fChain->SetBranchAddress("Kplus_TRACK_Type", &Kplus_TRACK_Type, &b_Kplus_TRACK_Type);
	fChain->SetBranchAddress("Kplus_TRACK_Key", &Kplus_TRACK_Key, &b_Kplus_TRACK_Key);
	fChain->SetBranchAddress("Kplus_TRACK_CHI2NDOF", &Kplus_TRACK_CHI2NDOF, &b_Kplus_TRACK_CHI2NDOF);
	fChain->SetBranchAddress("Kplus_TRACK_PCHI2", &Kplus_TRACK_PCHI2, &b_Kplus_TRACK_PCHI2);
	fChain->SetBranchAddress("Kplus_TRACK_MatchCHI2", &Kplus_TRACK_MatchCHI2, &b_Kplus_TRACK_MatchCHI2);
	fChain->SetBranchAddress("Kplus_TRACK_GhostProb", &Kplus_TRACK_GhostProb, &b_Kplus_TRACK_GhostProb);
	fChain->SetBranchAddress("Kplus_TRACK_CloneDist", &Kplus_TRACK_CloneDist, &b_Kplus_TRACK_CloneDist);
	fChain->SetBranchAddress("Kplus_TRACK_Likelihood", &Kplus_TRACK_Likelihood, &b_Kplus_TRACK_Likelihood);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNe", &Kminus_MC12TuneV2_ProbNNe, &b_Kminus_MC12TuneV2_ProbNNe);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNmu", &Kminus_MC12TuneV2_ProbNNmu, &b_Kminus_MC12TuneV2_ProbNNmu);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNpi", &Kminus_MC12TuneV2_ProbNNpi, &b_Kminus_MC12TuneV2_ProbNNpi);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNk", &Kminus_MC12TuneV2_ProbNNk, &b_Kminus_MC12TuneV2_ProbNNk);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNp", &Kminus_MC12TuneV2_ProbNNp, &b_Kminus_MC12TuneV2_ProbNNp);
	fChain->SetBranchAddress("Kminus_MC12TuneV2_ProbNNghost", &Kminus_MC12TuneV2_ProbNNghost, &b_Kminus_MC12TuneV2_ProbNNghost);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNe", &Kminus_MC12TuneV3_ProbNNe, &b_Kminus_MC12TuneV3_ProbNNe);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNmu", &Kminus_MC12TuneV3_ProbNNmu, &b_Kminus_MC12TuneV3_ProbNNmu);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNpi", &Kminus_MC12TuneV3_ProbNNpi, &b_Kminus_MC12TuneV3_ProbNNpi);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNk", &Kminus_MC12TuneV3_ProbNNk, &b_Kminus_MC12TuneV3_ProbNNk);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNp", &Kminus_MC12TuneV3_ProbNNp, &b_Kminus_MC12TuneV3_ProbNNp);
	fChain->SetBranchAddress("Kminus_MC12TuneV3_ProbNNghost", &Kminus_MC12TuneV3_ProbNNghost, &b_Kminus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNe", &Kminus_MC12TuneV4_ProbNNe, &b_Kminus_MC12TuneV4_ProbNNe);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNmu", &Kminus_MC12TuneV4_ProbNNmu, &b_Kminus_MC12TuneV4_ProbNNmu);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNpi", &Kminus_MC12TuneV4_ProbNNpi, &b_Kminus_MC12TuneV4_ProbNNpi);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNk", &Kminus_MC12TuneV4_ProbNNk, &b_Kminus_MC12TuneV4_ProbNNk);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNp", &Kminus_MC12TuneV4_ProbNNp, &b_Kminus_MC12TuneV4_ProbNNp);
	fChain->SetBranchAddress("Kminus_MC12TuneV4_ProbNNghost", &Kminus_MC12TuneV4_ProbNNghost, &b_Kminus_MC12TuneV4_ProbNNghost);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNe", &Kminus_MC15TuneV1_ProbNNe, &b_Kminus_MC15TuneV1_ProbNNe);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNmu", &Kminus_MC15TuneV1_ProbNNmu, &b_Kminus_MC15TuneV1_ProbNNmu);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNpi", &Kminus_MC15TuneV1_ProbNNpi, &b_Kminus_MC15TuneV1_ProbNNpi);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNk", &Kminus_MC15TuneV1_ProbNNk, &b_Kminus_MC15TuneV1_ProbNNk);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNp", &Kminus_MC15TuneV1_ProbNNp, &b_Kminus_MC15TuneV1_ProbNNp);
	fChain->SetBranchAddress("Kminus_MC15TuneV1_ProbNNghost", &Kminus_MC15TuneV1_ProbNNghost, &b_Kminus_MC15TuneV1_ProbNNghost);
	fChain->SetBranchAddress("Kminus_OWNPV_X", &Kminus_OWNPV_X, &b_Kminus_OWNPV_X);
	fChain->SetBranchAddress("Kminus_OWNPV_Y", &Kminus_OWNPV_Y, &b_Kminus_OWNPV_Y);
	fChain->SetBranchAddress("Kminus_OWNPV_Z", &Kminus_OWNPV_Z, &b_Kminus_OWNPV_Z);
	fChain->SetBranchAddress("Kminus_OWNPV_XERR", &Kminus_OWNPV_XERR, &b_Kminus_OWNPV_XERR);
	fChain->SetBranchAddress("Kminus_OWNPV_YERR", &Kminus_OWNPV_YERR, &b_Kminus_OWNPV_YERR);
	fChain->SetBranchAddress("Kminus_OWNPV_ZERR", &Kminus_OWNPV_ZERR, &b_Kminus_OWNPV_ZERR);
	fChain->SetBranchAddress("Kminus_OWNPV_CHI2", &Kminus_OWNPV_CHI2, &b_Kminus_OWNPV_CHI2);
	fChain->SetBranchAddress("Kminus_OWNPV_NDOF", &Kminus_OWNPV_NDOF, &b_Kminus_OWNPV_NDOF);
	fChain->SetBranchAddress("Kminus_OWNPV_COV_", Kminus_OWNPV_COV_, &b_Kminus_OWNPV_COV_);
	fChain->SetBranchAddress("Kminus_IP_OWNPV", &Kminus_IP_OWNPV, &b_Kminus_IP_OWNPV);
	fChain->SetBranchAddress("Kminus_IPCHI2_OWNPV", &Kminus_IPCHI2_OWNPV, &b_Kminus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("Kminus_ORIVX_X", &Kminus_ORIVX_X, &b_Kminus_ORIVX_X);
	fChain->SetBranchAddress("Kminus_ORIVX_Y", &Kminus_ORIVX_Y, &b_Kminus_ORIVX_Y);
	fChain->SetBranchAddress("Kminus_ORIVX_Z", &Kminus_ORIVX_Z, &b_Kminus_ORIVX_Z);
	fChain->SetBranchAddress("Kminus_ORIVX_XERR", &Kminus_ORIVX_XERR, &b_Kminus_ORIVX_XERR);
	fChain->SetBranchAddress("Kminus_ORIVX_YERR", &Kminus_ORIVX_YERR, &b_Kminus_ORIVX_YERR);
	fChain->SetBranchAddress("Kminus_ORIVX_ZERR", &Kminus_ORIVX_ZERR, &b_Kminus_ORIVX_ZERR);
	fChain->SetBranchAddress("Kminus_ORIVX_CHI2", &Kminus_ORIVX_CHI2, &b_Kminus_ORIVX_CHI2);
	fChain->SetBranchAddress("Kminus_ORIVX_NDOF", &Kminus_ORIVX_NDOF, &b_Kminus_ORIVX_NDOF);
	fChain->SetBranchAddress("Kminus_ORIVX_COV_", Kminus_ORIVX_COV_, &b_Kminus_ORIVX_COV_);
	fChain->SetBranchAddress("Kminus_P", &Kminus_P, &b_Kminus_P);
	fChain->SetBranchAddress("Kminus_PT", &Kminus_PT, &b_Kminus_PT);
	fChain->SetBranchAddress("Kminus_PE", &Kminus_PE, &b_Kminus_PE);
	fChain->SetBranchAddress("Kminus_PX", &Kminus_PX, &b_Kminus_PX);
	fChain->SetBranchAddress("Kminus_PY", &Kminus_PY, &b_Kminus_PY);
	fChain->SetBranchAddress("Kminus_PZ", &Kminus_PZ, &b_Kminus_PZ);
	fChain->SetBranchAddress("Kminus_M", &Kminus_M, &b_Kminus_M);
	fChain->SetBranchAddress("Kminus_ID", &Kminus_ID, &b_Kminus_ID);
	fChain->SetBranchAddress("Kminus_PIDe", &Kminus_PIDe, &b_Kminus_PIDe);
	fChain->SetBranchAddress("Kminus_PIDmu", &Kminus_PIDmu, &b_Kminus_PIDmu);
	fChain->SetBranchAddress("Kminus_PIDK", &Kminus_PIDK, &b_Kminus_PIDK);
	fChain->SetBranchAddress("Kminus_PIDp", &Kminus_PIDp, &b_Kminus_PIDp);
	fChain->SetBranchAddress("Kminus_ProbNNe", &Kminus_ProbNNe, &b_Kminus_ProbNNe);
	fChain->SetBranchAddress("Kminus_ProbNNk", &Kminus_ProbNNk, &b_Kminus_ProbNNk);
	fChain->SetBranchAddress("Kminus_ProbNNp", &Kminus_ProbNNp, &b_Kminus_ProbNNp);
	fChain->SetBranchAddress("Kminus_ProbNNpi", &Kminus_ProbNNpi, &b_Kminus_ProbNNpi);
	fChain->SetBranchAddress("Kminus_ProbNNmu", &Kminus_ProbNNmu, &b_Kminus_ProbNNmu);
	fChain->SetBranchAddress("Kminus_ProbNNghost", &Kminus_ProbNNghost, &b_Kminus_ProbNNghost);
	fChain->SetBranchAddress("Kminus_hasMuon", &Kminus_hasMuon, &b_Kminus_hasMuon);
	fChain->SetBranchAddress("Kminus_isMuon", &Kminus_isMuon, &b_Kminus_isMuon);
	fChain->SetBranchAddress("Kminus_hasRich", &Kminus_hasRich, &b_Kminus_hasRich);
	fChain->SetBranchAddress("Kminus_UsedRichAerogel", &Kminus_UsedRichAerogel, &b_Kminus_UsedRichAerogel);
	fChain->SetBranchAddress("Kminus_UsedRich1Gas", &Kminus_UsedRich1Gas, &b_Kminus_UsedRich1Gas);
	fChain->SetBranchAddress("Kminus_UsedRich2Gas", &Kminus_UsedRich2Gas, &b_Kminus_UsedRich2Gas);
	fChain->SetBranchAddress("Kminus_RichAboveElThres", &Kminus_RichAboveElThres, &b_Kminus_RichAboveElThres);
	fChain->SetBranchAddress("Kminus_RichAboveMuThres", &Kminus_RichAboveMuThres, &b_Kminus_RichAboveMuThres);
	fChain->SetBranchAddress("Kminus_RichAbovePiThres", &Kminus_RichAbovePiThres, &b_Kminus_RichAbovePiThres);
	fChain->SetBranchAddress("Kminus_RichAboveKaThres", &Kminus_RichAboveKaThres, &b_Kminus_RichAboveKaThres);
	fChain->SetBranchAddress("Kminus_RichAbovePrThres", &Kminus_RichAbovePrThres, &b_Kminus_RichAbovePrThres);
	fChain->SetBranchAddress("Kminus_hasCalo", &Kminus_hasCalo, &b_Kminus_hasCalo);
	fChain->SetBranchAddress("Kminus_PP_CombDLLe", &Kminus_PP_CombDLLe, &b_Kminus_PP_CombDLLe);
	fChain->SetBranchAddress("Kminus_PP_CombDLLmu", &Kminus_PP_CombDLLmu, &b_Kminus_PP_CombDLLmu);
	fChain->SetBranchAddress("Kminus_PP_CombDLLpi", &Kminus_PP_CombDLLpi, &b_Kminus_PP_CombDLLpi);
	fChain->SetBranchAddress("Kminus_PP_CombDLLk", &Kminus_PP_CombDLLk, &b_Kminus_PP_CombDLLk);
	fChain->SetBranchAddress("Kminus_PP_CombDLLp", &Kminus_PP_CombDLLp, &b_Kminus_PP_CombDLLp);
	fChain->SetBranchAddress("Kminus_PP_CombDLLd", &Kminus_PP_CombDLLd, &b_Kminus_PP_CombDLLd);
	fChain->SetBranchAddress("Kminus_PP_ProbNNe", &Kminus_PP_ProbNNe, &b_Kminus_PP_ProbNNe);
	fChain->SetBranchAddress("Kminus_PP_ProbNNmu", &Kminus_PP_ProbNNmu, &b_Kminus_PP_ProbNNmu);
	fChain->SetBranchAddress("Kminus_PP_ProbNNpi", &Kminus_PP_ProbNNpi, &b_Kminus_PP_ProbNNpi);
	fChain->SetBranchAddress("Kminus_PP_ProbNNk", &Kminus_PP_ProbNNk, &b_Kminus_PP_ProbNNk);
	fChain->SetBranchAddress("Kminus_PP_ProbNNp", &Kminus_PP_ProbNNp, &b_Kminus_PP_ProbNNp);
	fChain->SetBranchAddress("Kminus_PP_ProbNNghost", &Kminus_PP_ProbNNghost, &b_Kminus_PP_ProbNNghost);
	fChain->SetBranchAddress("Kminus_PP_ProbNNd", &Kminus_PP_ProbNNd, &b_Kminus_PP_ProbNNd);
	fChain->SetBranchAddress("Kminus_L0Global_Dec", &Kminus_L0Global_Dec, &b_Kminus_L0Global_Dec);
	fChain->SetBranchAddress("Kminus_L0Global_TIS", &Kminus_L0Global_TIS, &b_Kminus_L0Global_TIS);
	fChain->SetBranchAddress("Kminus_L0Global_TOS", &Kminus_L0Global_TOS, &b_Kminus_L0Global_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1Global_Dec", &Kminus_Hlt1Global_Dec, &b_Kminus_Hlt1Global_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1Global_TIS", &Kminus_Hlt1Global_TIS, &b_Kminus_Hlt1Global_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1Global_TOS", &Kminus_Hlt1Global_TOS, &b_Kminus_Hlt1Global_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1Phys_Dec", &Kminus_Hlt1Phys_Dec, &b_Kminus_Hlt1Phys_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1Phys_TIS", &Kminus_Hlt1Phys_TIS, &b_Kminus_Hlt1Phys_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1Phys_TOS", &Kminus_Hlt1Phys_TOS, &b_Kminus_Hlt1Phys_TOS);
	fChain->SetBranchAddress("Kminus_Hlt2Global_Dec", &Kminus_Hlt2Global_Dec, &b_Kminus_Hlt2Global_Dec);
	fChain->SetBranchAddress("Kminus_Hlt2Global_TIS", &Kminus_Hlt2Global_TIS, &b_Kminus_Hlt2Global_TIS);
	fChain->SetBranchAddress("Kminus_Hlt2Global_TOS", &Kminus_Hlt2Global_TOS, &b_Kminus_Hlt2Global_TOS);
	fChain->SetBranchAddress("Kminus_Hlt2Phys_Dec", &Kminus_Hlt2Phys_Dec, &b_Kminus_Hlt2Phys_Dec);
	fChain->SetBranchAddress("Kminus_Hlt2Phys_TIS", &Kminus_Hlt2Phys_TIS, &b_Kminus_Hlt2Phys_TIS);
	fChain->SetBranchAddress("Kminus_Hlt2Phys_TOS", &Kminus_Hlt2Phys_TOS, &b_Kminus_Hlt2Phys_TOS);
	fChain->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_Dec", &Kminus_L0DiHadron_lowMultDecision_Dec, &b_Kminus_L0DiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_TIS", &Kminus_L0DiHadron_lowMultDecision_TIS, &b_Kminus_L0DiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("Kminus_L0DiHadron,lowMultDecision_TOS", &Kminus_L0DiHadron_lowMultDecision_TOS, &b_Kminus_L0DiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_Dec", &Kminus_L0HRCDiHadron_lowMultDecision_Dec, &b_Kminus_L0HRCDiHadron_lowMultDecision_Dec);
	fChain->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_TIS", &Kminus_L0HRCDiHadron_lowMultDecision_TIS, &b_Kminus_L0HRCDiHadron_lowMultDecision_TIS);
	fChain->SetBranchAddress("Kminus_L0HRCDiHadron,lowMultDecision_TOS", &Kminus_L0HRCDiHadron_lowMultDecision_TOS, &b_Kminus_L0HRCDiHadron_lowMultDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS", &Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS, &b_Kminus_Hlt1LowMultPassThroughDecisionHlt1LowMultDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_Dec", &Kminus_Hlt1LowMultHerschelDecision_Dec, &b_Kminus_Hlt1LowMultHerschelDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_TIS", &Kminus_Hlt1LowMultHerschelDecision_TIS, &b_Kminus_Hlt1LowMultHerschelDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultHerschelDecision_TOS", &Kminus_Hlt1LowMultHerschelDecision_TOS, &b_Kminus_Hlt1LowMultHerschelDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS", &Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS, &b_Kminus_Hlt1LowMultVeloCut_HadronsDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS", &Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS, &b_Kminus_Hlt1LowMultVeloAndHerschel_HadronsDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_Dec", &Kminus_Hlt2LowMultLMR2HHDecision_Dec, &b_Kminus_Hlt2LowMultLMR2HHDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_TIS", &Kminus_Hlt2LowMultLMR2HHDecision_TIS, &b_Kminus_Hlt2LowMultLMR2HHDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHDecision_TOS", &Kminus_Hlt2LowMultLMR2HHDecision_TOS, &b_Kminus_Hlt2LowMultLMR2HHDecision_TOS);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_Dec", &Kminus_Hlt2LowMultLMR2HHWSDecision_Dec, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_Dec);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_TIS", &Kminus_Hlt2LowMultLMR2HHWSDecision_TIS, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_TIS);
	fChain->SetBranchAddress("Kminus_Hlt2LowMultLMR2HHWSDecision_TOS", &Kminus_Hlt2LowMultLMR2HHWSDecision_TOS, &b_Kminus_Hlt2LowMultLMR2HHWSDecision_TOS);
	fChain->SetBranchAddress("Kminus_TRACK_Type", &Kminus_TRACK_Type, &b_Kminus_TRACK_Type);
	fChain->SetBranchAddress("Kminus_TRACK_Key", &Kminus_TRACK_Key, &b_Kminus_TRACK_Key);
	fChain->SetBranchAddress("Kminus_TRACK_CHI2NDOF", &Kminus_TRACK_CHI2NDOF, &b_Kminus_TRACK_CHI2NDOF);
	fChain->SetBranchAddress("Kminus_TRACK_PCHI2", &Kminus_TRACK_PCHI2, &b_Kminus_TRACK_PCHI2);
	fChain->SetBranchAddress("Kminus_TRACK_MatchCHI2", &Kminus_TRACK_MatchCHI2, &b_Kminus_TRACK_MatchCHI2);
	fChain->SetBranchAddress("Kminus_TRACK_GhostProb", &Kminus_TRACK_GhostProb, &b_Kminus_TRACK_GhostProb);
	fChain->SetBranchAddress("Kminus_TRACK_CloneDist", &Kminus_TRACK_CloneDist, &b_Kminus_TRACK_CloneDist);
	fChain->SetBranchAddress("Kminus_TRACK_Likelihood", &Kminus_TRACK_Likelihood, &b_Kminus_TRACK_Likelihood);
	fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
	fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
	fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
	fChain->SetBranchAddress("LoKi_nCharged", &LoKi_nCharged, &b_LoKi_nCharged);
	fChain->SetBranchAddress("LoKi_nITClusters", &LoKi_nITClusters, &b_LoKi_nITClusters);
	fChain->SetBranchAddress("LoKi_nNeutrals", &LoKi_nNeutrals, &b_LoKi_nNeutrals);
	fChain->SetBranchAddress("LoKi_nOThits", &LoKi_nOThits, &b_LoKi_nOThits);
	fChain->SetBranchAddress("LoKi_nPVs", &LoKi_nPVs, &b_LoKi_nPVs);
	fChain->SetBranchAddress("LoKi_nSpdMult", &LoKi_nSpdMult, &b_LoKi_nSpdMult);
	fChain->SetBranchAddress("LoKi_nTTClusters", &LoKi_nTTClusters, &b_LoKi_nTTClusters);
	fChain->SetBranchAddress("LoKi_nVeloClusters", &LoKi_nVeloClusters, &b_LoKi_nVeloClusters);
	fChain->SetBranchAddress("LoKi_nVeloLiteClusters", &LoKi_nVeloLiteClusters, &b_LoKi_nVeloLiteClusters);
	fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
	fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
	fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
	fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
	fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
	fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
	fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
	fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
	fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
	fChain->SetBranchAddress("B00", &B00, &b_B00);
	fChain->SetBranchAddress("B01", &B01, &b_B01);
	fChain->SetBranchAddress("B02", &B02, &b_B02);
	fChain->SetBranchAddress("B03", &B03, &b_B03);
	fChain->SetBranchAddress("B10", &B10, &b_B10);
	fChain->SetBranchAddress("B11", &B11, &b_B11);
	fChain->SetBranchAddress("B12", &B12, &b_B12);
	fChain->SetBranchAddress("B13", &B13, &b_B13);
	fChain->SetBranchAddress("B20", &B20, &b_B20);
	fChain->SetBranchAddress("B21", &B21, &b_B21);
	fChain->SetBranchAddress("B22", &B22, &b_B22);
	fChain->SetBranchAddress("B23", &B23, &b_B23);
	fChain->SetBranchAddress("F10", &F10, &b_F10);
	fChain->SetBranchAddress("F11", &F11, &b_F11);
	fChain->SetBranchAddress("F12", &F12, &b_F12);
	fChain->SetBranchAddress("F13", &F13, &b_F13);
	fChain->SetBranchAddress("F20", &F20, &b_F20);
	fChain->SetBranchAddress("F21", &F21, &b_F21);
	fChain->SetBranchAddress("F22", &F22, &b_F22);
	fChain->SetBranchAddress("F23", &F23, &b_F23);
	fChain->SetBranchAddress("log_hrc_fom_v3", &log_hrc_fom_v3, &b_log_hrc_fom_v3);
	fChain->SetBranchAddress("log_hrc_fom_B_v3", &log_hrc_fom_B_v3, &b_log_hrc_fom_B_v3);
	fChain->SetBranchAddress("log_hrc_fom_F_v3", &log_hrc_fom_F_v3, &b_log_hrc_fom_F_v3);
	fChain->SetBranchAddress("nchB", &nchB, &b_nchB);
	fChain->SetBranchAddress("adc_B", adc_B, &b_adc_B);
	fChain->SetBranchAddress("nchF", &nchF, &b_nchF);
	fChain->SetBranchAddress("adc_F", adc_F, &b_adc_F);
	fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
	fChain->SetBranchAddress("PVX", PVX, &b_PVX);
	fChain->SetBranchAddress("PVY", PVY, &b_PVY);
	fChain->SetBranchAddress("PVZ", PVZ, &b_PVZ);
	fChain->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
	fChain->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
	fChain->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
	fChain->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
	fChain->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
	fChain->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
	fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
	fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	fChain->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
	fChain->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
	fChain->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
	fChain->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
	fChain->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
	fChain->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
	fChain->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
	fChain->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
	fChain->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
	fChain->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
	fChain->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
	fChain->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
	fChain->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
	fChain->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
	fChain->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
	fChain->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
	fChain->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
	fChain->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
	fChain->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
	Notify();
}

Bool_t skim::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void skim::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t skim::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef skim_cxx
