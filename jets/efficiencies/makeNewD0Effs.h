//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 24 17:44:11 2018 by ROOT version 6.06/00
// from TTree T/
// found on file: /tmp/dcraik/for_yandex_data_new_245X.root
//////////////////////////////////////////////////////////

#ifndef makeNewD0Effs_h
#define makeNewD0Effs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

class makeNewD0Effs {
	public :
		enum Binning {
			TwelveByThree,
			TwelveByFive,
			SixteenByFive,
			SixteenByEight,
			TwentyOneByFive,
			ThirtyFourByFive,
			ThirtySevenByFive,
			ThirtySevenByEight
		};

		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		Int_t           Evt;
		Int_t           Dec;
		Int_t           NPV;
		Double_t        JetPx;
		Double_t        JetPy;
		Double_t        JetPz;
		Double_t        JetE;
		Double_t        JetPT;
		Double_t        JetEta;
		Double_t        JetTruePx;
		Double_t        JetTruePy;
		Double_t        JetTruePz;
		Double_t        JetTrueE;
		Double_t        JetTruePT;
		Double_t        JetTrueEta;
		Double_t        JetTrueDR;
		Double_t        JetSigma1;
		Double_t        JetSigma2;
		Double_t        JetQ;
		Double_t        JetMult;
		Double_t        JetNChr;
		Double_t        JetNNeu;
		Double_t        JetPTD;
		Double_t        JetTRUEb;
		Double_t        JetTRUEc;
		Double_t        JetTRUED0;
		Double_t        JetTRUEDP;
		Double_t        JetTRUEDS;
		Double_t        JetTRUEJPSI;
		Double_t        JetTRUEDSV;
		Double_t        JetTRUEBSV;
		Double_t        JetTRUESV;
		Double_t        JetTRUEDPX;
		Double_t        JetTRUEDPY;
		Double_t        JetTRUEDPZ;
		Double_t        JetTRUEDE;
		Double_t        JetTRUEDDR;
		Double_t        TagPx;
		Double_t        TagPy;
		Double_t        TagPz;
		Double_t        TagE;
		Double_t        TagPT;
		Double_t        TagEta;
		Double_t        DeltaR;
		Double_t        DeltaPhi;
		std::vector<double>  *SVX;
		std::vector<double>  *SVY;
		std::vector<double>  *SVZ;
		std::vector<double>  *SVPerp;
		std::vector<double>  *SVPx;
		std::vector<double>  *SVPy;
		std::vector<double>  *SVPz;
		std::vector<double>  *SVE;
		std::vector<double>  *SVPT;
		std::vector<double>  *SVETA;
		std::vector<double>  *SVM;
		std::vector<double>  *SVMCor;
		std::vector<double>  *SVMCorErr;
		std::vector<double>  *SVMINPERP;
		std::vector<double>  *SVDRJ;
		std::vector<double>  *SVDRT;
		std::vector<double>  *SVN;
		std::vector<double>  *SVNTRUE;
		std::vector<double>  *SVNJ;
		std::vector<double>  *SVNT;
		std::vector<double>  *SVQ;
		std::vector<double>  *SVSumIPChi2;
		std::vector<double>  *SVTZ;
		std::vector<double>  *SVIPCHI2;
		std::vector<double>  *SVMINIPCHI2;
		std::vector<double>  *SVGhostMax;
		std::vector<double>  *SVISD0;
		std::vector<double>  *SVISDP;
		std::vector<double>  *SVD0M;
		std::vector<double>  *SVDPM;
		std::vector<double>  *SVTRUEIDX;
		std::vector<double>  *SVTRUETRK0IDX;
		std::vector<double>  *SVTRUETRK1IDX;
		std::vector<double>  *SVTRUETRK2IDX;
		std::vector<double>  *SVTRUETRK3IDX;
		std::vector<double>  *SVTRK0P;
		std::vector<double>  *SVTRK1P;
		std::vector<double>  *SVTRK2P;
		std::vector<double>  *SVTRK3P;
		std::vector<double>  *SVTRK0PT;
		std::vector<double>  *SVTRK1PT;
		std::vector<double>  *SVTRK2PT;
		std::vector<double>  *SVTRK3PT;
		std::vector<double>  *SVTRK0PX;
		std::vector<double>  *SVTRK1PX;
		std::vector<double>  *SVTRK2PX;
		std::vector<double>  *SVTRK3PX;
		std::vector<double>  *SVTRK0PY;
		std::vector<double>  *SVTRK1PY;
		std::vector<double>  *SVTRK2PY;
		std::vector<double>  *SVTRK3PY;
		std::vector<double>  *SVTRK0PZ;
		std::vector<double>  *SVTRK1PZ;
		std::vector<double>  *SVTRK2PZ;
		std::vector<double>  *SVTRK3PZ;
		std::vector<double>  *SVTRK0PNNPI;
		std::vector<double>  *SVTRK1PNNPI;
		std::vector<double>  *SVTRK2PNNPI;
		std::vector<double>  *SVTRK3PNNPI;
		std::vector<double>  *SVTRK0PNNK;
		std::vector<double>  *SVTRK1PNNK;
		std::vector<double>  *SVTRK2PNNK;
		std::vector<double>  *SVTRK3PNNK;
		std::vector<double>  *TSVX;
		std::vector<double>  *TSVY;
		std::vector<double>  *TSVZ;
		std::vector<double>  *TSVPerp;
		std::vector<double>  *TSVPx;
		std::vector<double>  *TSVPy;
		std::vector<double>  *TSVPz;
		std::vector<double>  *TSVE;
		std::vector<double>  *TSVPT;
		std::vector<double>  *TSVETA;
		std::vector<double>  *TSVM;
		std::vector<double>  *TSVMCor;
		std::vector<double>  *TSVMCorErr;
		std::vector<double>  *TSVMINPERP;
		std::vector<double>  *TSVDRJ;
		std::vector<double>  *TSVDRT;
		std::vector<double>  *TSVN;
		std::vector<double>  *TSVNTRUE;
		std::vector<double>  *TSVNJ;
		std::vector<double>  *TSVNT;
		std::vector<double>  *TSVQ;
		std::vector<double>  *TSVSumIPChi2;
		std::vector<double>  *TSVTZ;
		std::vector<double>  *TSVIPCHI2;
		std::vector<double>  *TSVMINIPCHI2;
		std::vector<double>  *TSVGhostMax;
		std::vector<double>  *TSVISD0;
		std::vector<double>  *TSVISDP;
		std::vector<double>  *TSVD0M;
		std::vector<double>  *TSVDPM;
		std::vector<double>  *TSVTRUEIDX;
		Double_t        NSV;
		Double_t        NTSV;
		Double_t        NTRK;
		Double_t        NNEU;
		Double_t        PVX;
		Double_t        PVY;
		Double_t        PVZ;
		Double_t        NDispl6;
		Double_t        NDispl9;
		Double_t        NDispl16;
		Double_t        MuPT;
		Double_t        MuIPChi2;
		Double_t        MuDR;
		Double_t        MuPNN;
		Double_t        NMu;
		Double_t        HardPT;
		Double_t        HardIPChi2;
		Double_t        HardDR;
		std::vector<double>  *TRUEDID;
		std::vector<double>  *TRUEDPX;
		std::vector<double>  *TRUEDPY;
		std::vector<double>  *TRUEDPZ;
		std::vector<double>  *TRUEDE;
		std::vector<double>  *TRUEDFROMB;
		std::vector<double>  *TRUEDSEL;
		std::vector<double>  *TRUEDTRK0IDX;
		std::vector<double>  *TRUEDTRK0ID;
		std::vector<double>  *TRUEDTRK0P;
		std::vector<double>  *TRUEDTRK0PT;
		std::vector<double>  *TRUEDTRK0INACC;
		std::vector<double>  *TRUEDTRK0RECO;
		std::vector<double>  *TRUEDTRK0PNNK;
		std::vector<double>  *TRUEDTRK0PNNPI;
		std::vector<double>  *TRUEDTRK1IDX;
		std::vector<double>  *TRUEDTRK1ID;
		std::vector<double>  *TRUEDTRK1P;
		std::vector<double>  *TRUEDTRK1PT;
		std::vector<double>  *TRUEDTRK1INACC;
		std::vector<double>  *TRUEDTRK1RECO;
		std::vector<double>  *TRUEDTRK1PNNK;
		std::vector<double>  *TRUEDTRK1PNNPI;
		std::vector<double>  *TRUEDTRK2IDX;
		std::vector<double>  *TRUEDTRK2ID;
		std::vector<double>  *TRUEDTRK2P;
		std::vector<double>  *TRUEDTRK2PT;
		std::vector<double>  *TRUEDTRK2INACC;
		std::vector<double>  *TRUEDTRK2RECO;
		std::vector<double>  *TRUEDTRK2PNNK;
		std::vector<double>  *TRUEDTRK2PNNPI;
		std::vector<double>  *TRUEDTRK3IDX;
		std::vector<double>  *TRUEDTRK3ID;
		std::vector<double>  *TRUEDTRK3P;
		std::vector<double>  *TRUEDTRK3PT;
		std::vector<double>  *TRUEDTRK3INACC;
		std::vector<double>  *TRUEDTRK3RECO;
		std::vector<double>  *TRUEDTRK3PNNK;
		std::vector<double>  *TRUEDTRK3PNNPI;
		std::vector<double>  *TRUEDTRK0PX;
		std::vector<double>  *TRUEDTRK1PX;
		std::vector<double>  *TRUEDTRK2PX;
		std::vector<double>  *TRUEDTRK3PX;
		std::vector<double>  *TRUEDTRK0PY;
		std::vector<double>  *TRUEDTRK1PY;
		std::vector<double>  *TRUEDTRK2PY;
		std::vector<double>  *TRUEDTRK3PY;
		std::vector<double>  *TRUEDTRK0PZ;
		std::vector<double>  *TRUEDTRK1PZ;
		std::vector<double>  *TRUEDTRK2PZ;
		std::vector<double>  *TRUEDTRK3PZ;
		std::vector<double>  *D0M;
		std::vector<double>  *D0PX;
		std::vector<double>  *D0PY;
		std::vector<double>  *D0PZ;
		std::vector<double>  *D0E;
		std::vector<double>  *D0X;
		std::vector<double>  *D0Y;
		std::vector<double>  *D0Z;
		std::vector<double>  *D0IP;
		std::vector<double>  *D0IPCHI2;
		std::vector<double>  *D0FD;
		std::vector<double>  *D0FDCHI2;
		std::vector<double>  *D0TAU;
		std::vector<double>  *D0DIRA;
		std::vector<double>  *D0VTXCHI2;
		std::vector<double>  *D0VTXNDOF;
		std::vector<double>  *D0DRJET;
		std::vector<double>  *D0DRTAG;
		std::vector<double>  *D0PIPNNK;
		std::vector<double>  *D0KPNNK;
		std::vector<double>  *D0PIPNNPI;
		std::vector<double>  *D0KPNNPI;
		std::vector<double>  *D0P;
		std::vector<double>  *D0PT;
		std::vector<double>  *D0ETA;
		std::vector<double>  *D0PIP;
		std::vector<double>  *D0KP;
		std::vector<double>  *D0PIPT;
		std::vector<double>  *D0KPT;
		std::vector<double>  *D0PIETA;
		std::vector<double>  *D0KETA;
		std::vector<double>  *D0PIPX;
		std::vector<double>  *D0KPX;
		std::vector<double>  *D0PIPY;
		std::vector<double>  *D0KPY;
		std::vector<double>  *D0PIPZ;
		std::vector<double>  *D0KPZ;
		std::vector<double>  *D0PIIPCHI2;
		std::vector<double>  *D0KIPCHI2;
		std::vector<double>  *D0TRK0;
		std::vector<double>  *D0TRK1;
		std::vector<double>  *D0TRUETRK0;
		std::vector<double>  *D0TRUETRK1;
		std::vector<double>  *D0NJ;
		std::vector<double>  *D0MAXDR;
		std::vector<double>  *D0TRUE;
		std::vector<double>  *D0TRUEDR;
		std::vector<double>  *D0TRUEIDX;
		std::vector<double>  *D0FROMB;
		std::vector<double>  *DM;
		std::vector<double>  *DPX;
		std::vector<double>  *DPY;
		std::vector<double>  *DPZ;
		std::vector<double>  *DE;
		std::vector<double>  *DX;
		std::vector<double>  *DY;
		std::vector<double>  *DZ;
		std::vector<double>  *DIP;
		std::vector<double>  *DIPCHI2;
		std::vector<double>  *DFD;
		std::vector<double>  *DFDCHI2;
		std::vector<double>  *DTAU;
		std::vector<double>  *DDIRA;
		std::vector<double>  *DVTXCHI2;
		std::vector<double>  *DVTXNDOF;
		std::vector<double>  *DDRJET;
		std::vector<double>  *DDRTAG;
		std::vector<double>  *DPI1PNN;
		std::vector<double>  *DPI2PNN;
		std::vector<double>  *DKPNN;
		std::vector<double>  *DP;
		std::vector<double>  *DPT;
		std::vector<double>  *DPI1PT;
		std::vector<double>  *DPI2PT;
		std::vector<double>  *DKPT;
		std::vector<double>  *DTRK0;
		std::vector<double>  *DTRK1;
		std::vector<double>  *DTRK2;
		std::vector<double>  *DNJ;
		std::vector<double>  *DMAXDR;
		std::vector<double>  *DTRUE;
		std::vector<double>  *DTRUEDR;
		std::vector<double>  *DTRUEIDX;
		std::vector<double>  *DFROMB;
		std::vector<double>  *DSM;
		std::vector<double>  *DSPX;
		std::vector<double>  *DSPY;
		std::vector<double>  *DSPZ;
		std::vector<double>  *DSE;
		std::vector<double>  *DSX;
		std::vector<double>  *DSY;
		std::vector<double>  *DSZ;
		std::vector<double>  *DSIP;
		std::vector<double>  *DSIPCHI2;
		std::vector<double>  *DSFD;
		std::vector<double>  *DSFDCHI2;
		std::vector<double>  *DSTAU;
		std::vector<double>  *DSDIRA;
		std::vector<double>  *DSVTXCHI2;
		std::vector<double>  *DSVTXNDOF;
		std::vector<double>  *DSDRJET;
		std::vector<double>  *DSDRTAG;
		std::vector<double>  *DSPIPNN;
		std::vector<double>  *DSK1PNN;
		std::vector<double>  *DSK2PNN;
		std::vector<double>  *DSPHIM;
		std::vector<double>  *DSP;
		std::vector<double>  *DSPT;
		std::vector<double>  *DSPIPT;
		std::vector<double>  *DSK1PT;
		std::vector<double>  *DSK2PT;
		std::vector<double>  *DSPHIPT;
		std::vector<double>  *DSTRK0;
		std::vector<double>  *DSTRK1;
		std::vector<double>  *DSTRK2;
		std::vector<double>  *DSNJ;
		std::vector<double>  *DSMAXDR;
		std::vector<double>  *DSTRUE;
		std::vector<double>  *DSTRUEDR;
		std::vector<double>  *DSTRUEIDX;
		std::vector<double>  *DSFROMB;
		std::vector<double>  *LCM;
		std::vector<double>  *LCPX;
		std::vector<double>  *LCPY;
		std::vector<double>  *LCPZ;
		std::vector<double>  *LCE;
		std::vector<double>  *LCX;
		std::vector<double>  *LCY;
		std::vector<double>  *LCZ;
		std::vector<double>  *LCIP;
		std::vector<double>  *LCIPCHI2;
		std::vector<double>  *LCFD;
		std::vector<double>  *LCFDCHI2;
		std::vector<double>  *LCTAU;
		std::vector<double>  *LCDIRA;
		std::vector<double>  *LCVTXCHI2;
		std::vector<double>  *LCVTXNDOF;
		std::vector<double>  *LCDRJET;
		std::vector<double>  *LCDRTAG;
		std::vector<double>  *LCPPNN;
		std::vector<double>  *LCKPNN;
		std::vector<double>  *LCPIPNN;
		std::vector<double>  *LCP;
		std::vector<double>  *LCPT;
		std::vector<double>  *LCPPT;
		std::vector<double>  *LCKPT;
		std::vector<double>  *LCPIPT;
		std::vector<double>  *LCTRK0;
		std::vector<double>  *LCTRK1;
		std::vector<double>  *LCTRK2;
		std::vector<double>  *LCNJ;
		std::vector<double>  *LCMAXDR;
		std::vector<double>  *K3PIM;
		std::vector<double>  *K3PIPX;
		std::vector<double>  *K3PIPY;
		std::vector<double>  *K3PIPZ;
		std::vector<double>  *K3PIE;
		std::vector<double>  *K3PIX;
		std::vector<double>  *K3PIY;
		std::vector<double>  *K3PIZ;
		std::vector<double>  *K3PIIP;
		std::vector<double>  *K3PIIPCHI2;
		std::vector<double>  *K3PIFD;
		std::vector<double>  *K3PIFDCHI2;
		std::vector<double>  *K3PITAU;
		std::vector<double>  *K3PIDIRA;
		std::vector<double>  *K3PIVTXCHI2;
		std::vector<double>  *K3PIVTXNDOF;
		std::vector<double>  *K3PIDRJET;
		std::vector<double>  *K3PIDRTAG;
		std::vector<double>  *K3PIPI1PNN;
		std::vector<double>  *K3PIPI2PNN;
		std::vector<double>  *K3PIPI3PNN;
		std::vector<double>  *K3PIKPNN;
		std::vector<double>  *K3PIP;
		std::vector<double>  *K3PIPT;
		std::vector<double>  *K3PIPI1PT;
		std::vector<double>  *K3PIPI2PT;
		std::vector<double>  *K3PIPI3PT;
		std::vector<double>  *K3PIKPT;
		std::vector<double>  *K3PITRK0;
		std::vector<double>  *K3PITRK1;
		std::vector<double>  *K3PITRK2;
		std::vector<double>  *K3PITRK3;
		std::vector<double>  *K3PINJ;
		std::vector<double>  *K3PIMAXDR;
		std::vector<double>  *JPSIM;
		std::vector<double>  *JPSIPX;
		std::vector<double>  *JPSIPY;
		std::vector<double>  *JPSIPZ;
		std::vector<double>  *JPSIE;
		std::vector<double>  *JPSIX;
		std::vector<double>  *JPSIY;
		std::vector<double>  *JPSIZ;
		std::vector<double>  *JPSIIP;
		std::vector<double>  *JPSIIPCHI2;
		std::vector<double>  *JPSIFD;
		std::vector<double>  *JPSIFDCHI2;
		std::vector<double>  *JPSITAU;
		std::vector<double>  *JPSIDIRA;
		std::vector<double>  *JPSIVTXCHI2;
		std::vector<double>  *JPSIVTXNDOF;
		std::vector<double>  *JPSIDRJET;
		std::vector<double>  *JPSIDRTAG;
		std::vector<double>  *JPSIP;
		std::vector<double>  *JPSIPT;
		std::vector<double>  *JPSIPIPT;
		std::vector<double>  *JPSIKPT;
		std::vector<double>  *JPSITRK0;
		std::vector<double>  *JPSITRK1;
		std::vector<double>  *JPSINJ;
		std::vector<double>  *JPSIMAXDR;
		std::vector<double>  *JPSITRUE;
		std::vector<double>  *JPSITRUEDR;
		std::vector<double>  *JPSITRUEIDX;
		std::vector<double>  *JPSIFROMB;

		// List of branches
		TBranch        *b_Evt;   //!
		TBranch        *b_Dec;   //!
		TBranch        *b_NPV;   //!
		TBranch        *b_JetPx;   //!
		TBranch        *b_JetPy;   //!
		TBranch        *b_JetPz;   //!
		TBranch        *b_JetE;   //!
		TBranch        *b_JetPT;   //!
		TBranch        *b_JetEta;   //!
		TBranch        *b_JetTruePx;   //!
		TBranch        *b_JetTruePy;   //!
		TBranch        *b_JetTruePz;   //!
		TBranch        *b_JetTrueE;   //!
		TBranch        *b_JetTruePT;   //!
		TBranch        *b_JetTrueEta;   //!
		TBranch        *b_JetTrueDR;   //!
		TBranch        *b_JetSigma1;   //!
		TBranch        *b_JetSigma2;   //!
		TBranch        *b_JetQ;   //!
		TBranch        *b_JetMult;   //!
		TBranch        *b_JetNChr;   //!
		TBranch        *b_JetNNeu;   //!
		TBranch        *b_JetPTD;   //!
		TBranch        *b_JetTRUEb;   //!
		TBranch        *b_JetTRUEc;   //!
		TBranch        *b_JetTRUED0;   //!
		TBranch        *b_JetTRUEDP;   //!
		TBranch        *b_JetTRUEDS;   //!
		TBranch        *b_JetTRUEJPSI;   //!
		TBranch        *b_JetTRUEDSV;   //!
		TBranch        *b_JetTRUEBSV;   //!
		TBranch        *b_JetTRUESV;   //!
		TBranch        *b_JetTRUEDPX;   //!
		TBranch        *b_JetTRUEDPY;   //!
		TBranch        *b_JetTRUEDPZ;   //!
		TBranch        *b_JetTRUEDE;   //!
		TBranch        *b_JetTRUEDDR;   //!
		TBranch        *b_TagPx;   //!
		TBranch        *b_TagPy;   //!
		TBranch        *b_TagPz;   //!
		TBranch        *b_TagE;   //!
		TBranch        *b_TagPT;   //!
		TBranch        *b_TagEta;   //!
		TBranch        *b_DeltaR;   //!
		TBranch        *b_DeltaPhi;   //!
		TBranch        *b_SVX;   //!
		TBranch        *b_SVY;   //!
		TBranch        *b_SVZ;   //!
		TBranch        *b_SVPerp;   //!
		TBranch        *b_SVPx;   //!
		TBranch        *b_SVPy;   //!
		TBranch        *b_SVPz;   //!
		TBranch        *b_SVE;   //!
		TBranch        *b_SVPT;   //!
		TBranch        *b_SVETA;   //!
		TBranch        *b_SVM;   //!
		TBranch        *b_SVMCor;   //!
		TBranch        *b_SVMCorErr;   //!
		TBranch        *b_SVMINPERP;   //!
		TBranch        *b_SVDRJ;   //!
		TBranch        *b_SVDRT;   //!
		TBranch        *b_SVN;   //!
		TBranch        *b_SVNTRUE;   //!
		TBranch        *b_SVNJ;   //!
		TBranch        *b_SVNT;   //!
		TBranch        *b_SVQ;   //!
		TBranch        *b_SVSumIPChi2;   //!
		TBranch        *b_SVTZ;   //!
		TBranch        *b_SVIPCHI2;   //!
		TBranch        *b_SVMINIPCHI2;   //!
		TBranch        *b_SVGhostMax;   //!
		TBranch        *b_SVISD0;   //!
		TBranch        *b_SVISDP;   //!
		TBranch        *b_SVD0M;   //!
		TBranch        *b_SVDPM;   //!
		TBranch        *b_SVTRUEIDX;   //!
		TBranch        *b_SVTRUETRK0IDX;   //!
		TBranch        *b_SVTRUETRK1IDX;   //!
		TBranch        *b_SVTRUETRK2IDX;   //!
		TBranch        *b_SVTRUETRK3IDX;   //!
		TBranch        *b_SVTRK0P;   //!
		TBranch        *b_SVTRK1P;   //!
		TBranch        *b_SVTRK2P;   //!
		TBranch        *b_SVTRK3P;   //!
		TBranch        *b_SVTRK0PT;   //!
		TBranch        *b_SVTRK1PT;   //!
		TBranch        *b_SVTRK2PT;   //!
		TBranch        *b_SVTRK3PT;   //!
		TBranch        *b_SVTRK0PX;   //!
		TBranch        *b_SVTRK1PX;   //!
		TBranch        *b_SVTRK2PX;   //!
		TBranch        *b_SVTRK3PX;   //!
		TBranch        *b_SVTRK0PY;   //!
		TBranch        *b_SVTRK1PY;   //!
		TBranch        *b_SVTRK2PY;   //!
		TBranch        *b_SVTRK3PY;   //!
		TBranch        *b_SVTRK0PZ;   //!
		TBranch        *b_SVTRK1PZ;   //!
		TBranch        *b_SVTRK2PZ;   //!
		TBranch        *b_SVTRK3PZ;   //!
		TBranch        *b_SVTRK0PNNPI;   //!
		TBranch        *b_SVTRK1PNNPI;   //!
		TBranch        *b_SVTRK2PNNPI;   //!
		TBranch        *b_SVTRK3PNNPI;   //!
		TBranch        *b_SVTRK0PNNK;   //!
		TBranch        *b_SVTRK1PNNK;   //!
		TBranch        *b_SVTRK2PNNK;   //!
		TBranch        *b_SVTRK3PNNK;   //!
		TBranch        *b_TSVX;   //!
		TBranch        *b_TSVY;   //!
		TBranch        *b_TSVZ;   //!
		TBranch        *b_TSVPerp;   //!
		TBranch        *b_TSVPx;   //!
		TBranch        *b_TSVPy;   //!
		TBranch        *b_TSVPz;   //!
		TBranch        *b_TSVE;   //!
		TBranch        *b_TSVPT;   //!
		TBranch        *b_TSVETA;   //!
		TBranch        *b_TSVM;   //!
		TBranch        *b_TSVMCor;   //!
		TBranch        *b_TSVMCorErr;   //!
		TBranch        *b_TSVMINPERP;   //!
		TBranch        *b_TSVDRJ;   //!
		TBranch        *b_TSVDRT;   //!
		TBranch        *b_TSVN;   //!
		TBranch        *b_TSVNTRUE;   //!
		TBranch        *b_TSVNJ;   //!
		TBranch        *b_TSVNT;   //!
		TBranch        *b_TSVQ;   //!
		TBranch        *b_TSVSumIPChi2;   //!
		TBranch        *b_TSVTZ;   //!
		TBranch        *b_TSVIPCHI2;   //!
		TBranch        *b_TSVMINIPCHI2;   //!
		TBranch        *b_TSVGhostMax;   //!
		TBranch        *b_TSVISD0;   //!
		TBranch        *b_TSVISDP;   //!
		TBranch        *b_TSVD0M;   //!
		TBranch        *b_TSVDPM;   //!
		TBranch        *b_TSVTRUEIDX;   //!
		TBranch        *b_NSV;   //!
		TBranch        *b_NTSV;   //!
		TBranch        *b_NTRK;   //!
		TBranch        *b_NNEU;   //!
		TBranch        *b_PVX;   //!
		TBranch        *b_PVY;   //!
		TBranch        *b_PVZ;   //!
		TBranch        *b_NDispl6;   //!
		TBranch        *b_NDispl9;   //!
		TBranch        *b_NDispl16;   //!
		TBranch        *b_MuPT;   //!
		TBranch        *b_MuIPChi2;   //!
		TBranch        *b_MuDR;   //!
		TBranch        *b_MuPNN;   //!
		TBranch        *b_NMu;   //!
		TBranch        *b_HardPT;   //!
		TBranch        *b_HardIPChi2;   //!
		TBranch        *b_HardDR;   //!
		TBranch        *b_TRUEDID;   //!
		TBranch        *b_TRUEDPX;   //!
		TBranch        *b_TRUEDPY;   //!
		TBranch        *b_TRUEDPZ;   //!
		TBranch        *b_TRUEDE;   //!
		TBranch        *b_TRUEDFROMB;   //!
		TBranch        *b_TRUEDSEL;   //!
		TBranch        *b_TRUEDTRK0IDX;   //!
		TBranch        *b_TRUEDTRK0ID;   //!
		TBranch        *b_TRUEDTRK0P;   //!
		TBranch        *b_TRUEDTRK0PT;   //!
		TBranch        *b_TRUEDTRK0INACC;   //!
		TBranch        *b_TRUEDTRK0RECO;   //!
		TBranch        *b_TRUEDTRK0PNNK;   //!
		TBranch        *b_TRUEDTRK0PNNPI;   //!
		TBranch        *b_TRUEDTRK1IDX;   //!
		TBranch        *b_TRUEDTRK1ID;   //!
		TBranch        *b_TRUEDTRK1P;   //!
		TBranch        *b_TRUEDTRK1PT;   //!
		TBranch        *b_TRUEDTRK1INACC;   //!
		TBranch        *b_TRUEDTRK1RECO;   //!
		TBranch        *b_TRUEDTRK1PNNK;   //!
		TBranch        *b_TRUEDTRK1PNNPI;   //!
		TBranch        *b_TRUEDTRK2IDX;   //!
		TBranch        *b_TRUEDTRK2ID;   //!
		TBranch        *b_TRUEDTRK2P;   //!
		TBranch        *b_TRUEDTRK2PT;   //!
		TBranch        *b_TRUEDTRK2INACC;   //!
		TBranch        *b_TRUEDTRK2RECO;   //!
		TBranch        *b_TRUEDTRK2PNNK;   //!
		TBranch        *b_TRUEDTRK2PNNPI;   //!
		TBranch        *b_TRUEDTRK3IDX;   //!
		TBranch        *b_TRUEDTRK3ID;   //!
		TBranch        *b_TRUEDTRK3P;   //!
		TBranch        *b_TRUEDTRK3PT;   //!
		TBranch        *b_TRUEDTRK3INACC;   //!
		TBranch        *b_TRUEDTRK3RECO;   //!
		TBranch        *b_TRUEDTRK3PNNK;   //!
		TBranch        *b_TRUEDTRK3PNNPI;   //!
		TBranch        *b_TRUEDTRK0PX;   //!
		TBranch        *b_TRUEDTRK1PX;   //!
		TBranch        *b_TRUEDTRK2PX;   //!
		TBranch        *b_TRUEDTRK3PX;   //!
		TBranch        *b_TRUEDTRK0PY;   //!
		TBranch        *b_TRUEDTRK1PY;   //!
		TBranch        *b_TRUEDTRK2PY;   //!
		TBranch        *b_TRUEDTRK3PY;   //!
		TBranch        *b_TRUEDTRK0PZ;   //!
		TBranch        *b_TRUEDTRK1PZ;   //!
		TBranch        *b_TRUEDTRK2PZ;   //!
		TBranch        *b_TRUEDTRK3PZ;   //!
		TBranch        *b_D0M;   //!
		TBranch        *b_D0PX;   //!
		TBranch        *b_D0PY;   //!
		TBranch        *b_D0PZ;   //!
		TBranch        *b_D0E;   //!
		TBranch        *b_D0X;   //!
		TBranch        *b_D0Y;   //!
		TBranch        *b_D0Z;   //!
		TBranch        *b_D0IP;   //!
		TBranch        *b_D0IPCHI2;   //!
		TBranch        *b_D0FD;   //!
		TBranch        *b_D0FDCHI2;   //!
		TBranch        *b_D0TAU;   //!
		TBranch        *b_D0DIRA;   //!
		TBranch        *b_D0VTXCHI2;   //!
		TBranch        *b_D0VTXNDOF;   //!
		TBranch        *b_D0DRJET;   //!
		TBranch        *b_D0DRTAG;   //!
		TBranch        *b_D0PIPNNK;   //!
		TBranch        *b_D0KPNNK;   //!
		TBranch        *b_D0PIPNNPI;   //!
		TBranch        *b_D0KPNNPI;   //!
		TBranch        *b_D0P;   //!
		TBranch        *b_D0PT;   //!
		TBranch        *b_D0ETA;   //!
		TBranch        *b_D0PIP;   //!
		TBranch        *b_D0KP;   //!
		TBranch        *b_D0PIPT;   //!
		TBranch        *b_D0KPT;   //!
		TBranch        *b_D0PIETA;   //!
		TBranch        *b_D0KETA;   //!
		TBranch        *b_D0PIPX;   //!
		TBranch        *b_D0KPX;   //!
		TBranch        *b_D0PIPY;   //!
		TBranch        *b_D0KPY;   //!
		TBranch        *b_D0PIPZ;   //!
		TBranch        *b_D0KPZ;   //!
		TBranch        *b_D0PIIPCHI2;   //!
		TBranch        *b_D0KIPCHI2;   //!
		TBranch        *b_D0TRK0;   //!
		TBranch        *b_D0TRK1;   //!
		TBranch        *b_D0TRUETRK0;   //!
		TBranch        *b_D0TRUETRK1;   //!
		TBranch        *b_D0NJ;   //!
		TBranch        *b_D0MAXDR;   //!
		TBranch        *b_D0TRUE;   //!
		TBranch        *b_D0TRUEDR;   //!
		TBranch        *b_D0TRUEIDX;   //!
		TBranch        *b_D0FROMB;   //!
		TBranch        *b_DM;   //!
		TBranch        *b_DPX;   //!
		TBranch        *b_DPY;   //!
		TBranch        *b_DPZ;   //!
		TBranch        *b_DE;   //!
		TBranch        *b_DX;   //!
		TBranch        *b_DY;   //!
		TBranch        *b_DZ;   //!
		TBranch        *b_DIP;   //!
		TBranch        *b_DIPCHI2;   //!
		TBranch        *b_DFD;   //!
		TBranch        *b_DFDCHI2;   //!
		TBranch        *b_DTAU;   //!
		TBranch        *b_DDIRA;   //!
		TBranch        *b_DVTXCHI2;   //!
		TBranch        *b_DVTXNDOF;   //!
		TBranch        *b_DDRJET;   //!
		TBranch        *b_DDRTAG;   //!
		TBranch        *b_DPI1PNN;   //!
		TBranch        *b_DPI2PNN;   //!
		TBranch        *b_DKPNN;   //!
		TBranch        *b_DP;   //!
		TBranch        *b_DPT;   //!
		TBranch        *b_DPI1PT;   //!
		TBranch        *b_DPI2PT;   //!
		TBranch        *b_DKPT;   //!
		TBranch        *b_DTRK0;   //!
		TBranch        *b_DTRK1;   //!
		TBranch        *b_DTRK2;   //!
		TBranch        *b_DNJ;   //!
		TBranch        *b_DMAXDR;   //!
		TBranch        *b_DTRUE;   //!
		TBranch        *b_DTRUEDR;   //!
		TBranch        *b_DTRUEIDX;   //!
		TBranch        *b_DFROMB;   //!
		TBranch        *b_DSM;   //!
		TBranch        *b_DSPX;   //!
		TBranch        *b_DSPY;   //!
		TBranch        *b_DSPZ;   //!
		TBranch        *b_DSE;   //!
		TBranch        *b_DSX;   //!
		TBranch        *b_DSY;   //!
		TBranch        *b_DSZ;   //!
		TBranch        *b_DSIP;   //!
		TBranch        *b_DSIPCHI2;   //!
		TBranch        *b_DSFD;   //!
		TBranch        *b_DSFDCHI2;   //!
		TBranch        *b_DSTAU;   //!
		TBranch        *b_DSDIRA;   //!
		TBranch        *b_DSVTXCHI2;   //!
		TBranch        *b_DSVTXNDOF;   //!
		TBranch        *b_DSDRJET;   //!
		TBranch        *b_DSDRTAG;   //!
		TBranch        *b_DSPIPNN;   //!
		TBranch        *b_DSK1PNN;   //!
		TBranch        *b_DSK2PNN;   //!
		TBranch        *b_DSPHIM;   //!
		TBranch        *b_DSP;   //!
		TBranch        *b_DSPT;   //!
		TBranch        *b_DSPIPT;   //!
		TBranch        *b_DSK1PT;   //!
		TBranch        *b_DSK2PT;   //!
		TBranch        *b_DSPHIPT;   //!
		TBranch        *b_DSTRK0;   //!
		TBranch        *b_DSTRK1;   //!
		TBranch        *b_DSTRK2;   //!
		TBranch        *b_DSNJ;   //!
		TBranch        *b_DSMAXDR;   //!
		TBranch        *b_DSTRUE;   //!
		TBranch        *b_DSTRUEDR;   //!
		TBranch        *b_DSTRUEIDX;   //!
		TBranch        *b_DSFROMB;   //!
		TBranch        *b_LCM;   //!
		TBranch        *b_LCPX;   //!
		TBranch        *b_LCPY;   //!
		TBranch        *b_LCPZ;   //!
		TBranch        *b_LCE;   //!
		TBranch        *b_LCX;   //!
		TBranch        *b_LCY;   //!
		TBranch        *b_LCZ;   //!
		TBranch        *b_LCIP;   //!
		TBranch        *b_LCIPCHI2;   //!
		TBranch        *b_LCFD;   //!
		TBranch        *b_LCFDCHI2;   //!
		TBranch        *b_LCTAU;   //!
		TBranch        *b_LCDIRA;   //!
		TBranch        *b_LCVTXCHI2;   //!
		TBranch        *b_LCVTXNDOF;   //!
		TBranch        *b_LCDRJET;   //!
		TBranch        *b_LCDRTAG;   //!
		TBranch        *b_LCPPNN;   //!
		TBranch        *b_LCKPNN;   //!
		TBranch        *b_LCPIPNN;   //!
		TBranch        *b_LCP;   //!
		TBranch        *b_LCPT;   //!
		TBranch        *b_LCPPT;   //!
		TBranch        *b_LCKPT;   //!
		TBranch        *b_LCPIPT;   //!
		TBranch        *b_LCTRK0;   //!
		TBranch        *b_LCTRK1;   //!
		TBranch        *b_LCTRK2;   //!
		TBranch        *b_LCNJ;   //!
		TBranch        *b_LCMAXDR;   //!
		TBranch        *b_K3PIM;   //!
		TBranch        *b_K3PIPX;   //!
		TBranch        *b_K3PIPY;   //!
		TBranch        *b_K3PIPZ;   //!
		TBranch        *b_K3PIE;   //!
		TBranch        *b_K3PIX;   //!
		TBranch        *b_K3PIY;   //!
		TBranch        *b_K3PIZ;   //!
		TBranch        *b_K3PIIP;   //!
		TBranch        *b_K3PIIPCHI2;   //!
		TBranch        *b_K3PIFD;   //!
		TBranch        *b_K3PIFDCHI2;   //!
		TBranch        *b_K3PITAU;   //!
		TBranch        *b_K3PIDIRA;   //!
		TBranch        *b_K3PIVTXCHI2;   //!
		TBranch        *b_K3PIVTXNDOF;   //!
		TBranch        *b_K3PIDRJET;   //!
		TBranch        *b_K3PIDRTAG;   //!
		TBranch        *b_K3PIPI1PNN;   //!
		TBranch        *b_K3PIPI2PNN;   //!
		TBranch        *b_K3PIPI3PNN;   //!
		TBranch        *b_K3PIKPNN;   //!
		TBranch        *b_K3PIP;   //!
		TBranch        *b_K3PIPT;   //!
		TBranch        *b_K3PIPI1PT;   //!
		TBranch        *b_K3PIPI2PT;   //!
		TBranch        *b_K3PIPI3PT;   //!
		TBranch        *b_K3PIKPT;   //!
		TBranch        *b_K3PITRK0;   //!
		TBranch        *b_K3PITRK1;   //!
		TBranch        *b_K3PITRK2;   //!
		TBranch        *b_K3PITRK3;   //!
		TBranch        *b_K3PINJ;   //!
		TBranch        *b_K3PIMAXDR;   //!
		TBranch        *b_JPSIM;   //!
		TBranch        *b_JPSIPX;   //!
		TBranch        *b_JPSIPY;   //!
		TBranch        *b_JPSIPZ;   //!
		TBranch        *b_JPSIE;   //!
		TBranch        *b_JPSIX;   //!
		TBranch        *b_JPSIY;   //!
		TBranch        *b_JPSIZ;   //!
		TBranch        *b_JPSIIP;   //!
		TBranch        *b_JPSIIPCHI2;   //!
		TBranch        *b_JPSIFD;   //!
		TBranch        *b_JPSIFDCHI2;   //!
		TBranch        *b_JPSITAU;   //!
		TBranch        *b_JPSIDIRA;   //!
		TBranch        *b_JPSIVTXCHI2;   //!
		TBranch        *b_JPSIVTXNDOF;   //!
		TBranch        *b_JPSIDRJET;   //!
		TBranch        *b_JPSIDRTAG;   //!
		TBranch        *b_JPSIP;   //!
		TBranch        *b_JPSIPT;   //!
		TBranch        *b_JPSIPIPT;   //!
		TBranch        *b_JPSIKPT;   //!
		TBranch        *b_JPSITRK0;   //!
		TBranch        *b_JPSITRK1;   //!
		TBranch        *b_JPSINJ;   //!
		TBranch        *b_JPSIMAXDR;   //!
		TBranch        *b_JPSITRUE;   //!
		TBranch        *b_JPSITRUEDR;   //!
		TBranch        *b_JPSITRUEIDX;   //!
		TBranch        *b_JPSIFROMB;   //!

		makeNewD0Effs(TString sample="2XX", Binning binning=ThirtySevenByEight);
		virtual ~makeNewD0Effs();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		Binning binning_;
		TString sample_;

		//TH2F* pidPi;
		TH2F* pidK;
};

#endif

#ifdef makeNewD0Effs_cxx
makeNewD0Effs::makeNewD0Effs(TString sample, Binning binning) : fChain(0), binning_(binning), sample_(sample)
{
	TTree* tree(0);
	TFile* f = new TFile("/eos/user/d/dcraik/jets-tuples-new-181203/for_yandex_data_new_"+sample_+".root");
	f->GetObject("T",tree);
	Init(tree);

	TFile* fpidk = TFile::Open("pidcalib/PerfHists_K_Turbo16_MagDown_kaons4_Brunel_P_Brunel_PT.root");
	//TFile* fpidpi= TFile::Open("pidcalib/PerfHists_Pi_Turbo16_MagDown_kaons4_Brunel_P_Brunel_PT.root");
	pidK = dynamic_cast<TH2F*>(fpidk->Get("K_Brunel_MC15TuneV1_ProbNNK > 0.2_All"));
	//pidPi= dynamic_cast<TH2F*>(fpidpi->Get("Pi_Brunel_MC15TuneV1_ProbNNpi > 0.1_All"));
}

makeNewD0Effs::~makeNewD0Effs()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t makeNewD0Effs::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t makeNewD0Effs::LoadTree(Long64_t entry)
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

void makeNewD0Effs::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	SVX = 0;
	SVY = 0;
	SVZ = 0;
	SVPerp = 0;
	SVPx = 0;
	SVPy = 0;
	SVPz = 0;
	SVE = 0;
	SVPT = 0;
	SVETA = 0;
	SVM = 0;
	SVMCor = 0;
	SVMCorErr = 0;
	SVMINPERP = 0;
	SVDRJ = 0;
	SVDRT = 0;
	SVN = 0;
	SVNTRUE = 0;
	SVNJ = 0;
	SVNT = 0;
	SVQ = 0;
	SVSumIPChi2 = 0;
	SVTZ = 0;
	SVIPCHI2 = 0;
	SVMINIPCHI2 = 0;
	SVGhostMax = 0;
	SVISD0 = 0;
	SVISDP = 0;
	SVD0M = 0;
	SVDPM = 0;
	SVTRUEIDX = 0;
	SVTRUETRK0IDX = 0;
	SVTRUETRK1IDX = 0;
	SVTRUETRK2IDX = 0;
	SVTRUETRK3IDX = 0;
	SVTRK0P = 0;
	SVTRK1P = 0;
	SVTRK2P = 0;
	SVTRK3P = 0;
	SVTRK0PT = 0;
	SVTRK1PT = 0;
	SVTRK2PT = 0;
	SVTRK3PT = 0;
	SVTRK0PX = 0;
	SVTRK1PX = 0;
	SVTRK2PX = 0;
	SVTRK3PX = 0;
	SVTRK0PY = 0;
	SVTRK1PY = 0;
	SVTRK2PY = 0;
	SVTRK3PY = 0;
	SVTRK0PZ = 0;
	SVTRK1PZ = 0;
	SVTRK2PZ = 0;
	SVTRK3PZ = 0;
	SVTRK0PNNPI = 0;
	SVTRK1PNNPI = 0;
	SVTRK2PNNPI = 0;
	SVTRK3PNNPI = 0;
	SVTRK0PNNK = 0;
	SVTRK1PNNK = 0;
	SVTRK2PNNK = 0;
	SVTRK3PNNK = 0;
	TSVX = 0;
	TSVY = 0;
	TSVZ = 0;
	TSVPerp = 0;
	TSVPx = 0;
	TSVPy = 0;
	TSVPz = 0;
	TSVE = 0;
	TSVPT = 0;
	TSVETA = 0;
	TSVM = 0;
	TSVMCor = 0;
	TSVMCorErr = 0;
	TSVMINPERP = 0;
	TSVDRJ = 0;
	TSVDRT = 0;
	TSVN = 0;
	TSVNTRUE = 0;
	TSVNJ = 0;
	TSVNT = 0;
	TSVQ = 0;
	TSVSumIPChi2 = 0;
	TSVTZ = 0;
	TSVIPCHI2 = 0;
	TSVMINIPCHI2 = 0;
	TSVGhostMax = 0;
	TSVISD0 = 0;
	TSVISDP = 0;
	TSVD0M = 0;
	TSVDPM = 0;
	TSVTRUEIDX = 0;
	TRUEDID = 0;
	TRUEDPX = 0;
	TRUEDPY = 0;
	TRUEDPZ = 0;
	TRUEDE = 0;
	TRUEDFROMB = 0;
	TRUEDSEL = 0;
	TRUEDTRK0IDX = 0;
	TRUEDTRK0ID = 0;
	TRUEDTRK0P = 0;
	TRUEDTRK0PT = 0;
	TRUEDTRK0INACC = 0;
	TRUEDTRK0RECO = 0;
	TRUEDTRK0PNNK = 0;
	TRUEDTRK0PNNPI = 0;
	TRUEDTRK1IDX = 0;
	TRUEDTRK1ID = 0;
	TRUEDTRK1P = 0;
	TRUEDTRK1PT = 0;
	TRUEDTRK1INACC = 0;
	TRUEDTRK1RECO = 0;
	TRUEDTRK1PNNK = 0;
	TRUEDTRK1PNNPI = 0;
	TRUEDTRK2IDX = 0;
	TRUEDTRK2ID = 0;
	TRUEDTRK2P = 0;
	TRUEDTRK2PT = 0;
	TRUEDTRK2INACC = 0;
	TRUEDTRK2RECO = 0;
	TRUEDTRK2PNNK = 0;
	TRUEDTRK2PNNPI = 0;
	TRUEDTRK3IDX = 0;
	TRUEDTRK3ID = 0;
	TRUEDTRK3P = 0;
	TRUEDTRK3PT = 0;
	TRUEDTRK3INACC = 0;
	TRUEDTRK3RECO = 0;
	TRUEDTRK3PNNK = 0;
	TRUEDTRK3PNNPI = 0;
	TRUEDTRK0PX = 0;
	TRUEDTRK1PX = 0;
	TRUEDTRK2PX = 0;
	TRUEDTRK3PX = 0;
	TRUEDTRK0PY = 0;
	TRUEDTRK1PY = 0;
	TRUEDTRK2PY = 0;
	TRUEDTRK3PY = 0;
	TRUEDTRK0PZ = 0;
	TRUEDTRK1PZ = 0;
	TRUEDTRK2PZ = 0;
	TRUEDTRK3PZ = 0;
	D0M = 0;
	D0PX = 0;
	D0PY = 0;
	D0PZ = 0;
	D0E = 0;
	D0X = 0;
	D0Y = 0;
	D0Z = 0;
	D0IP = 0;
	D0IPCHI2 = 0;
	D0FD = 0;
	D0FDCHI2 = 0;
	D0TAU = 0;
	D0DIRA = 0;
	D0VTXCHI2 = 0;
	D0VTXNDOF = 0;
	D0DRJET = 0;
	D0DRTAG = 0;
	D0PIPNNK = 0;
	D0KPNNK = 0;
	D0PIPNNPI = 0;
	D0KPNNPI = 0;
	D0P = 0;
	D0PT = 0;
	D0ETA = 0;
	D0PIP = 0;
	D0KP = 0;
	D0PIPT = 0;
	D0KPT = 0;
	D0PIETA = 0;
	D0KETA = 0;
	D0PIPX = 0;
	D0KPX = 0;
	D0PIPY = 0;
	D0KPY = 0;
	D0PIPZ = 0;
	D0KPZ = 0;
	D0PIIPCHI2 = 0;
	D0KIPCHI2 = 0;
	D0TRK0 = 0;
	D0TRK1 = 0;
	D0TRUETRK0 = 0;
	D0TRUETRK1 = 0;
	D0NJ = 0;
	D0MAXDR = 0;
	D0TRUE = 0;
	D0TRUEDR = 0;
	D0TRUEIDX = 0;
	D0FROMB = 0;
	DM = 0;
	DPX = 0;
	DPY = 0;
	DPZ = 0;
	DE = 0;
	DX = 0;
	DY = 0;
	DZ = 0;
	DIP = 0;
	DIPCHI2 = 0;
	DFD = 0;
	DFDCHI2 = 0;
	DTAU = 0;
	DDIRA = 0;
	DVTXCHI2 = 0;
	DVTXNDOF = 0;
	DDRJET = 0;
	DDRTAG = 0;
	DPI1PNN = 0;
	DPI2PNN = 0;
	DKPNN = 0;
	DP = 0;
	DPT = 0;
	DPI1PT = 0;
	DPI2PT = 0;
	DKPT = 0;
	DTRK0 = 0;
	DTRK1 = 0;
	DTRK2 = 0;
	DNJ = 0;
	DMAXDR = 0;
	DTRUE = 0;
	DTRUEDR = 0;
	DTRUEIDX = 0;
	DFROMB = 0;
	DSM = 0;
	DSPX = 0;
	DSPY = 0;
	DSPZ = 0;
	DSE = 0;
	DSX = 0;
	DSY = 0;
	DSZ = 0;
	DSIP = 0;
	DSIPCHI2 = 0;
	DSFD = 0;
	DSFDCHI2 = 0;
	DSTAU = 0;
	DSDIRA = 0;
	DSVTXCHI2 = 0;
	DSVTXNDOF = 0;
	DSDRJET = 0;
	DSDRTAG = 0;
	DSPIPNN = 0;
	DSK1PNN = 0;
	DSK2PNN = 0;
	DSPHIM = 0;
	DSP = 0;
	DSPT = 0;
	DSPIPT = 0;
	DSK1PT = 0;
	DSK2PT = 0;
	DSPHIPT = 0;
	DSTRK0 = 0;
	DSTRK1 = 0;
	DSTRK2 = 0;
	DSNJ = 0;
	DSMAXDR = 0;
	DSTRUE = 0;
	DSTRUEDR = 0;
	DSTRUEIDX = 0;
	DSFROMB = 0;
	LCM = 0;
	LCPX = 0;
	LCPY = 0;
	LCPZ = 0;
	LCE = 0;
	LCX = 0;
	LCY = 0;
	LCZ = 0;
	LCIP = 0;
	LCIPCHI2 = 0;
	LCFD = 0;
	LCFDCHI2 = 0;
	LCTAU = 0;
	LCDIRA = 0;
	LCVTXCHI2 = 0;
	LCVTXNDOF = 0;
	LCDRJET = 0;
	LCDRTAG = 0;
	LCPPNN = 0;
	LCKPNN = 0;
	LCPIPNN = 0;
	LCP = 0;
	LCPT = 0;
	LCPPT = 0;
	LCKPT = 0;
	LCPIPT = 0;
	LCTRK0 = 0;
	LCTRK1 = 0;
	LCTRK2 = 0;
	LCNJ = 0;
	LCMAXDR = 0;
	K3PIM = 0;
	K3PIPX = 0;
	K3PIPY = 0;
	K3PIPZ = 0;
	K3PIE = 0;
	K3PIX = 0;
	K3PIY = 0;
	K3PIZ = 0;
	K3PIIP = 0;
	K3PIIPCHI2 = 0;
	K3PIFD = 0;
	K3PIFDCHI2 = 0;
	K3PITAU = 0;
	K3PIDIRA = 0;
	K3PIVTXCHI2 = 0;
	K3PIVTXNDOF = 0;
	K3PIDRJET = 0;
	K3PIDRTAG = 0;
	K3PIPI1PNN = 0;
	K3PIPI2PNN = 0;
	K3PIPI3PNN = 0;
	K3PIKPNN = 0;
	K3PIP = 0;
	K3PIPT = 0;
	K3PIPI1PT = 0;
	K3PIPI2PT = 0;
	K3PIPI3PT = 0;
	K3PIKPT = 0;
	K3PITRK0 = 0;
	K3PITRK1 = 0;
	K3PITRK2 = 0;
	K3PITRK3 = 0;
	K3PINJ = 0;
	K3PIMAXDR = 0;
	JPSIM = 0;
	JPSIPX = 0;
	JPSIPY = 0;
	JPSIPZ = 0;
	JPSIE = 0;
	JPSIX = 0;
	JPSIY = 0;
	JPSIZ = 0;
	JPSIIP = 0;
	JPSIIPCHI2 = 0;
	JPSIFD = 0;
	JPSIFDCHI2 = 0;
	JPSITAU = 0;
	JPSIDIRA = 0;
	JPSIVTXCHI2 = 0;
	JPSIVTXNDOF = 0;
	JPSIDRJET = 0;
	JPSIDRTAG = 0;
	JPSIP = 0;
	JPSIPT = 0;
	JPSIPIPT = 0;
	JPSIKPT = 0;
	JPSITRK0 = 0;
	JPSITRK1 = 0;
	JPSINJ = 0;
	JPSIMAXDR = 0;
	JPSITRUE = 0;
	JPSITRUEDR = 0;
	JPSITRUEIDX = 0;
	JPSIFROMB = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
	fChain->SetBranchAddress("Dec", &Dec, &b_Dec);
	fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
	fChain->SetBranchAddress("JetPx", &JetPx, &b_JetPx);
	fChain->SetBranchAddress("JetPy", &JetPy, &b_JetPy);
	fChain->SetBranchAddress("JetPz", &JetPz, &b_JetPz);
	fChain->SetBranchAddress("JetE", &JetE, &b_JetE);
	fChain->SetBranchAddress("JetPT", &JetPT, &b_JetPT);
	fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
	fChain->SetBranchAddress("JetTruePx", &JetTruePx, &b_JetTruePx);
	fChain->SetBranchAddress("JetTruePy", &JetTruePy, &b_JetTruePy);
	fChain->SetBranchAddress("JetTruePz", &JetTruePz, &b_JetTruePz);
	fChain->SetBranchAddress("JetTrueE", &JetTrueE, &b_JetTrueE);
	fChain->SetBranchAddress("JetTruePT", &JetTruePT, &b_JetTruePT);
	fChain->SetBranchAddress("JetTrueEta", &JetTrueEta, &b_JetTrueEta);
	fChain->SetBranchAddress("JetTrueDR", &JetTrueDR, &b_JetTrueDR);
	fChain->SetBranchAddress("JetSigma1", &JetSigma1, &b_JetSigma1);
	fChain->SetBranchAddress("JetSigma2", &JetSigma2, &b_JetSigma2);
	fChain->SetBranchAddress("JetQ", &JetQ, &b_JetQ);
	fChain->SetBranchAddress("JetMult", &JetMult, &b_JetMult);
	fChain->SetBranchAddress("JetNChr", &JetNChr, &b_JetNChr);
	fChain->SetBranchAddress("JetNNeu", &JetNNeu, &b_JetNNeu);
	fChain->SetBranchAddress("JetPTD", &JetPTD, &b_JetPTD);
	fChain->SetBranchAddress("JetTRUEb", &JetTRUEb, &b_JetTRUEb);
	fChain->SetBranchAddress("JetTRUEc", &JetTRUEc, &b_JetTRUEc);
	fChain->SetBranchAddress("JetTRUED0", &JetTRUED0, &b_JetTRUED0);
	fChain->SetBranchAddress("JetTRUEDP", &JetTRUEDP, &b_JetTRUEDP);
	fChain->SetBranchAddress("JetTRUEDS", &JetTRUEDS, &b_JetTRUEDS);
	fChain->SetBranchAddress("JetTRUEJPSI", &JetTRUEJPSI, &b_JetTRUEJPSI);
	fChain->SetBranchAddress("JetTRUEDSV", &JetTRUEDSV, &b_JetTRUEDSV);
	fChain->SetBranchAddress("JetTRUEBSV", &JetTRUEBSV, &b_JetTRUEBSV);
	fChain->SetBranchAddress("JetTRUESV", &JetTRUESV, &b_JetTRUESV);
	fChain->SetBranchAddress("JetTRUEDPX", &JetTRUEDPX, &b_JetTRUEDPX);
	fChain->SetBranchAddress("JetTRUEDPY", &JetTRUEDPY, &b_JetTRUEDPY);
	fChain->SetBranchAddress("JetTRUEDPZ", &JetTRUEDPZ, &b_JetTRUEDPZ);
	fChain->SetBranchAddress("JetTRUEDE", &JetTRUEDE, &b_JetTRUEDE);
	fChain->SetBranchAddress("JetTRUEDDR", &JetTRUEDDR, &b_JetTRUEDDR);
	fChain->SetBranchAddress("TagPx", &TagPx, &b_TagPx);
	fChain->SetBranchAddress("TagPy", &TagPy, &b_TagPy);
	fChain->SetBranchAddress("TagPz", &TagPz, &b_TagPz);
	fChain->SetBranchAddress("TagE", &TagE, &b_TagE);
	fChain->SetBranchAddress("TagPT", &TagPT, &b_TagPT);
	fChain->SetBranchAddress("TagEta", &TagEta, &b_TagEta);
	fChain->SetBranchAddress("DeltaR", &DeltaR, &b_DeltaR);
	fChain->SetBranchAddress("DeltaPhi", &DeltaPhi, &b_DeltaPhi);
	fChain->SetBranchAddress("SVX", &SVX, &b_SVX);
	fChain->SetBranchAddress("SVY", &SVY, &b_SVY);
	fChain->SetBranchAddress("SVZ", &SVZ, &b_SVZ);
	fChain->SetBranchAddress("SVPerp", &SVPerp, &b_SVPerp);
	fChain->SetBranchAddress("SVPx", &SVPx, &b_SVPx);
	fChain->SetBranchAddress("SVPy", &SVPy, &b_SVPy);
	fChain->SetBranchAddress("SVPz", &SVPz, &b_SVPz);
	fChain->SetBranchAddress("SVE", &SVE, &b_SVE);
	fChain->SetBranchAddress("SVPT", &SVPT, &b_SVPT);
	fChain->SetBranchAddress("SVETA", &SVETA, &b_SVETA);
	fChain->SetBranchAddress("SVM", &SVM, &b_SVM);
	fChain->SetBranchAddress("SVMCor", &SVMCor, &b_SVMCor);
	fChain->SetBranchAddress("SVMCorErr", &SVMCorErr, &b_SVMCorErr);
	fChain->SetBranchAddress("SVMINPERP", &SVMINPERP, &b_SVMINPERP);
	fChain->SetBranchAddress("SVDRJ", &SVDRJ, &b_SVDRJ);
	fChain->SetBranchAddress("SVDRT", &SVDRT, &b_SVDRT);
	fChain->SetBranchAddress("SVN", &SVN, &b_SVN);
	fChain->SetBranchAddress("SVNTRUE", &SVNTRUE, &b_SVNTRUE);
	fChain->SetBranchAddress("SVNJ", &SVNJ, &b_SVNJ);
	fChain->SetBranchAddress("SVNT", &SVNT, &b_SVNT);
	fChain->SetBranchAddress("SVQ", &SVQ, &b_SVQ);
	fChain->SetBranchAddress("SVSumIPChi2", &SVSumIPChi2, &b_SVSumIPChi2);
	fChain->SetBranchAddress("SVTZ", &SVTZ, &b_SVTZ);
	fChain->SetBranchAddress("SVIPCHI2", &SVIPCHI2, &b_SVIPCHI2);
	fChain->SetBranchAddress("SVMINIPCHI2", &SVMINIPCHI2, &b_SVMINIPCHI2);
	fChain->SetBranchAddress("SVGhostMax", &SVGhostMax, &b_SVGhostMax);
	fChain->SetBranchAddress("SVISD0", &SVISD0, &b_SVISD0);
	fChain->SetBranchAddress("SVISDP", &SVISDP, &b_SVISDP);
	fChain->SetBranchAddress("SVD0M", &SVD0M, &b_SVD0M);
	fChain->SetBranchAddress("SVDPM", &SVDPM, &b_SVDPM);
	fChain->SetBranchAddress("SVTRUEIDX", &SVTRUEIDX, &b_SVTRUEIDX);
	fChain->SetBranchAddress("SVTRUETRK0IDX", &SVTRUETRK0IDX, &b_SVTRUETRK0IDX);
	fChain->SetBranchAddress("SVTRUETRK1IDX", &SVTRUETRK1IDX, &b_SVTRUETRK1IDX);
	fChain->SetBranchAddress("SVTRUETRK2IDX", &SVTRUETRK2IDX, &b_SVTRUETRK2IDX);
	fChain->SetBranchAddress("SVTRUETRK3IDX", &SVTRUETRK3IDX, &b_SVTRUETRK3IDX);
	fChain->SetBranchAddress("SVTRK0P", &SVTRK0P, &b_SVTRK0P);
	fChain->SetBranchAddress("SVTRK1P", &SVTRK1P, &b_SVTRK1P);
	fChain->SetBranchAddress("SVTRK2P", &SVTRK2P, &b_SVTRK2P);
	fChain->SetBranchAddress("SVTRK3P", &SVTRK3P, &b_SVTRK3P);
	fChain->SetBranchAddress("SVTRK0PT", &SVTRK0PT, &b_SVTRK0PT);
	fChain->SetBranchAddress("SVTRK1PT", &SVTRK1PT, &b_SVTRK1PT);
	fChain->SetBranchAddress("SVTRK2PT", &SVTRK2PT, &b_SVTRK2PT);
	fChain->SetBranchAddress("SVTRK3PT", &SVTRK3PT, &b_SVTRK3PT);
	fChain->SetBranchAddress("SVTRK0PX", &SVTRK0PX, &b_SVTRK0PX);
	fChain->SetBranchAddress("SVTRK1PX", &SVTRK1PX, &b_SVTRK1PX);
	fChain->SetBranchAddress("SVTRK2PX", &SVTRK2PX, &b_SVTRK2PX);
	fChain->SetBranchAddress("SVTRK3PX", &SVTRK3PX, &b_SVTRK3PX);
	fChain->SetBranchAddress("SVTRK0PY", &SVTRK0PY, &b_SVTRK0PY);
	fChain->SetBranchAddress("SVTRK1PY", &SVTRK1PY, &b_SVTRK1PY);
	fChain->SetBranchAddress("SVTRK2PY", &SVTRK2PY, &b_SVTRK2PY);
	fChain->SetBranchAddress("SVTRK3PY", &SVTRK3PY, &b_SVTRK3PY);
	fChain->SetBranchAddress("SVTRK0PZ", &SVTRK0PZ, &b_SVTRK0PZ);
	fChain->SetBranchAddress("SVTRK1PZ", &SVTRK1PZ, &b_SVTRK1PZ);
	fChain->SetBranchAddress("SVTRK2PZ", &SVTRK2PZ, &b_SVTRK2PZ);
	fChain->SetBranchAddress("SVTRK3PZ", &SVTRK3PZ, &b_SVTRK3PZ);
	fChain->SetBranchAddress("SVTRK0PNNPI", &SVTRK0PNNPI, &b_SVTRK0PNNPI);
	fChain->SetBranchAddress("SVTRK1PNNPI", &SVTRK1PNNPI, &b_SVTRK1PNNPI);
	fChain->SetBranchAddress("SVTRK2PNNPI", &SVTRK2PNNPI, &b_SVTRK2PNNPI);
	fChain->SetBranchAddress("SVTRK3PNNPI", &SVTRK3PNNPI, &b_SVTRK3PNNPI);
	fChain->SetBranchAddress("SVTRK0PNNK", &SVTRK0PNNK, &b_SVTRK0PNNK);
	fChain->SetBranchAddress("SVTRK1PNNK", &SVTRK1PNNK, &b_SVTRK1PNNK);
	fChain->SetBranchAddress("SVTRK2PNNK", &SVTRK2PNNK, &b_SVTRK2PNNK);
	fChain->SetBranchAddress("SVTRK3PNNK", &SVTRK3PNNK, &b_SVTRK3PNNK);
	fChain->SetBranchAddress("TSVX", &TSVX, &b_TSVX);
	fChain->SetBranchAddress("TSVY", &TSVY, &b_TSVY);
	fChain->SetBranchAddress("TSVZ", &TSVZ, &b_TSVZ);
	fChain->SetBranchAddress("TSVPerp", &TSVPerp, &b_TSVPerp);
	fChain->SetBranchAddress("TSVPx", &TSVPx, &b_TSVPx);
	fChain->SetBranchAddress("TSVPy", &TSVPy, &b_TSVPy);
	fChain->SetBranchAddress("TSVPz", &TSVPz, &b_TSVPz);
	fChain->SetBranchAddress("TSVE", &TSVE, &b_TSVE);
	fChain->SetBranchAddress("TSVPT", &TSVPT, &b_TSVPT);
	fChain->SetBranchAddress("TSVETA", &TSVETA, &b_TSVETA);
	fChain->SetBranchAddress("TSVM", &TSVM, &b_TSVM);
	fChain->SetBranchAddress("TSVMCor", &TSVMCor, &b_TSVMCor);
	fChain->SetBranchAddress("TSVMCorErr", &TSVMCorErr, &b_TSVMCorErr);
	fChain->SetBranchAddress("TSVMINPERP", &TSVMINPERP, &b_TSVMINPERP);
	fChain->SetBranchAddress("TSVDRJ", &TSVDRJ, &b_TSVDRJ);
	fChain->SetBranchAddress("TSVDRT", &TSVDRT, &b_TSVDRT);
	fChain->SetBranchAddress("TSVN", &TSVN, &b_TSVN);
	fChain->SetBranchAddress("TSVNTRUE", &TSVNTRUE, &b_TSVNTRUE);
	fChain->SetBranchAddress("TSVNJ", &TSVNJ, &b_TSVNJ);
	fChain->SetBranchAddress("TSVNT", &TSVNT, &b_TSVNT);
	fChain->SetBranchAddress("TSVQ", &TSVQ, &b_TSVQ);
	fChain->SetBranchAddress("TSVSumIPChi2", &TSVSumIPChi2, &b_TSVSumIPChi2);
	fChain->SetBranchAddress("TSVTZ", &TSVTZ, &b_TSVTZ);
	fChain->SetBranchAddress("TSVIPCHI2", &TSVIPCHI2, &b_TSVIPCHI2);
	fChain->SetBranchAddress("TSVMINIPCHI2", &TSVMINIPCHI2, &b_TSVMINIPCHI2);
	fChain->SetBranchAddress("TSVGhostMax", &TSVGhostMax, &b_TSVGhostMax);
	fChain->SetBranchAddress("TSVISD0", &TSVISD0, &b_TSVISD0);
	fChain->SetBranchAddress("TSVISDP", &TSVISDP, &b_TSVISDP);
	fChain->SetBranchAddress("TSVD0M", &TSVD0M, &b_TSVD0M);
	fChain->SetBranchAddress("TSVDPM", &TSVDPM, &b_TSVDPM);
	fChain->SetBranchAddress("TSVTRUEIDX", &TSVTRUEIDX, &b_TSVTRUEIDX);
	fChain->SetBranchAddress("NSV", &NSV, &b_NSV);
	fChain->SetBranchAddress("NTSV", &NTSV, &b_NTSV);
	fChain->SetBranchAddress("NTRK", &NTRK, &b_NTRK);
	fChain->SetBranchAddress("NNEU", &NNEU, &b_NNEU);
	fChain->SetBranchAddress("PVX", &PVX, &b_PVX);
	fChain->SetBranchAddress("PVY", &PVY, &b_PVY);
	fChain->SetBranchAddress("PVZ", &PVZ, &b_PVZ);
	fChain->SetBranchAddress("NDispl6", &NDispl6, &b_NDispl6);
	fChain->SetBranchAddress("NDispl9", &NDispl9, &b_NDispl9);
	fChain->SetBranchAddress("NDispl16", &NDispl16, &b_NDispl16);
	fChain->SetBranchAddress("MuPT", &MuPT, &b_MuPT);
	fChain->SetBranchAddress("MuIPChi2", &MuIPChi2, &b_MuIPChi2);
	fChain->SetBranchAddress("MuDR", &MuDR, &b_MuDR);
	fChain->SetBranchAddress("MuPNN", &MuPNN, &b_MuPNN);
	fChain->SetBranchAddress("NMu", &NMu, &b_NMu);
	fChain->SetBranchAddress("HardPT", &HardPT, &b_HardPT);
	fChain->SetBranchAddress("HardIPChi2", &HardIPChi2, &b_HardIPChi2);
	fChain->SetBranchAddress("HardDR", &HardDR, &b_HardDR);
	fChain->SetBranchAddress("TRUEDID", &TRUEDID, &b_TRUEDID);
	fChain->SetBranchAddress("TRUEDPX", &TRUEDPX, &b_TRUEDPX);
	fChain->SetBranchAddress("TRUEDPY", &TRUEDPY, &b_TRUEDPY);
	fChain->SetBranchAddress("TRUEDPZ", &TRUEDPZ, &b_TRUEDPZ);
	fChain->SetBranchAddress("TRUEDE", &TRUEDE, &b_TRUEDE);
	fChain->SetBranchAddress("TRUEDFROMB", &TRUEDFROMB, &b_TRUEDFROMB);
	fChain->SetBranchAddress("TRUEDSEL", &TRUEDSEL, &b_TRUEDSEL);
	fChain->SetBranchAddress("TRUEDTRK0IDX", &TRUEDTRK0IDX, &b_TRUEDTRK0IDX);
	fChain->SetBranchAddress("TRUEDTRK0ID", &TRUEDTRK0ID, &b_TRUEDTRK0ID);
	fChain->SetBranchAddress("TRUEDTRK0P", &TRUEDTRK0P, &b_TRUEDTRK0P);
	fChain->SetBranchAddress("TRUEDTRK0PT", &TRUEDTRK0PT, &b_TRUEDTRK0PT);
	fChain->SetBranchAddress("TRUEDTRK0INACC", &TRUEDTRK0INACC, &b_TRUEDTRK0INACC);
	fChain->SetBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO, &b_TRUEDTRK0RECO);
	fChain->SetBranchAddress("TRUEDTRK0PNNK", &TRUEDTRK0PNNK, &b_TRUEDTRK0PNNK);
	fChain->SetBranchAddress("TRUEDTRK0PNNPI", &TRUEDTRK0PNNPI, &b_TRUEDTRK0PNNPI);
	fChain->SetBranchAddress("TRUEDTRK1IDX", &TRUEDTRK1IDX, &b_TRUEDTRK1IDX);
	fChain->SetBranchAddress("TRUEDTRK1ID", &TRUEDTRK1ID, &b_TRUEDTRK1ID);
	fChain->SetBranchAddress("TRUEDTRK1P", &TRUEDTRK1P, &b_TRUEDTRK1P);
	fChain->SetBranchAddress("TRUEDTRK1PT", &TRUEDTRK1PT, &b_TRUEDTRK1PT);
	fChain->SetBranchAddress("TRUEDTRK1INACC", &TRUEDTRK1INACC, &b_TRUEDTRK1INACC);
	fChain->SetBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO, &b_TRUEDTRK1RECO);
	fChain->SetBranchAddress("TRUEDTRK1PNNK", &TRUEDTRK1PNNK, &b_TRUEDTRK1PNNK);
	fChain->SetBranchAddress("TRUEDTRK1PNNPI", &TRUEDTRK1PNNPI, &b_TRUEDTRK1PNNPI);
	fChain->SetBranchAddress("TRUEDTRK2IDX", &TRUEDTRK2IDX, &b_TRUEDTRK2IDX);
	fChain->SetBranchAddress("TRUEDTRK2ID", &TRUEDTRK2ID, &b_TRUEDTRK2ID);
	fChain->SetBranchAddress("TRUEDTRK2P", &TRUEDTRK2P, &b_TRUEDTRK2P);
	fChain->SetBranchAddress("TRUEDTRK2PT", &TRUEDTRK2PT, &b_TRUEDTRK2PT);
	fChain->SetBranchAddress("TRUEDTRK2INACC", &TRUEDTRK2INACC, &b_TRUEDTRK2INACC);
	fChain->SetBranchAddress("TRUEDTRK2RECO", &TRUEDTRK2RECO, &b_TRUEDTRK2RECO);
	fChain->SetBranchAddress("TRUEDTRK2PNNK", &TRUEDTRK2PNNK, &b_TRUEDTRK2PNNK);
	fChain->SetBranchAddress("TRUEDTRK2PNNPI", &TRUEDTRK2PNNPI, &b_TRUEDTRK2PNNPI);
	fChain->SetBranchAddress("TRUEDTRK3IDX", &TRUEDTRK3IDX, &b_TRUEDTRK3IDX);
	fChain->SetBranchAddress("TRUEDTRK3ID", &TRUEDTRK3ID, &b_TRUEDTRK3ID);
	fChain->SetBranchAddress("TRUEDTRK3P", &TRUEDTRK3P, &b_TRUEDTRK3P);
	fChain->SetBranchAddress("TRUEDTRK3PT", &TRUEDTRK3PT, &b_TRUEDTRK3PT);
	fChain->SetBranchAddress("TRUEDTRK3INACC", &TRUEDTRK3INACC, &b_TRUEDTRK3INACC);
	fChain->SetBranchAddress("TRUEDTRK3RECO", &TRUEDTRK3RECO, &b_TRUEDTRK3RECO);
	fChain->SetBranchAddress("TRUEDTRK3PNNK", &TRUEDTRK3PNNK, &b_TRUEDTRK3PNNK);
	fChain->SetBranchAddress("TRUEDTRK3PNNPI", &TRUEDTRK3PNNPI, &b_TRUEDTRK3PNNPI);
	fChain->SetBranchAddress("TRUEDTRK0PX", &TRUEDTRK0PX, &b_TRUEDTRK0PX);
	fChain->SetBranchAddress("TRUEDTRK1PX", &TRUEDTRK1PX, &b_TRUEDTRK1PX);
	fChain->SetBranchAddress("TRUEDTRK2PX", &TRUEDTRK2PX, &b_TRUEDTRK2PX);
	fChain->SetBranchAddress("TRUEDTRK3PX", &TRUEDTRK3PX, &b_TRUEDTRK3PX);
	fChain->SetBranchAddress("TRUEDTRK0PY", &TRUEDTRK0PY, &b_TRUEDTRK0PY);
	fChain->SetBranchAddress("TRUEDTRK1PY", &TRUEDTRK1PY, &b_TRUEDTRK1PY);
	fChain->SetBranchAddress("TRUEDTRK2PY", &TRUEDTRK2PY, &b_TRUEDTRK2PY);
	fChain->SetBranchAddress("TRUEDTRK3PY", &TRUEDTRK3PY, &b_TRUEDTRK3PY);
	fChain->SetBranchAddress("TRUEDTRK0PZ", &TRUEDTRK0PZ, &b_TRUEDTRK0PZ);
	fChain->SetBranchAddress("TRUEDTRK1PZ", &TRUEDTRK1PZ, &b_TRUEDTRK1PZ);
	fChain->SetBranchAddress("TRUEDTRK2PZ", &TRUEDTRK2PZ, &b_TRUEDTRK2PZ);
	fChain->SetBranchAddress("TRUEDTRK3PZ", &TRUEDTRK3PZ, &b_TRUEDTRK3PZ);
	fChain->SetBranchAddress("D0M", &D0M, &b_D0M);
	fChain->SetBranchAddress("D0PX", &D0PX, &b_D0PX);
	fChain->SetBranchAddress("D0PY", &D0PY, &b_D0PY);
	fChain->SetBranchAddress("D0PZ", &D0PZ, &b_D0PZ);
	fChain->SetBranchAddress("D0E", &D0E, &b_D0E);
	fChain->SetBranchAddress("D0X", &D0X, &b_D0X);
	fChain->SetBranchAddress("D0Y", &D0Y, &b_D0Y);
	fChain->SetBranchAddress("D0Z", &D0Z, &b_D0Z);
	fChain->SetBranchAddress("D0IP", &D0IP, &b_D0IP);
	fChain->SetBranchAddress("D0IPCHI2", &D0IPCHI2, &b_D0IPCHI2);
	fChain->SetBranchAddress("D0FD", &D0FD, &b_D0FD);
	fChain->SetBranchAddress("D0FDCHI2", &D0FDCHI2, &b_D0FDCHI2);
	fChain->SetBranchAddress("D0TAU", &D0TAU, &b_D0TAU);
	fChain->SetBranchAddress("D0DIRA", &D0DIRA, &b_D0DIRA);
	fChain->SetBranchAddress("D0VTXCHI2", &D0VTXCHI2, &b_D0VTXCHI2);
	fChain->SetBranchAddress("D0VTXNDOF", &D0VTXNDOF, &b_D0VTXNDOF);
	fChain->SetBranchAddress("D0DRJET", &D0DRJET, &b_D0DRJET);
	fChain->SetBranchAddress("D0DRTAG", &D0DRTAG, &b_D0DRTAG);
	fChain->SetBranchAddress("D0PIPNNK", &D0PIPNNK, &b_D0PIPNNK);
	fChain->SetBranchAddress("D0KPNNK", &D0KPNNK, &b_D0KPNNK);
	fChain->SetBranchAddress("D0PIPNNPI", &D0PIPNNPI, &b_D0PIPNNPI);
	fChain->SetBranchAddress("D0KPNNPI", &D0KPNNPI, &b_D0KPNNPI);
	fChain->SetBranchAddress("D0P", &D0P, &b_D0P);
	fChain->SetBranchAddress("D0PT", &D0PT, &b_D0PT);
	fChain->SetBranchAddress("D0ETA", &D0ETA, &b_D0ETA);
	fChain->SetBranchAddress("D0PIP", &D0PIP, &b_D0PIP);
	fChain->SetBranchAddress("D0KP", &D0KP, &b_D0KP);
	fChain->SetBranchAddress("D0PIPT", &D0PIPT, &b_D0PIPT);
	fChain->SetBranchAddress("D0KPT", &D0KPT, &b_D0KPT);
	fChain->SetBranchAddress("D0PIETA", &D0PIETA, &b_D0PIETA);
	fChain->SetBranchAddress("D0KETA", &D0KETA, &b_D0KETA);
	fChain->SetBranchAddress("D0PIPX", &D0PIPX, &b_D0PIPX);
	fChain->SetBranchAddress("D0KPX", &D0KPX, &b_D0KPX);
	fChain->SetBranchAddress("D0PIPY", &D0PIPY, &b_D0PIPY);
	fChain->SetBranchAddress("D0KPY", &D0KPY, &b_D0KPY);
	fChain->SetBranchAddress("D0PIPZ", &D0PIPZ, &b_D0PIPZ);
	fChain->SetBranchAddress("D0KPZ", &D0KPZ, &b_D0KPZ);
	fChain->SetBranchAddress("D0PIIPCHI2", &D0PIIPCHI2, &b_D0PIIPCHI2);
	fChain->SetBranchAddress("D0KIPCHI2", &D0KIPCHI2, &b_D0KIPCHI2);
	fChain->SetBranchAddress("D0TRK0", &D0TRK0, &b_D0TRK0);
	fChain->SetBranchAddress("D0TRK1", &D0TRK1, &b_D0TRK1);
	fChain->SetBranchAddress("D0TRUETRK0", &D0TRUETRK0, &b_D0TRUETRK0);
	fChain->SetBranchAddress("D0TRUETRK1", &D0TRUETRK1, &b_D0TRUETRK1);
	fChain->SetBranchAddress("D0NJ", &D0NJ, &b_D0NJ);
	fChain->SetBranchAddress("D0MAXDR", &D0MAXDR, &b_D0MAXDR);
	fChain->SetBranchAddress("D0TRUE", &D0TRUE, &b_D0TRUE);
	fChain->SetBranchAddress("D0TRUEDR", &D0TRUEDR, &b_D0TRUEDR);
	fChain->SetBranchAddress("D0TRUEIDX", &D0TRUEIDX, &b_D0TRUEIDX);
	fChain->SetBranchAddress("D0FROMB", &D0FROMB, &b_D0FROMB);
	fChain->SetBranchAddress("DM", &DM, &b_DM);
	fChain->SetBranchAddress("DPX", &DPX, &b_DPX);
	fChain->SetBranchAddress("DPY", &DPY, &b_DPY);
	fChain->SetBranchAddress("DPZ", &DPZ, &b_DPZ);
	fChain->SetBranchAddress("DE", &DE, &b_DE);
	fChain->SetBranchAddress("DX", &DX, &b_DX);
	fChain->SetBranchAddress("DY", &DY, &b_DY);
	fChain->SetBranchAddress("DZ", &DZ, &b_DZ);
	fChain->SetBranchAddress("DIP", &DIP, &b_DIP);
	fChain->SetBranchAddress("DIPCHI2", &DIPCHI2, &b_DIPCHI2);
	fChain->SetBranchAddress("DFD", &DFD, &b_DFD);
	fChain->SetBranchAddress("DFDCHI2", &DFDCHI2, &b_DFDCHI2);
	fChain->SetBranchAddress("DTAU", &DTAU, &b_DTAU);
	fChain->SetBranchAddress("DDIRA", &DDIRA, &b_DDIRA);
	fChain->SetBranchAddress("DVTXCHI2", &DVTXCHI2, &b_DVTXCHI2);
	fChain->SetBranchAddress("DVTXNDOF", &DVTXNDOF, &b_DVTXNDOF);
	fChain->SetBranchAddress("DDRJET", &DDRJET, &b_DDRJET);
	fChain->SetBranchAddress("DDRTAG", &DDRTAG, &b_DDRTAG);
	fChain->SetBranchAddress("DPI1PNN", &DPI1PNN, &b_DPI1PNN);
	fChain->SetBranchAddress("DPI2PNN", &DPI2PNN, &b_DPI2PNN);
	fChain->SetBranchAddress("DKPNN", &DKPNN, &b_DKPNN);
	fChain->SetBranchAddress("DP", &DP, &b_DP);
	fChain->SetBranchAddress("DPT", &DPT, &b_DPT);
	fChain->SetBranchAddress("DPI1PT", &DPI1PT, &b_DPI1PT);
	fChain->SetBranchAddress("DPI2PT", &DPI2PT, &b_DPI2PT);
	fChain->SetBranchAddress("DKPT", &DKPT, &b_DKPT);
	fChain->SetBranchAddress("DTRK0", &DTRK0, &b_DTRK0);
	fChain->SetBranchAddress("DTRK1", &DTRK1, &b_DTRK1);
	fChain->SetBranchAddress("DTRK2", &DTRK2, &b_DTRK2);
	fChain->SetBranchAddress("DNJ", &DNJ, &b_DNJ);
	fChain->SetBranchAddress("DMAXDR", &DMAXDR, &b_DMAXDR);
	fChain->SetBranchAddress("DTRUE", &DTRUE, &b_DTRUE);
	fChain->SetBranchAddress("DTRUEDR", &DTRUEDR, &b_DTRUEDR);
	fChain->SetBranchAddress("DTRUEIDX", &DTRUEIDX, &b_DTRUEIDX);
	fChain->SetBranchAddress("DFROMB", &DFROMB, &b_DFROMB);
	fChain->SetBranchAddress("DSM", &DSM, &b_DSM);
	fChain->SetBranchAddress("DSPX", &DSPX, &b_DSPX);
	fChain->SetBranchAddress("DSPY", &DSPY, &b_DSPY);
	fChain->SetBranchAddress("DSPZ", &DSPZ, &b_DSPZ);
	fChain->SetBranchAddress("DSE", &DSE, &b_DSE);
	fChain->SetBranchAddress("DSX", &DSX, &b_DSX);
	fChain->SetBranchAddress("DSY", &DSY, &b_DSY);
	fChain->SetBranchAddress("DSZ", &DSZ, &b_DSZ);
	fChain->SetBranchAddress("DSIP", &DSIP, &b_DSIP);
	fChain->SetBranchAddress("DSIPCHI2", &DSIPCHI2, &b_DSIPCHI2);
	fChain->SetBranchAddress("DSFD", &DSFD, &b_DSFD);
	fChain->SetBranchAddress("DSFDCHI2", &DSFDCHI2, &b_DSFDCHI2);
	fChain->SetBranchAddress("DSTAU", &DSTAU, &b_DSTAU);
	fChain->SetBranchAddress("DSDIRA", &DSDIRA, &b_DSDIRA);
	fChain->SetBranchAddress("DSVTXCHI2", &DSVTXCHI2, &b_DSVTXCHI2);
	fChain->SetBranchAddress("DSVTXNDOF", &DSVTXNDOF, &b_DSVTXNDOF);
	fChain->SetBranchAddress("DSDRJET", &DSDRJET, &b_DSDRJET);
	fChain->SetBranchAddress("DSDRTAG", &DSDRTAG, &b_DSDRTAG);
	fChain->SetBranchAddress("DSPIPNN", &DSPIPNN, &b_DSPIPNN);
	fChain->SetBranchAddress("DSK1PNN", &DSK1PNN, &b_DSK1PNN);
	fChain->SetBranchAddress("DSK2PNN", &DSK2PNN, &b_DSK2PNN);
	fChain->SetBranchAddress("DSPHIM", &DSPHIM, &b_DSPHIM);
	fChain->SetBranchAddress("DSP", &DSP, &b_DSP);
	fChain->SetBranchAddress("DSPT", &DSPT, &b_DSPT);
	fChain->SetBranchAddress("DSPIPT", &DSPIPT, &b_DSPIPT);
	fChain->SetBranchAddress("DSK1PT", &DSK1PT, &b_DSK1PT);
	fChain->SetBranchAddress("DSK2PT", &DSK2PT, &b_DSK2PT);
	fChain->SetBranchAddress("DSPHIPT", &DSPHIPT, &b_DSPHIPT);
	fChain->SetBranchAddress("DSTRK0", &DSTRK0, &b_DSTRK0);
	fChain->SetBranchAddress("DSTRK1", &DSTRK1, &b_DSTRK1);
	fChain->SetBranchAddress("DSTRK2", &DSTRK2, &b_DSTRK2);
	fChain->SetBranchAddress("DSNJ", &DSNJ, &b_DSNJ);
	fChain->SetBranchAddress("DSMAXDR", &DSMAXDR, &b_DSMAXDR);
	fChain->SetBranchAddress("DSTRUE", &DSTRUE, &b_DSTRUE);
	fChain->SetBranchAddress("DSTRUEDR", &DSTRUEDR, &b_DSTRUEDR);
	fChain->SetBranchAddress("DSTRUEIDX", &DSTRUEIDX, &b_DSTRUEIDX);
	fChain->SetBranchAddress("DSFROMB", &DSFROMB, &b_DSFROMB);
	fChain->SetBranchAddress("LCM", &LCM, &b_LCM);
	fChain->SetBranchAddress("LCPX", &LCPX, &b_LCPX);
	fChain->SetBranchAddress("LCPY", &LCPY, &b_LCPY);
	fChain->SetBranchAddress("LCPZ", &LCPZ, &b_LCPZ);
	fChain->SetBranchAddress("LCE", &LCE, &b_LCE);
	fChain->SetBranchAddress("LCX", &LCX, &b_LCX);
	fChain->SetBranchAddress("LCY", &LCY, &b_LCY);
	fChain->SetBranchAddress("LCZ", &LCZ, &b_LCZ);
	fChain->SetBranchAddress("LCIP", &LCIP, &b_LCIP);
	fChain->SetBranchAddress("LCIPCHI2", &LCIPCHI2, &b_LCIPCHI2);
	fChain->SetBranchAddress("LCFD", &LCFD, &b_LCFD);
	fChain->SetBranchAddress("LCFDCHI2", &LCFDCHI2, &b_LCFDCHI2);
	fChain->SetBranchAddress("LCTAU", &LCTAU, &b_LCTAU);
	fChain->SetBranchAddress("LCDIRA", &LCDIRA, &b_LCDIRA);
	fChain->SetBranchAddress("LCVTXCHI2", &LCVTXCHI2, &b_LCVTXCHI2);
	fChain->SetBranchAddress("LCVTXNDOF", &LCVTXNDOF, &b_LCVTXNDOF);
	fChain->SetBranchAddress("LCDRJET", &LCDRJET, &b_LCDRJET);
	fChain->SetBranchAddress("LCDRTAG", &LCDRTAG, &b_LCDRTAG);
	fChain->SetBranchAddress("LCPPNN", &LCPPNN, &b_LCPPNN);
	fChain->SetBranchAddress("LCKPNN", &LCKPNN, &b_LCKPNN);
	fChain->SetBranchAddress("LCPIPNN", &LCPIPNN, &b_LCPIPNN);
	fChain->SetBranchAddress("LCP", &LCP, &b_LCP);
	fChain->SetBranchAddress("LCPT", &LCPT, &b_LCPT);
	fChain->SetBranchAddress("LCPPT", &LCPPT, &b_LCPPT);
	fChain->SetBranchAddress("LCKPT", &LCKPT, &b_LCKPT);
	fChain->SetBranchAddress("LCPIPT", &LCPIPT, &b_LCPIPT);
	fChain->SetBranchAddress("LCTRK0", &LCTRK0, &b_LCTRK0);
	fChain->SetBranchAddress("LCTRK1", &LCTRK1, &b_LCTRK1);
	fChain->SetBranchAddress("LCTRK2", &LCTRK2, &b_LCTRK2);
	fChain->SetBranchAddress("LCNJ", &LCNJ, &b_LCNJ);
	fChain->SetBranchAddress("LCMAXDR", &LCMAXDR, &b_LCMAXDR);
	fChain->SetBranchAddress("K3PIM", &K3PIM, &b_K3PIM);
	fChain->SetBranchAddress("K3PIPX", &K3PIPX, &b_K3PIPX);
	fChain->SetBranchAddress("K3PIPY", &K3PIPY, &b_K3PIPY);
	fChain->SetBranchAddress("K3PIPZ", &K3PIPZ, &b_K3PIPZ);
	fChain->SetBranchAddress("K3PIE", &K3PIE, &b_K3PIE);
	fChain->SetBranchAddress("K3PIX", &K3PIX, &b_K3PIX);
	fChain->SetBranchAddress("K3PIY", &K3PIY, &b_K3PIY);
	fChain->SetBranchAddress("K3PIZ", &K3PIZ, &b_K3PIZ);
	fChain->SetBranchAddress("K3PIIP", &K3PIIP, &b_K3PIIP);
	fChain->SetBranchAddress("K3PIIPCHI2", &K3PIIPCHI2, &b_K3PIIPCHI2);
	fChain->SetBranchAddress("K3PIFD", &K3PIFD, &b_K3PIFD);
	fChain->SetBranchAddress("K3PIFDCHI2", &K3PIFDCHI2, &b_K3PIFDCHI2);
	fChain->SetBranchAddress("K3PITAU", &K3PITAU, &b_K3PITAU);
	fChain->SetBranchAddress("K3PIDIRA", &K3PIDIRA, &b_K3PIDIRA);
	fChain->SetBranchAddress("K3PIVTXCHI2", &K3PIVTXCHI2, &b_K3PIVTXCHI2);
	fChain->SetBranchAddress("K3PIVTXNDOF", &K3PIVTXNDOF, &b_K3PIVTXNDOF);
	fChain->SetBranchAddress("K3PIDRJET", &K3PIDRJET, &b_K3PIDRJET);
	fChain->SetBranchAddress("K3PIDRTAG", &K3PIDRTAG, &b_K3PIDRTAG);
	fChain->SetBranchAddress("K3PIPI1PNN", &K3PIPI1PNN, &b_K3PIPI1PNN);
	fChain->SetBranchAddress("K3PIPI2PNN", &K3PIPI2PNN, &b_K3PIPI2PNN);
	fChain->SetBranchAddress("K3PIPI3PNN", &K3PIPI3PNN, &b_K3PIPI3PNN);
	fChain->SetBranchAddress("K3PIKPNN", &K3PIKPNN, &b_K3PIKPNN);
	fChain->SetBranchAddress("K3PIP", &K3PIP, &b_K3PIP);
	fChain->SetBranchAddress("K3PIPT", &K3PIPT, &b_K3PIPT);
	fChain->SetBranchAddress("K3PIPI1PT", &K3PIPI1PT, &b_K3PIPI1PT);
	fChain->SetBranchAddress("K3PIPI2PT", &K3PIPI2PT, &b_K3PIPI2PT);
	fChain->SetBranchAddress("K3PIPI3PT", &K3PIPI3PT, &b_K3PIPI3PT);
	fChain->SetBranchAddress("K3PIKPT", &K3PIKPT, &b_K3PIKPT);
	fChain->SetBranchAddress("K3PITRK0", &K3PITRK0, &b_K3PITRK0);
	fChain->SetBranchAddress("K3PITRK1", &K3PITRK1, &b_K3PITRK1);
	fChain->SetBranchAddress("K3PITRK2", &K3PITRK2, &b_K3PITRK2);
	fChain->SetBranchAddress("K3PITRK3", &K3PITRK3, &b_K3PITRK3);
	fChain->SetBranchAddress("K3PINJ", &K3PINJ, &b_K3PINJ);
	fChain->SetBranchAddress("K3PIMAXDR", &K3PIMAXDR, &b_K3PIMAXDR);
	fChain->SetBranchAddress("JPSIM", &JPSIM, &b_JPSIM);
	fChain->SetBranchAddress("JPSIPX", &JPSIPX, &b_JPSIPX);
	fChain->SetBranchAddress("JPSIPY", &JPSIPY, &b_JPSIPY);
	fChain->SetBranchAddress("JPSIPZ", &JPSIPZ, &b_JPSIPZ);
	fChain->SetBranchAddress("JPSIE", &JPSIE, &b_JPSIE);
	fChain->SetBranchAddress("JPSIX", &JPSIX, &b_JPSIX);
	fChain->SetBranchAddress("JPSIY", &JPSIY, &b_JPSIY);
	fChain->SetBranchAddress("JPSIZ", &JPSIZ, &b_JPSIZ);
	fChain->SetBranchAddress("JPSIIP", &JPSIIP, &b_JPSIIP);
	fChain->SetBranchAddress("JPSIIPCHI2", &JPSIIPCHI2, &b_JPSIIPCHI2);
	fChain->SetBranchAddress("JPSIFD", &JPSIFD, &b_JPSIFD);
	fChain->SetBranchAddress("JPSIFDCHI2", &JPSIFDCHI2, &b_JPSIFDCHI2);
	fChain->SetBranchAddress("JPSITAU", &JPSITAU, &b_JPSITAU);
	fChain->SetBranchAddress("JPSIDIRA", &JPSIDIRA, &b_JPSIDIRA);
	fChain->SetBranchAddress("JPSIVTXCHI2", &JPSIVTXCHI2, &b_JPSIVTXCHI2);
	fChain->SetBranchAddress("JPSIVTXNDOF", &JPSIVTXNDOF, &b_JPSIVTXNDOF);
	fChain->SetBranchAddress("JPSIDRJET", &JPSIDRJET, &b_JPSIDRJET);
	fChain->SetBranchAddress("JPSIDRTAG", &JPSIDRTAG, &b_JPSIDRTAG);
	fChain->SetBranchAddress("JPSIP", &JPSIP, &b_JPSIP);
	fChain->SetBranchAddress("JPSIPT", &JPSIPT, &b_JPSIPT);
	fChain->SetBranchAddress("JPSIPIPT", &JPSIPIPT, &b_JPSIPIPT);
	fChain->SetBranchAddress("JPSIKPT", &JPSIKPT, &b_JPSIKPT);
	fChain->SetBranchAddress("JPSITRK0", &JPSITRK0, &b_JPSITRK0);
	fChain->SetBranchAddress("JPSITRK1", &JPSITRK1, &b_JPSITRK1);
	fChain->SetBranchAddress("JPSINJ", &JPSINJ, &b_JPSINJ);
	fChain->SetBranchAddress("JPSIMAXDR", &JPSIMAXDR, &b_JPSIMAXDR);
	fChain->SetBranchAddress("JPSITRUE", &JPSITRUE, &b_JPSITRUE);
	fChain->SetBranchAddress("JPSITRUEDR", &JPSITRUEDR, &b_JPSITRUEDR);
	fChain->SetBranchAddress("JPSITRUEIDX", &JPSITRUEIDX, &b_JPSITRUEIDX);
	fChain->SetBranchAddress("JPSIFROMB", &JPSIFROMB, &b_JPSIFROMB);
	Notify();
}

Bool_t makeNewD0Effs::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void makeNewD0Effs::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t makeNewD0Effs::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef makeNewD0Effs_cxx
