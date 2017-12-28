#include <iostream>
#include <vector>
#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVectorT.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TSystem.h"

using namespace std;

bool pidCut(double pnn, double px, double py) {
	return true;
	//if(pnn>0.3) return true;
	//double pt2 = px*px + py*py;
	////if(pt2> 9000000. && pnn>0.1) return true;
	//if(pt2>36000000. && pnn>0.1) return true;
	//if(pt2>81000000.) return true;
	//return false;
}

//Calculate distance of closest approach between two lines A and B (each given by a point and a direction vector)
// Use this one! Verified against the second method and documented here:
// http://geomalgorithms.com/a07-_distance.html
double calcDoca(TVector3 &v, const TVector3 &pa, const TVector3 &da, 
		const TVector3 &pb, const TVector3 &db) {

	//lines defined by pa + s*da and pb + t*db
	//vector connecting pa and pb
	TVector3 w0 = pa-pb;

	//define some dot products to simplify the maths
	double a = da.Mag2();
	double b = da.Dot(db);
	double c = db.Mag2();
	double d = da.Dot(w0);
	double e = db.Dot(w0);

	//calculate coefficients of closest point
	double sc(0), tc(0);

	if(a*c - b*b == 0) {
		//lines are parallel - set one coefficient to zero and solve for t'other
		sc = 0;
		tc = e/c;
	} else {
		//general case - see http://geomalgorithms.com/a07-_distance.html
		sc=(b*e - c*d)/(a*c - b*b);
		tc=(a*e - b*d)/(a*c - b*b);
	}

	//points on lines with shortest distance
	TVector3 Pc = pa + sc*da;
	TVector3 Qc = pb + tc*db;

	//give vertex at centre point of the connecting line
	v = Pc+Qc;
	v *= 0.5;

	//return separation at closest point
	return (Pc-Qc).Mag();
}

//Calculate distance of closest approach of line given by point A and direction A to point B
//Equation from Wolfram Alpha
double calcDocaPoint(const TVector3 &pa, const TVector3 &da, const TVector3 &pb) {
	TVector3 x0 = pb;
	TVector3 x1 = pa;
	TVector3 x2 = pa + da;

	return (x0-x1).Cross(x0-x2).Mag()/(x2-x1).Mag();
}

int main(int argc, char *argv[]){

	TString type = argv[1];
	TRandom3 rand;

	char str[1000];

	TFile* fout = TFile::Open("/tmp/dcraik/for_yandex_data_"+type+".root","recreate");
	TTree *tout = new TTree("T","");
	int e(0);
	double JPX,JPY,JPZ,JE,JPT,JETA,JS1,JS2,JQ,JN,JNQ,JNN,JPTD,JDIJETDEC,JDIJETSVDEC, JDIJETSVSVDEC, JDIJETSVMUDEC, JDIJETMUMUDEC;
	double PARTON,BMAXPT,CMAXPT;
	double SVM,SVMCOR,SVMINPERP,SVPT,SVDR,SVN,SVNJ,SVQ,SVFDCHI2,SVPERP,SVETA,SVTZ,
	       SVMINIPCHI2,SVPX,SVPY,SVPZ,SVE,SVMAXGHOST,SVX,SVY,SVZ,SVSUMIPCHI2;
	double NSV,PVX,PVY,PVZ,NDISPL6,NDISPL9,NDISPL16;
	double MUPT,MUIPCHI2,MUDR,MUPNN,NMU;
	double HPT,HIPCHI2,HDR;
	double NTRK, NNEU;
	//double D0M, D0PX, D0PY, D0PZ, D0E, D0X, D0Y, D0Z, D0FD, D0DIRA, D0DOCA, D0IPCHI2MIN, D0DOCAKPI, D0VTXCHI2;
	//double DPMM, DPMPX, DPMPY, DPMPZ, DPME, DPMX, DPMY, DPMZ, DPMFD, DPMDIRA, DPMDOCA, DPMIPCHI2MIN, DPMDOCAMAX;
	//double DSM, DSPX, DSPY, DSPZ, DSE, DSX, DSY, DSZ, DSFD, DSDIRA, DSDOCA, DSIPCHI2MIN, DSDOCAMAX;
	//double LCM, LCPX, LCPY, LCPZ, LCE, LCX, LCY, LCZ, LCFD, LCDIRA, LCDOCA, LCIPCHI2MIN, LCDOCAMAX;
	//double D2K3PIM, D2K3PIPX, D2K3PIPY, D2K3PIPZ, D2K3PIE, D2K3PIX, D2K3PIY, D2K3PIZ, D2K3PIFD, D2K3PIDIRA, D2K3PIDOCA, D2K3PIIPCHI2MIN, D2K3PIDOCAMAX;
	double     D0M,     D0PX,     D0PY,     D0PZ,     D0E,     D0X,     D0Y,     D0Z,     D0FD,     D0FDCHI2,     D0IP,     D0IPCHI2,     D0VTXCHI2,     D0VTXNDOF,     D0TAU;
	double    DPMM,    DPMPX,    DPMPY,    DPMPZ,    DPME,    DPMX,    DPMY,    DPMZ,    DPMFD,    DPMFDCHI2,    DPMIP,    DPMIPCHI2,    DPMVTXCHI2,    DPMVTXNDOF,    DPMTAU;
	double     DSM,     DSPX,     DSPY,     DSPZ,     DSE,     DSX,     DSY,     DSZ,     DSFD,     DSFDCHI2,     DSIP,     DSIPCHI2,     DSVTXCHI2,     DSVTXNDOF,     DSTAU;
	double     LCM,     LCPX,     LCPY,     LCPZ,     LCE,     LCX,     LCY,     LCZ,     LCFD,     LCFDCHI2,     LCIP,     LCIPCHI2,     LCVTXCHI2,     LCVTXNDOF,     LCTAU;
	double D2K3PIM, D2K3PIPX, D2K3PIPY, D2K3PIPZ, D2K3PIE, D2K3PIX, D2K3PIY, D2K3PIZ, D2K3PIFD, D2K3PIFDCHI2, D2K3PIIP, D2K3PIIPCHI2, D2K3PIVTXCHI2, D2K3PIVTXNDOF, D2K3PITAU;
	double D0KPNN, D0PIPNN;
	double DPMKPNN, DPMPI1PNN, DPMPI2PNN;
	double DSK1PNN, DSK2PNN, DSPIPNN, DSPHIM;
	double LCPPNN, LCKPNN, LCPIPNN;
	double D2K3PIKPNN, D2K3PIPI1PNN, D2K3PIPI2PNN, D2K3PIPI3PNN;
	//double D0Mb, DPMMb, DSMb;

	tout->Branch("EVT",&e);
	tout->Branch("TrueParton",&PARTON);
	tout->Branch("TrueMaxBPT",&BMAXPT);
	tout->Branch("TrueMaxCPT",&CMAXPT);
	tout->Branch("TrueJetPx",&JPX);
	tout->Branch("TrueJetPy",&JPY);
	tout->Branch("TrueJetPz",&JPZ);
	tout->Branch("TrueJetE",&JE);
	tout->Branch("JetPT",&JPT);
	tout->Branch("JetEta",&JETA);
	tout->Branch("JetSigma1",&JS1);
	tout->Branch("JetSigma2",&JS2);
	tout->Branch("JetQ",&JQ);
	tout->Branch("JetMult",&JN);
	tout->Branch("JetNChr",&JNQ);
	tout->Branch("JetNNeu",&JNN);
	tout->Branch("JetPTD",&JPTD);
	tout->Branch("JetDijetDec",&JDIJETDEC);
	tout->Branch("JetDijetSVDec",&JDIJETSVDEC);
	tout->Branch("JetDijetSVSVDec",&JDIJETSVSVDEC);
	tout->Branch("JetDijetSVMuDec",&JDIJETSVMUDEC);
	tout->Branch("JetDijetMuMuDec",&JDIJETMUMUDEC);
	tout->Branch("SVX",&SVX);
	tout->Branch("SVY",&SVY);
	tout->Branch("SVZ",&SVZ);
	tout->Branch("SVPerp",&SVPERP);
	tout->Branch("SVPx",&SVPX);
	tout->Branch("SVPy",&SVPY);
	tout->Branch("SVPz",&SVPZ);
	tout->Branch("SVE",&SVE);
	tout->Branch("SVPT",&SVPT);
	tout->Branch("SVETA",&SVETA);
	tout->Branch("SVM",&SVM);
	tout->Branch("SVMCor",&SVMCOR);
	tout->Branch("SVMINPERP",&SVMINPERP);
	tout->Branch("SVDR",&SVDR);
	tout->Branch("SVN",&SVN);
	tout->Branch("SVNJ",&SVNJ);
	tout->Branch("SVQ",&SVQ);
	//tout->Branch("SVFDChi2",&SVFDCHI2);
	tout->Branch("SVSumIPChi2",&SVSUMIPCHI2);
	tout->Branch("SVTZ",&SVTZ);
	tout->Branch("SVMINIPCHI2",&SVMINIPCHI2);
	tout->Branch("SVGhostMax",&SVMAXGHOST);
	tout->Branch("NSV",&NSV);
	tout->Branch("NTRK",&NTRK);
	tout->Branch("NNEU",&NNEU);
	tout->Branch("PVX",&PVX);
	tout->Branch("PVY",&PVY);
	tout->Branch("PVZ",&PVZ);
	tout->Branch("NDispl6",&NDISPL6);
	tout->Branch("NDispl9",&NDISPL9);
	tout->Branch("NDispl16",&NDISPL16);
	tout->Branch("MuPT",&MUPT);
	tout->Branch("MuIPChi2",&MUIPCHI2);
	tout->Branch("MuDR",&MUDR);
	tout->Branch("MuPNN",&MUPNN);
	tout->Branch("NMu",&NMU);
	tout->Branch("HardPT",&HPT);
	tout->Branch("HardIPChi2",&HIPCHI2);
	tout->Branch("HardDR",&HDR);
	tout->Branch("D0M",        &D0M);
	tout->Branch("D0PX",       &D0PX);
	tout->Branch("D0PY",       &D0PY);
	tout->Branch("D0PZ",       &D0PZ);
	tout->Branch("D0E",        &D0E);
	tout->Branch("D0X",        &D0X);
	tout->Branch("D0Y",        &D0Y);
	tout->Branch("D0Z",        &D0Z);
	tout->Branch("D0IP",       &D0IP);
	tout->Branch("D0IPCHI2",   &D0IPCHI2);
	tout->Branch("D0FD",       &D0FD);
	tout->Branch("D0FDCHI2",   &D0FDCHI2);
	tout->Branch("D0TAU",      &D0TAU);
	tout->Branch("D0VTXCHI2",  &D0VTXCHI2);
	tout->Branch("D0VTXNDOF",  &D0VTXNDOF);
	tout->Branch("D0PIPNN",    &D0PIPNN);
	tout->Branch("D0KPNN",     &D0KPNN);
	//tout->Branch("D0Mb",        &D0Mb);
	//tout->Branch("D0PX",       &D0PX);
	//tout->Branch("D0PY",       &D0PY);
	//tout->Branch("D0PZ",       &D0PZ);
	//tout->Branch("D0E",        &D0E);
	//tout->Branch("D0X",        &D0X);
	//tout->Branch("D0Y",        &D0Y);
	//tout->Branch("D0Z",        &D0Z);
	//tout->Branch("D0FD",       &D0FD);
	//tout->Branch("D0DIRA",     &D0DIRA);
	//tout->Branch("D0IP",     &D0DOCA);
	//tout->Branch("D0IPCHI2MIN",    &D0IPCHI2MIN);
	//tout->Branch("D0DOCAKPI",  &D0DOCAKPI);
	//tout->Branch("D0VTXCHI2",  &D0VTXCHI2);
	tout->Branch("DM",        &DPMM);
	tout->Branch("DPX",       &DPMPX);
	tout->Branch("DPY",       &DPMPY);
	tout->Branch("DPZ",       &DPMPZ);
	tout->Branch("DE",        &DPME);
	tout->Branch("DX",        &DPMX);
	tout->Branch("DY",        &DPMY);
	tout->Branch("DZ",        &DPMZ);
	tout->Branch("DIP",       &DPMIP);
	tout->Branch("DIPCHI2",   &DPMIPCHI2);
	tout->Branch("DFD",       &DPMFD);
	tout->Branch("DFDCHI2",   &DPMFDCHI2);
	tout->Branch("DTAU",      &DPMTAU);
	tout->Branch("DVTXCHI2",  &DPMVTXCHI2);
	tout->Branch("DVTXNDOF",  &DPMVTXNDOF);
	tout->Branch("DPI1PNN",   &DPMPI1PNN);
	tout->Branch("DPI2PNN",   &DPMPI2PNN);
	tout->Branch("DKPNN",     &DPMKPNN);
	//tout->Branch("DMb",        &DPMMb);
	//tout->Branch("DPX",       &DPMPX);
	//tout->Branch("DPY",       &DPMPY);
	//tout->Branch("DPZ",       &DPMPZ);
	//tout->Branch("DE",        &DPME);
	//tout->Branch("DX",        &DPMX);
	//tout->Branch("DY",        &DPMY);
	//tout->Branch("DZ",        &DPMZ);
	//tout->Branch("DFD",       &DPMFD);
	//tout->Branch("DDIRA",     &DPMDIRA);
	//tout->Branch("DIP",     &DPMDOCA);
	//tout->Branch("DIPCHI2MIN",    &DPMIPCHI2MIN);
	//tout->Branch("DDOCAMAX",  &DPMDOCAMAX);
	tout->Branch("DSM",        &DSM);
	tout->Branch("DSPX",       &DSPX);
	tout->Branch("DSPY",       &DSPY);
	tout->Branch("DSPZ",       &DSPZ);
	tout->Branch("DSE",        &DSE);
	tout->Branch("DSX",        &DSX);
	tout->Branch("DSY",        &DSY);
	tout->Branch("DSZ",        &DSZ);
	tout->Branch("DSIP",       &DSIP);
	tout->Branch("DSIPCHI2",   &DSIPCHI2);
	tout->Branch("DSFD",       &DSFD);
	tout->Branch("DSFDCHI2",   &DSFDCHI2);
	tout->Branch("DSTAU",      &DSTAU);
	tout->Branch("DSVTXCHI2",  &DSVTXCHI2);
	tout->Branch("DSVTXNDOF",  &DSVTXNDOF);
	tout->Branch("DSPIPNN",    &DSPIPNN);
	tout->Branch("DSK1PNN",    &DSK1PNN);
	tout->Branch("DSK2PNN",    &DSK2PNN);
	tout->Branch("DSPHIM",     &DSPHIM);
	//tout->Branch("DSMb",        &DSMb);
	//tout->Branch("DSPX",       &DSPX);
	//tout->Branch("DSPY",       &DSPY);
	//tout->Branch("DSPZ",       &DSPZ);
	//tout->Branch("DSE",        &DSE);
	//tout->Branch("DSX",        &DSX);
	//tout->Branch("DSY",        &DSY);
	//tout->Branch("DSZ",        &DSZ);
	//tout->Branch("DSFD",       &DSFD);
	//tout->Branch("DSDIRA",     &DSDIRA);
	//tout->Branch("DSIP",     &DSDOCA);
	//tout->Branch("DSIPCHI2MIN",    &DSIPCHI2MIN);
	//tout->Branch("DSDOCAMAX",  &DSDOCAMAX);
	tout->Branch("LCM",        &LCM);
	tout->Branch("LCPX",       &LCPX);
	tout->Branch("LCPY",       &LCPY);
	tout->Branch("LCPZ",       &LCPZ);
	tout->Branch("LCE",        &LCE);
	tout->Branch("LCX",        &LCX);
	tout->Branch("LCY",        &LCY);
	tout->Branch("LCZ",        &LCZ);
	tout->Branch("LCIP",       &LCIP);
	tout->Branch("LCIPCHI2",   &LCIPCHI2);
	tout->Branch("LCFD",       &LCFD);
	tout->Branch("LCFDCHI2",   &LCFDCHI2);
	tout->Branch("LCTAU",      &LCTAU);
	tout->Branch("LCVTXCHI2",  &LCVTXCHI2);
	tout->Branch("LCVTXNDOF",  &LCVTXNDOF);
	tout->Branch("LCPIPNN",    &LCPIPNN);
	tout->Branch("LCKPNN",    &LCKPNN);
	tout->Branch("LCPPNN",    &LCPPNN);
	//tout->Branch("LCM",        &LCM);
	//tout->Branch("LCPX",       &LCPX);
	//tout->Branch("LCPY",       &LCPY);
	//tout->Branch("LCPZ",       &LCPZ);
	//tout->Branch("LCE",        &LCE);
	//tout->Branch("LCX",        &LCX);
	//tout->Branch("LCY",        &LCY);
	//tout->Branch("LCZ",        &LCZ);
	//tout->Branch("LCFD",       &LCFD);
	//tout->Branch("LCDIRA",     &LCDIRA);
	//tout->Branch("LCIP",     &LCDOCA);
	//tout->Branch("LCIPCHI2MIN",    &LCIPCHI2MIN);
	//tout->Branch("LCDOCAMAX",  &LCDOCAMAX);
	tout->Branch("D2K3PIM",   &D2K3PIM);
	tout->Branch("D2K3PIPX",       &D2K3PIPX);
	tout->Branch("D2K3PIPY",       &D2K3PIPY);
	tout->Branch("D2K3PIPZ",       &D2K3PIPZ);
	tout->Branch("D2K3PIE",        &D2K3PIE);
	tout->Branch("D2K3PIX",        &D2K3PIX);
	tout->Branch("D2K3PIY",        &D2K3PIY);
	tout->Branch("D2K3PIZ",        &D2K3PIZ);
	tout->Branch("D2K3PIIP",       &D2K3PIIP);
	tout->Branch("D2K3PIIPCHI2",   &D2K3PIIPCHI2);
	tout->Branch("D2K3PIFD",       &D2K3PIFD);
	tout->Branch("D2K3PIFDCHI2",   &D2K3PIFDCHI2);
	tout->Branch("D2K3PITAU",      &D2K3PITAU);
	tout->Branch("D2K3PIVTXCHI2",  &D2K3PIVTXCHI2);
	tout->Branch("D2K3PIVTXNDOF",  &D2K3PIVTXNDOF);
	tout->Branch("D2K3PIPI1PNN",   &D2K3PIPI1PNN);
	tout->Branch("D2K3PIPI2PNN",   &D2K3PIPI2PNN);
	tout->Branch("D2K3PIPI3PNN",   &D2K3PIPI3PNN);
	tout->Branch("D2K3PIKPNN",     &D2K3PIKPNN);
	//tout->Branch("D2K3PIPX",       &D2K3PIPX);
	//tout->Branch("D2K3PIPY",       &D2K3PIPY);
	//tout->Branch("D2K3PIPZ",       &D2K3PIPZ);
	//tout->Branch("D2K3PIE",        &D2K3PIE);
	//tout->Branch("D2K3PIX",        &D2K3PIX);
	//tout->Branch("D2K3PIY",        &D2K3PIY);
	//tout->Branch("D2K3PIZ",        &D2K3PIZ);
	//tout->Branch("D2K3PIFD",       &D2K3PIFD);
	//tout->Branch("D2K3PIDIRA",     &D2K3PIDIRA);
	//tout->Branch("D2K3PIIP",     &D2K3PIDOCA);
	//tout->Branch("D2K3PIIPCHI2MIN",    &D2K3PIIPCHI2MIN);
	//tout->Branch("D2K3PIDOCAMAX",  &D2K3PIDOCAMAX);

	TChain *t = new TChain("data");
//	t->Add("/afs/cern.ch/user/d/dcraik/workspace/gangadir-jets/workspace/dcraik/LocalXML/100/0/output/output.root");
	//boost::progress_display show_addfile_progress( 270 );
	//for(int i=0; i<270; ++i) {
	//      ++show_addfile_progress;
	//	sprintf(str,"/eos/lhcb/user/d/dcraik/jets/158/%d/output.root",i);
	//	TFile* f = TFile::Open(str);
	//	if(f) {
	//	      delete f;
	//	      t->Add(str);
	//	}
	//}
	//boost::progress_display show_addfile_progress2( 953 );
	//for(int i=0; i<953; ++i) {
	//      ++show_addfile_progress2;
	//	sprintf(str,"/eos/lhcb/user/d/dcraik/jets/193/%d/output.root",i);
	//	TFile* f = TFile::Open(str);
	//	if(f) {
	//	      delete f;
	//	      t->Add(str);
	//	}
	//}
	//boost::progress_display show_addfile_progress( 949 );
	//for(int i=0; i<949; ++i) {
	//      ++show_addfile_progress;
	//	sprintf(str,"/eos/lhcb/user/d/dcraik/jets/159/%d/output.root",i);
	//	TFile* f = TFile::Open(str);
	//	if(f) {
	//	      delete f;
	//	      t->Add(str);
	//	}
	//}
	//boost::progress_display show_addfile_progress( 1049 );
	//for(int i=0; i<1049; ++i) {
	//	++show_addfile_progress;
	//	sprintf(str,"/eos/lhcb/user/d/dcraik/jets/227/%d/output.root",i);
	//	if(gSystem->AccessPathName(str)) continue;
	//	t->Add(str);
	//}
	//job 289: found 948 of 1050
	//job 290: found 784 of 950
	if(type=="testA") {
		boost::progress_display show_addfile_progress( 1 );
		for(int i=0; i<1; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/289/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="testA1k") {
		boost::progress_display show_addfile_progress( 1050 );
		for(int i=0; i<1050; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/289/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="testB") {
		boost::progress_display show_addfile_progress( 1 );
		for(int i=0; i<1; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/361/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="testB1k") {
		boost::progress_display show_addfile_progress( 1050 );
		for(int i=0; i<1050; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/361/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="testE1k") {
		boost::progress_display show_addfile_progress( 1050 );
		for(int i=0; i<1050; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/379/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="MD") {
		boost::progress_display show_addfile_progress( 1050 );
		for(int i=0; i<1050; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/289/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}
	if(type=="MU") {
		boost::progress_display show_addfile_progress( 950 );
		for(int i=0; i<950; ++i) {
			++show_addfile_progress;
			sprintf(str,"/eos/lhcb/user/d/dcraik/jets/290/%d/output.root",i);
			if(gSystem->AccessPathName(str)) continue;
			t->Add(str);
		}
	}

	//  vector<double> *gen_px = new vector<double>();
	//  vector<double> *gen_py = new vector<double>();
	//  vector<double> *gen_pz = new vector<double>();
	//  vector<double> *gen_e = new vector<double>();
	//  vector<double> *gen_pid = new vector<double>();
	//  t->SetBranchAddress("gen_px",&gen_px);  
	//  t->SetBranchAddress("gen_py",&gen_py);  
	//  t->SetBranchAddress("gen_pz",&gen_pz);  
	//  t->SetBranchAddress("gen_e",&gen_e);  
	//  t->SetBranchAddress("gen_pid",&gen_pid);  

	vector<double> *jet_px = new vector<double>();
	vector<double> *jet_py = new vector<double>();
	vector<double> *jet_pz = new vector<double>();
	vector<double> *jet_e = new vector<double>();
	vector<double> *jet_pv = new vector<double>();

	vector<double> *evt_dec = new vector<double>();

	vector<double> *evt_j1_idx = new vector<double>();
	vector<double> *evt_j2_idx = new vector<double>();

	vector<double> *evt_j1_nsv = new vector<double>();
	vector<double> *evt_j2_nsv = new vector<double>();

	vector<double> *evt_j1_nmu = new vector<double>();
	vector<double> *evt_j2_nmu = new vector<double>();

	vector<double> *evt_j1_px = new vector<double>();
	vector<double> *evt_j1_py = new vector<double>();
	vector<double> *evt_j1_pz = new vector<double>();
	vector<double> *evt_j2_px = new vector<double>();
	vector<double> *evt_j2_py = new vector<double>();
	vector<double> *evt_j2_pz = new vector<double>();

	vector<double> *d0_mass     = new vector<double>();
	vector<double> *d_mass      = new vector<double>();
	vector<double> *ds_mass     = new vector<double>();
	vector<double> *lc_mass     = new vector<double>();
	vector<double> *d2k3pi_mass = new vector<double>();
	vector<double> *d0_j = new vector<double>();
	vector<double> *d_j = new vector<double>();
	vector<double> *ds_j = new vector<double>();
	vector<double> *lc_j = new vector<double>();
	vector<double> *d2k3pi_j = new vector<double>();
	vector<double> *d0_nj = new vector<double>();
	vector<double> *d_nj = new vector<double>();
	vector<double> *ds_nj = new vector<double>();
	vector<double> *lc_nj = new vector<double>();
	vector<double> *d2k3pi_nj = new vector<double>();
	vector<double> *d0_trk0 = new vector<double>();
	vector<double> *d_trk0 = new vector<double>();
	vector<double> *ds_trk0 = new vector<double>();
	vector<double> *lc_trk0 = new vector<double>();
	vector<double> *d2k3pi_trk0 = new vector<double>();
	vector<double> *d0_trk1 = new vector<double>();
	vector<double> *d_trk1 = new vector<double>();
	vector<double> *ds_trk1 = new vector<double>();
	vector<double> *lc_trk1 = new vector<double>();
	vector<double> *d2k3pi_trk1 = new vector<double>();
	vector<double> *d_trk2 = new vector<double>();
	vector<double> *ds_trk2 = new vector<double>();
	vector<double> *lc_trk2 = new vector<double>();
	vector<double> *d2k3pi_trk2 = new vector<double>();
	vector<double> *d2k3pi_trk3 = new vector<double>();
	vector<double> *d0_px     = new vector<double>();
	vector<double> *d_px      = new vector<double>();
	vector<double> *ds_px     = new vector<double>();
	vector<double> *lc_px     = new vector<double>();
	vector<double> *d2k3pi_px = new vector<double>();
	vector<double> *d0_py     = new vector<double>();
	vector<double> *d_py      = new vector<double>();
	vector<double> *ds_py     = new vector<double>();
	vector<double> *lc_py     = new vector<double>();
	vector<double> *d2k3pi_py = new vector<double>();
	vector<double> *d0_pz     = new vector<double>();
	vector<double> *d_pz      = new vector<double>();
	vector<double> *ds_pz     = new vector<double>();
	vector<double> *lc_pz     = new vector<double>();
	vector<double> *d2k3pi_pz = new vector<double>();
	vector<double> *d0_e     = new vector<double>();
	vector<double> *d_e      = new vector<double>();
	vector<double> *ds_e     = new vector<double>();
	vector<double> *lc_e     = new vector<double>();
	vector<double> *d2k3pi_e = new vector<double>();
	vector<double> *d0_x     = new vector<double>();
	vector<double> *d_x      = new vector<double>();
	vector<double> *ds_x     = new vector<double>();
	vector<double> *lc_x     = new vector<double>();
	vector<double> *d2k3pi_x = new vector<double>();
	vector<double> *d0_y     = new vector<double>();
	vector<double> *d_y      = new vector<double>();
	vector<double> *ds_y     = new vector<double>();
	vector<double> *lc_y     = new vector<double>();
	vector<double> *d2k3pi_y = new vector<double>();
	vector<double> *d0_z     = new vector<double>();
	vector<double> *d_z      = new vector<double>();
	vector<double> *ds_z     = new vector<double>();
	vector<double> *lc_z     = new vector<double>();
	vector<double> *d2k3pi_z = new vector<double>();
	vector<double> *d0_ip     = new vector<double>();
	vector<double> *d_ip      = new vector<double>();
	vector<double> *ds_ip     = new vector<double>();
	vector<double> *lc_ip     = new vector<double>();
	vector<double> *d2k3pi_ip = new vector<double>();
	vector<double> *d0_ipchi2     = new vector<double>();
	vector<double> *d_ipchi2      = new vector<double>();
	vector<double> *ds_ipchi2     = new vector<double>();
	vector<double> *lc_ipchi2     = new vector<double>();
	vector<double> *d2k3pi_ipchi2 = new vector<double>();
	vector<double> *d0_fd     = new vector<double>();
	vector<double> *d_fd      = new vector<double>();
	vector<double> *ds_fd     = new vector<double>();
	vector<double> *lc_fd     = new vector<double>();
	vector<double> *d2k3pi_fd = new vector<double>();
	vector<double> *d0_fdchi2     = new vector<double>();
	vector<double> *d_fdchi2      = new vector<double>();
	vector<double> *ds_fdchi2     = new vector<double>();
	vector<double> *lc_fdchi2     = new vector<double>();
	vector<double> *d2k3pi_fdchi2 = new vector<double>();
	vector<double> *d0_tau     = new vector<double>();
	vector<double> *d_tau      = new vector<double>();
	vector<double> *ds_tau     = new vector<double>();
	vector<double> *lc_tau     = new vector<double>();
	vector<double> *d2k3pi_tau = new vector<double>();
	vector<double> *d0_vtxchi2     = new vector<double>();
	vector<double> *d_vtxchi2      = new vector<double>();
	vector<double> *ds_vtxchi2     = new vector<double>();
	vector<double> *lc_vtxchi2     = new vector<double>();
	vector<double> *d2k3pi_vtxchi2 = new vector<double>();
	vector<double> *d0_vtxndof     = new vector<double>();
	vector<double> *d_vtxndof      = new vector<double>();
	vector<double> *ds_vtxndof     = new vector<double>();
	vector<double> *lc_vtxndof     = new vector<double>();
	vector<double> *d2k3pi_vtxndof = new vector<double>();

	t->SetBranchAddress("jet_px",&jet_px);  
	t->SetBranchAddress("jet_py",&jet_py);  
	t->SetBranchAddress("jet_pz",&jet_pz);  
	t->SetBranchAddress("jet_e",&jet_e);  
	t->SetBranchAddress("jet_idx_pvr",&jet_pv);  

	t->SetBranchAddress("evt_dec",&evt_dec);  

	t->SetBranchAddress("evt_j1_idx",&evt_j1_idx);  
	t->SetBranchAddress("evt_j2_idx",&evt_j2_idx);  

	t->SetBranchAddress("evt_j1_nsv",&evt_j1_nsv);  
	t->SetBranchAddress("evt_j2_nsv",&evt_j2_nsv);  

	t->SetBranchAddress("evt_j1_nmu",&evt_j1_nmu);  
	t->SetBranchAddress("evt_j2_nmu",&evt_j2_nmu);  

	t->SetBranchAddress("evt_j1_px",&evt_j1_px);  
	t->SetBranchAddress("evt_j1_py",&evt_j1_py);  
	t->SetBranchAddress("evt_j1_pz",&evt_j1_pz);  

	t->SetBranchAddress("evt_j2_px",&evt_j2_px);  
	t->SetBranchAddress("evt_j2_py",&evt_j2_py);  
	t->SetBranchAddress("evt_j2_pz",&evt_j2_pz);  

	t->SetBranchAddress("d0_m",&d0_mass);  
	t->SetBranchAddress("dp_m",&d_mass);  
	t->SetBranchAddress("ds_m",&ds_mass);  
	t->SetBranchAddress("lc_m",&lc_mass);  
	t->SetBranchAddress("k3pi_m",&d2k3pi_mass);  

	t->SetBranchAddress("d0_idx_jet",&d0_j);  
	t->SetBranchAddress("dp_idx_jet",&d_j);  
	t->SetBranchAddress("ds_idx_jet",&ds_j);  
	t->SetBranchAddress("lc_idx_jet",&lc_j);  
	t->SetBranchAddress("k3pi_idx_jet",&d2k3pi_j);  

	t->SetBranchAddress("d0_ntrk_jet",&d0_nj);  
	t->SetBranchAddress("dp_ntrk_jet",&d_nj);  
	t->SetBranchAddress("ds_ntrk_jet",&ds_nj);  
	t->SetBranchAddress("lc_ntrk_jet",&lc_nj);  
	t->SetBranchAddress("k3pi_ntrk_jet",&d2k3pi_nj);  

	t->SetBranchAddress("d0_idx_trk0",&d0_trk0);  
	t->SetBranchAddress("dp_idx_trk0",&d_trk0);  
	t->SetBranchAddress("ds_idx_trk0",&ds_trk0);  
	t->SetBranchAddress("lc_idx_trk0",&lc_trk0);  
	t->SetBranchAddress("k3pi_idx_trk0",&d2k3pi_trk0);  

	t->SetBranchAddress("d0_idx_trk1",&d0_trk1);  
	t->SetBranchAddress("dp_idx_trk1",&d_trk1);  
	t->SetBranchAddress("ds_idx_trk1",&ds_trk1);  
	t->SetBranchAddress("lc_idx_trk1",&lc_trk1);  
	t->SetBranchAddress("k3pi_idx_trk1",&d2k3pi_trk1);  

	t->SetBranchAddress("dp_idx_trk2",&d_trk2);  
	t->SetBranchAddress("ds_idx_trk2",&ds_trk2);  
	t->SetBranchAddress("lc_idx_trk2",&lc_trk2);  
	t->SetBranchAddress("k3pi_idx_trk2",&d2k3pi_trk2);  

	t->SetBranchAddress("k3pi_idx_trk3",&d2k3pi_trk3);  

	t->SetBranchAddress(     "d0_px",     &d0_px);  
	t->SetBranchAddress(     "dp_px",      &d_px);  
	t->SetBranchAddress(     "ds_px",     &ds_px);  
	t->SetBranchAddress(     "lc_px",     &lc_px);  
	t->SetBranchAddress("k3pi_px", &d2k3pi_px);  

	t->SetBranchAddress(     "d0_py",     &d0_py);  
	t->SetBranchAddress(     "dp_py",      &d_py);  
	t->SetBranchAddress(     "ds_py",     &ds_py);  
	t->SetBranchAddress(     "lc_py",     &lc_py);  
	t->SetBranchAddress("k3pi_py", &d2k3pi_py);  

	t->SetBranchAddress(     "d0_pz",     &d0_pz);  
	t->SetBranchAddress(     "dp_pz",      &d_pz);  
	t->SetBranchAddress(     "ds_pz",     &ds_pz);  
	t->SetBranchAddress(     "lc_pz",     &lc_pz);  
	t->SetBranchAddress("k3pi_pz", &d2k3pi_pz);  

	t->SetBranchAddress(     "d0_e",     &d0_e);  
	t->SetBranchAddress(     "dp_e",      &d_e);  
	t->SetBranchAddress(     "ds_e",     &ds_e);  
	t->SetBranchAddress(     "lc_e",     &lc_e);  
	t->SetBranchAddress("k3pi_e", &d2k3pi_e);  

	t->SetBranchAddress(     "d0_x",     &d0_x);  
	t->SetBranchAddress(     "dp_x",      &d_x);  
	t->SetBranchAddress(     "ds_x",     &ds_x);  
	t->SetBranchAddress(     "lc_x",     &lc_x);  
	t->SetBranchAddress("k3pi_x", &d2k3pi_x);  

	t->SetBranchAddress(     "d0_y",     &d0_y);  
	t->SetBranchAddress(     "dp_y",      &d_y);  
	t->SetBranchAddress(     "ds_y",     &ds_y);  
	t->SetBranchAddress(     "lc_y",     &lc_y);  
	t->SetBranchAddress("k3pi_y", &d2k3pi_y);  

	t->SetBranchAddress(     "d0_z",     &d0_z);  
	t->SetBranchAddress(     "dp_z",      &d_z);  
	t->SetBranchAddress(     "ds_z",     &ds_z);  
	t->SetBranchAddress(     "lc_z",     &lc_z);  
	t->SetBranchAddress("k3pi_z", &d2k3pi_z);  

	t->SetBranchAddress(     "d0_ip",     &d0_ip);  
	t->SetBranchAddress(     "dp_ip",      &d_ip);  
	t->SetBranchAddress(     "ds_ip",     &ds_ip);  
	t->SetBranchAddress(     "lc_ip",     &lc_ip);  
	t->SetBranchAddress("k3pi_ip", &d2k3pi_ip);  

	t->SetBranchAddress(     "d0_ip_chi2",     &d0_ipchi2);  
	t->SetBranchAddress(     "dp_ip_chi2",      &d_ipchi2);  
	t->SetBranchAddress(     "ds_ip_chi2",     &ds_ipchi2);  
	t->SetBranchAddress(     "lc_ip_chi2",     &lc_ipchi2);  
	t->SetBranchAddress("k3pi_ip_chi2", &d2k3pi_ipchi2);  

	t->SetBranchAddress(     "d0_fd",     &d0_fd);  
	t->SetBranchAddress(     "dp_fd",      &d_fd);  
	t->SetBranchAddress(     "ds_fd",     &ds_fd);  
	t->SetBranchAddress(     "lc_fd",     &lc_fd);  
	t->SetBranchAddress("k3pi_fd", &d2k3pi_fd);  

	t->SetBranchAddress(     "d0_fd_chi2",     &d0_fdchi2);  
	t->SetBranchAddress(     "dp_fd_chi2",      &d_fdchi2);  
	t->SetBranchAddress(     "ds_fd_chi2",     &ds_fdchi2);  
	t->SetBranchAddress(     "lc_fd_chi2",     &lc_fdchi2);  
	t->SetBranchAddress("k3pi_fd_chi2", &d2k3pi_fdchi2);  

	t->SetBranchAddress(     "d0_tau",     &d0_tau);  
	t->SetBranchAddress(     "dp_tau",      &d_tau);  
	t->SetBranchAddress(     "ds_tau",     &ds_tau);  
	t->SetBranchAddress(     "lc_tau",     &lc_tau);  
	t->SetBranchAddress("k3pi_tau", &d2k3pi_tau);  

	t->SetBranchAddress(     "d0_vtx_chi2",     &d0_vtxchi2);  
	t->SetBranchAddress(     "dp_vtx_chi2",      &d_vtxchi2);  
	t->SetBranchAddress(     "ds_vtx_chi2",     &ds_vtxchi2);  
	t->SetBranchAddress(     "lc_vtx_chi2",     &lc_vtxchi2);  
	t->SetBranchAddress("k3pi_vtx_chi2", &d2k3pi_vtxchi2);  

	t->SetBranchAddress(     "d0_vtx_ndof",     &d0_vtxndof);  
	t->SetBranchAddress(     "dp_vtx_ndof",      &d_vtxndof);  
	t->SetBranchAddress(     "ds_vtx_ndof",     &ds_vtxndof);  
	t->SetBranchAddress(     "lc_vtx_ndof",     &lc_vtxndof);  
	t->SetBranchAddress("k3pi_vtx_ndof", &d2k3pi_vtxndof);  

	double npv;
	t->SetBranchAddress("evt_pvr_n",&npv); 

	vector<double> *pvx = new vector<double>();
	vector<double> *pvy = new vector<double>();
	vector<double> *pvz = new vector<double>();
	t->SetBranchAddress("pvr_x",&pvx);  
	t->SetBranchAddress("pvr_y",&pvy);  
	t->SetBranchAddress("pvr_z",&pvz);  

	vector<double> *svx = new vector<double>();
	vector<double> *svy = new vector<double>();
	vector<double> *svz = new vector<double>();
	vector<double> *svpx = new vector<double>();
	vector<double> *svpy = new vector<double>();
	vector<double> *svpz = new vector<double>();  
	vector<double> *sve = new vector<double>();  
	vector<double> *svfdchi2 = new vector<double>();  
	vector<double> *svpv = new vector<double>();  
	vector<double> *svmcor = new vector<double>();  
	vector<double> *svfdmin = new vector<double>();  
	vector<double> *svtrk[10];
	t->SetBranchAddress("svr_x",&svx);  
	t->SetBranchAddress("svr_y",&svy);  
	t->SetBranchAddress("svr_z",&svz);  
	t->SetBranchAddress("svr_px",&svpx);  
	t->SetBranchAddress("svr_py",&svpy);  
	t->SetBranchAddress("svr_pz",&svpz);  
	t->SetBranchAddress("svr_e",&sve);  
	//  t->SetBranchAddress("svr_fd_chi2",&svfdchi2);  
	t->SetBranchAddress("svr_idx_pvr",&svpv);  
	t->SetBranchAddress("svr_m_cor",&svmcor);  
	t->SetBranchAddress("svr_fd_min",&svfdmin);  
	for(int i=0; i<10; i++){
		svtrk[i] = new vector<double>();
		sprintf(str,"svr_idx_trk%d",i);
		t->SetBranchAddress(str,&svtrk[i]);
	}

	vector<double> *trk_q = new vector<double>();
	vector<double> *trk_ghost = new vector<double>();
	vector<double> *trk_x = new vector<double>();
	vector<double> *trk_y = new vector<double>();
	vector<double> *trk_z = new vector<double>();
	vector<double> *trk_px = new vector<double>();
	vector<double> *trk_py = new vector<double>();
	vector<double> *trk_pz = new vector<double>();
	vector<double> *trk_e = new vector<double>();
	vector<double> *trk_ip = new vector<double>();
	vector<double> *trk_ipchi2 = new vector<double>();
	vector<double> *trk_pnnk = new vector<double>();
	vector<double> *trk_pnnpi = new vector<double>();
	vector<double> *trk_pnnp = new vector<double>();
	vector<double> *trk_pnnmu = new vector<double>();
	vector<double> *trk_ismu = new vector<double>();
	vector<double> *trk_j = new vector<double>();
	vector<double>  *trk_pid = new vector<double>();
	vector<double>  *trk_type = new vector<double>();
	t->SetBranchAddress("trk_q",&trk_q);  
	t->SetBranchAddress("trk_prb_ghost",&trk_ghost);  
	t->SetBranchAddress("trk_x",&trk_x);  
	t->SetBranchAddress("trk_y",&trk_y);  
	t->SetBranchAddress("trk_z",&trk_z);  
	t->SetBranchAddress("trk_px",&trk_px);  
	t->SetBranchAddress("trk_py",&trk_py);  
	t->SetBranchAddress("trk_pz",&trk_pz);  
	t->SetBranchAddress("trk_e",&trk_e);    
	t->SetBranchAddress("trk_ip",&trk_ip);    
	t->SetBranchAddress("trk_ip_chi2",&trk_ipchi2);    
	t->SetBranchAddress("trk_pnn_k",&trk_pnnk);    
	t->SetBranchAddress("trk_pnn_pi",&trk_pnnpi);    
	t->SetBranchAddress("trk_pnn_p",&trk_pnnp);    
	t->SetBranchAddress("trk_pnn_mu",&trk_pnnmu);    
	t->SetBranchAddress("trk_is_mu",&trk_ismu);    
	t->SetBranchAddress("trk_idx_jet",&trk_j);    
	t->SetBranchAddress("trk_pid",&trk_pid);    
	t->SetBranchAddress("trk_type",&trk_type);    

	vector<double> *neu_px = new vector<double>();
	vector<double> *neu_py = new vector<double>();
	vector<double> *neu_pz = new vector<double>();
	vector<double> *neu_e = new vector<double>();
	vector<double> *neu_j = new vector<double>();
	t->SetBranchAddress("neu_px",&neu_px);  
	t->SetBranchAddress("neu_py",&neu_py);  
	t->SetBranchAddress("neu_pz",&neu_pz);  
	t->SetBranchAddress("neu_e",&neu_e);   
	t->SetBranchAddress("neu_idx_jet",&neu_j);    

	int nent = t->GetEntries();
	//nent = 1e6;
	cout << nent << endl;
	boost::progress_display show_progress( nent );

	//  int noTOSFoundA(0), oneTOSFoundA(0), sameTOSFoundA(0), okTOSA(0);
	//  int noTOSFoundB(0), oneTOSFoundB(0), sameTOSFoundB(0), okTOSB(0);
	//  int methodsMatch(0), methodsDiffer(0);

	int nFoundSVTOS(0), nFoundSameSVTOS(0), nFoundOneSVTOS(0), nFoundBothSVTOS(0), nFoundSVForSVOnline(0), nFoundSVForSVOffline(0), nFoundSVForSVOnOff(0);
	int nFoundSVSVTOS(0), nFoundSameSVSVTOS(0), nFoundOneSVSVTOS(0), nFoundBothSVSVTOS(0), nFoundSV1ForSVSVOnline(0), nFoundSV1ForSVSVOffline(0), nFoundSV1ForSVSVOnOff(0), nFoundSV2ForSVSVOnline(0), nFoundSV2ForSVSVOffline(0), nFoundSV2ForSVSVOnOff(0), nFound1SVForSVSVOnline(0), nFound1SVForSVSVOffline(0), nFound1SVForSVSVOnOff(0), nFound2SVForSVSVOnline(0), nFound2SVForSVSVOffline(0), nFound2SVForSVSVOnOff(0);
	int nFoundSVMuTOS(0), nFoundSameSVMuTOS(0), nFoundOneSVMuTOS(0), nFoundBothSVMuTOS(0), nFoundSVForSVMuOnline(0), nFoundSVForSVMuOffline(0), nFoundSVForSVMuOnOff(0), nFoundMuForSVMuOnline(0), nFoundMuForSVMuOffline(0), nFoundMuForSVMuOnOff(0), nFoundSVMuForSVMuOnline(0), nFoundSVMuForSVMuOffline(0), nFoundSVMuForSVMuOnOff(0); 
	int nFoundMuMuTOS(0), nFoundSameMuMuTOS(0), nFoundOneMuMuTOS(0), nFoundBothMuMuTOS(0), nFoundMu1ForMuMuOnline(0), nFoundMu1ForMuMuOffline(0), nFoundMu1ForMuMuOnOff(0), nFoundMu2ForMuMuOnline(0), nFoundMu2ForMuMuOffline(0), nFoundMu2ForMuMuOnOff(0), nFound1MuForMuMuOnline(0), nFound1MuForMuMuOffline(0), nFound1MuForMuMuOnOff(0), nFound2MuForMuMuOnline(0), nFound2MuForMuMuOffline(0), nFound2MuForMuMuOnOff(0);

	int njtot = 0, njsv = 0;
	for(e=0; e<nent; e++){
		++show_progress;

		bool foundSVTOS(false), foundSameSVTOS(false), foundOneSVTOS(false), foundBothSVTOS(false), foundSVForSVOnline(false), foundSVForSVOffline(false), foundSVForSVOnOff(false);
		bool foundSVSVTOS(false), foundSameSVSVTOS(false), foundOneSVSVTOS(false), foundBothSVSVTOS(false), foundSV1ForSVSVOnline(false), foundSV1ForSVSVOffline(false), foundSV1ForSVSVOnOff(false), foundSV2ForSVSVOnline(false), foundSV2ForSVSVOffline(false), foundSV2ForSVSVOnOff(false);
		bool foundSVMuTOS(false), foundSameSVMuTOS(false), foundOneSVMuTOS(false), foundBothSVMuTOS(false), foundSVForSVMuOnline(false), foundSVForSVMuOffline(false), foundSVForSVMuOnOff(false), foundMuForSVMuOnline(false), foundMuForSVMuOffline(false), foundMuForSVMuOnOff(false), foundSVMuForSVMuOnline(false), foundSVMuForSVMuOffline(false), foundSVMuForSVMuOnOff(false);
		bool foundMuMuTOS(false), foundSameMuMuTOS(false), foundOneMuMuTOS(false), foundBothMuMuTOS(false), foundMu1ForMuMuOnline(false), foundMu1ForMuMuOffline(false), foundMu1ForMuMuOnOff(false), foundMu2ForMuMuOnline(false), foundMu2ForMuMuOffline(false), foundMu2ForMuMuOnOff(false);

		t->GetEntry(e);
		if(npv != 1) continue;

		//    int ng = gen_pid->size();
		//    if(ng < 2) continue;
		int nj = jet_e->size();
		int nsv = svz->size();
		int ntrk = trk_e->size();
		int nneu = neu_e->size();

		NTRK=ntrk;
		NNEU=nneu;

		//TODO
		std::vector<TLorentzVector> p4j1;
		std::vector<TLorentzVector> p4j2;
		int ntrig = evt_dec->size();
		//std::cout << std::endl;
		for(int itrig=0; itrig<ntrig; ++itrig) {
			p4j1.push_back(TLorentzVector(evt_j1_px->at(itrig),evt_j1_py->at(itrig),evt_j1_pz->at(itrig),0));
			p4j2.push_back(TLorentzVector(evt_j2_px->at(itrig),evt_j2_py->at(itrig),evt_j2_pz->at(itrig),0));
			//std::cout << evt_dec->at(itrig) << std::endl << evt_j1_idx->at(itrig) << "\t" << evt_j2_idx->at(itrig) << std::endl;
		}

		for(int j=0; j<nj; j++){
			TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),
					jet_e->at(j));
			if(p4j.Pt() < 10e3) continue;

			njtot++;
			int ipv = jet_pv->at(j);
			TVector3 pv(pvx->at(ipv),pvy->at(ipv),pvz->at(ipv)); 

			int ihard = -1, imu = -1, nmu = 0, jnchr = 0, jnneu = 0, ndispl6 = 0, ndispl9 = 0, ndispl16 = 0;
			double ptd = 0, jetq = 0, ry = 0, rp = 0, m11 = 0, m12 = 0, m22 = 0, sumpt2 = 0, pnnmu_best = 0;
			TLorentzVector p4mu,p4hard;
			for(int i=0; i<ntrk; i++){
				if(trk_j->at(i) != j) continue;
				jnchr++;
				TLorentzVector p4trk(trk_px->at(i),trk_py->at(i),trk_pz->at(i),trk_e->at(i));
				int q = trk_q->at(i);
				jetq += q*p4trk.Pt();
				ptd += pow(p4trk.Pt(),2);
				if(p4trk.Pt() > p4hard.Pt()) {p4hard = p4trk; ihard = i;}
				if(trk_ismu->at(i) > 0 && trk_pnnmu->at(i) > 0.5 && p4trk.Pt() > 500){
					nmu++;
					//if(trk_pnnmu->at(i) > pnnmu_best) 
					if(p4trk.Pt() > p4mu.Pt()) {p4mu = p4trk; imu=i; pnnmu_best = trk_pnnmu->at(i);}
				}
				double dy = p4trk.Rapidity()-p4j.Rapidity();
				double dp = p4trk.DeltaPhi(p4j);
				double r = sqrt(dy*dy+dp*dp);
				ry += r*dy*p4trk.Pt();
				rp += r*dp*p4trk.Pt();
				m11 += pow(p4trk.Pt()*dy,2);
				m22 += pow(p4trk.Pt()*dp,2);
				m12 += pow(p4trk.Pt(),2)*dy*dp;
				sumpt2 += pow(p4trk.Pt(),2);
				if(p4trk.Pt() > 500){
					if(trk_ipchi2->at(i) > 6) ndispl6++;
					if(trk_ipchi2->at(i) > 9) ndispl9++;
					if(trk_ipchi2->at(i) > 16) ndispl16++;
				}
			}
			D0M = -1000.;
			D0PX      = -1000.;
			D0PY      = -1000.;
			D0PZ      = -1000.;
			D0E       = -1000.;
			D0X       = -1000.;
			D0Y       = -1000.;
			D0Z       = -1000.;
			D0IP      = -1000.;
			D0IPCHI2  = -1000;
			D0FD      = -1000.;
			D0FDCHI2  = -1000;
			D0TAU     = -1000;
			D0VTXCHI2 = -1000;
			D0VTXNDOF = -1000;
			D0PIPNN   = -1000;
			D0KPNN    = -1000;
			for(uint iD=0; iD<d0_mass->size(); ++iD) {
				TLorentzVector p4d(d0_px->at(iD),d0_py->at(iD),d0_pz->at(iD),d0_e->at(iD));
				if(p4d.DeltaR(p4j)<0.5) {
				//if(d0_j->at(iD) == j && d0_nj->at(iD)==2) {
					if(!(trk_pnnpi->at(d0_trk1->at(iD))>0.2 && trk_pnnk->at( d0_trk0->at(iD))>0.3)) continue;
					//if(!pidCut(trk_pnnpi->at(d0_trk1->at(iD)),trk_px->at(d0_trk1->at(iD)),trk_py->at(d0_trk1->at(iD))) ||
					//   !pidCut(trk_pnnk->at( d0_trk0->at(iD)),trk_px->at(d0_trk0->at(iD)),trk_py->at(d0_trk0->at(iD))) )
					//	continue;
					D0M = d0_mass->at(iD);
					D0PX      = d0_px->at(iD);
					D0PY      = d0_py->at(iD);
					D0PZ      = d0_pz->at(iD);
					D0E       = d0_e->at(iD);
					D0X       = d0_x->at(iD);
					D0Y       = d0_y->at(iD);
					D0Z       = d0_z->at(iD);
					D0IP      = d0_ip->at(iD);
					D0IPCHI2  = d0_ipchi2->at(iD);
					D0FD      = d0_fd->at(iD);
					D0FDCHI2  = d0_fdchi2->at(iD);
					D0TAU     = d0_tau->at(iD);
					D0VTXCHI2 = d0_vtxchi2->at(iD);
					D0VTXNDOF = d0_vtxndof->at(iD);
					D0PIPNN   = trk_pnnpi->at(d0_trk1->at(iD));
					D0KPNN    = trk_pnnk->at( d0_trk0->at(iD));
					break;
				}
			}
			DPMM = -1000.;
			DPMPX      = -1000.;
			DPMPY      = -1000.;
			DPMPZ      = -1000.;
			DPME       = -1000.;
			DPMX       = -1000.;
			DPMY       = -1000.;
			DPMZ       = -1000.;
			DPMIP      = -1000.;
			DPMIPCHI2  = -1000;
			DPMFD      = -1000.;
			DPMFDCHI2  = -1000;
			DPMTAU     = -1000;
			DPMVTXCHI2 = -1000;
			DPMVTXNDOF = -1000;
			DPMPI1PNN   = -1000;
			DPMPI2PNN   = -1000;
			DPMKPNN    = -1000;
			for(uint iD=0; iD<d_mass->size(); ++iD) {
				TLorentzVector p4d(d_px->at(iD),d_py->at(iD),d_pz->at(iD),d_e->at(iD));
				if(p4d.DeltaR(p4j)<0.5) {
				//if(d_j->at(iD) == j && d_nj->at(iD)==3) {
					if(!(trk_pnnpi->at(d_trk1->at(iD))>0.2 && trk_pnnpi->at(d_trk2->at(iD))>0.2 && trk_pnnk->at( d_trk0->at(iD))>0.3)) continue;
					//if(!pidCut(trk_pnnpi->at(d_trk2->at(iD)),trk_px->at(d_trk2->at(iD)),trk_py->at(d_trk2->at(iD))) ||
					//   !pidCut(trk_pnnpi->at(d_trk1->at(iD)),trk_px->at(d_trk1->at(iD)),trk_py->at(d_trk1->at(iD))) ||
					//   !pidCut(trk_pnnk->at( d_trk0->at(iD)),trk_px->at(d_trk0->at(iD)),trk_py->at(d_trk0->at(iD))) )
					//	continue;
					DPMM = d_mass->at(iD);
					DPMPX      = d_px->at(iD);
					DPMPY      = d_py->at(iD);
					DPMPZ      = d_pz->at(iD);
					DPME       = d_e->at(iD);
					DPMX       = d_x->at(iD);
					DPMY       = d_y->at(iD);
					DPMZ       = d_z->at(iD);
					DPMIP      = d_ip->at(iD);
					DPMIPCHI2  = d_ipchi2->at(iD);
					DPMFD      = d_fd->at(iD);
					DPMFDCHI2  = d_fdchi2->at(iD);
					DPMTAU     = d_tau->at(iD);
					DPMVTXCHI2 = d_vtxchi2->at(iD);
					DPMVTXNDOF = d_vtxndof->at(iD);
					DPMPI1PNN   = trk_pnnpi->at(d_trk1->at(iD));
					DPMPI2PNN   = trk_pnnpi->at(d_trk2->at(iD));
					DPMKPNN    = trk_pnnk->at( d_trk0->at(iD));
					break;
				}
			}
			DSM       = -1000.;
			DSPX      = -1000.;
			DSPY      = -1000.;
			DSPZ      = -1000.;
			DSE       = -1000.;
			DSX       = -1000.;
			DSY       = -1000.;
			DSZ       = -1000.;
			DSIP      = -1000.;
			DSIPCHI2  = -1000;
			DSFD      = -1000.;
			DSFDCHI2  = -1000;
			DSTAU     = -1000;
			DSVTXCHI2 = -1000;
			DSVTXNDOF = -1000;
			DSPIPNN   = -1000;
			DSK1PNN   = -1000;
			DSK2PNN   = -1000;
			DSPHIM    = -1000.;
			for(uint iD=0; iD<ds_mass->size(); ++iD) {
				TLorentzVector p4d(ds_px->at(iD),ds_py->at(iD),ds_pz->at(iD),ds_e->at(iD));
				if(p4d.DeltaR(p4j)<0.5) {
				//if(ds_j->at(iD) == j && ds_nj->at(iD)==3) {
					if(!(trk_pnnk->at(ds_trk0->at(iD))>0.3 && trk_pnnk->at(ds_trk1->at(iD))>0.3 && trk_pnnpi->at( ds_trk2->at(iD))>0.2)) continue;
					//if(!pidCut(trk_pnnpi->at(ds_trk2->at(iD)),trk_px->at(ds_trk2->at(iD)),trk_py->at(ds_trk2->at(iD))) ||
					//   !pidCut(trk_pnnk->at( ds_trk1->at(iD)),trk_px->at(ds_trk1->at(iD)),trk_py->at(ds_trk1->at(iD))) ||
					//   !pidCut(trk_pnnk->at( ds_trk0->at(iD)),trk_px->at(ds_trk0->at(iD)),trk_py->at(ds_trk0->at(iD))) )
					//	continue;
					TLorentzVector phi(trk_px->at(ds_trk0->at(iD))+trk_px->at(ds_trk1->at(iD)),
							   trk_py->at(ds_trk0->at(iD))+trk_py->at(ds_trk1->at(iD)),
							   trk_pz->at(ds_trk0->at(iD))+trk_pz->at(ds_trk1->at(iD)),
							   trk_e->at( ds_trk0->at(iD))+trk_e->at( ds_trk1->at(iD)));
					DSM = ds_mass->at(iD);
					DSPX      = ds_px->at(iD);
					DSPY      = ds_py->at(iD);
					DSPZ      = ds_pz->at(iD);
					DSE       = ds_e->at(iD);
					DSX       = ds_x->at(iD);
					DSY       = ds_y->at(iD);
					DSZ       = ds_z->at(iD);
					DSIP      = ds_ip->at(iD);
					DSIPCHI2  = ds_ipchi2->at(iD);
					DSFD      = ds_fd->at(iD);
					DSFDCHI2  = ds_fdchi2->at(iD);
					DSTAU     = ds_tau->at(iD);
					DSVTXCHI2 = ds_vtxchi2->at(iD);
					DSVTXNDOF = ds_vtxndof->at(iD);
					DSPIPNN   = trk_pnnpi->at(ds_trk2->at(iD));
					DSK1PNN   = trk_pnnk->at( ds_trk0->at(iD));
					DSK2PNN   = trk_pnnk->at( ds_trk1->at(iD));
					DSPHIM = phi.M();
					break;
				}
			}
			LCM       = -1000.;
			LCPX      = -1000.;
			LCPY      = -1000.;
			LCPZ      = -1000.;
			LCE       = -1000.;
			LCX       = -1000.;
			LCY       = -1000.;
			LCZ       = -1000.;
			LCIP      = -1000.;
			LCIPCHI2  = -1000;
			LCFD      = -1000.;
			LCFDCHI2  = -1000;
			LCTAU     = -1000;
			LCVTXCHI2 = -1000;
			LCVTXNDOF = -1000;
			LCPIPNN   = -1000;
			LCKPNN    = -1000;
			LCPPNN    = -1000;
			for(uint iD=0; iD<lc_mass->size(); ++iD) {
				TLorentzVector p4d(lc_px->at(iD),lc_py->at(iD),lc_pz->at(iD),lc_e->at(iD));
				if(p4d.DeltaR(p4j)<0.5) {
				//if(lc_j->at(iD) == j && lc_nj->at(iD)==3) {
					if(!(trk_pnnp->at(lc_trk0->at(iD))>0.3 && trk_pnnk->at(lc_trk1->at(iD))>0.3 && trk_pnnpi->at( lc_trk2->at(iD))>0.2)) continue;
					//if(!pidCut(trk_pnnpi->at(lc_trk2->at(iD)),trk_px->at(lc_trk2->at(iD)),trk_py->at(lc_trk2->at(iD))) ||
					//   !pidCut(trk_pnnk->at( lc_trk1->at(iD)),trk_px->at(lc_trk1->at(iD)),trk_py->at(lc_trk1->at(iD))) ||
					//   !pidCut(trk_pnnp->at( lc_trk0->at(iD)),trk_px->at(lc_trk0->at(iD)),trk_py->at(lc_trk0->at(iD))) )
					//	continue;
					LCM       = lc_mass->at(iD);
					LCPX      = lc_px->at(iD);
					LCPY      = lc_py->at(iD);
					LCPZ      = lc_pz->at(iD);
					LCE       = lc_e->at(iD);
					LCX       = lc_x->at(iD);
					LCY       = lc_y->at(iD);
					LCZ       = lc_z->at(iD);
					LCIP      = lc_ip->at(iD);
					LCIPCHI2  = lc_ipchi2->at(iD);
					LCFD      = lc_fd->at(iD);
					LCFDCHI2  = lc_fdchi2->at(iD);
					LCTAU     = lc_tau->at(iD);
					LCVTXCHI2 = lc_vtxchi2->at(iD);
					LCVTXNDOF = lc_vtxndof->at(iD);
					LCPIPNN   = trk_pnnpi->at(lc_trk2->at(iD));
					LCPPNN    = trk_pnnp->at( lc_trk0->at(iD));
					LCKPNN    = trk_pnnk->at( lc_trk1->at(iD));
					break;
				}
			}
			D2K3PIM       = -1000.;
			D2K3PIPX      = -1000.;
			D2K3PIPY      = -1000.;
			D2K3PIPZ      = -1000.;
			D2K3PIE       = -1000.;
			D2K3PIX       = -1000.;
			D2K3PIY       = -1000.;
			D2K3PIZ       = -1000.;
			D2K3PIIP      = -1000.;
			D2K3PIIPCHI2  = -1000;
			D2K3PIFD      = -1000.;
			D2K3PIFDCHI2  = -1000;
			D2K3PITAU     = -1000;
			D2K3PIVTXCHI2 = -1000;
			D2K3PIVTXNDOF = -1000;
			D2K3PIPI1PNN  = -1000;
			D2K3PIPI2PNN  = -1000;
			D2K3PIPI3PNN  = -1000;
			D2K3PIKPNN   = -1000;
			for(uint iD=0; iD<d2k3pi_mass->size(); ++iD) {
				TLorentzVector p4d(d2k3pi_px->at(iD),d2k3pi_py->at(iD),d2k3pi_pz->at(iD),d2k3pi_e->at(iD));
				if(p4d.DeltaR(p4j)<0.5) {
				//if(d2k3pi_j->at(iD) == j && d2k3pi_nj->at(iD)==4) {
					if(!(trk_pnnk->at(d2k3pi_trk0->at(iD))>0.3 && trk_pnnpi->at(d2k3pi_trk1->at(iD))>0.2 && trk_pnnpi->at(d2k3pi_trk2->at(iD))>0.2 && trk_pnnpi->at( d2k3pi_trk3->at(iD))>0.2)) continue;
					//if(!pidCut(trk_pnnpi->at(d2k3pi_trk3->at(iD)),trk_px->at(d2k3pi_trk3->at(iD)),trk_py->at(d2k3pi_trk3->at(iD))) ||
					//   !pidCut(trk_pnnpi->at(d2k3pi_trk2->at(iD)),trk_px->at(d2k3pi_trk2->at(iD)),trk_py->at(d2k3pi_trk2->at(iD))) ||
					//   !pidCut(trk_pnnpi->at(d2k3pi_trk1->at(iD)),trk_px->at(d2k3pi_trk1->at(iD)),trk_py->at(d2k3pi_trk1->at(iD))) ||
					//   !pidCut(trk_pnnk->at( d2k3pi_trk0->at(iD)),trk_px->at(d2k3pi_trk0->at(iD)),trk_py->at(d2k3pi_trk0->at(iD))) )
					//	continue;
					D2K3PIM = d2k3pi_mass->at(iD);
					D2K3PIPX      = d2k3pi_px->at(iD);
					D2K3PIPY      = d2k3pi_py->at(iD);
					D2K3PIPZ      = d2k3pi_pz->at(iD);
					D2K3PIE       = d2k3pi_e->at(iD);
					D2K3PIX       = d2k3pi_x->at(iD);
					D2K3PIY       = d2k3pi_y->at(iD);
					D2K3PIZ       = d2k3pi_z->at(iD);
					D2K3PIIP      = d2k3pi_ip->at(iD);
					D2K3PIIPCHI2  = d2k3pi_ipchi2->at(iD);
					D2K3PIFD      = d2k3pi_fd->at(iD);
					D2K3PIFDCHI2  = d2k3pi_fdchi2->at(iD);
					D2K3PITAU     = d2k3pi_tau->at(iD);
					D2K3PIVTXCHI2 = d2k3pi_vtxchi2->at(iD);
					D2K3PIVTXNDOF = d2k3pi_vtxndof->at(iD);
					D2K3PIPI1PNN   = trk_pnnpi->at(d2k3pi_trk1->at(iD));
					D2K3PIPI2PNN   = trk_pnnpi->at(d2k3pi_trk2->at(iD));
					D2K3PIPI3PNN   = trk_pnnpi->at(d2k3pi_trk3->at(iD));
					D2K3PIKPNN     = trk_pnnk->at( d2k3pi_trk0->at(iD));
					break;
				}
			}

			jetq /= p4j.Pt();
			ptd = sqrt(ptd) / p4j.Pt();      
			for(int i=0; i<nneu; i++){
				if(neu_j->at(i) != j) continue; 
				jnneu++;
				TLorentzVector p4neu(neu_px->at(i),neu_py->at(i),neu_pz->at(i),neu_e->at(i));
				double dy = p4neu.Rapidity()-p4j.Rapidity();
				double dp = p4neu.DeltaPhi(p4j);
				double r = sqrt(dy*dy+dp*dp);
				ry += r*dy*p4neu.Pt();
				rp += r*dp*p4neu.Pt();
				m11 += pow(p4neu.Pt()*dy,2);
				m22 += pow(p4neu.Pt()*dp,2);
				m12 += pow(p4neu.Pt(),2)*dy*dp;
				sumpt2 += pow(p4neu.Pt(),2);
			}
			ry /= p4j.Pt();
			rp /= p4j.Pt();
			TMatrixD m(2,2);
			m(0,0)=m11;
			m(1,1)=m22;
			m(1,0)=m12;
			m(0,1)=m12;
			TVectorD eig;
			m.EigenVectors(eig);
			double js1=sqrt(eig[0]/sumpt2);
			double js2=sqrt(eig[1]/sumpt2);

			PARTON = -1; //parton;
			BMAXPT = -1; //bpt;
			CMAXPT = -1; //cpt;
			JPX = p4j.Px();
			JPY = p4j.Py();
			JPZ = p4j.Pz();
			JE = p4j.E();
			JPT = p4j.Pt();
			JETA = p4j.Eta();
			JS1 = js1; 
			JS2 = js2; 
			JQ = jetq; 
			JN = jnchr + jnneu; 
			JNQ = jnchr;
			JNN = jnneu;
			JPTD = ptd;

			JDIJETDEC=0;
			JDIJETSVDEC=0;
			JDIJETSVSVDEC=0;
			JDIJETSVMUDEC=0;
			JDIJETMUMUDEC=0;

			for(int itrig=0; itrig<ntrig; ++itrig) {
				if(evt_dec->at(itrig) == 0) {
					if(evt_j1_idx->at(itrig) == j || evt_j2_idx->at(itrig) == j) {
						JDIJETDEC=2;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nsv->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nsv->at(itrig) > 0))
							JDIJETDEC+=4;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nmu->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nmu->at(itrig) > 0))
							JDIJETDEC+=8;
						if(evt_j1_idx->at(itrig) == j && evt_j2_idx->at(itrig) == j)
							JDIJETDEC+=16;
					} else {
						JDIJETDEC=1;
					}
				} else if(evt_dec->at(itrig) == 1) {
					if(evt_j1_idx->at(itrig) == j || evt_j2_idx->at(itrig) == j) {
						JDIJETSVDEC=2;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nsv->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nsv->at(itrig) > 0))
							JDIJETSVDEC+=4;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nmu->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nmu->at(itrig) > 0))
							JDIJETSVDEC+=8;
						if(evt_j1_idx->at(itrig) == j && evt_j2_idx->at(itrig) == j)
							JDIJETSVDEC+=16;
					} else {
						JDIJETSVDEC=1;
					}
				} else if(evt_dec->at(itrig) == 2) {
					if(evt_j1_idx->at(itrig) == j || evt_j2_idx->at(itrig) == j) {
						JDIJETSVSVDEC=2;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nsv->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nsv->at(itrig) > 0))
							JDIJETSVSVDEC+=4;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nmu->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nmu->at(itrig) > 0))
							JDIJETSVSVDEC+=8;
						if(evt_j1_idx->at(itrig) == j && evt_j2_idx->at(itrig) == j)
							JDIJETSVSVDEC+=16;
					} else {
						JDIJETSVSVDEC=1;
					}
				} else if(evt_dec->at(itrig) == 3) {
					if(evt_j1_idx->at(itrig) == j || evt_j2_idx->at(itrig) == j) {
						JDIJETSVMUDEC=2;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nsv->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nsv->at(itrig) > 0))
							JDIJETSVMUDEC+=4;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nmu->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nmu->at(itrig) > 0))
							JDIJETSVMUDEC+=8;
						if(evt_j1_idx->at(itrig) == j && evt_j2_idx->at(itrig) == j)
							JDIJETSVMUDEC+=16;
					} else {
						JDIJETSVMUDEC=1;
					}
				} else if(evt_dec->at(itrig) == 4) {
					if(evt_j1_idx->at(itrig) == j || evt_j2_idx->at(itrig) == j) {
						JDIJETMUMUDEC=2;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nsv->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nsv->at(itrig) > 0))
							JDIJETMUMUDEC+=4;
						if((evt_j1_idx->at(itrig) == j && evt_j1_nmu->at(itrig) > 0) ||
						   (evt_j2_idx->at(itrig) == j && evt_j2_nmu->at(itrig) > 0))
							JDIJETMUMUDEC+=8;
						if(evt_j1_idx->at(itrig) == j && evt_j2_idx->at(itrig) == j)
							JDIJETMUMUDEC+=16;
					} else {
						JDIJETMUMUDEC=1;
					}
				}
			}
			PVX = pv.X();
			PVY = pv.Y();
			PVZ = pv.Z();
			NDISPL6 = ndispl6;
			NDISPL9 = ndispl9;
			NDISPL16 = ndispl16;
			NSV = 0; 
			int s = -1; double svptmax = -1;
			for(int ss = 0; ss < nsv; ss++){
				TVector3 sv(svx->at(ss),svy->at(ss),svz->at(ss));
				TVector3 fly = sv-pv;
				if(p4j.Vect().DeltaR(fly) < 0.5){
					TLorentzVector p4sv(svpx->at(ss),svpy->at(ss),svpz->at(ss),sve->at(ss));
					if(fly.Z()*p4sv.M()/p4sv.Pz()/(3e11)*(1e12) > 10) continue;	  
					NSV++;
					if(p4sv.Pt() > svptmax) {svptmax = p4sv.Pt(); s=ss;}
				}
			}
			if(svptmax < 1000) NSV = 0;
			if(NSV > 0){
				TVector3 sv(svx->at(s),svy->at(s),svz->at(s));
				TVector3 fly = sv-pv;
				SVX = sv.X();
				SVY = sv.Y();
				SVZ = sv.Z();
				SVPERP = sv.Perp();
				SVMCOR = svmcor->at(s); 
				SVMINPERP = svfdmin->at(s);
				SVDR = p4j.Vect().DeltaR(fly);
				SVN = 0; 
				SVNJ = 0;
				SVQ = 0; 
				SVSUMIPCHI2 = 0;
				SVMINIPCHI2 = 1e9;
				SVMAXGHOST = 0; 
				//	TLorentzVector p4sv(svpx->at(s),svpy->at(s),svpz->at(s),sve->at(s));
				TLorentzVector p4sv;
				bool veto=false;
				for(int i=0; i<10; i++){
					if(svtrk[i]->at(s) < 0) break;
					SVN++;
					int ii = svtrk[i]->at(s);
					TVector3 hit = TVector3(trk_x->at(ii),trk_y->at(ii),trk_z->at(ii));
					if(hit.Z() < sv.Z()) veto=true;
					TLorentzVector p4trk(trk_px->at(ii),trk_py->at(ii),trk_pz->at(ii),trk_e->at(ii));
					p4sv += p4trk;
					if(p4trk.DeltaR(p4j) < 0.5) SVNJ++;
					SVQ += trk_q->at(ii);
					double ipchi2 = trk_ipchi2->at(ii);
					//cout << "sv track IP chi2: " << ipchi2 << endl;
					SVSUMIPCHI2 += ipchi2;
					if(ipchi2 < SVMINIPCHI2) SVMINIPCHI2 = ipchi2;
					double ghost = trk_ghost->at(ii);
					if(ghost > SVMAXGHOST) SVMAXGHOST = ghost;
				}	
				SVPX = p4sv.Px();
				SVPY = p4sv.Py();
				SVPZ = p4sv.Pz();
				SVE = p4sv.E();
				SVPT = p4sv.Pt();
				SVETA = p4sv.Eta();
				SVM = p4sv.M();
				//	SVFDCHI2 = svfdchi2->at(s); 
				SVTZ = fly.Z()*p4sv.M()/p4sv.Pz()/(3e11)*(1e12);
				if(SVNJ < 1 || veto) NSV=0;
			}
			if(NSV==0){
				SVX = -1000;
				SVY = -1000;
				SVZ = -1000;
				SVPERP = -1000;
				SVPX = -1000;
				SVPY = -1000;
				SVPZ = -1000;
				SVE = -1000;
				SVPT = -1000;
				SVETA = -1000;
				SVM = -1000;
				SVMCOR = -1000;
				SVMINPERP = -1000;
				SVDR = -1000; 
				SVN = -1000; 
				SVNJ = -1000; 
				SVQ = -1000; 
				SVFDCHI2 = -1000; 
				SVSUMIPCHI2 = -1000; 
				SVTZ = -1000; 
				SVMINIPCHI2 = -1000; 
				SVMAXGHOST = -1000; 
			}
			else njsv++;
			NMU = nmu;
			if(NMU > 0){	
				MUPT = p4mu.Pt();
				MUIPCHI2 = trk_ipchi2->at(imu);
				MUDR = p4j.DeltaR(p4mu);
				MUPNN = trk_pnnmu->at(imu);
			}
			else{
				MUPT = -1000;
				MUIPCHI2 = -1000;
				MUDR = -1000;
				MUPNN = -1000;
			}
			if(p4hard.Pt() > 500){
				HPT = p4hard.Pt();
				HIPCHI2 = trk_ipchi2->at(ihard);
				HDR = p4hard.DeltaR(p4j);
			}
			else{
				HPT = -1000;
				HIPCHI2 = -1000;
				HDR = -1000;
			}


	//int nSV(0), nFoundOneSVTOS(0), nFOundBothSVTOS(0), nFoundSVForSV(0);
	//int nSVMu(0), nFoundOneSVMuTOS(0), nFoundBothSVMuTOS(0), nFoundMuForSVMu(0), nFoundSVForSVMu(0), nFoundSVMuForSVMu(0);
			for(int itrig=0; itrig<ntrig; ++itrig) {
				if(evt_dec->at(itrig) == 1) {
					foundSVTOS=true;
					if(evt_j1_idx->at(itrig) == evt_j2_idx->at(itrig)) {
						foundSameSVTOS=true;
					} else {
						if(evt_j1_idx->at(itrig) >= 0 || evt_j2_idx->at(itrig) >= 0) foundOneSVTOS=true;
						if(evt_j1_idx->at(itrig) >= 0 && evt_j2_idx->at(itrig) >= 0) {
							foundBothSVTOS=true;
							if(evt_j1_idx->at(itrig) == j) {
								if(evt_j1_nsv->at(itrig)>0  && NSV==0) foundSVForSVOnline=true;
								if(evt_j1_nsv->at(itrig)==0 && NSV>0) foundSVForSVOffline=true;
								if(evt_j1_nsv->at(itrig)>0  && NSV>0) foundSVForSVOnOff=true;
							} else if(evt_j2_idx->at(itrig) == j) {
								if(evt_j2_nsv->at(itrig)>0  && NSV==0) foundSVForSVOnline=true;
								if(evt_j2_nsv->at(itrig)==0 && NSV>0) foundSVForSVOffline=true;
								if(evt_j2_nsv->at(itrig)>0  && NSV>0) foundSVForSVOnOff=true;
							}
						}
					}
				} else if(evt_dec->at(itrig) == 2) {
					foundSVSVTOS=true;
					if(evt_j1_idx->at(itrig) == evt_j2_idx->at(itrig)) {
						foundSameSVSVTOS=true;
					} else {
						if(evt_j1_idx->at(itrig) >= 0 || evt_j2_idx->at(itrig) >= 0) foundOneSVSVTOS=true;
						if(evt_j1_idx->at(itrig) >= 0 && evt_j2_idx->at(itrig) >= 0) {
							foundBothSVSVTOS=true;
							if(evt_j1_idx->at(itrig) == j) {
								if(evt_j1_nsv->at(itrig)>0  && NSV==0) foundSV1ForSVSVOnline=true;
								if(evt_j1_nsv->at(itrig)==0 && NSV>0) foundSV1ForSVSVOffline=true;
								if(evt_j1_nsv->at(itrig)>0  && NSV>0) foundSV1ForSVSVOnOff=true;
							} else if(evt_j2_idx->at(itrig) == j) {
								if(evt_j2_nsv->at(itrig)>0  && NSV==0) foundSV2ForSVSVOnline=true;
								if(evt_j2_nsv->at(itrig)==0 && NSV>0) foundSV2ForSVSVOffline=true;
								if(evt_j2_nsv->at(itrig)>0  && NSV>0) foundSV2ForSVSVOnOff=true;
							}
						}
					}
					
				} else if(evt_dec->at(itrig) == 3) {
					foundSVMuTOS=true;
					if(evt_j1_idx->at(itrig) == evt_j2_idx->at(itrig)) {
						foundSameSVMuTOS=true;
					} else {
						if(evt_j1_idx->at(itrig) >= 0 || evt_j2_idx->at(itrig) >= 0) foundOneSVMuTOS=true;
						if(evt_j1_idx->at(itrig) >= 0 && evt_j2_idx->at(itrig) >= 0) {
							foundBothSVMuTOS=true;
							if(evt_j1_idx->at(itrig) == j) {
								//first check the ones that require multiple different jets
								if(foundSVForSVMuOnline  && evt_j1_nmu->at(itrig)>0  && NMU==0) foundSVMuForSVMuOnline=true;
								if(foundSVForSVMuOffline && evt_j1_nmu->at(itrig)==0 && NMU>0)  foundSVMuForSVMuOffline=true;
								if(foundSVForSVMuOnOff   && evt_j1_nmu->at(itrig)>0  && NMU>0)  foundSVMuForSVMuOnOff=true;
								if(foundMuForSVMuOnline  && evt_j1_nsv->at(itrig)>0  && NSV==0) foundSVMuForSVMuOnline=true;
								if(foundMuForSVMuOffline && evt_j1_nsv->at(itrig)==0 && NSV>0)  foundSVMuForSVMuOffline=true;
								if(foundMuForSVMuOnOff   && evt_j1_nsv->at(itrig)>0  && NSV>0)  foundSVMuForSVMuOnOff=true;
								//now the rest
								if(evt_j1_nsv->at(itrig)>0  && NSV==0) foundSVForSVMuOnline=true;
								if(evt_j1_nsv->at(itrig)==0 && NSV>0) foundSVForSVMuOffline=true;
								if(evt_j1_nsv->at(itrig)>0  && NSV>0) foundSVForSVMuOnOff=true;
								if(evt_j1_nmu->at(itrig)>0  && NMU==0) foundMuForSVMuOnline=true;
								if(evt_j1_nmu->at(itrig)==0 && NMU>0) foundMuForSVMuOffline=true;
								if(evt_j1_nmu->at(itrig)>0  && NMU>0) foundMuForSVMuOnOff=true;
							} else if(evt_j2_idx->at(itrig) == j) {
								//first check the ones that require multiple different jets
								if(foundSVForSVMuOnline  && evt_j1_nmu->at(itrig)>0  && NMU==0) foundSVMuForSVMuOnline=true;
								if(foundSVForSVMuOffline && evt_j1_nmu->at(itrig)==0 && NMU>0)  foundSVMuForSVMuOffline=true;
								if(foundSVForSVMuOnOff   && evt_j1_nmu->at(itrig)>0  && NMU>0)  foundSVMuForSVMuOnOff=true;
								if(foundMuForSVMuOnline  && evt_j1_nsv->at(itrig)>0  && NSV==0) foundSVMuForSVMuOnline=true;
								if(foundMuForSVMuOffline && evt_j1_nsv->at(itrig)==0 && NSV>0)  foundSVMuForSVMuOffline=true;
								if(foundMuForSVMuOnOff   && evt_j1_nsv->at(itrig)>0  && NSV>0)  foundSVMuForSVMuOnOff=true;
								//now the rest
								if(evt_j2_nsv->at(itrig)>0  && NSV==0) foundSVForSVMuOnline=true;
								if(evt_j2_nsv->at(itrig)==0 && NSV>0) foundSVForSVMuOffline=true;
								if(evt_j2_nsv->at(itrig)>0  && NSV>0) foundSVForSVMuOnOff=true;
								if(evt_j2_nmu->at(itrig)>0  && NMU==0) foundMuForSVMuOnline=true;
								if(evt_j2_nmu->at(itrig)==0 && NMU>0) foundMuForSVMuOffline=true;
								if(evt_j2_nmu->at(itrig)>0  && NMU>0) foundMuForSVMuOnOff=true;
							}
						}
					}
					
				} else if(evt_dec->at(itrig) == 4) {
					foundMuMuTOS=true;
					if(evt_j1_idx->at(itrig) == evt_j2_idx->at(itrig)) {
						foundSameMuMuTOS=true;
					} else {
						if(evt_j1_idx->at(itrig) >= 0 || evt_j2_idx->at(itrig) >= 0) foundOneMuMuTOS=true;
						if(evt_j1_idx->at(itrig) >= 0 && evt_j2_idx->at(itrig) >= 0) {
							foundBothMuMuTOS=true;
							if(evt_j1_idx->at(itrig) == j) {
								if(evt_j1_nmu->at(itrig)>0  && NMU==0) foundMu1ForMuMuOnline=true;
								if(evt_j1_nmu->at(itrig)==0 && NMU>0) foundMu1ForMuMuOffline=true;
								if(evt_j1_nmu->at(itrig)>0  && NMU>0) foundMu1ForMuMuOnOff=true;
							} else if(evt_j2_idx->at(itrig) == j) {
								if(evt_j2_nmu->at(itrig)>0  && NMU==0) foundMu2ForMuMuOnline=true;
								if(evt_j2_nmu->at(itrig)==0 && NMU>0) foundMu2ForMuMuOffline=true;
								if(evt_j2_nmu->at(itrig)>0  && NMU>0) foundMu2ForMuMuOnOff=true;
							}
						}
					}
					
				}
			}

			tout->Fill();
		}
		if(foundSVTOS) ++nFoundSVTOS;
		if(foundSameSVTOS) ++nFoundSameSVTOS;
		if(foundOneSVTOS) ++nFoundOneSVTOS;
		if(foundBothSVTOS) ++nFoundBothSVTOS;
		if(foundSVForSVOnline || foundSVForSVOnOff) ++nFoundSVForSVOnline;
		if(foundSVForSVOffline || foundSVForSVOnOff) ++nFoundSVForSVOffline;
		if(foundSVForSVOnOff) ++nFoundSVForSVOnOff;

		if(foundSVSVTOS) ++nFoundSVSVTOS;
		if(foundSameSVSVTOS) ++nFoundSameSVSVTOS;
		if(foundOneSVSVTOS) ++nFoundOneSVSVTOS;
		if(foundBothSVSVTOS) ++nFoundBothSVSVTOS;
		if((foundSV1ForSVSVOnline || foundSV1ForSVSVOnOff)) ++nFoundSV1ForSVSVOnline;
		if((foundSV1ForSVSVOffline || foundSV1ForSVSVOnOff)) ++nFoundSV1ForSVSVOffline;
		if(foundSV1ForSVSVOnOff) ++nFoundSV1ForSVSVOnOff;
		if((foundSV2ForSVSVOnline || foundSV2ForSVSVOnOff)) ++nFoundSV2ForSVSVOnline;
		if((foundSV2ForSVSVOffline || foundSV2ForSVSVOnOff)) ++nFoundSV2ForSVSVOffline;
		if(foundSV2ForSVSVOnOff) ++nFoundSV2ForSVSVOnOff;
		if((foundSV1ForSVSVOnline || foundSV1ForSVSVOnOff) || (foundSV2ForSVSVOnline || foundSV2ForSVSVOnOff)) ++nFound1SVForSVSVOnline;
		if((foundSV1ForSVSVOffline || foundSV1ForSVSVOnOff) || (foundSV2ForSVSVOffline || foundSV2ForSVSVOnOff)) ++nFound1SVForSVSVOffline;
		if(foundSV1ForSVSVOnOff || foundSV2ForSVSVOnOff) ++nFound1SVForSVSVOnOff;
		if((foundSV1ForSVSVOnline || foundSV1ForSVSVOnOff) && (foundSV2ForSVSVOnline || foundSV2ForSVSVOnOff)) ++nFound2SVForSVSVOnline;
		if((foundSV1ForSVSVOffline || foundSV1ForSVSVOnOff) && (foundSV2ForSVSVOffline || foundSV2ForSVSVOnOff)) ++nFound2SVForSVSVOffline;
		if(foundSV1ForSVSVOnOff && foundSV2ForSVSVOnOff) ++nFound2SVForSVSVOnOff;

		if(foundSVMuTOS) ++nFoundSVMuTOS;
		if(foundSameSVMuTOS) ++nFoundSameSVMuTOS;
		if(foundOneSVMuTOS) ++nFoundOneSVMuTOS;
		if(foundBothSVMuTOS) ++nFoundBothSVMuTOS;
		if(foundSVForSVMuOnline || foundSVForSVMuOnOff) ++nFoundSVForSVMuOnline;
		if(foundSVForSVMuOffline || foundSVForSVMuOnOff) ++nFoundSVForSVMuOffline;
		if(foundSVForSVMuOnOff) ++nFoundSVForSVMuOnOff;
		if(foundMuForSVMuOnline || foundMuForSVMuOnOff) ++nFoundMuForSVMuOnline;
		if(foundMuForSVMuOffline || foundMuForSVMuOnOff) ++nFoundMuForSVMuOffline;
		if(foundMuForSVMuOnOff) ++nFoundMuForSVMuOnOff;
		if(foundSVMuForSVMuOnline  || foundSVMuForSVMuOnOff) ++nFoundSVMuForSVMuOnline;
		if(foundSVMuForSVMuOffline || foundSVMuForSVMuOnOff) ++nFoundSVMuForSVMuOffline;
		if(foundSVMuForSVMuOnOff) ++nFoundSVMuForSVMuOnOff;

		if(foundMuMuTOS) ++nFoundMuMuTOS;
		if(foundSameMuMuTOS) ++nFoundSameMuMuTOS;
		if(foundOneMuMuTOS) ++nFoundOneMuMuTOS;
		if(foundBothMuMuTOS) ++nFoundBothMuMuTOS;
		if((foundMu1ForMuMuOnline || foundMu1ForMuMuOnOff)) ++nFoundMu1ForMuMuOnline;
		if((foundMu1ForMuMuOffline || foundMu1ForMuMuOnOff)) ++nFoundMu1ForMuMuOffline;
		if(foundMu1ForMuMuOnOff) ++nFoundMu1ForMuMuOnOff;
		if((foundMu2ForMuMuOnline || foundMu2ForMuMuOnOff)) ++nFoundMu2ForMuMuOnline;
		if((foundMu2ForMuMuOffline || foundMu2ForMuMuOnOff)) ++nFoundMu2ForMuMuOffline;
		if(foundMu2ForMuMuOnOff) ++nFoundMu2ForMuMuOnOff;
		if((foundMu1ForMuMuOnline || foundMu1ForMuMuOnOff) || (foundMu2ForMuMuOnline || foundMu2ForMuMuOnOff)) ++nFound1MuForMuMuOnline;
		if((foundMu1ForMuMuOffline || foundMu1ForMuMuOnOff) || (foundMu2ForMuMuOffline || foundMu2ForMuMuOnOff)) ++nFound1MuForMuMuOffline;
		if(foundMu1ForMuMuOnOff || foundMu2ForMuMuOnOff) ++nFound1MuForMuMuOnOff;
		if((foundMu1ForMuMuOnline || foundMu1ForMuMuOnOff) && (foundMu2ForMuMuOnline || foundMu2ForMuMuOnOff)) ++nFound2MuForMuMuOnline;
		if((foundMu1ForMuMuOffline || foundMu1ForMuMuOnOff) && (foundMu2ForMuMuOffline || foundMu2ForMuMuOnOff)) ++nFound2MuForMuMuOffline;
		if(foundMu1ForMuMuOnOff && foundMu2ForMuMuOnOff) ++nFound2MuForMuMuOnOff;
	}
	cout << njsv / (double) njtot << endl;
	//cout << noTOSFoundA << "\t" << oneTOSFoundA << "\t" << sameTOSFoundA << "\t" << okTOSA << std::endl;
	//  cout << noTOSFoundB << "\t" << oneTOSFoundB << "\t" << sameTOSFoundB << "\t" << okTOSB << std::endl;
	//  cout << methodsMatch << "\t" << methodsDiffer << std::endl;
	std::cout << std::endl;
	std::cout << nFoundSVTOS
	  << "\t" << nFoundSameSVTOS << "\t" << (double)nFoundSameSVTOS/nFoundSVTOS << std::endl;
	std::cout << nFoundOneSVTOS
	  << "\t" << nFoundBothSVTOS << "\t" << (double)nFoundBothSVTOS/nFoundSVTOS << std::endl;
	std::cout << nFoundSVForSVOnline
	  << "\t" << nFoundSVForSVOffline
	  << "\t" << nFoundSVForSVOnOff << "\t" << (double)nFoundSVForSVOnline/nFoundBothSVTOS << "\t" << (double)nFoundSVForSVOffline/nFoundBothSVTOS << std::endl;

	std::cout << std::endl;
	std::cout << nFoundSVSVTOS
	  << "\t" << nFoundSameSVSVTOS << "\t" << (double)nFoundSameSVSVTOS/nFoundSVSVTOS << std::endl;
	std::cout << nFoundOneSVSVTOS
	  << "\t" << nFoundBothSVSVTOS << "\t" << (double)nFoundBothSVSVTOS/nFoundSVSVTOS << std::endl;
	std::cout << nFoundSV1ForSVSVOnline
	  << "\t" << nFoundSV1ForSVSVOffline
	  << "\t" << nFoundSV1ForSVSVOnOff << "\t" << (double)nFoundSV1ForSVSVOnline/nFoundBothSVSVTOS << "\t" << (double)nFoundSV1ForSVSVOffline/nFoundBothSVSVTOS << std::endl;
	std::cout << nFoundSV2ForSVSVOnline
	  << "\t" << nFoundSV2ForSVSVOffline
	  << "\t" << nFoundSV2ForSVSVOnOff << "\t" << (double)nFoundSV2ForSVSVOnline/nFoundBothSVSVTOS << "\t" << (double)nFoundSV2ForSVSVOffline/nFoundBothSVSVTOS << std::endl;
	std::cout << nFound1SVForSVSVOnline
	  << "\t" << nFound1SVForSVSVOffline
	  << "\t" << nFound1SVForSVSVOnOff << "\t" << (double)nFound1SVForSVSVOnline/nFoundBothSVSVTOS << "\t" << (double)nFound1SVForSVSVOffline/nFoundBothSVSVTOS << std::endl;
	std::cout << nFound2SVForSVSVOnline
	  << "\t" << nFound2SVForSVSVOffline
	  << "\t" << nFound2SVForSVSVOnOff << "\t" << (double)nFound2SVForSVSVOnline/nFoundBothSVSVTOS << "\t" << (double)nFound2SVForSVSVOffline/nFoundBothSVSVTOS << std::endl;

	std::cout << std::endl;
	std::cout << nFoundSVMuTOS
	  << "\t" << nFoundSameSVMuTOS << "\t" << (double)nFoundSameSVMuTOS/nFoundSVMuTOS << std::endl;
	std::cout << nFoundOneSVMuTOS
	  << "\t" << nFoundBothSVMuTOS << "\t" << (double)nFoundBothSVMuTOS/nFoundSVMuTOS << std::endl;
	std::cout << nFoundSVForSVMuOnline
	  << "\t" << nFoundSVForSVMuOffline
	  << "\t" << nFoundSVForSVMuOnOff << "\t" << (double)nFoundSVForSVMuOnline/nFoundBothSVMuTOS << "\t" << (double)nFoundSVForSVMuOffline/nFoundBothSVMuTOS << std::endl;
	std::cout << nFoundMuForSVMuOnline
	  << "\t" << nFoundMuForSVMuOffline
	  << "\t" << nFoundMuForSVMuOnOff << "\t" << (double)nFoundMuForSVMuOnline/nFoundBothSVMuTOS << "\t" << (double)nFoundMuForSVMuOffline/nFoundBothSVMuTOS << std::endl;
	std::cout << nFoundSVMuForSVMuOnline
	  << "\t" << nFoundSVMuForSVMuOffline
	  << "\t" << nFoundSVMuForSVMuOnOff << "\t" << (double)nFoundSVMuForSVMuOnline/nFoundBothSVMuTOS << "\t" << (double)nFoundSVMuForSVMuOffline/nFoundBothSVMuTOS << std::endl;

	std::cout << std::endl;
	std::cout << nFoundMuMuTOS
	  << "\t" << nFoundSameMuMuTOS << "\t" << (double)nFoundSameMuMuTOS/nFoundMuMuTOS << std::endl;
	std::cout << nFoundOneMuMuTOS
	  << "\t" << nFoundBothMuMuTOS << "\t" << (double)nFoundBothMuMuTOS/nFoundMuMuTOS << std::endl;
	std::cout << nFoundMu1ForMuMuOnline
	  << "\t" << nFoundMu1ForMuMuOffline
	  << "\t" << nFoundMu1ForMuMuOnOff << "\t" << (double)nFoundMu1ForMuMuOnline/nFoundBothMuMuTOS << "\t" << (double)nFoundMu1ForMuMuOffline/nFoundBothMuMuTOS << std::endl;
	std::cout << nFoundMu2ForMuMuOnline
	  << "\t" << nFoundMu2ForMuMuOffline
	  << "\t" << nFoundMu2ForMuMuOnOff << "\t" << (double)nFoundMu2ForMuMuOnline/nFoundBothMuMuTOS << "\t" << (double)nFoundMu2ForMuMuOffline/nFoundBothMuMuTOS << std::endl;
	std::cout << nFound1MuForMuMuOnline
	  << "\t" << nFound1MuForMuMuOffline
	  << "\t" << nFound1MuForMuMuOnOff << "\t" << (double)nFound1MuForMuMuOnline/nFoundBothMuMuTOS << "\t" << (double)nFound1MuForMuMuOffline/nFoundBothMuMuTOS << std::endl;
	std::cout << nFound2MuForMuMuOnline
	  << "\t" << nFound2MuForMuMuOffline
	  << "\t" << nFound2MuForMuMuOnOff << "\t" << (double)nFound2MuForMuMuOnline/nFoundBothMuMuTOS << "\t" << (double)nFound2MuForMuMuOffline/nFoundBothMuMuTOS << std::endl;

	//fout.cd();
	//tout->Write("T",TObject::kOverwrite);
	tout->AutoSave();
	fout->Close();

	return 0;
}
