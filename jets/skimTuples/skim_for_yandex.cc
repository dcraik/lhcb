#include <iostream>
#include <vector>
#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorT.h"

using namespace std;

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

	int type = atoi(argv[1]);
	TRandom3 rand(9000+type);

	char str[1000];
	sprintf(str,"/tmp/dcraik/for_yandex_%d.root",type);
	//sprintf(str,"/tmp/dcraik/for_yandex_%d_resampled.root",type);
	TFile fout(str,"recreate");

	TTree *tout = new TTree("T","");
	int e(0);
	double JPX,JPY,JPZ,JE,JPT,JETA,JS1,JS2,JQ,JN,JNQ,JNN,JPTD;
	double PARTON,BMAXPT,CMAXPT;
	double SVM,SVMCOR,SVMCORERR,SVMCORERRFULL,SVMINPERP,SVPT,SVDR,SVN,SVNJ,SVQ,SVFDCHI2,SVPERP,SVETA,SVTZ,
	       SVMINIPCHI2,SVPX,SVPY,SVPZ,SVE,SVMAXGHOST,SVX,SVY,SVZ,SVSUMIPCHI2;
	double NSV,PVX,PVY,PVZ,NDISPL6,NDISPL9,NDISPL16,NADDDISPL6,NADDDISPL9,NADDDISPL16;
	double MUPT,MUIPCHI2,MUDR,MUPNN,NMU;
	double HPT,HIPCHI2,HDR;
	//double D0M, D0PX, D0PY, D0PZ, D0E, D0X, D0Y, D0Z, D0FD, D0DIRA, D0DOCA, D0IPCHI2MIN, D0DOCAKPI, D0VTXCHI2;
	//double DPMM, DPMPX, DPMPY, DPMPZ, DPME, DPMX, DPMY, DPMZ, DPMFD, DPMDIRA, DPMDOCA, DPMIPCHI2MIN, DPMDOCAMAX;
	//double DSM, DSPX, DSPY, DSPZ, DSE, DSX, DSY, DSZ, DSFD, DSDIRA, DSDOCA, DSIPCHI2MIN, DSDOCAMAX;
	//double LCM, LCPX, LCPY, LCPZ, LCE, LCX, LCY, LCZ, LCFD, LCDIRA, LCDOCA, LCIPCHI2MIN, LCDOCAMAX;
	//double D2K3PIM, D2K3PIPX, D2K3PIPY, D2K3PIPZ, D2K3PIE, D2K3PIX, D2K3PIY, D2K3PIZ, D2K3PIFD, D2K3PIDIRA, D2K3PIDOCA, D2K3PIIPCHI2MIN, D2K3PIDOCAMAX;
	double     D0M,     D0PX,     D0PY,     D0PZ,     D0E,     D0X,     D0Y,     D0Z,     D0FD,     D0FDCHI2,     D0IP,     D0IPCHI2,     D0VTXCHI2,     D0VTXNDOF,     D0TAU;
	double    DPMM,    DPMPX,    DPMPY,    DPMPZ,    DPME,    DPMX,    DPMY,    DPMZ,    DPMFD,    DPMFDCHI2,    DPMIP,    DPMIPCHI2,    DPMVTXCHI2,    DPMVTXNDOF,    DPMTAU;
	double     DSM,     DSPX,     DSPY,     DSPZ,     DSE,     DSX,     DSY,     DSZ,     DSFD,     DSFDCHI2,     DSIP,     DSIPCHI2,     DSVTXCHI2,     DSVTXNDOF,     DSTAU;
	double D2K3PIM, D2K3PIPX, D2K3PIPY, D2K3PIPZ, D2K3PIE, D2K3PIX, D2K3PIY, D2K3PIZ, D2K3PIFD, D2K3PIFDCHI2, D2K3PIIP, D2K3PIIPCHI2, D2K3PIVTXCHI2, D2K3PIVTXNDOF, D2K3PITAU;
	double NTRK, NNEU;

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
	tout->Branch("SVMCorErr",&SVMCORERR);
	tout->Branch("SVMCorErrFull",&SVMCORERRFULL);
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
	tout->Branch("NAddDispl6",&NADDDISPL6);
	tout->Branch("NAddDispl9",&NADDDISPL9);
	tout->Branch("NAddDispl16",&NADDDISPL16);
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
	//tout->Branch("D0M",        &D0M);
	//tout->Branch("D0PX",       &D0PX);
	//tout->Branch("D0PY",       &D0PY);
	//tout->Branch("D0PZ",       &D0PZ);
	//tout->Branch("D0E",        &D0E);
	//tout->Branch("D0X",        &D0X);
	//tout->Branch("D0Y",        &D0Y);
	//tout->Branch("D0Z",        &D0Z);
	//tout->Branch("D0FD",       &D0FD);
	//tout->Branch("D0DIRA",     &D0DIRA);
	//tout->Branch("D0IP",       &D0DOCA);
	//tout->Branch("D0IPCHI2MIN",    &D0IPCHI2MIN);
	//tout->Branch("D0DOCAKPI",  &D0DOCAKPI);
	//tout->Branch("D0VTXCHI2",  &D0VTXCHI2);
	//tout->Branch("DM",        &DPMM);
	//tout->Branch("DPX",       &DPMPX);
	//tout->Branch("DPY",       &DPMPY);
	//tout->Branch("DPZ",       &DPMPZ);
	//tout->Branch("DE",        &DPME);
	//tout->Branch("DX",        &DPMX);
	//tout->Branch("DY",        &DPMY);
	//tout->Branch("DZ",        &DPMZ);
	//tout->Branch("DFD",       &DPMFD);
	//tout->Branch("DDIRA",     &DPMDIRA);
	//tout->Branch("DIP",       &DPMDOCA);
	//tout->Branch("DIPCHI2MIN",    &DPMIPCHI2MIN);
	//tout->Branch("DDOCAMAX",  &DPMDOCAMAX);
	//tout->Branch("DSM",        &DSM);
	//tout->Branch("DSPX",       &DSPX);
	//tout->Branch("DSPY",       &DSPY);
	//tout->Branch("DSPZ",       &DSPZ);
	//tout->Branch("DSE",        &DSE);
	//tout->Branch("DSX",        &DSX);
	//tout->Branch("DSY",        &DSY);
	//tout->Branch("DSZ",        &DSZ);
	//tout->Branch("DSFD",       &DSFD);
	//tout->Branch("DSDIRA",     &DSDIRA);
	//tout->Branch("DSIP",       &DSDOCA);
	//tout->Branch("DSIPCHI2MIN",    &DSIPCHI2MIN);
	//tout->Branch("DSDOCAMAX",  &DSDOCAMAX);
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
	//tout->Branch("LCIP",       &LCDOCA);
	//tout->Branch("LCIPCHI2MIN",    &LCIPCHI2MIN);
	//tout->Branch("LCDOCAMAX",  &LCDOCAMAX);
	//tout->Branch("D2K3PIM",        &D2K3PIM);
	//tout->Branch("D2K3PIPX",       &D2K3PIPX);
	//tout->Branch("D2K3PIPY",       &D2K3PIPY);
	//tout->Branch("D2K3PIPZ",       &D2K3PIPZ);
	//tout->Branch("D2K3PIE",        &D2K3PIE);
	//tout->Branch("D2K3PIX",        &D2K3PIX);
	//tout->Branch("D2K3PIY",        &D2K3PIY);
	//tout->Branch("D2K3PIZ",        &D2K3PIZ);
	//tout->Branch("D2K3PIFD",       &D2K3PIFD);
	//tout->Branch("D2K3PIDIRA",     &D2K3PIDIRA);
	//tout->Branch("D2K3PIIP",       &D2K3PIDOCA);
	//tout->Branch("D2K3PIIPCHI2MIN",    &D2K3PIIPCHI2MIN);
	//tout->Branch("D2K3PIDOCAMAX",  &D2K3PIDOCAMAX);

	TChain *t = new TChain("data");
	//for(int i=1; i<41; ++i) {
	//        t->Add(TString::Format("/tmp/dcraik/lightjets_filtered_resampled%d.root",i));
	//}
	//  t->Add("/tmp/dcraik/jets_sim_mix.root");
	//  160 9 RIIJ_testMagDown49000000
	//  161 9 RIIJ_testMagDown49000001
	//  162 9 RIIJ_testMagDown49000002
	//  163 9 RIIJ_testMagDown49000003
	//  164 11 RIIJ_testMagDown49000004
	//  165 9 RIIJ_testMagDown49000040
	//  166 9 RIIJ_testMagDown49000041
	//  167 9 RIIJ_testMagDown49000042
	//  168 9 RIIJ_testMagDown49000043
	//  169 10 RIIJ_testMagDown49000044
	//  170 9 RIIJ_testMagDown49000050
	//  171 9 RIIJ_testMagDown49000051
	//  172 10 RIIJ_testMagDown49000052
	//  173 10 RIIJ_testMagDown49000053
	//  174 10 RIIJ_testMagDown49000054
	//  175 8 RIIJ_testMagUp49000000
	//  176 8 RIIJ_testMagUp49000001
	//  177 9 RIIJ_testMagUp49000002
	//  178 9 RIIJ_testMagUp49000003
	//  179 9 RIIJ_testMagUp49000004
	//  180 8 RIIJ_testMagUp49000040
	//  182 9 RIIJ_testMagUp49000041
	//  183 9 RIIJ_testMagUp49000042
	//  184 10 RIIJ_testMagUp49000043
	//  185 9 RIIJ_testMagUp49000044
	//  186 9 RIIJ_testMagUp49000050
	//  187 9 RIIJ_testMagUp49000051
	//  188 10 RIIJ_testMagUp49000052
	//  189 9 RIIJ_testMagUp49000053
	//  190 9 RIIJ_testMagUp49000054
	//int jobs0[] = {160,161,162,163,164,175,176,177,178,179};
	//int jobs4[] = {165,166,167,168,169,180,182,183,184,185};
	//int jobs5[] = {170,171,172,173,174,186,187,188,189,190};
	//int nsubjobs0[] = { 9, 9, 9, 9,11, 8, 8, 9, 9, 9};
	//int nsubjobs4[] = { 9, 9, 9, 9,10, 8, 9, 9,10, 9};
	//int nsubjobs5[] = { 9, 9,10,10,10, 9, 9,10, 9, 9};
	//256 RIIJ_testMagDown49000050 9
	//257 RIIJ_testMagDown49000051 9
	//258 RIIJ_testMagDown49000052 10
	//259 RIIJ_testMagDown49000053 10
	//260 RIIJ_testMagDown49000054 10
	//261 RIIJ_testMagDown49000040 9
	//262 RIIJ_testMagDown49000041 9
	//263 RIIJ_testMagDown49000042 9
	//264 RIIJ_testMagDown49000043 9
	//265 RIIJ_testMagDown49000044 10
	//266 RIIJ_testMagDown49000000 9
	//267 RIIJ_testMagDown49000001 9
	//268 RIIJ_testMagDown49000002 9
	//269 RIIJ_testMagDown49000003 9
	//270 RIIJ_testMagDown49000004 11
	//271 RIIJ_testMagUp49000050 9
	//272 RIIJ_testMagUp49000051 9
	//273 RIIJ_testMagUp49000052 10
	//274 RIIJ_testMagUp49000053 9
	//275 RIIJ_testMagUp49000054 9
	//276 RIIJ_testMagUp49000040 8
	//277 RIIJ_testMagUp49000041 9
	//278 RIIJ_testMagUp49000042 9
	//279 RIIJ_testMagUp49000043 10
	//280 RIIJ_testMagUp49000044 9
	//281 RIIJ_testMagUp49000000 8
	//282 RIIJ_testMagUp49000001 8
	//283 RIIJ_testMagUp49000002 9
	//284 RIIJ_testMagUp49000003 9
	//285 RIIJ_testMagUp49000004 9
	////TODO////int jobs0[] = {310,311,312,313,314};//{266,267,268,269,270,281,282,283,284,285};
	////TODO////int jobs4[] = {300,301,302,303,304};//{261,262,263,264,265,276,277,278,279,280};
	////TODO////int jobs5[] = {305,306,307,308,309};//{256,257,258,259,260,271,272,273,274,275};
	////TODO////int nsubjobs0[] = { 9, 9, 9, 9,11};//, 8, 8, 9, 9, 9};
	////TODO////int nsubjobs4[] = { 9, 9, 9, 9,10};//, 8, 9, 9,10, 9};
	////TODO////int nsubjobs5[] = { 9, 9,10,10,10};//, 9, 9,10, 9, 9};

	////TODO////boost::progress_display show_addfile_progress( 10 );
	////TODO////switch(type) {
	////TODO////	case 0:
	////TODO////		for(int i=0; i<5; ++i) {
	////TODO////			for(int j=0; j<nsubjobs0[i]; ++j) {
	////TODO////				sprintf(str,"/eos/lhcb/user/d/dcraik/jets/%d/%d/output.root",jobs0[i],j);
	////TODO////				if(gSystem->AccessPathName(str)) continue;
	////TODO////				t->Add(str);
	////TODO////			}
	////TODO////			++show_addfile_progress;
	////TODO////		}
	////TODO////		break;
	////TODO////	case 4:
	////TODO////		for(int i=0; i<5; ++i) {
	////TODO////			for(int j=0; j<nsubjobs4[i]; ++j) {
	////TODO////				sprintf(str,"/eos/lhcb/user/d/dcraik/jets/%d/%d/output.root",jobs4[i],j);
	////TODO////				if(gSystem->AccessPathName(str)) continue;
	////TODO////				t->Add(str);
	////TODO////			}
	////TODO////			++show_addfile_progress;
	////TODO////		}
	////TODO////		break;
	////TODO////	case 5:
	////TODO////		for(int i=0; i<5; ++i) {
	////TODO////			for(int j=0; j<nsubjobs5[i]; ++j) {
	////TODO////				sprintf(str,"/eos/lhcb/user/d/dcraik/jets/%d/%d/output.root",jobs5[i],j);
	////TODO////				if(gSystem->AccessPathName(str)) continue;
	////TODO////				t->Add(str);
	////TODO////			}
	////TODO////			++show_addfile_progress;
	////TODO////		}
	////TODO////		break;
	////TODO////}
	//t->Add("lightjets_filtered_addVars_resampled1.root");
	//for(int i=1; i<=4; i++){
	//  sprintf(str,"/Users/philten/data/Run2Jets.MC15.MU.490000%d%d.0.160925.00.root",type,i);
	//  sprintf(str,"/Users/philten/data/Run2Jets.MC15.MD.490000%d%d.0.160925.00.root",type,i);
	//  t->Add(str);
	//}

	switch(type) {
		case 0:
			t->Add("../davinci/light.root");
			break;
		case 4:
			t->Add("../davinci/charm.root");
			break;
		case 5:
			t->Add("../davinci/beauty.root");
			break;
	}

	vector<double> *gen_px = new vector<double>();
	vector<double> *gen_py = new vector<double>();
	vector<double> *gen_pz = new vector<double>();
	vector<double> *gen_e = new vector<double>();
	vector<double> *gen_pid = new vector<double>();
	t->SetBranchAddress("gen_px",&gen_px);  
	t->SetBranchAddress("gen_py",&gen_py);  
	t->SetBranchAddress("gen_pz",&gen_pz);  
	t->SetBranchAddress("gen_e",&gen_e);  
	t->SetBranchAddress("gen_pid",&gen_pid);  

	vector<double> *jet_px = new vector<double>();
	vector<double> *jet_py = new vector<double>();
	vector<double> *jet_pz = new vector<double>();
	vector<double> *jet_e = new vector<double>();
	vector<double> *jet_pv = new vector<double>();
	t->SetBranchAddress("jet_px",&jet_px);  
	t->SetBranchAddress("jet_py",&jet_py);  
	t->SetBranchAddress("jet_pz",&jet_pz);  
	t->SetBranchAddress("jet_e",&jet_e);  
	t->SetBranchAddress("jet_idx_pvr",&jet_pv);  

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
	vector<double> *svmcorerr = new vector<double>();  
	vector<double> *svmcorerrfull = new vector<double>();  
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
	t->SetBranchAddress("svr_m_cor_err",&svmcorerr);  
	t->SetBranchAddress("svr_m_cor_err_full",&svmcorerrfull);  
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
	vector<double> *trk_g = new vector<double>();
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
	t->SetBranchAddress("trk_idx_gen",&trk_g);    
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

	vector<double> *d0_mass     = new vector<double>();
	vector<double> *d_mass      = new vector<double>();
	vector<double> *ds_mass     = new vector<double>();
	vector<double> *d2k3pi_mass = new vector<double>();
	vector<double> *d0_j = new vector<double>();
	vector<double> *d_j = new vector<double>();
	vector<double> *ds_j = new vector<double>();
	vector<double> *d2k3pi_j = new vector<double>();
	vector<double> *d0_nj = new vector<double>();
	vector<double> *d_nj = new vector<double>();
	vector<double> *ds_nj = new vector<double>();
	vector<double> *d2k3pi_nj = new vector<double>();
	vector<double> *d0_trk0 = new vector<double>();
	vector<double> *d_trk0 = new vector<double>();
	vector<double> *ds_trk0 = new vector<double>();
	vector<double> *d2k3pi_trk0 = new vector<double>();
	vector<double> *d0_trk1 = new vector<double>();
	vector<double> *d_trk1 = new vector<double>();
	vector<double> *ds_trk1 = new vector<double>();
	vector<double> *d2k3pi_trk1 = new vector<double>();
	vector<double> *d_trk2 = new vector<double>();
	vector<double> *ds_trk2 = new vector<double>();
	vector<double> *d2k3pi_trk2 = new vector<double>();
	vector<double> *d2k3pi_trk3 = new vector<double>();
	vector<double> *d0_px     = new vector<double>();
	vector<double> *d_px      = new vector<double>();
	vector<double> *ds_px     = new vector<double>();
	vector<double> *d2k3pi_px = new vector<double>();
	vector<double> *d0_py     = new vector<double>();
	vector<double> *d_py      = new vector<double>();
	vector<double> *ds_py     = new vector<double>();
	vector<double> *d2k3pi_py = new vector<double>();
	vector<double> *d0_pz     = new vector<double>();
	vector<double> *d_pz      = new vector<double>();
	vector<double> *ds_pz     = new vector<double>();
	vector<double> *d2k3pi_pz = new vector<double>();
	vector<double> *d0_e     = new vector<double>();
	vector<double> *d_e      = new vector<double>();
	vector<double> *ds_e     = new vector<double>();
	vector<double> *d2k3pi_e = new vector<double>();
	vector<double> *d0_x     = new vector<double>();
	vector<double> *d_x      = new vector<double>();
	vector<double> *ds_x     = new vector<double>();
	vector<double> *d2k3pi_x = new vector<double>();
	vector<double> *d0_y     = new vector<double>();
	vector<double> *d_y      = new vector<double>();
	vector<double> *ds_y     = new vector<double>();
	vector<double> *d2k3pi_y = new vector<double>();
	vector<double> *d0_z     = new vector<double>();
	vector<double> *d_z      = new vector<double>();
	vector<double> *ds_z     = new vector<double>();
	vector<double> *d2k3pi_z = new vector<double>();
	vector<double> *d0_ip     = new vector<double>();
	vector<double> *d_ip      = new vector<double>();
	vector<double> *ds_ip     = new vector<double>();
	vector<double> *d2k3pi_ip = new vector<double>();
	vector<double> *d0_ipchi2     = new vector<double>();
	vector<double> *d_ipchi2      = new vector<double>();
	vector<double> *ds_ipchi2     = new vector<double>();
	vector<double> *d2k3pi_ipchi2 = new vector<double>();
	vector<double> *d0_fd     = new vector<double>();
	vector<double> *d_fd      = new vector<double>();
	vector<double> *ds_fd     = new vector<double>();
	vector<double> *d2k3pi_fd = new vector<double>();
	vector<double> *d0_fdchi2     = new vector<double>();
	vector<double> *d_fdchi2      = new vector<double>();
	vector<double> *ds_fdchi2     = new vector<double>();
	vector<double> *d2k3pi_fdchi2 = new vector<double>();
	vector<double> *d0_tau     = new vector<double>();
	vector<double> *d_tau      = new vector<double>();
	vector<double> *ds_tau     = new vector<double>();
	vector<double> *d2k3pi_tau = new vector<double>();
	vector<double> *d0_vtxchi2     = new vector<double>();
	vector<double> *d_vtxchi2      = new vector<double>();
	vector<double> *ds_vtxchi2     = new vector<double>();
	vector<double> *d2k3pi_vtxchi2 = new vector<double>();
	vector<double> *d0_vtxndof     = new vector<double>();
	vector<double> *d_vtxndof      = new vector<double>();
	vector<double> *ds_vtxndof     = new vector<double>();
	vector<double> *d2k3pi_vtxndof = new vector<double>();

	t->SetBranchAddress("d0_m",&d0_mass);  
	t->SetBranchAddress("dp_m",&d_mass);  
	t->SetBranchAddress("ds_m",&ds_mass);  
	t->SetBranchAddress("d02k3pi_m",&d2k3pi_mass);  

	t->SetBranchAddress("d0_idx_jet",&d0_j);  
	t->SetBranchAddress("dp_idx_jet",&d_j);  
	t->SetBranchAddress("ds_idx_jet",&ds_j);  
	t->SetBranchAddress("d02k3pi_idx_jet",&d2k3pi_j);  

	t->SetBranchAddress("d0_ntrk_jet",&d0_nj);  
	t->SetBranchAddress("dp_ntrk_jet",&d_nj);  
	t->SetBranchAddress("ds_ntrk_jet",&ds_nj);  
	t->SetBranchAddress("d02k3pi_ntrk_jet",&d2k3pi_nj);  

	t->SetBranchAddress("d0_idx_trk0",&d0_trk0);  
	t->SetBranchAddress("dp_idx_trk0",&d_trk0);  
	t->SetBranchAddress("ds_idx_trk0",&ds_trk0);  
	t->SetBranchAddress("d02k3pi_idx_trk0",&d2k3pi_trk0);  

	t->SetBranchAddress("d0_idx_trk1",&d0_trk1);  
	t->SetBranchAddress("dp_idx_trk1",&d_trk1);  
	t->SetBranchAddress("ds_idx_trk1",&ds_trk1);  
	t->SetBranchAddress("d02k3pi_idx_trk1",&d2k3pi_trk1);  

	t->SetBranchAddress("dp_idx_trk2",&d_trk2);  
	t->SetBranchAddress("ds_idx_trk2",&ds_trk2);  
	t->SetBranchAddress("d02k3pi_idx_trk2",&d2k3pi_trk2);  

	t->SetBranchAddress("d02k3pi_idx_trk3",&d2k3pi_trk3);  

	t->SetBranchAddress(     "d0_px",     &d0_px);  
	t->SetBranchAddress(     "dp_px",      &d_px);  
	t->SetBranchAddress(     "ds_px",     &ds_px);  
	t->SetBranchAddress("d02k3pi_px", &d2k3pi_px);  

	t->SetBranchAddress(     "d0_py",     &d0_py);  
	t->SetBranchAddress(     "dp_py",      &d_py);  
	t->SetBranchAddress(     "ds_py",     &ds_py);  
	t->SetBranchAddress("d02k3pi_py", &d2k3pi_py);  

	t->SetBranchAddress(     "d0_pz",     &d0_pz);  
	t->SetBranchAddress(     "dp_pz",      &d_pz);  
	t->SetBranchAddress(     "ds_pz",     &ds_pz);  
	t->SetBranchAddress("d02k3pi_pz", &d2k3pi_pz);  

	t->SetBranchAddress(     "d0_e",     &d0_e);  
	t->SetBranchAddress(     "dp_e",      &d_e);  
	t->SetBranchAddress(     "ds_e",     &ds_e);  
	t->SetBranchAddress("d02k3pi_e", &d2k3pi_e);  

	t->SetBranchAddress(     "d0_x",     &d0_x);  
	t->SetBranchAddress(     "dp_x",      &d_x);  
	t->SetBranchAddress(     "ds_x",     &ds_x);  
	t->SetBranchAddress("d02k3pi_x", &d2k3pi_x);  

	t->SetBranchAddress(     "d0_y",     &d0_y);  
	t->SetBranchAddress(     "dp_y",      &d_y);  
	t->SetBranchAddress(     "ds_y",     &ds_y);  
	t->SetBranchAddress("d02k3pi_y", &d2k3pi_y);  

	t->SetBranchAddress(     "d0_z",     &d0_z);  
	t->SetBranchAddress(     "dp_z",      &d_z);  
	t->SetBranchAddress(     "ds_z",     &ds_z);  
	t->SetBranchAddress("d02k3pi_z", &d2k3pi_z);  

	t->SetBranchAddress(     "d0_ip",     &d0_ip);  
	t->SetBranchAddress(     "dp_ip",      &d_ip);  
	t->SetBranchAddress(     "ds_ip",     &ds_ip);  
	t->SetBranchAddress("d02k3pi_ip", &d2k3pi_ip);  

	t->SetBranchAddress(     "d0_ip_chi2",     &d0_ipchi2);  
	t->SetBranchAddress(     "dp_ip_chi2",      &d_ipchi2);  
	t->SetBranchAddress(     "ds_ip_chi2",     &ds_ipchi2);  
	t->SetBranchAddress("d02k3pi_ip_chi2", &d2k3pi_ipchi2);  

	t->SetBranchAddress(     "d0_fd",     &d0_fd);  
	t->SetBranchAddress(     "dp_fd",      &d_fd);  
	t->SetBranchAddress(     "ds_fd",     &ds_fd);  
	t->SetBranchAddress("d02k3pi_fd", &d2k3pi_fd);  

	t->SetBranchAddress(     "d0_fd_chi2",     &d0_fdchi2);  
	t->SetBranchAddress(     "dp_fd_chi2",      &d_fdchi2);  
	t->SetBranchAddress(     "ds_fd_chi2",     &ds_fdchi2);  
	t->SetBranchAddress("d02k3pi_fd_chi2", &d2k3pi_fdchi2);  

	t->SetBranchAddress(     "d0_tau",     &d0_tau);  
	t->SetBranchAddress(     "dp_tau",      &d_tau);  
	t->SetBranchAddress(     "ds_tau",     &ds_tau);  
	t->SetBranchAddress("d02k3pi_tau", &d2k3pi_tau);  

	t->SetBranchAddress(     "d0_vtx_chi2",     &d0_vtxchi2);  
	t->SetBranchAddress(     "dp_vtx_chi2",      &d_vtxchi2);  
	t->SetBranchAddress(     "ds_vtx_chi2",     &ds_vtxchi2);  
	t->SetBranchAddress("d02k3pi_vtx_chi2", &d2k3pi_vtxchi2);  

	t->SetBranchAddress(     "d0_vtx_ndof",     &d0_vtxndof);  
	t->SetBranchAddress(     "dp_vtx_ndof",      &d_vtxndof);  
	t->SetBranchAddress(     "ds_vtx_ndof",     &ds_vtxndof);  
	t->SetBranchAddress("d02k3pi_vtx_ndof", &d2k3pi_vtxndof);  

	int nent = t->GetEntries();
	//nent = 1e5;
	cout << nent << endl;
	boost::progress_display show_progress( nent );

	//Setup PID histograms
	TFile* fE = TFile::Open("../pidcalib/muons/PerfHists_Pi_Turbo16_MagDown_misIDAsMu_Brunel_P_Brunel_PT.root");//TODO
	TH2F* histE = dynamic_cast<TH2F*>(fE->Get("Pi_IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5_0"));//TODO
	TFile* fMu = TFile::Open("../pidcalib/muons/PerfHists_Mu_Turbo16_MagDown_muons5_Brunel_P_Brunel_PT.root");
	TH2F* histMu = dynamic_cast<TH2F*>(fMu->Get("Mu_IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5_0"));
	TFile* fPi = TFile::Open("../pidcalib/muons/PerfHists_Pi_Turbo16_MagDown_misIDAsMu_Brunel_P_Brunel_PT.root");
	TH2F* histPi = dynamic_cast<TH2F*>(fPi->Get("Pi_IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5_0"));
	TFile* fK = TFile::Open("../pidcalib/muons/PerfHists_K_Turbo16_MagDown_misIDAsMu_Brunel_P_Brunel_PT.root");
	TH2F* histK = dynamic_cast<TH2F*>(fK->Get("K_IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5_0"));
	TFile* fP = TFile::Open("../pidcalib/muons/PerfHists_P_Turbo16_MagDown_misIDAsMu_Brunel_P_Brunel_PT.root");
	TH2F* histP = dynamic_cast<TH2F*>(fP->Get("P_IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5_0"));

	int countE(0), countMu(0), countPi(0), countK(0), countP(0);
	double njmu(0.), njmuErr(0.), njmuFromMu(0.);
	int njmu2(0);

	int nInJet(0), nOutOfJet(0), nSVOutOfJet(0);
	int njtot = 0, njsv = 0;
	for(e=0; e<nent; e++){
		++show_progress;

		t->GetEntry(e);
		if(npv != 1) continue;

		int ng = gen_pid->size();
		if(ng < 2) continue;
		int nj = jet_e->size();
		int nsv = svz->size();
		int ntrk = trk_e->size();
		int nneu = neu_e->size();

		NTRK=ntrk;
		NNEU=nneu;

		for(int j=0; j<nj; j++){
			TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),
					jet_e->at(j));
			if(p4j.Pt() < 10e3) continue;
			TLorentzVector p4mcj;

			// match to true jet
			for(int g=0; g<ng; g++){
				if(fabs(gen_pid->at(g)) != 98) continue;
				TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
				if(p4g.DeltaR(p4j) < 0.4 && p4g.Pt() > p4mcj.Pt())
					p4mcj = p4g;
			}
			if(p4mcj.Pt() < 10e3) continue;
			if(p4mcj.Eta() < 2.2 || p4mcj.Eta() > 4.2) continue;      

			double parton=-1000;
			TLorentzVector p4g0(gen_px->at(0),gen_py->at(0),gen_pz->at(0),gen_e->at(0));
			if(p4g0.DeltaR(p4j) < 0.5){
				parton = gen_pid->at(0);
			}
			for(int g=1; g<ng; g++){
				TLorentzVector p4g1(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
				if(p4g1.DeltaR(p4j) < 0.5 && (parton == -1000 || p4g1.Pt() > p4g0.Pt())){
					parton = gen_pid->at(1);
					p4g0 = p4g1;
				}
			}

			if(type == 4 && gen_pid->at(0) == 5) continue;
			if(type == 4 && gen_pid->at(1) == 5) continue;
			if(type == 4 && fabs(parton) == 5) continue;
			if(type == 5 && fabs(parton) != 5) continue;
			if(type == 0 && (gen_pid->at(0) == 4 || gen_pid->at(0) == 5)) continue;
			if(type == 0 && (gen_pid->at(1) == 4 || gen_pid->at(1) == 5)) continue;
			if(type == 0 && fabs(parton) == 4) continue;
			if(type == 0 && fabs(parton) == 5) continue;

			double bpt=0,cpt=0;
			for(int g=0; g<ng; g++){
				double pid = fabs(gen_pid->at(g));
				bool isb = (pid==511 || pid==521 || pid==531 || pid==5122);
				bool isc = (pid==411 || pid==421 || pid==431 || pid==4122);
				if(!isb && !isc) continue;
				TLorentzVector p4g(gen_px->at(g),gen_py->at(g),gen_pz->at(g),gen_e->at(g));
				if(p4g.DeltaR(p4j) > 0.5) continue;
				if(isb && p4g.Pt() > bpt) bpt = p4g.Pt();
				if(isc && p4g.Pt() > cpt) cpt = p4g.Pt();
			}

			if(type == 5 && bpt < 2000) continue;
			if(type == 4 && bpt > 0) continue;
			if(type == 4 && cpt < 2000) continue;
			if(type == 0 && (bpt > 0 || cpt > 0)) continue;

			//cout << parton << " " << ng << " " << bpt << " " << cpt << endl;

			njtot++;
			int ipv = jet_pv->at(j);
			TVector3 pv(pvx->at(ipv),pvy->at(ipv),pvz->at(ipv)); 

			//double D0ptmax(0), Dptmax(0), Dsptmax(0), Lcptmax(0), D2K3piptmax(0);
			//int bestD0(0), bestD(0), bestDs(0), bestLc(0), bestD2K3pi(0);

			//int countD0(0);
			//vector<double> d0_m;
			//vector<double> d0_px;
			//vector<double> d0_py;
			//vector<double> d0_pz;
			//vector<double> d0_e;
			//vector<double> d0_x;
			//vector<double> d0_y;
			//vector<double> d0_z;
			//vector<double> d0_fd;
			//vector<double> d0_dira;
			//vector<double> d0_doca;
			//vector<double> d0_ipmin;
			//vector<double> d0_docaKpi;
			//vector<double> d0_vtxchi2;

			//int countDpm(0);
			//vector<double> d_m;
			//vector<double> d_px;
			//vector<double> d_py;
			//vector<double> d_pz;
			//vector<double> d_e;
			//vector<double> d_x;
			//vector<double> d_y;
			//vector<double> d_z;
			//vector<double> d_fd;
			//vector<double> d_dira;
			//vector<double> d_doca;
			//vector<double> d_ipmin;
			//vector<double> d_docamax;

			//int countDs(0);
			//vector<double> ds_m;
			//vector<double> ds_px;
			//vector<double> ds_py;
			//vector<double> ds_pz;
			//vector<double> ds_e;
			//vector<double> ds_x;
			//vector<double> ds_y;
			//vector<double> ds_z;
			//vector<double> ds_fd;
			//vector<double> ds_dira;
			//vector<double> ds_doca;
			//vector<double> ds_ipmin;
			//vector<double> ds_docamax;

			//int countLc(0);
			//vector<double> lc_m;
			//vector<double> lc_px;
			//vector<double> lc_py;
			//vector<double> lc_pz;
			//vector<double> lc_e;
			//vector<double> lc_x;
			//vector<double> lc_y;
			//vector<double> lc_z;
			//vector<double> lc_fd;
			//vector<double> lc_dira;
			//vector<double> lc_doca;
			//vector<double> lc_ipmin;
			//vector<double> lc_docamax;

			//int countD2K3pi(0);
			//vector<double> d2k3pi_m;
			//vector<double> d2k3pi_px;
			//vector<double> d2k3pi_py;
			//vector<double> d2k3pi_pz;
			//vector<double> d2k3pi_e;
			//vector<double> d2k3pi_x;
			//vector<double> d2k3pi_y;
			//vector<double> d2k3pi_z;
			//vector<double> d2k3pi_fd;
			//vector<double> d2k3pi_dira;
			//vector<double> d2k3pi_doca;
			//vector<double> d2k3pi_ipmin;
			//vector<double> d2k3pi_docamax;

			int ihard = -1, imu = -1, nmu = 0, jnchr = 0, jnneu = 0, ndispl6 = 0, ndispl9 = 0, ndispl16 = 0;
			double ptd = 0, jetq = 0, ry = 0, rp = 0, m11 = 0, m12 = 0, m22 = 0, sumpt2 = 0, pnnmu_best = 0;
			TLorentzVector p4mu,p4hard;

			double probMu(0.), probMuErr(0.);
			double nMuExpected(0.), probNoMu(1.), probNoMuErr(0.), probNoMuErrSq(0.);
			int pidBin(0);
			int countE2(0), countMu2(0), countPi2(0), countK2(0), countP2(0);

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
				if(p4trk.Pt() > 500) {//probabilistic muon count
					//get true PID from gen - if no gen associated then use PID from trk
					int pid = TMath::Abs(trk_pid->at(i));
					if(trk_g->at(i) != -1) pid = TMath::Abs(gen_pid->at(trk_g->at(i)));

					switch(pid) {
						case 11:
							pidBin = 0.;
							probMu = 0.;
							probMuErr = 0.;
							++countE;
							++countE2;
							break;
						case 13:
							//std::cout << TMath::Abs(trk_pid->at(i));
							//if(trk_g->at(i) != -1) std::cout << "\t" << TMath::Abs(gen_pid->at(trk_g->at(i)));
							//std::cout << std::endl;
							pidBin = histMu->FindBin(p4trk.P(), p4trk.Pt());
							probMu = histMu->GetBinContent(pidBin);
							probMuErr = histMu->GetBinError(pidBin);
							++countMu;
							++countMu2;
							break;
						case 211:
							pidBin = histPi->FindBin(p4trk.P(), p4trk.Pt());
							probMu = histPi->GetBinContent(pidBin);
							probMuErr = histPi->GetBinError(pidBin);
							++countPi;
							++countPi2;
							break;
						case 321:
							pidBin = histK->FindBin(p4trk.P(), p4trk.Pt());
							probMu = histK->GetBinContent(pidBin);
							probMuErr = histK->GetBinError(pidBin);
							++countK;
							++countK2;
							break;
						case 2212:
							pidBin = histP->FindBin(p4trk.P(), p4trk.Pt());
							probMu = histP->GetBinContent(pidBin);
							probMuErr = histP->GetBinError(pidBin);
							++countP;
							++countP2;
							break;
						default:
							probMu = 0.;
							probMuErr = 0.;
							//std::cout << "Unknown track PID: " << pid << std::endl;
					}

					if(probMu>1.-probMuErr) probMu=1.-probMuErr;
					nMuExpected += probMu;
					probNoMu *= (1. - probMu);
					if(probMu!=1) {
						probNoMuErrSq += TMath::Power(probMuErr/(1.-probMu),2.);
					}
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
				//all of the combinations we're going to try have opp charged Kpi pair so start with that
				//if(trk_pnnk->at(i)>0.3 /*&& TMath::Abs(trk_pid->at(i))==321*/ && trk_ipchi2->at(i)>16. && trk_type->at(i)==3) {
				//  	TVector3 xtrk1 = TVector3(trk_x->at(i),trk_y->at(i),trk_z->at(i));
				//  	TLorentzVector p4trk1(trk_px->at(i),trk_py->at(i),trk_pz->at(i),trk_e->at(i));
				//	p4trk1.SetE(TMath::Sqrt(p4trk1.P()*p4trk1.P() + 493.7*493.7));
				//	//make D0 candidates
				//	for(int ii=0; ii<ntrk; ++ii) {//try every pair twice as first picked is the kaon
				//		if(trk_j->at(ii) != j) continue;
				//		if(ii==i) continue;
				//		if(trk_pid->at(i)*trk_pid->at(ii) > 0) continue;
				//		if(trk_pnnpi->at(ii)>0.3 /*&& TMath::Abs(trk_pid->at(ii))==211*/ && trk_ipchi2->at(ii)>16. && trk_type->at(ii)==3) {
				//  			TVector3 xtrk2 = TVector3(trk_x->at(ii),trk_y->at(ii),trk_z->at(ii));
				//  			TLorentzVector p4trk2(trk_px->at(ii),trk_py->at(ii),trk_pz->at(ii),trk_e->at(ii));
				//			p4trk2.SetE(TMath::Sqrt(p4trk2.P()*p4trk2.P() + 139.6*139.6));

				//			//also try to make D+ candidates
				//			for(int iii=ii+1; iii<ntrk; ++iii) {
				//				if(trk_j->at(iii) != j) continue;
				//				if(iii==i) continue;
				//				if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//same charge pions
				//				if(trk_pnnpi->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==211*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
				//  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
				//  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
				//					p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 139.6*139.6));
				//					TLorentzVector p4Dpm = p4trk1 + p4trk2 + p4trk3;
				//					if(TMath::Abs(p4Dpm.M()-1870.) > 160.) continue;
				//					TVector3 sv12, sv13, sv23, sv123;
				//					double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				//					double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
				//					double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
				//					double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
				//					double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
				//					if(docamax>0.1) continue;
				//					sv123 = sv12 + sv13 + sv23;
				//					sv123 *= (1./3.);
				//					double docaDpm = calcDocaPoint(sv123, p4Dpm.Vect(), pv);
				//					TVector3 fvDpm = sv123-pv;
				//					double diraDpm = fvDpm.Unit().Dot(p4Dpm.Vect().Unit());
				//					//printf("found D+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Dpm.M(), docamax, diraDpm);
				//					++countDpm;
				//					d_m.push_back(p4Dpm.M());
				//					d_px.push_back(p4Dpm.Px());
				//					d_py.push_back(p4Dpm.Py());
				//					d_pz.push_back(p4Dpm.Pz());
				//					d_e.push_back(p4Dpm.E());
				//					d_x.push_back(sv123.X());
				//					d_y.push_back(sv123.Y());
				//					d_z.push_back(sv123.Z());
				//					d_fd.push_back(fvDpm.Mag());
				//					d_dira.push_back(diraDpm);
				//					d_doca.push_back(docaDpm);
				//					d_ipmin.push_back(ipmin);
				//					d_docamax.push_back(docamax);

				//					if(p4Dpm.Pt() > Dptmax) {
				//						Dptmax = p4Dpm.Pt();
				//						bestD = d_m.size()-1;
				//					}
				//				}
				//			}

				//			//also try to make Ds+ candidates
				//			for(int iii=0; iii<ntrk; ++iii) {
				//				if(trk_j->at(iii) != j) continue;
				//				if(iii==i || iii==ii) continue;
				//				if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//opp charged kaons
				//				if(trk_pnnk->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
				//  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
				//  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
				//					p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 493.7*493.7));
				//					TLorentzVector p4Ds = p4trk1 + p4trk2 + p4trk3;
				//					if(TMath::Abs(p4Ds.M()-1968.) > 160.) continue;
				//					TVector3 sv12, sv13, sv23, sv123;
				//					double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				//					double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
				//					double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
				//					double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
				//					double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
				//					if(docamax>0.1) continue;
				//					sv123 = sv12 + sv13 + sv23;
				//					sv123 *= (1./3.);
				//					double docaDs = calcDocaPoint(sv123, p4Ds.Vect(), pv);
				//					TVector3 fvDs = sv123-pv;
				//					double diraDs = fvDs.Unit().Dot(p4Ds.Vect().Unit());
				//					//printf("found Ds+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Ds.M(), docamax, diraDs);
				//					++countDs;
				//					ds_m.push_back(p4Ds.M());
				//					ds_px.push_back(p4Ds.Px());
				//					ds_py.push_back(p4Ds.Py());
				//					ds_pz.push_back(p4Ds.Pz());
				//					ds_e.push_back(p4Ds.E());
				//					ds_x.push_back(sv123.X());
				//					ds_y.push_back(sv123.Y());
				//					ds_z.push_back(sv123.Z());
				//					ds_fd.push_back(fvDs.Mag());
				//					ds_dira.push_back(diraDs);
				//					ds_doca.push_back(docaDs);
				//					ds_ipmin.push_back(ipmin);
				//					ds_docamax.push_back(docamax);

				//					if(p4Ds.Pt() > Dsptmax) {
				//						Dsptmax = p4Ds.Pt();
				//						bestDs = ds_m.size()-1;
				//					}
				//				}
				//			}

				//			//also try to make Lc+ candidates
				//			for(int iii=0; iii<ntrk; ++iii) {
				//				if(trk_j->at(iii) != j) continue;
				//				if(iii==i || iii==ii) continue;
				//				if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//opp charged kaon and proton
				//				if(trk_pnnp->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==2212*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
				//  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
				//  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
				//					p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 938.3*938.3));
				//					TLorentzVector p4Lc = p4trk1 + p4trk2 + p4trk3;
				//					if(TMath::Abs(p4Lc.M()-2286.) > 160.) continue;
				//					TVector3 sv12, sv13, sv23, sv123;
				//					double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				//					double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
				//					double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
				//					double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
				//					double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
				//					if(docamax>0.1) continue;
				//					sv123 = sv12 + sv13 + sv23;
				//					sv123 *= (1./3.);
				//					double docaLc = calcDocaPoint(sv123, p4Lc.Vect(), pv);
				//					TVector3 fvLc = sv123-pv;
				//					double diraLc = fvLc.Unit().Dot(p4Lc.Vect().Unit());
				//					//printf("found Lc+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Lc.M(), docamax, diraLc);
				//					++countLc;
				//					lc_m.push_back(p4Lc.M());
				//					lc_px.push_back(p4Lc.Px());
				//					lc_py.push_back(p4Lc.Py());
				//					lc_pz.push_back(p4Lc.Pz());
				//					lc_e.push_back(p4Lc.E());
				//					lc_x.push_back(sv123.X());
				//					lc_y.push_back(sv123.Y());
				//					lc_z.push_back(sv123.Z());
				//					lc_fd.push_back(fvLc.Mag());
				//					lc_dira.push_back(diraLc);
				//					lc_doca.push_back(docaLc);
				//					lc_ipmin.push_back(ipmin);
				//					lc_docamax.push_back(docamax);

				//					if(p4Lc.Pt() > Lcptmax) {
				//						Lcptmax = p4Lc.Pt();
				//						bestLc = lc_m.size()-1;
				//					}
				//				}
				//			}

				//			//also try to make D0->Kpipipi candidates
				//			for(int iii=0; iii<ntrk; ++iii) {
				//				if(trk_j->at(iii) != j) continue;
				//				if(iii==i || iii==ii) continue;
				//				if(trk_pnnpi->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
				//  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
				//  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
				//					p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 139.6*139.6));
				//					for(int iv=iii+1; iv<ntrk; ++iv) {
				//						if(trk_j->at(iv) != j) continue;
				//						if(iv==i || iv==ii) continue;
				//						if(trk_pid->at(iv)*trk_pid->at(iii) > 0) continue;//adding an opp charged pion pair
				//						if(trk_pnnpi->at(iv)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iv)>16. && trk_type->at(iv)==3) {
				//  							TVector3 xtrk4 = TVector3(trk_x->at(iv),trk_y->at(iv),trk_z->at(iv));
				//  							TLorentzVector p4trk4(trk_px->at(iv),trk_py->at(iv),trk_pz->at(iv),trk_e->at(iv));
				//							p4trk4.SetE(TMath::Sqrt(p4trk4.P()*p4trk4.P() + 139.6*139.6));
				//							TLorentzVector p4D2K3pi = p4trk1 + p4trk2 + p4trk3 + p4trk4;
				//							if(TMath::Abs(p4D2K3pi.M()-1864.) > 160.) continue;
				//							TVector3 sv12, sv13, sv14, sv23, sv24, sv34, sv1234;
				//							double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				//							double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
				//							double doca14 = calcDoca(sv14, xtrk1, p4trk1.Vect(), xtrk4, p4trk4.Vect());
				//							double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
				//							double doca24 = calcDoca(sv24, xtrk2, p4trk2.Vect(), xtrk4, p4trk4.Vect());
				//							double doca34 = calcDoca(sv34, xtrk3, p4trk3.Vect(), xtrk4, p4trk4.Vect());
				//							double docamax = TMath::Max(TMath::Max(doca12, TMath::Max(doca13,doca14)),TMath::Max(doca23, TMath::Max(doca24,doca34)));
				//							double ipmin = TMath::Min(TMath::Min(trk_ipchi2->at(i),trk_ipchi2->at(ii)), TMath::Min(trk_ipchi2->at(iii),trk_ipchi2->at(iv)));
				//							if(docamax>0.1) continue;
				//							sv1234 = sv12 + sv13 + sv14 + sv23 + sv24 + sv34;
				//							sv1234 *= (1./6.);
				//							double docaD2K3pi = calcDocaPoint(sv1234, p4D2K3pi.Vect(), pv);
				//							TVector3 fvD2K3pi = sv1234-pv;
				//							double diraD2K3pi = fvD2K3pi.Unit().Dot(p4D2K3pi.Vect().Unit());
				//							//printf("found D0->K3pi: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4D2K3pi.M(), docamax, diraD2K3pi);
				//							++countD2K3pi;
				//							d2k3pi_m.push_back(p4D2K3pi.M());
				//							d2k3pi_px.push_back(p4D2K3pi.Px());
				//							d2k3pi_py.push_back(p4D2K3pi.Py());
				//							d2k3pi_pz.push_back(p4D2K3pi.Pz());
				//							d2k3pi_e.push_back(p4D2K3pi.E());
				//							d2k3pi_x.push_back(sv1234.X());
				//							d2k3pi_y.push_back(sv1234.Y());
				//							d2k3pi_z.push_back(sv1234.Z());
				//							d2k3pi_fd.push_back(fvD2K3pi.Mag());
				//							d2k3pi_dira.push_back(diraD2K3pi);
				//							d2k3pi_doca.push_back(docaD2K3pi);
				//							d2k3pi_ipmin.push_back(ipmin);
				//							d2k3pi_docamax.push_back(docamax);

				//							if(p4D2K3pi.Pt() > D2K3piptmax) {
				//								D2K3piptmax = p4D2K3pi.Pt();
				//								bestD2K3pi = d2k3pi_m.size()-1;
				//							}
				//						}
				//					}
				//				}
				//			}

				//			//D0 candidates
				//			TLorentzVector p4D0 = p4trk1 + p4trk2;
				//			if(TMath::Abs(p4D0.M()-1864.) > 160.) continue;
				//			TVector3 sv;
				//			double docaKpi = calcDoca(sv, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				//			double ipmin = TMath::Min(trk_ipchi2->at(i), trk_ipchi2->at(ii));
				//			if(docaKpi>0.1) continue;
				//			double sigma2i = trk_ip->at(i)*trk_ip->at(i) / trk_ipchi2->at(i);
				//			double sigma2j = trk_ip->at(ii)*trk_ip->at(ii) / trk_ipchi2->at(ii);
				//			double vtxchi2 = docaKpi*docaKpi/(sigma2i+sigma2j);
				//			if(vtxchi2>10.) continue;
				//			double docaD0 = calcDocaPoint(sv, p4D0.Vect(), pv);
				//			TVector3 fv = sv-pv;
				//			double dira = fv.Unit().Dot(p4D0.Vect().Unit());
				//			//printf("found D0: mass=%.2f, doca=%.3f, dira=%.4f, vtxchi2=%.1f, KID=%.0f, pNNk=%.2f, piID=%.0f, pNNpi=%.2f\n", p4D0.M(), docaKpi, dira, vtxchi2, trk_pid->at(i), trk_pnnk->at(i), trk_pid->at(ii), trk_pnnpi->at(ii));
				//			++countD0;
				//			d0_m.push_back(p4D0.M());
				//			d0_px.push_back(p4D0.Px());
				//			d0_py.push_back(p4D0.Py());
				//			d0_pz.push_back(p4D0.Pz());
				//			d0_e.push_back(p4D0.E());
				//			d0_x.push_back(sv.X());
				//			d0_y.push_back(sv.Y());
				//			d0_z.push_back(sv.Z());
				//			d0_fd.push_back(fv.Mag());
				//			d0_dira.push_back(dira);
				//			d0_doca.push_back(docaD0);
				//			d0_ipmin.push_back(ipmin);
				//			d0_docaKpi.push_back(docaKpi);
				//			d0_vtxchi2.push_back(vtxchi2);

				//			if(p4D0.Pt() > D0ptmax) {
				//				D0ptmax = p4D0.Pt();
				//				bestD0 = d0_m.size()-1;
				//			}
				//		}
				//	}
				//}
			}

			//if(probNoMu<0.9) {
			//	std::cout << probNoMu << "\t" << nMuExpected << "\t" << nmu << std::endl;
			//	std::cout << countE2 << "\t" << countMu2 << "\t" << countPi2 << "\t" << countK2 << "\t" << countP2 << std::endl;
			//}

			probNoMuErr = probNoMu*TMath::Sqrt(probNoMuErrSq);
			njmu += 1.-probNoMu;
			if(countMu2>0) njmuFromMu += (1.-probNoMu);
			njmuErr += probNoMuErr;
			if(nmu>0) ++njmu2;
			//std::cout << 1.-probNoMu << "+/-" << probNoMuErr << "\t" << nMuExpected << std::endl;
			//std::cout << countE << "\t" << countMu << "\t" << countPi << "\t" << countK << "\t" << countP << std::endl;

			//TODO
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
			for(uint iD=0; iD<d0_mass->size(); ++iD) {
				if(d0_j->at(iD) == j && d0_nj->at(iD)==2) {
					if(!(trk_pnnpi->at(d0_trk1->at(iD))>0.2 && trk_pnnk->at( d0_trk0->at(iD))>0.3)) continue;
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
			for(uint iD=0; iD<d_mass->size(); ++iD) {
				if(d_j->at(iD) == j && d_nj->at(iD)==3) {
					if(!(trk_pnnpi->at(d_trk1->at(iD))>0.2 && trk_pnnpi->at(d_trk2->at(iD))>0.2 && trk_pnnk->at( d_trk0->at(iD))>0.3)) continue;
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
			for(uint iD=0; iD<ds_mass->size(); ++iD) {
				if(ds_j->at(iD) == j && ds_nj->at(iD)==3) {
					if(!(trk_pnnk->at(ds_trk0->at(iD))>0.3 && trk_pnnk->at(ds_trk1->at(iD))>0.3 && trk_pnnpi->at( ds_trk2->at(iD))>0.2)) continue;
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
			for(uint iD=0; iD<d2k3pi_mass->size(); ++iD) {
				if(d2k3pi_j->at(iD) == j && d2k3pi_nj->at(iD)==4) {
					if(!(trk_pnnk->at(d2k3pi_trk0->at(iD))>0.3 && trk_pnnk->at(d2k3pi_trk1->at(iD))>0.2 && trk_pnnk->at(d2k3pi_trk2->at(iD))>0.2 && trk_pnnpi->at( d2k3pi_trk3->at(iD))>0.2)) continue;
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
					break;
				}
			}

			//pick a random D0
			//if(!d0_px.empty()) {
			//  int whichD0 = bestD0;//rand.Integer(d0_px.size());
			//  D0M       = d0_m[whichD0];
			//  D0PX      = d0_px[whichD0];
			//  D0PY      = d0_py[whichD0];
			//  D0PZ      = d0_pz[whichD0];
			//  D0E       = d0_e[whichD0];
			//  D0X       = d0_x[whichD0];
			//  D0Y       = d0_y[whichD0];
			//  D0Z       = d0_z[whichD0];
			//  D0FD      = d0_fd[whichD0];
			//  D0DIRA    = d0_dira[whichD0];
			//  D0DOCA    = d0_doca[whichD0];
			//  D0IPCHI2MIN= d0_ipmin[whichD0];
			//  D0DOCAKPI = d0_docaKpi[whichD0];
			//  D0VTXCHI2 = d0_vtxchi2[whichD0];
			//  //if(D0PX==0.) std::cout << D0M << "\t" << D0PX << "\t" << D0PY << "\t" << D0PZ << std::endl;
			//} else {
			//  D0M       = -1000.;
			//  D0PX      = -1000.;
			//  D0PY      = -1000.;
			//  D0PZ      = -1000.;
			//  D0E       = -1000.;
			//  D0X       = -1000.;
			//  D0Y       = -1000.;
			//  D0Z       = -1000.;
			//  D0FD      = -1000.;
			//  D0DIRA    = -1000.;
			//  D0DOCA    = -1000.;
			//  D0IPCHI2MIN= -1000;
			//  D0DOCAKPI = -1000.;
			//  D0VTXCHI2 = -1000.;
			//}
			////pick a random D+
			//if(!d_px.empty()) {
			//  int whichD = bestD;//rand.Integer(d_px.size());
			//  DPMM       = d_m[whichD];
			//  DPMPX      = d_px[whichD];
			//  DPMPY      = d_py[whichD];
			//  DPMPZ      = d_pz[whichD];
			//  DPME       = d_e[whichD];
			//  DPMX       = d_x[whichD];
			//  DPMY       = d_y[whichD];
			//  DPMZ       = d_z[whichD];
			//  DPMFD      = d_fd[whichD];
			//  DPMDIRA    = d_dira[whichD];
			//  DPMIPCHI2MIN= d_ipmin[whichD];
			//  DPMDOCA    = d_doca[whichD];
			//  DPMDOCAMAX = d_docamax[whichD];
			//} else {
			//  DPMM       = -1000.;
			//  DPMPX      = -1000.;
			//  DPMPY      = -1000.;
			//  DPMPZ      = -1000.;
			//  DPME       = -1000.;
			//  DPMX       = -1000.;
			//  DPMY       = -1000.;
			//  DPMZ       = -1000.;
			//  DPMFD      = -1000.;
			//  DPMDIRA    = -1000.;
			//  DPMIPCHI2MIN= -1000;
			//  DPMDOCA    = -1000.;
			//  DPMDOCAMAX = -1000.;
			//}
			////pick a random Ds+
			//if(!ds_px.empty()) {
			//  int whichDs = bestDs;//rand.Integer(ds_px.size());
			//  DSM       = ds_m[whichDs];
			//  DSPX      = ds_px[whichDs];
			//  DSPY      = ds_py[whichDs];
			//  DSPZ      = ds_pz[whichDs];
			//  DSE       = ds_e[whichDs];
			//  DSX       = ds_x[whichDs];
			//  DSY       = ds_y[whichDs];
			//  DSZ       = ds_z[whichDs];
			//  DSFD      = ds_fd[whichDs];
			//  DSDIRA    = ds_dira[whichDs];
			//  DSIPCHI2MIN= ds_ipmin[whichDs];
			//  DSDOCA    = ds_doca[whichDs];
			//  DSDOCAMAX = ds_docamax[whichDs];
			//} else {
			//  DSM       = -1000.;
			//  DSPX      = -1000.;
			//  DSPY      = -1000.;
			//  DSPZ      = -1000.;
			//  DSE       = -1000.;
			//  DSX       = -1000.;
			//  DSY       = -1000.;
			//  DSZ       = -1000.;
			//  DSFD      = -1000.;
			//  DSDIRA    = -1000.;
			//  DSIPCHI2MIN= -1000;
			//  DSDOCA    = -1000.;
			//  DSDOCAMAX = -1000.;
			//}
			////pick a random Lc+
			//if(!lc_px.empty()) {
			//  int whichLc = bestLc;//rand.Integer(lc_px.size());
			//  LCM       = lc_m[whichLc];
			//  LCPX      = lc_px[whichLc];
			//  LCPY      = lc_py[whichLc];
			//  LCPZ      = lc_pz[whichLc];
			//  LCE       = lc_e[whichLc];
			//  LCX       = lc_x[whichLc];
			//  LCY       = lc_y[whichLc];
			//  LCZ       = lc_z[whichLc];
			//  LCFD      = lc_fd[whichLc];
			//  LCDIRA    = lc_dira[whichLc];
			//  LCIPCHI2MIN= lc_ipmin[whichLc];
			//  LCDOCA    = lc_doca[whichLc];
			//  LCDOCAMAX = lc_docamax[whichLc];
			//} else {
			//  LCM       = -1000.;
			//  LCPX      = -1000.;
			//  LCPY      = -1000.;
			//  LCPZ      = -1000.;
			//  LCE       = -1000.;
			//  LCX       = -1000.;
			//  LCY       = -1000.;
			//  LCZ       = -1000.;
			//  LCFD      = -1000.;
			//  LCDIRA    = -1000.;
			//  LCIPCHI2MIN= -1000;
			//  LCDOCA    = -1000.;
			//  LCDOCAMAX = -1000.;
			//}
			////pick a random D0->K3pi
			//if(!d2k3pi_px.empty()) {
			//  int whichD0 = bestD2K3pi;//rand.Integer(d2k3pi_px.size());
			//  D2K3PIM       = d2k3pi_m[whichD0];
			//  D2K3PIPX      = d2k3pi_px[whichD0];
			//  D2K3PIPY      = d2k3pi_py[whichD0];
			//  D2K3PIPZ      = d2k3pi_pz[whichD0];
			//  D2K3PIE       = d2k3pi_e[whichD0];
			//  D2K3PIX       = d2k3pi_x[whichD0];
			//  D2K3PIY       = d2k3pi_y[whichD0];
			//  D2K3PIZ       = d2k3pi_z[whichD0];
			//  D2K3PIFD      = d2k3pi_fd[whichD0];
			//  D2K3PIDIRA    = d2k3pi_dira[whichD0];
			//  D2K3PIIPCHI2MIN= d2k3pi_ipmin[whichD0];
			//  D2K3PIDOCA    = d2k3pi_doca[whichD0];
			//  D2K3PIDOCAMAX = d2k3pi_docamax[whichD0];
			//} else {
			//  D2K3PIM       = -1000.;
			//  D2K3PIPX      = -1000.;
			//  D2K3PIPY      = -1000.;
			//  D2K3PIPZ      = -1000.;
			//  D2K3PIE       = -1000.;
			//  D2K3PIX       = -1000.;
			//  D2K3PIY       = -1000.;
			//  D2K3PIZ       = -1000.;
			//  D2K3PIFD      = -1000.;
			//  D2K3PIDIRA    = -1000.;
			//  D2K3PIIPCHI2MIN= -1000;
			//  D2K3PIDOCA    = -1000.;
			//  D2K3PIDOCAMAX = -1000.;
			//}

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

			PARTON = parton;
			BMAXPT = bpt;
			CMAXPT = cpt;
			JPX = p4mcj.Px();
			JPY = p4mcj.Py();
			JPZ = p4mcj.Pz();
			JE = p4mcj.E();
			JPT = p4j.Pt();
			JETA = p4j.Eta();
			JS1 = js1; 
			JS2 = js2; 
			JQ = jetq; 
			JN = jnchr + jnneu; 
			JNQ = jnchr;
			JNN = jnneu;
			JPTD = ptd;
			PVX = pv.X();
			PVY = pv.Y();
			PVZ = pv.Z();
			NDISPL6 = ndispl6;
			NDISPL9 = ndispl9;
			NDISPL16 = ndispl16;
			NADDDISPL6 = ndispl6;
			NADDDISPL9 = ndispl9;
			NADDDISPL16 = ndispl16;
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
				SVMCORERR = svmcorerr->at(s); 
				SVMCORERRFULL = svmcorerrfull->at(s); 
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
				bool foundOutOfJet=false;
				for(int i=0; i<10; i++){
					if(svtrk[i]->at(s) < 0) break;
					SVN++;
					int ii = svtrk[i]->at(s);
					TVector3 hit = TVector3(trk_x->at(ii),trk_y->at(ii),trk_z->at(ii));
					if(hit.Z() < sv.Z()) veto=true;
					TLorentzVector p4trk(trk_px->at(ii),trk_py->at(ii),trk_pz->at(ii),trk_e->at(ii));
					p4sv += p4trk;
					if(p4trk.DeltaR(p4j) < 0.5) SVNJ++;
					if(trk_j->at(ii) == j) ++nInJet;
					else {
						++nOutOfJet;
						foundOutOfJet=true;
					}
					SVQ += trk_q->at(ii);
					double ipchi2 = trk_ipchi2->at(ii);
					//cout << "sv track IP chi2: " << ipchi2 << endl;
					SVSUMIPCHI2 += ipchi2;
					if(ipchi2 < SVMINIPCHI2) SVMINIPCHI2 = ipchi2;
					double ghost = trk_ghost->at(ii);
					if(ghost > SVMAXGHOST) SVMAXGHOST = ghost;
					if(p4trk.Pt() > 500){//subtract off the tracks from our SV from the additional displaced tracks
						if(trk_ipchi2->at(i) > 6) NADDDISPL6--;
						if(trk_ipchi2->at(i) > 9) NADDDISPL9--;
						if(trk_ipchi2->at(i) > 16) NADDDISPL16--;
					}
				}	
				if(foundOutOfJet) ++nSVOutOfJet;
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

			tout->Fill();
		}
	}
	cout << nInJet << "\t" << nOutOfJet << "\t" << nSVOutOfJet << std::endl;
	cout << njsv / (double) njtot << endl;

	cout << njmu << " +/- " << njmuErr << "\t" << njmuFromMu << "\t" << njmu2 << endl;
	cout << njmu / (double) njtot << " +/- " << njmuErr / (double) njtot << endl;
	std::cout << countE << "\t" << countMu << "\t" << countPi << "\t" << countK << "\t" << countP << std::endl;

	fout.cd();
	tout->Write("T",TObject::kOverwrite);

	return 0;
}
