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
  TFile fout(str,"recreate");

  TTree *tout = new TTree("T","");
  int e(0);
  double JPX,JPY,JPZ,JE,JPT,JETA,JS1,JS2,JQ,JN,JNQ,JNN,JPTD;
  double PARTON,BMAXPT,CMAXPT;
  double SVM,SVMCOR,SVMINPERP,SVPT,SVDR,SVN,SVNJ,SVQ,SVFDCHI2,SVPERP,SVETA,SVTZ,
    SVMINIPCHI2,SVPX,SVPY,SVPZ,SVE,SVMAXGHOST,SVX,SVY,SVZ,SVSUMIPCHI2;
  double NSV,PVX,PVY,PVZ,NDISPL6,NDISPL9,NDISPL16;
  double MUPT,MUIPCHI2,MUDR,MUPNN,NMU;
  double HPT,HIPCHI2,HDR;
  double D0M, D0PX, D0PY, D0PZ, D0E, D0X, D0Y, D0Z, D0FD, D0DIRA, D0DOCA, D0IPCHI2MIN, D0DOCAKPI, D0VTXCHI2;
  double DPMM, DPMPX, DPMPY, DPMPZ, DPME, DPMX, DPMY, DPMZ, DPMFD, DPMDIRA, DPMDOCA, DPMIPCHI2MIN, DPMDOCAMAX;
  double DSM, DSPX, DSPY, DSPZ, DSE, DSX, DSY, DSZ, DSFD, DSDIRA, DSDOCA, DSIPCHI2MIN, DSDOCAMAX;
  double LCM, LCPX, LCPY, LCPZ, LCE, LCX, LCY, LCZ, LCFD, LCDIRA, LCDOCA, LCIPCHI2MIN, LCDOCAMAX;
  double D2K3PIM, D2K3PIPX, D2K3PIPY, D2K3PIPZ, D2K3PIE, D2K3PIX, D2K3PIY, D2K3PIZ, D2K3PIFD, D2K3PIDIRA, D2K3PIDOCA, D2K3PIIPCHI2MIN, D2K3PIDOCAMAX;

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
  tout->Branch("D0FD",       &D0FD);
  tout->Branch("D0DIRA",     &D0DIRA);
  tout->Branch("D0IP",       &D0DOCA);
  tout->Branch("D0IPCHI2MIN",    &D0IPCHI2MIN);
  tout->Branch("D0DOCAKPI",  &D0DOCAKPI);
  tout->Branch("D0VTXCHI2",  &D0VTXCHI2);
  tout->Branch("DM",        &DPMM);
  tout->Branch("DPX",       &DPMPX);
  tout->Branch("DPY",       &DPMPY);
  tout->Branch("DPZ",       &DPMPZ);
  tout->Branch("DE",        &DPME);
  tout->Branch("DX",        &DPMX);
  tout->Branch("DY",        &DPMY);
  tout->Branch("DZ",        &DPMZ);
  tout->Branch("DFD",       &DPMFD);
  tout->Branch("DDIRA",     &DPMDIRA);
  tout->Branch("DIP",       &DPMDOCA);
  tout->Branch("DIPCHI2MIN",    &DPMIPCHI2MIN);
  tout->Branch("DDOCAMAX",  &DPMDOCAMAX);
  tout->Branch("DSM",        &DSM);
  tout->Branch("DSPX",       &DSPX);
  tout->Branch("DSPY",       &DSPY);
  tout->Branch("DSPZ",       &DSPZ);
  tout->Branch("DSE",        &DSE);
  tout->Branch("DSX",        &DSX);
  tout->Branch("DSY",        &DSY);
  tout->Branch("DSZ",        &DSZ);
  tout->Branch("DSFD",       &DSFD);
  tout->Branch("DSDIRA",     &DSDIRA);
  tout->Branch("DSIP",       &DSDOCA);
  tout->Branch("DSIPCHI2MIN",    &DSIPCHI2MIN);
  tout->Branch("DSDOCAMAX",  &DSDOCAMAX);
  tout->Branch("LCM",        &LCM);
  tout->Branch("LCPX",       &LCPX);
  tout->Branch("LCPY",       &LCPY);
  tout->Branch("LCPZ",       &LCPZ);
  tout->Branch("LCE",        &LCE);
  tout->Branch("LCX",        &LCX);
  tout->Branch("LCY",        &LCY);
  tout->Branch("LCZ",        &LCZ);
  tout->Branch("LCFD",       &LCFD);
  tout->Branch("LCDIRA",     &LCDIRA);
  tout->Branch("LCIP",       &LCDOCA);
  tout->Branch("LCIPCHI2MIN",    &LCIPCHI2MIN);
  tout->Branch("LCDOCAMAX",  &LCDOCAMAX);
  tout->Branch("D2K3PIM",        &D2K3PIM);
  tout->Branch("D2K3PIPX",       &D2K3PIPX);
  tout->Branch("D2K3PIPY",       &D2K3PIPY);
  tout->Branch("D2K3PIPZ",       &D2K3PIPZ);
  tout->Branch("D2K3PIE",        &D2K3PIE);
  tout->Branch("D2K3PIX",        &D2K3PIX);
  tout->Branch("D2K3PIY",        &D2K3PIY);
  tout->Branch("D2K3PIZ",        &D2K3PIZ);
  tout->Branch("D2K3PIFD",       &D2K3PIFD);
  tout->Branch("D2K3PIDIRA",     &D2K3PIDIRA);
  tout->Branch("D2K3PIIP",       &D2K3PIDOCA);
  tout->Branch("D2K3PIIPCHI2MIN",    &D2K3PIIPCHI2MIN);
  tout->Branch("D2K3PIDOCAMAX",  &D2K3PIDOCAMAX);

  TChain *t = new TChain("data");
  for(int i=1; i<101; ++i) {
	  t->Add(TString::Format("/tmp/dcraik/lightjets_filtered_addVars_resampled%d.root",i));
  }
  //t->Add("/tmp/dcraik/jets_sim_mix.root");
  //t->Add("/tmp/dcraik/light.root");
  //t->Add("lightjets_filtered_addVars_resampled1.root");
  //for(int i=1; i<=4; i++){
  //  sprintf(str,"/Users/philten/data/Run2Jets.MC15.MU.490000%d%d.0.160925.00.root",type,i);
  //  sprintf(str,"/Users/philten/data/Run2Jets.MC15.MD.490000%d%d.0.160925.00.root",type,i);
  //  t->Add(str);
  //}

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
  //nent = 1e5;
  cout << nent << endl;

  int njtot = 0, njsv = 0;
  for(e=0; e<nent; e++){

    t->GetEntry(e);
    if(npv != 1) continue;

    int ng = gen_pid->size();
    if(ng < 2) continue;
    int nj = jet_e->size();
    int nsv = svz->size();
    int ntrk = trk_e->size();
    int nneu = neu_e->size();

    for(int j=0; j<nj; j++){
      TLorentzVector p4j(jet_px->at(j),jet_py->at(j),jet_pz->at(j),
			 jet_e->at(j));
      if(p4j.Pt() < 20e3) continue;
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

      double D0ptmax(0), Dptmax(0), Dsptmax(0), Lcptmax(0), D2K3piptmax(0);
      int bestD0(0), bestD(0), bestDs(0), bestLc(0), bestD2K3pi(0);

      int countD0(0);
      vector<double> d0_m;
      vector<double> d0_px;
      vector<double> d0_py;
      vector<double> d0_pz;
      vector<double> d0_e;
      vector<double> d0_x;
      vector<double> d0_y;
      vector<double> d0_z;
      vector<double> d0_fd;
      vector<double> d0_dira;
      vector<double> d0_doca;
      vector<double> d0_ipmin;
      vector<double> d0_docaKpi;
      vector<double> d0_vtxchi2;

      int countDpm(0);
      vector<double> d_m;
      vector<double> d_px;
      vector<double> d_py;
      vector<double> d_pz;
      vector<double> d_e;
      vector<double> d_x;
      vector<double> d_y;
      vector<double> d_z;
      vector<double> d_fd;
      vector<double> d_dira;
      vector<double> d_doca;
      vector<double> d_ipmin;
      vector<double> d_docamax;

      int countDs(0);
      vector<double> ds_m;
      vector<double> ds_px;
      vector<double> ds_py;
      vector<double> ds_pz;
      vector<double> ds_e;
      vector<double> ds_x;
      vector<double> ds_y;
      vector<double> ds_z;
      vector<double> ds_fd;
      vector<double> ds_dira;
      vector<double> ds_doca;
      vector<double> ds_ipmin;
      vector<double> ds_docamax;

      int countLc(0);
      vector<double> lc_m;
      vector<double> lc_px;
      vector<double> lc_py;
      vector<double> lc_pz;
      vector<double> lc_e;
      vector<double> lc_x;
      vector<double> lc_y;
      vector<double> lc_z;
      vector<double> lc_fd;
      vector<double> lc_dira;
      vector<double> lc_doca;
      vector<double> lc_ipmin;
      vector<double> lc_docamax;

      int countD2K3pi(0);
      vector<double> d2k3pi_m;
      vector<double> d2k3pi_px;
      vector<double> d2k3pi_py;
      vector<double> d2k3pi_pz;
      vector<double> d2k3pi_e;
      vector<double> d2k3pi_x;
      vector<double> d2k3pi_y;
      vector<double> d2k3pi_z;
      vector<double> d2k3pi_fd;
      vector<double> d2k3pi_dira;
      vector<double> d2k3pi_doca;
      vector<double> d2k3pi_ipmin;
      vector<double> d2k3pi_docamax;

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
	//all of the combinations we're going to try have opp charged Kpi pair so start with that
	if(trk_pnnk->at(i)>0.3 /*&& TMath::Abs(trk_pid->at(i))==321*/ && trk_ipchi2->at(i)>16. && trk_type->at(i)==3) {
	  	TVector3 xtrk1 = TVector3(trk_x->at(i),trk_y->at(i),trk_z->at(i));
	  	TLorentzVector p4trk1(trk_px->at(i),trk_py->at(i),trk_pz->at(i),trk_e->at(i));
		p4trk1.SetE(TMath::Sqrt(p4trk1.P()*p4trk1.P() + 493.7*493.7));
		//make D0 candidates
		for(int ii=0; ii<ntrk; ++ii) {//try every pair twice as first picked is the kaon
			if(trk_j->at(ii) != j) continue;
			if(ii==i) continue;
			if(trk_pid->at(i)*trk_pid->at(ii) > 0) continue;
			if(trk_pnnpi->at(ii)>0.3 /*&& TMath::Abs(trk_pid->at(ii))==211*/ && trk_ipchi2->at(ii)>16. && trk_type->at(ii)==3) {
	  			TVector3 xtrk2 = TVector3(trk_x->at(ii),trk_y->at(ii),trk_z->at(ii));
	  			TLorentzVector p4trk2(trk_px->at(ii),trk_py->at(ii),trk_pz->at(ii),trk_e->at(ii));
				p4trk2.SetE(TMath::Sqrt(p4trk2.P()*p4trk2.P() + 139.6*139.6));

				//also try to make D+ candidates
				for(int iii=ii+1; iii<ntrk; ++iii) {
					if(trk_j->at(iii) != j) continue;
					if(iii==i) continue;
					if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//same charge pions
					if(trk_pnnpi->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==211*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
	  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
	  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
						p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 139.6*139.6));
						TLorentzVector p4Dpm = p4trk1 + p4trk2 + p4trk3;
						if(TMath::Abs(p4Dpm.M()-1870.) > 160.) continue;
						TVector3 sv12, sv13, sv23, sv123;
						double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
						double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
						double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
						double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
						double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
						if(docamax>0.1) continue;
						sv123 = sv12 + sv13 + sv23;
						sv123 *= (1./3.);
						double docaDpm = calcDocaPoint(sv123, p4Dpm.Vect(), pv);
						TVector3 fvDpm = sv123-pv;
						double diraDpm = fvDpm.Unit().Dot(p4Dpm.Vect().Unit());
						//printf("found D+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Dpm.M(), docamax, diraDpm);
						++countDpm;
						d_m.push_back(p4Dpm.M());
						d_px.push_back(p4Dpm.Px());
  						d_py.push_back(p4Dpm.Py());
  						d_pz.push_back(p4Dpm.Pz());
  						d_e.push_back(p4Dpm.E());
  						d_x.push_back(sv123.X());
  						d_y.push_back(sv123.Y());
  						d_z.push_back(sv123.Z());
  						d_fd.push_back(fvDpm.Mag());
  						d_dira.push_back(diraDpm);
  						d_doca.push_back(docaDpm);
  						d_ipmin.push_back(ipmin);
  						d_docamax.push_back(docamax);

						if(p4Dpm.Pt() > Dptmax) {
							Dptmax = p4Dpm.Pt();
							bestD = d_m.size()-1;
						}
					}
				}

				//also try to make Ds+ candidates
				for(int iii=0; iii<ntrk; ++iii) {
					if(trk_j->at(iii) != j) continue;
					if(iii==i || iii==ii) continue;
					if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//opp charged kaons
					if(trk_pnnk->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
	  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
	  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
						p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 493.7*493.7));
						TLorentzVector p4Ds = p4trk1 + p4trk2 + p4trk3;
						if(TMath::Abs(p4Ds.M()-1968.) > 160.) continue;
						TVector3 sv12, sv13, sv23, sv123;
						double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
						double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
						double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
						double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
						double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
						if(docamax>0.1) continue;
						sv123 = sv12 + sv13 + sv23;
						sv123 *= (1./3.);
						double docaDs = calcDocaPoint(sv123, p4Ds.Vect(), pv);
						TVector3 fvDs = sv123-pv;
						double diraDs = fvDs.Unit().Dot(p4Ds.Vect().Unit());
						//printf("found Ds+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Ds.M(), docamax, diraDs);
						++countDs;
						ds_m.push_back(p4Ds.M());
						ds_px.push_back(p4Ds.Px());
  						ds_py.push_back(p4Ds.Py());
  						ds_pz.push_back(p4Ds.Pz());
  						ds_e.push_back(p4Ds.E());
  						ds_x.push_back(sv123.X());
  						ds_y.push_back(sv123.Y());
  						ds_z.push_back(sv123.Z());
  						ds_fd.push_back(fvDs.Mag());
  						ds_dira.push_back(diraDs);
  						ds_doca.push_back(docaDs);
  						ds_ipmin.push_back(ipmin);
  						ds_docamax.push_back(docamax);

						if(p4Ds.Pt() > Dsptmax) {
							Dsptmax = p4Ds.Pt();
							bestDs = ds_m.size()-1;
						}
					}
				}

				//also try to make Lc+ candidates
				for(int iii=0; iii<ntrk; ++iii) {
					if(trk_j->at(iii) != j) continue;
					if(iii==i || iii==ii) continue;
					if(trk_pid->at(i)*trk_pid->at(iii) > 0) continue;//opp charged kaon and proton
					if(trk_pnnp->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==2212*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
	  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
	  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
						p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 938.3*938.3));
						TLorentzVector p4Lc = p4trk1 + p4trk2 + p4trk3;
						if(TMath::Abs(p4Lc.M()-2286.) > 160.) continue;
						TVector3 sv12, sv13, sv23, sv123;
						double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
						double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
						double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
						double docamax = TMath::Max(doca12, TMath::Max(doca13,doca23));
						double ipmin = TMath::Min(trk_ipchi2->at(i), TMath::Min(trk_ipchi2->at(ii),trk_ipchi2->at(iii)));
						if(docamax>0.1) continue;
						sv123 = sv12 + sv13 + sv23;
						sv123 *= (1./3.);
						double docaLc = calcDocaPoint(sv123, p4Lc.Vect(), pv);
						TVector3 fvLc = sv123-pv;
						double diraLc = fvLc.Unit().Dot(p4Lc.Vect().Unit());
						//printf("found Lc+: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4Lc.M(), docamax, diraLc);
						++countLc;
						lc_m.push_back(p4Lc.M());
						lc_px.push_back(p4Lc.Px());
  						lc_py.push_back(p4Lc.Py());
  						lc_pz.push_back(p4Lc.Pz());
  						lc_e.push_back(p4Lc.E());
  						lc_x.push_back(sv123.X());
  						lc_y.push_back(sv123.Y());
  						lc_z.push_back(sv123.Z());
  						lc_fd.push_back(fvLc.Mag());
  						lc_dira.push_back(diraLc);
  						lc_doca.push_back(docaLc);
  						lc_ipmin.push_back(ipmin);
  						lc_docamax.push_back(docamax);

						if(p4Lc.Pt() > Lcptmax) {
							Lcptmax = p4Lc.Pt();
							bestLc = lc_m.size()-1;
						}
					}
				}

				//also try to make D0->Kpipipi candidates
				for(int iii=0; iii<ntrk; ++iii) {
					if(trk_j->at(iii) != j) continue;
					if(iii==i || iii==ii) continue;
					if(trk_pnnpi->at(iii)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iii)>16. && trk_type->at(iii)==3) {
	  					TVector3 xtrk3 = TVector3(trk_x->at(iii),trk_y->at(iii),trk_z->at(iii));
	  					TLorentzVector p4trk3(trk_px->at(iii),trk_py->at(iii),trk_pz->at(iii),trk_e->at(iii));
						p4trk3.SetE(TMath::Sqrt(p4trk3.P()*p4trk3.P() + 139.6*139.6));
						for(int iv=iii+1; iv<ntrk; ++iv) {
							if(trk_j->at(iv) != j) continue;
							if(iv==i || iv==ii) continue;
							if(trk_pid->at(iv)*trk_pid->at(iii) > 0) continue;//adding an opp charged pion pair
							if(trk_pnnpi->at(iv)>0.3 /*&& TMath::Abs(trk_pid->at(iii))==321*/ && trk_ipchi2->at(iv)>16. && trk_type->at(iv)==3) {
	  							TVector3 xtrk4 = TVector3(trk_x->at(iv),trk_y->at(iv),trk_z->at(iv));
	  							TLorentzVector p4trk4(trk_px->at(iv),trk_py->at(iv),trk_pz->at(iv),trk_e->at(iv));
								p4trk4.SetE(TMath::Sqrt(p4trk4.P()*p4trk4.P() + 139.6*139.6));
								TLorentzVector p4D2K3pi = p4trk1 + p4trk2 + p4trk3 + p4trk4;
								if(TMath::Abs(p4D2K3pi.M()-1864.) > 160.) continue;
								TVector3 sv12, sv13, sv14, sv23, sv24, sv34, sv1234;
								double doca12 = calcDoca(sv12, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
								double doca13 = calcDoca(sv13, xtrk1, p4trk1.Vect(), xtrk3, p4trk3.Vect());
								double doca14 = calcDoca(sv14, xtrk1, p4trk1.Vect(), xtrk4, p4trk4.Vect());
								double doca23 = calcDoca(sv23, xtrk2, p4trk2.Vect(), xtrk3, p4trk3.Vect());
								double doca24 = calcDoca(sv24, xtrk2, p4trk2.Vect(), xtrk4, p4trk4.Vect());
								double doca34 = calcDoca(sv34, xtrk3, p4trk3.Vect(), xtrk4, p4trk4.Vect());
								double docamax = TMath::Max(TMath::Max(doca12, TMath::Max(doca13,doca14)),TMath::Max(doca23, TMath::Max(doca24,doca34)));
								double ipmin = TMath::Min(TMath::Min(trk_ipchi2->at(i),trk_ipchi2->at(ii)), TMath::Min(trk_ipchi2->at(iii),trk_ipchi2->at(iv)));
								if(docamax>0.1) continue;
								sv1234 = sv12 + sv13 + sv14 + sv23 + sv24 + sv34;
								sv1234 *= (1./6.);
								double docaD2K3pi = calcDocaPoint(sv1234, p4D2K3pi.Vect(), pv);
								TVector3 fvD2K3pi = sv1234-pv;
								double diraD2K3pi = fvD2K3pi.Unit().Dot(p4D2K3pi.Vect().Unit());
								//printf("found D0->K3pi: mass=%.2f, docamax=%.3f, dira=%.4f\n", p4D2K3pi.M(), docamax, diraD2K3pi);
								++countD2K3pi;
								d2k3pi_m.push_back(p4D2K3pi.M());
								d2k3pi_px.push_back(p4D2K3pi.Px());
  								d2k3pi_py.push_back(p4D2K3pi.Py());
  								d2k3pi_pz.push_back(p4D2K3pi.Pz());
  								d2k3pi_e.push_back(p4D2K3pi.E());
  								d2k3pi_x.push_back(sv1234.X());
  								d2k3pi_y.push_back(sv1234.Y());
  								d2k3pi_z.push_back(sv1234.Z());
  								d2k3pi_fd.push_back(fvD2K3pi.Mag());
  								d2k3pi_dira.push_back(diraD2K3pi);
  								d2k3pi_doca.push_back(docaD2K3pi);
  								d2k3pi_ipmin.push_back(ipmin);
  								d2k3pi_docamax.push_back(docamax);

								if(p4D2K3pi.Pt() > D2K3piptmax) {
									D2K3piptmax = p4D2K3pi.Pt();
									bestD2K3pi = d2k3pi_m.size()-1;
								}
							}
						}
					}
				}

				//D0 candidates
				TLorentzVector p4D0 = p4trk1 + p4trk2;
				if(TMath::Abs(p4D0.M()-1864.) > 160.) continue;
				TVector3 sv;
				double docaKpi = calcDoca(sv, xtrk1, p4trk1.Vect(), xtrk2, p4trk2.Vect());
				double ipmin = TMath::Min(trk_ipchi2->at(i), trk_ipchi2->at(ii));
				if(docaKpi>0.1) continue;
				double sigma2i = trk_ip->at(i)*trk_ip->at(i) / trk_ipchi2->at(i);
				double sigma2j = trk_ip->at(ii)*trk_ip->at(ii) / trk_ipchi2->at(ii);
				double vtxchi2 = docaKpi*docaKpi/(sigma2i+sigma2j);
				if(vtxchi2>10.) continue;
				double docaD0 = calcDocaPoint(sv, p4D0.Vect(), pv);
				TVector3 fv = sv-pv;
				double dira = fv.Unit().Dot(p4D0.Vect().Unit());
				//printf("found D0: mass=%.2f, doca=%.3f, dira=%.4f, vtxchi2=%.1f, KID=%.0f, pNNk=%.2f, piID=%.0f, pNNpi=%.2f\n", p4D0.M(), docaKpi, dira, vtxchi2, trk_pid->at(i), trk_pnnk->at(i), trk_pid->at(ii), trk_pnnpi->at(ii));
				++countD0;
				d0_m.push_back(p4D0.M());
				d0_px.push_back(p4D0.Px());
  				d0_py.push_back(p4D0.Py());
  				d0_pz.push_back(p4D0.Pz());
  				d0_e.push_back(p4D0.E());
  				d0_x.push_back(sv.X());
  				d0_y.push_back(sv.Y());
  				d0_z.push_back(sv.Z());
  				d0_fd.push_back(fv.Mag());
  				d0_dira.push_back(dira);
  				d0_doca.push_back(docaD0);
  				d0_ipmin.push_back(ipmin);
  				d0_docaKpi.push_back(docaKpi);
  				d0_vtxchi2.push_back(vtxchi2);

				if(p4D0.Pt() > D0ptmax) {
					D0ptmax = p4D0.Pt();
					bestD0 = d0_m.size()-1;
				}
			}
		}
	}
      }
      //pick a random D0
      if(!d0_px.empty()) {
        int whichD0 = bestD0;//rand.Integer(d0_px.size());
        D0M       = d0_m[whichD0];
        D0PX      = d0_px[whichD0];
        D0PY      = d0_py[whichD0];
        D0PZ      = d0_pz[whichD0];
        D0E       = d0_e[whichD0];
        D0X       = d0_x[whichD0];
        D0Y       = d0_y[whichD0];
        D0Z       = d0_z[whichD0];
        D0FD      = d0_fd[whichD0];
        D0DIRA    = d0_dira[whichD0];
        D0DOCA    = d0_doca[whichD0];
        D0IPCHI2MIN= d0_ipmin[whichD0];
        D0DOCAKPI = d0_docaKpi[whichD0];
        D0VTXCHI2 = d0_vtxchi2[whichD0];
	//if(D0PX==0.) std::cout << D0M << "\t" << D0PX << "\t" << D0PY << "\t" << D0PZ << std::endl;
      } else {
        D0M       = -1000.;
        D0PX      = -1000.;
        D0PY      = -1000.;
        D0PZ      = -1000.;
        D0E       = -1000.;
        D0X       = -1000.;
        D0Y       = -1000.;
        D0Z       = -1000.;
        D0FD      = -1000.;
        D0DIRA    = -1000.;
        D0DOCA    = -1000.;
        D0IPCHI2MIN= -1000;
        D0DOCAKPI = -1000.;
        D0VTXCHI2 = -1000.;
      }
      //pick a random D+
      if(!d_px.empty()) {
        int whichD = bestD;//rand.Integer(d_px.size());
        DPMM       = d_m[whichD];
        DPMPX      = d_px[whichD];
        DPMPY      = d_py[whichD];
        DPMPZ      = d_pz[whichD];
        DPME       = d_e[whichD];
        DPMX       = d_x[whichD];
        DPMY       = d_y[whichD];
        DPMZ       = d_z[whichD];
        DPMFD      = d_fd[whichD];
        DPMDIRA    = d_dira[whichD];
        DPMIPCHI2MIN= d_ipmin[whichD];
        DPMDOCA    = d_doca[whichD];
        DPMDOCAMAX = d_docamax[whichD];
      } else {
        DPMM       = -1000.;
        DPMPX      = -1000.;
        DPMPY      = -1000.;
        DPMPZ      = -1000.;
        DPME       = -1000.;
        DPMX       = -1000.;
        DPMY       = -1000.;
        DPMZ       = -1000.;
        DPMFD      = -1000.;
        DPMDIRA    = -1000.;
        DPMIPCHI2MIN= -1000;
        DPMDOCA    = -1000.;
        DPMDOCAMAX = -1000.;
      }
      //pick a random Ds+
      if(!ds_px.empty()) {
        int whichDs = bestDs;//rand.Integer(ds_px.size());
        DSM       = ds_m[whichDs];
        DSPX      = ds_px[whichDs];
        DSPY      = ds_py[whichDs];
        DSPZ      = ds_pz[whichDs];
        DSE       = ds_e[whichDs];
        DSX       = ds_x[whichDs];
        DSY       = ds_y[whichDs];
        DSZ       = ds_z[whichDs];
        DSFD      = ds_fd[whichDs];
        DSDIRA    = ds_dira[whichDs];
        DSIPCHI2MIN= ds_ipmin[whichDs];
        DSDOCA    = ds_doca[whichDs];
        DSDOCAMAX = ds_docamax[whichDs];
      } else {
        DSM       = -1000.;
        DSPX      = -1000.;
        DSPY      = -1000.;
        DSPZ      = -1000.;
        DSE       = -1000.;
        DSX       = -1000.;
        DSY       = -1000.;
        DSZ       = -1000.;
        DSFD      = -1000.;
        DSDIRA    = -1000.;
        DSIPCHI2MIN= -1000;
        DSDOCA    = -1000.;
        DSDOCAMAX = -1000.;
      }
      //pick a random Lc+
      if(!lc_px.empty()) {
        int whichLc = bestLc;//rand.Integer(lc_px.size());
        LCM       = lc_m[whichLc];
        LCPX      = lc_px[whichLc];
        LCPY      = lc_py[whichLc];
        LCPZ      = lc_pz[whichLc];
        LCE       = lc_e[whichLc];
        LCX       = lc_x[whichLc];
        LCY       = lc_y[whichLc];
        LCZ       = lc_z[whichLc];
        LCFD      = lc_fd[whichLc];
        LCDIRA    = lc_dira[whichLc];
        LCIPCHI2MIN= lc_ipmin[whichLc];
        LCDOCA    = lc_doca[whichLc];
        LCDOCAMAX = lc_docamax[whichLc];
      } else {
        LCM       = -1000.;
        LCPX      = -1000.;
        LCPY      = -1000.;
        LCPZ      = -1000.;
        LCE       = -1000.;
        LCX       = -1000.;
        LCY       = -1000.;
        LCZ       = -1000.;
        LCFD      = -1000.;
        LCDIRA    = -1000.;
        LCIPCHI2MIN= -1000;
        LCDOCA    = -1000.;
        LCDOCAMAX = -1000.;
      }
      //pick a random D0->K3pi
      if(!d2k3pi_px.empty()) {
        int whichD0 = bestD2K3pi;//rand.Integer(d2k3pi_px.size());
        D2K3PIM       = d2k3pi_m[whichD0];
        D2K3PIPX      = d2k3pi_px[whichD0];
        D2K3PIPY      = d2k3pi_py[whichD0];
        D2K3PIPZ      = d2k3pi_pz[whichD0];
        D2K3PIE       = d2k3pi_e[whichD0];
        D2K3PIX       = d2k3pi_x[whichD0];
        D2K3PIY       = d2k3pi_y[whichD0];
        D2K3PIZ       = d2k3pi_z[whichD0];
        D2K3PIFD      = d2k3pi_fd[whichD0];
        D2K3PIDIRA    = d2k3pi_dira[whichD0];
        D2K3PIIPCHI2MIN= d2k3pi_ipmin[whichD0];
        D2K3PIDOCA    = d2k3pi_doca[whichD0];
        D2K3PIDOCAMAX = d2k3pi_docamax[whichD0];
      } else {
        D2K3PIM       = -1000.;
        D2K3PIPX      = -1000.;
        D2K3PIPY      = -1000.;
        D2K3PIPZ      = -1000.;
        D2K3PIE       = -1000.;
        D2K3PIX       = -1000.;
        D2K3PIY       = -1000.;
        D2K3PIZ       = -1000.;
        D2K3PIFD      = -1000.;
        D2K3PIDIRA    = -1000.;
        D2K3PIIPCHI2MIN= -1000;
        D2K3PIDOCA    = -1000.;
        D2K3PIDOCAMAX = -1000.;
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

      tout->Fill();
    }
  }
  cout << njsv / (double) njtot << endl;

  fout.cd();
  tout->Write("T",TObject::kOverwrite);

  return 0;
}
