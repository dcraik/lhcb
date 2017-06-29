#include <iostream>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVectorT.h"
#include "TMatrixD.h"
#include "TRandom.h"

using namespace std;

int main(int argc, char *argv[]){

  int type = atoi(argv[1]);

  char str[1000];
  sprintf(str,"for_yandex_%d.root",type);
  TFile fout(str,"recreate");

  TTree *tout = new TTree("T","");
  double JPX,JPY,JPZ,JE,JPT,JETA,JS1,JS2,JQ,JN,JNQ,JNN,JPTD;
  double PARTON,BMAXPT,CMAXPT;
  double SVM,SVMCOR,SVMINPERP,SVPT,SVDR,SVN,SVNJ,SVQ,SVFDCHI2,SVPERP,SVETA,SVTZ,
    SVMINIPCHI2,SVPX,SVPY,SVPZ,SVE,SVMAXGHOST,SVX,SVY,SVZ,SVSUMIPCHI2;
  double NSV,PVX,PVY,PVZ,NDISPL6,NDISPL9,NDISPL16;
  double MUPT,MUIPCHI2,MUDR,MUPNN,NMU;
  double HPT,HIPCHI2,HDR;

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

  TChain *t = new TChain("data");
  for(int i=1; i<=4; i++){
    sprintf(str,"/Users/philten/data/Run2Jets.MC15.MU.490000%d%d.0.160925.00.root",type,i);
    sprintf(str,"/Users/philten/data/Run2Jets.MC15.MD.490000%d%d.0.160925.00.root",type,i);
    t->Add(str);
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
  vector<double> *trk_ipchi2 = new vector<double>();
  vector<double> *trk_pnnmu = new vector<double>();
  vector<double> *trk_ismu = new vector<double>();
  vector<double> *trk_j = new vector<double>();
  t->SetBranchAddress("trk_q",&trk_q);  
  t->SetBranchAddress("trk_prb_ghost",&trk_ghost);  
  t->SetBranchAddress("trk_x",&trk_x);  
  t->SetBranchAddress("trk_y",&trk_y);  
  t->SetBranchAddress("trk_z",&trk_z);  
  t->SetBranchAddress("trk_px",&trk_px);  
  t->SetBranchAddress("trk_py",&trk_py);  
  t->SetBranchAddress("trk_pz",&trk_pz);  
  t->SetBranchAddress("trk_e",&trk_e);    
  t->SetBranchAddress("trk_ip_chi2",&trk_ipchi2);    
  t->SetBranchAddress("trk_pnn_mu",&trk_pnnmu);    
  t->SetBranchAddress("trk_is_mu",&trk_ismu);    
  t->SetBranchAddress("trk_idx_jet",&trk_j);    

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
  for(int e=0; e<nent; e++){

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
