#define skimZs_cxx
#include "skimZs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>

void skimZs::Loop()
{
	if (fChain == 0) return;

	TH1D hZ("hZ","",100,40000.,120000.);
	TH1D hZsel("hZsel","",100,40000.,120000.);

	TH2D hMuIso("hMuIso","",20,40000.,120000.,20,0.,1.2);

	TH1D hD0("hD0","",80,1780.,1940.);
	TH2D hD02D("hD02D","",8,1780.,1940.,6,-2.,10.);
	TH1D hD0sel("hD0sel","",80,1780.,1940.);
	TH2D hD02Dsel("hD02Dsel","",8,1780.,1940.,6,-2.,10.);
	TH2D hD0PID("hD0PID","",5,0.,1.,5,0.,1.);
	TH2D hD0PID2("hD0PID2","",5,0.,1.,5,0.,1.);

	TH1D hDp("hDp","",80,1780.,1940.);
	TH2D hDp2D("hDp2D","",8,1780.,1940.,6,-2.,10.);
	TH1D hDpsel("hDpsel","",80,1780.,1940.);
	TH2D hDp2Dsel("hDp2Dsel","",8,1780.,1940.,6,-2.,10.);

	TH1D hDs("hDs","",80,1880.,2040.);
	TH2D hDs2D("hDs2D","",8,1880.,2040.,6,-2.,10.);
	TH1D hDssel("hDssel","",80,1880.,2040.);
	TH2D hDs2Dsel("hDs2Dsel","",8,1880.,2040.,6,-2.,10.);
	TH1D hPhi("hPhi","",100,700.,2000.);

	TFile* fout = TFile::Open("/tmp/dcraik/ZjetSkimmed.root","RECREATE");
	TTree* tout = new TTree("T","");

	double DJetPX,DJetPY,DJetPZ,DJetE,DJetPT,DJetETA;//,DJetS1,DJetS2,DJetQ,DJetN,DJetNQ,DJetNN,DJetPTD;
	double ZJet0PX,ZJet0PY,ZJet0PZ,ZJet0E,ZJet0PT,ZJet0ETA;//,ZJet0S1,ZJet0S2,ZJet0Q,ZJet0N,ZJet0NQ,ZJet0NN,ZJet0PTD;
	double ZJet1PX,ZJet1PY,ZJet1PZ,ZJet1E,ZJet1PT,ZJet1ETA;//,ZJet1S1,ZJet1S2,ZJet1Q,ZJet1N,ZJet1NQ,ZJet1NN,ZJet1PTD;
	double D0M, D0PX, D0PY, D0PZ, D0PT, D0E, D0X, D0Y, D0Z, D0FD, D0FDCHI2, D0IP, D0IPCHI2, D0LOGIPCHI2, D0VTXCHI2, D0VTXNDOF, D0TAU;
	double Z0M, Z0PX, Z0PY, Z0PZ, Z0PT, Z0E, Z0X, Z0Y, Z0Z, Z0IP, Z0IPCHI2;

	tout->Branch("DJetPx",&DJetPX);
	tout->Branch("DJetPy",&DJetPY);
	tout->Branch("DJetPz",&DJetPZ);
	tout->Branch("DJetE",&DJetE);
	tout->Branch("DJetPT",&DJetPT);
	tout->Branch("DJetEta",&DJetETA);
	//tout->Branch("DJetSigma1",&DJetS1);
	//tout->Branch("DJetSigma2",&DJetS2);
	//tout->Branch("DJetQ",&DJetQ);
	//tout->Branch("DJetMult",&DJetN);
	//tout->Branch("DJetNChr",&DJetNQ);
	//tout->Branch("DJetNNeu",&DJetNN);
	//tout->Branch("DJetPTD",&DJetPTD);

	tout->Branch("ZJet0Px",&ZJet0PX);
	tout->Branch("ZJet0Py",&ZJet0PY);
	tout->Branch("ZJet0Pz",&ZJet0PZ);
	tout->Branch("ZJet0E",&ZJet0E);
	tout->Branch("ZJet0PT",&ZJet0PT);
	tout->Branch("ZJet0Eta",&ZJet0ETA);
	//tout->Branch("ZJet0Sigma1",&ZJet0S1);
	//tout->Branch("ZJet0Sigma2",&ZJet0S2);
	//tout->Branch("ZJet0Q",&ZJet0Q);
	//tout->Branch("ZJet0Mult",&ZJet0N);
	//tout->Branch("ZJet0NChr",&ZJet0NQ);
	//tout->Branch("ZJet0NNeu",&ZJet0NN);
	//tout->Branch("ZJet0PTD",&ZJet0PTD);

	tout->Branch("ZJet1Px",&ZJet1PX);
	tout->Branch("ZJet1Py",&ZJet1PY);
	tout->Branch("ZJet1Pz",&ZJet1PZ);
	tout->Branch("ZJet1E",&ZJet1E);
	tout->Branch("ZJet1PT",&ZJet1PT);
	tout->Branch("ZJet1Eta",&ZJet1ETA);
	//tout->Branch("ZJet1Sigma1",&ZJet1S1);
	//tout->Branch("ZJet1Sigma2",&ZJet1S2);
	//tout->Branch("ZJet1Q",&ZJet1Q);
	//tout->Branch("ZJet1Mult",&ZJet1N);
	//tout->Branch("ZJet1NChr",&ZJet1NQ);
	//tout->Branch("ZJet1NNeu",&ZJet1NN);
	//tout->Branch("ZJet1PTD",&ZJet1PTD);

	tout->Branch("D0M",         &D0M);
	tout->Branch("D0PX",        &D0PX);
	tout->Branch("D0PY",        &D0PY);
	tout->Branch("D0PZ",        &D0PZ);
	tout->Branch("D0PT",        &D0PT);
	tout->Branch("D0E",         &D0E);
	tout->Branch("D0X",         &D0X);
	tout->Branch("D0Y",         &D0Y);
	tout->Branch("D0Z",         &D0Z);
	tout->Branch("D0IP",        &D0IP);
	tout->Branch("D0IPCHI2",    &D0IPCHI2);
	tout->Branch("D0LOGIPCHI2", &D0LOGIPCHI2);
	tout->Branch("D0FD",        &D0FD);
	tout->Branch("D0FDCHI2",    &D0FDCHI2);
	tout->Branch("D0TAU",       &D0TAU);
	tout->Branch("D0VTXCHI2",   &D0VTXCHI2);
	tout->Branch("D0VTXNDOF",   &D0VTXNDOF);

	tout->Branch("Z0M",        &Z0M);
	tout->Branch("Z0PX",       &Z0PX);
	tout->Branch("Z0PY",       &Z0PY);
	tout->Branch("Z0PZ",       &Z0PZ);
	tout->Branch("Z0PT",       &Z0PT);
	tout->Branch("Z0E",        &Z0E);
	tout->Branch("Z0X",        &Z0X);
	tout->Branch("Z0Y",        &Z0Y);
	tout->Branch("Z0Z",        &Z0Z);
	tout->Branch("Z0IP",       &Z0IP);
	tout->Branch("Z0IPCHI2",   &Z0IPCHI2);
	
	Long64_t nentries = fChain->GetEntries();

	int nZ(0), nZJet(0), nZD0(0), nZDp(0), nZDs(0);
	int nMultZ(0), nMultZJet(0);
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		fChain->GetEntry(jentry);

		int nZInEvent(0);

		for(int iZ=0; iZ<z0_m->size(); ++iZ) {
			if(z0_m->at(iZ)<40000. || z0_m->at(iZ)>120000.) continue;
			hZ.Fill(z0_m->at(iZ));
			int iMu0 = z0_idx_trk0->at(iZ);
			int iMu1 = z0_idx_trk1->at(iZ);
			int iJ0  = z0_idx_jet_trk0->at(iZ);
			int iJ1  = z0_idx_jet_trk1->at(iZ);
			int iPV  = z0_idx_pvr->at(iZ);
			if(iMu0<0 || iMu1<0 || iJ0<0 || iJ1<0) continue;

			double mu0_pt = TMath::Sqrt(TMath::Power(trk_px->at(iMu0),2.)+TMath::Power(trk_py->at(iMu0),2.));
			double mu1_pt = TMath::Sqrt(TMath::Power(trk_px->at(iMu1),2.)+TMath::Power(trk_py->at(iMu1),2.));
			double j0_pt(0.), j1_pt(0.);
			for(int iTrk=0; iTrk<trk_px->size(); ++iTrk) {
				if(trk_idx_jet->at(iTrk) == iJ0) j0_pt += TMath::Sqrt(TMath::Power(trk_px->at(iTrk ),2.)+TMath::Power(trk_py->at(iTrk ),2.));
				if(trk_idx_jet->at(iTrk) == iJ1) j1_pt += TMath::Sqrt(TMath::Power(trk_px->at(iTrk ),2.)+TMath::Power(trk_py->at(iTrk ),2.));
			}
			hMuIso.Fill(z0_m->at(iZ),mu0_pt/j0_pt);
			hMuIso.Fill(z0_m->at(iZ),mu1_pt/j1_pt);
			if(mu0_pt/j0_pt < 0.80) continue;
			if(mu1_pt/j1_pt < 0.80) continue;
			if(trk_idx_pvr->at(iMu0) != trk_idx_pvr->at(iMu1)) continue;
			if(trk_idx_pvr->at(iMu0) != iPV) continue;
			TLorentzVector pMu0(trk_px->at(iMu0),
					    trk_py->at(iMu0),
					    trk_pz->at(iMu0),
					    trk_e->at( iMu0));

			TLorentzVector pMu1(trk_px->at(iMu1),
					    trk_py->at(iMu1),
					    trk_pz->at(iMu1),
					    trk_e->at( iMu1));

			TLorentzVector pJMu0(jet_px->at(iJ0),
					     jet_py->at(iJ0),
					     jet_pz->at(iJ0),
					     jet_e->at( iJ0));

			TLorentzVector pJMu1(jet_px->at(iJ1),
					     jet_py->at(iJ1),
					     jet_pz->at(iJ1),
					     jet_e->at( iJ1));

			Z0M       = z0_m->at(iZ);
			Z0PX      = z0_px->at(iZ);
			Z0PY      = z0_py->at(iZ);
			Z0PZ      = z0_pz->at(iZ);
			Z0PT      = TMath::Sqrt(Z0PX*Z0PX + Z0PY*Z0PY);
			Z0E       = z0_e->at(iZ);
			Z0X       = z0_x->at(iZ);
			Z0Y       = z0_y->at(iZ);
			Z0Z       = z0_z->at(iZ);
			Z0IP      = z0_ip->at(iZ);
			Z0IPCHI2  = z0_ip_chi2->at(iZ);

			ZJet0PX = pJMu0.Px();
			ZJet0PY = pJMu0.Py();
			ZJet0PZ = pJMu0.Pz();
			ZJet0E = pJMu0.E();
			ZJet0PT = pJMu0.Pt();
			ZJet0ETA = pJMu0.Eta();

			ZJet1PX = pJMu1.Px();
			ZJet1PY = pJMu1.Py();
			ZJet1PZ = pJMu1.Pz();
			ZJet1E = pJMu1.E();
			ZJet1PT = pJMu1.Pt();
			ZJet1ETA = pJMu1.Eta();
			//ZJetSigma1;
			//ZJetSigma2;
			//ZJetQ;
			//ZJetMult;
			//ZJetNChr;
			//ZJetNNeu;
			//ZJetPTD;
		
			DJetPX = -999.;
			DJetPY = -999.;
			DJetPZ = -999.;
			DJetE = -999.;
			DJetPT = -999.;
			DJetETA = -999.;
			//DJetS1 = -999.;
			//DJetS2 = -999.;
			//DJetQ = -999.;
			//DJetN = -999.;
			//DJetNQ = -999.;
			//DJetNN = -999.;
			//DJetPTD = -999.;

			D0M = -999.;
			D0PX = -999.;
			D0PY = -999.;
			D0PZ = -999.;
			D0PT = -999.;
			D0E = -999.;
			D0X = -999.;
			D0Y = -999.;
			D0Z = -999.;
			D0FD = -999.;
			D0FDCHI2 = -999.;
			D0IP = -999.;
			D0IPCHI2 = -999.;
			D0LOGIPCHI2 = -999.;
			D0VTXCHI2 = -999.;
			D0VTXNDOF = -999.;
			D0TAU = -999.;

			++nZ;
			++nZInEvent;
			hZsel.Fill(z0_m->at(iZ));
			int nJetInZ(0);
			for(int iJet=0; iJet<jet_px->size(); ++iJet) {
				if(iJet==z0_idx_jet_trk0->at(iZ) || iJet==z0_idx_jet_trk1->at(iZ)) continue;
				int nInPV(0);
				for(int iTrk=0; iTrk<trk_px->size(); ++iTrk) {
					if(trk_idx_jet->at(iTrk) == iJet && trk_idx_pvr->at(iTrk) == iPV) ++nInPV;
				}
				if(nInPV<2) continue;
				TLorentzVector pJet(jet_px->at(iJet),
						    jet_py->at(iJet),
						    jet_pz->at(iJet),
						    jet_e->at( iJet));
				if(pJet.Pt()<15000.) continue;
				if(pMu0.DeltaR(pJet)<0.5) continue;
				if(pMu1.DeltaR(pJet)<0.5) continue;
				++nZJet;
				++nJetInZ;
				if(DJetPT < pJet.Pt() && D0M<0.) {
					DJetPX = pJet.Px();
					DJetPY = pJet.Py();
					DJetPZ = pJet.Pz();
					DJetE = pJet.E();
					DJetPT = pJet.Pt();
					DJetETA = pJet.Eta();
				}
				
				for(int iD=0; iD<d0_px->size(); ++iD) {
					if(!(trk_pnn_pi->at(d0_idx_trk1->at(iD))>0.2 && trk_pnn_k->at( d0_idx_trk0->at(iD))>0.3)) continue;

					if(d0_idx_jet->at(iD) == iJet && d0_ntrk_jet->at(iD)>0.) {
						++nZD0;
						hD0sel.Fill(d0_m->at(iD));
						hD02Dsel.Fill(d0_m->at(iD),TMath::Log(d0_ip_chi2->at(iD)));

						if(DJetPT < pJet.Pt() || D0M<0.) {
							D0M         = d0_m->at(iD);
							D0PX        = d0_px->at(iD);
							D0PY        = d0_py->at(iD);
							D0PZ        = d0_pz->at(iD);
							D0PT        = TMath::Sqrt(D0PX*D0PX + D0PY*D0PY);
							D0E         = d0_e->at(iD);
							D0X         = d0_x->at(iD);
							D0Y         = d0_y->at(iD);
							D0Z         = d0_z->at(iD);
							D0IP        = d0_ip->at(iD);
							D0IPCHI2    = d0_ip_chi2->at(iD);
							D0LOGIPCHI2 = TMath::Log(D0IPCHI2);
							D0FD        = d0_fd->at(iD);
							D0FDCHI2    = d0_fd_chi2->at(iD);
							D0TAU       = d0_tau->at(iD);
							D0VTXCHI2   = d0_vtx_chi2->at(iD);
							D0VTXNDOF   = d0_vtx_ndof->at(iD);

							DJetPX = pJet.Px();
							DJetPY = pJet.Py();
							DJetPZ = pJet.Pz();
							DJetE = pJet.E();
							DJetPT = pJet.Pt();
							DJetETA = pJet.Eta();
						}
						//DJetSigma1;
						//DJetSigma2;
						//DJetQ;
						//DJetMult;
						//DJetNChr;
						//DJetNNeu;
						//DJetPTD;
					}
				}
				for(int iD=0; iD<dp_px->size(); ++iD) {
					if(!(trk_pnn_pi->at(dp_idx_trk2->at(iD))>0.2 && trk_pnn_pi->at(dp_idx_trk1->at(iD))>0.2 && trk_pnn_k->at( dp_idx_trk0->at(iD))>0.3)) continue;

					if(dp_idx_jet->at(iD) == iJet && dp_ntrk_jet->at(iD)>1.) {
						++nZDp;
						hDpsel.Fill(dp_m->at(iD));
						hDp2Dsel.Fill(dp_m->at(iD),TMath::Log(dp_ip_chi2->at(iD)));
					}
				}
				for(int iD=0; iD<ds_px->size(); ++iD) {
					if(!(trk_pnn_pi->at(ds_idx_trk2->at(iD))>0.2 && trk_pnn_k->at(ds_idx_trk1->at(iD))>0.3 && trk_pnn_k->at( ds_idx_trk0->at(iD))>0.3)) continue;
					TLorentzVector phi(trk_px->at(ds_idx_trk0->at(iD))+trk_px->at(ds_idx_trk1->at(iD)),
							   trk_py->at(ds_idx_trk0->at(iD))+trk_py->at(ds_idx_trk1->at(iD)),
							   trk_pz->at(ds_idx_trk0->at(iD))+trk_pz->at(ds_idx_trk1->at(iD)),
							   trk_e->at( ds_idx_trk0->at(iD))+trk_e->at( ds_idx_trk1->at(iD)));
					if(phi.M()<990.||phi.M()>1050.) continue;

					if(ds_idx_jet->at(iD) == iJet && ds_ntrk_jet->at(iD)>1.) {
						++nZDs;
						hDssel.Fill(ds_m->at(iD));
						hDs2Dsel.Fill(ds_m->at(iD),TMath::Log(ds_ip_chi2->at(iD)));
					}
				}
			}
			if(nJetInZ>1) ++nMultZJet;

			if(DJetPT>0.) tout->Fill();
		}
		if(nZInEvent>1) ++nMultZ;

		for(int iD=0; iD<d0_px->size(); ++iD) {
			hD0PID.Fill(trk_pnn_pi->at(d0_idx_trk0->at(iD)), trk_pnn_k->at( d0_idx_trk0->at(iD)));
			hD0PID2.Fill(trk_pnn_pi->at(d0_idx_trk1->at(iD)), trk_pnn_k->at( d0_idx_trk1->at(iD)));
			if(!(trk_pnn_pi->at(d0_idx_trk1->at(iD))>0.2 && trk_pnn_k->at( d0_idx_trk0->at(iD))>0.3)) continue;
			if(d0_idx_jet<0) continue;
			hD0.Fill(d0_m->at(iD));
			hD02D.Fill(d0_m->at(iD),TMath::Log(d0_ip_chi2->at(iD)));
		}
		for(int iD=0; iD<dp_px->size(); ++iD) {
			if(!(trk_pnn_pi->at(dp_idx_trk2->at(iD))>0.2 && trk_pnn_pi->at(dp_idx_trk1->at(iD))>0.2 && trk_pnn_k->at( dp_idx_trk0->at(iD))>0.3)) continue;
			if(dp_idx_jet<0) continue;
			hDp.Fill(dp_m->at(iD));
			hDp2D.Fill(dp_m->at(iD),TMath::Log(dp_ip_chi2->at(iD)));
		}
		for(int iD=0; iD<ds_px->size(); ++iD) {
			if(!(trk_pnn_pi->at(ds_idx_trk2->at(iD))>0.2 && trk_pnn_k->at(ds_idx_trk1->at(iD))>0.3 && trk_pnn_k->at( ds_idx_trk0->at(iD))>0.3)) continue;
			TLorentzVector phi(trk_px->at(ds_idx_trk0->at(iD))+trk_px->at(ds_idx_trk1->at(iD)),
					   trk_py->at(ds_idx_trk0->at(iD))+trk_py->at(ds_idx_trk1->at(iD)),
					   trk_pz->at(ds_idx_trk0->at(iD))+trk_pz->at(ds_idx_trk1->at(iD)),
					   trk_e->at( ds_idx_trk0->at(iD))+trk_e->at( ds_idx_trk1->at(iD)));
			if(ds_idx_jet<0) continue;
			hPhi.Fill(phi.M());
			if(phi.M()<990.||phi.M()>1050.) continue;
			hDs.Fill(ds_m->at(iD));
			hDs2D.Fill(ds_m->at(iD),TMath::Log(ds_ip_chi2->at(iD)));
		}

		//tout->Fill();
	}
	TCanvas c;
	hZ.Draw();
	hZsel.SetLineColor(kRed);
	hZsel.Draw("same");
	c.SaveAs("mZ.pdf");

	hMuIso.Draw("colz");
	c.SaveAs("muonIsolation.pdf");

	hD0.Draw();
	hD0sel.SetLineColor(kRed);
	hD0sel.Draw("same");
	c.SaveAs("mD0.pdf");
	hD02D.Draw("colz");
	c.SaveAs("logipchi2D0_mD0_all.pdf");
	hD02Dsel.Draw("colz");
	c.SaveAs("logipchi2D0_mD0_jet.pdf");

	hDp.Draw();
	hDpsel.SetLineColor(kRed);
	hDpsel.Draw("same");
	c.SaveAs("mDp.pdf");
	hDp2D.Draw("colz");
	c.SaveAs("logipchi2Dp_mDp_all.pdf");
	hDp2Dsel.Draw("colz");
	c.SaveAs("logipchi2Dp_mDp_jet.pdf");

	hDs.Draw();
	hDssel.SetLineColor(kRed);
	hDssel.Draw("same");
	c.SaveAs("mDs.pdf");
	hDs2D.Draw("colz");
	c.SaveAs("logipchi2Ds_mDs_all.pdf");
	hDs2Dsel.Draw("colz");
	c.SaveAs("logipchi2Ds_mDs_jet.pdf");
	hPhi.Draw();
	c.SaveAs("mPhi.pdf");

	c.SetLogz();
	hD0PID.Draw("colz");
	c.SaveAs("D0K_PID.pdf");
	hD0PID2.Draw("colz");
	c.SaveAs("D0pi_PID.pdf");

	tout->AutoSave();
	fout->Close();

	std::cout << nentries << "\t" << nZ << "\t" << nZJet << "\t" << nZD0 << "\t" << nZDp << "\t" << nZDs << std::endl;
	std::cout << "\t" << nMultZ << "\t" << nMultZJet <<std::endl;
}

int main() {
	skimZs a;
	a.Loop();
	return 0;
}
