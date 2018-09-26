#include <iostream>
#include <map>
#include <set>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"

void tagEvents(TString dir="./") {
	gStyle->SetOptStat(0);

	//TFile* f = TFile::Open(dir+"for_yandex_data_SV.root");
	//TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	TChain* t = new TChain("T");
	t->Add("/tmp/dcraik/for_yandex_data_testE1k.root");
	//t->Add("/tmp/dcraik/for_yandex_data_MD.root");
	//t->Add("/tmp/dcraik/for_yandex_data_MU.root");

	TFile* fX = TFile::Open(dir+"for_yandex_data_SV_Xtag_testE1k.root","RECREATE");
	TTree* tX = t->CloneTree(0);

	//TFile* f0 = TFile::Open(dir+"for_yandex_data_SV_0tag_testE1k.root","RECREATE");
	//TTree* t0 = t->CloneTree(0);

	//TFile* f4 = TFile::Open(dir+"for_yandex_data_SV_4tag_testE1k.root","RECREATE");
	//TTree* t4 = t->CloneTree(0);

	//TFile* f5 = TFile::Open(dir+"for_yandex_data_SV_5tag_testE1k.root","RECREATE");
	//TTree* t5 = t->CloneTree(0);

	//	TFile* fU = TFile::Open(dir+"for_yandex_data_SV_notag.root","RECREATE");
	//	TTree* tU = t->CloneTree(0);

	int Evt;
	double D0M, D0IP;
	double D2K3piM, D2K3piIP;
	double DM, DIP;
	double DsM, DsIP;
	double LcM, LcIP;
	double SVMCor, SVMINIPCHI2, SVN;

	double JPX,JPY, JPZ, JE;
	double JPT;

	double JetDijetSVDec;

	double NSV;

	t->SetBranchAddress("TrueJetPx", &JPX);
	t->SetBranchAddress("TrueJetPy", &JPY);
	t->SetBranchAddress("TrueJetPz", &JPZ);
	t->SetBranchAddress("TrueJetE", &JE);
	t->SetBranchAddress("JetPT", &JPT);

	t->SetBranchAddress("D0M", &D0M);
	t->SetBranchAddress("D0IP", &D0IP);
	t->SetBranchAddress("D2K3PIM", &D2K3piM);
	t->SetBranchAddress("D2K3PIIP", &D2K3piIP);
	t->SetBranchAddress("DM", &DM);
	t->SetBranchAddress("DIP", &DIP);
	t->SetBranchAddress("DSM", &DsM);
	t->SetBranchAddress("DSIP", &DsIP);
	//	t->SetBranchAddress("LCM", &LcM);
	//	t->SetBranchAddress("LCIP", &LcIP);
	t->SetBranchAddress("EVT", &Evt);
	t->SetBranchAddress("SVMCor", &SVMCor);
	t->SetBranchAddress("SVMINIPCHI2", &SVMINIPCHI2);
	t->SetBranchAddress("SVN", &SVN);
	t->SetBranchAddress("JetDijetSVDec", &JetDijetSVDec);
	t->SetBranchAddress("NSV", &NSV);

	int n = t->GetEntries();
//	n/=30.;

	int prevEvt(0), firstInEvt(0);

	bool /*beauty(false), charm(false), light(false),*/ tagged(false); //flags to keep track of whether the current event has been tagged
	//std::set<int> foundB, foundC, foundQ; //store the events where we found beauty or charm
	std::map<int, std::vector<int>*> taggedJets;
	std::map<int, std::vector<TLorentzVector>*> taggedJetP4;

	int nPT20(0), nTag(0), nQ(0), nC(0), nB(0), nU(0), nQT(0), nCT(0), nBT(0), nUT(0);

	std::vector<int>* evtTaggedJets = new std::vector<int>();
	std::vector<TLorentzVector>* evtTaggedJetP4 = new std::vector<TLorentzVector>();

	//first loop - do the tagging
	for(int i=0; i<n; ++i) {
		if(!((int)i % (int)(n/20.))) std::cout << i << " of " << n << std::endl;
		t->GetEntry(i);

		if(JetDijetSVDec!=6 && JetDijetSVDec!=14) continue;
		//if(JPT<15000.) continue;
		++nPT20;

		if(Evt != prevEvt) {//this jet is from a new event so finish off the previous one
			if(tagged) {
				taggedJets[prevEvt] = evtTaggedJets;
				taggedJetP4[prevEvt] = evtTaggedJetP4;

				evtTaggedJets = new std::vector<int>();
				evtTaggedJetP4 = new std::vector<TLorentzVector>();

				//if(beauty) {
				//	foundB.insert(prevEvt);
				//} else if(charm) {
				//	foundC.insert(prevEvt);
				//} else if(light) {
				//	foundQ.insert(prevEvt);
				//}
			}

			//reset for next event
			//beauty=false;
			//charm=false;
			//light=false;
			tagged=false;
			prevEvt = Evt;
		}

		//check this jet for beauty or charm tag
		//if(!tagged) {
			++nTag;
			tagged=true;
			//if(SVMCor>500. && SVMCor<2500) {// && SVN>2) {
			//	charm=true;
			//	++nCT;
			//} else if(SVMCor>3000. && SVMCor<5500. && SVN>2) {
			//	beauty=true;
			//	++nBT;
			//} else if(SVMCor<500 && SVN<=2) {
			//	light = true;
			//	++nQT;
			//}
			//if(charm || beauty || light) {
			//} else {
			//	++nUT;
			//}
			evtTaggedJets->push_back(i);
			evtTaggedJetP4->push_back(TLorentzVector(JPX,JPY,JPZ,JE));
		//	} else {
		//		if(beauty) ++nB;
		//		else if(charm) ++nC;
		//		else if(light) ++nQ;
		//		else ++nU;
		//	}
		}

		//now deal with the last event
		if(tagged) {
			taggedJets[prevEvt] = evtTaggedJets;
			taggedJetP4[prevEvt] = evtTaggedJetP4;
		}
		//if(beauty) foundB.insert(prevEvt);
		//if(charm) foundC.insert(prevEvt);
		//if(light) foundQ.insert(prevEvt);

		//std::cout << foundB.size() << "\t" << foundC.size() << "\t" << foundQ.size() << "\t" << n << std::endl;
		std::cout << taggedJets.size() << "\t" << taggedJetP4.size() << std::endl;
		std::cout << nPT20 << "\t" << nTag << std::endl;
		//std::cout << nBT << "\t" << nCT << "\t" << nQT << "\t" << nUT << std::endl;
		//std::cout << nB << "\t" << nC << "\t" << nQ << "\t" << nU << std::endl;

		//histograms to plot SVMCor
		//TH1D hBT("hBT","",40,0.,10000.);
		//TH1D hCT("hCT","",40,0.,10000.);
		TH1D hQT("hQT","",40,0.,10000.);
		//TH1D hB("hB","",40,0.,10000.);
		//TH1D hC("hC","",40,0.,10000.);
		TH1D hQ("hQ","",40,0.,10000.);

		//count how many jets we find
		//int njetsB(0);
		//int njetsC(0);
		int njetsQ(0);
		//int njetsBT(0);
		//int njetsCT(0);
		int njetsQT(0);
		//int njetsBNoSV(0);
		//int njetsCNoSV(0);
		int njetsQNoSV(0);
		int njetsU(0);
		int njetsT(0);

		int nMultPair(0);
		std::map<int, std::vector<std::pair<int,int>>*> goodPairs;
		std::vector<std::pair<int,int>>* evtGoodPairs = new std::vector<std::pair<int,int>>();

		//second loop - find good probes
		for(int i=0; i<n; ++i) {
			if(!((int)i % (int)(n/20.))) std::cout << i << " of " << n << std::endl;
			t->GetEntry(i);

			//if(NSV==0) continue;
			if(JPT<10000.) continue;
			//if(JetDijetSVDec<2) continue;

			if(Evt != prevEvt) {//this jet is from a new event so finish off the previous one
				if(evtGoodPairs->size()>0) {
					goodPairs[prevEvt] = evtGoodPairs;
					evtGoodPairs = new std::vector<std::pair<int,int>>();
				}
				prevEvt=Evt;
			}

			//if(JetDijetSVDec==0 || JetDijetSVDec>15) continue;
			if(taggedJets.count(Evt)==0) {
				++njetsU;
				//			tU->Fill();
			} else {
				//if(foundB.count(Evt)) {
				//	if(taggedJet[Evt] == i) {
				//		hBT.Fill(SVMCor);
				//		++njetsBT;
				//	} else {
				//		//untagged beauty jet
				//		TLorentzVector p4(JPX,JPY,JPZ,JE);
				//		if(TMath::Abs(p4.DeltaPhi(taggedJetP4[Evt]))<2.6) {
				//			//hQ.Fill(SVMCor);
				//			continue;
				//		}
				//		if(SVMCor<0.) ++njetsBNoSV;
				//		++njetsB;
				//		hB.Fill(SVMCor);
				//		t5->Fill();
				//	}
				//} else if(foundC.count(Evt)) {
				//	if(taggedJet[Evt] == i) {
				//		hCT.Fill(SVMCor);
				//		++njetsCT;
				//	} else {
				//		//untagged charm jet
				//		TLorentzVector p4(JPX,JPY,JPZ,JE);
				//		if(TMath::Abs(p4.DeltaPhi(taggedJetP4[Evt]))<2.6) {
				//			//hQ.Fill(SVMCor);
				//			continue;
				//		}
				//		if(SVMCor<0.) ++njetsCNoSV;
				//		++njetsC;
				//		hC.Fill(SVMCor);
				//		t4->Fill();
				//	}
				//} //else if(foundQ.count(Evt)) {
				//	//event not tagged
				//
				//	Now put ALL tagged events in the "light" sample
				//if(taggedJet[Evt] == i) {
				//	hQT.Fill(SVMCor);
				//	++njetsQT;
				//} else {
				//	TLorentzVector p4(JPX,JPY,JPZ,JE);
				//	if(TMath::Abs(p4.DeltaPhi(taggedJetP4[Evt]))<2.6) {
				//		//hQ.Fill(SVMCor);
				//		continue;
				//	}
				//	if(SVMCor<0.) ++njetsQNoSV;
				//	++njetsQ;
				//	hQ.Fill(SVMCor);
				//	t0->Fill();
				//}
				//}

				TLorentzVector p4(JPX,JPY,JPZ,JE);
				if(p4.Pt()<10000. || p4.Eta()<2.5 || p4.Eta()>4.0) continue;

				evtTaggedJets = taggedJets[Evt];
				evtTaggedJetP4 = taggedJetP4[Evt];

				for(int itag=0; itag<evtTaggedJets->size(); ++itag) {
					double ptAsym = (p4.Pt()-evtTaggedJetP4->at(itag).Pt()) /
							(p4.Pt()+evtTaggedJetP4->at(itag).Pt());
					if(TMath::Abs(p4.DeltaPhi(evtTaggedJetP4->at(itag)))>2.6) {
						if(TMath::Abs(ptAsym)<0.25) {
							++njetsQ;
							evtGoodPairs->push_back(std::pair<int,int>(i,evtTaggedJets->at(itag)));
						}
					}

				}
			}
		}
		if(evtGoodPairs->size()>0) {
			if(evtGoodPairs->size() > 1) ++nMultPair;
			goodPairs[prevEvt] = evtGoodPairs;
		}


//		std::cout << njetsB << "\t" << njetsC << "\t" << njetsQ << std::endl;
//		std::cout << njetsBNoSV << "\t" << njetsCNoSV << "\t" << njetsQNoSV << std::endl;
//		std::cout << njetsBT << "\t" << njetsCT << "\t" << njetsQT << std::endl;
		std::cout << goodPairs.size() << "\t" << nMultPair << "\t" << njetsQ << "\t" << njetsU << std::endl;

		int itag(0);
		int iprobe(0);
		int whichTag(0);

		//third loop - fill histograms
		for(int i=0; i<n; ++i) {
			if(!((int)i % (int)(n/20.))) std::cout << i << " of " << n << std::endl;
			t->GetEntry(i);

			if(goodPairs.count(Evt)==0) {
				continue;
			}

			if(prevEvt!=Evt) {
				prevEvt=Evt;

				evtGoodPairs = goodPairs[Evt];
				if(evtGoodPairs->size()>1) {
					//TODO
				}
				iprobe=evtGoodPairs->at(whichTag).first;
				itag=evtGoodPairs->at(whichTag).second;
				//std::cout << evtGoodPairs->size() << "\t" << whichTag << "\t" << iprobe << "\t" << itag << std::endl;
			}

			if(iprobe!=i) continue;
			tX->Fill();
		}

//		TCanvas c;
//
//		hCT.GetXaxis()->SetTitle("SV M_{cor}");
//
//		hB.Scale(1./static_cast<double>(hB.GetEntries() - hB.GetBinContent(0) - hB.GetBinContent(hB.GetNbinsX())));
//		hC.Scale(1./static_cast<double>(hC.GetEntries() - hC.GetBinContent(0) - hC.GetBinContent(hC.GetNbinsX())));
//		hQ.Scale(1./static_cast<double>(hQ.GetEntries() - hQ.GetBinContent(0) - hQ.GetBinContent(hQ.GetNbinsX())));
//
//		hBT.Scale(0.2/static_cast<double>(hBT.GetEntries() - hBT.GetBinContent(0) - hBT.GetBinContent(hBT.GetNbinsX())));
//		hCT.Scale(0.2/static_cast<double>(hCT.GetEntries() - hCT.GetBinContent(0) - hCT.GetBinContent(hCT.GetNbinsX())));
//		hQT.Scale(0.2/static_cast<double>(hQT.GetEntries() - hQT.GetBinContent(0) - hQT.GetBinContent(hQT.GetNbinsX())));
//
//		hC.SetLineColor(kRed);
//		hQ.SetLineColor(kGreen+2);
//
//		hCT.SetLineColor(kRed);
//		hQT.SetLineColor(kGreen+2);
//
//		hQT.SetLineStyle(kDashed);
//		hCT.SetLineStyle(kDashed);
//		hBT.SetLineStyle(kDashed);
//
//		hC.Draw();
//		hQT.Draw("same");
//		hB.Draw("same");
//
//		hCT.Draw("same");
//		hBT.Draw("same");
//		hQ.Draw("same");
//
//		c.SaveAs("MCorTagged.pdf");
//
//		t0->AutoSave();
//		f0->Close();
//
//		t4->AutoSave();
//		f4->Close();
//
//		t5->AutoSave();
//		f5->Close();

		tX->AutoSave();
		fX->Close();

		//	tU->AutoSave();
		//	fU->Close();
	}

	int main() {
		tagEvents();
		return 0;
	}
