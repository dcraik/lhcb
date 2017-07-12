{
	gStyle->SetOptStat(0);

	TFile* f = TFile::Open("for_yandex_data_SV.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

	TFile* f0 = TFile::Open("for_yandex_data_SV_0tag.root","RECREATE");
	TTree* t0 = t->CloneTree(0);

	TFile* f4 = TFile::Open("for_yandex_data_SV_4tag.root","RECREATE");
	TTree* t4 = t->CloneTree(0);

	TFile* f5 = TFile::Open("for_yandex_data_SV_4tag.root","RECREATE");
	TTree* t5 = t->CloneTree(0);

	int Evt;
	double D0M, D0IP;
	double D2K3piM, D2K3piIP;
	double DM, DIP;
	double DsM, DsIP;
	double LcM, LcIP;
	double SVMCor;

	double JPX,JPY, JPZ;

	double JetDijetSVDec;

	t->SetBranchAddress("TrueJetPx", &JPX);
	t->SetBranchAddress("TrueJetPy", &JPY);
	t->SetBranchAddress("TrueJetPz", &JPZ);

	t->SetBranchAddress("D0M", &D0M);
	t->SetBranchAddress("D0IP", &D0IP);
	t->SetBranchAddress("D2K3PIM", &D2K3piM);
	t->SetBranchAddress("D2K3PIIP", &D2K3piIP);
	t->SetBranchAddress("DM", &DM);
	t->SetBranchAddress("DIP", &DIP);
	t->SetBranchAddress("DSM", &DsM);
	t->SetBranchAddress("DSIP", &DsIP);
	t->SetBranchAddress("LCM", &LcM);
	t->SetBranchAddress("LCIP", &LcIP);
	t->SetBranchAddress("EVT", &Evt);
	t->SetBranchAddress("SVMCor", &SVMCor);
	t->SetBranchAddress("JetDijetSVDec", &JetDijetSVDec);

	int n = t->GetEntries();

	int prevEvt(0), firstInEvt(0);

	bool beauty(false), charm(false); //flags to keep track of whether the current event has been tagged
	std::set<int> foundB, foundC; //store the events where we found beauty or charm

	//first loop - do the tagging
	for(int i=0; i<n; ++i) {
		if(!((int)i % (int)(n/20.))) std::cout << i << " of " << n << std::endl;
		t->GetEntry(i);

		if(JetDijetSVDec!=6 && JetDijetSVDec!=14) continue;

		if(Evt != prevEvt) {//this jet is from a new event so finish off the previous one
			if(beauty) {
				foundB.insert(prevEvt);
			} else if(charm) {
				foundC.insert(prevEvt);
			}

			//reset for next event
			beauty=false;
			charm=false;
			prevEvt = Evt;
			whichJetTagged = -1;
		}

		//check this jet for beauty or charm tag
		if(!beauty && !charm) {
			if(TMath::Abs(D0M-1864.)<25.) {
				if(D0IP>0.05) beauty=true;
				else charm=true;
			} else if(TMath::Abs(D2K3piM-1864.)<25.) {
				if(D2K3piIP>0.05) beauty=true;
				else charm=true;
			} else if(TMath::Abs(DM-1870.)<25.) {
				if(DIP>0.05) beauty=true;
				else charm=true;
			} else if(TMath::Abs(DsM-1968.)<25.) {
				if(DsIP>0.05) beauty=true;
				else charm=true;
			} else if(TMath::Abs(LcM-2286.)<25.) {
				if(LcIP>0.05) beauty=true;
				else charm=true;
			}
		}
	}

	//now deal with the last event
	if(beauty) foundB.insert(prevEvt);
	if(charm) foundC.insert(prevEvt);

	std::cout << foundB.size() << "\t" << foundC.size() << "\t" << n << std::endl;

	//histograms to plot SVMCor
	TH1D hBT("hBT","",20,0.,10000.);
	TH1D hCT("hCT","",20,0.,10000.);
	TH1D hB("hB","",20,0.,10000.);
	TH1D hC("hC","",20,0.,10000.);
	TH1D hQ("hQ","",20,0.,10000.);

	//count how many jets we find
	int njetsB(0);
	int njetsC(0);
	int njetsQ(0);

	//second loop - fill the histograms
	for(int i=0; i<n; ++i) {
		if(!((int)i % (int)(n/20.))) std::cout << i << " of " << n << std::endl;
		t->GetEntry(i);

		if(JetDijetSVDec==0) continue;

		if(foundB.count(Evt)) {
			if((TMath::Abs(D0M-1864.)<25. && D0IP>0.05) ||
			   (TMath::Abs(D2K3piM-1864.)<25. && D2K3piIP>0.05) || 
			   (TMath::Abs(DM-1870.)<25. && DIP>0.05) ||
			   (TMath::Abs(DsM-1968.)<25. && DsIP>0.05) ||
			   (TMath::Abs(LcM-2286.)<25. && LcIP>0.05)) {
				//tagged beauty jet
				hBT.Fill(SVMCor);
			} else {
				//untagged beauty jet
				++njetsB;
				hB.Fill(SVMCor);
				t5->Fill();
			}
		} else if(foundC.count(Evt)) {
			if((TMath::Abs(D0M-1864.)<25. && D0IP<0.05)  ||
			   (TMath::Abs(D2K3piM-1864.)<25. && D2K3piIP<0.05) || 
			   (TMath::Abs(DM-1870.)<25. && DIP<0.05) ||
			   (TMath::Abs(DsM-1968.)<25. && DsIP<0.05) ||
			   (TMath::Abs(LcM-2286.)<25. && LcIP<0.05)) {
				//tagged charm jet
				hCT.Fill(SVMCor);
			} else {
				//untagged charm jet
				++njetsC;
				hC.Fill(SVMCor);
				t4->Fill();
			}
		} else {
			//event not tagged
			++njetsQ;
			hQ.Fill(SVMCor);
			t0->Fill();
		}
	}


	std::cout << njetsB << "\t" << njetsC << "\t" << njetsQ << std::endl;

	TCanvas c;

	hCT.GetXaxis()->SetTitle("SV M_{cor}");
	
	hB.Scale(1./static_cast<double>(hB.GetEntries() - hB.GetBinContent(0) - hB.GetBinContent(hB.GetNbinsX())));
	hC.Scale(1./static_cast<double>(hC.GetEntries() - hC.GetBinContent(0) - hC.GetBinContent(hC.GetNbinsX())));
	hQ.Scale(1./static_cast<double>(hQ.GetEntries() - hQ.GetBinContent(0) - hQ.GetBinContent(hQ.GetNbinsX())));

	hBT.Scale(1./static_cast<double>(hBT.GetEntries() - hBT.GetBinContent(0) - hBT.GetBinContent(hBT.GetNbinsX())));
	hCT.Scale(1./static_cast<double>(hCT.GetEntries() - hCT.GetBinContent(0) - hCT.GetBinContent(hCT.GetNbinsX())));
	
	hC.SetLineColor(kRed);
	hQ.SetLineColor(kGreen+2);
	
	hCT.SetLineColor(kRed);
	
	hCT.SetLineStyle(kDashed);
	hBT.SetLineStyle(kDashed);
	
	hCT.Draw();
	hBT.Draw("same");

	hC.Draw("same");
	hB.Draw("same");
	hQ.Draw("same");
	
	c.SaveAs("MCorTagged.pdf");

	t0->AutoSave();
	f0->Close();

	t4->AutoSave();
	f4->Close();

	t5->AutoSave();
	f5->Close();
}
