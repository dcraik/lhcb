void makeTemplates(TString dir="./") {
	TFile* f0 = TFile::Open(dir+"for_yandex_0.root");
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f4 = TFile::Open(dir+"for_yandex_4.root");
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));

	TFile* f5 = TFile::Open(dir+"for_yandex_5.root");
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	if(!t0 || !t4 || !t5) return;

	TFile* fout = TFile::Open(dir+"templates.root","RECREATE");
	TH1D SVM_0("SVM_0","",30,0.,10000.);
	TH1D SVM_4("SVM_4","",30,0.,10000.);
	TH1D SVM_5("SVM_5","",30,0.,10000.);
	TH1D SVMCor_0("SVMCor_0","",30,0.,10000.);
	TH1D SVMCor_4("SVMCor_4","",30,0.,10000.);
	TH1D SVMCor_5("SVMCor_5","",30,0.,10000.);
	TH1D SVSumIPChi2_0("SVSumIPChi2_0","",20,0.,20.);
	TH1D SVSumIPChi2_4("SVSumIPChi2_4","",20,0.,20.);
	TH1D SVSumIPChi2_5("SVSumIPChi2_5","",20,0.,20.);
	TH1D MuIPChi2_0("MuIPChi2_0","",20,0.,20.);
	TH1D MuIPChi2_4("MuIPChi2_4","",20,0.,20.);
	TH1D MuIPChi2_5("MuIPChi2_5","",20,0.,20.);
	TH1D HardIPChi2_0("HardIPChi2_0","",20,0.,20.);
	TH1D HardIPChi2_4("HardIPChi2_4","",20,0.,20.);
	TH1D HardIPChi2_5("HardIPChi2_5","",20,0.,20.);
	TH1D NDispl6_0("NDispl6_0","",12,0.,12.);
	TH1D NDispl6_4("NDispl6_4","",12,0.,12.);
	TH1D NDispl6_5("NDispl6_5","",12,0.,12.);
	TH2D twod_0("twod_0","",20,0.,10000.,10.,2.,12.);
	TH2D twod_4("twod_4","",20,0.,10000.,10.,2.,12.);
	TH2D twod_5("twod_5","",20,0.,10000.,10.,2.,12.);

	t0->Draw("SVM>>SVM_0","JetPT>20000 && JetPT<30000");
	t4->Draw("SVM>>SVM_4","JetPT>20000 && JetPT<30000");
	t5->Draw("SVM>>SVM_5","JetPT>20000 && JetPT<30000");

	t0->Draw("SVMCor>>SVMCor_0","JetPT>20000 && JetPT<30000");
	t4->Draw("SVMCor>>SVMCor_4","JetPT>20000 && JetPT<30000");
	t5->Draw("SVMCor>>SVMCor_5","JetPT>20000 && JetPT<30000");

	t0->Draw("SVN:SVMCor>>twod_0","JetPT>20000 && JetPT<30000");
	t4->Draw("SVN:SVMCor>>twod_4","JetPT>20000 && JetPT<30000");
	t5->Draw("SVN:SVMCor>>twod_5","JetPT>20000 && JetPT<30000");

	t0->Draw("TMath::Log(SVSumIPChi2)>>SVSumIPChi2_0","JetPT>20000 && JetPT<30000");
	t4->Draw("TMath::Log(SVSumIPChi2)>>SVSumIPChi2_4","JetPT>20000 && JetPT<30000");
	t5->Draw("TMath::Log(SVSumIPChi2)>>SVSumIPChi2_5","JetPT>20000 && JetPT<30000");

	t0->Draw("TMath::Log(MuIPChi2)>>MuIPChi2_0","JetPT>20000 && JetPT<30000");
	t4->Draw("TMath::Log(MuIPChi2)>>MuIPChi2_4","JetPT>20000 && JetPT<30000");
	t5->Draw("TMath::Log(MuIPChi2)>>MuIPChi2_5","JetPT>20000 && JetPT<30000");

	t0->Draw("TMath::Log(HardIPChi2)>>HardIPChi2_0","JetPT>20000 && JetPT<30000");
	t4->Draw("TMath::Log(HardIPChi2)>>HardIPChi2_4","JetPT>20000 && JetPT<30000");
	t5->Draw("TMath::Log(HardIPChi2)>>HardIPChi2_5","JetPT>20000 && JetPT<30000");

	t0->Draw("NDispl6>>NDispl6_0","JetPT>20000 && JetPT<30000");
	t4->Draw("NDispl6>>NDispl6_4","JetPT>20000 && JetPT<30000");
	t5->Draw("NDispl6>>NDispl6_5","JetPT>20000 && JetPT<30000");

	SVM_0.Write();
	SVM_4.Write();
	SVM_5.Write();
	SVMCor_0.Write();
	SVMCor_4.Write();
	SVMCor_5.Write();
	SVSumIPChi2_0.Write();
	SVSumIPChi2_4.Write();
	SVSumIPChi2_5.Write();
	MuIPChi2_0.Write();
	MuIPChi2_4.Write();
	MuIPChi2_5.Write();
	HardIPChi2_0.Write();
	HardIPChi2_4.Write();
	HardIPChi2_5.Write();
	NDispl6_0.Write();
	NDispl6_4.Write();
	NDispl6_5.Write();
	twod_0.Write();
	twod_4.Write();
	twod_5.Write();
	fout->Close();

}
