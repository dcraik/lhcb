void makePlot(TString name, TTree* t0, TTree* t4, TTree* t5) {
	TCanvas c;

	double maxX = t0->GetMaximum(name);
	double minX = t0->GetMinimum(name);
	if(t4->GetMaximum(name)>maxX) maxX = t4->GetMaximum(name);
	if(t4->GetMinimum(name)<minX) minX = t4->GetMinimum(name);
	if(t5->GetMaximum(name)>maxX) maxX = t5->GetMaximum(name);
	if(t5->GetMinimum(name)<minX) minX = t5->GetMinimum(name);

	TH1D* h0= new TH1D("h0","",50,minX,maxX);
	TH1D* h4= new TH1D("h4","",50,minX,maxX);
	TH1D* h5= new TH1D("h5","",50,minX,maxX);

	t0->Draw(name+">>h0",name+"!=-1000");
	t4->Draw(name+">>h4",name+"!=-1000");
	t5->Draw(name+">>h5",name+"!=-1000");

	while(h0->Integral(1,10)/h0->Integral() <0.01 && h4->Integral(1,10)/h4->Integral() <0.01 && h5->Integral(1,10)/h5->Integral() <0.01) {
		//std::cout << minX << "," << maxX << std::endl;
		//std::cout << "1-10\t" << h0->Integral(1,10) <<"/"<< h0->Integral() <<"\t"<< h4->Integral(1,10) <<"/"<< h4->Integral() <<"\t"<< h5->Integral(1,10) <<"/"<< h5->Integral() << std::endl;
		minX = h0->GetBinLowEdge(11);

		delete h0;
		delete h4;
		delete h5;

		h0= new TH1D("h0","",50,minX,maxX);
		h4= new TH1D("h4","",50,minX,maxX);
		h5= new TH1D("h5","",50,minX,maxX);

		t0->Draw(name+">>h0",name+"!=-1000");
		t4->Draw(name+">>h4",name+"!=-1000");
		t5->Draw(name+">>h5",name+"!=-1000");
	}
	while(h0->Integral(41,50)/h0->Integral() <0.01 && h4->Integral(41,50)/h4->Integral() <0.01 && h5->Integral(41,50)/h5->Integral() <0.01) {
		//std::cout << minX << "," << maxX << std::endl;
		//std::cout << "41-50\t" << h0->Integral(41,50) <<"/"<< h0->Integral() <<"\t"<< h4->Integral(41,50) <<"/"<< h4->Integral() <<"\t"<< h5->Integral(41,50) <<"/"<< h5->Integral() << std::endl;
		maxX = h0->GetBinLowEdge(41);

		delete h0;
		delete h4;
		delete h5;

		h0= new TH1D("h0","",50,minX,maxX);
		h4= new TH1D("h4","",50,minX,maxX);
		h5= new TH1D("h5","",50,minX,maxX);

		t0->Draw(name+">>h0",name+"!=-1000");
		t4->Draw(name+">>h4",name+"!=-1000");
		t5->Draw(name+">>h5",name+"!=-1000");
	}

	h0->SetLineColor(kGreen+2);
	h5->SetLineColor(kRed);

	h0->Scale(1./h0->Integral());
	h4->Scale(1./h4->Integral());
	h5->Scale(1./h5->Integral());

	h0->SetMaximum(1.1*h0->GetMaximum());
	if(1.1*h4->GetMaximum() > h0->GetMaximum()) h0->SetMaximum(1.1*h4->GetMaximum());
	if(1.1*h5->GetMaximum() > h0->GetMaximum()) h0->SetMaximum(1.1*h5->GetMaximum());

	h0->GetXaxis()->SetTitle(name);

	h0->Draw("");
	h4->Draw("same");
	h5->Draw("same");
	c.SaveAs("plots-MC/"+name+".png");

	delete h0;
	delete h4;
	delete h5;
}

void makePlots() {
	gStyle->SetOptStat(0);

	TFile* f0 = TFile::Open("/tmp/dcraik/for_yandex_0.root");//"../davinci/light.root");
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f4 = TFile::Open("/tmp/dcraik/for_yandex_4.root");//"../davinci/charm.root");
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));

	TFile* f5 = TFile::Open("/tmp/dcraik/for_yandex_5.root");//"../davinci/beauty.root");
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	if(!t0 || !t4 || !t5) return;

	TObjArray* branches = t0->GetListOfBranches();

	for(int i=0; i<branches->GetSize(); ++i) {
		TBranch* b = dynamic_cast<TBranch*>(branches->At(i));
		if(!b) continue;

		TString name = b->GetName();

		makePlot(name,t0,t4,t5);
	}
}
