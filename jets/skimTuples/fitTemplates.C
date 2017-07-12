double fcn(double* x, double* pars, TString varName) {
	static std::map<TString,TH1D*> hist_0;
	static std::map<TString,TH1D*> hist_4;
	static std::map<TString,TH1D*> hist_5;

	if (hist_0.count(varName)==0) {
		std::cout << "foo" << std::endl;
		TFile* ft = TFile::Open("templates.root");
		hist_0[varName] = dynamic_cast<TH1D*>(ft->Get(varName+"_0"));
		hist_4[varName] = dynamic_cast<TH1D*>(ft->Get(varName+"_4"));
		hist_5[varName] = dynamic_cast<TH1D*>(ft->Get(varName+"_5"));
		hist_0[varName]->Scale(1./hist_0[varName]->GetSumOfWeights());
		hist_4[varName]->Scale(1./hist_4[varName]->GetSumOfWeights());
		hist_5[varName]->Scale(1./hist_5[varName]->GetSumOfWeights());
	}

	const double xx = x[0];
	const double N  = pars[0];
	const double w1 = pars[1];
	const double w2 = pars[2];
	const double w3 = 1 - pars[1] - pars[2];
	
	const double y1 = hist_4[varName]->GetBinContent(hist_0[varName]->GetXaxis()->FindFixBin(xx));
	const double y2 = hist_5[varName]->GetBinContent(hist_4[varName]->GetXaxis()->FindFixBin(xx));
	const double y3 = hist_0[varName]->GetBinContent(hist_5[varName]->GetXaxis()->FindFixBin(xx));

	return N*(w1*y1 + w2*y2 + w3*y3);
}

double svm(double* x, double* pars) {
	return fcn(x, pars, "SVM");
}

double svmcor(double* x, double* pars) {
	return fcn(x, pars, "SVMCor");
}

double svsumipchi2(double* x, double* pars) {
	return fcn(x, pars, "SVSumIPChi2");
}

double muipchi2(double* x, double* pars) {
	return fcn(x, pars, "MuIPChi2");
}

double hardipchi2(double* x, double* pars) {
	return fcn(x, pars, "HardIPChi2");
}

double ndispl6(double* x, double* pars) {
	return fcn(x, pars, "NDispl6");
}

void fitTemplate(TString hname, TString varName) {
	TFile* ft = TFile::Open("templates.root");
	TH1D* sim_0 = dynamic_cast<TH1D*>(ft->Get(hname+"_0"));
	TH1D* sim_4 = dynamic_cast<TH1D*>(ft->Get(hname+"_4"));
	TH1D* sim_5 = dynamic_cast<TH1D*>(ft->Get(hname+"_5"));

	int n= sim_0->GetXaxis()->GetNbins();
	double min = sim_0->GetXaxis()->GetBinLowEdge(1);
	double max = sim_0->GetXaxis()->GetBinUpEdge(n);

	TFile* f0 = TFile::Open("for_yandex_data_SV_0tag.root");
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));
	TFile* f4 = TFile::Open("for_yandex_data_SV_4tag.root");
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));
	TFile* f5 = TFile::Open("for_yandex_data_SV_5tag.root");
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	TH1D data0(hname+"data0","",n,min,max);
	TH1D data4(hname+"data4","",n,min,max);
	TH1D data5(hname+"data5","",n,min,max);

	t0->Draw(varName+">>"+hname+"data0","JetPT>20000. && JetPT<30000.");
	t4->Draw(varName+">>"+hname+"data4","JetPT>20000. && JetPT<30000.");
	t5->Draw(varName+">>"+hname+"data5","JetPT>20000. && JetPT<30000.");

	double N[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double S[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

	TF1 f;
	if(hname=="SVM") f = TF1("fcn", svm, min, max, 3);
	else if(hname=="SVMCor") f = TF1("fcn", svmcor, min, max, 3);
	else if(hname=="SVSumIPChi2") f = TF1("fcn", svsumipchi2, min, max, 3);
	else if(hname=="MuIPChi2") f = TF1("fcn", muipchi2, min, max, 3);
	else if(hname=="HardIPChi2") f = TF1("fcn", hardipchi2, min, max, 3);
	else if(hname=="NDispl6") f = TF1("fcn", ndispl6, min, max, 3);
	else {
		std::cout << "unknown variable " << hname << std::endl;
		return;
	}
	f.SetParLimits(1,0.,1.);
	f.SetParLimits(2,0.,1.);

	f.SetParameter(0,data0.GetEntries());
	f.SetParameter(1,0.05);
	f.SetParameter(2,0.36);
	data0.Fit(&f, "L");
	N[0] = f.GetParameter(0) * f.GetParameter(1);
	N[1] = f.GetParameter(0) * f.GetParameter(2);
	N[2] = f.GetParameter(0) - N[1] - N[2];

	std::cout << f.GetParameter(1) << "+/-" << f.GetParError(1) 
	  << "\t" << f.GetParameter(2) << "+/-" << f.GetParError(2) 
	  << "\t" << 1. - f.GetParameter(1) - f.GetParameter(2) << "+/-" << TMath::Sqrt(TMath::Power(f.GetParError(1),2) + TMath::Power(f.GetParError(1),2)) << std::endl;

	f.SetParameter(0,data4.GetEntries());
	f.SetParameter(1,0.10);
	f.SetParameter(2,0.62);
	data4.Fit(&f, "L");
	N[3] = f.GetParameter(0) * f.GetParameter(1);
	N[4] = f.GetParameter(0) * f.GetParameter(2);
	N[5] = f.GetParameter(0) - N[4] - N[5];

	std::cout << f.GetParameter(1) << "+/-" << f.GetParError(1) 
	  << "\t" << f.GetParameter(2) << "+/-" << f.GetParError(2) 
	  << "\t" << 1. - f.GetParameter(1) - f.GetParameter(2) << "+/-" << TMath::Sqrt(TMath::Power(f.GetParError(1),2) + TMath::Power(f.GetParError(1),2)) << std::endl;

	f.SetParameter(0,data5.GetEntries());
	f.SetParameter(1,0.03);
	f.SetParameter(2,0.21);
	data5.Fit(&f, "L");
	N[6] = f.GetParameter(0) * f.GetParameter(1);
	N[7] = f.GetParameter(0) * f.GetParameter(2);
	N[8] = f.GetParameter(0) - N[7] - N[8];

	std::cout << f.GetParameter(1) << "+/-" << f.GetParError(1) 
	  << "\t" << f.GetParameter(2) << "+/-" << f.GetParError(2) 
	  << "\t" << 1. - f.GetParameter(1) - f.GetParameter(2) << "+/-" << TMath::Sqrt(TMath::Power(f.GetParError(1),2) + TMath::Power(f.GetParError(1),2)) << std::endl;

	TCanvas c;
	sim_0->SetLineColor(kRed);
	sim_4->SetLineColor(kGreen+2);
	sim_5->SetLineColor(kMagenta);

	sim_4->Scale(N[0]/sim_4->GetSumOfWeights());
	sim_5->Scale(N[1]/sim_5->GetSumOfWeights());
	sim_0->Scale(N[2]/sim_0->GetSumOfWeights());
	data0.Draw();
	sim_0->Draw("same");
	sim_4->Draw("same");
	sim_5->Draw("same");
	c.SaveAs(hname+"fit0.pdf");

	sim_4->Scale(N[3]/N[0]);
	sim_5->Scale(N[4]/N[1]);
	sim_0->Scale(N[5]/N[2]);
	data4.Draw();
	sim_0->Draw("same");
	sim_4->Draw("same");
	sim_5->Draw("same");
	c.SaveAs(hname+"fit4.pdf");

	sim_4->Scale(N[6]/N[3]);
	sim_5->Scale(N[7]/N[4]);
	sim_0->Scale(N[8]/N[5]);
	data5.Draw();
	sim_0->Draw("same");
	sim_4->Draw("same");
	sim_5->Draw("same");
	c.SaveAs(hname+"fit5.pdf");
}

void fitTemplates() {
	fitTemplate("SVM","SVM");
	fitTemplate("SVMCor","SVMCor");
	fitTemplate("SVSumIPChi2","TMath::Log(SVSumIPChi2)");
	fitTemplate("MuIPChi2","TMath::Log(MuIPChi2)");
	fitTemplate("HardIPChi2","TMath::Log(HardIPChi2)");
	fitTemplate("NDispl6","NDispl6");
}
