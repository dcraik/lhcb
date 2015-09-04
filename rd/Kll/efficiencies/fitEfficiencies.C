Double_t fcn(Double_t * abscissa, Double_t * parameter)
{
        Double_t eff(0.0);

        Double_t c = abscissa[0];

        Double_t N = parameter[0];
        Double_t A = parameter[1];
        Double_t B = parameter[2];

	eff = N * ( 1 + A*c*c + B*c*c*c*c );

        return eff;
}



void fitEfficiencies(Int_t weighting=0, Int_t nbins=50, Bool_t fixB=false) {
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();
	gStyle->SetOptStat(0);

	TString weightStr("");
	TString binStr(""); binStr += nbins;
	TString fixStr(""); if(fixB) fixStr+="_fixB";

	switch(weighting) {
		case 0:
			weightStr="PIDonly";
			break;
		case 1:
			weightStr="PID_noMuNN";
			break;
		case 2:
			weightStr="PID_noMuPID";
			break;
		case 3:
			weightStr="PID_trackMult";
			break;
		case 4:
			weightStr="PID_Bkin";
			break;
		case 5:
			weightStr="PID_trackMult_Bkin";
			break;
		default:
			std::cout << "Unknown weighting scheme." << std::endl;
			return;
	}

	std::cout << "Fitting efficiency histograms with " << nbins << " bins and " << weightStr << " weighting..." << std::endl;

	TFile * paramsFile = new TFile("effParams_"+weightStr+fixStr+"_"+binStr+".root","RECREATE");
	TTree * paramsTree = new TTree("params","params");

	Double_t     N(0.),     A(0.),     B(0.);
	Double_t err_N(0.), err_A(0.), err_B(0.);
	Int_t i(0);

	paramsTree->Branch("qSqBin", &i,     "qSqBin/I");
	paramsTree->Branch("N",      &N,     "N/D");
	paramsTree->Branch("A",      &A,     "A/D");
	paramsTree->Branch("B",      &B,     "B/D");
	paramsTree->Branch("err_N",  &err_N, "err_N/D");
	paramsTree->Branch("err_A",  &err_A, "err_A/D");
	paramsTree->Branch("err_B",  &err_B, "err_B/D");
	
	TFile * fHists = TFile::Open("efficiencies_"+weightStr+"_"+binStr+".root");

	TH1D * hist(0);
	TF1  * func = new TF1("func", fcn, -1.0, 1.0, 3);
	TCanvas c;

	for(i=0; i<19; ++i) {
		TString qStr; qStr +=i;
		TString nameE("efficiency_");  nameE+=i;
		hist = dynamic_cast<TH1D*>(fHists->Get(nameE));
		func->SetParNames("N", "A", "B");
		func->SetParameter(0, hist->GetMaximum());
		func->SetParameter(1,-1.0);
		if(fixB) func->FixParameter(2, 0.0);
		else func->SetParameter(2, 0.0);

		TFitResultPtr r = hist->Fit(func,"S");

		hist->GetFunction("func")->SetLineColor(kRed);
		hist->GetXaxis()->SetTitle("cos #theta_{l}");
		hist->GetYaxis()->SetTitle("efficiency");
		hist->Draw();
		c.SaveAs("plots/fit/pdf/effFit"+weightStr+fixStr+"_"+binStr+"_"+qStr+".pdf");
		c.SaveAs("plots/fit/png/effFit"+weightStr+fixStr+"_"+binStr+"_"+qStr+".png");

		N = func->GetParameter(0);
		A = func->GetParameter(1);
		B = func->GetParameter(2);

		TMatrixDSym cm = r->GetCovarianceMatrix();
		err_N = TMath::Sqrt(cm(0,0));
		err_A = TMath::Sqrt(cm(1,1));
		err_B = TMath::Sqrt(cm(2,2));

		paramsTree->Fill();

	}

	paramsTree->AutoSave();
	paramsFile->Close();
}
