{
	std::ofstream fout("eff.dat");

	fout << "qSq bin\tmoment\tvalue\t\tstat\t\tsyst\t\t\t\tbias1\t\tbias2" << std::endl;

	for(Int_t q=0; q<19; ++q) {
		TString toyName("eff/"); toyName+=q; toyName+="_results.root";
		TString dataName("results_"); dataName+=q; dataName+="_P.root";

		TFile* toyFile  = TFile::Open(toyName);
		TFile* dataFile = TFile::Open(dataName);

		TTree* toy   = dynamic_cast<TTree*>(toyFile->Get("total"));
		TTree* data  = dynamic_cast<TTree*>(dataFile->Get("total"));
		TTree* error = dynamic_cast<TTree*>(dataFile->Get("error"));

		for(Int_t g=1; g<5; ++g) {
			TString gName("G00"); gName+=g;
			
			//Double_t min = TMath::Min(toy->GetMinimum(gName),data->GetMinimum(gName)-error->GetMaximum(gName));
			//Double_t max = TMath::Max(toy->GetMaximum(gName),data->GetMaximum(gName)+error->GetMaximum(gName));
			Double_t min = toy->GetMinimum(gName);
			Double_t max = toy->GetMaximum(gName);
			Double_t range = max-min;
			if(range<=0) continue;

			min-=0.10*range;
			max+=0.10*range;

			TH1D h("h", "", 30, min, max);
			toy->Draw(gName+">>h");

			TF1 gaus("gaus","gaus(0)",min,max);
			gaus.SetParameter(0,100);
			gaus.SetParameter(1,h.GetBinCenter(h.GetMaximumBin()));
			gaus.SetParameter(2,h.GetStdDev());
			gaus.SetParLimits(1,min+0.2*range,max-0.2*range);
			gaus.SetParLimits(2,h.GetStdDev()/3.,h.GetStdDev()*3.);
			h->Fit(&gaus,"S");

			TCanvas c;
			h->Draw();
			TString pName(""); pName+="effFits/Q"; pName+=q; pName+="_"; pName+=gName; pName+=".pdf";
			c.SaveAs(pName);

			fout << q << "\t" << g << "\t" << data->GetMinimum(gName) << "\t" << error->GetMinimum(gName) << "\t" << gaus.GetParameter(2) << "+/-" << gaus.GetParError(2) << "\t" 
			     << data->GetMinimum(gName)-gaus.GetParameter(1) << "\t" << gaus.GetParError(1) << std::endl;

			//std::cout << data->GetMinimum(gName)  << "\t" << gaus.GetParameter(1) << "+/-" << gaus.GetParError(1) << std::endl;
			//std::cout << error->GetMinimum(gName) << "\t" << gaus.GetParameter(2) << "+/-" << gaus.GetParError(2) << std::endl;
		}

		fout << std::endl;
	}

	fout.close();

}
