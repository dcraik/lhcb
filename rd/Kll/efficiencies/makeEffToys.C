TRandom* r;

Bool_t doCorrelatedBinFluctuation(TEfficiency* hin, TH1* hout, std::vector<Int_t> bins, Double_t corrScale)
{
	Int_t nBinsX = hout->GetNbinsX();
	Int_t nPar = bins.size();

	TMatrix cm(nPar,nPar);

	Float_t* cmVals = new Float_t[nPar*nPar];

	for(Int_t i=0; i<nPar; ++i) {
			
		Int_t ii = bins[i];
		Double_t errI = (hin->GetEfficiencyErrorUp(ii) + hin->GetEfficiencyErrorLow(ii)) / 2.;

		for(Int_t k=0; k<nPar; ++k) {
			Int_t kk = bins[k];
			Double_t errK = (hin->GetEfficiencyErrorUp(kk) + hin->GetEfficiencyErrorLow(kk)) / 2.;

			Double_t xTerm = (ii-kk)*(kk-ii)/(corrScale*corrScale*nBinsX*nBinsX);

			Double_t cov = TMath::Exp(xTerm) * errI * errK;
			cmVals[k*nPar + i] = cov;
		}
	}

	cm.SetMatrixArray(cmVals);
	delete[] cmVals;

	// calculate the elements of the upper-triangular matrix L that gives
	// Lt*L = C where Lt is the transpose of L (the "square-root method")
	TMatrix L(nPar,nPar);
	for(UInt_t iPar(0); iPar < nPar; ++iPar) {
		// calculate the diagonal term first
		L(iPar,iPar) = cm(iPar,iPar);
		for(UInt_t kPar(0); kPar < iPar; ++kPar) {
			Double_t tmp = L(kPar,iPar);
			L(iPar,iPar) -= tmp*tmp;
		}
		if(L(iPar,iPar)<=0.) {
			std::cerr << "ERROR : Numerical imprecision has produced a bad covariance matrix." << std::endl;
			std::cerr << "      : This can be caused by combining a large correlation scale with a finely binned histogram." << std::endl;
			std::cerr << "      : Correlation scale is " << corrScale << " and histogram has " << nBinsX << " bins." << std::endl;
			std::cerr << "      : Reducing the correlation scale by a factor of 2." << std::endl;

			return kFALSE;
		}
		L(iPar,iPar) = TMath::Sqrt( L(iPar,iPar) );
		// then the off-diagonal terms
		for(UInt_t jPar(iPar+1); jPar < nPar; ++jPar) {
			L(iPar,jPar) = cm(iPar,jPar);
			for(UInt_t kPar(0); kPar < iPar; ++kPar) {
				L(iPar,jPar) -= L(kPar,iPar)*L(kPar,jPar);
			}
			L(iPar,jPar) /= L(iPar,iPar);
		}
	}
	// transpose to get Lt
	TMatrix * Lt = new TMatrix(TMatrix::kTransposed,L);

	// create a vector of random unit Gaussian variables 
	TVector g(nPar);
	for (UInt_t j(0); j<nPar; ++j) {
		g[j] = r->Gaus(); 
	}

	// multiply this vector by Lt to introduce the appropriate correlations
	g *= (*Lt);

	for (Int_t i(0); i<nPar; ++i) {
		Double_t eff = hin->GetEfficiency(bins[i])+g[i];
		if(eff > 1.0) eff = 1.0;
		hout->SetBinContent(bins[i], eff);
	}

	return kTRUE;
}

void makeEffToys(Int_t seed, TString veto="D") {
//	r.SetSeed(seed);
	r = RooRandom::randomGenerator();
	r->SetSeed(seed);

	TFile* fin = TFile::Open(veto+"veto_200.root");

	TString fName("toys/"); fName+=seed; fName+="/"+veto+"veto_200.root";

	TFile* fout = new TFile(fName,"RECREATE");

	for(Int_t j=0; j<21; ++j) {
		TString hName  = "efficiency_"; hName+=j;
		TString hName2 = "efficiencyHist_"; hName2+=j;
		TEfficiency* hin = dynamic_cast<TEfficiency*>(fin->Get(hName));

		TH1* hout = dynamic_cast<TH1*>(hin->GetTotalHistogram()->Clone(hName2));

		std::vector<Int_t> corrBins;

		Int_t n = hout->GetNbinsX();
		for(Int_t i=0; i<n; ++i) {
			Double_t eff = hin->GetEfficiency(i+1);
			Double_t erm = hin->GetEfficiencyErrorLow(i+1);
			Double_t erp = hin->GetEfficiencyErrorUp(i+1);

			Bool_t fluctuate = kTRUE;

			// don't fluctuate if the veto hasn't affected this bin
			// also ignore the odd missing entry - not sure what causes these but they don't seem reasonable
			if((eff > 0.99 && eff + erp > 0.999) //efficiency close to 1 and not significantly different
			|| ((i<1 || hin->GetEfficiency(i) == 1) && (i>n-1 || hin->GetEfficiency(i+2) == 1))) { //single bin dip (careful with this one)
				eff = 1;
				fluctuate = kFALSE;
			}
			if(eff < 0.01 && eff - erm < 0.001) {//efficiency close to 0 and not significantly different
				eff = 0;
				fluctuate = kFALSE;
			}

			//otherwise we're fluctuating the bin
			if(fluctuate) {
				//if the errors are roughly symmetric then we can symmetrise them and introduce some correlation between neighbouring bins
				//this is difficult to do with asymmetric errors so if asymmetry > 10% lets just ignore correlations
				//note that a large asymmetry in neighbouring bins will also lead to same-sign fluctuations anyway
				if((erm - erp) / (erm + erp) < 0.1) {
					//correlation is more important than asymmetry

					//add bin to the list to be fluctuated later
					corrBins.push_back(i+1);

				} else {
					//asymmetry is more important than correlation

					//first catch any cases on a limit (the previous checks for eff > 0.99 and eff < 0.01 should catch these but play it safe)
					if(erm <= 0) {
						//vary with a half Gaussian
						eff += TMath::Abs(r->Gaus(0.,erp));
					} else if(erp <= 0) {
						//vary with a half Gaussian
						eff -= TMath::Abs(r->Gaus(0.,erm));
					} else {
						//vary with a bifurcated Gaussian
						RooRealVar effVar( "effVar", "",-1.,2.);
						RooRealVar muVar(  "muVar",  "",eff);
						RooRealVar sigmVar("sigmVar","",erm);
						RooRealVar sigpVar("sigpVar","",erp);

						RooBifurGauss pdf("pdf","",effVar,muVar,sigmVar,sigpVar);
						RooDataSet* ds = pdf.generate(RooArgSet(effVar),1);
						eff = ds->get(0)->getRealValue("effVar");
						delete ds;
					}
				}
			}
			if(eff > 1.0) eff = 1.0; //std::cout << i << "\t" << eff << "\t" << erp << "\t" << erm << std::endl;
//			std::cout << hin->GetEfficiency(i+1) << "\t" << eff << std::endl;

			hout->SetBinContent(i+1, eff);
		}

		//now deal with the correlated efficiencies
		Double_t corrFactor(0.01);

		while(!doCorrelatedBinFluctuation(hin,hout,corrBins,corrFactor)) {
			corrFactor /= 2.;
		}

//		std::cout << std::endl;

		TCanvas c;
		hin->Draw();
		hout->SetMarkerColor(kRed);
		hout->SetMarkerStyle(4);
		hout->Draw("Psame");
		TString pName = "plots/toys/"; pName+=seed; pName+="/"+veto+"veto_Q"; pName+=j; pName+=".pdf";
		c.SaveAs(pName);
	}

	hout->Write();
	fout->Close();
}
