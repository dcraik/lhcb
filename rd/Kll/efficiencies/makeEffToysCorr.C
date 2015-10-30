void doCorrelatedBinFluctuation(TH1* hist, Double_t corrScale)
{
	TRandom* random = LauRandom::randomFun();

	Int_t nPar = static_cast<Int_t>(hist->GetNbinsX());

	TMatrix cm(nPar,nPar);

	Float_t* cmVals = new Float_t[nPar*nPar];

	for(Int_t i=0; i<nPar; ++i) {
		for(Int_t k=0; k<nPar; ++k) {
			Double_t xTerm = (i-k)*(k-i)/(corrScale*corrScale*nBinsX*nBinsX);

			Double_t cov = TMath::Exp(xTerm) * hist->GetBinError(i+1) * hist->GetBinError(k+1);
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
		for(UInt_t k(0); k < iPar; ++k) {
			Double_t tmp = L(k,iPar);
			L(iPar,iPar) -= tmp*tmp;
		}
		if(L(iPar,iPar)<=0.) {
			std::cerr << "ERROR : Numerical imprecision has produced a bad covariance matrix." << std::endl;
			std::cerr << "      : This can be caused by combining a large correlation scale with a finely binned histogram." << std::endl;
			std::cerr << "      : Correlation scale is " << corrScale << " and histogram has " << nBinsX << "x" << nBinsY << " bins." << std::endl;
			std::cerr << "      : Try reducing the correlation scale." << std::endl;
			gSystem->Exit(EXIT_FAILURE);
		}
		L(iPar,iPar) = TMath::Sqrt( L(iPar,iPar) );
		// then the off-diagonal terms
		for(UInt_t jPar(iPar+1); jPar < nPar; ++jPar) {
			L(iPar,jPar) = cm(iPar,jPar);
			for(UInt_t k(0); k < iPar; ++k) {
				L(iPar,jPar) -= L(k,iPar)*L(k,jPar);
			}
			L(iPar,jPar) /= L(iPar,iPar);
		}
	}
	// transpose to get Lt
	TMatrix * Lt = new TMatrix(TMatrix::kTransposed,L);

	LauFitData genData;

	// create a vector of random unit Gaussian variables 
	TVector g(nPar);
	for (UInt_t i(0); i<nPar; ++i) {
		g[i] = random->Gaus(); 
	}

	// multiply this vector by Lt to introduce the appropriate correlations
	g *= (*Lt);

	for (Int_t i(0); i<nBinsX; ++i) {
		for (Int_t j(0); j<nBinsY; ++j) {
			hist->SetBinContent(i+1,j+1,hist->GetBinContent(i+1,j+1)+g[j*nBinsX + i]);
		}
	}
}

void makeEffToys() {
	TFile* fin = TFile::Open("");

	TFile* fout = new TFile("test.root","RECREATE");

	for(Int_t i=0; i<19; ++i) {
		TString hName = "efficiency_"; hName+=i;
		TH1* hin = dynamic_cast<TH1*>(fin->Get(hName));

		TH1* hout = tin->Clone();
		doCorrelatedBinFluctuations(hout,0.1);
	}

	hout->Write();
	fout->Close();
}
