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

void getEffsNew(Int_t weighting=0, Int_t nbins=50) {

	TString weightStr("");
	TString binStr(""); binStr += nbins;

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

	std::cout << "Generating efficiency histograms with " << nbins << " bins and " << weightStr << " weighting..." << std::endl;

	TH2D * passed[19];
	TH1D * total[19];
	TEfficiency * eff[19];

	TH1D * tempPassed[50];
	TEfficiency * tempEff[50];

	Int_t ngenInQ[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		             0, 0, 0, 0, 0, 0, 0, 0, 0};
	Int_t nselInQ[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		             0, 0, 0, 0, 0, 0, 0, 0, 0};

	Double_t qSq(0.), costhetal(0.), q1k(0.);

	Double_t ProbNNKweight(0.), ProbNNmuweight(0.), MuonPIDweight(0.), Bplus_P_weight(0.), Bplus_PT_weight(0.), Bplus_ENDVERTEX_CHI2_weight(0.), nTracks_weight(0.);

	for(Int_t i=0; i<19; ++i) {
		TString nameP("passed_");      nameP+=i;
		TString nameT("total_");       nameT+=i;
//		TString nameE("efficiency_");  nameE+=i;

		passed[i] = new TH2D(nameP, "", nbins, -1.0, 1.0, 50, 0., 1.);
		total[i]  = new TH1D(nameT, "", nbins, -1.0, 1.0);
//		eff[i]    = new TEfficiency(nameE, "", nbins, -1.0, 1.0);
	}

	TFile * filegen = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/kmumu_gen_addVars.root");
	TTree * treegen = dynamic_cast<TTree*>(filegen->Get("DecayTuple"));

	TFile * filesel = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromPatrick/2012_mumu.root");
	TTree * treesel = dynamic_cast<TTree*>(filesel->Get("finalTree_KMuMu"));

	treegen->SetBranchAddress("qSq",       &qSq);
	treegen->SetBranchAddress("costhetal", &costhetal);

	treesel->SetBranchAddress("Jpsi_M",    &q1k);
	treesel->SetBranchAddress("costhetal", &costhetal);
	treesel->SetBranchAddress("costhetal", &costhetal);

	treesel->SetBranchAddress("ProbNNKweight",               &ProbNNKweight);
	treesel->SetBranchAddress("ProbNNmuweight",              &ProbNNmuweight);
	treesel->SetBranchAddress("MuonPIDweight",               &MuonPIDweight);
	treesel->SetBranchAddress("Bplus_P_weight",              &Bplus_P_weight);
	treesel->SetBranchAddress("Bplus_PT_weight",             &Bplus_PT_weight);
	treesel->SetBranchAddress("Bplus_ENDVERTEX_CHI2_weight", &Bplus_ENDVERTEX_CHI2_weight);
	treesel->SetBranchAddress("nTracks_weight",              &nTracks_weight);

	Int_t ngen = treegen->GetEntries();
	Int_t nsel = treesel->GetEntries();

	Int_t infoGen = ngen/10;
	Int_t infoSel = nsel/10;

	Double_t weightMin( 9999.);
	Double_t weightMax(-9999.);

	for(Int_t i=0; i<nsel; ++i) {
		if(i%infoSel == 0) std::cout << "Processing entry " << i << " of " << nsel << " selected..." << std::endl;

		Int_t binA(-1), binB(-1);
		Double_t weight(1.0);

		treesel->GetEntry(i);

		switch(weighting) {
			case 0:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight;
				break;
			case 1:
				weight = ProbNNKweight*MuonPIDweight;
				break;
			case 2:
				weight = ProbNNKweight*ProbNNmuweight;
				break;
			case 3:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*nTracks_weight;
				break;
			case 4:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*Bplus_P_weight*Bplus_PT_weight*Bplus_ENDVERTEX_CHI2_weight;
				break;
			case 5:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*Bplus_P_weight*Bplus_PT_weight*Bplus_ENDVERTEX_CHI2_weight*nTracks_weight;
				break;
		}
		if(weight>weightMax) weightMax=weight;
		if(weight<weightMin) weightMin=weight;
	}

	Double_t weightRange = weightMax - weightMin;

	for(Int_t i=0; i<ngen; ++i) {
		if(i%infoGen == 0) std::cout << "Processing entry " << i << " of " << ngen << " generated..." << std::endl;
		Int_t binA(-1), binB(-1);

		treegen->GetEntry(i);

		if(       qSq >  0.10 && qSq <=  0.98) {	binA=0;
		} else if(qSq >  1.10 && qSq <=  2.00) {	binA=1;
		} else if(               qSq <=  3.00) {	binA=2;
		} else if(               qSq <=  4.00) {	binA=3;
		} else if(               qSq <=  5.00) {	binA=4;
		} else if(               qSq <=  6.00) {	binA=5;
		} else if(               qSq <=  7.00) {	binA=6;
		} else if(               qSq <=  8.00) {	binA=7;
		} else if(qSq > 11.00 && qSq <= 11.75) {	binA=8;
		} else if(               qSq <= 12.50) {	binA=9;
		} else if(qSq > 15.00 && qSq <= 16.00) {	binA=10;
		} else if(               qSq <= 17.00) {	binA=11;
		} else if(               qSq <= 18.00) {	binA=12;
		} else if(               qSq <= 19.00) {	binA=13;
		} else if(               qSq <= 20.00) {	binA=14;
		} else if(               qSq <= 21.00) {	binA=15;
		} else if(               qSq <= 22.00) {	binA=16;
		}

		if(       qSq >  1.10 && qSq <=  6.00) {	binB=17;
		} else if(qSq > 15.00 && qSq <= 22.00) {	binB=18;
		}

		if(binA>=0) {
			total[binA]->Fill(costhetal);
			++ngenInQ[binA];
//			eff[binA]->Fill(kTRUE,costhetal);
//		} else {
//			eff[binA]->Fill(kFALSE,costhetal);
		}
		if(binB>=0) {
			total[binB]->Fill(costhetal);
			++ngenInQ[binB];
//			eff[binB]->Fill(kTRUE,costhetal);
//		} else {
//			eff[binB]->Fill(kFALSE,costhetal);
		}
	}

	for(Int_t i=0; i<nsel; ++i) {
		if(i%infoSel == 0) std::cout << "Processing entry " << i << " of " << nsel << " selected..." << std::endl;

		Int_t binA(-1), binB(-1);
		Double_t weight(1.0);

		treesel->GetEntry(i);

		switch(weighting) {
			case 0:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight;
				break;
			case 1:
				weight = ProbNNKweight*MuonPIDweight;
				break;
			case 2:
				weight = ProbNNKweight*ProbNNmuweight;
				break;
			case 3:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*nTracks_weight;
				break;
			case 4:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*Bplus_P_weight*Bplus_PT_weight*Bplus_ENDVERTEX_CHI2_weight;
				break;
			case 5:
				weight = ProbNNKweight*ProbNNmuweight*MuonPIDweight*Bplus_P_weight*Bplus_PT_weight*Bplus_ENDVERTEX_CHI2_weight*nTracks_weight;
				break;
		}

		qSq = q1k*q1k/1.0e6;

		if(       qSq >  0.10 && qSq <=  0.98) {	binA=0;
		} else if(qSq >  1.10 && qSq <=  2.00) {	binA=1;
		} else if(               qSq <=  3.00) {	binA=2;
		} else if(               qSq <=  4.00) {	binA=3;
		} else if(               qSq <=  5.00) {	binA=4;
		} else if(               qSq <=  6.00) {	binA=5;
		} else if(               qSq <=  7.00) {	binA=6;
		} else if(               qSq <=  8.00) {	binA=7;
		} else if(qSq > 11.00 && qSq <= 11.75) {	binA=8;
		} else if(               qSq <= 12.50) {	binA=9;
		} else if(qSq > 15.00 && qSq <= 16.00) {	binA=10;
		} else if(               qSq <= 17.00) {	binA=11;
		} else if(               qSq <= 18.00) {	binA=12;
		} else if(               qSq <= 19.00) {	binA=13;
		} else if(               qSq <= 20.00) {	binA=14;
		} else if(               qSq <= 21.00) {	binA=15;
		} else if(               qSq <= 22.00) {	binA=16;
		}

		if(       qSq >  1.10 && qSq <=  6.00) {	binB=17;
		} else if(qSq > 15.00 && qSq <= 22.00) {	binB=18;
		}

		if(binA>=0) {
			passed[binA]->Fill(costhetal,(weight-weightMin)/weightRange);
			++nselInQ[binA];
		}
		if(binB>=0) {
			passed[binB]->Fill(costhetal,(weight-weightMin)/weightRange);
			++nselInQ[binB];
		}
	}

	printf("\nAverage efficiencies...\n");
	printf("Total:\t%7d\t%7d\n",nsel,ngen);
	TFile * out = new TFile("efficiencies_"+weightStr+"_"+binStr+".root","RECREATE");
	TF1  * func = new TF1("func", fcn, -1.0, 1.0, 3);
	for(Int_t i=0; i<19; ++i) {

		for(Int_t j=0; j<50; ++j) {
			tempPassed[j] = passed[i]->ProjectionX(name,j+1,j+1);
			tempEff[j] = new TEfficiency(*tempPassed[i], *total[i]);
			tempEff[j].SetWeight(weightMin + (j+0.5)*weightRange);
		}

		TString nameE("efficiency_");  nameE+=i;
		eff[i]    = new TEfficiency(*passed[i], *total[i]);
		eff[i]->SetName(nameE);
		func->SetParNames("N", "A", "B");
		func->SetParameter(0, eff[i]->GetEfficiency(nbins/2));
		func->SetParameter(1,-1.0);
		func->FixParameter(2, 0.0);
		TFitResultPtr r = eff[i]->Fit(func,"S");
		eff[i]->Write();
		eff2[i]->Add(passed[i]);
		eff2[i]->Divide(total[i]);
		eff2[i]->Write();
		passed[i]->Write();
		total[i]->Write();

		printf("Bin %2d:\t%7d\t%7d\n",i,nselInQ[i],ngenInQ[i]);
	}

	out->Close();
	filesel->Close();
	filegen->Close();
}
