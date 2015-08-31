void getDVetoCorrections(Int_t nbins=200) {

	TString binStr(""); binStr += nbins;

	std::cout << "Generating D veto correction histograms with " << nbins << " bins..." << std::endl;

	TH1D * passed[19];
	TH1D * total[19];
	TH1D * eff[19];

	Int_t ngenInQ[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		             0, 0, 0, 0, 0, 0, 0, 0, 0};
	Int_t nselInQ[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		             0, 0, 0, 0, 0, 0, 0, 0, 0};

	Double_t qSq(0.), costhetal(0.), mKmu_D;

	for(Int_t i=0; i<19; ++i) {
		TString nameP("passed_");      nameP+=i;
		TString nameT("total_");       nameT+=i;
		TString nameE("efficiency_");  nameE+=i;

		passed[i] = new TH1D(nameP, "", nbins, -1.0, 1.0);
		total[i]  = new TH1D(nameT, "", nbins, -1.0, 1.0);
		eff[i]    = new TH1D(nameE, "", nbins, -1.0, 1.0);

		passed[i]->Sumw2();
		total[i]->Sumw2();
		eff[i]->Sumw2();
	}

	TFile * filegen = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/kmumu_gen_addVars.root");
	TTree * treegen = dynamic_cast<TTree*>(filegen->Get("DecayTuple"));

	treegen->SetBranchAddress("qSq",       &qSq);
	treegen->SetBranchAddress("costhetal", &costhetal);
	treegen->SetBranchAddress("mKmu_D",    &mKmu_D);

	Int_t ngen = treegen->GetEntries();
	Int_t nsel(0);

	Int_t infoGen = ngen/10;

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
			if(mKmu_D < 1850. || mKmu_D > 1880.) {
				passed[binA]->Fill(costhetal);
				++nselInQ[binA];
			}
		}
		if(binB>=0) {
			total[binB]->Fill(costhetal);
			++ngenInQ[binB];
			if(mKmu_D < 1850. || mKmu_D > 1880.) {
				passed[binB]->Fill(costhetal);
				++nselInQ[binB];
			}
		}
		if(mKmu_D < 1850. || mKmu_D > 1880.) {
			++nsel;
		}
	}

	printf("\nAverage efficiencies...\n");
	printf("Total:\t%7d\t%7d\n",nsel,ngen);
	TFile * out = new TFile("Dveto_"+binStr+".root","RECREATE");
	for(Int_t i=0; i<19; ++i) {
		TString qStr; qStr+=i;

		eff[i]->Add(passed[i]);
		eff[i]->Divide(total[i]);
		eff[i]->Write();
		passed[i]->Write();
		total[i]->Write();

		TCanvas c;
		eff[i]->Draw();
		c.SaveAs("plots/Dveto/Q"+qStr+"_"+binStr+".pdf");
		c.SaveAs("plots/Dveto/Q"+qStr+"_"+binStr+".png");

		printf("Bin %2d:\t%7d\t%7d\n",i,nselInQ[i],ngenInQ[i]);
	}

	out->Close();
	filegen->Close();
}
