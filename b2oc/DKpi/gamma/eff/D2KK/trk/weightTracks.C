{
  Double_t mB(5278.34), sB(11.94);

  Int_t nBinsP(5), nBinsEta(2), nBinsMult(5);

  Float_t pBins[nBinsP+1], etaBins[nBinsEta+1], multBins[nBinsMult+1];

  TFile* dataFile    = TFile::Open("/data/lhcb/phrkbf/B2DKpi/Job2012_B2DKpi_D2KK_Scaled/B2D0Kpi_D02KK_selBd_Dsignal_Dst25_vetoes_PID3_NND2KK_addIMs.root");
  TFile* mcFile      = TFile::Open("/data/lhcb/phrkbf/B2DKpi/JobSim08aMC/D2KK/Bd2D0Kpi/D2KK_Bd2D0Kpi_selBd_Dsignal_vetoes_noDPV_NND2KK_addIMs_addMisIDIMs.root");

  TTree* dataTree    = dataFile->Get("DecayTree");
  TTree*   mcTree    =   mcFile->Get("DecayTree");

  Int_t nEntriesData = dataTree->GetEntries();
  Int_t nEntriesMC   =   mcTree->GetEntries();

  Double_t B_D0_B_CM_M(0.);
  Double_t B_D0_p_CM_P(0.),   B_D0_p_CM_PX(0.),   B_D0_p_CM_PY(0.),   B_D0_p_CM_PZ(0.);
  Double_t B_D0_m_CM_P(0.),   B_D0_m_CM_PX(0.),   B_D0_m_CM_PY(0.),   B_D0_m_CM_PZ(0.);
  Double_t B_D0_D0p_CM_P(0.), B_D0_D0p_CM_PX(0.), B_D0_D0p_CM_PY(0.), B_D0_D0p_CM_PZ(0.);
  Double_t B_D0_D0m_CM_P(0.), B_D0_D0m_CM_PX(0.), B_D0_D0m_CM_PY(0.), B_D0_D0m_CM_PZ(0.);
  Int_t pi_ID(0.), D0m_ID(0.), nTracks(0.);
  Int_t B_TRUEID,D_TRUEID,K_TRUEID,pi_TRUEID,D0p_TRUEID,D0m_TRUEID,K_MC_MOTHER_KEY,pi_MC_MOTHER_KEY,D0p_MC_GD_MOTHER_KEY,D0m_MC_GD_MOTHER_KEY;

  dataTree->SetBranchAddress("B_D0_B_CM_M",    &B_D0_B_CM_M   );
  dataTree->SetBranchAddress("B_D0_p_CM_PX",   &B_D0_p_CM_PX  );
  dataTree->SetBranchAddress("B_D0_p_CM_PY",   &B_D0_p_CM_PY  );
  dataTree->SetBranchAddress("B_D0_p_CM_PZ",   &B_D0_p_CM_PZ  );
  dataTree->SetBranchAddress("B_D0_m_CM_PX",   &B_D0_m_CM_PX  );
  dataTree->SetBranchAddress("B_D0_m_CM_PY",   &B_D0_m_CM_PY  );
  dataTree->SetBranchAddress("B_D0_m_CM_PZ",   &B_D0_m_CM_PZ  );
  dataTree->SetBranchAddress("B_D0_D0p_CM_PX", &B_D0_D0p_CM_PX);
  dataTree->SetBranchAddress("B_D0_D0p_CM_PY", &B_D0_D0p_CM_PY);
  dataTree->SetBranchAddress("B_D0_D0p_CM_PZ", &B_D0_D0p_CM_PZ);
  dataTree->SetBranchAddress("B_D0_D0m_CM_PX", &B_D0_D0m_CM_PX);
  dataTree->SetBranchAddress("B_D0_D0m_CM_PY", &B_D0_D0m_CM_PY);
  dataTree->SetBranchAddress("B_D0_D0m_CM_PZ", &B_D0_D0m_CM_PZ);
  dataTree->SetBranchAddress("pi_ID",          &pi_ID         );
  dataTree->SetBranchAddress("D0m_ID",        &D0m_ID       );
  dataTree->SetBranchAddress("nTracks",        &nTracks       );
  mcTree->SetBranchAddress(  "B_D0_p_CM_PX",   &B_D0_p_CM_PX  );
  mcTree->SetBranchAddress(  "B_D0_p_CM_PY",   &B_D0_p_CM_PY  );
  mcTree->SetBranchAddress(  "B_D0_p_CM_PZ",   &B_D0_p_CM_PZ  );
  mcTree->SetBranchAddress(  "B_D0_m_CM_PX",   &B_D0_m_CM_PX  );
  mcTree->SetBranchAddress(  "B_D0_m_CM_PY",   &B_D0_m_CM_PY  );
  mcTree->SetBranchAddress(  "B_D0_m_CM_PZ",   &B_D0_m_CM_PZ  );
  mcTree->SetBranchAddress(  "B_D0_D0p_CM_PX", &B_D0_D0p_CM_PX);
  mcTree->SetBranchAddress(  "B_D0_D0p_CM_PY", &B_D0_D0p_CM_PY);
  mcTree->SetBranchAddress(  "B_D0_D0p_CM_PZ", &B_D0_D0p_CM_PZ);
  mcTree->SetBranchAddress(  "B_D0_D0m_CM_PX", &B_D0_D0m_CM_PX);
  mcTree->SetBranchAddress(  "B_D0_D0m_CM_PY", &B_D0_D0m_CM_PY);
  mcTree->SetBranchAddress(  "B_D0_D0m_CM_PZ", &B_D0_D0m_CM_PZ);
  mcTree->SetBranchAddress(  "pi_ID",          &pi_ID         );
  mcTree->SetBranchAddress(  "D0m_ID",        &D0m_ID       );
  mcTree->SetBranchAddress(  "nTracks",        &nTracks       );

  mcTree->SetBranchAddress("B_TRUEID",                &B_TRUEID);
  mcTree->SetBranchAddress("D_TRUEID",                &D_TRUEID);
  mcTree->SetBranchAddress("D0p_TRUEID",              &D0p_TRUEID);
  mcTree->SetBranchAddress("D0m_TRUEID",             &D0m_TRUEID);
  mcTree->SetBranchAddress("K_TRUEID",                &K_TRUEID);
  mcTree->SetBranchAddress("pi_TRUEID",               &pi_TRUEID);
  mcTree->SetBranchAddress("D0p_MC_GD_MOTHER_KEY",    &D0p_MC_GD_MOTHER_KEY);
  mcTree->SetBranchAddress("D0m_MC_GD_MOTHER_KEY",   &D0m_MC_GD_MOTHER_KEY);
  mcTree->SetBranchAddress("K_MC_MOTHER_KEY",         &K_MC_MOTHER_KEY);
  mcTree->SetBranchAddress("pi_MC_MOTHER_KEY",        &pi_MC_MOTHER_KEY);

  pBins[0] = 0;
  pBins[nBinsP] = 200000.;

  Int_t perBin   = 4*(nEntriesMC/nBinsP);
  Int_t perBinLo = 4*(nEntriesMC/nBinsP)*0.99;
  Int_t perBinHi = 4*(nEntriesMC/nBinsP)*1.01;
  Float_t cutVal(0.);
  Float_t minStep=100.;
  for(Int_t i=0; i<nBinsP-1; ++i) {
    Float_t stepSize = 100.;
    cout << i << endl;
    while(true) {
      TString cut1="(B_D0_p_CM_PX**2+B_D0_p_CM_PY**2+B_D0_p_CM_PZ**2)**.5<";
      cut1 += cutVal;
      TString cut2="(B_D0_m_CM_PX**2+B_D0_m_CM_PY**2+B_D0_m_CM_PZ**2)**.5<";
      cut2 += cutVal;
      TString cut3="(B_D0_D0p_CM_PX**2+B_D0_D0p_CM_PY**2+B_D0_D0p_CM_PZ**2)**.5<";
      cut3 += cutVal;
      TString cut4="(B_D0_D0m_CM_PX**2+B_D0_D0m_CM_PY**2+B_D0_D0m_CM_PZ**2)**.5<";
      cut4 += cutVal;
      
      Int_t nBelow = mcTree->GetEntries(cut1);
      nBelow += mcTree->GetEntries(cut2);
      nBelow += mcTree->GetEntries(cut3);
      nBelow += mcTree->GetEntries(cut4);
      
      if(nBelow < perBin*i + perBinLo) {
        cutVal += stepSize;
	stepSize *= 1.2;
	continue;
      }
      if(nBelow > perBin*i + perBinHi) {
        cutVal -= stepSize / 2.;
	stepSize /= 2.;
	if(stepSize>minStep/5.) continue;
      }

      pBins[i+1] = TMath::Nint(cutVal/minStep)*minStep;
      break;

    }
  }

  etaBins[0] = 1.5;
  etaBins[nBinsEta] = 5.5;

  Int_t perBin   = 4*(nEntriesMC/nBinsEta);
  Int_t perBinLo = 4*(nEntriesMC/nBinsEta)*0.99;
  Int_t perBinHi = 4*(nEntriesMC/nBinsEta)*1.01;
  Float_t cutVal(1.5);
  Float_t minStep=0.01;
  for(Int_t i=0; i<nBinsEta-1; ++i) {
    Float_t stepSize = 0.01;
    cout << i << endl;
    while(true) {
      TString cut1="0.5*log(((B_D0_p_CM_PX**2+B_D0_p_CM_PY**2+B_D0_p_CM_PZ**2)**.5+B_D0_p_CM_PZ)/((B_D0_p_CM_PX**2+B_D0_p_CM_PY**2+B_D0_p_CM_PZ**2)**.5-B_D0_p_CM_PZ))<";
      cut1 += cutVal;
      TString cut2="0.5*log(((B_D0_m_CM_PX**2+B_D0_m_CM_PY**2+B_D0_m_CM_PZ**2)**.5+B_D0_m_CM_PZ)/((B_D0_m_CM_PX**2+B_D0_m_CM_PY**2+B_D0_m_CM_PZ**2)**.5-B_D0_m_CM_PZ))<";
      cut2 += cutVal;
      TString cut3="0.5*log(((B_D0_D0p_CM_PX**2+B_D0_D0p_CM_PY**2+B_D0_D0p_CM_PZ**2)**.5+B_D0_D0p_CM_PZ)/((B_D0_D0p_CM_PX**2+B_D0_D0p_CM_PY**2+B_D0_D0p_CM_PZ**2)**.5-B_D0_D0p_CM_PZ))<";
      cut3 += cutVal;
      TString cut4="0.5*log(((B_D0_D0m_CM_PX**2+B_D0_D0m_CM_PY**2+B_D0_D0m_CM_PZ**2)**.5+B_D0_D0m_CM_PZ)/((B_D0_D0m_CM_PX**2+B_D0_D0m_CM_PY**2+B_D0_D0m_CM_PZ**2)**.5-B_D0_D0m_CM_PZ))<";
      cut4 += cutVal;
      
      Int_t nBelow = mcTree->GetEntries(cut1);
      nBelow += mcTree->GetEntries(cut2);
      nBelow += mcTree->GetEntries(cut3);
      nBelow += mcTree->GetEntries(cut4);
      
      if(nBelow < perBin*i + perBinLo) {
        cutVal += stepSize;
	stepSize *= 1.2;
	continue;
      }
      if(nBelow > perBin*i + perBinHi) {
        cutVal -= stepSize / 2.;
	stepSize /= 2.;
	if(stepSize>minStep/5.) continue;
      }
      
      etaBins[i+1] = TMath::Nint(cutVal/minStep)*minStep;
      break;
      
    }
  }

  multBins[0] = 0;
  multBins[nBinsMult] = 1000.;

  Int_t perBin   = (nEntriesMC/nBinsMult);
  Int_t perBinLo = (nEntriesMC/nBinsMult)*0.99;
  Int_t perBinHi = (nEntriesMC/nBinsMult)*1.01;
  Float_t cutVal(0.);
  Float_t minStep=1.;
  for(Int_t i=0; i<nBinsMult-1; ++i) {
    Float_t stepSize = 10.;
    cout << i << endl;
    while(true) {
//      std::cout << stepSize << "\t" << cutVal << "\t" << nBelow << "\t" << perBin*(i+1) << endl;
      TString cut="nTracks<";
      cut += cutVal;
      
      Int_t nBelow = mcTree->GetEntries(cut);
      
      if(nBelow < perBin*i + perBinLo) {
        cutVal += stepSize;
	stepSize *= 1.2;
	continue;
      }
      if(nBelow > perBin*i + perBinHi) {
        cutVal -= stepSize;
	stepSize /= 2.;
	if(stepSize>minStep/5.) continue;
      }

      multBins[i+1] = TMath::Nint(cutVal/minStep)*minStep;
      break;

    }
  }
  cout << "p bins:" << endl;
  for(Int_t i=0; i<nBinsP+1; ++i) {
    cout << pBins[i] << endl;
  }
  cout << "eta bins:" << endl;
  for(Int_t i=0; i<nBinsEta+1; ++i) {
    cout << etaBins[i] << endl;
  }
  cout << "nTrack bins:" << endl;
  for(Int_t i=0; i<nBinsMult+1; ++i) {
    cout << multBins[i] << endl;
  }

  TH3F*   dataHist1 = new TH3F("dataHist1",   "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*     mcHist1 = new TH3F("mcHist1",     "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F* weightHist1 = new TH3F("weightHist1", "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*   dataHist2 = new TH3F("dataHist2",   "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*     mcHist2 = new TH3F("mcHist2",     "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F* weightHist2 = new TH3F("weightHist2", "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*   dataHist3 = new TH3F("dataHist3",   "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*     mcHist3 = new TH3F("mcHist3",     "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F* weightHist3 = new TH3F("weightHist3", "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*   dataHist4 = new TH3F("dataHist4",   "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F*     mcHist4 = new TH3F("mcHist4",     "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);
  TH3F* weightHist4 = new TH3F("weightHist4", "", nBinsP, pBins, nBinsEta, etaBins, nBinsMult, multBins);

  for(Int_t i=0; i<nEntriesData; ++i) {
    dataTree->GetEntry(i);

    if(B_D0_B_CM_M<mB-2.5*sB || B_D0_B_CM_M>mB+2.5*sB) continue;

    B_D0_p_CM_P   = sqrt(B_D0_p_CM_PX*B_D0_p_CM_PX+B_D0_p_CM_PY*B_D0_p_CM_PY+B_D0_p_CM_PZ*B_D0_p_CM_PZ);
    B_D0_m_CM_P   = sqrt(B_D0_m_CM_PX*B_D0_m_CM_PX+B_D0_m_CM_PY*B_D0_m_CM_PY+B_D0_m_CM_PZ*B_D0_m_CM_PZ);
    B_D0_D0p_CM_P = sqrt(B_D0_D0p_CM_PX*B_D0_D0p_CM_PX+B_D0_D0p_CM_PY*B_D0_D0p_CM_PY+B_D0_D0p_CM_PZ*B_D0_D0p_CM_PZ);
    B_D0_D0m_CM_P = sqrt(B_D0_D0m_CM_PX*B_D0_D0m_CM_PX+B_D0_D0m_CM_PY*B_D0_D0m_CM_PY+B_D0_D0m_CM_PZ*B_D0_D0m_CM_PZ);

    if(pi_ID>0) {
      dataHist1->Fill(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      dataHist2->Fill(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else if(pi_ID<0) {
      dataHist2->Fill(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      dataHist1->Fill(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else {
      cout << "Bad pion ID" << endl;
    }

    dataHist4->Fill(B_D0_D0p_CM_P,0.5*log((B_D0_D0p_CM_P+B_D0_D0p_CM_PZ)/(B_D0_D0p_CM_P-B_D0_D0p_CM_PZ)),nTracks);
    dataHist3->Fill(B_D0_D0m_CM_P,0.5*log((B_D0_D0m_CM_P+B_D0_D0m_CM_PZ)/(B_D0_D0m_CM_P-B_D0_D0m_CM_PZ)),nTracks);
  }
  
  for(Int_t i=0; i<nEntriesMC; ++i) {
    mcTree->GetEntry(i);

    if(abs(B_TRUEID)!=511||abs(D_TRUEID)!=421||abs(K_TRUEID)!=321||abs(pi_TRUEID)!=211||abs(D0p_TRUEID)!=321||abs(D0m_TRUEID)!=321) continue;
    if(K_MC_MOTHER_KEY!=pi_MC_MOTHER_KEY||K_MC_MOTHER_KEY!=D0p_MC_GD_MOTHER_KEY||K_MC_MOTHER_KEY!=D0m_MC_GD_MOTHER_KEY) continue;

    B_D0_p_CM_P   = sqrt(B_D0_p_CM_PX*B_D0_p_CM_PX+B_D0_p_CM_PY*B_D0_p_CM_PY+B_D0_p_CM_PZ*B_D0_p_CM_PZ);
    B_D0_m_CM_P   = sqrt(B_D0_m_CM_PX*B_D0_m_CM_PX+B_D0_m_CM_PY*B_D0_m_CM_PY+B_D0_m_CM_PZ*B_D0_m_CM_PZ);
    B_D0_D0p_CM_P = sqrt(B_D0_D0p_CM_PX*B_D0_D0p_CM_PX+B_D0_D0p_CM_PY*B_D0_D0p_CM_PY+B_D0_D0p_CM_PZ*B_D0_D0p_CM_PZ);
    B_D0_D0m_CM_P = sqrt(B_D0_D0m_CM_PX*B_D0_D0m_CM_PX+B_D0_D0m_CM_PY*B_D0_D0m_CM_PY+B_D0_D0m_CM_PZ*B_D0_D0m_CM_PZ);

    if(pi_ID>0) {
      mcHist1->Fill(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      mcHist2->Fill(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else if(pi_ID<0) {
      mcHist2->Fill(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      mcHist1->Fill(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else {
      cout << "Bad pion ID" << endl;
    }

    mcHist4->Fill(B_D0_D0p_CM_P,0.5*log((B_D0_D0p_CM_P+B_D0_D0p_CM_PZ)/(B_D0_D0p_CM_P-B_D0_D0p_CM_PZ)),nTracks);
    mcHist3->Fill(B_D0_D0m_CM_P,0.5*log((B_D0_D0m_CM_P+B_D0_D0m_CM_PZ)/(B_D0_D0m_CM_P-B_D0_D0m_CM_PZ)),nTracks);
  }

  weightHist1->Divide(dataHist1,mcHist1);
  weightHist2->Divide(dataHist2,mcHist2);
  weightHist3->Divide(dataHist3,mcHist3);
  weightHist4->Divide(dataHist4,mcHist4);

  Double_t piWeight(0.), KWeight(0.), D0mWeight(0.), D0pWeight(0.);
  Double_t piWeightErr(0.), KWeightErr(0.), D0mWeightErr(0.), D0pWeightErr(0.);
  Int_t piBin(0.), KBin(0.), D0mBin(0.), D0pBin(0.);

  TFile * newFile = TFile::Open("D2KK_Bd2D0Kpi_trackWeight.root","recreate");
  TTree * newTree = mcTree->CloneTree(0);

  newTree->Branch("piWeight",      &piWeight,      "piWeight/D");
  newTree->Branch("KWeight",       &KWeight,       "KWeight/D");
  newTree->Branch("D0mWeight",    &D0mWeight,    "D0mWeight/D");
  newTree->Branch("D0pWeight",     &D0pWeight,     "D0pWeight/D");
  newTree->Branch("piWeightErr",   &piWeightErr,   "piWeightErr/D");
  newTree->Branch("KWeightErr",    &KWeightErr,    "KWeightErr/D");
  newTree->Branch("D0mWeightErr", &D0mWeightErr, "D0mWeightErr/D");
  newTree->Branch("D0pWeightErr",  &D0pWeightErr,  "D0pWeightErr/D");
 
  for(Int_t i=0; i<nEntriesMC; ++i) {
    mcTree->GetEntry(i);

    B_D0_p_CM_P   = sqrt(B_D0_p_CM_PX*B_D0_p_CM_PX+B_D0_p_CM_PY*B_D0_p_CM_PY+B_D0_p_CM_PZ*B_D0_p_CM_PZ);
    B_D0_m_CM_P   = sqrt(B_D0_m_CM_PX*B_D0_m_CM_PX+B_D0_m_CM_PY*B_D0_m_CM_PY+B_D0_m_CM_PZ*B_D0_m_CM_PZ);
    B_D0_D0p_CM_P = sqrt(B_D0_D0p_CM_PX*B_D0_D0p_CM_PX+B_D0_D0p_CM_PY*B_D0_D0p_CM_PY+B_D0_D0p_CM_PZ*B_D0_D0p_CM_PZ);
    B_D0_D0m_CM_P = sqrt(B_D0_D0m_CM_PX*B_D0_D0m_CM_PX+B_D0_D0m_CM_PY*B_D0_D0m_CM_PY+B_D0_D0m_CM_PZ*B_D0_D0m_CM_PZ);

    if(pi_ID>0) {
      piBin = mcHist1->FindBin(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      KBin  = mcHist2->FindBin(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else if(pi_ID<0) {
      KBin  = mcHist2->FindBin(B_D0_p_CM_P,0.5*log((B_D0_p_CM_P+B_D0_p_CM_PZ)/(B_D0_p_CM_P-B_D0_p_CM_PZ)),nTracks);
      piBin = mcHist1->FindBin(B_D0_m_CM_P,0.5*log((B_D0_m_CM_P+B_D0_m_CM_PZ)/(B_D0_m_CM_P-B_D0_m_CM_PZ)),nTracks);

    } else {
      cout << "Bad pion ID" << endl;
    }

    D0pBin = mcHist4->FindBin(B_D0_D0p_CM_P,0.5*log((B_D0_D0p_CM_P+B_D0_D0p_CM_PZ)/(B_D0_D0p_CM_P-B_D0_D0p_CM_PZ)),nTracks);
    D0mBin = mcHist3->FindBin(B_D0_D0m_CM_P,0.5*log((B_D0_D0m_CM_P+B_D0_D0m_CM_PZ)/(B_D0_D0m_CM_P-B_D0_D0m_CM_PZ)),nTracks);

    piWeight   = weightHist1->GetBinContent(piBin  );
    KWeight    = weightHist2->GetBinContent(KBin   );
    D0mWeight = weightHist3->GetBinContent(D0mBin);
    D0pWeight  = weightHist4->GetBinContent(D0pBin );

    piWeightErr   = weightHist1->GetBinError(piBin  );
    KWeightErr    = weightHist2->GetBinError(KBin   );
    D0mWeightErr = weightHist3->GetBinError(D0mBin);
    D0pWeightErr  = weightHist4->GetBinError(D0pBin );

    newTree->Fill();
  }

  newTree->AutoSave();

//  TCanvas c1;
//  weightHist1->Draw("LEGO");
//  c1.SaveAs("weights_pi_p_eta_nTracks.pdf");
}
