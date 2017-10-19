void applyBsWindow()
{
	TFile* file = TFile::Open("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs.root");
	TTree* tree = (TTree*)file->Get("DecayTree");

	TFile* newfile = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25.root","recreate");
	TTree* newtree = new TTree("DecayTree","DecayTree");

	TFile* newfileA = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_a.root","recreate");
	TTree* newtreeA = new TTree("DecayTree","DecayTree");

	TFile* newfileB = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_b.root","recreate");
	TTree* newtreeB = new TTree("DecayTree","DecayTree");

	TFile* newfileC = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_c.root","recreate");
	TTree* newtreeC = new TTree("DecayTree","DecayTree");

	TFile* newfileD = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_d.root","recreate");
	TTree* newtreeD = new TTree("DecayTree","DecayTree");

	TFile* newfileE = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_e.root","recreate");
	TTree* newtreeE = new TTree("DecayTree","DecayTree");

	TFile* newfileZ = new TFile("B2D0Kpi_D02piK_selBd_Dsignal_Dst25_vetoes_PID3_NND2Kpi_addIMs_Bs25_overlapTest.root","recreate");
	TTree* newtreeZ = new TTree("DecayTree","DecayTree");

   Double_t mBs(5278.34+87.19), sBs(11.94);

   Double_t Bs_CM_m23Sq(0.), Bs_CM_m12Sq(0.), Bs_CM_m13Sq(0.), B_D0_B_CM_M(0.);
   Int_t K_ID(0);
   Float_t NN(0.);

   //load Laura++ library
   gSystem->Load( "libTree" );
   gSystem->Load( "libHist" );
   gSystem->Load( "libMatrix" );
   gSystem->Load( "libEG" );
   gSystem->Load( "/data/lhcb/phrkbf/Laura-gamma/Laura++/lib/libLaura++.so" );

   LauDaughters* daughters = new LauDaughters( "B_s0", "D0_bar", "K-", "pi+", true );
   LauKinematics* kinematics = daughters->getKinematics();

   Double_t m12Sq(0.);
   Double_t m13(0.);
   Double_t m23(0.);
   Double_t m12(0.);
   Double_t c13(0.);
   Double_t c23(0.);
   Double_t c12(0.);
   Double_t mPrime(0.);
   Double_t thetaPrime(0.);
   Int_t charge(0);

   UInt_t run(0);
   ULong64_t evt(0);
   
   Int_t iExpt(0);
   //these were numbered in KpiD before - from here we fix them to DKpi
   TBranch *newbranch1a = newtree->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch1b = newtree->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch1c = newtree->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch1d = newtree->Branch("m13",        &m13);//Dpi
   TBranch *newbranch1e = newtree->Branch("m23",        &m23);//Kpi
   TBranch *newbranch1f = newtree->Branch("m12",        &m12);//DK
   TBranch *newbranch1g = newtree->Branch("c13",        &c13);
   TBranch *newbranch1h = newtree->Branch("c23",        &c23);
   TBranch *newbranch1i = newtree->Branch("c12",        &c12);
   TBranch *newbranch1j = newtree->Branch("mPrime",     &mPrime);
   TBranch *newbranch1k = newtree->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch1l = newtree->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch1m = newtree->Branch("iExpt",      &iExpt);
   TBranch *newbranch1n = newtree->Branch("charge",     &charge);

   TBranch *newbranch2a = newtreeA->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch2b = newtreeA->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch2c = newtreeA->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch2d = newtreeA->Branch("m13",        &m13);//Dpi
   TBranch *newbranch2e = newtreeA->Branch("m23",        &m23);//Kpi
   TBranch *newbranch2f = newtreeA->Branch("m12",        &m12);//DK
   TBranch *newbranch2g = newtreeA->Branch("c13",        &c13);
   TBranch *newbranch2h = newtreeA->Branch("c23",        &c23);
   TBranch *newbranch2i = newtreeA->Branch("c12",        &c12);
   TBranch *newbranch2j = newtreeA->Branch("mPrime",     &mPrime);
   TBranch *newbranch2k = newtreeA->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch2l = newtreeA->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch2m = newtreeA->Branch("iExpt",      &iExpt);
   TBranch *newbranch2n = newtreeA->Branch("charge",     &charge);

   TBranch *newbranch3a = newtreeB->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch3b = newtreeB->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch3c = newtreeB->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch3d = newtreeB->Branch("m13",        &m13);//Dpi
   TBranch *newbranch3e = newtreeB->Branch("m23",        &m23);//Kpi
   TBranch *newbranch3f = newtreeB->Branch("m12",        &m12);//DK
   TBranch *newbranch3g = newtreeB->Branch("c13",        &c13);
   TBranch *newbranch3h = newtreeB->Branch("c23",        &c23);
   TBranch *newbranch3i = newtreeB->Branch("c12",        &c12);
   TBranch *newbranch3j = newtreeB->Branch("mPrime",     &mPrime);
   TBranch *newbranch3k = newtreeB->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch3l = newtreeB->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch3m = newtreeB->Branch("iExpt",      &iExpt);
   TBranch *newbranch3n = newtreeB->Branch("charge",     &charge);

   TBranch *newbranch4a = newtreeC->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch4b = newtreeC->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch4c = newtreeC->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch4d = newtreeC->Branch("m13",        &m13);//Dpi
   TBranch *newbranch4e = newtreeC->Branch("m23",        &m23);//Kpi
   TBranch *newbranch4f = newtreeC->Branch("m12",        &m12);//DK
   TBranch *newbranch4g = newtreeC->Branch("c13",        &c13);
   TBranch *newbranch4h = newtreeC->Branch("c23",        &c23);
   TBranch *newbranch4i = newtreeC->Branch("c12",        &c12);
   TBranch *newbranch4j = newtreeC->Branch("mPrime",     &mPrime);
   TBranch *newbranch4k = newtreeC->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch4l = newtreeC->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch4m = newtreeC->Branch("iExpt",      &iExpt);
   TBranch *newbranch4n = newtreeC->Branch("charge",     &charge);

   TBranch *newbranch5a = newtreeD->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch5b = newtreeD->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch5c = newtreeD->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch5d = newtreeD->Branch("m13",        &m13);//Dpi
   TBranch *newbranch5e = newtreeD->Branch("m23",        &m23);//Kpi
   TBranch *newbranch5f = newtreeD->Branch("m12",        &m12);//DK
   TBranch *newbranch5g = newtreeD->Branch("c13",        &c13);
   TBranch *newbranch5h = newtreeD->Branch("c23",        &c23);
   TBranch *newbranch5i = newtreeD->Branch("c12",        &c12);
   TBranch *newbranch5j = newtreeD->Branch("mPrime",     &mPrime);
   TBranch *newbranch5k = newtreeD->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch5l = newtreeD->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch5m = newtreeD->Branch("iExpt",      &iExpt);
   TBranch *newbranch5n = newtreeD->Branch("charge",     &charge);

   TBranch *newbranch6a = newtreeE->Branch("m13Sq",      &Bs_CM_m23Sq);//Dpi
   TBranch *newbranch6b = newtreeE->Branch("m23Sq",      &Bs_CM_m12Sq);//Kpi
   TBranch *newbranch6c = newtreeE->Branch("m12Sq",      &m12Sq);//DK
   TBranch *newbranch6d = newtreeE->Branch("m13",        &m13);//Dpi
   TBranch *newbranch6e = newtreeE->Branch("m23",        &m23);//Kpi
   TBranch *newbranch6f = newtreeE->Branch("m12",        &m12);//DK
   TBranch *newbranch6g = newtreeE->Branch("c13",        &c13);
   TBranch *newbranch6h = newtreeE->Branch("c23",        &c23);
   TBranch *newbranch6i = newtreeE->Branch("c12",        &c12);
   TBranch *newbranch6j = newtreeE->Branch("mPrime",     &mPrime);
   TBranch *newbranch6k = newtreeE->Branch("thetaPrime", &thetaPrime);
   TBranch *newbranch6l = newtreeE->Branch("B_M",        &B_D0_B_CM_M);
   TBranch *newbranch6m = newtreeE->Branch("iExpt",      &iExpt);
   TBranch *newbranch6n = newtreeE->Branch("charge",     &charge);

   TBranch *newbranch7m = newtreeZ->Branch("runNumber",     &run);
   TBranch *newbranch7n = newtreeZ->Branch("eventNumber",   &evt);
   TBranch *newbranch7l = newtreeZ->Branch("B_M",           &B_D0_B_CM_M);
   TBranch *newbranch7n = newtreeZ->Branch("charge",        &charge);

   tree->SetBranchAddress("Bs_CM_m23Sq", &Bs_CM_m23Sq);
   tree->SetBranchAddress("Bs_CM_m12Sq", &Bs_CM_m12Sq);
   tree->SetBranchAddress("Bs_CM_m13Sq", &Bs_CM_m13Sq);
   tree->SetBranchAddress("B_D0_B_CM_M", &B_D0_B_CM_M);
   tree->SetBranchAddress("K_ID", &K_ID);
   tree->SetBranchAddress("NN", &NN);
   tree->SetBranchAddress("runNumber", &run);
   tree->SetBranchAddress("eventNumber", &evt);

   Long64_t nentries = tree->GetEntries();

   for (Int_t i=0; i<nentries;++i) {
      tree->GetEntry(i);

      if(B_D0_B_CM_M > mBs - 2.5*sBs && B_D0_B_CM_M < mBs + 2.5*sBs) {
         kinematics->updateKinematics(Bs_CM_m23Sq,Bs_CM_m12Sq);//fixed to DKpi here

	 charge = K_ID/321;

         m12Sq = kinematics->getm12Sq();

         m13 = kinematics->getm13();
         m23 = kinematics->getm23();
         m12 = kinematics->getm12();

         c13 = kinematics->getc13();
         c23 = kinematics->getc23();
         c12 = kinematics->getc12();

         mPrime = kinematics->getmPrime();
         thetaPrime = kinematics->getThetaPrime();

	 if(NN>-0.8) newtreeZ->Fill();

         if(m23 > 1.835 && m23 < 1.880) {//Kpi is now 23
            cout << m23 << endl;
            continue;
         }
         if(m13 > 2.00778 && m13 < 2.01278) {//Dpi is now 13
           cout << m13 << endl;
           continue;
         }

         if(NN>0.23) newtree->Fill();
         if(NN>-0.80 && NN<= 0.00) newtreeA->Fill();
         if(NN> 0.00 && NN<= 0.50) newtreeB->Fill();
         if(NN> 0.50 && NN<= 0.72) newtreeC->Fill();
         if(NN> 0.72 && NN<= 0.81) newtreeD->Fill();
         if(NN> 0.81 && NN<= 1.00) newtreeE->Fill();
      }
   }
   newtree->AutoSave();
   newfile->Close();
   newtreeA->AutoSave();
   newfileA->Close();
   newtreeB->AutoSave();
   newfileB->Close();
   newtreeC->AutoSave();
   newfileC->Close();
   newtreeD->AutoSave();
   newfileD->Close();
   newtreeE->AutoSave();
   newfileE->Close();

   newtreeZ->AutoSave();
   newfileZ->Close();
}
