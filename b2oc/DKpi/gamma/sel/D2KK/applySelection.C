void applySelection(TString year, TString mag, TString part="")
{
   TString fname("B2D0Kpi_D02KK_");
   fname += year; fname += "_"; fname += mag;
   if(part!="") fname += "."; fname += part;
   fname += ".bdt.root";

   TFile* file = TFile::Open(fname);
   if(!file) {
      std::cerr << "File " << fname << " does not exist!" << std::endl;
      return;
   }
   TTree* tree = (TTree*)file->Get("DecayTree");
   if(!tree) {
      std::cerr << "File " << fname << " does not contain the tree DecayTree" << std::endl;
      return;
   }

   Float_t  Bd_K_E(0.),   Bd_K_PX(0.),   Bd_K_PY(0.),   Bd_K_PZ(0.);
   Float_t  Bd_pi_E(0.),  Bd_pi_PX(0.),  Bd_pi_PY(0.),  Bd_pi_PZ(0.);
   Float_t  Bd_Dp_E(0.),  Bd_Dp_PX(0.),  Bd_Dp_PY(0.),  Bd_Dp_PZ(0.);
   Float_t  Bd_Dm_E(0.),  Bd_Dm_PX(0.),  Bd_Dm_PY(0.),  Bd_Dm_PZ(0.);
                           
   Float_t  Bs_K_E(0.),   Bs_K_PX(0.),   Bs_K_PY(0.),   Bs_K_PZ(0.);
   Float_t  Bs_pi_E(0.),  Bs_pi_PX(0.),  Bs_pi_PY(0.),  Bs_pi_PZ(0.);
   Float_t  Bs_Dp_E(0.),  Bs_Dp_PX(0.),  Bs_Dp_PY(0.),  Bs_Dp_PZ(0.);
   Float_t  Bs_Dm_E(0.),  Bs_Dm_PX(0.),  Bs_Dm_PY(0.),  Bs_Dm_PZ(0.);

   Double_t D_MM(0.), B_D0_B_CM_M(0.), B_D0_B_CM_DIRA_OWNPV(0.), B_D0_B_CM_MINIPCHI2(0.);
   Double_t B_D0_B_CM_ENDVERTEX_CHI2(0.), B_D0_B_CM_ENDVERTEX_NDOF(0.);
   Double_t B_ENDVERTEX_CHI2(0.);
   Int_t B_ENDVERTEX_NDOF(0.);

   Float_t lab1_bag(0.);
   Bool_t B_Hlt2Topo2BodyBBDTDecision_TOS(0.), B_Hlt2Topo3BodyBBDTDecision_TOS(0.), B_Hlt2Topo4BodyBBDTDecision_TOS(0.);
   Bool_t B_L0HadronDecision_TOS(0.), B_L0Global_TIS(0.);

   Double_t D_ENDVERTEX_X(0.), D_ENDVERTEX_Y(0.), D_ENDVERTEX_Z(0.);
   Double_t B_ENDVERTEX_X(0.), B_ENDVERTEX_Y(0.), B_ENDVERTEX_Z(0.);

   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PE"         ,   &Bd_Dp_E ); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PX"         ,   &Bd_Dp_PX); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PY"         ,   &Bd_Dp_PY); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_PZ"         ,   &Bd_Dp_PZ); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_0_PE"       ,   &Bd_Dm_E ); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_0_PX"       ,   &Bd_Dm_PX); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_0_PY"       ,   &Bd_Dm_PY); 
   tree->SetBranchAddress("B_ConstB0Fit_D0_Kplus_0_PZ"       ,   &Bd_Dm_PZ); 

   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PE"  ,   &Bd_K_E  ); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PX"  ,   &Bd_K_PX ); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PY"  ,   &Bd_K_PY ); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_Kplus_PZ"  ,   &Bd_K_PZ ); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PE" ,   &Bd_pi_E ); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PX" ,   &Bd_pi_PX); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PY" ,   &Bd_pi_PY); 
   tree->SetBranchAddress("B_ConstB0Fit_Kst_892_0_piplus_PZ" ,   &Bd_pi_PZ);

   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_PE"         ,   &Bs_Dp_E ); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_PX"         ,   &Bs_Dp_PX); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_PY"         ,   &Bs_Dp_PY); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_PZ"         ,   &Bs_Dp_PZ); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_0_PE"       ,   &Bs_Dm_E ); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_0_PX"       ,   &Bs_Dm_PX); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_0_PY"       ,   &Bs_Dm_PY); 
   tree->SetBranchAddress("B_ConstBsFit_D0_Kplus_0_PZ"       ,   &Bs_Dm_PZ); 

   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_Kplus_PE"  ,   &Bs_K_E  ); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_Kplus_PX"  ,   &Bs_K_PX ); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_Kplus_PY"  ,   &Bs_K_PY ); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_Kplus_PZ"  ,   &Bs_K_PZ ); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_piplus_PE" ,   &Bs_pi_E ); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_piplus_PX" ,   &Bs_pi_PX); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_piplus_PY" ,   &Bs_pi_PY); 
   tree->SetBranchAddress("B_ConstBsFit_Kst_892_0_piplus_PZ" ,   &Bs_pi_PZ);

   tree->SetBranchAddress("D_MM"                             ,   &D_MM);
   tree->SetBranchAddress("B_D0_B_CM_M"                      ,   &B_D0_B_CM_M);
   tree->SetBranchAddress("B_D0_B_CM_DIRA_OWNPV"             ,   &B_D0_B_CM_DIRA_OWNPV);
   tree->SetBranchAddress("B_D0_B_CM_MINIPCHI2"              ,   &B_D0_B_CM_MINIPCHI2);
   tree->SetBranchAddress("B_D0_B_CM_ENDVERTEX_CHI2"         ,   &B_D0_B_CM_ENDVERTEX_CHI2);
   tree->SetBranchAddress("B_D0_B_CM_ENDVERTEX_NDOF"         ,   &B_D0_B_CM_ENDVERTEX_NDOF);
   tree->SetBranchAddress("B_ENDVERTEX_CHI2"                 ,   &B_ENDVERTEX_CHI2);
   tree->SetBranchAddress("B_ENDVERTEX_NDOF"                 ,   &B_ENDVERTEX_NDOF);
   tree->SetBranchAddress("lab1_bag"                         ,   &lab1_bag);

   tree->SetBranchAddress("B_Hlt2Topo2BodyBBDTDecision_TOS"  , &B_Hlt2Topo2BodyBBDTDecision_TOS);
   tree->SetBranchAddress("B_Hlt2Topo3BodyBBDTDecision_TOS"  , &B_Hlt2Topo3BodyBBDTDecision_TOS);
   tree->SetBranchAddress("B_Hlt2Topo4BodyBBDTDecision_TOS"  , &B_Hlt2Topo4BodyBBDTDecision_TOS);

   tree->SetBranchAddress("B_L0HadronDecision_TOS"           ,   &B_L0HadronDecision_TOS);
   tree->SetBranchAddress("B_L0Global_TIS"                   ,   &B_L0Global_TIS);

   tree->SetBranchAddress("D_ENDVERTEX_X"                    ,   &D_ENDVERTEX_X);
   tree->SetBranchAddress("D_ENDVERTEX_Y"                    ,   &D_ENDVERTEX_Y);
   tree->SetBranchAddress("D_ENDVERTEX_Z"                    ,   &D_ENDVERTEX_Z);
   tree->SetBranchAddress("B_ENDVERTEX_X"                    ,   &B_ENDVERTEX_X);
   tree->SetBranchAddress("B_ENDVERTEX_Y"                    ,   &B_ENDVERTEX_Y);
   tree->SetBranchAddress("B_ENDVERTEX_Z"                    ,   &B_ENDVERTEX_Z);

   TString newname("B2D0Kpi_D02KK_");
   newname += year; newname += "_"; newname += mag;
   newname += "_selBd_Dsidebands";
   if(part!="") newname += "."; newname += part;
   newname += ".root";

   TFile* newfile = new TFile(newname,"recreate");
   TTree* newtree = tree->CloneTree(0);

   Int_t n(tree->GetEntries());

   TLorentzVector Bd_CM_D0p_4mom;
   TLorentzVector Bd_CM_D0m_4mom;
   TLorentzVector Bd_CM_D0_4mom;
   TLorentzVector Bd_CM_K_4mom;
   TLorentzVector Bd_CM_Pi_4mom;
   TLorentzVector Bd_CM_DK;
   TLorentzVector Bd_CM_DPi;
   TLorentzVector Bd_CM_KPi;

   TLorentzVector Bs_CM_D0p_4mom;
   TLorentzVector Bs_CM_D0m_4mom;
   TLorentzVector Bs_CM_D0_4mom;
   TLorentzVector Bs_CM_K_4mom;
   TLorentzVector Bs_CM_Pi_4mom;
   TLorentzVector Bs_CM_DK;
   TLorentzVector Bs_CM_DPi;
   TLorentzVector Bs_CM_KPi;

   Double_t Bd_CM_m13Sq(-1);
   Double_t Bd_CM_m23Sq(-1);
   Double_t Bd_CM_m12Sq(-1);
   Double_t Bd_CM_m13(-1);
   Double_t Bd_CM_m23(-1);
   Double_t Bd_CM_m12(-1);

   Double_t Bs_CM_m13Sq(-1);
   Double_t Bs_CM_m23Sq(-1);
   Double_t Bs_CM_m12Sq(-1);
   Double_t Bs_CM_m13(-1);
   Double_t Bs_CM_m23(-1);
   Double_t Bs_CM_m12(-1);

   Double_t dataset = year.Atoi();

   TBranch *newbranch1a = newtree->Branch("Bd_CM_m13Sq",&Bd_CM_m13Sq);
   TBranch *newbranch1b = newtree->Branch("Bd_CM_m23Sq",&Bd_CM_m23Sq);
   TBranch *newbranch1c = newtree->Branch("Bd_CM_m12Sq",&Bd_CM_m12Sq);
   TBranch *newbranch1d = newtree->Branch("Bd_CM_m13",  &Bd_CM_m13);
   TBranch *newbranch1e = newtree->Branch("Bd_CM_m23",  &Bd_CM_m23);
   TBranch *newbranch1f = newtree->Branch("Bd_CM_m12",  &Bd_CM_m12);

   TBranch *newbranch1g = newtree->Branch("Bs_CM_m13Sq",&Bs_CM_m13Sq);
   TBranch *newbranch1h = newtree->Branch("Bs_CM_m23Sq",&Bs_CM_m23Sq);
   TBranch *newbranch1i = newtree->Branch("Bs_CM_m12Sq",&Bs_CM_m12Sq);
   TBranch *newbranch1j = newtree->Branch("Bs_CM_m13",  &Bs_CM_m13);
   TBranch *newbranch1k = newtree->Branch("Bs_CM_m23",  &Bs_CM_m23);
   TBranch *newbranch1l = newtree->Branch("Bs_CM_m12",  &Bs_CM_m12);

   TBranch *newbranch1m = newtree->Branch("dataset",    &dataset);

   for (Long64_t j=0; j<n;++j) {
      if(j%10000==0) {std::cout<<"Processing entry "<<j<<"..."<<std::endl;}
      tree->GetEntry(j);

//      if( D_MM>=1814 && D_MM<= 1914 ){//D0 mass
         if( B_D0_B_CM_M >= 5000 && B_D0_B_CM_M <= 6000 ){//B0 mass 
            if( B_D0_B_CM_DIRA_OWNPV > 0.99995 && B_D0_B_CM_MINIPCHI2 < 9 && B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF < 4){// && B_D0_B_CM_ENDVERTEX_CHI2/B_D0_B_CM_ENDVERTEX_NDOF < 4){//B0 
               if(lab1_bag>0.2) {//DO
                  if(B_Hlt2Topo2BodyBBDTDecision_TOS==1 || B_Hlt2Topo3BodyBBDTDecision_TOS==1 || B_Hlt2Topo4BodyBBDTDecision_TOS==1){//Trigger HLT
                     if(B_L0HadronDecision_TOS==1 || B_L0Global_TIS==1) {//Trigger L0
               	        if(sqrt((D_ENDVERTEX_X-B_ENDVERTEX_X)**2+(D_ENDVERTEX_Y-B_ENDVERTEX_Y)**2+(D_ENDVERTEX_Z-B_ENDVERTEX_Z)**2)>1.) {//D FD
                           //Bd and D0 mass constraints
               	           Bd_CM_D0p_4mom.SetPxPyPzE(    Bd_Dp_PX/1000, Bd_Dp_PY/1000, Bd_Dp_PZ/1000, Bd_Dp_E/1000);
               	           Bd_CM_D0m_4mom.SetPxPyPzE(    Bd_Dm_PX/1000, Bd_Dm_PY/1000, Bd_Dm_PZ/1000, Bd_Dm_E/1000);
               	           Bd_CM_K_4mom.SetPxPyPzE(      Bd_K_PX/1000,  Bd_K_PY/1000,  Bd_K_PZ/1000,  Bd_K_E/1000);
               	           Bd_CM_Pi_4mom.SetPxPyPzE(     Bd_pi_PX/1000, Bd_pi_PY/1000, Bd_pi_PZ/1000, Bd_pi_E/1000);

		           Bd_CM_D0_4mom = Bd_CM_D0p_4mom + Bd_CM_D0m_4mom;

               	           Bd_CM_DK  = Bd_CM_D0_4mom + Bd_CM_K_4mom;
               	           Bd_CM_DPi = Bd_CM_D0_4mom + Bd_CM_Pi_4mom;
               	           Bd_CM_KPi = Bd_CM_K_4mom  + Bd_CM_Pi_4mom;

                           //Bs and D0 mass constraints
               	           Bs_CM_D0p_4mom.SetPxPyPzE(    Bs_Dp_PX/1000, Bs_Dp_PY/1000, Bs_Dp_PZ/1000, Bs_Dp_E/1000);
               	           Bs_CM_D0m_4mom.SetPxPyPzE(    Bs_Dm_PX/1000, Bs_Dm_PY/1000, Bs_Dm_PZ/1000, Bs_Dm_E/1000);
               	           Bs_CM_K_4mom.SetPxPyPzE(       Bs_K_PX/1000,  Bs_K_PY/1000,  Bs_K_PZ/1000,  Bs_K_E/1000);
               	           Bs_CM_Pi_4mom.SetPxPyPzE(      Bs_pi_PX/1000, Bs_pi_PY/1000, Bs_pi_PZ/1000, Bs_pi_E/1000);

		           Bs_CM_D0_4mom = Bs_CM_D0p_4mom + Bs_CM_D0m_4mom;

               	           Bs_CM_DK  = Bs_CM_D0_4mom + Bs_CM_K_4mom;
               	           Bs_CM_DPi = Bs_CM_D0_4mom + Bs_CM_Pi_4mom;
               	           Bs_CM_KPi = Bs_CM_K_4mom  + Bs_CM_Pi_4mom;

              	           //Bd and D0 mass contraints			
              	           Bd_CM_m13Sq = Bd_CM_DK.M2();
              	           Bd_CM_m23Sq = Bd_CM_DPi.M2();
              	           Bd_CM_m12Sq = Bd_CM_KPi.M2();
              	           Bd_CM_m13 = Bd_CM_DK.M();
              	           Bd_CM_m23 = Bd_CM_DPi.M();
              	           Bd_CM_m12 = Bd_CM_KPi.M();

              	           //Bs and D0 mass contraints			
              	           Bs_CM_m13Sq = Bs_CM_DK.M2();
              	           Bs_CM_m23Sq = Bs_CM_DPi.M2();
              	           Bs_CM_m12Sq = Bs_CM_KPi.M2();
              	           Bs_CM_m13 = Bs_CM_DK.M();
              	           Bs_CM_m23 = Bs_CM_DPi.M();
              	           Bs_CM_m12 = Bs_CM_KPi.M();

			   //Bd specific section
             	           if((Bd_K_PX**2+Bd_K_PY**2+Bd_K_PZ**2)**.5<100000 && (Bd_pi_PX**2+Bd_pi_PY**2+Bd_pi_PZ**2)**.5<100000) {//max p
             	              if((Bd_Dp_PX**2+Bd_Dp_PY**2+Bd_Dp_PZ**2)**.5<100000 && (Bd_Dm_PX**2+Bd_Dm_PY**2+Bd_Dm_PZ**2)**.5<100000) {//max p
                                 newtree->Fill();
			      }
			   }
			   //End specific sections
			   //Note cuts on PID, D0 daughter PID, NN and vetoes are still to come
			}
		     }
               	  }
               }
	    }
//         }
      }
   }
   newtree->AutoSave();
   newfile->Close();
}
