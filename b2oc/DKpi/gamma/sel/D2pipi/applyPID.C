void applyPID()
{
   TString fname("B2D0Kpi_D02pipi_selBd_NND2pipi.root");

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

   Double_t D0p_ProbNNk(0.), D0p_ProbNNpi(0.), D0m_ProbNNk(0.), D0m_ProbNNpi(0.);
   Double_t K_ProbNNk(0.),   K_ProbNNpi(0.),   pi_ProbNNk(0.),  pi_ProbNNpi(0.);
   Double_t B_D0_B_CM_M(0.), Bd_CM_m23(0.), B_D0_D_CM_M(0.);
                           
   tree->SetBranchAddress("D0p_ProbNNk"        ,   &D0p_ProbNNk ); 
   tree->SetBranchAddress("D0p_ProbNNpi"       ,   &D0p_ProbNNpi); 
   tree->SetBranchAddress("D0m_ProbNNk"        ,   &D0m_ProbNNk ); 
   tree->SetBranchAddress("D0m_ProbNNpi"       ,   &D0m_ProbNNpi); 

   tree->SetBranchAddress("K_ProbNNk"          ,   &K_ProbNNk ); 
   tree->SetBranchAddress("K_ProbNNpi"         ,   &K_ProbNNpi); 
   tree->SetBranchAddress("pi_ProbNNk"         ,   &pi_ProbNNk ); 
   tree->SetBranchAddress("pi_ProbNNpi"        ,   &pi_ProbNNpi);

   tree->SetBranchAddress("B_D0_B_CM_M"        ,   &B_D0_B_CM_M);

   tree->SetBranchAddress("Bd_CM_m23"          ,   &Bd_CM_m23);
   tree->SetBranchAddress("B_D0_D_CM_M"        ,   &B_D0_D_CM_M);

   TString newname("B2D0Kpi_D02pipi_selBd_Dst25_PID3_NND2pipi.root");

   TFile* newfile = new TFile(newname,"recreate");
   TTree* newtree = tree->CloneTree(0);

   Int_t n(tree->GetEntries());

   for (Long64_t j=0; j<n;++j) {
      if(j%10000==0) {std::cout<<"Processing entry "<<j<<"..."<<std::endl;}
      tree->GetEntry(j);

      if( K_ProbNNk*(1-K_ProbNNpi)>0.3 && pi_ProbNNpi*(1-pi_ProbNNk)>0.2 ){//bachelor PID
         if( D0p_ProbNNpi*(1-D0p_ProbNNk)>0.1 && D0m_ProbNNpi*(1-D0m_ProbNNk)>0.1 ) {//D0 PID
	    if( Bd_CM_m23 - B_D0_D_CM_M/1000. < (145.4-2.5)/1000. || Bd_CM_m23 - B_D0_D_CM_M/1000. > (145.4+2.5)/1000.) {//Dst veto
               newtree->Fill();
	    }
         }
      }
   }
   newtree->AutoSave();
   newfile->Close();
}
