{
TFile* f0 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/BuKee_offline_selected_030915.root");
//TFile* f1 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT.root");
TFile* f1 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel.root");
//TFile* f2 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/Kee/B2Kee_Strip21_data_folded.root");
TFile* f2 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/Kee/B2Kee_Strip21_data_folded_presel.root");

TCanvas c;

((TTree*)f0->Get("DecayTree"))->Draw("B_M","Psi_M*Psi_M>1.e6 && Psi_M*Psi_M<6.e6 && Kplus_ProbNNghost<0.3 && eminus_ProbNNghost<0.3 && eplus_ProbNNghost<0.3 && Kplus_ProbNNk>0.2 && Kplus_PIDe<0 && eminus_PIDe>3 && eplus_PIDe>3 &&  B_M>4880 && B_M<5700 &&  B_L0ElectronDecision_TOS==1 && B_Hlt1TrackAllL0Decision_TOS==1 && (B_Hlt2TopoE2BodyBBDTDecision_TOS == 1 || B_Hlt2TopoE3BodyBBDTDecision_TOS == 1 || B_Hlt2Topo2BodyBBDTDecision_TOS == 1 || B_Hlt2Topo3BodyBBDTDecision_TOS == 1)");
c.SaveAs("dataset0.pdf");
//
((TTree*)f1->Get("DecayTree"))->Draw("B_M","JPs_M*JPs_M>1.e6 && JPs_M*JPs_M<6.e6 && E1_TRACK_GhostProb < 0.3 && E2_TRACK_GhostProb < 0.3 && K_TRACK_GhostProb < 0.3&& E1_PIDe > 3.0 && E2_PIDe > 3.0 && K_ProbNNk > 0.2 && K_PIDe < 0 && B_M > 4880 && B_M < 5700 && (B_L0ElectronDecision_TOS==1 || B_L0HadronDecision_TOS==1 || ((B_L0ElectronDecision_TIS==1) || (B_L0HadronDecision_TIS==1) || (B_L0MuonDecision_TIS==1) || (B_L0PhotonDecision_TIS))) && B_Hlt1TrackAllL0Decision_TOS==1 && (B_Hlt2TopoE2BodyBBDTDecision_TOS == 1 || B_Hlt2TopoE3BodyBBDTDecision_TOS == 1 || B_Hlt2Topo2BodyBBDTDecision_TOS == 1 || B_Hlt2Topo3BodyBBDTDecision_TOS == 1)");
c.SaveAs("dataset1c.pdf");
//
((TTree*)f2->Get("DecayTree"))->Draw("B_plus_M","J_psi_1S_M*J_psi_1S_M > 1.0e6 && J_psi_1S_M*J_psi_1S_M < 6.0e6 && K_Kst_TRACK_GhostProb<0.3 && e_plus_TRACK_GhostProb<0.3 && e_minus_TRACK_GhostProb<0.3 && K_Kst_ProbNNk>0.2 && K_Kst_PIDe<0 && e_plus_PIDe>3 && e_minus_PIDe>3 && B_plus_M>4880 && B_plus_M<5700 && (B_plus_L0ElectronDecision_TOS==1 || B_plus_L0HadronDecision_TOS==1 || ((B_plus_L0ElectronDecision_TIS==1) || (B_plus_L0HadronDecision_TIS==1) || (B_plus_L0MuonDecision_TIS==1) || (B_plus_L0PhotonDecision_TIS))) && B_plus_Hlt1TrackAllL0Decision_TOS==1 && ((B_plus_Hlt2TopoE3BodyBBDTDecision_TOS==1) || (B_plus_Hlt2TopoE2BodyBBDTDecision_TOS==1) || (B_plus_Hlt2Topo3BodyBBDTDecision_TOS==1) || (B_plus_Hlt2Topo2BodyBBDTDecision_TOS==1))");
c.SaveAs("dataset2c.pdf");
//
}
