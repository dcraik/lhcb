void makeBootstrappedSamples() {

	/Disk/ecdf-nfs-ppe/lhcb/dcraik/kmumu_gen_addVars.root

	treegen->SetBranchAddress("qSq",       &qSq);
	treegen->SetBranchAddress("costhetal", &costhetal);
	treegen->SetBranchAddress("mKmu_D",    &mKmu_D);
	treegen->SetBranchAddress("mKmu_Psi",  &mKmu_Psi);
}
