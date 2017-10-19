void makeTableHist(TString file) {
  TFile* f = TFile::Open("tables2012L0HadronTOS_S20realET/"+file+".root");
  TTree* t = (TTree*)f->Get("eff");

  Int_t n = t->GetEntries();

  Double_t et(0.), effIn(0.), errIn(0.), effOut(0.), errOut(0.);

  t->SetBranchAddress("binStart", &et);
  t->SetBranchAddress("effIn",    &effIn);
  t->SetBranchAddress("errIn",    &errIn);
  t->SetBranchAddress("effOut",   &effOut);
  t->SetBranchAddress("errOut",   &errOut);

  Float_t* etBins = new Float_t[n+1];

  for(Int_t i=0; i<n; ++i) {
    t->GetEntry(i);
    etBins[i] = et;
  }
  etBins[n] = 100000;

  for(Int_t i=0; i<n+1; ++i) {
//    std::cout << etBins[i] << std::endl;
  }

  TH1F* hIn = new TH1F("hIn","",n,etBins);
  TH1F* hOut = new TH1F("hOut","",n,etBins);

  for(Int_t i=0; i<n; ++i) {
    t->GetEntry(i);
    hIn->SetBinContent(i+1,effIn);
    hIn->SetBinError(i+1,errIn);
    hOut->SetBinContent(i+1,effOut);
    hOut->SetBinError(i+1,errOut);
  }

  outfile = new TFile("tables/"+file+".root","RECREATE");
  hIn->SetName("inner");
  hIn->Write();
  hOut->SetName("outer");
  hOut->Write();
  outfile->Close();
}
