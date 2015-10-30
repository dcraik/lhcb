//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 16:43:29 2015 by ROOT version 5.34/10
// from TTree dataNTuple/All data
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef get_moments_B2Kstll_h
#define get_moments_B2Kstll_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TNtuple.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class get_moments_B2Kstll {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        cosThetaL;
   Double_t        cosThetaK;
   Double_t        phi;
   Double_t        sWeight;
   Double_t        B_M;

   // List of branches
   TBranch        *b_cosThetaL;   //!
   TBranch        *b_cosThetaK;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_sWeight;   //!
   TBranch        *b_B_M;   //!

   get_moments_B2Kstll(TTree *tree=0, Double_t effA=-0.5, Double_t effB=0.0, TH1* dVetoHist=0, TH1* psiVetoHist=0,
		      Bool_t useSWeights=true, Double_t sigWeight=1.0, Double_t bkgWeight=-1.0, Double_t sigMin=5170., Double_t sigMax=5970., Double_t bkgMin=5170., Double_t bkgMax=5970., Bool_t calcErrors=true);
   virtual ~get_moments_B2Kstll();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString filename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

private:
   virtual void     GetMoments(bool bootstrap=false);
   virtual void     SetIndicesForBootstrapping(int seed=123456);
   virtual void     GetUncertainties();

   Long64_t nentries;
   std::vector<Long64_t> indices;

   TFile *outFile;
   TNtuple *tup;
   TNtuple *tup_total;
   TNtuple *tup_error;
   TNtuple *tup_bootstrapped;

   //efficiency coefficients
   Double_t a, b;

   //veto histograms
   TH1 *dVeto, *psiVeto;

   //if useSWeights is false then we perform a background subtraction
   Bool_t useSWeights_;
   Double_t sigWeight_;
   Double_t bkgWeight_;
   Double_t sigMin_;
   Double_t sigMax_;
   Double_t bkgMin_;
   Double_t bkgMax_;

   Bool_t calcErrors_;
};

#endif

#ifdef get_moments_B2Kstll_cxx
get_moments_B2Kstll::get_moments_B2Kstll(TTree *tree, Double_t effA, Double_t effB, TH1* dVetoHist, TH1* psiVetoHist, 
		                         Bool_t useSWeights, Double_t sigWeight, Double_t bkgWeight, Double_t sigMin, Double_t sigMax, Double_t bkgMin, Double_t bkgMax, Bool_t calcErrors) 
    : fChain(0), cosThetaL(1.), cosThetaK(1.), phi(0.), sWeight(1.), a(effA), b(effB), dVeto(dVetoHist), psiVeto(psiVetoHist), useSWeights_(useSWeights), sigWeight_(sigWeight), bkgWeight_(bkgWeight), sigMin_(sigMin), sigMax_(sigMax), bkgMin_(bkgMin), bkgMax_(bkgMax), calcErrors_(calcErrors)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("dataNTuple",tree);

   }
   Init(tree);
}

get_moments_B2Kstll::~get_moments_B2Kstll()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t get_moments_B2Kstll::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t get_moments_B2Kstll::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void get_moments_B2Kstll::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   Int_t n = fChain->GetListOfBranches()->GetEntries();
   TBranch* branch(0);

   for(Int_t i=0; i<n; ++i) {
       branch = dynamic_cast<TBranch*>(fChain->GetListOfBranches()->At(i));
       TString name = branch->GetName();

       if( !name.CompareTo("cosThetaL") ) {
           fChain->SetBranchAddress("cosThetaL", &cosThetaL, &b_cosThetaL);
       } else if( !name.CompareTo("cosThetaK") ) {
           fChain->SetBranchAddress("cosThetaK", &cosThetaK, &b_cosThetaK);
       } else if( !name.CompareTo("phi") ) {
           fChain->SetBranchAddress("phi", &phi, &b_phi);
       } else if( !name.CompareTo("sWeight") ) {
           fChain->SetBranchAddress("sWeight", &sWeight, &b_sWeight);
       } else if( !name.CompareTo("Bplus_M") ) {
           fChain->SetBranchAddress("Bplus_M", &B_M, &b_B_M);
       }
   }


   //setup the indices for bootstrapping
   nentries = fChain->GetEntriesFast();
   indices.reserve(nentries);

   Notify();
}

Bool_t get_moments_B2Kstll::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void get_moments_B2Kstll::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t get_moments_B2Kstll::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef get_moments_B2Kstll_cxx
