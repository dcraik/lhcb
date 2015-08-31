#define get_moments_cxx
#include "get_moments.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>

using namespace std;

// These are just the Legendre polynomials
//
double D000(double x) { return 1.; }
double D001(double x) { return x; }
double D002(double x) { return (3.*x*x - 1.)/2.; }
double D003(double x) { return (5.*x*x*x - 3.*x)/2.; }
double D004(double x) { return (35.*x*x*x*x - 30.*x*x + 3.)/8.; }

double Omega000(double cosTheta) { return 1.; }
double Omega001(double cosThetaL){ return D000(0.) * D001(cosThetaL); }
double Omega002(double cosThetaL){ return D000(0.) * D002(cosThetaL); }
double Omega003(double cosThetaL){ return D000(0.) * D003(cosThetaL); }
double Omega004(double cosThetaL){ return D000(0.) * D004(cosThetaL); }

void get_moments::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    double G000(0.);
    double G001(0.);
    double G002(0.);
    double G003(0.);
    double G004(0.);
    double eff(0.);

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        eff = 1. - cosThetaL*cosThetaL/2.;
        G000 += sWeight*Omega000(cosThetaL)/eff;
        G001 += sWeight*Omega001(cosThetaL)/eff;
        G002 += sWeight*Omega002(cosThetaL)/eff;
        G003 += sWeight*Omega003(cosThetaL)/eff;
        G004 += sWeight*Omega004(cosThetaL)/eff;
    }
    cout << "G000 \t 3*G001 \t 5*G002 \t 7*G003 \t 9*G004" << endl;
    cout << G000/G000 << "\t" << 3*G001/G000 << "\t" << 5*G002/G000 << "\t" << 7*G003/G000 << "\t" << 9*G004/G000 << endl;
    TFile * f = TFile::Open("moments.root", "RECREATE");
    TNtuple * tup = new TNtuple("coeffs", "coeffs", "G000:3G001:5G002:7G003:9G004");
    tup->Fill(G000/G000, 3.*G001/G000, 5.*G002/G000, 7.*G003/G000, 9.*G004/G000);
    tup->Write();
    f->Close();
}

int main( int argv, char * argc[] )
{
    TFile * f = TFile::Open(argc[1]);
    TTree * tree = (TTree*)f->Get("DecayTree");
    get_moments t(tree);
    t.Loop();
    return 1;
}

