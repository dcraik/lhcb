#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"
#include "RooRandom.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "RooBDecay.h"
#include "RooAddModel.h"
#include "RooGaussModel.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooB2Kll.h"
#include "RooFitResult.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace RooFit;
using namespace std;

TRandom3 r;

TFile * file(0);
TTree * DecayTree(0);
RooRealVar * cosTheta(0);
RooRealVar * sWeight(0);
RooDataSet * data(0);

char namestr[30];

void fit(Int_t i, Double_t va=-0.43, Double_t vb=0.){
    sprintf(namestr,"%d_%.2f_%.2f",i,va,vb);

    Double_t g1, g2;
    g1 = 0.1*r.Rndm() - 0.05;
    g2 =     r.Rndm() - 0.5;

    RooRealVar G000("G000", "G000",  0.5);
    RooRealVar G001("G001", "G001",  g1, -5., 5.);
    RooRealVar G002("G002", "G002",  g2, -5., 5.);
    RooRealVar G003("G003", "G003", 0.0);//, -1., 1.);
    RooRealVar G004("G004", "G004", 0.0);//, -1., 1.);
    RooRealVar    a("a",    "a",    va);
    RooRealVar    b("b",    "b",    vb);
    RooRealVar    n("n",    "n",    1000., -100., 5000.);

    RooB2Kll pdf("pdf", "pdf", *cosTheta, G000, G001, G002, G003, G004, a, b, n);

    RooFitResult * fitresult = pdf.fitTo(*data,Save(kTRUE), Minos(kFALSE), NumCPU(4), SumW2Error(kTRUE));
    
    cout << i << "\t" << a.getVal() << "\t" << b.getVal() << "\t" << fitresult->minNll() << "\t" << fitresult->status() << "\t" 
         << G000.getVal() << "\t" << G001.getVal() << "\t" << G002.getVal() << "\t" << G003.getVal() << "\t" << G004.getVal() << endl;
    
    RooPlot * frame = cosTheta->frame();
    data->plotOn(frame);
    pdf.plotOn(frame);

    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TCanvas c("fit","fit", 800, 800);
    frame->Draw();
    c.SaveAs("singleFits/fit"+TString(namestr)+".png");
    c.SaveAs("singleFits/fit"+TString(namestr)+".pdf");
}

int main(int argc, char** argv)
{
    file = TFile::Open("~/lhcb/rd/Kll/tuples/fromPatrick/Kmm_Q17_reduced.root");
    DecayTree = dynamic_cast<TTree*>(file->Get("finalTree_KMuMu"));
    cosTheta = new RooRealVar("cosThetaL", "cosThetaL", -1., 1.);
    sWeight  = new RooRealVar("sWeight",   "sWeight",   -5., 5.);
    data = new RooDataSet("data", "", DecayTree, RooArgList(*cosTheta,*sWeight), 0, "sWeight");

    Double_t a(-0.43), b(0.);
    Int_t i(0);

    if(argc>1) {
        i = atoi(argv[1]);
        if(argc>3) {
            a = atof(argv[2]);
            b = atof(argv[3]);
	}
    }
    
    r.SetSeed(i);
    fit(i,a,b);

    return 0;
}

