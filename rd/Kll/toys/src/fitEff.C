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
#include "TCut.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "RooBDecay.h"
#include "RooAddModel.h"
#include "RooGaussModel.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooB2Kll.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace RooFit;
using namespace std;

void fit(Int_t i, Int_t n, Double_t genG1=0.5, Double_t genG2=0.5, Double_t genA=-0.5, Double_t genB=0.0, Double_t fitG1=0, Double_t fitG2=0){

    RooRealVar * cosTheta = new RooRealVar("cosTheta", "cosTheta", -1., 1.);
    RooRealVar * G000 = new RooRealVar("G000", "G000", 0.5);
    RooRealVar * G001 = new RooRealVar("G001", "G001", genG1);
    RooRealVar * G002 = new RooRealVar("G002", "G002", genG2);
    RooRealVar * G003 = new RooRealVar("G003", "G003", 0.);
    RooRealVar * G004 = new RooRealVar("G004", "G004", 0.);
    RooRealVar * a    = new RooRealVar("a",    "a",    genA, -2.0, 2.0);
    RooRealVar * b    = new RooRealVar("b",    "b",    genB, -2.0, 2.0);

    RooB2Kll * pdf = new RooB2Kll("pdf", "pdf", *cosTheta, *G000, *G001, *G002, *G003, *G004, *a, *b);

    RooRandom::randomGenerator()->SetSeed(i);

    RooDataSet * data = pdf->generate(RooArgSet(*cosTheta), n);

    G001->setVal(fitG1);
    G002->setVal(fitG2);

    RooFitResult * fitresult = pdf->fitTo(*data,Save(kTRUE), Minos(kFALSE), NumCPU(4));
    
    RooPlot * frame = cosTheta->frame();
    data->plotOn(frame);
    pdf->plotOn(frame);

    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TCanvas* c = new TCanvas("fit","fit", 800, 800);
    frame->Draw();
    c->SaveAs("fit.png");
    c->SaveAs("fit.pdf");

    //data->write("toys/toy_"+nameStr+".dat");
}

int main(int argc, char** argv)
{
    Int_t i(1), n(1000);
    Double_t genG1(0.0), genG2(-0.5), genA(-0.5), genB(0.), fitG1(0.0), fitG2(-0.5);

    if(argc>1) {
        i = atoi(argv[1]);
        if(argc>2) {
            n = atoi(argv[2]);
            if(argc>8) {
                genG1 = atof(argv[3]);
                genG2 = atof(argv[4]);
                genA  = atof(argv[5]);
                genB  = atof(argv[6]);
                fitG1 = atof(argv[7]);
                fitG2 = atof(argv[8]);
            }
        }
    }
    fit(i,n,genG1,genG2,genA,genB,fitG1,fitG2);
    return 0;
}

