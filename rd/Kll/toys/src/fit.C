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

void fit(Int_t i, Int_t n, Double_t vG000=0.5, Double_t vG001=0.5, Double_t vG002=0.5, Double_t vG003=0, Double_t vG004=0){

    TString nameStr=""; nameStr+=i; nameStr+="_"; nameStr+=n; nameStr+="_"; nameStr+=vG000; nameStr+="_"; nameStr+=vG001; nameStr+="_"; nameStr+=vG002; nameStr+="_"; nameStr+=vG003; nameStr+="_"; nameStr+=vG004;

    RooRealVar * cosTheta = new RooRealVar("cosTheta", "cosTheta", -1., 1.);
    RooRealVar * G000 = new RooRealVar("G000", "G000", vG000, -1., 1.);
    RooRealVar * G001 = new RooRealVar("G001", "G001", vG001, -1., 1.);
    RooRealVar * G002 = new RooRealVar("G002", "G002", vG002, -1., 1.);
    RooRealVar * G003 = new RooRealVar("G003", "G003", vG003);// -1., 1.);
    RooRealVar * G004 = new RooRealVar("G004", "G004", vG004);//, -1., 1.);
    RooRealVar * a    = new RooRealVar("a",    "a",    -0.5);
    RooRealVar * b    = new RooRealVar("b",    "b",     0.0);

    RooB2Kll * pdf = new RooB2Kll("pdf", "pdf", *cosTheta, *G000, *G001, *G002, *G003, *G004, *a, *b);

    //RooFormulaVar eff("eff","1 + @0*@0*@1 + @0*@0*@0*@0*@2",RooArgList(cosTheta,a,b));

    RooRandom::randomGenerator()->SetSeed(i);

    RooDataSet * data = pdf->generate(RooArgSet(*cosTheta), n);
    RooFitResult * fitresult = pdf->fitTo(*data,Save(kTRUE), Minos(kFALSE), NumCPU(4));
    
    RooPlot * frame = cosTheta->frame();
    data->plotOn(frame);
    pdf->plotOn(frame);

    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TCanvas* c = new TCanvas("fit","fit", 800, 800);
    frame->Draw();
    c->SaveAs("fits/"+nameStr+".png");
    c->SaveAs("fits/"+nameStr+".pdf");

    data->write("toys/toy_"+nameStr+".dat");
}

int main(int argc, char** argv)
{
    Int_t i(1), n(1000);
    Double_t g0(0.5), g1(0.5), g2(0.5), g3(0.), g4(0.);

    if(argc>1) {
        i = atoi(argv[1]);
        if(argc>2) {
            n = atoi(argv[2]);
            if(argc>5) {
                g0 = atof(argv[3]);
                g1 = atof(argv[4]);
                g2 = atof(argv[5]);
                if(argc>6) {
                    g3 = atof(argv[6]);
                    if(argc>7) {
		        g4 = atof(argv[7]);
                    }
                }
            }
        }
    }
    fit(i,n,g0,g1,g2,g3,g4);
    return 0;
}

