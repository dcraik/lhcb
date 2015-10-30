#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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
#include "RooBDecay.h"
#include "RooAddModel.h"
#include "RooGaussModel.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooB2Kll.h"
#include "RooFitResult.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace RooFit;
using namespace std;


void fit(Int_t i, Double_t va=-0.43, Double_t vb=0.){
    char namestr[30];

    Double_t sigMin(5242.4), sigMax(5324.9), bkgMin(5400), bkgMax(5970);

    TString sigStr="Bplus_M > "; sigStr+=sigMin; sigStr+=" && Bplus_M < "; sigStr+=sigMax;
    TString bkgStr="Bplus_M > "; bkgStr+=bkgMin; bkgStr+=" && Bplus_M < "; bkgStr+=bkgMax;

    std::cout << sigStr << std::endl << bkgStr << std::endl;

    TFile * file(0);
    TTree * DecayTree(0);
    RooRealVar * cosTheta(0);
    RooRealVar * Bplus_M(0);
    RooDataSet * data(0);
    
    file = TFile::Open("~/lhcb/rd/Kll/tuples/fromPatrick/Kmm_Q17_reduced.root");
    DecayTree = dynamic_cast<TTree*>(file->Get("finalTree_KMuMu"));
    cosTheta = new RooRealVar("cosThetaL", "cosThetaL", -1., 1.);
    Bplus_M  = new RooRealVar("Bplus_M", "Bplus_M", 5170., 5970.);
    data = new RooDataSet("data", "", DecayTree, RooArgList(*cosTheta,*Bplus_M));

    RooDataSet* data1 = dynamic_cast<RooDataSet*>(data->reduce(sigStr));
    RooDataSet* data2 = dynamic_cast<RooDataSet*>(data->reduce(bkgStr));

    RooKeysPdf bkg("bkg","bkg",*cosTheta,*data2);

    sprintf(namestr,"%d_%.2f_%.2f",i,va,vb);

    Double_t g1, g2;
    g1 = 0.0;
    g2 = -0.5;

    RooRealVar G000("G000", "G000",  0.5);
    RooRealVar G001("G001", "G001",  g1, -5., 5.);
    RooRealVar G002("G002", "G002",  g2, -5., 5.);
    RooRealVar G003("G003", "G003", 0.0);//, -1., 1.);
    RooRealVar G004("G004", "G004", 0.0);//, -1., 1.);
    RooRealVar    a("a",    "a",    va);
    RooRealVar    b("b",    "b",    vb);
    RooRealVar    n("n",    "n",    1.);

    RooRealVar nsig("nsig", "nsig", 996.); //1000., -100., 5000.);
    RooRealVar nbkg("nbkg", "nbkg", 217.); // 200., -100., 5000.);

    RooB2Kll sig("sig", "sig", *cosTheta, G000, G001, G002, G003, G004, a, b, n);

    TFile * fileK1 = TFile::Open("../efficiencies/Dveto_200.root");
    TH1D  * histK1  = dynamic_cast<TH1D*>(fileK1->Get("efficiency_17"));
    RooDataHist * dataK1 = new RooDataHist("dataK1", "", RooArgSet(*cosTheta), histK1);
    RooHistPdf DVeto("DVeto","", RooArgSet(*cosTheta), *dataK1, 2);

    TFile * fileK2 = TFile::Open("../efficiencies/Psiveto_200.root");
    TH1D  * histK2  = dynamic_cast<TH1D*>(fileK2->Get("efficiency_17"));
    RooDataHist * dataK2 = new RooDataHist("dataK2", "", RooArgSet(*cosTheta), histK2);
    RooHistPdf PsiVeto("PsiVeto","", RooArgSet(*cosTheta), *dataK2, 2);

    RooProdPdf sigeff("sigeff","sigeff", RooArgList(sig,DVeto,PsiVeto));

    RooAddPdf pdf("pdf","pdf",RooArgList(sig,bkg), RooArgList(nsig,nbkg));
    RooAddPdf pdfeff("pdfeff","pdfeff",RooArgList(sigeff,bkg), RooArgList(nsig,nbkg));

    RooFitResult * fitresult = pdfeff.fitTo(*data1,Save(kTRUE), Minos(kFALSE), NumCPU(4), SumW2Error(kTRUE));
    
    cout << i << "\t" << a.getVal() << "\t" << b.getVal() << "\t" << fitresult->minNll() << "\t" << fitresult->status() << "\t" 
         << G000.getVal() << "\t" << G001.getVal() << "\t" << G002.getVal() << "\t" << G003.getVal() << "\t" << G004.getVal() << endl;
    
    RooPlot * frame = cosTheta->frame(50);
    data1->plotOn(frame);
    pdfeff.plotOn(frame);
    pdfeff.plotOn(frame, RooFit::Components("sigeff"),RooFit::LineStyle(kDashed));
    pdfeff.plotOn(frame, RooFit::Components("bkg"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
   
    //std::cout << sig.getNorm(RooArgSet(*cosTheta)) << "\t" << nsig.getVal() << std::endl;
    //std::cout << bkg.getNorm(RooArgSet(*cosTheta)) << "\t" << nbkg.getVal() << std::endl;
    
    //std::cout << sig.createIntegral(RooArgSet(*cosTheta),RooFit::NormSet(*cosTheta))->getVal() << "\t" << nsig.getVal() << std::endl;
    //std::cout << bkg.createIntegral(RooArgSet(*cosTheta),RooFit::NormSet(*cosTheta))->getVal() << "\t" << nbkg.getVal() << std::endl;

    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TCanvas c("fit","fit", 800, 800);
    frame->Draw();
    c.SaveAs("bkgFits/fit"+TString(namestr)+".png");
    c.SaveAs("bkgFits/fit"+TString(namestr)+".pdf");
}

int main(int argc, char** argv)
{
    Double_t a(-0.43), b(0.);
    Int_t i(0);

    if(argc>1) {
        i = atoi(argv[1]);
        if(argc>3) {
            a = atof(argv[2]);
            b = atof(argv[3]);
	}
    }
    
    fit(i,a,b);

    return 0;
}

