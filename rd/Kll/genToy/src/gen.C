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


void gen(Int_t iToy, Int_t qSqBin){
    TString namestr("toy/toy_"); namestr+=iToy; namestr+="_"; namestr+=qSqBin; namestr+=".dat";

    Double_t bkgMin(5400), bkgMax(5970);

    TString bkgStr="Bplus_M > "; bkgStr+=bkgMin; bkgStr+=" && Bplus_M < "; bkgStr+=bkgMax;

    TFile * file(0);
    TTree * DecayTree(0);
    RooRealVar * cosTheta(0);
    RooRealVar * Bplus_M(0);
    RooDataSet * data(0);
   
    TString fname("~/lhcb/rd/Kll/tuples/fromPatrick/Kmm_Q"); fname+=qSqBin; fname+="_reduced.root";
    file = TFile::Open(fname);
    DecayTree = dynamic_cast<TTree*>(file->Get("finalTree_KMuMu"));
    cosTheta = new RooRealVar("cosThetaL", "cosThetaL", -1., 1.);
    Bplus_M  = new RooRealVar("Bplus_M", "Bplus_M", 5170., 5970.);
    data = new RooDataSet("data", "", DecayTree, RooArgList(*cosTheta,*Bplus_M));

    RooDataSet* data2 = dynamic_cast<RooDataSet*>(data->reduce(bkgStr));
    RooKeysPdf bkg("bkg","bkg",*cosTheta,*data2);

    Double_t g0, g1, g2, g3, g4, g5, g6, va, vb, ignore;
    Int_t nsig, nbkg;

    std::ifstream fin;

    //load efficiency parameters
    fin.open("../efficiencies/effParams_PIDonly_fixB_50.txt");
    Int_t row(0);
    fin >> row;
    while(row!=qSqBin) {
	    fin.ignore(256,'\n');
	    fin >> row;
    }
    fin >> ignore >> ignore >> va >> ignore >> vb >> ignore;
    std::cout << va << "\t" << vb << std::endl;

    fin.close();

    //load moments parameters
    fin.open("../MoM/moments.dat");
    fin >> row;
    while(row!=qSqBin) {
	    fin.ignore(256,'\n');
	    fin >> row;
    }
    fin >> g0 >> g1 >> g2 >> g3 >> g4 >> g5 >> g6;
    std::cout << g0 << "\t" << g1 << "\t" << g2 << "\t" << g3 << "\t" << g4 << "\t" << g5 << "\t" << g6 << std::endl;

    fin.close();

    //load sig/bkg yields
    TString yieldsname("../massFit/yields/"); yieldsname+= qSqBin; yieldsname+= ".dat";
    fin.open(yieldsname);
    fin >> nsig >> ignore >> nbkg;
    std::cout << nsig << "\t" << nbkg << std::endl;

    fin.close();

    RooRealVar G000("G000", "G000",  g0);
    RooRealVar G001("G001", "G001",  g1);
    RooRealVar G002("G002", "G002",  g2);
    RooRealVar G003("G003", "G003",  g3);
    RooRealVar G004("G004", "G004",  g4);
    RooRealVar G005("G005", "G005",  g5);
    RooRealVar G006("G006", "G006",  g6);
    RooRealVar    a("a",    "a",     va);
    RooRealVar    b("b",    "b",     vb);
    RooRealVar    n("n",    "n",    1.0);

    RooB2Kll sig("sig", "sig", *cosTheta, G000, G001, G002, G003, G004, G005, G006, a, b, n);

    TString vetoHistName("efficiencyHist_");
    vetoHistName+=qSqBin;
    TFile * fileK1 = TFile::Open("../efficiencies/Dveto_200.root");
    TH1D  * histK1  = dynamic_cast<TH1D*>(fileK1->Get(vetoHistName));
    RooDataHist * dataK1 = new RooDataHist("dataK1", "", RooArgSet(*cosTheta), histK1);
    RooHistPdf DVeto("DVeto","", RooArgSet(*cosTheta), *dataK1, 2);

    TFile * fileK2 = TFile::Open("../efficiencies/Psiveto_200.root");
    TH1D  * histK2  = dynamic_cast<TH1D*>(fileK2->Get(vetoHistName));
    RooDataHist * dataK2 = new RooDataHist("dataK2", "", RooArgSet(*cosTheta), histK2);
    RooHistPdf PsiVeto("PsiVeto","", RooArgSet(*cosTheta), *dataK2, 2);

    RooProdPdf sigeff("sigeff","sigeff", RooArgList(sig,DVeto,PsiVeto));

//    RooAddPdf pdf("pdf","pdf",RooArgList(sig,bkg), RooArgList(nsig,nbkg));
//    RooAddPdf pdfeff("pdfeff","pdfeff",RooArgList(sigeff,bkg), RooArgList(nsig,nbkg));

    RooRandom::randomGenerator()->SetSeed(iToy);

    RooDataSet * genData  = sigeff.generate(RooArgSet(*cosTheta), nsig);
//    RooDataSet * genData2 = bkg.generate(RooArgSet(*cosTheta), nsig);
//    genData->append(*genData2);
    genData->write(namestr);

    file->Close();
    fileK1->Close();
    fileK2->Close();

    delete genData;
//    delete genData2;

//    RooFitResult * fitresult = pdfeff.fitTo(*data1,Save(kTRUE), Minos(kFALSE), NumCPU(4), SumW2Error(kTRUE));
//    
//    cout << i << "\t" << a.getVal() << "\t" << b.getVal() << "\t" << fitresult->minNll() << "\t" << fitresult->status() << "\t" 
//         << G000.getVal() << "\t" << G001.getVal() << "\t" << G002.getVal() << "\t" << G003.getVal() << "\t" << G004.getVal() << endl;
//    
//    RooPlot * frame = cosTheta->frame(50);
//    data1->plotOn(frame);
//    pdfeff.plotOn(frame);
//    pdfeff.plotOn(frame, RooFit::Components("sigeff"),RooFit::LineStyle(kDashed));
//    pdfeff.plotOn(frame, RooFit::Components("bkg"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
//   
//    //std::cout << sig.getNorm(RooArgSet(*cosTheta)) << "\t" << nsig.getVal() << std::endl;
//    //std::cout << bkg.getNorm(RooArgSet(*cosTheta)) << "\t" << nbkg.getVal() << std::endl;
//    
//    //std::cout << sig.createIntegral(RooArgSet(*cosTheta),RooFit::NormSet(*cosTheta))->getVal() << "\t" << nsig.getVal() << std::endl;
//    //std::cout << bkg.createIntegral(RooArgSet(*cosTheta),RooFit::NormSet(*cosTheta))->getVal() << "\t" << nbkg.getVal() << std::endl;
//
//    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
//    TCanvas c("fit","fit", 800, 800);
//    frame->Draw();
//    c.SaveAs("bkgFits/fit"+TString(namestr)+".png");
//    c.SaveAs("bkgFits/fit"+TString(namestr)+".pdf");
}

int main(int argc, char** argv)
{
    Int_t i(0), q(0);
    if(argc>1) {
        i = atoi(argv[1]);
        if(argc>2) {
            q = atof(argv[2]);
	}
    }
    
    gen(i,q);

    return 0;
}

