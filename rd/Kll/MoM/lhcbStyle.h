#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>

#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

TStyle* lhcbStyle; // general lhcb style
TPaveText* lhcbName; // standard lhcb text for plot
TText* lhcbLabel; // style for Ttext
TLatex *lhcbLatex; //style for TLatex;

void setLHCbStyle() {

  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.07); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);
  
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  // add LHCb label
  lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);

  lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;  
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;
  std::cout << "-------------------------" << std::endl;  
  
}

void printLHCb(TString optLR="L", TString optPrelim="Final", TString optText="")
{
//////////////////////////////////////////////////////////////////////////
// routine to print 'LHCb', 'LHCb Preliminary' on plots 
// options: optLR=L (top left) / R (top right) of plots
//          optPrelim= Final (LHCb), Prelim (LHCb Preliminary), Other
//          optText= text printed if 'Other' specified
////////////////////////////////////////////////////////////////////
  if (optLR=="R"){    
    lhcbName = new TPaveText(0.65 - lhcbStyle->GetPadRightMargin(),
					   0.75 - lhcbStyle->GetPadTopMargin(),
					   0.85 - lhcbStyle->GetPadRightMargin(),
					   0.85 - lhcbStyle->GetPadTopMargin(),
					   "BRNDC");
  }
  else if (optLR=="FTR"){
    lhcbName = new TPaveText(0.88 - lhcbStyle->GetPadRightMargin(),
					   0.87 - lhcbStyle->GetPadTopMargin(),
					   0.98 - lhcbStyle->GetPadRightMargin(),
					   0.98 - lhcbStyle->GetPadTopMargin(),
                                        "BRNDC");
  }
  else if (optLR=="R2"){
    lhcbName = new TPaveText(0.60 - lhcbStyle->GetPadRightMargin(),
					   0.88 - lhcbStyle->GetPadTopMargin(),
					   0.85 - lhcbStyle->GetPadRightMargin(),
					   0.96 - lhcbStyle->GetPadTopMargin(),
                                        "BRNDC");
  }
  else if (optLR=="RW"){
    lhcbName = new TPaveText(0.40 - lhcbStyle->GetPadRightMargin(),
					   0.78 - lhcbStyle->GetPadTopMargin(),
					   0.90 - lhcbStyle->GetPadRightMargin(),
					   0.85 - lhcbStyle->GetPadTopMargin(),
					   "BRNDC");
  }
  else if (optLR=="TR"){
    lhcbName = new TPaveText(0.55 - lhcbStyle->GetPadRightMargin(),
					   0.90, //0.91,
					   1.00 - lhcbStyle->GetPadRightMargin(),
					   1.00, //0.97,
					   "BRNDC");
  }
  else if (optLR=="L"){
    lhcbName = new TPaveText(lhcbStyle->GetPadLeftMargin() + 0.05,
				       0.88 - lhcbStyle->GetPadTopMargin(),//87
				       lhcbStyle->GetPadLeftMargin() + 0.30,
				       0.96 - lhcbStyle->GetPadTopMargin(),//95
                                        "BRNDC");
  }
  else if (optLR=="LW"){
    lhcbName = new TPaveText(lhcbStyle->GetPadLeftMargin() + 0.05,
				       0.78 - lhcbStyle->GetPadTopMargin(),
				       lhcbStyle->GetPadLeftMargin() + 0.55,
				       0.85 - lhcbStyle->GetPadTopMargin(),
                                        "BRNDC");
  }
  else{
   std::cout << "printLHCb: option unknown" << optLR << std::endl;   
  }
  if (optPrelim=="Final"){
    lhcbName->AddText("LHCb");
  }
  else if (optPrelim=="MC"){
    lhcbName->AddText("#splitline{LHCb}{#scale[1.0]{Monte Carlo}}");  
  }
  else if (optPrelim=="S"){
    lhcbName->AddText("#splitline{LHCb}{#scale[1.0]{Simulation}}");  
  }
  else if (optPrelim=="Prelim"){
    lhcbName->AddText("#splitline{LHCb}{#scale[1.0]{Preliminary}}");  
  }
  else if (optPrelim=="Other"){
    lhcbName->AddText(optText);
  }
  else{
    std::cout << "printLHCb: option unknown " << optPrelim << std::endl;   
  }

  lhcbName->SetFillColor(0);
  lhcbName->SetFillStyle(0); //26);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
  lhcbName->SetLineColor(0);
  lhcbName->SetLineWidth(0);
  lhcbName->Draw();
}

