#include "outputFunctions.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"

#include "RooAbsData.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"

TString gSaveDir=".";

// function to make plot of 1D projection of the fit
void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf,
             std::vector<std::string>& sig_pdfs,
	     std::vector<std::string>& bkg_pdfs,
	     TString name, TString title, TString label,
	     int typeIdx,
	     RooCmdArg /*extraArg0*/,
	     RooCmdArg /*extraArg1*/,
	     RooCmdArg /*extraArg2*/) {
	if(gSaveDir=="") {
		std::cout << "WARNING in plotFit: gSaveDir not set" << std::endl;
		std::cout << "                    setting to \".\"" << std::endl;
		gSaveDir=".";
	}
	TCanvas c1;
	c1.SetBottomMargin(0.19);
	RooPlot* plot = var.frame(min, max, nbins);
	int iCol(0);
	int colours[6] = {kBlue, kRed, kGreen+2, kMagenta, kOrange, kGray};
	int fills[6] = {3245, 3454, 3644, 3205, 3495, 3690};

	RooCmdArg normCmd = RooFit::NormRange("FIT");
	//RooCmdArg normCmd = RooFit::Normalization(1.);

	//if we have a simultaneous PDF then use dataset to project the index category
	RooCmdArg projCmd   = RooCmdArg::none();
	RooCmdArg cutCmd    = RooCmdArg::none();
	RooCmdArg sliceCmd  = RooCmdArg::none();
	if(pdf.InheritsFrom(RooSimultaneous::Class())) {
		RooSimultaneous* sim = dynamic_cast<RooSimultaneous*>(&pdf);
		const RooCategory* cat = dynamic_cast<const RooCategory*>(&sim->indexCat());
		projCmd = RooFit::ProjWData(*cat, *dh);
		//make plots of the slices
		if(typeIdx>-1) {// && typeIdx<cat->numTypes()) {
			if(!cat->lookupType(typeIdx)) {
				delete plot;
				return;
			}
			TString cutStr=cat->GetName();
			cutStr+="==";
			cutStr+=cat->GetName();
			cutStr+="::";
			cutStr+=cat->lookupType(typeIdx)->GetName();
			//skip if no entries in this slice
			if(dh->sumEntries(cutStr)==0) {
				delete plot;
				return;
			}
			if(!sim->getPdf(cat->lookupType(typeIdx)->GetName())) {
				delete plot;
				return;
			}
			cutCmd = RooFit::Cut(cutStr);
			sliceCmd = RooFit::Slice(*const_cast<RooCategory*>(cat), cat->lookupType(typeIdx)->GetName());
		//moved to the end
		//} else {
		//	for(int iType=0; iType<cat->numTypes(); ++iType) {
		//		//TODO//std::cout << iType << ": making plot for " << cat->lookupType(iType)->GetName() << std::endl;
		//		plotFit(var,min,max,nbins,dh,pdf,sig_pdfs,bkg_pdfs,name+cat->lookupType(iType)->GetName(),title,label,iType);
		//	}
		}
	}

	//plot data
	dh->plotOn( plot, cutCmd );
	//pdf.plotOn( plot, RooFit::Components("*"), normCmd, sliceCmd, projCmd, RooFit::LineColor(kBlack));
	
	//TODO skip background components for now - one fit freezes at plotOn
//	if(typeIdx<0) {
		//plot selected components
		for (std::vector<std::string>::iterator it = bkg_pdfs.begin(); it != bkg_pdfs.end(); ++it){
			if(iCol>=6) iCol=0;
			pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, sliceCmd, projCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(9));//, extraArg0, extraArg1, extraArg2 );//TODO
			++iCol;
		}

		for (std::vector<std::string>::iterator it = sig_pdfs.begin(); it != sig_pdfs.end(); ++it){
			if(iCol>=6) iCol=0;
			pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, sliceCmd, projCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(kDashed));
			pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, sliceCmd, projCmd, RooFit::LineColor( colours[iCol] ), RooFit::VLines(), RooFit::FillStyle(fills[iCol]), RooFit::FillColor( colours[iCol] ), RooFit::DrawOption("LF"));
			++iCol;
		}
	//}

	//plot full PDF
	pdf.plotOn( plot, RooFit::Components("*"), normCmd, sliceCmd, projCmd, RooFit::LineColor(kBlack));

	//re-plot data
	dh->plotOn( plot, cutCmd );

	plot->SetXTitle(title);
	plot->SetTitle("");
	plot->GetXaxis()->SetLabelOffset(0.02);
	plot->GetXaxis()->SetTitleOffset(1.18);
	plot->Draw();
	TPaveText txt(0.65,0.8,0.9,0.9,"BRNDC");
	txt.AddText(label);
	txt.SetFillColor(0);
	txt.SetTextAlign(12);
	txt.SetBorderSize(0);
	if(label!="") txt.Draw();
	//std::cout << "!" << label << std::endl;//TODO

	c1.SaveAs(gSaveDir+"/"+name+"_fit.pdf");
	//TODO//c1.SaveAs(gSaveDir+"/"+name+"_fit.png");
	//TODO//c1.SetLogy();
	//TODO//c1.SaveAs(gSaveDir+"/"+name+"_fit_log.pdf");

	delete plot;

	//make subplots
	if(pdf.InheritsFrom(RooSimultaneous::Class())) {
		if(typeIdx==-1) {
			const RooCategory* cat = dynamic_cast<const RooCategory*>(&dynamic_cast<RooSimultaneous*>(&pdf)->indexCat());
//			for(int iType=0; iType<cat->numTypes(); ++iType) {
//			for(auto catType=cat->begin(); catType!=cat->end(); ++catType) {
			TIterator* it = cat->typeIterator();
			RooCatType* catType(0);
			while(catType = static_cast<RooCatType*>(it->Next())) {
				plotFit(var,min,max,nbins,dh,pdf,sig_pdfs,bkg_pdfs,name+catType->GetName(),title,label,catType->getVal());
			}
		}
	}
}

void plotComparison(TString plotName, TString yLabel, std::vector<TH1D*> hists) {
	int colours[7] = {kBlue, kRed, kGreen+2, kMagenta, kOrange, kGray, kTeal};
	std::vector<TH1D*> plotHists;
	std::vector<TLine> plotDivs;
	std::vector<TLine> plotDivsRatio;

	if(hists.empty()) return;

	size_t nHists = hists.size();
	size_t nBins = hists.front()->GetNbinsX();
	TH1D basePlot("basePlot","",nBins,0.,nBins);

	for(auto it=hists.begin(); it!=hists.end(); ++it) {
		if((*it)->GetNbinsX() != nBins) {
			return;
		}
	}

	double min(0.), max(0.);
	double plotmin(0.), plotmax(0.);

	for(size_t ihist=0; ihist<nHists; ++ihist) {
		TString name = "plot";
		TH1D* hist = hists.at(ihist);
		name += hist->GetName();
		plotHists.push_back(new TH1D(name,"",nBins*nHists,0.,nBins));
		for(size_t ibin=0; ibin<nBins; ++ibin) {
			plotHists.back()->SetBinContent(ibin*nHists + ihist+1, hist->GetBinContent(ibin+1));
			plotHists.back()->SetBinError(ibin*nHists + ihist+1, hist->GetBinError(ibin+1));
			if(hist->GetBinContent(ibin+1)+hist->GetBinError(ibin+1) > max) max = hist->GetBinContent(ibin+1)+hist->GetBinError(ibin+1);
			if(hist->GetBinContent(ibin+1)-hist->GetBinError(ibin+1) < min) min = hist->GetBinContent(ibin+1)-hist->GetBinError(ibin+1);
		}
		plotHists.back()->SetMarkerColor(colours[ihist%7]);
		plotHists.back()->SetLineColor(colours[ihist%7]);
		plotHists.back()->SetMarkerStyle(20+ihist);
	}

	if(min!=0) plotmin = min - 0.05*(max-min);
	else plotmin = min;
	plotmax= max + 0.05*(max-min);

	for(size_t iline=1; iline<nBins; ++iline) {
		plotDivs.push_back(TLine(iline,plotmin,iline,plotmax));
		plotDivs.back().SetLineColor(kBlack);
		plotDivs.back().SetLineStyle(kDashed);
	}
	for(size_t ibin=1; ibin<=nBins; ++ibin) {
		TString binName = TString::Format("%3.0f-%3.0f GeV",hists.front()->GetBinLowEdge(ibin)/1e3,hists.front()->GetBinLowEdge(ibin+1)/1e3);
		basePlot.GetXaxis()->SetBinLabel(ibin,binName);
	}

	basePlot.SetMinimum(plotmin);
	basePlot.SetMaximum(plotmax);
	basePlot.GetYaxis()->SetTitle(yLabel);

	TCanvas c;
	basePlot.Draw("EP");
	for(auto it=plotHists.begin(); it!=plotHists.end(); ++it) {
		(*it)->Draw("EP SAME");
	}
	for(auto it=plotDivs.begin(); it!=plotDivs.end(); ++it) {
		(*it).Draw();
	}
	c.SaveAs(gSaveDir+"/"+plotName+".pdf");

	//now make the ratio plot
	min=0.9;
	max=1.1;
	for(size_t ihist=0; ihist<plotHists.size()-1; ihist+=2) {
		plotHists[ihist]->Rebin();
		plotHists[ihist+1]->Rebin();
		plotHists[ihist+1]->Divide(plotHists[ihist]);
		for(size_t ibin=0; ibin<plotHists[ihist+1]->GetNbinsX(); ++ibin) {
			if(plotHists[ihist+1]->GetBinContent(ibin+1)==0) continue;
			if(plotHists[ihist+1]->GetBinContent(ibin+1)+plotHists[ihist+1]->GetBinError(ibin+1) > max) max = plotHists[ihist+1]->GetBinContent(ibin+1)+plotHists[ihist+1]->GetBinError(ibin+1);
			if(plotHists[ihist+1]->GetBinContent(ibin+1)-plotHists[ihist+1]->GetBinError(ibin+1) < min) min = plotHists[ihist+1]->GetBinContent(ibin+1)-plotHists[ihist+1]->GetBinError(ibin+1);
		}
	}

	if(min!=0) plotmin = min - 0.05*(max-min);
	else plotmin = min;
	plotmax= max + 0.05*(max-min);

	for(size_t iline=1; iline<nBins; ++iline) {
		plotDivsRatio.push_back(TLine(iline,plotmin,iline,plotmax));
		plotDivsRatio.back().SetLineColor(kBlack);
		plotDivsRatio.back().SetLineStyle(kDashed);
	}

	basePlot.SetMinimum(plotmin);
	basePlot.SetMaximum(plotmax);
	basePlot.GetYaxis()->SetTitle("ratio "+yLabel);

	basePlot.Draw("EP");
	for(size_t ihist=0; ihist<plotHists.size()-1; ihist+=2) {
		plotHists[ihist+1]->Draw("EP SAME");
	}
	for(auto it=plotDivsRatio.begin(); it!=plotDivsRatio.end(); ++it) {
		(*it).Draw();
	}

	c.SaveAs(gSaveDir+"/"+plotName+"-ratio.pdf");

	//cleanup
	for(auto it=plotHists.begin(); it!=plotHists.end(); ++it) {
		delete (*it);
	}
}

//function to set precision and exponent for printing a parameter
void getPrecision(double value, double error, int& precision, int& exponent) {
	while(TMath::Abs(value) < TMath::Power(10,exponent)) exponent-=1;
	while(TMath::Abs(value) > TMath::Power(10,exponent+1)) exponent+=1;

	if(error==0) {
		precision = 6 - exponent;
		return;
	}

	int exponentErr(0);
	while(TMath::Abs(error) < TMath::Power(10,exponentErr)) exponentErr-=1;
	while(TMath::Abs(error) > TMath::Power(10,exponentErr+1)) exponentErr+=1;

	if(error < 3.5*TMath::Power(10,exponentErr)) precision = 1-exponentErr;
	else precision = -exponentErr;

	if(precision > 6 - exponent) precision = 6 - exponent;
}

// function to print parameters to a file
void printParams(TString file, const RooArgList& params) {
	if(gSaveDir=="") {
		std::cout << "WARNING in plotFit: gSaveDir not set" << std::endl;
		std::cout << "                    setting to \".\"" << std::endl;
		gSaveDir=".";
	}
	FILE * pFile = fopen((gSaveDir+"/"+file).Data(), "w");

	int nPar = params.getSize();

	for ( int i=0; i<nPar; ++i) {
		RooRealVar* par = dynamic_cast<RooRealVar*>(params.at(i));

		TString title = par->getTitle();
		float value = par->getValV();
		float error = par->getError();
		if(error==0) continue;

		int exponent(0);
		int exponentErr(0);
		int precision(0);

		while(TMath::Abs(value) < TMath::Power(10,exponent)) exponent-=1;
		while(TMath::Abs(value) > TMath::Power(10,exponent+1)) exponent+=1;
		while(TMath::Abs(error) < TMath::Power(10,exponentErr)) exponentErr-=1;
		while(TMath::Abs(error) > TMath::Power(10,exponentErr+1)) exponentErr+=1;

		if(error < 3.5*TMath::Power(10,exponentErr)) precision = 1-exponentErr;
		else precision = -exponentErr;

		if(exponent<-2) {
			title+=" ($10^{";
			title+=exponent;
			title+="}$)";

			precision+=exponent;

			value/=TMath::Power(10,exponent);
			error/=TMath::Power(10,exponent);
		}

		if(precision<0) precision=0;

		TString format = "%-30s";
		format+="\t& $% 8.";
		format+=precision;
		format+="f \\pm % 8.";
		format+=precision;
		format+="f$ \\\\\n";

		fprintf(pFile, format, title.Data(), value, error);
	}
	fclose(pFile);
}

// function to print parameters to a file with split parameters on a single line
void printParamsSim(TString file, const RooArgList& params) {
	int nCol(0);
	int maxTitleLen(40);//start at 40 but expand if we find any long titles
	//store names and titles of rows
	std::vector<TString> rowLookup;
	std::vector<TString> titles;
	//store values and errors of all parameters
	std::map<int, std::map<int, double> > values;
	std::map<int, std::map<int, double> > errors;
	//store the precision and exponent for printing each row
	std::vector<int> precisions;
	std::vector<int> exponents;

	int nPar = params.getSize();

	for ( int i=0; i<nPar; ++i) {
		RooRealVar* par = dynamic_cast<RooRealVar*>(params.at(i));
		TString name = par->GetName();
		TString baseName;
		TString binName;
		int col(-1);

		if(name.Last('_') != -1) {
			baseName = name(0,name.Last('_'));
			binName  = name(name.Last('_')+1,name.Length());
			col = TString(binName(TRegexp("[0-9]+"))).Atoi();
			if(col>=nCol) nCol = col+1;
		} else {
			baseName = name;
		}

		uint row;
		for(row=0; row<rowLookup.size(); ++row) {
			if(rowLookup.at(row) == baseName) {
				break;
			}
		}
		if(row==rowLookup.size()) {
			rowLookup.push_back(baseName);
			TString title = par->getTitle();
			TString baseTitle;
			if(title.Last('(')!=-1) baseTitle = TString(title(0,title.Last('(')));
			else baseTitle = title;
			titles.push_back(baseTitle);
			if(baseTitle.Length() > maxTitleLen) maxTitleLen=baseTitle.Length();
		}

		values[row][col] = par->getValV();
		errors[row][col] = par->getError();

		int prec(0), exp(0);
		getPrecision(par->getValV(),par->getError(),prec,exp);
		if(row < precisions.size()) {
			//store highest precision but largest exponent
			if(prec > precisions[row]) precisions[row] = prec;
			if(exp > exponents[row]) exponents[row] = exp;

			if(precisions[row] > 6 - exponents[row]) precisions[row] = 6 - exponents[row];
		} else {
			precisions.push_back(prec);
			exponents.push_back(exp);
		}
	}

	//allow for an exponent
	maxTitleLen += 10;
	TString formatT = TString::Format("%%-%ds",maxTitleLen);

	if(gSaveDir=="") {
		std::cout << "WARNING in plotFit: gSaveDir not set" << std::endl;
		std::cout << "                    setting to \".\"" << std::endl;
		gSaveDir=".";
	}

	FILE * pFile = fopen((gSaveDir+"/"+file).Data(), "w");

	for(uint irow=0; irow<titles.size(); ++irow) {
		if(values.find(irow) != values.end()) {
			double scale = 1.;
			if(exponents[irow]<-2) {
				titles[irow]+=TString::Format(" ($10^{%d}$)",exponents[irow]);
				precisions[irow]+=exponents[irow];
				scale = TMath::Power(0.1,exponents[irow]);
			}
			if(precisions[irow]<0) precisions[irow]=0;

			TString format = TString::Format("%% 8.%df",precisions[irow]);

			fprintf(pFile,formatT,titles[irow].Data());
			if(values[irow].find(-1) != values[irow].end()) {
				fprintf(pFile,"\t& \\multicolumn{%d}{c}{$"+format+" \\pm "+format+"$}\\\\\n",nCol,scale*values[irow][-1],scale*errors[irow][-1]);
			} else {
				for(int icol=0; icol<nCol; ++icol) {
					if(values[irow].find(icol)!=values[irow].end()) {
						fprintf(pFile,"\t& $"+format+" \\pm "+format+"$",scale*values[irow][icol],scale*errors[irow][icol]);
					} else {
						fprintf(pFile,"\t&           --           ");
					}
				}
				fprintf(pFile,"\\\\\n");
			}
		} else {
			fprintf(pFile,formatT+"\t& \\multicolumn{%d}{c}{--}\\\\\n",titles[irow].Data(),nCol);
		}
	}

	fclose(pFile);
}

//alternative form that returns all non-constant parameters of the pdf
//dataset is used to remove abscissa
void printParams(TString file, RooAbsData const* data, RooAbsPdf const* pdf) {
	RooArgList* floatPars = static_cast<RooArgList*>(pdf->getParameters(*data)->selectByAttrib("Constant",kFALSE));
	if(pdf->InheritsFrom(RooSimultaneous::Class())) {
		printParamsSim(file, *floatPars);
	} else {
		printParams(file, *floatPars);
	}
}

void printAllParams(TString file, RooAbsData const* data, RooAbsPdf const* pdf) {
	RooArgList* floatPars = static_cast<RooArgList*>(pdf->getParameters(*data)->selectByName("*"));
	if(pdf->InheritsFrom(RooSimultaneous::Class())) {
		printParamsSim(file, *floatPars);
	} else {
		printParams(file, *floatPars);
	}
}

void setLHCbStyle() {
  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
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
  lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
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
  //lhcbStyle->SetOptStat(0);  
  lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
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
/*
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
*/
  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  //std::cout << "-------------------------" << std::endl;  
  //std::cout << "Set LHCb Style - Feb 2012" << std::endl;
  //std::cout << "-------------------------" << std::endl;  
  
}


