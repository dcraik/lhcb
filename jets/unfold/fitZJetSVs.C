#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <boost/progress.hpp>

#include "outputFunctions.h"
#include "SVFitter.h"
#include "ZFitter.h"
#include "MCJets.h"

//globals to save passing these around
TString savedir("output-fitZj-statOnly-24bins-newFit-new-new-fixSVSel-fixJetPVMatch2-ptCorr");

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
//TString charmSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_14X.root";
//TString beautySimFile = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_15X.root";
//TString lightSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_101.root";
//TString dataFile      = "/tmp/dcraik/skimmed.root";
TString charmSimFile  = "/data/zjet/zjet_sim4.root";
TString beautySimFile = "/data/zjet/zjet_sim5.root";
TString lightSimFile  = "/data/zjet/zjet_sim0.root";
TString lightDataFile  = "/data/zjet/zjet_back_201X.root";
TString dataFile      = "/data/zjet/zjet_201X.root";
TString ssDataFile    = "/data/zjet/zjet_ss_201X.root";

//TString lightHistFile = "svFitHists0_zj.root";
//TString charmHistFile = "svFitHists4_zj.root";
//TString beautyHistFile = "svFitHists5_zj.root";
//TString dataHistFile = "svFitHistsD_zj.root";

void fillJetsHist(TH2D* h) {
	TFile* f = TFile::Open(dataFile);
	if(!f) return;
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return;

	double JetPT, ZPZ, ZE;

	int NPV;//TODO
	t->SetBranchAddress("NPV", &NPV);//TODO
	t->SetBranchAddress("JetPT", &JetPT);
	t->SetBranchAddress("ZPZ",   &ZPZ);
	t->SetBranchAddress("ZE",    &ZE);

	std::cout << t->GetEntries() << std::endl;
	for(int i=0; i<t->GetEntries(); ++i) {
		t->GetEntry(i);
		//if(NPV==1) continue;//TODO
		//std::cout << JetPT << "\t" << 0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)) << "\t" << h->FindBin(JetPT,0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ))) << std::endl;
		h->Fill(JetPT,0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)));
	}
	//TString cutStr="0.5*TMath::Log((ZE+ZPZ)/(ZE-ZPZ)):JetPT>>";
	//cutStr+=h->GetName();
	//t->Draw(cutStr,"","goff");
}

int main() {
	gSaveDir = savedir;
	gSystem->Exec("mkdir -p "+savedir);
	gStyle->SetOptStat(0);

	double nB(0.), eB(0.), nC(0.), eC(0.), nQ(0.), eQ(0.), nZ(0.), eZ(0.);

	const int npt(3);
	double ptBounds[npt+1] = {15000.,20000.,30000.,50000.};//TODO 100000.
	//const int nzy(2);
	//double zyBounds[nzy+1] = {2.0,3.0,4.5};
	const int nzy(5);
	double zyBounds[nzy+1] = {2.0,2.5,3.0,3.5,4.0,4.5};

	TH1D ptScheme("ptScheme","",npt,ptBounds);
	TH1D zyScheme("zyScheme","",nzy,zyBounds);

	TH2D hC("hC","",npt,ptBounds,nzy,zyBounds);
	TH2D hJ("hJ","",npt,ptBounds,nzy,zyBounds);
	TH2D hE("hE","",npt,ptBounds,nzy,zyBounds);
	TH2D hE2("hE2","",npt,ptBounds,1,zyBounds[0],zyBounds[nzy]);
	TH2D hR("hR","",npt,ptBounds,nzy,zyBounds);
	TH2D hUR("hUR","",npt,ptBounds,nzy,zyBounds);
	TH2D hCE("hCE","",npt,ptBounds,nzy,zyBounds);
	TH2D hUCE("hUCE","",npt,ptBounds,nzy,zyBounds);
	TH2D hZ("hZ","",npt,ptBounds,nzy,zyBounds);
	TH2D hCSS("hCSS","",npt,ptBounds,nzy,zyBounds);

	TH1D hR20up("hR20up","",nzy,zyBounds);
	TH1D hUR20up("hUR20up","",nzy,zyBounds);

	hC.Sumw2();
	hJ.Sumw2();
	hE.Sumw2();
	hE2.Sumw2();
	hR.Sumw2();
	hUR.Sumw2();
	hCE.Sumw2();
	hUCE.Sumw2();
	hZ.Sumw2();
	hCSS.Sumw2();

	hR20up.Sumw2();
	hUR20up.Sumw2();

	//unsigned int nmcpt=4;
	//double* mcPtInputWeightBins  = new double[nmcpt +1]{10000.,15000.,20000.,50000.,100000.};

	//TH1D* jetTruePtWeights = new TH1D("jetTruePtWeights","",nmcpt,mcPtInputWeightBins);
	//jetTruePtWeights->SetBinContent(1,1.);
	//jetTruePtWeights->SetBinContent(2,0.6);
	//jetTruePtWeights->SetBinContent(3,0.5);
	//jetTruePtWeights->SetBinContent(4,0.04);

	MCJets wmc("Zj");
	wmc.setInputs(lightSimFile,charmSimFile,beautySimFile);
	//wmc.setInputTruePtWeights(MCJets::jetRecoSV4,jetTruePtWeights);
	//wmc.setInputTruePtWeights(MCJets::jetRecoSV5,jetTruePtWeights);
	if(!wmc.weightMC(MCJets::jetRecoSV4,&hE2)) return 1;
	if(!wmc.weightMC(MCJets::jetRecoSV5)) return 1;

	//for(int j=1; j<=nzy; ++j) {
	//	hE.SetBinContent(1,j,0.191);
	//	hE.SetBinError(1,j,0.);//0.021);
	//	hE.SetBinContent(2,j,0.237);
	//	hE.SetBinError(2,j,0.);//0.016);
	//	hE.SetBinContent(3,j,0.241);
	//	hE.SetBinError(3,j,0.);//0.022);
	//}

	fillJetsHist(&hJ);

	//ZFitter zfit("ZFit", &ptScheme, &zyScheme);
	//zfit.setInputs(dataFile);
	//zfit.fit(nZ,eZ);
	//zfit.fixShape();

	SVFitter svfit("ZjFit", &ptScheme, &zyScheme);
	svfit.setInputs(lightDataFile,wmc.outputName(MCJets::jetRecoSV4),wmc.outputName(MCJets::jetRecoSV5),dataFile);
	//svfit.setSVBinning(18,500.,5000.,3);
	svfit.makeSVFitHists(0);
	svfit.makeSVFitHists(4);
	svfit.makeSVFitHists(5);
	svfit.makeSVFitHists(7);

	//SVFitter sssvfit("SSFit", &ptScheme, &zyScheme);
	//sssvfit.setInputs(lightDataFile,wmc.outputName(MCJets::jetRecoSV4),wmc.outputName(MCJets::jetRecoSV5),ssDataFile);
	////svfit.setSVBinning(18,500.,5000.,3);
	//sssvfit.makeSVFitHists(0);
	//sssvfit.makeSVFitHists(4);
	//sssvfit.makeSVFitHists(5);
	//sssvfit.makeSVFitHists(7);

	for(int i=1; i<=npt; ++i) {
		for(int j=1; j<=nzy; ++j) {
			hE.SetBinContent(i,j,hE2.GetBinContent(i,1));
			hE.SetBinError  (i,j,hE2.GetBinError  (i,1));
			if(svfit.fitSV(nB,eB,nC,eC,nQ,eQ,i,j)) {
				double corr = wmc.getPtCorrFactor(MCJets::jetRecoSV4,ptBounds[i-1],ptBounds[i]);
				hC.SetBinContent(i,j,corr*nC);
				hC.SetBinError(  i,j,corr*eC);
			}
			//if(sssvfit.fitSV(nB,eB,nC,eC,nQ,eQ,i,j)) {
			//	hCSS.SetBinContent(i,j,nC);
			//	hCSS.SetBinError(  i,j,eC);
			//}
			//if(zfit.fit(nZ,eZ,i,j)) {
			//	hZ.SetBinContent(i,j,nZ);
			//	hZ.SetBinError(  i,j,eZ);
			//}
		}
	}

	TH2D* hUC = wmc.unfold(&hC, MCJets::jetRecoSV4);//, true);
	TH2D* hUJ = wmc.unfold(&hJ, MCJets::jetAll0);//, true);

	hZ.Divide(&hJ);
	hCSS.Divide(&hC);

	hR.Divide(&hC,&hJ);
	hR.Divide(&hE);
	hUR.Divide(hUC,hUJ);
	hUR.Divide(&hE);

	hCE.Divide(&hC,&hE);
	hUCE.Divide(hUC,&hE);
	TH1D* hCE20up  = hCE .ProjectionY("hCE20up", 2,3);
	TH1D* hJ20up   = hJ  .ProjectionY("hJ20up",  2,3);
	TH1D* hUCE20up = hUCE.ProjectionY("hUCE20up",2,3);
	TH1D* hUJ20up  = hUJ->ProjectionY("hUJ20up", 2,3);

	TH1D* hPtCE  = hCE .ProjectionX("hPtCE");
	TH1D* hPtJ   = hJ  .ProjectionX("hPtJ");
	TH1D* hPtUCE = hUCE.ProjectionX("hPtUCE");
	TH1D* hPtUJ  = hUJ->ProjectionX("hPtUJ");

	hR20up.Divide(hCE20up,hJ20up);
	hUR20up.Divide(hUCE20up,hUJ20up);

	hUJ20up->SetLineColor(kRed);
	hUJ20up->SetMarkerColor(kRed);
	hUCE20up->SetLineColor(kRed);
	hUCE20up->SetMarkerColor(kRed);
	hUCE20up->SetLineStyle(kDashed);
	hCE20up->SetLineStyle(kDashed);
	hUR20up.SetLineColor(kRed);
	hUR20up.SetMarkerColor(kRed);

	hPtUJ->SetLineColor(kRed);
	hPtUJ->SetMarkerColor(kRed);
	hPtUCE->SetLineColor(kRed);
	hPtUCE->SetMarkerColor(kRed);
	hPtUCE->SetLineStyle(kDashed);
	hPtCE->SetLineStyle(kDashed);

	TCanvas c;
	hR20up.SetMaximum(1.1*TMath::Max(hR20up.GetMaximum(),hUR20up.GetMaximum()));
	hR20up.SetMinimum(0.);
	hR20up.GetXaxis()->SetTitle("#it{y}(#it{Z})");
	hR20up.GetYaxis()->SetTitle("#it{f}_{#it{c}}");
	hR20up.Draw();
	hUR20up.Draw("same");
	c.SaveAs(savedir+"/fracC_20up.pdf");
	hJ20up->SetMaximum(1.1*TMath::Max(hJ20up->GetMaximum(),hUJ20up->GetMaximum()));
	hJ20up->SetMinimum(0.);
	hJ20up->GetXaxis()->SetTitle("#it{y}(#it{Z})");
	hJ20up->GetYaxis()->SetTitle("#it{N}_{jet}");
	hJ20up->Draw();
	hUJ20up->Draw("same");
	hCE20up->Draw("same");
	hUCE20up->Draw("same");
	c.SaveAs(savedir+"/NJNC_20up.pdf");
	hPtJ->SetMaximum(1.1*TMath::Max(hPtJ->GetMaximum(),hPtUJ->GetMaximum()));
	hPtJ->SetMinimum(0.);
	hPtJ->GetXaxis()->SetTitle("#it{p}_{T}(#it{j})");
	hPtJ->GetYaxis()->SetTitle("#it{N}_{jet}");
	hPtJ->Draw();
	hPtUJ->Draw("same");
	hPtCE->Draw("same");
	hPtUCE->Draw("same");
	c.SaveAs(savedir+"/NJNC_pT.pdf");

	TH1D* h1520 = hR.ProjectionY("h1520",1,1);
	TH1D* h2030 = hR.ProjectionY("h2030",2,2);
	TH1D* h30up = hR.ProjectionY("h30up",3,3);

	h30up->SetLineColor(kRed);
	h30up->SetMarkerColor(kRed);
	h2030->SetLineColor(kGreen+2);
	h2030->SetMarkerColor(kGreen+2);
	h1520->SetLineColor(kBlue);
	h1520->SetMarkerColor(kBlue);
	h30up->SetMinimum(0.);
	h30up->SetMaximum(0.15);
	gStyle->SetOptStat(0);
	h30up->GetXaxis()->SetTitle("#it{y}_{#it{Z}^{0}}");
	h30up->GetYaxis()->SetTitle("#it{#sigma}_{#it{Z}^{0}+#it{c}}/#it{#sigma}_{#it{Z}^{0}+#it{j}}");
	h30up->Draw("E1 P");
	h2030->Draw("E1 P same");
	h1520->Draw("E1 P same");
	c.SaveAs(savedir+"/fracCharmZj.pdf");

	gStyle->SetOptStat(0);
	hE.Draw("colz text45");
	c.SaveAs(savedir+"/charmEffs.pdf");
	hC.Draw("colz text45");
	c.SaveAs(savedir+"/charmYields.pdf");
	hJ.Draw("colz text45");
	c.SaveAs(savedir+"/totalYields.pdf");
	hR.Draw("colz text45");
	c.SaveAs(savedir+"/charmFractions.pdf");
	//hZ.Draw("colz text45");
	//c.SaveAs(savedir+"/fracZ.pdf");
	hCSS.Draw("colz text45");
	c.SaveAs(savedir+"/fracCharmBkg.pdf");

	for(int i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(int j=1; j<=nzy; ++j) {
			//std::cout << ptBounds[i-1] << "-" << ptBounds[i] << "," << zyBounds[j-1] << "-" << zyBounds[j] << ":\t" << hJ.GetBinContent(i,j) << "+/-" << hJ.GetBinError(i,j) << "\t" << hC.GetBinContent(i,j) << "+/-" << hC.GetBinError(i,j) << "\t" << hE.GetBinContent(i,j) << "+/-" << hE.GetBinError(i,j) << "\t" << hR.GetBinContent(i,j) << "+/-" << hR.GetBinError(i,j) << std::endl;
			printf("%5.0f--%5.0f & %.1f--%.1f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ \\\\\n", ptBounds[i-1], ptBounds[i], zyBounds[j-1], zyBounds[j], hJ.GetBinContent(i,j), hJ.GetBinError(i,j), hC.GetBinContent(i,j), hC.GetBinError(i,j), hE.GetBinContent(i,j), hE.GetBinError(i,j), hR.GetBinContent(i,j), hR.GetBinError(i,j));
		}
	}
	std::cout << "unfolded" << std::endl;
	for(int i=1; i<=npt; ++i) {
		printf("\\midrule\n");
		for(int j=1; j<=nzy; ++j) {
			//std::cout << ptBounds[i-1] << "-" << ptBounds[i] << "," << zyBounds[j-1] << "-" << zyBounds[j] << ":\t" << hUJ->GetBinContent(i,j) << "+/-" << hUJ->GetBinError(i,j) << "\t" << hUC->GetBinContent(i,j) << "+/-" << hUC->GetBinError(i,j) << "\t" << hE.GetBinContent(i,j) << "+/-" << hE.GetBinError(i,j) << "\t" << hUR.GetBinContent(i,j) << "+/-" << hUR.GetBinError(i,j) << std::endl;
			printf("%5.0f--%5.0f & %.1f--%.1f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ \\\\\n", ptBounds[i-1], ptBounds[i], zyBounds[j-1], zyBounds[j], hUJ->GetBinContent(i,j), hUJ->GetBinError(i,j), hUC->GetBinContent(i,j), hUC->GetBinError(i,j), hE.GetBinContent(i,j), hE.GetBinError(i,j), hUR.GetBinContent(i,j), hUR.GetBinError(i,j));
		}
	}

	FILE* rFile = fopen(savedir+"/results.log","w");
	for(int i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(int j=1; j<=nzy; ++j) {
			fprintf(rFile,"%5.0f--%5.0f & %.1f--%.1f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ \\\\\n", ptBounds[i-1], ptBounds[i], zyBounds[j-1], zyBounds[j], hJ.GetBinContent(i,j), hJ.GetBinError(i,j), hC.GetBinContent(i,j), hC.GetBinError(i,j), hE.GetBinContent(i,j), hE.GetBinError(i,j), hR.GetBinContent(i,j), hR.GetBinError(i,j));
		}
	}
	fprintf(rFile,"unfolded\n");
	for(int i=1; i<=npt; ++i) {
		fprintf(rFile,"\\midrule\n");
		for(int j=1; j<=nzy; ++j) {
			fprintf(rFile,"%5.0f--%5.0f & %.1f--%.1f & $%5.0f\\pm%3.0f$ & $%5.0f\\pm%3.0f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ \\\\\n", ptBounds[i-1], ptBounds[i], zyBounds[j-1], zyBounds[j], hUJ->GetBinContent(i,j), hUJ->GetBinError(i,j), hUC->GetBinContent(i,j), hUC->GetBinError(i,j), hE.GetBinContent(i,j), hE.GetBinError(i,j), hUR.GetBinContent(i,j), hUR.GetBinError(i,j));
		}
	}
	fclose(rFile);




	TFile* fout = TFile::Open(savedir+"/hists.root","RECREATE");
	hC.Write();
	hJ.Write();
	hE.Write();
	hR.Write();
	hUC->Write();
	hUJ->Write();
	hUR.Write();
	fout->Close();

	return 0;
}
