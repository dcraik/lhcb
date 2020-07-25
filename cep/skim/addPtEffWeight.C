#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TDecompChol.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"
#include "TSystem.h"

//#include <boost/progress.hpp>

double func(double* abs, double* par) {
	double x = abs[0];
	double x0 = par[0];
	double y0 = par[1];
	double a1 = par[2];
	double a2 = par[3];

	if(x<x0) return y0 + (x - x0)*a1;
	else     return y0 + (x - x0)*a2;
}

double invfunc(double* abs, double* par) {
	return 1./func(abs,par);
}

TVectorD applyFluctuations(TVectorD v, TMatrixD& L) {
	int N = v.GetNrows();
	if(N!=L.GetNrows() || N!=L.GetNcols()) {
		std::cout << "Nrows mis-match" << std::endl;
		return v;
	}
	//v.Print();
	TVectorD s(N);
	for(int i=0; i<N; ++i) {
		s[i] = gRandom->Gaus();
	}
	//s.Print();
	s*=L;
	//s.Print();
	s+=v;
	//s.Print();
	return s;
}

TF1* fitEff() {
	TFile* f = TFile::Open("/data/cep-phi/phi-sim.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

	TH1D den("den","",100,0.,2.);
	TH1D num("num","",100,0.,2.);
	TH1D eff("eff","",100,0.,2.);
	TH1D cor("cor","",100,0.,2.);
	TH1D min("min","",100,0.,2.);
	TH1D max("max","",100,0.,2.);

	den.Sumw2();
	num.Sumw2();
	eff.Sumw2();
	cor.Sumw2();

	t->Draw("phi_TRUEPT**2/1e6>>den");
	t->Draw("phi_PT**2/1e6>>num","phi_PT>0. && Kplus_P>5e3 && Kminus_P>5e3");

	eff.Divide(&num,&den);
	cor.Divide(&den,&num);

	TF1 *f1 = new TF1("f1",func,0.,2.,4);//used to fit efficiency
	f1->SetParameters(0.3,0.1,0.1,0.01);
	f1->SetParNames("x0","y0","a1","a2");
	TFitResultPtr r = eff.Fit(f1,"S");
	r->GetCovarianceMatrix().Print();
	r->GetCorrelationMatrix().Print();

	//produce random sampling of parameters with correlations present
	TDecompChol dc(r->GetCovarianceMatrix());
	dc.Decompose();
	TMatrixD L = dc.GetU();
	L.Print();
	TMatrixD LT(4,4);
	LT.Transpose(L);
	LT.Print();

	TH1D hpar0("hpar0", "", 100, f1->GetParameter(0)+5*f1->GetParError(0), f1->GetParameter(0)-5*f1->GetParError(0));
	TH1D hpar1("hpar1", "", 100, f1->GetParameter(1)+5*f1->GetParError(1), f1->GetParameter(1)-5*f1->GetParError(1));
	TH1D hpar2("hpar2", "", 100, f1->GetParameter(2)+5*f1->GetParError(2), f1->GetParameter(2)-5*f1->GetParError(2));
	TH1D hpar3("hpar3", "", 100, f1->GetParameter(3)+5*f1->GetParError(3), f1->GetParameter(3)-5*f1->GetParError(3));

	TF1 *f2 = new TF1("f2",invfunc,0.,2.,4);//used to calculate fluctuations
	for(int i=0; i<1000; ++i) {
		//initialise vector to fit results
		TVectorD p(4);
		for(int j=0; j<4; ++j) {
			p[j] = f1->GetParameter(j);
		}
		//apply correlated fluctuations to parameters array
		//p.Print();
		p = applyFluctuations(p, LT);
		//p.Print();
		//fill histograms
		hpar0.Fill(p[0]);
		hpar1.Fill(p[1]);
		hpar2.Fill(p[2]);
		hpar3.Fill(p[3]);
		for(int j=0; j<4; ++j) {
			f2->SetParameter(j,p[j]);
		}
		for(int j=1; j<=min.GetNbinsX(); ++j) {
			double val = f2->Eval(min.GetBinCenter(j));
			if(i==0 || val<min.GetBinContent(j)) min.SetBinContent(j,val);
			if(i==0 || val>max.GetBinContent(j)) max.SetBinContent(j,val);
		}
	}
	TCanvas c;
	c.Divide(2,2);
	c.cd(1);
	hpar0.Draw();
	c.cd(2);
	hpar1.Draw();
	c.cd(3);
	hpar2.Draw();
	c.cd(4);
	hpar3.Draw();
	c.SaveAs("hpars.pdf");
	c.Clear();

	//setup the inverse function to plot with correction
	TF1 *f3 = new TF1("f3",invfunc,0.,2.,4);//used to show fit to efficiency on correction
	TF1 *f4 = new TF1("f4",invfunc,0.,2.,4);//used to fit to correction data
	for(int i=0; i<4; ++i) {
		f3->FixParameter(i,f1->GetParameter(i));
		f4->SetParameter(i,f1->GetParameter(i));
	}
	f3->SetLineColor(kBlue);

	TFitResultPtr r2 = cor.Fit(f4,"S");
	r2->GetCovarianceMatrix().Print();
	r2->GetCorrelationMatrix().Print();

	eff.Draw();
	c.SaveAs("fitPhiEffPtSq.pdf");
	cor.Draw();
	f3->Draw("same");
	min.Draw("same");
	max.Draw("same");
	c.SaveAs("fitPhiEffPtSqCorFunc.pdf");

	//for(int j=1; j<=min.GetNbinsX(); ++j) {
	//	printf("%.1f-%.1f\t",min.GetBinContent(j), max.GetBinContent(j));
	//}
	
	return f4;
}

void addPtEff(TString dataset, TString mode, std::vector<TString>& inputs) {
	TChain* c = new TChain("DecayTree");
	for(std::vector<TString>::iterator it = inputs.begin(); it!=inputs.end(); ++it) {
		for(int i=0; i<3000; ++i) {
			TString name="/data/cep-phi/";
			name+=(*it);
			name+="/";
			name+=i;
			name+="/skimmed";
			if(mode!="phi") name+="-"+mode;
			name+=".root";
			if(!gSystem->AccessPathName(name)) {
				c->AddFile(name);
			}
		}
	}
	c->SetBranchStatus("*",0);
	c->SetBranchStatus("phi_M",1);
	c->SetBranchStatus("phi_PT",1);
	c->SetBranchStatus("log_hrc_fom_v4",1);
	c->SetBranchStatus("L0DUTCK",1);
	c->SetBranchStatus("HLT1TCK",1);
	c->SetBranchStatus("HLT2TCK",1);
	c->SetBranchStatus("Kplus_PT",1);
	c->SetBranchStatus("Kminus_PT",1);
	c->SetBranchStatus("nLongTracks",1);
	c->SetBranchStatus("nVeloTracks",1);
	if(dataset=="pPb" || dataset=="Pbp") {
		c->SetBranchStatus("pPb_l0_flags",1);
		c->SetBranchStatus("pPb_hlt1_flags",1);
		c->SetBranchStatus("pPb_hlt2_flags",1);
	}
	if(dataset=="PbPb") {
		c->SetBranchStatus("PbPb_l0_flags",1);
		c->SetBranchStatus("PbPb_hlt1_flags",1);
		c->SetBranchStatus("PbPb_hlt2_flags",1);
	}

	double phi_PT(0.);

	c->SetBranchAddress("phi_PT", &phi_PT);

	TFile* fout = TFile::Open("/data/cep-phi/"+dataset+"2"+mode+".root","RECREATE");
	TTree* tout = dynamic_cast<TTree*>(c->CloneTree(0));

	double phi_PTSq(0.), weight(1.);

	tout->Branch("phi_PTSq", &phi_PTSq);
	tout->Branch("weight",   &weight);

	//eff-weight function
	TF1* func = fitEff();

	//boost::progress_display progress( c->GetEntries() );
	for(int i=0; i<c->GetEntries(); ++i) {
		//++progress;
		c->GetEntry(i);
		phi_PTSq = phi_PT*phi_PT/1e6;
		weight = func->Eval(phi_PTSq);
		tout->Fill();
	}
	tout->AutoSave();
	fout->Close();
}

int main(int argc, char** argv) {
	TString dataset, mode;
	if(argc<4) {
		std::cout << "usage: " << argv[0] << " <outputName> <mode> <jobNo> [<jobNo> ...]" << std::endl;
		return 1;
	}

	dataset = argv[1];
	mode = argv[2];
	std::vector<TString> inputs;

	for(int i=3; i<argc; ++i) {
		inputs.push_back(argv[i]);
	}

	addPtEff(dataset, mode, inputs);

	return 0;
}
