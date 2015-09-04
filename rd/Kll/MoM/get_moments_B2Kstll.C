#define get_moments_B2Kstll_cxx
#include "get_moments_B2Kstll.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TComplex.h>
#include <TRandom3.h>
#include <TF1.h>
#include <iostream>
#include <fstream>

using namespace std;

double D000(double phi, double cosTheta) { return 1.; }
double D001(double phi, double cosTheta) { return cosTheta; }
double D002(double phi, double cosTheta) { return (3. *pow(cosTheta,2) - 1.)/2.; }
double D003(double phi, double cosTheta) { return (5. *pow(cosTheta,3) - 3. *cosTheta)/2.; }
double D004(double phi, double cosTheta) { return (35.*pow(cosTheta,4) - 30.*pow(cosTheta,2) + 3.)/8.; }
double D101(double phi, double cosTheta)
{
    double sinTheta = sqrt(1.-pow(cosTheta,2));
    //return -1./sqrt(2.)*TComplex::Exp(-TComplex::I()*phi)*sinTheta;
    return -1./sqrt(2.)*TComplex::Exp(+TComplex::I()*phi)*sinTheta;
}
double D202(double phi, double cosTheta)
{
    //return sqrt(3./8.)*TComplex::Exp(-2.*TComplex::I()*phi)*(1.-pow(cosTheta,2));
    return sqrt(3./8.)*TComplex::Exp(+2.*TComplex::I()*phi)*(1.-pow(cosTheta,2));
}
double D102(double phi, double cosTheta)
{
    double sinTheta = sqrt(1.-pow(cosTheta,2));
    //return -sqrt(3./8.)*TComplex::Exp(-TComplex::I()*phi)*2.*cosTheta*sinTheta;
    return -sqrt(3./8.)*TComplex::Exp(+TComplex::I()*phi)*2.*cosTheta*sinTheta;
}

double Omega000(double phi, double cosThetaK, double cosThetaL) { return 1.; }
double Omega001(double phi, double cosThetaK, double cosThetaL){ return D000(phi, cosThetaK) * D001(0., cosThetaL); }
double Omega002(double phi, double cosThetaK, double cosThetaL){ return D000(phi, cosThetaK) * D002(0., cosThetaL); }
double Omega003(double phi, double cosThetaK, double cosThetaL){ return D000(phi, cosThetaK) * D003(0., cosThetaL); }
double Omega004(double phi, double cosThetaK, double cosThetaL){ return D000(phi, cosThetaK) * D004(0., cosThetaL); }
double Omega020(double phi, double cosThetaK, double cosThetaL){ return D002(phi, cosThetaK) * D000(0., cosThetaL); }
double Omega021(double phi, double cosThetaK, double cosThetaL){ return D002(phi, cosThetaK) * D001(0., cosThetaL); }
double Omega121(double phi, double cosThetaK, double cosThetaL){ return D102(phi, cosThetaK) * D101(0., cosThetaL); }
double Omega022(double phi, double cosThetaK, double cosThetaL){ return D002(phi, cosThetaK) * D002(0., cosThetaL); }
double Omega122(double phi, double cosThetaK, double cosThetaL){ return D102(phi, cosThetaK) * D102(0., cosThetaL); }
double Omega222(double phi, double cosThetaK, double cosThetaL){ return D202(phi, cosThetaK) * D202(0., cosThetaL); }

int delta(int m) {if (m == 0) return 1; else return 0;}
double coeff(int m, int lK, int ll) {return (1 + delta(m))/(2.*(2*lK+1)*(2*ll+1));}

void get_moments_B2Kstll::Loop(TString filename) {
    if (fChain == 0) return;

    outFile = TFile::Open(filename, "RECREATE");
    tup = new TNtuple("partial", "coeffs", "cosThetaL:phi:cosThetaK:D000:D002:D102:D202");
    tup_total = new TNtuple("total", "coeffs",               "G000:G001:G002:G003:G004:G020:G021:G121:G022:G122:G222");
    tup_error = new TNtuple("error", "coeffs",               "G000:G001:G002:G003:G004:G020:G021:G121:G022:G122:G222");
    tup_bootstrapped = new TNtuple("bootstrapped", "coeffs", "G000:G001:G002:G003:G004:G020:G021:G121:G022:G122:G222");

    this->GetMoments();

    for(Int_t i=0; i<1000; ++i) {
        std::cout << "bootstrapped sample " << i << " of 1000..." << std::endl;
        SetIndicesForBootstrapping(i+9000);
        this->GetMoments(true);
    }

    this->GetUncertainties();

    tup->Write();
    tup_total->Write();
    tup_error->Write();
    tup_bootstrapped->Write();
    outFile->Close();
}

void get_moments_B2Kstll::SetIndicesForBootstrapping(int seed)
{
    TRandom3 r(seed);
    for(Int_t i=0; i<nentries; ++i) {
       indices[i] = (Int_t)(r.Rndm()*nentries);
    }

    std::sort(indices.begin(),indices.end());
}

void get_moments_B2Kstll::GetMoments(bool bootstrap)
{
    double G000(0.);
    double G001(0.);
    double G002(0.);
    double G003(0.);
    double G004(0.);
    double G020(0.);
    double G021(0.);
    double G121(0.);
    double G022(0.);
    double G122(0.);
    double G222(0.);

    double D000_v(0.);
    double D002_v(0.);
    double D102_v(0.);
    double D202_v(0.);
    double eff(0.);

    Long64_t nbytes = 0, nb = 0;
    Long64_t ientry(0);
    for (Long64_t jentry=0; jentry<nentries; ++jentry) {
        if(bootstrap) {
            ientry = LoadTree(indices[jentry]);
            if (ientry < 0) break;
            nb = fChain->GetEntry(indices[jentry]);   nbytes += nb;
        } else {
            ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;
        }
        eff = 1. + a*cosThetaL*cosThetaL + b*cosThetaL*cosThetaL*cosThetaL*cosThetaL;
	if(dVeto)   eff *= dVeto->Interpolate(cosThetaL);
	if(psiVeto) eff *= psiVeto->Interpolate(cosThetaL);
        
        //partial moments
        double norm_p = 1./4./TMath::Pi();
        D000_v = D000(phi, cosThetaK)/norm_p;
        D002_v = D002(phi, cosThetaK)/norm_p;
        D102_v = D102(phi, cosThetaK)/norm_p;
        D202_v = D202(phi, cosThetaK)/norm_p;
        if(!bootstrap) tup->Fill(cosThetaL, phi, cosThetaK, D000_v, D002_v, D102_v, D202_v);

	if(!useSWeights_) {
		if(B_M > sigMin_ && B_M < sigMax_) {
			sWeight=sigWeight_;
		} else if(B_M > bkgMin_ && B_M < bkgMax_) {
			sWeight=-1.*bkgWeight_;
		} else {
			continue;
		}
	}
    
        //total moments
        double norm_t = 1./8./TMath::Pi();
        G000 += sWeight*Omega000(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G001 += sWeight*Omega001(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G002 += sWeight*Omega002(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G003 += sWeight*Omega003(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G004 += sWeight*Omega004(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G020 += sWeight*Omega020(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G021 += sWeight*Omega021(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G121 += sWeight*Omega121(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G022 += sWeight*Omega022(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G122 += sWeight*Omega122(phi, cosThetaK, cosThetaL)/norm_t/eff;
        G222 += sWeight*Omega222(phi, cosThetaK, cosThetaL)/norm_t/eff;
    }

    cout << "All normalised to G000 " << G000/coeff(0,0,0) << endl;
    cout << "G000 \t " << G000/G000/2./coeff(0,0,0) << endl;
    cout << "G001 \t " << G001/G000/2./coeff(0,0,1) << endl;
    cout << "G002 \t " << G002/G000/2./coeff(0,0,2) << endl;
    cout << "G003 \t " << G003/G000/2./coeff(0,0,3) << endl;
    cout << "G004 \t " << G004/G000/2./coeff(0,0,4) << endl;
    cout << "G020 \t " << G020/G000/2./coeff(0,2,0) << endl;
    cout << "G021 \t " << G021/G000/2./coeff(0,2,1) << endl;
    cout << "G121 \t " << G121/G000/2./coeff(1,2,1) << endl;
    cout << "G022 \t " << G022/G000/2./coeff(0,2,2) << endl;
    cout << "G122 \t " << G122/G000/2./coeff(1,2,2) << endl;
    cout << "G222 \t " << G222/G000/2./coeff(2,2,2) << endl;

    if(bootstrap) {
        tup_bootstrapped->Fill(G000/G000/2./coeff(0,0,0), 
                               G001/G000/2./coeff(0,0,1),
                               G002/G000/2./coeff(0,0,2),
                               G003/G000/2./coeff(0,0,3),
                               G004/G000/2./coeff(0,0,4),
                               G020/G000/2./coeff(0,2,0),
                               G021/G000/2./coeff(0,2,1),
                               G121/G000/2./coeff(1,2,1),
                               G022/G000/2./coeff(0,2,2),
                               G122/G000/2./coeff(1,2,2),
                               G222/G000/2./coeff(2,2,2)
                               );
    } else {
        tup_total->Fill(G000/G000/2./coeff(0,0,0), 
                        G001/G000/2./coeff(0,0,1),
                        G002/G000/2./coeff(0,0,2),
                        G003/G000/2./coeff(0,0,3),
                        G004/G000/2./coeff(0,0,4),
                        G020/G000/2./coeff(0,2,0),
                        G021/G000/2./coeff(0,2,1),
                        G121/G000/2./coeff(1,2,1),
                        G022/G000/2./coeff(0,2,2),
                        G122/G000/2./coeff(1,2,2),
                        G222/G000/2./coeff(2,2,2)
                        );
    }

}

void get_moments_B2Kstll::GetUncertainties()
{
    double errors[11];

    Int_t nVars = 11;

    for(Int_t i=0; i<nVars; ++i) {
        TString varName = tup_total->GetListOfBranches()->At(i)->GetName();
        Double_t min = tup_bootstrapped->GetMinimum(varName);
        Double_t max = tup_bootstrapped->GetMaximum(varName);
	Double_t range = max-min;

	if(range>0) {
            min -= 0.1*range;
	    max += 0.1*range;
	    range *= 1.2;

	    TH1D h("h","",100,min,max);
	    tup_bootstrapped->Draw(varName+">>h");

	    while(range > 10*h.GetRMS()) {
                range = 10*h.GetRMS();
		min = h.GetBinCenter(h.GetMaximumBin()) - range/2.;
		max = h.GetBinCenter(h.GetMaximumBin()) + range/2.;

		h = TH1D("h","",100,min,max);
	    	tup_bootstrapped->Draw(varName+">>h");
	    }

            TF1 gaus("gaus","gaus(0)",min,max);
            gaus.SetParameter(0,1000);
            gaus.SetParameter(1,h.GetBinCenter(h.GetMaximumBin()));
            gaus.SetParameter(2,h.GetStdDev());
            gaus.SetParLimits(1,min+0.2*range,max-0.2*range);
            gaus.SetParLimits(2,h.GetStdDev()/3.,h.GetStdDev()*3.);
     
	    TCanvas c;
            TFitResultPtr r = h.Fit(&gaus,"S");
            h.GetFunction("gaus")->SetLineColor(kRed);
            h.GetXaxis()->SetTitle(varName);
            h.GetYaxis()->SetTitle("N_{toys}");
            h.Draw();
            c.SaveAs("plots/error_"+varName+".pdf");
      
            errors[i] = h.GetFunction("gaus")->GetParameter(2);
	} else {
	    errors[i] = 0;
	}

    }

    tup_error->Fill(errors[ 0], 
                    errors[ 1],
                    errors[ 2],
                    errors[ 3],
                    errors[ 4],
                    errors[ 5],
                    errors[ 6],
                    errors[ 7],
                    errors[ 8],
                    errors[ 9],
                    errors[10]
                    );
    
}

int main( int argc, char * argv[] )
{
    TString inputFile("");
    TString outputFile("results.root");
    TString treeName("DecayTree");
    TString effFile(""), dVetoFile(""), psiVetoFile("");
    Int_t effEntry(0);
    Double_t a(-0.5), b(0.0);
    TH1 *hDVeto(0), *hPsiVeto(0);
    TString bkgFile("");
    Bool_t useSWeights(true);
    Double_t sigWeight(1.);
    Double_t bkgWeight(-1.);
    Double_t sigMin(5170.);
    Double_t sigMax(5400.);
    Double_t bkgMin(5400.);
    Double_t bkgMax(5970.);
	
    if(argc < 2 ) {
        std::cout << "Usage: " << argv[0] << " inputFile [outputFile = \"results.root\"] [inputTree = \"DecayTree\"]" << std::endl;
    } else {
        inputFile = argv[1];
        if(argc > 3) {
	    outputFile = argv[2];
            if(argc > 3) {
	        treeName = argv[3];
            	if(argc > 4) {
	    	    effFile = argv[4];
            		if(argc > 5) {
	    		    effEntry = atoi( argv[5] );
            			if(argc > 7) {
	    	    			dVetoFile = argv[6];
	    	    			psiVetoFile = argv[7];
            				if(argc > 8) {
						bkgFile = argv[8];
					}
				}
			}
		}
	    }
	}
    }
    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TFile * f = TFile::Open(inputFile);
    TTree * tree = dynamic_cast<TTree*>(f->Get(treeName));
    if(effFile != "") {
        if(effFile=="NOEFF") {
	    a=0.0;
	    b=0.0;
	} else {
    	    TFile * fEff = TFile::Open(effFile);
    	    TTree * tEff = dynamic_cast<TTree*>(fEff->Get("params"));
	    tEff->SetBranchAddress("A",&a);
	    tEff->SetBranchAddress("B",&b);
	    tEff->GetEntry(effEntry);
	}
    }
    if(dVetoFile != "") {
	    TString dVetoHistName("efficiency_"); dVetoHistName+=effEntry;
	    TFile * fDVeto = TFile::Open(dVetoFile);
	    hDVeto = dynamic_cast<TH1*>(fDVeto->Get(dVetoHistName));
    }
    if(psiVetoFile != "") {
	    TString psiVetoHistName("efficiency_"); psiVetoHistName+=effEntry;
	    TFile * fPsiVeto = TFile::Open(psiVetoFile);
	    hPsiVeto = dynamic_cast<TH1*>(fPsiVeto->Get(psiVetoHistName));
    }
    if(bkgFile != "") {
	    useSWeights=false;

	    std::ifstream fin;
	    fin.open(bkgFile);
	    fin >> sigMin >> sigMax >> bkgWeight;
	    fin.close();
    }
    get_moments_B2Kstll t(tree,a,b,hDVeto,hPsiVeto,useSWeights,sigWeight,bkgWeight,sigMin,sigMax,bkgMin,bkgMax);
    t.Loop(outputFile);
    return 1;
}

