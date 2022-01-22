#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <string>
#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "dataset/weights/BXSignal_NN.class.C"

void applyTMVA(){
    // Load the library
    TMVA::Tools::Instance();
    
    // Initialize the reader
    TMVA::Reader* reader = new TMVA::Reader();

    // Create local variables to store updated inputs
    Float_t b_fd(0.);
    Float_t logAcosb_dira(0.);
    Float_t b_vtx_chi2(0.);
    Float_t b_doca_kpi(0.);
    Float_t b_p(0.);
    Float_t logb_pt(0.);
    Float_t b_ip(0.);
    Float_t logb_ip_chi2(0.);
    Float_t b_tau(0.);
    Float_t b_tau_chi2(0.);
    Float_t b_maxdoca(0.);
    Float_t x_fd(0.);
    Float_t x_vtx_chi2(0.);
    Float_t x_p(0.);
    Float_t logx_pt(0.);
    Float_t x_ip(0.);
    Float_t logx_ip_chi2(0.);
    Float_t x_tau(0.);
    Float_t x_tau_chi2(0.);
    Float_t x_maxdoca(0.);
    Float_t logk_ip_chi2(0.);
    Float_t logk_pt(0.);
    Float_t logpi_ip_chi2(0.);
    Float_t logpi_pt(0.);
    Float_t logxpip_ip_chi2(0.);
    Float_t logxpip_pt(0.);
    Float_t logxpim_ip_chi2(0.);
    Float_t logxpim_pt(0.);
    Float_t logxeta_pt(0.);

    // Add variables used to make cuts
    reader->AddVariable("b_fd", &b_fd);
    reader->AddVariable("log(acos(b_dira))", &logAcosb_dira);
    reader->AddVariable("b_vtx_chi2", &b_vtx_chi2);
    reader->AddVariable("b_doca_kpi", &b_doca_kpi);
    reader->AddVariable("b_p", &b_p);
    reader->AddVariable("log(b_pt)", &logb_pt);
    reader->AddVariable("b_ip", &b_ip);
    reader->AddVariable("log(b_ip_chi2)", &logb_ip_chi2);
    reader->AddVariable("b_tau", &b_tau);
    reader->AddVariable("b_tau_chi2", &b_tau_chi2);
    reader->AddVariable("b_maxdoca", &b_maxdoca);
    reader->AddVariable("x_fd", &x_fd);
    reader->AddVariable("x_vtx_chi2", &x_vtx_chi2);
    reader->AddVariable("x_p", &x_p);
    reader->AddVariable("log(x_pt)", &logx_pt);
    reader->AddVariable("x_ip", &x_ip);
    reader->AddVariable("log(x_ip_chi2)", &logx_ip_chi2);
    reader->AddVariable("x_tau", &x_tau);
    reader->AddVariable("x_tau_chi2", &x_tau_chi2);
    reader->AddVariable("x_maxdoca", &x_maxdoca);
    reader->AddVariable("log(k_ip_chi2)", &logk_ip_chi2);
    reader->AddVariable("log(k_pt)", &logk_pt);
    reader->AddVariable("log(pi_ip_chi2)", &logpi_ip_chi2);
    reader->AddVariable("log(pi_pt)", &logpi_pt);
    reader->AddVariable("log(xpip_ip_chi2)", &logxpip_ip_chi2);
    reader->AddVariable("log(xpip_pt)", &logxpip_pt);
    reader->AddVariable("log(xpim_ip_chi2)", &logxpim_ip_chi2);
    reader->AddVariable("log(xpim_pt)", &logxpim_pt);
    reader->AddVariable("log(xeta_pt)", &logxeta_pt);

    // Book method
    reader->BookMVA("NN","dataset/weights/BXSignal_NN.weights.xml");

    // Get classifier response and its error
    /*cout << "Booked Method\n";
    Float_t cValue = reader->EvaluateMVA("NN");
    cout << "Evaluated Method\n";
    Float_t cError = reader->GetMVAError();
    cout << "Got MVA Error\n";*/

    // Create output histogram
    TH1F *cutHisto = new TH1F( "NN", "NN", 1000, 0.,1.);

    //// Open input file
    //TFile *input = TFile::Open("/data/alps/4/kpietapipi.root");

    //// Prepare event tree
    //TTree* applyTree;
    //input->GetObject("T;l",applyTree);

    TChain c("T");
    c.AddFile("/data/alps/4/kpietapipi.root");
    c.AddFile("/data/alps/5/kpietapipi.root");
    c.AddFile("/data/alps/6/kpietapipi.root");
    c.AddFile("/data/alps/7/kpietapipi.root");
    c.AddFile("/data/alps/8/kpietapipi.root");
    c.AddFile("/data/alps/9/kpietapipi.root");
    c.AddFile("/data/alps/10/kpietapipi.root");
    c.AddFile("/data/alps/11/kpietapipi.root");

    Double_t tree_b_fd(0.);
    Double_t tree_b_dira(0.);
    Double_t tree_b_vtx_chi2(0.);
    Double_t tree_b_doca_kpi(0.);
    Double_t tree_b_p(0.);
    Double_t tree_b_pt(0.);
    Double_t tree_b_ip(0.);
    Double_t tree_b_ip_chi2(0.);
    Double_t tree_b_tau(0.);
    Double_t tree_b_tau_chi2(0.);
    Double_t tree_b_maxdoca(0.);
    Double_t tree_x_fd(0.);
    Double_t tree_x_vtx_chi2(0.);
    Double_t tree_x_p(0.);
    Double_t tree_x_pt(0.);
    Double_t tree_x_ip(0.);
    Double_t tree_x_ip_chi2(0.);
    Double_t tree_x_tau(0.);
    Double_t tree_x_tau_chi2(0.);
    Double_t tree_x_maxdoca(0.);
    Double_t tree_k_ip_chi2(0.);
    Double_t tree_k_pt(0.);
    Double_t tree_pi_ip_chi2(0.);
    Double_t tree_pi_pt(0.);
    Double_t tree_xpip_ip_chi2(0.);
    Double_t tree_xpip_pt(0.);
    Double_t tree_xpim_ip_chi2(0.);
    Double_t tree_xpim_pt(0.);
    Double_t tree_xeta_pt(0.);

    c.SetBranchAddress("b_fd",        &tree_b_fd);
    c.SetBranchAddress("b_dira",      &tree_b_dira);
    c.SetBranchAddress("b_vtx_chi2",  &tree_b_vtx_chi2);
    c.SetBranchAddress("b_doca_kpi",  &tree_b_doca_kpi);
    c.SetBranchAddress("b_p",         &tree_b_p);
    c.SetBranchAddress("b_pt",        &tree_b_pt);
    c.SetBranchAddress("b_ip",        &tree_b_ip);
    c.SetBranchAddress("b_ip_chi2",   &tree_b_ip_chi2);
    c.SetBranchAddress("b_tau",       &tree_b_tau);
    c.SetBranchAddress("b_tau_chi2",  &tree_b_tau_chi2);
    c.SetBranchAddress("b_maxdoca",   &tree_b_maxdoca);
    c.SetBranchAddress("x_fd",        &tree_x_fd);
    c.SetBranchAddress("x_vtx_chi2",  &tree_x_vtx_chi2);
    c.SetBranchAddress("x_p",         &tree_x_p);
    c.SetBranchAddress("x_pt",        &tree_x_pt);
    c.SetBranchAddress("x_ip",        &tree_x_ip);
    c.SetBranchAddress("x_ip_chi2",   &tree_x_ip_chi2);
    c.SetBranchAddress("x_tau",       &tree_x_tau);
    c.SetBranchAddress("x_tau_chi2",  &tree_x_tau_chi2);
    c.SetBranchAddress("x_maxdoca",   &tree_x_maxdoca);
    c.SetBranchAddress("k_ip_chi2",   &tree_k_ip_chi2);
    c.SetBranchAddress("k_pt",        &tree_k_pt);
    c.SetBranchAddress("pi_ip_chi2",  &tree_pi_ip_chi2);
    c.SetBranchAddress("pi_pt",       &tree_pi_pt);
    c.SetBranchAddress("xpip_ip_chi2",&tree_xpip_ip_chi2);
    c.SetBranchAddress("xpip_pt",     &tree_xpip_pt);
    c.SetBranchAddress("xpim_ip_chi2",&tree_xpim_ip_chi2);
    c.SetBranchAddress("xpim_pt",     &tree_xpim_pt);
    c.SetBranchAddress("xeta_pt",     &tree_xeta_pt);

    TFile* fout = TFile::Open("/data/alps/kpietapipi_addNN.root","RECREATE");
    TTree* tout = c.CloneTree(0);
    double NN(0.);
    tout->Branch("NN", &NN);

    // Loop over events
    int n=c.GetEntries();
    boost::progress_display progress(n);
    for (int i=0; i<n; i++)
    {
	++progress;
        c.GetEntry(i);

    	b_fd =           tree_b_fd;
    	logAcosb_dira =  log(acos(tree_b_dira));
    	b_vtx_chi2 =     tree_b_vtx_chi2;
    	b_doca_kpi =     tree_b_doca_kpi;
    	b_p =            tree_b_p;
    	logb_pt =        log(tree_b_pt);
    	b_ip =           tree_b_ip;
    	logb_ip_chi2 =   log(TMath::Max(tree_b_ip_chi2,0.));
    	b_tau =          tree_b_tau;
    	b_tau_chi2 =     tree_b_tau_chi2;
    	b_maxdoca =      tree_b_maxdoca;
    	x_fd =           tree_x_fd;
    	x_vtx_chi2 =     tree_x_vtx_chi2;
    	x_p =            tree_x_p;
    	logx_pt =        log(tree_x_pt);
    	x_ip =           tree_x_ip;
    	logx_ip_chi2 =   log(TMath::Max(tree_x_ip_chi2,0.));
    	x_tau =          tree_x_tau;
    	x_tau_chi2 =     tree_x_tau_chi2;
    	x_maxdoca =      tree_x_maxdoca;
    	logk_ip_chi2 =   log(TMath::Max(tree_k_ip_chi2,0.));
    	logk_pt =        log(tree_k_pt);
    	logpi_ip_chi2 =  log(TMath::Max(tree_pi_ip_chi2,0.));
    	logpi_pt =       log(tree_pi_pt);
    	logxpip_ip_chi2 =log(TMath::Max(tree_xpip_ip_chi2,0.));
    	logxpip_pt =     log(tree_xpip_pt);
    	logxpim_ip_chi2 =log(TMath::Max(tree_xpim_ip_chi2,0.));
    	logxpim_pt =     log(tree_xpim_pt);
    	logxeta_pt =     log(tree_xeta_pt);

	//if(i>639) {
	//std::cout << b_fd <<"	"<<           tree_b_fd<<std::endl;
    	//std::cout << logAcosb_dira <<"	"<<  tree_b_dira<<std::endl;
    	//std::cout << b_vtx_chi2 <<"	"<<     tree_b_vtx_chi2<<std::endl;
    	//std::cout << b_doca_kpi <<"	"<<     tree_b_doca_kpi<<std::endl;
    	//std::cout << b_p <<"	"<<            tree_b_p<<std::endl;
    	//std::cout << logb_pt <<"	"<<        tree_b_pt<<std::endl;
    	//std::cout << b_ip <<"	"<<           tree_b_ip<<std::endl;
    	//std::cout << logb_ip_chi2 <<"	"<<   tree_b_ip_chi2<<std::endl;
    	//std::cout << b_tau <<"	"<<          tree_b_tau<<std::endl;
    	//std::cout << b_tau_chi2 <<"	"<<     tree_b_tau_chi2<<std::endl;
    	//std::cout << b_maxdoca <<"	"<<      tree_b_maxdoca<<std::endl;
    	//std::cout << x_fd <<"	"<<           tree_x_fd<<std::endl;
    	//std::cout << x_vtx_chi2 <<"	"<<     tree_x_vtx_chi2<<std::endl;
    	//std::cout << x_p <<"	"<<            tree_x_p<<std::endl;
    	//std::cout << logx_pt <<"	"<<        tree_x_pt<<std::endl;
    	//std::cout << x_ip <<"	"<<           tree_x_ip<<std::endl;
    	//std::cout << logx_ip_chi2 <<"	"<<   tree_x_ip_chi2<<std::endl;
    	//std::cout << x_tau <<"	"<<          tree_x_tau<<std::endl;
    	//std::cout << x_tau_chi2 <<"	"<<     tree_x_tau_chi2<<std::endl;
    	//std::cout << x_maxdoca <<"	"<<      tree_x_maxdoca<<std::endl;
    	//std::cout << logk_ip_chi2 <<"	"<<   tree_k_ip_chi2<<std::endl;
    	//std::cout << logk_pt <<"	"<<        tree_k_pt<<std::endl;
    	//std::cout << logpi_ip_chi2 <<"	"<<  tree_pi_ip_chi2<<std::endl;
    	//std::cout << logpi_pt <<"	"<<       tree_pi_pt<<std::endl;
    	//std::cout << logxpip_ip_chi2 <<"	"<<tree_xpip_ip_chi2<<std::endl;
    	//std::cout << logxpip_pt <<"	"<<     tree_xpip_pt<<std::endl;
    	//std::cout << logxpim_ip_chi2 <<"	"<<tree_xpim_ip_chi2<<std::endl;
    	//std::cout << logxpim_pt <<"	"<<     tree_xpim_pt<<std::endl;
    	//std::cout << logxeta_pt <<"	"<<     tree_xeta_pt<<std::endl;
	//}

        //Fill histogram
        NN = reader->EvaluateMVA("NN");
        cutHisto->Fill(NN);
	if(NN>0.01) tout->Fill();
    }

    TCanvas can;
    cutHisto->Draw();
    can.SetLogy();
    can.SaveAs("NNdata.pdf");

    tout->Write();
    fout->Close();

    delete reader;
}

int main() {
	applyTMVA();
}
