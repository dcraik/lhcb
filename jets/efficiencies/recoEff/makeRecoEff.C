#define makeRecoEff_cxx
#include "makeRecoEff.h"

#include <iostream>

#include <TCanvas.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "TMatrixD.h"
#include <TRandom3.h>
#include <TStyle.h>
#include "TVectorT.h"

void makeRecoEff::Loop(int nmax)
{
	if (fChain == 0) return;

	int np(50), neta(25);

	double* pbins   = new double[np+1  ]{5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,21000.,22000.,23000.,24000.,25000.,26000.,27000.,28000.,29000.,30000.,32500.,35000.,37500.,40000.,45000.,50000.,60000.,70000.,80000.,90000.,100000.,150000.,200000.,250000.,300000.,400000.,750000.};
	double* etabins = new double[neta+1]{2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5};

	TH2D denomK  ("denomK"  , "", np, pbins, neta, etabins);
	TH2D denomKp ("denomKp" , "", np, pbins, neta, etabins);
	TH2D denomKm ("denomKm" , "", np, pbins, neta, etabins);

	TH2D denomPi ("denomPi" , "", np, pbins, neta, etabins);
	TH2D denomPip("denomPip", "", np, pbins, neta, etabins);
	TH2D denomPim("denomPim", "", np, pbins, neta, etabins);

	TH2D denomP  ("denomP" , "", np, pbins, neta, etabins);
	TH2D denomPp ("denomPp", "", np, pbins, neta, etabins);
	TH2D denomPm ("denomPm", "", np, pbins, neta, etabins);

	TH2D denomMu ("denomMu" , "", np, pbins, neta, etabins);
	TH2D denomMup("denomMup", "", np, pbins, neta, etabins);
	TH2D denomMum("denomMum", "", np, pbins, neta, etabins);

	TH2D denomE  ("denomE" ,  "", np, pbins, neta, etabins);
	TH2D denomEp ("denomEp",  "", np, pbins, neta, etabins);
	TH2D denomEm ("denomEm",  "", np, pbins, neta, etabins);

	TH2D numK    ("numK"    , "", np, pbins, neta, etabins);
	TH2D numKp   ("numKp"   , "", np, pbins, neta, etabins);
	TH2D numKm   ("numKm"   , "", np, pbins, neta, etabins);

	TH2D numPi   ("numPi"   , "", np, pbins, neta, etabins);
	TH2D numPip  ("numPip"  , "", np, pbins, neta, etabins);
	TH2D numPim  ("numPim"  , "", np, pbins, neta, etabins);

	TH2D numP    ("numP"    , "", np, pbins, neta, etabins);
	TH2D numPp   ("numPp"   , "", np, pbins, neta, etabins);
	TH2D numPm   ("numPm"   , "", np, pbins, neta, etabins);

	TH2D numMu   ("numMu"   , "", np, pbins, neta, etabins);
	TH2D numMup  ("numMup"  , "", np, pbins, neta, etabins);
	TH2D numMum  ("numMum"  , "", np, pbins, neta, etabins);

	TH2D numE    ("numE"    , "", np, pbins, neta, etabins);
	TH2D numEp   ("numEp"   , "", np, pbins, neta, etabins);
	TH2D numEm   ("numEm"   , "", np, pbins, neta, etabins);

	TH2D effK    ("effK"    , "", np, pbins, neta, etabins);
	TH2D effKp   ("effKp"   , "", np, pbins, neta, etabins);
	TH2D effKm   ("effKm"   , "", np, pbins, neta, etabins);

	TH2D effPi   ("effPi"   , "", np, pbins, neta, etabins);
	TH2D effPip  ("effPip"  , "", np, pbins, neta, etabins);
	TH2D effPim  ("effPim"  , "", np, pbins, neta, etabins);

	TH2D effP    ("effP"    , "", np, pbins, neta, etabins);
	TH2D effPp   ("effPp"   , "", np, pbins, neta, etabins);
	TH2D effPm   ("effPm"   , "", np, pbins, neta, etabins);

	TH2D effMu   ("effMu"   , "", np, pbins, neta, etabins);
	TH2D effMup  ("effMup"  , "", np, pbins, neta, etabins);
	TH2D effMum  ("effMum"  , "", np, pbins, neta, etabins);

	TH2D effE    ("effE"    , "", np, pbins, neta, etabins);
	TH2D effEp   ("effEp"   , "", np, pbins, neta, etabins);
	TH2D effEm   ("effEm"   , "", np, pbins, neta, etabins);

	numK.Sumw2();
	numKp.Sumw2();
	numKm.Sumw2();
	numPi.Sumw2();
	numPip.Sumw2();
	numPim.Sumw2();
	numP.Sumw2();
	numPp.Sumw2();
	numPm.Sumw2();
	numMu.Sumw2();
	numMup.Sumw2();
	numMum.Sumw2();
	numE.Sumw2();
	numEp.Sumw2();
	numEm.Sumw2();

	Long64_t nentries = fChain->GetEntries();
	std::cout << nentries << " events..." << std::endl;
	if(nmax>=0 && nmax<nentries) {
		nentries=nmax;
	}
	std::cout << "Will process first " << nentries << std::endl;

	boost::progress_display progress( nentries );
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		++progress;
		fChain->GetEntry(jentry);

		double ptot(0.), pt(0.), eta(0.);

		for(unsigned int igen=0; igen<gen_pid->size(); ++igen) {
			int pid = gen_pid->at(igen);
			if(TMath::Abs(pid)!=211 && TMath::Abs(pid)!=321 && TMath::Abs(pid)!=2212 && TMath::Abs(pid)!=11 && TMath::Abs(pid)!=13) continue;

			//remove anything far from the PV (material, KS->pipi, etc)
			if(TMath::Abs(gen_z->at(igen))>200) continue;
			if(TMath::Abs(gen_x->at(igen))>3||TMath::Abs(gen_y->at(igen))>3) continue;

			TVector3 p3(gen_px->at(igen),gen_py->at(igen),gen_pz->at(igen));
			ptot = p3.Mag();
			pt   = p3.Pt();
			eta  = p3.Eta();

			//geometric acceptance
			if(ptot<5000. || pt<500 || eta<2. || eta>4.5) continue;

			switch(pid) {
				case -2212:
					denomP .Fill(ptot,eta);
					denomPm.Fill(ptot,eta);
					break;
				case -321:
					denomK .Fill(ptot,eta);
					denomKm.Fill(ptot,eta);
					break;
				case -211:
					denomPi .Fill(ptot,eta);
					denomPim.Fill(ptot,eta);
					break;
				case -13:
					denomMu .Fill(ptot,eta);
					denomMum.Fill(ptot,eta);
					break;
				case -11:
					denomE .Fill(ptot,eta);
					denomEm.Fill(ptot,eta);
					break;
				case 11:
					denomE .Fill(ptot,eta);
					denomEp.Fill(ptot,eta);
					break;
				case 13:
					denomMu .Fill(ptot,eta);
					denomMup.Fill(ptot,eta);
					break;
				case 211:
					denomPi .Fill(ptot,eta);
					denomPip.Fill(ptot,eta);
					break;
				case 321:
					denomK .Fill(ptot,eta);
					denomKp.Fill(ptot,eta);
					break;
				case 2212:
					denomP .Fill(ptot,eta);
					denomPp.Fill(ptot,eta);
					break;
			}

			for(unsigned int itrk=0; itrk<trk_idx_gen->size(); ++itrk) {
				if(trk_idx_gen->at(itrk)==igen) {
					switch(pid) {
						case -2212:
							numP .Fill(ptot,eta);
							numPm.Fill(ptot,eta);
							break;
						case -321:
							numK .Fill(ptot,eta);
							numKm.Fill(ptot,eta);
							break;
						case -211:
							numPi .Fill(ptot,eta);
							numPim.Fill(ptot,eta);
							break;
						case -13:
							numMu .Fill(ptot,eta);
							numMum.Fill(ptot,eta);
							break;
						case -11:
							numE .Fill(ptot,eta);
							numEm.Fill(ptot,eta);
							break;
						case 11:
							numE .Fill(ptot,eta);
							numEp.Fill(ptot,eta);
							break;
						case 13:
							numMu .Fill(ptot,eta);
							numMup.Fill(ptot,eta);
							break;
						case 211:
							numPi .Fill(ptot,eta);
							numPip.Fill(ptot,eta);
							break;
						case 321:
							numK .Fill(ptot,eta);
							numKp.Fill(ptot,eta);
							break;
						case 2212:
							numP .Fill(ptot,eta);
							numPp.Fill(ptot,eta);
							break;
					}
					break;
				}
			}
		}


	}//end loop over entries

	effK    .Divide(&numK  , &denomK );
	effKp   .Divide(&numKp , &denomKp);
	effKm   .Divide(&numKm , &denomKm);

	effPi   .Divide(&numPi , &denomPi );
	effPip  .Divide(&numPip, &denomPip);
	effPim  .Divide(&numPim, &denomPim);

	effP    .Divide(&numP  , &denomP );
	effPp   .Divide(&numPp , &denomPp);
	effPm   .Divide(&numPm , &denomPm);

	effMu   .Divide(&numMu , &denomMu );
	effMup  .Divide(&numMup, &denomMup);
	effMum  .Divide(&numMum, &denomMum);

	effE    .Divide(&numE  , &denomE );
	effEp   .Divide(&numEp , &denomEp);
	effEm   .Divide(&numEm , &denomEm);

	TCanvas c;
	gStyle->SetOptStat(0);
	c.SetLogx();
	effK    .Draw("colz"); c.SaveAs("effK.pdf");
	effPi   .Draw("colz"); c.SaveAs("effPi.pdf");
	effP    .Draw("colz"); c.SaveAs("effP.pdf");
	effMu   .Draw("colz"); c.SaveAs("effMu.pdf");
	effE    .Draw("colz"); c.SaveAs("effE.pdf");
	effKp   .Draw("colz"); c.SaveAs("effKp.pdf");
	effPip  .Draw("colz"); c.SaveAs("effPip.pdf");
	effPp   .Draw("colz"); c.SaveAs("effPp.pdf");
	effMup  .Draw("colz"); c.SaveAs("effMup.pdf");
	effEp   .Draw("colz"); c.SaveAs("effEp.pdf");
	effKm   .Draw("colz"); c.SaveAs("effKm.pdf");
	effPim  .Draw("colz"); c.SaveAs("effPim.pdf");
	effPm   .Draw("colz"); c.SaveAs("effPm.pdf");
	effMum  .Draw("colz"); c.SaveAs("effMum.pdf");
	effEm   .Draw("colz"); c.SaveAs("effEm.pdf");

	TFile* fout = TFile::Open("trkRecoEffs_fix.root","RECREATE");

	denomK  .Write();
	denomKp .Write();
	denomKm .Write();

	denomPi .Write();
	denomPip.Write();
	denomPim.Write();

	denomP   .Write();
	denomPp  .Write();
	denomPm  .Write();

	denomMu   .Write();
	denomMup  .Write();
	denomMum  .Write();

	denomE   .Write();
	denomEp  .Write();
	denomEm  .Write();

	numK    .Write();
	numKp   .Write();
	numKm   .Write();

	numPi   .Write();
	numPip  .Write();
	numPim  .Write();

	numP   .Write();
	numPp  .Write();
	numPm  .Write();

	numMu   .Write();
	numMup  .Write();
	numMum  .Write();

	numE   .Write();
	numEp  .Write();
	numEm  .Write();

	effK    .Write();
	effKp   .Write();
	effKm   .Write();

	effPi   .Write();
	effPip  .Write();
	effPim  .Write();

	effP   .Write();
	effPp  .Write();
	effPm  .Write();

	effMu   .Write();
	effMup  .Write();
	effMum  .Write();

	effE   .Write();
	effEp  .Write();
	effEm  .Write();

	fout->Close();

}

int main(int argc, char** argv) {
	int nmax(-1);
	if(argc>1) nmax = atoi(argv[1]);

	makeRecoEff a;
	a.Loop(nmax);
	return 0;
}
