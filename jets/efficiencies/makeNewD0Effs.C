#define makeNewD0Effs_cxx
#include "makeNewD0Effs.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <boost/progress.hpp>

void makeNewD0Effs::Loop()
{
	if (fChain == 0) return;

	//int npt(9), neta(3);
	//double ptBins[] = {/*0.,2000.,2500.,*/3000.,3500.,4500.,5500.,7000.,9000.,11000.,15000.,20000.,100000.};
	//double etaBins[] = {2.,3.,4.,5.};

	//int npt(24-5), neta(6);
	//double ptBins[] = {/*0.,2000.,2250.,2500.,2750.,*/3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,7000.,8000.,9000.,10000.,11000.,12000.,13000.,14000.,15000.,17500.,20000.,100000.};
	//double etaBins[] = {2.,2.5,3.,3.5,4.,4.5,5.};

	int npt, neta;
	double *ptBins, *etaBins;

	TString name = sample_;

	switch(binning_) {
		case TwelveByThree:
			npt=12;
			neta=3;
			ptBins  = new double[npt +1]{3000.,4000.,5000.,6000.,7000.,8000.,9000.,10000.,12500.,15000.,20000.,30000.,100000.};
			etaBins = new double[neta+1]{2.,2.75,3.5,4.5};
			name+="_12x3bins";
			break;
		case TwelveByFive:
			npt=12;
			neta=5;
			ptBins  = new double[npt +1]{3000.,4000.,5000.,6000.,7000.,8000.,9000.,10000.,12500.,15000.,20000.,30000.,100000.};
			etaBins = new double[neta+1]{2.,2.5,3.,3.5,4.,4.5};
			name+="_12x5bins";
			break;
		case SixteenByFive:
			npt=16;
			neta=5;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,7000.,8000.,9000.,10000.,12500.,15000.,20000.,30000.,100000.};
			etaBins = new double[neta+1]{2.,2.5,3.,3.5,4.,4.5};
			name+="_16x5bins";
			break;
		case SixteenByEight:
			npt=16;
			neta=8;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,7000.,8000.,9000.,10000.,12500.,15000.,20000.,30000.,100000.};
			etaBins = new double[neta+1]{2.,2.2,2.5,2.75,3.,3.25,3.5,4.,4.5};
			name+="_16x8bins";
			break;
		case TwentyOneByFive:
			npt=21;
			neta=5;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,7000.,8000.,9000.,10000.,11000.,12000.,13000.,14000.,15000.,17500.,20000.,30000.,50000.,100000.};
			etaBins = new double[neta+1]{2.,2.5,3.,3.5,4.,4.5};
			name+="_21x5bins";
			break;
		case ThirtyFourByFive:
			npt=34;
			neta=5;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,12000.,13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,22500.,25000.,27500.,30000.,35000.,40000.,50000.,100000.};
			etaBins = new double[neta+1]{2.,2.5,3.,3.5,4.,4.5};
			name+="_34x5bins";
			break;
		case ThirtySevenByFive:
			npt=37;
			neta=5;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,12000.,13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,22500.,25000.,27500.,30000.,32500.,35000.,37500.,40000.,45000.,50000.,100000.};
			etaBins = new double[neta+1]{2.,2.5,3.,3.5,4.,4.5};
			name+="_37x5bins";
			break;
		case ThirtySevenByEight:
			npt=37;
			neta=8;
			ptBins  = new double[npt +1]{3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,12000.,13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,22500.,25000.,27500.,30000.,32500.,35000.,37500.,40000.,45000.,50000.,100000.};
			etaBins = new double[neta+1]{2.,2.2,2.5,2.75,3.,3.25,3.5,4.,4.5};
			name+="_37x8bins";
			break;
		default:
			std::cout << "Unknown binning" << std::endl;
			return;
	}

	gStyle->SetOptStat(0);

	TFile* fout = TFile::Open("D0Effs_"+name+"_up181204.root","RECREATE");
	TH2D* inaccD04 = new TH2D("inaccD04","",npt,ptBins,neta,etaBins);
	TH2D* recoD04 = new TH2D("recoD04","",npt,ptBins,neta,etaBins);
	TH2D* selD04 = new TH2D("selD04","",npt,ptBins,neta,etaBins);
	TH2D* numD04 = new TH2D("numD04","",npt,ptBins,neta,etaBins);
	TH2D* pidD04 = new TH2D("pidD04","",npt,ptBins,neta,etaBins);
	TH2D* pidmcD04 = new TH2D("pidmcD04","",npt,ptBins,neta,etaBins);
	TH2D* acceffD04 = new TH2D("acceffD04","",npt,ptBins,neta,etaBins);
	TH2D* accerrD04 = new TH2D("accerrD04","",npt,ptBins,neta,etaBins);
	TH2D* receffD04 = new TH2D("receffD04","",npt,ptBins,neta,etaBins);
	TH2D* recerrD04 = new TH2D("recerrD04","",npt,ptBins,neta,etaBins);
	TH2D* seleffD04 = new TH2D("seleffD04","",npt,ptBins,neta,etaBins);
	TH2D* selerrD04 = new TH2D("selerrD04","",npt,ptBins,neta,etaBins);
	TH2D* seleffRD04 = new TH2D("seleffRD04","",npt,ptBins,neta,etaBins);
	TH2D* selerrRD04 = new TH2D("selerrRD04","",npt,ptBins,neta,etaBins);
	TH2D* pideffD04 = new TH2D("pideffD04","",npt,ptBins,neta,etaBins);
	TH2D* piderrD04 = new TH2D("piderrD04","",npt,ptBins,neta,etaBins);
	TH2D* pidmceffD04 = new TH2D("pidmceffD04","",npt,ptBins,neta,etaBins);
	TH2D* pidmcerrD04 = new TH2D("pidmcerrD04","",npt,ptBins,neta,etaBins);
//	TH2D* diffD04 = new TH2D("diffD04","",npt,ptBins,neta,etaBins);
//	TH2D* ratioD04 = new TH2D("ratioD04","",npt,ptBins,neta,etaBins);

	TH2D* inaccD05 = new TH2D("inaccD05","",npt,ptBins,neta,etaBins);
	TH2D* recoD05 = new TH2D("recoD05","",npt,ptBins,neta,etaBins);
	TH2D* selD05 = new TH2D("selD05","",npt,ptBins,neta,etaBins);
	TH2D* numD05 = new TH2D("numD05","",npt,ptBins,neta,etaBins);
	TH2D* pidD05 = new TH2D("pidD05","",npt,ptBins,neta,etaBins);
	TH2D* pidmcD05 = new TH2D("pidmcD05","",npt,ptBins,neta,etaBins);
	TH2D* acceffD05 = new TH2D("acceffD05","",npt,ptBins,neta,etaBins);
	TH2D* accerrD05 = new TH2D("accerrD05","",npt,ptBins,neta,etaBins);
	TH2D* receffD05 = new TH2D("receffD05","",npt,ptBins,neta,etaBins);
	TH2D* recerrD05 = new TH2D("recerrD05","",npt,ptBins,neta,etaBins);
	TH2D* seleffD05 = new TH2D("seleffD05","",npt,ptBins,neta,etaBins);
	TH2D* selerrD05 = new TH2D("selerrD05","",npt,ptBins,neta,etaBins);
	TH2D* seleffRD05 = new TH2D("seleffRD05","",npt,ptBins,neta,etaBins);
	TH2D* selerrRD05 = new TH2D("selerrRD05","",npt,ptBins,neta,etaBins);
	TH2D* pideffD05 = new TH2D("pideffD05","",npt,ptBins,neta,etaBins);
	TH2D* piderrD05 = new TH2D("piderrD05","",npt,ptBins,neta,etaBins);
	TH2D* pidmceffD05 = new TH2D("pidmceffD05","",npt,ptBins,neta,etaBins);
	TH2D* pidmcerrD05 = new TH2D("pidmcerrD05","",npt,ptBins,neta,etaBins);
//	TH2D* diffD05 = new TH2D("diffD05","",npt,ptBins,neta,etaBins);
//	TH2D* ratioD05 = new TH2D("ratioD05","",npt,ptBins,neta,etaBins);

	//TH2D* numD04 = new TH2D("numD04","",npt,ptBins,neta,etaBins);
	TH2D* denomD04 = new TH2D("denomD04","",npt,ptBins,neta,etaBins);
	//TH2D* effD04 = new TH2D("effD04","",npt,ptBins,neta,etaBins);
	//TH2D* errD04 = new TH2D("errD04","",npt,ptBins,neta,etaBins);
	//TH2D* diffD04 = new TH2D("diffD04","",npt,ptBins,neta,etaBins);
	//TH2D* ratioD04 = new TH2D("ratioD04","",npt,ptBins,neta,etaBins);

	//TH2D* numD05 = new TH2D("numD05","",npt,ptBins,neta,etaBins);
	TH2D* denomD05 = new TH2D("denomD05","",npt,ptBins,neta,etaBins);
	//TH2D* effD05 = new TH2D("effD05","",npt,ptBins,neta,etaBins);
	//TH2D* errD05 = new TH2D("errD05","",npt,ptBins,neta,etaBins);
	//TH2D* diffD05 = new TH2D("diffD05","",npt,ptBins,neta,etaBins);
	//TH2D* ratioD05 = new TH2D("ratioD05","",npt,ptBins,neta,etaBins);

	TH2D* denomD045 = new TH2D("denomD045","",npt,ptBins,neta,etaBins);
	TH2D* inaccD045 = new TH2D("inaccD045","",npt,ptBins,neta,etaBins);
	TH2D* recoD045 = new TH2D("recoD045","",npt,ptBins,neta,etaBins);
	TH2D* selD045 = new TH2D("selD045","",npt,ptBins,neta,etaBins);
	TH2D* numD045 = new TH2D("numD045","",npt,ptBins,neta,etaBins);
	TH2D* pidD045 = new TH2D("pidD045","",npt,ptBins,neta,etaBins);
	TH2D* pidmcD045 = new TH2D("pidmcD045","",npt,ptBins,neta,etaBins);
	TH2D* acceffD045 = new TH2D("acceffD045","",npt,ptBins,neta,etaBins);
	TH2D* accerrD045 = new TH2D("accerrD045","",npt,ptBins,neta,etaBins);
	TH2D* receffD045 = new TH2D("receffD045","",npt,ptBins,neta,etaBins);
	TH2D* recerrD045 = new TH2D("recerrD045","",npt,ptBins,neta,etaBins);
	TH2D* seleffD045 = new TH2D("seleffD045","",npt,ptBins,neta,etaBins);
	TH2D* selerrD045 = new TH2D("selerrD045","",npt,ptBins,neta,etaBins);
	TH2D* seleffRD045 = new TH2D("seleffRD045","",npt,ptBins,neta,etaBins);
	TH2D* selerrRD045 = new TH2D("selerrRD045","",npt,ptBins,neta,etaBins);
	TH2D* pideffD045 = new TH2D("pideffD045","",npt,ptBins,neta,etaBins);
	TH2D* piderrD045 = new TH2D("piderrD045","",npt,ptBins,neta,etaBins);
	TH2D* pidmceffD045 = new TH2D("pidmceffD045","",npt,ptBins,neta,etaBins);
	TH2D* pidmcerrD045 = new TH2D("pidmcerrD045","",npt,ptBins,neta,etaBins);

	inaccD04->Sumw2();
	recoD04->Sumw2();
	selD04->Sumw2();
	numD04->Sumw2();
	pidD04->Sumw2();
	pidmcD04->Sumw2();
	acceffD04->Sumw2();
	accerrD04->Sumw2();
	receffD04->Sumw2();
	recerrD04->Sumw2();
	seleffD04->Sumw2();
	selerrD04->Sumw2();
	seleffRD04->Sumw2();
	selerrRD04->Sumw2();
	pideffD04->Sumw2();
	piderrD04->Sumw2();
	pidmceffD04->Sumw2();
	pidmcerrD04->Sumw2();
//	diffD04->Sumw2();
//	ratioD04->Sumw2();

	inaccD05->Sumw2();
	recoD05->Sumw2();
	selD05->Sumw2();
	numD05->Sumw2();
	pidD05->Sumw2();
	pidmcD05->Sumw2();
	acceffD05->Sumw2();
	accerrD05->Sumw2();
	receffD05->Sumw2();
	recerrD05->Sumw2();
	seleffD05->Sumw2();
	selerrD05->Sumw2();
	seleffRD05->Sumw2();
	selerrRD05->Sumw2();
	pideffD05->Sumw2();
	piderrD05->Sumw2();
	pidmceffD05->Sumw2();
	pidmcerrD05->Sumw2();
//	diffD05->Sumw2();
//	ratioD05->Sumw2();

	//numD04->Sumw2();
	denomD04->Sumw2();
	//effD04->Sumw2();
	//errD04->Sumw2();
	//diffD04->Sumw2();
	//ratioD04->Sumw2();

	//numD05->Sumw2();
	denomD05->Sumw2();
	//effD05->Sumw2();
	//errD05->Sumw2();
	//diffD05->Sumw2();
	//ratioD05->Sumw2();

	denomD045->Sumw2();
	inaccD045->Sumw2();
	recoD045->Sumw2();
	selD045->Sumw2();
	numD045->Sumw2();
	pidD045->Sumw2();
	pidmcD045->Sumw2();
	acceffD045->Sumw2();
	accerrD045->Sumw2();
	receffD045->Sumw2();
	recerrD045->Sumw2();
	seleffD045->Sumw2();
	selerrD045->Sumw2();
	seleffRD045->Sumw2();
	selerrRD045->Sumw2();
	pideffD045->Sumw2();
	piderrD045->Sumw2();
	pidmceffD045->Sumw2();
	pidmcerrD045->Sumw2();

	Long64_t nentries = fChain->GetEntries();

	Long64_t nbytes = 0, nb = 0;
	boost::progress_display progress( nentries );
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		++progress;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		for(unsigned int d=0; d<TRUEDID->size(); ++d) {
			if(TMath::Abs(TRUEDID->at(d))!=421) continue; //not D0
			//fixed - previously allowed through other D0 decays (where the children are not saved)
			if(TRUEDTRK1IDX->at(d)==-1 || TRUEDTRK2IDX->at(d)!=-1) continue; //not 2-body
			bool fromb = TRUEDFROMB->at(d);
			TLorentzVector p(TRUEDPX->at(d),TRUEDPY->at(d),TRUEDPZ->at(d),TRUEDE->at(d));

			//denominator
			if(fromb) denomD05->Fill(p.Pt(),p.Eta());
			else denomD04->Fill(p.Pt(),p.Eta());
			denomD045->Fill(p.Pt(),p.Eta());

			//tracks in acceptance (DecProdCut: ~1.6<eta<5.3)
			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
			//tighter cuts now included in inacc
			////also apply tighter acceptance cuts
			//if(TRUEDTRK0P->at(d) <5. || TRUEDTRK1P->at(d) <5.) continue;
			//if(TRUEDTRK0PT->at(d)<.5 || TRUEDTRK1PT->at(d)<.5) continue;
			//double trk0PL = TMath::Sqrt(TRUEDTRK0P->at(d)*TRUEDTRK0P->at(d)-TRUEDTRK0PT->at(d)*TRUEDTRK0PT->at(d));
			//double trk1PL = TMath::Sqrt(TRUEDTRK1P->at(d)*TRUEDTRK1P->at(d)-TRUEDTRK1PT->at(d)*TRUEDTRK1PT->at(d));
			//double trk0ETA = 0.5*TMath::Log((TRUEDTRK0P->at(d)+trk0PL)/(TRUEDTRK0P->at(d)-trk0PL));
			//double trk1ETA = 0.5*TMath::Log((TRUEDTRK1P->at(d)+trk1PL)/(TRUEDTRK1P->at(d)-trk1PL));
			//if(trk0ETA<2. || trk0ETA>4.5 || trk1ETA<2. || trk1ETA>4.5) continue;
			if(fromb) inaccD05->Fill(p.Pt(),p.Eta());
			else inaccD04->Fill(p.Pt(),p.Eta());
			inaccD045->Fill(p.Pt(),p.Eta());
			
			//tracks reconstructed
			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;
			if(fromb) recoD05->Fill(p.Pt(),p.Eta());
			else recoD04->Fill(p.Pt(),p.Eta());
			recoD045->Fill(p.Pt(),p.Eta());

			//for(unsigned int s=0; s<D0TRUEIDX->size(); ++s) {
			//	if(D0TRUEIDX->at(s)==d) {
			//		if(fromb) numD05->Fill(p.Pt(),p.Eta());
			//		else numD04->Fill(p.Pt(),p.Eta());
			//		break;
			//	}
			//}

			//D0 selected
			for(unsigned int s=0; s<D0TRUEIDX->size(); ++s) {
				if(D0TRUEIDX->at(s)==d) {
					//get reconstructed 4-momentum to include bin migration
					TVector3 prec  (D0PX->at(s),  D0PY->at(s),  D0PZ->at(s));
					TVector3 precK (D0KPX->at(s), D0KPY->at(s), D0KPZ->at(s));
					TVector3 precPi(D0PIPX->at(s),D0PIPY->at(s),D0PIPZ->at(s));

					double pK  = precK.Mag();
					//double pPi = precPi.Mag();
					double ptK = precK.Pt();
					//double ptPi= precPi.Pt();

					//if(pK  >500000.) pK  =499000.;
					//if(pPi >500000.) pPi =499000.;
					//if(ptK > 50000.) ptK = 49000.;
					//if(ptPi> 25000.) ptPi= 24000.;

					double effK = pidK->GetBinContent(pidK->FindBin(pK,ptK));
					double effPi = 1.;//pidPi->GetBinContent(pidPi->FindBin(pPi,ptPi)); //pion PID turned off

					if(pK>500000. || ptK>25000.) effK=1.; //no PID cut at high P or PT

					if(fromb) {
						numD05->Fill(p.Pt(),p.Eta());
						pidD05->Fill(p.Pt(),p.Eta(),effK*effPi);
						selD05->Fill(prec.Pt(),prec.Eta());
						if(D0KPNNK->at(s)>0.2 || pK>500000. || ptK>25000.) //&&D0PIPNNPI->at(s)>0.1)
							pidmcD05->Fill(p.Pt(),p.Eta());
					} else {
						numD04->Fill(p.Pt(),p.Eta());
						pidD04->Fill(p.Pt(),p.Eta(),effK*effPi);
						selD04->Fill(prec.Pt(),prec.Eta());
						if(D0KPNNK->at(s)>0.2 || pK>500000. || ptK>25000.) //&&D0PIPNNPI->at(s)>0.1)
							pidmcD04->Fill(p.Pt(),p.Eta());
					}
					numD045->Fill(p.Pt(),p.Eta());
					pidD045->Fill(p.Pt(),p.Eta(),effK*effPi);
					selD045->Fill(prec.Pt(),prec.Eta());
						if(D0KPNNK->at(s)>0.2 || pK>500000. || ptK>25000.) //&&D0PIPNNPI->at(s)>0.1)
							pidmcD045->Fill(p.Pt(),p.Eta());
					break;
				}
			}
		}
	}

	acceffD04->Divide(inaccD04,denomD04);
	acceffD05->Divide(inaccD05,denomD05);
	acceffD045->Divide(inaccD045,denomD045);
	receffD04->Divide(recoD04,inaccD04);
	receffD05->Divide(recoD05,inaccD05);
	receffD045->Divide(recoD045,inaccD045);
	seleffD04->Divide(numD04,recoD04);
	seleffD05->Divide(numD05,recoD05);
	seleffD045->Divide(numD045,recoD045);
	seleffRD04->Divide(selD04,recoD04);
	seleffRD05->Divide(selD05,recoD05);
	seleffRD045->Divide(selD045,recoD045);
	pideffD04->Divide(pidD04,numD04);
	pideffD05->Divide(pidD05,numD05);
	pideffD045->Divide(pidD045,numD045);
	pidmceffD04->Divide(pidmcD04,numD04);
	pidmceffD05->Divide(pidmcD05,numD05);
	pidmceffD045->Divide(pidmcD045,numD045);

	//effD04->Divide(numD04,denomD04);
	//effD05->Divide(numD05,denomD05);

	for(int i=1; i<=npt; ++i) {
		for(int j=1; j<=neta; ++j) {
			accerrD04->SetBinContent(i,j,TMath::Sqrt(inaccD04->GetBinContent(i,j))/  denomD04->GetBinContent(i,j));
			accerrD05->SetBinContent(i,j,TMath::Sqrt(inaccD05->GetBinContent(i,j))/  denomD05->GetBinContent(i,j));
			accerrD045->SetBinContent(i,j,TMath::Sqrt(inaccD045->GetBinContent(i,j))/  denomD045->GetBinContent(i,j));
			recerrD04->SetBinContent(i,j,TMath::Sqrt( recoD04->GetBinContent(i,j))/inaccD04->GetBinContent(i,j));
			recerrD05->SetBinContent(i,j,TMath::Sqrt( recoD05->GetBinContent(i,j))/inaccD05->GetBinContent(i,j));
			recerrD045->SetBinContent(i,j,TMath::Sqrt( recoD045->GetBinContent(i,j))/inaccD045->GetBinContent(i,j));
			selerrD04->SetBinContent(i,j,TMath::Sqrt(  numD04->GetBinContent(i,j))/ recoD04->GetBinContent(i,j));
			selerrD05->SetBinContent(i,j,TMath::Sqrt(  numD05->GetBinContent(i,j))/ recoD05->GetBinContent(i,j));
			selerrD045->SetBinContent(i,j,TMath::Sqrt(  numD045->GetBinContent(i,j))/ recoD045->GetBinContent(i,j));
			selerrRD04->SetBinContent(i,j,TMath::Sqrt(  selD04->GetBinContent(i,j))/ recoD04->GetBinContent(i,j));
			selerrRD05->SetBinContent(i,j,TMath::Sqrt(  selD05->GetBinContent(i,j))/ recoD05->GetBinContent(i,j));
			selerrRD045->SetBinContent(i,j,TMath::Sqrt(  selD045->GetBinContent(i,j))/ recoD045->GetBinContent(i,j));
			//TODO piderrD04->SetBinContent(i,j,TMath::Sqrt(  numD04->GetBinContent(i,j))/ recoD04->GetBinContent(i,j));
			//TODO piderrD05->SetBinContent(i,j,TMath::Sqrt(  numD05->GetBinContent(i,j))/ recoD05->GetBinContent(i,j));
			//TODO piderrD045->SetBinContent(i,j,TMath::Sqrt(  numD045->GetBinContent(i,j))/ recoD045->GetBinContent(i,j));
			//errD04->SetBinContent(i,j,TMath::Sqrt(numD04->GetBinContent(i,j))/denomD04->GetBinContent(i,j));
			//errD05->SetBinContent(i,j,TMath::Sqrt(numD05->GetBinContent(i,j))/denomD05->GetBinContent(i,j));

//			//differences from neighbouring bins
//			double diffS4(0.), diffS5(0.);//, diffD4(0.), diffD5(0.);
//			int nTerms(0);
//			if(i>1) {
//				diffS4+= TMath::Abs(effD04->GetBinContent(i,j) - effD04->GetBinContent(i-1,j));
//				diffS5+= TMath::Abs(effD05->GetBinContent(i,j) - effD05->GetBinContent(i-1,j));
//				//diffD4+= TMath::Abs(effD04->GetBinContent(i,j)   - effD04->GetBinContent(i-1,j));
//				//diffD5+= TMath::Abs(effD05->GetBinContent(i,j)   - effD05->GetBinContent(i-1,j));
//				++nTerms;
//			}
//			if(i<npt) {
//				diffS4+= TMath::Abs(effD04->GetBinContent(i,j) - effD04->GetBinContent(i+1,j));
//				diffS5+= TMath::Abs(effD05->GetBinContent(i,j) - effD05->GetBinContent(i+1,j));
//				//diffD4+= TMath::Abs(effD04->GetBinContent(i,j)   - effD04->GetBinContent(i+1,j));
//				//diffD5+= TMath::Abs(effD05->GetBinContent(i,j)   - effD05->GetBinContent(i+1,j));
//				++nTerms;
//			}
//			if(j>1) {
//				diffS4+= TMath::Abs(effD04->GetBinContent(i,j) - effD04->GetBinContent(i,j-1));
//				diffS5+= TMath::Abs(effD05->GetBinContent(i,j) - effD05->GetBinContent(i,j-1));
//				//diffD4+= TMath::Abs(effD04->GetBinContent(i,j)   - effD04->GetBinContent(i,j-1));
//				//diffD5+= TMath::Abs(effD05->GetBinContent(i,j)   - effD05->GetBinContent(i,j-1));
//				++nTerms;
//			}
//			if(j<neta) {
//				diffS4+= TMath::Abs(effD04->GetBinContent(i,j) - effD04->GetBinContent(i,j+1));
//				diffS5+= TMath::Abs(effD05->GetBinContent(i,j) - effD05->GetBinContent(i,j+1));
//				//diffD4+= TMath::Abs(effD04->GetBinContent(i,j)   - effD04->GetBinContent(i,j+1));
//				//diffD5+= TMath::Abs(effD05->GetBinContent(i,j)   - effD05->GetBinContent(i,j+1));
//				++nTerms;
//			}
//
//			diffS4/=nTerms;
//			diffS5/=nTerms;
//			//diffD4/=nTerms;
//			//diffD5/=nTerms;
//
//			diffD04->SetBinContent(i,j,diffS4);
//			diffD05->SetBinContent(i,j,diffS5);
//			//diffD04->SetBinContent(i,j,diffD4);
//			//diffD05->SetBinContent(i,j,diffD5);

		}
	}

//	ratioD04->Divide(diffD04,errD04);
//	ratioD05->Divide(diffD05,errD05);
//	//ratioD04->Divide(diffD04,errD04);
//	//ratioD05->Divide(diffD05,errD05);


	inaccD04->Write();
	recoD04->Write();
	selD04->Write();
	numD04->Write();
	pidD04->Write();
	pidmcD04->Write();
	//numD04->Write();
	denomD04->Write();
	acceffD04->Write();
	receffD04->Write();
	seleffD04->Write();
	seleffRD04->Write();
	pideffD04->Write();
	pidmceffD04->Write();
	//effD04->Write();
	accerrD04->Write();
	recerrD04->Write();
	selerrD04->Write();
	selerrRD04->Write();
	piderrD04->Write();
	pidmcerrD04->Write();
	//errD04->Write();
//	diffD04->Write();
	//diffD04->Write();
//	ratioD04->Write();
	//ratioD04->Write();
	
	inaccD05->Write();
	recoD05->Write();
	selD05->Write();
	numD05->Write();
	pidD05->Write();
	pidmcD05->Write();
	//numD05->Write();
	denomD05->Write();
	acceffD05->Write();
	receffD05->Write();
	seleffD05->Write();
	seleffRD05->Write();
	pideffD05->Write();
	pidmceffD05->Write();
	//effD05->Write();
	accerrD05->Write();
	recerrD05->Write();
	selerrD05->Write();
	selerrRD05->Write();
	piderrD05->Write();
	pidmcerrD05->Write();
	//errD05->Write();
//	diffD05->Write();
	//diffD05->Write();
//	ratioD05->Write();
	//ratioD05->Write();
	
	inaccD045->Write();
	recoD045->Write();
	selD045->Write();
	pidD045->Write();
	pidmcD045->Write();
	numD045->Write();
	denomD045->Write();
	acceffD045->Write();
	receffD045->Write();
	seleffD045->Write();
	seleffRD045->Write();
	pideffD045->Write();
	pidmceffD045->Write();
	accerrD045->Write();
	recerrD045->Write();
	selerrD045->Write();
	selerrRD045->Write();
	piderrD045->Write();
	pidmcerrD045->Write();
	fout->Close();
}

int main(int argc, char** argv) {
	TString sample="2XX";
	if(argc>1) sample=argv[1];
	makeNewD0Effs::Binning binning=makeNewD0Effs::SixteenByEight;
	if(argc>2) {
		TString binStr = argv[2];
		if(binStr=="12x3") {
			binning = makeNewD0Effs::TwelveByThree;
		} else if(binStr=="12x5") {
			binning = makeNewD0Effs::TwelveByFive;
		} else if(binStr=="16x5") {
			binning = makeNewD0Effs::SixteenByFive;
		} else if(binStr=="16x8") {
			binning = makeNewD0Effs::SixteenByEight;
		} else if(binStr=="21x5") {
			binning = makeNewD0Effs::TwentyOneByFive;
		} else if(binStr=="34x5") {
			binning = makeNewD0Effs::ThirtyFourByFive;
		} else if(binStr=="37x5") {
			binning = makeNewD0Effs::ThirtySevenByFive;
		} else if(binStr=="37x8") {
			binning = makeNewD0Effs::ThirtySevenByEight;
		} else {
			std::cout << "Unknown binning" << std::endl;
			return 1;
		}
	}
	makeNewD0Effs a(sample, binning);
	a.Loop();
	return 0;
}
