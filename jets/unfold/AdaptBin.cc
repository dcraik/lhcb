// Class to produce a 2D histogram with equally populated bins based on some input data
// Author : dcraik

#include "AdaptBin.h"

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

bool AdaptBin::loadDataFromTree(TString fname, TString tname, TString xname, TString yname) {
	if(loaded_) {
		std::cout << "ERROR in AdaptBin::loadDataFromTree : Histogram already loaded." << std::endl;
		return false;
	}
	if(usingLists_) {
		std::cout << "ERROR in AdaptBin::loadDataFromTree : Entries already added to lists." << std::endl;
		std::cout << "                                      Unable to also load from tree." << std::endl;
		return false;
	}

	TFile* file = TFile::Open(fname);
	if(!file) {
		std::cout << "ERROR in AdaptBin::loadDataFromTree : file " << fname << " does not exist." << std::endl;
		return false;
	}

	TTree* tree = dynamic_cast<TTree*>(file->Get(tname));
	if(!tree) {
		std::cout << "ERROR in AdaptBin::loadDataFromTree : tree " << tname << " does not exist in file " << fname << "." << std::endl;
		return false;
	}

	double x;
	double y;

	tree->SetBranchAddress(xname.Data(),&x);
	tree->SetBranchAddress(yname.Data(),&y);

	nentries_ = tree->GetEntries();

	xs_ = new double[nentries_];
	ys_ = new double[nentries_];

	for ( int i=0; i < nentries_; ++i ) {
		tree->GetEntry( i );
		xs_[i] = x;
		ys_[i] = y;
	}

	file->Close();
	loaded_=true;

	if(verbose_) std::cout << "INFO in AdaptBin::loadDataFromTree : Histogram successfully loaded from tree." << std::endl;
	return true;
}

bool AdaptBin::loadDataFromHist(TString fname, TString hname) {
	if(loaded_) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : Histogram already loaded." << std::endl;
		return false;
	}
	if(usingLists_) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : Entries already added to lists." << std::endl;
		std::cout << "                                      Unable to also load from histogram." << std::endl;
		return false;
	}

	TFile* file = TFile::Open(fname);
	if(!file) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : file " << fname << " does not exist." << std::endl;
		return false;
	}

	TH2* hist = dynamic_cast<TH2*>(file->Get(hname));
	if(!hist) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : histogram " << hname << " does not exist in file " << fname << "." << std::endl;
		return false;
	}

	bool ret = loadDataFromHist(hist);
	file->Close();

	return ret;
}

bool AdaptBin::loadDataFromHist(TH2* hist) {
	if(loaded_) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : Histogram already loaded." << std::endl;
		return false;
	}
	if(usingLists_) {
		std::cout << "ERROR in AdaptBin::loadDataFromHist : Entries already added to lists." << std::endl;
		std::cout << "                                      Unable to also load from histogram." << std::endl;
		return false;
	}

	nentries_ = hist->GetSumOfWeights();

	xs_ = new double[nentries_];
	ys_ = new double[nentries_];

	int ientry(0);
	int skipped(0);
	for ( int i=1; i <= hist->GetNbinsX(); ++i ) {
		for ( int j=1; j < hist->GetNbinsY(); ++j ) {
			double x = hist->GetXaxis()->GetBinCenter(i);
			double y = hist->GetYaxis()->GetBinCenter(j);
			int n = TMath::Nint(hist->GetBinContent( i, j ));
			for ( int k=0; k < n; ++k ) {
				if(ientry==nentries_) {
					++skipped;
					break;
				}
				xs_[ientry] = x;
				ys_[ientry] = y;
				++ientry;
			}
		}
	}
	if(skipped>0) {
		std::cout << "WARNING in AdaptBin::loadDataFromHist : input histogram had too many entries" << std::endl;
		std::cout << "                                        this might happen if the bin contents were non-integer so we had to round them" << std::endl;
		std::cout << "                                        last " << skipped << " of " << nentries_+skipped << " were ignored" << std::endl;
	}

	loaded_=true;

	if(verbose_) std::cout << "INFO in AdaptBin::loadDataFromHist : Histogram successfully loaded from histogram." << std::endl;
	return true;
}

bool AdaptBin::loadDataFromLists() {
	if(loaded_) {
		std::cout << "ERROR in AdaptBin::loadDataFromLists : Histogram already loaded." << std::endl;
		return false;
	}
	if(!usingLists_) {
		std::cout << "ERROR in AdaptBin::loadDataFromLists : No entries in lists." << std::endl;
		std::cout << "                                       Unable to also load from lists." << std::endl;
		return false;
	}
	if(xList_.size()!=yList_.size()) {
		std::cout << "ERROR in AdaptBin::loadDataFromLists : Different numbers of entries in lists." << std::endl;
		return false;
	}
	nentries_ = xList_.size();

	xs_ = new double[nentries_];
	ys_ = new double[nentries_];

	std::list<double>::iterator itX, itY;
	itX=xList_.begin();
	itY=yList_.begin();
	for(int i=0; i<nentries_; ++i) {
		if(itX==xList_.end() || itY==yList_.end()) {
			std::cout << "ERROR in AdaptBin::loadDataFromLists : Reached end of list sooner than expected." << std::endl;
			return false;
		}
		xs_[i] = (*itX);
		ys_[i] = (*itY);
		++itX;
		++itY;
	}

	xList_.clear();
	yList_.clear();
	loaded_=true;

	if(verbose_) std::cout << "INFO in AdaptBin::loadDataFromLists : Histogram successfully loaded from lists." << std::endl;
	return true;
}

bool AdaptBin::addEntry(double x, double y) {
	if(loaded_) {
		std::cout << "ERROR in AdaptBin::addEntry : Histogram already loaded." << std::endl;
		std::cout << "                              Unable to add further entries to lists." << std::endl;
		return false;
	}
	usingLists_=true;
	if(x<xmin_ || x>xmax_ || y<ymin_ || y>ymax_) return false;//ignore out of range entries

	xList_.push_back(x);
	yList_.push_back(y);

	return true;
}

TH2Poly* AdaptBin::getHisto(TString name) {
	if(theHisto_) {
		if(name) {
			return static_cast<TH2Poly*>(theHisto_->Clone(name));
		} else {
			return theHisto_;
		}
	}

	if(usingLists_) {
		if(!loadDataFromLists()) {
			std::cout << "ERROR in AdaptBin::getHisto : Failed to load histogram from lists." << std::endl;
			return 0;
		}
	}

	if(!loaded_) {
		std::cout << "ERROR in AdaptBin::getHisto : Histogram not loaded." << std::endl;
		std::cout << "                              Either load from a tree or add entries to the lists to load from lists." << std::endl;
		return 0;
	}

	theHisto_ = new TH2Poly(name_, "", xmin_, xmax_, ymin_, ymax_);
	initHisto(xmin_,xmax_,ymin_,ymax_,0);

	if(verbose_) std::cout << "INFO in AdaptBin::getHisto : Histogram successfully initialised." << std::endl;
	if(name) {
		return static_cast<TH2Poly*>(theHisto_->Clone(name));
	} else {
		return theHisto_;
	}
}

void AdaptBin::calcDivisions(int ntarget) {
	std::cout << "INFO in AdaptBin::calcDivisions : Target is " << ntarget << " bins." << std::endl;

	//we will iteratively sub-divide histogram bins into either 4, 9, 25, 49 or 121
	//here we figure out how many 4s, 9s, 25s, 49s and 121s to use to best match our target without exceeding it
	int nDivisions(0), nNines(0), nTwentyFives(0), nFortyNines(0), nEleventyElevens(0), nBins(1);
	int nDivisionsBest(0), nNinesBest(0), nTwentyFivesBest(0), nFortyNinesBest(0), nEleventyElevensBest(0), nBinsBest(1);

	do {
		++nDivisions;
		for(nNines=0; nNines<=nDivisions; ++nNines) {
			for(nTwentyFives=0; nTwentyFives<=nDivisions-nNines; ++nTwentyFives) {
				for(nFortyNines=0; nFortyNines<=nDivisions-nNines-nTwentyFives; ++nFortyNines) {
					for(nEleventyElevens=0; nEleventyElevens<=nDivisions-nNines-nTwentyFives-nFortyNines; ++nEleventyElevens) {
						nBins =  TMath::Power(4,nDivisions-nNines-nTwentyFives-nFortyNines-nEleventyElevens)
							*TMath::Power(9,nNines)*TMath::Power(25,nTwentyFives)
							*TMath::Power(49,nFortyNines)*TMath::Power(121,nEleventyElevens);
						if(nBins <= ntarget && nBins > nBinsBest) {
							//keep track of the best number of bins so far
							nBinsBest = nBins;
							nDivisionsBest = nDivisions;
							nNinesBest = nNines;
							nTwentyFivesBest = nTwentyFives;
							nFortyNinesBest = nFortyNines;
							nEleventyElevensBest = nEleventyElevens;
						}
					}
				}
			}
		}
	} while(TMath::Power(4,nDivisions+1) <= ntarget);//if 4^n > target then we've gone far enough

	std::cout << "INFO in AdaptBin::calcDivisions : Using " << nBinsBest << " bins." << std::endl;

	//fill the vector with the divisions that we want to make
	for(int i=0; i<nEleventyElevensBest; ++i) {
		divisions_.push_back(11);
	}
	for(int i=0; i<nFortyNinesBest; ++i) {
		divisions_.push_back(7);
	}
	for(int i=0; i<nTwentyFivesBest; ++i) {
		divisions_.push_back(5);
	}
	for(int i=0; i<nNinesBest; ++i) {
		divisions_.push_back(3);
	}
	for(int i=0; i<nDivisionsBest-nNinesBest-nTwentyFivesBest-nFortyNinesBest-nEleventyElevensBest; ++i) {
		divisions_.push_back(2);
	}
}

void AdaptBin::initHisto(double xmin, double xmax, double ymin, double ymax, uint iter) {

	//If it's the last iteration create the bin and return
	if(iter == divisions_.size()) {
		double * x_new = new double[5];
		double * y_new = new double[5];
		x_new[0] = xmin; x_new[1] = xmin; x_new[2] = xmax; x_new[3] = xmax; x_new[4] = xmin;
		y_new[0] = ymin; y_new[1] = ymax; y_new[2] = ymax; y_new[3] = ymin; y_new[4] = ymin;
		theHisto_->AddBin(5, x_new, y_new);
		if(verbose_) std::cout << "INFO in AdaptBin::initHisto : Adding bin from (" << xmin << "," << ymin << ") to (" << xmax << "," << ymax << ")." << std::endl;
		return;
	}

	//If not the last iteration then divide the bin
	int n_divx=divisions_[iter];
	int n_divy=divisions_[iter];

	if(verbose_) std::cout << "INFO in AdaptBin::initHisto : Dividing bin from (" << xmin << "," << ymin << ") to (" << xmax << "," << ymax << ") into " << n_divx << " by " << n_divy << " subbins" << std::endl;

	double *xIn = new double[nentries_]; 
	double *yIn = new double[nentries_]; 
	int *xIndex = new int [nentries_+2];
	int *yIndex = new int [nentries_+2];

	int xCountIn = 0; 
	for(int i = 0; i<nentries_; ++i) {
		if ((xs_[i]<xmin)||(xs_[i]>xmax)||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
		xIn[xCountIn] = xs_[i];
		++xCountIn;
	}

	//find the delimitting x and y values for the sub bins
	double xLimits[n_divx + 1];
	double yLimits[n_divx][n_divy + 1];

	//first sort entries in x and divide bin into equally populated bins in x
	TMath::Sort(xCountIn, xIn, xIndex, false);

	xLimits[0] = xmin;
	xLimits[n_divx] = xmax;
	for (int nDivx = 0; nDivx < n_divx; ++nDivx){
		if (nDivx  < (n_divx-1)){
			if(xCountIn>0) {
				xLimits[nDivx+1] = xIn[xIndex[xCountIn*(nDivx+1)/n_divx]];
			} else {
				//if no entries then use equal bin widths
				xLimits[nDivx+1] = xmin + (xmax-xmin)*(nDivx+1)/n_divx;
			}
		}

		//for each bin in x divide into equally populated bins in y
		yLimits[nDivx][0] = ymin;
		yLimits[nDivx][n_divy] = ymax;    
		int yCountIn = 0;

		for(int i = 0; i<nentries_; ++i) {
			if ((xs_[i]<xmin)||(xs_[i]>xmax)||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
			if ((xs_[i]<xLimits[nDivx])||(xs_[i]>=xLimits[nDivx+1])||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
			yIn[yCountIn] = ys_[i];
			++yCountIn;
		}

		TMath::Sort(yCountIn, yIn, yIndex, false);

		for (int nDivy = 1; nDivy < n_divy; ++nDivy){
			if(yCountIn>0) {
				yLimits[nDivx][nDivy] = yIn[yIndex[yCountIn*nDivy/n_divy]];
			} else {
				//if no entries then use equal bin widths
				yLimits[nDivx][nDivy] = ymin + (ymax-ymin)*nDivy/n_divy;
			}
		}
	}

	delete[] xIn;
	delete[] yIn;
	delete[] xIndex;
	delete[] yIndex;

	//call for each sub bin
	for (int nDivx = 0; nDivx < n_divx; ++nDivx){
		for (int nDivy = 0; nDivy < n_divy; ++nDivy){
			initHisto(xLimits[nDivx], xLimits[nDivx + 1], yLimits[nDivx][nDivy], yLimits[nDivx][nDivy + 1],iter+1);
		}
	}
}

