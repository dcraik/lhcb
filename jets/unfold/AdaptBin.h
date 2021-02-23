// Class to produce a 2D histogram with equally populated bins based on some input data
// Author : dcraik

#include <iostream>
#include <list>
#include <vector>

#include "TH2Poly.h"
#include "TString.h"

class AdaptBin {
	public:
		AdaptBin(TString name, int nbin, double xmin, double xmax, double ymin, double ymax)
			: name_(name), xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
			xs_(0), ys_(0), nentries_(0),
			loaded_(false), usingLists_(false), verbose_(false),
			theHisto_(0) {calcDivisions(nbin);}

		~AdaptBin() {
			delete[] xs_;
			delete[] ys_;
			if(theHisto_) delete theHisto_;
		}

		bool loadDataFromTree(TString fname, TString tname, TString xname, TString yname);
		bool loadDataFromHist(TString fname, TString hname);
		bool loadDataFromHist(TH2* hist);

		bool addEntry(double x, double y);

		TH2Poly* getHisto(TString name="");

		void setVerbose(bool verbose=true) {verbose_ = verbose;}

	private:
		void calcDivisions(int ntarget);
		bool loadDataFromLists();
		void initHisto(double xmin, double xmax, double ymin, double ymax, uint iter=0);

		TString name_;
		double xmin_, xmax_, ymin_, ymax_;

		std::list<double> xList_, yList_;
		double *xs_, *ys_;
		std::vector<int> divisions_;
		int nentries_;
		bool loaded_, usingLists_, verbose_;

		TH2Poly* theHisto_;
};
