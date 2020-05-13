#ifndef DATASETMANAGER_H
#define DATASETMANAGER_H

#include <map>
#include <set>
#include <vector>

#include "TString.h"

class TChain;
class TH1;

//class to manager read-only datasets
//datasets may contain multiple files and may be random subsets or bootstraps of the data
class DatasetManager {
	public:
		//load a new dataset
		bool loadDataset(TString name, TString treeName, std::vector<TString> files,
				 unsigned long entries=0, bool bootstrap=false);
		//build a dataset from existing datasets
		bool buildDataset(TString name, TString treeName, std::vector<TString> datasets);
		//sample a dataset from an existing dataset
		bool sampleDataset(TString name, TString dataset, unsigned long entries=0, bool bootstrap=false);

		//get the TChain or list of entries for a named dataset
		//TChain* getDataset(TString name);
		//std::multiset<unsigned int>* getEntriesSet(TString name);

		//the following functions can be used to select and access a dataset
		bool setDataset(TString name, bool reset=true);
		bool setBranchAddress(TString bname, void* obj);
		TString getBranchType(TString bname);
		double draw(TString var, TString sel, TH1* h);
		unsigned long getEntries();
		double getEntries(TString sel);
		std::vector<double> getEntries(std::vector<TString> sels);
		bool getNext();
		void rewind();
		void reset();

		//get singleton instance
		static DatasetManager* getInstance() {
			if(!_instance) {
				_instance = new DatasetManager();
			}
			return _instance;
		}

	protected:
		DatasetManager() {}
	private:
		int findDataset(TString name, bool* flip=0);

		static DatasetManager* _instance;

		std::map<TString, int> _datasetLookup;
		std::vector<TChain*> _datasets;
		std::vector<std::multiset<unsigned long> > _entries;

		TChain* _currentSet;
		std::multiset<unsigned long>* _currentEntryList;
		std::multiset<unsigned long>::iterator _nextEntryIter;
		unsigned long _nextEntry;
		unsigned long _nEntries;
		bool _flipped;
};

#endif
