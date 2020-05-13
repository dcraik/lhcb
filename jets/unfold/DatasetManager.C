#include "DatasetManager.h"

#include <iostream>

//#include "TCanvas.h"//TODO
#include "TChain.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTreeFormula.h"

DatasetManager* DatasetManager::_instance=0;

bool DatasetManager::loadDataset(TString name, TString treeName, std::vector<TString> files,
		                 unsigned long entries, bool bootstrap) {
	std::cout << "INFO in DatasetManager::loadDataset : loading " << name << " from " << files.size() << " files..." << std::endl;
	TChain* c = new TChain(treeName);
	for(std::vector<TString>::iterator it=files.begin(); it!=files.end(); ++it) {
		c->AddFile(*it);
	}
	unsigned long maxEntries = c->GetEntries();

	std::multiset<unsigned long> s;

	//need to construct a set of entries to use
	if(entries>0) {
		if(bootstrap) {
			//allow repeats
			while(s.size()<entries) {
				s.insert(gRandom->Integer(maxEntries));
			}
		} else {
			//no repeats - use set and copy to multiset
			if(entries<maxEntries) {
				std::set<unsigned long> s2;
				while(s2.size()<entries) {
					s2.insert(gRandom->Integer(maxEntries));
				}
				s = std::multiset<unsigned long>(s2.begin(), s2.end());
			}
		}
	}

	_datasetLookup[name] = _datasets.size();
	_datasets.push_back(c);
	_entries.push_back(s);
	std::cout << "INFO in DatasetManager::loadDataset : " << name << " loaded successfully (" << (s.size()==0?maxEntries:s.size()) << " events)" << std::endl;

	return true;
}

bool DatasetManager::buildDataset(TString name, TString treeName, std::vector<TString> datasets) {
	std::cout << "INFO in DatasetManager::buildDataset : building " << name << " from " << datasets.size() << " existing datasets..." << std::endl;
	//first check all datasets exist
	for(unsigned int i=0; i<datasets.size(); ++i) {
		if(findDataset(datasets[i])<0) return false;
	}

	TChain* c = new TChain(treeName);
	std::multiset<unsigned long> s;

	for(std::vector<TString>::iterator it=datasets.begin(); it!=datasets.end(); ++it) {
		bool flip(false);
		int idx = findDataset(*it,&flip);
		TChain* ds = _datasets[idx];
		std::multiset<unsigned long> e = _entries[idx];
		unsigned long offset = c->GetEntries();
		//TODO//std::cout << s.size() << " " << offset << std::endl;//TODO
		c->Add(ds);
		if(flip) {
			for(int evt=0; evt<ds->GetEntries(); ++evt) {
				if(e.find(evt)==e.end()) {
					s.insert(s.end(),evt+offset);
				}
			}
		} else {
			for(std::multiset<unsigned long>::iterator it2=e.begin(); it2!=e.end(); ++it2) {
				auto r = s.insert(s.end(),(*it2)+offset);
				//TODO//if(offset!=0) std::cout << offset << " " << (*it2) << " " << offset+(*it2) << " " << e.size() << " " << *(r) << std::endl;//TODO
			}
		}
	}

	_datasetLookup[name] = _datasets.size();
	_datasets.push_back(c);
	_entries.push_back(s);
	std::cout << "INFO in DatasetManager::buildDataset : " << name << " built successfully (" << s.size() << " events)" << std::endl;

	return true;
}

//construct a sampled dataset using the TChain in an existing dataset
//note that this ignores any entrylist in the sampled dataset
//TODO need some sort of iterative accept/reject to allow pT sculpting
bool DatasetManager::sampleDataset(TString name, TString dataset, unsigned long entries, bool bootstrap) {
	std::cout << "INFO in DatasetManager::sampleDataset : loading " << name << " from " << dataset << "..." << std::endl;

	int input = findDataset(dataset);
	if(input<0) return false;

	TChain* c = _datasets[input];
	unsigned long maxEntries = c->GetEntries();

	std::multiset<unsigned long> s;

	//need to construct a set of entries to use
	if(entries>0) {
		if(bootstrap) {
			//allow repeats
			while(s.size()<entries) {
				s.insert(gRandom->Integer(maxEntries));
			}
		} else {
			//no repeats - use set and copy to multiset
			if(entries<maxEntries) {
				std::set<unsigned long> s2;
				while(s2.size()<entries) {
					s2.insert(gRandom->Integer(maxEntries));
				}
				s = std::multiset<unsigned long>(s2.begin(), s2.end());
			}
		}
	}

	_datasetLookup[name] = _datasets.size();
	_datasets.push_back(c);
	_entries.push_back(s);
	std::cout << "INFO in DatasetManager::sampleDataset : " << name << " loaded successfully (" << (s.size()==0?maxEntries:s.size()) << " events)" << std::endl;

	return true;
}

//get the TChain or list of entries for a named dataset
//TChain* DatasetManager::getDataset(TString name) {
//	if(_datasets.find(name)==_datasets.end()) return 0;
//	return _datasets[name];
//}
//
//std::multiset<unsigned int>* DatasetManager::getEntriesSet(TString name) {
//	if(_entries.find(name)==_entries.end()) return 0;
//	return &_entries[name];
//}

bool DatasetManager::setDataset(TString name, bool reset) {
	bool flip(false);
	int idx = findDataset(name, &flip);
	if(idx<0) return false;

	_currentSet = _datasets[idx];
	_nEntries = _currentSet->GetEntries();
	_flipped = flip;

	_currentEntryList = &_entries[idx];
	rewind();

	//optionally remove branch addresses
	if(reset) {
		_currentSet->ResetBranchAddresses();
	}
	return true;
}

bool DatasetManager::setBranchAddress(TString bname, void* obj) {
	if(!_currentSet) return false;
	_currentSet->SetBranchAddress(bname, obj);
	return true;
}

TString DatasetManager::getBranchType(TString bname) {
	if(!_currentSet) return "";
	TBranch* b = _currentSet->FindBranch("SVMCor");
	if(!b) return "";
	return b->GetClassName();
}

unsigned long DatasetManager::getEntries() {
	if(!_currentSet) return 0;

	if(_currentEntryList->empty()) {
		return _nEntries;
	} else if(_flipped) {
		return _nEntries-_currentEntryList->size();
	} else {
		return _currentEntryList->size();
	}
}

double DatasetManager::draw(TString var, TString sel, TH1* h) {
	if(!h) return 0.;
	h->Reset();
	rewind();
	TTreeFormula tfs("tfs",sel,_currentSet);
	TTreeFormula tfv("tfv",var,_currentSet);
	while(getNext()) {
		h->Fill(tfv.EvalInstance(),tfs.EvalInstance());
	}
	rewind();
	return h->Integral();
}

double DatasetManager::getEntries(TString sel) {
	double count(0);
	rewind();
	TTreeFormula tf("tf",sel,_currentSet);
	while(getNext()) {
		count += tf.EvalInstance();
	}
	rewind();
	return count;
}

std::vector<double> DatasetManager::getEntries(std::vector<TString> sels) {
	std::vector<double> counts;
	std::vector<TTreeFormula*> tfs;
	for(std::vector<TString>::iterator it=sels.begin(); it!=sels.end(); ++it) {
		tfs.push_back(new TTreeFormula((*it),(*it),_currentSet));
		counts.push_back(0);
	}
	rewind();
	while(getNext()) {
		std::vector<double>::iterator itC=counts.begin();
		std::vector<TTreeFormula*>::iterator itF=tfs.begin();
		for( ; itC!=counts.end(); ++itC,++itF) {
			(*itC) += (*itF)->EvalInstance();
		}
	}
	rewind();
	for(std::vector<TTreeFormula*>::iterator it=tfs.begin(); it!=tfs.end(); ++it) {
		delete(*it);
	}
	return counts;
}

bool DatasetManager::getNext() {
	if(!_currentSet) return false;

	//load entry before incrementing nextEntry
	//this allows the first entry to be loaded on first call
	if(_nextEntry>=_nEntries) return false;
	//TODO//std::cout << _nextEntry << "\t";//TODO
	_currentSet->GetEntry(_nextEntry);

	if(!_currentEntryList->empty()) {
		if(_flipped) {
			//in the flipped case we treat the entry list as a "skip list"
			do {
				++_nextEntry;
			} while(_currentEntryList->find(_nextEntry)!=_currentEntryList->end());
		} else {
			++_nextEntryIter;
			if(_nextEntryIter==_currentEntryList->end()) {
				_nextEntry=_nEntries;
			} else {
				_nextEntry = (*_nextEntryIter);
			}
		}
	} else {
		++_nextEntry;
	}
	return true;
}

void DatasetManager::rewind() {
	_nextEntry=0;
	if(!_currentEntryList->empty()) {
		_nextEntryIter = _currentEntryList->begin();
		if(_flipped) {
			while(_currentEntryList->find(_nextEntry)!=_currentEntryList->end()) {
				++_nextEntry;
			}
		} else {
			_nextEntry = (*_nextEntryIter);
		}
	}
	_currentSet->GetEntry(0);
}

void DatasetManager::reset() {
	if(_currentSet) _currentSet->ResetBranchAddresses();
}

int DatasetManager::findDataset(TString name, bool* flip) {
	TString notName;
	if(name.BeginsWith('!')) {
		notName=name.Strip(TString::kLeading,'!');
	} else {
		notName="!"+name;
	}

	if(_datasetLookup.find(name)!=_datasetLookup.end()) {
		if(flip) *flip=false;
		return _datasetLookup[name];
	}
	if(_datasetLookup.find(notName)!=_datasetLookup.end()) {
		if(flip) *flip=true;
		return _datasetLookup[notName];
	}
	return -1;
}

//int main() {
//	gRandom->SetSeed(124);
//	DatasetManager* dm = DatasetManager::getInstance();
//	dm->loadDataset("!test4", "T", std::vector<TString>{TString("/data/dijets/144/1/skimmed.root")}, 1800);
//	dm->loadDataset("test5", "T", std::vector<TString>{TString("/data/dijets/154/1/skimmed.root")}, 1600);
//	dm->buildDataset("test45", "T", std::vector<TString>{TString("test4"),TString("test5")});
//	//dm->loadDataset("test", "T", std::vector<TString>{TString("/data/dijets/144/1/skimmed.root"),TString("/data/dijets/154/1/skimmed.root")}, 2000);
//	//double JPT;
//	//dm->setDataset("test45");
//	//dm->setBranchAddress("JetPT",&JPT);
//	//for(int i=0; i<10; ++i) {
//	//	dm->getNext();
//	//	std::cout << JPT << std::endl;
//	//}
//	dm->setDataset("test4");
//	std::cout << dm->getEntries() << "\t" << dm->getEntries("JetPT>0") << std::endl;
//	std::vector<double> res = dm->getEntries(std::vector<TString>{"1","JetPT>0"});
//	std::cout << res[0] << "\t" << res[1] << std::endl;
//	TH1D h("h","",10,10e3,50e3);
//	dm->draw("JetPT","JetPT>0.",&h);
//	TCanvas c;
//	h.Draw();
//	c.SaveAs("testing.pdf");
//	dm->setDataset("test5");
//	std::cout << dm->getEntries() << std::endl;
//	dm->setDataset("!test4");
//	std::cout << dm->getEntries() << std::endl;
//	dm->setDataset("!test5");
//	std::cout << dm->getEntries() << std::endl;
//	dm->setDataset("test45");
//	std::cout << dm->getEntries() << std::endl;
//	dm->setDataset("!test45");
//	std::cout << dm->getEntries() << std::endl;
//	//dm->setBranchAddress("JetPT",&JPT);
//	//for(int i=0; i<100; ++i) {
//	//	if(!dm->getNext()) break;
//	//	std::cout << JPT << std::endl;
//	//}
//
//	return 0;
//}
