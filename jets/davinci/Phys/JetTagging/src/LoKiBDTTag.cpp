// Implementation file for class BDTTag.
// Written by Philip Ilten and Mike Williams, 2015-01-30.

// Event.
#include "Event/MCParticle.h"

// Local.
#include "LoKiBDTTag.h"

// Declaration of tool factory.
DECLARE_NAMESPACE_TOOL_FACTORY(LoKi, BDTTag)

//==========================================================================

// Jet tagging via the secondary vertex (SVR) method of LHCb-ANA-2014-074.

//--------------------------------------------------------------------------

// Standard constructor.

LoKi::BDTTag::BDTTag(const std::string& type, const std::string& name,
		     const IInterface* parent) : GaudiTool(type, name, parent) {
  declareInterface<IJetTagTool>(this);
  m_velStrc = 0; m_fitTool = 0; m_dstTool = 0; m_reader  = 0;
  declareProperty("Bdt0Weights", m_bdt0Weights =
		  "data/bdt_configs/LHCb_ANA_2014_076_BDT0.weights.xml",
		  "BDT0 TMVA weights filename and path from $JETTAGGINGROOT.");
  declareProperty("Bdt1Weights", m_bdt1Weights =
		  "data/bdt_configs/LHCb_ANA_2014_076_BDT1.weights.xml",
		  "BDT1 TMVA weights filename and path from $JETTAGGINGROOT.");
  declareProperty("TmvaOptions", m_tmvaOptions = "Silent",
		  "Options to pass to the TMVA reader.");
  declareProperty("FitName", m_fitName = "LoKi::VertexFitter",
		  "Name of the vertex fitter tool used to create SVRs.");
  declareProperty("DstName", m_dstName = "LoKi::DistanceCalculator:PUBLIC",
		  "Name of the distance calculator tool used to create SVRs.");
  declareProperty("DR", m_dr = 0.5,
		  "The maximum dR(SVR flight direction, jet momentum) "
		  "for linking two-body SVRs.");
  declareProperty("Backwards", m_backwards = false,
		  "If true, build backwards SVRs by reversing the SVR flight "
		  "direction.");
  declareProperty("PrtSelect", m_prtSelect = true,
		  "If true, apply the default selection to the particles.");
  declareProperty("NbvSelect", m_nbvSelect = true,
		  "If true, apply the default selection to the n-body SVRs.");
  declareProperty("NbvSort", m_nbvSort = "pt",
		  "Sort the n-body SVRs by \"pt\", \"bdt0\", or \"bdt1\".");
  declareProperty("TbvLocation", m_tbvLocation = "",
		  "Optional TES location of two-body SVRs. If not set, the "
		  "two-body SVRs will be automatically built.");
  declareProperty("PrtLocation", m_prtLocation =
		  "Phys/StdAllNoPIDsPions/Particles",
		  "TES location of particles used to build the SVRs.");
  declareProperty("PvrLocation", m_pvrLocation = RecVertexLocation::Primary,
		  "TES location of PVRs used when building the SVRs.");
  declareProperty("OdnLocation", m_odnLocation = ODINLocation::Default,
		  "TES location of ODIN.");
}

//--------------------------------------------------------------------------

// Standard destructor.

LoKi::BDTTag::~BDTTag() {}

//--------------------------------------------------------------------------

// Initialize the tagger.

StatusCode LoKi::BDTTag::initialize() {

  // Set the run and event number.
  m_runNumber = -1; m_evtNumber = -1;

  // Get the tools.
  if (!GaudiTool::initialize()) {
    error() << "Failed to initialize the GaudiTool." << endmsg;
    return StatusCode::FAILURE;
  }
  m_velStrc = getDetIfExists<DeVelo>(DeVeloLocation::Default);
  if (!m_velStrc) {
    error() << "Failed to retrieve the VELO structure." << endmsg;
    return StatusCode::FAILURE;
  }
  m_fitTool = tool<IVertexFit>(m_fitName, this);
  if (!m_fitTool) {
    error() << "Failed to create the vertex fitter." << endmsg;
    return StatusCode::FAILURE;
  }
  m_dstTool = tool<IDistanceCalculator>(m_dstName, this);
  if (!m_dstTool) {//TODO
    error() << "Failed to create the distance calculator." << endmsg;
    return StatusCode::FAILURE;
  }

  // Create the TMVA reader.
  if (m_reader) delete m_reader;
  m_reader = new TMVA::Reader(m_tmvaOptions);
  if (!m_reader) {
    error() << "Failed to create the TMVA reader." << endmsg;
    return StatusCode::FAILURE;
  }

  // Add the variables.
  m_reader->AddVariable("FDR",   &m_vars.m_fdrMin);
  m_reader->AddVariable("PTR",   &m_vars.m_ptSvrJet);
  m_reader->AddVariable("NTRK",  &m_vars.m_nTrk);
  m_reader->AddVariable("NTRKJ", &m_vars.m_nTrkJet);
  m_reader->AddVariable("DR",    &m_vars.m_drSvrJet);
  m_reader->AddVariable("Q",     &m_vars.m_absQSum);
  m_reader->AddVariable("M",     &m_vars.m_m);
  m_reader->AddVariable("MCOR",  &m_vars.m_mCor);
  m_reader->AddVariable("FD",    &m_vars.m_fdChi2);
  m_reader->AddVariable("SIP",   &m_vars.m_ipChi2Sum);

  // Book the two BDT methods.
  string root = getenv("JETTAGGINGROOT");
  string bdt0Weights = root + "/" + m_bdt0Weights;
  string bdt1Weights = root + "/" + m_bdt1Weights;
  ifstream bdt0File(bdt0Weights), bdt1File(bdt1Weights);
  if (!bdt0File.good()) {
    error() << "Cannot read the BDT0 file, " + bdt0Weights + "." << endmsg;
    return StatusCode::FAILURE;
  }
  if (!bdt1File.good()) {
    error() << "Cannot read the BDT1 file, " + bdt1Weights + "." << endmsg;
    return StatusCode::FAILURE;
  }
  m_reader->BookMVA("BDT0", bdt0Weights);
  m_reader->BookMVA("BDT1", bdt1Weights);
  return StatusCode::SUCCESS;
}

//--------------------------------------------------------------------------

// Finalize the tagger.

StatusCode LoKi::BDTTag::finalize() {
  if (m_reader) delete m_reader;
  return GaudiTool::finalize();
}

//--------------------------------------------------------------------------

// Calculate the tagger properties for the jet.

bool LoKi::BDTTag::calculateJetProperty(const Particle *jet,
					map<string, double> &props) {

  // Create the n-body SVRs.
  if (!jet) {
    error() << "The passed jet is not valid." << endmsg; return false;
  }
  if (!nbvs(jet)) {
   error() << "The n-body SVRs are not valid." << endmsg; return false;
  }

  // Fill the jet property information.
  int n(0);
  if (m_nbvSelect)
  {if (m_nbvs.size() > 0 && m_nbvs[0].info(n, jet, props)) ++n;}
  else
      for (int nbv = 0; nbv < (int)m_nbvs.size(); ++nbv)
          if (m_nbvs[nbv].info(n, jet, props)) ++n;

  props["Tag"] = n;
  props["extraInfo"] = 5200;
  if (n > 0) return true;
  else return false;
}

//--------------------------------------------------------------------------

// Calculate dR.

double LoKi::BDTTag::deltaR(const Particle *prt1, const Particle *prt2) {
  double phi1(prt1->momentum().Phi()), phi2(prt2->momentum().Phi());
  double dPhi(phi1 - phi2);
  while (dPhi > M_PI) dPhi -= 2*M_PI;
  while (dPhi <= -M_PI) dPhi += 2*M_PI;
  double dEta(prt1->momentum().Eta()- prt2->momentum().Eta());
  return sqrt(dPhi * dPhi + dEta * dEta);
}

double LoKi::BDTTag::deltaR(const Svr *svr, const Particle *prt, double sgn) {
  Gaudi::LorentzVector fv(sgn * svr->m_fv.X(), sgn * svr->m_fv.Y(),
			  sgn * svr->m_fv.Z(), 0);
  double phi1(fv.Phi()), phi2(prt->momentum().Phi());
  double dPhi(phi1 - phi2);
  while (dPhi > M_PI) dPhi -= 2*M_PI;
  while (dPhi <= -M_PI) dPhi += 2*M_PI;
  double dEta(fv.Eta()- prt->momentum().Eta());
  return sqrt(dPhi * dPhi + dEta * dEta);
}

//--------------------------------------------------------------------------

// Static comparison functions.

bool LoKi::BDTTag::comparePrtPt(const Particle *prt1, const Particle *prt2) {
  return prt1->momentum().Pt() > prt2->momentum().Pt();}
bool LoKi::BDTTag::compareSvrPt(const Svr &svr1, const Svr &svr2) {
  return svr1.m_prt.momentum().Pt() > svr2.m_prt.momentum().Pt();}
bool LoKi::BDTTag::compareSvrBdt0(const Svr &svr1, const Svr &svr2) {
  return svr1.m_bdt0 > svr2.m_bdt0;}
bool LoKi::BDTTag::compareSvrBdt1(const Svr &svr1, const Svr &svr2) {
  return svr1.m_bdt1 > svr2.m_bdt1;}

//--------------------------------------------------------------------------

// Return the vector of two-body SVRs.

const vector<LoKi::BDTTag::Svr> *LoKi::BDTTag::tbvs(bool force) {

  // Check if new event and if SVRs need re-building.
  const ODIN *odin = getIfExists<ODIN>(m_odnLocation);
  if (!odin)
    {error() << "Could not retrieve ODIN." << endmsg; m_tbvs.clear(); return 0;}
  if (!force && m_runNumber == (int)odin->runNumber() &&
      m_evtNumber == odin->eventNumber()) return &m_tbvs;
  m_runNumber = odin->runNumber();
  m_evtNumber = odin->eventNumber();
  m_tbvs.clear();

  // Retrieve the primary vertices.
  m_pvrs = getIfExists<RecVertices>(m_pvrLocation);
  if (!m_pvrs) {
    error() << "PVR location " + m_prtLocation + " not found." << endmsg;
    return 0;
  }

  // Create the two-body SVRs from the pre-built SVR TES location.
  if (m_tbvLocation != "") {
    Particles *svrs = getIfExists<Particles>(m_tbvLocation);
    if (!svrs) {
      error() << "SVR location " + m_tbvLocation + " not found." << endmsg;
      return 0;
    }
    for (Particles::iterator svr = svrs->begin(); svr != svrs->begin(); ++svr) {
      SmartRefVector<Particle> prts = (*svr)->daughters();
      if (prts.size() == 2) m_tbvs.push_back(Svr(this, prts[0], prts[1]));
    }
    return &m_tbvs;
  }

  // Retrieve the particles.
  Particles *prts = getIfExists<Particles>(m_prtLocation);
  if (!prts) {
    error() << "Particle location " + m_prtLocation + " not found." << endmsg;
    return 0;
  }

  // Select and sort the particles.
  // A lower minimum IP chi-squared requirement is applied for MC.
  MCParticles *gens   = getIfExists<MCParticles>(MCParticleLocation::Default);
  double ipChi2MinCut = !gens ? 6 : 6;//TODO lowered from 9
  m_prts.clear();
  for (Particles::iterator prt = prts->begin(); prt != prts->end();
       prt++) {
    const Track *trk((*prt)->proto() ? (*prt)->proto()->track() : 0);
    if (m_prtSelect) {
      if (!trk)                             continue;
      if (!(trk->pt() > 500))               continue;
      if (!(trk->type() == Track::Long))    continue;
      if (!(trk->ghostProbability() < 0.3)) continue;
    
      // Calculate track minimum IP chi-squared.
      double ip(-1), ipChi2(-1), ipChi2Min(-1);
      for (RecVertices::iterator pvr = m_pvrs->begin();
	   pvr != m_pvrs->end(); ++pvr) {
	m_dstTool->distance(*prt, *pvr, ip, ipChi2);
	if (ipChi2 >= 0 && (ipChi2Min == -1 || ipChi2 < ipChi2Min))
	  ipChi2Min = ipChi2;
      }
      if (!(ipChi2Min > ipChi2MinCut)) continue;
    }
    m_prts.push_back(*prt);
  }

//  //TODO
//  Particles::iterator prt = prts->begin();
//  Gaudi::LorentzVector p = (*prt)->momentum();
//  Gaudi::XYZPoint x = (*prt)->referencePoint();
//  double ip(-1), ipChi2(-1), ipChi2Min(-1), ipBest(-1);
//  for (RecVertices::iterator pvr = m_pvrs->begin();
//       pvr != m_pvrs->end(); ++pvr) {
//    m_dstTool->distance(*prt, *pvr, ip, ipChi2);
//    if (ipChi2 >= 0 && (ipChi2Min == -1 || ipChi2 < ipChi2Min)) {
//      ipChi2Min = ipChi2;
//      ipBest = ip;
//    }
//  }
//  info() << "pt is " << p.Pt() << ", (x, y, z) is (" << x.x() << ", " << x.y() << ", " << x.z() << ") and best IP is " << ipBest << endmsg;
//
//  for(int i=0; i<2; ++i) {//TODO
////	  p.SetPx(1.1*p.Px());
////	  p.SetPy(1.1*p.Py());
////	  (*prt)->setMomentum(p);
////	  x.SetX(1.1*x.x());
////	  x.SetY(1.1*x.y());
//	  x.SetZ(1.1*x.z());
//	  (*prt)->setReferencePoint(x);
//    
//          // Calculate track minimum IP chi-squared.
//	  ipChi2Min=-1;
//	  ipBest=-1;
//          for (RecVertices::iterator pvr = m_pvrs->begin();
//               pvr != m_pvrs->end(); ++pvr) {
//            m_dstTool->distance(*prt, *pvr, ip, ipChi2);
//            if (ipChi2 >= 0 && (ipChi2Min == -1 || ipChi2 < ipChi2Min)) {
//              ipChi2Min = ipChi2;
//	      ipBest = ip;
//	    }
//          }
//  	  p = (*prt)->momentum();
//  	  x = (*prt)->referencePoint();
//  	  info() << "pt is " << p.Pt() << ", (x, y, z) is (" << x.x() << ", " << x.y() << ", " << x.z() << ") and best IP is " << ipBest << endmsg;
//  }//TODO

//  info() << "Nparts " << m_prts.size() << endmsg; //TODO
//  int nRerolls(0);
  //reroll IPs if needed
//  while(m_tbvs.empty()) {
  	sort(m_prts.begin(), m_prts.end(), comparePrtPt);

//	info() << "***" << endmsg;
  	// Create and select the two-body SVRs.
  	for (Particles::iterator prt1 = m_prts.begin(); prt1 != m_prts.end(); ++prt1) {
  	  for (Particles::iterator prt2 = prt1 + 1; prt2 != m_prts.end(); ++prt2) {
		  //TODO added 3-body vertices without linking (for charm)
  	  	for (Particles::iterator prt3 = prt2 + 1; prt3 != m_prts.end(); ++prt3) {
		    Svr tbv(this, *prt1, *prt2, *prt3);
		    if (tbv.pass()) m_tbvs.push_back(tbv);
		}
		    Svr tbv(this, *prt1, *prt2);
		    if (tbv.pass()) m_tbvs.push_back(tbv);
	    }
	}
//  	    {Svr tbv(this, *prt1, *prt2); if (tbv.pass()) m_tbvs.push_back(tbv); info() << "VRTCHI2" << tbv.vrtchi2() << endmsg;}//TODO

//	//limit to 99 rerolls
//	if(nRerolls>=9) break;
//  	//reroll IPs if needed
//  	if(m_tbvs.empty()) {
//		++nRerolls;
//  		for (Particles::iterator prt = m_prts.begin(); prt != m_prts.end(); ++prt) {
//  			Gaudi::LorentzVector p = (*prt)->momentum();
//			Gaudi::XYZPoint x = (*prt)->referencePoint();
//			x.SetZ(1.1*x.z());
//	  		p.SetPx(1.1*p.Px());
//	  		p.SetPy(1.1*p.Py());
//	  		(*prt)->setMomentum(p);
//			(*prt)->setReferencePoint(x);
//  	  		info() << "pT is " << p.Pt() << " and (x, y, z) is (" << x.x() << ", " << x.y() << ", " << x.z() << ")" << endmsg;
//		}
//	}
//  }
//  info() << "Rerolled particle IPs " << nRerolls << " times." << endmsg;

  return &m_tbvs;
}

//--------------------------------------------------------------------------

// Return the vector of n-body SVRs.

const vector<LoKi::BDTTag::Svr> *LoKi::BDTTag::nbvs(const Particle *jet) {

  // Create the two-body SVRs.
  m_nbvs.clear();
  if (!tbvs()) {
    error() << "The two-body SVRs are not valid." << endmsg; return 0;
  }
  if (!jet) {error() << "The jet is not valid." << endmsg; return 0;}

  // Link the two-body SVRs.
  for (int tbv = 0; tbv < (int)m_tbvs.size(); ++tbv) {
    if (deltaR(&m_tbvs[tbv], jet, m_backwards ? -1 : 1) > m_dr) continue;
    for (int trk = 0; trk < (int)m_tbvs[tbv].m_trks.size(); ++trk)
      if (deltaR(m_tbvs[tbv].m_trks[trk], jet) <= m_dr)
	{m_nbvs.push_back(Svr(&m_tbvs[tbv])); break;}
  }
  //TODO linking is turned off on but originals also kept
  double n(0);
  while (m_nbvs.size() != n) {
    n = m_nbvs.size();
    for (int nbv1 = 0; nbv1 < (int)m_nbvs.size(); ++nbv1)
      for (int nbv2 = nbv1 + 1; nbv2 < (int)m_nbvs.size();)
        if (m_nbvs[nbv1].link(m_nbvs[nbv2]))
          m_nbvs.erase(m_nbvs.begin() + nbv2);
        else ++nbv2;
  }
  //Add in the originals
  for (int tbv = 0; tbv < (int)m_tbvs.size(); ++tbv) {
    if (deltaR(&m_tbvs[tbv], jet, m_backwards ? -1 : 1) > m_dr) continue;
    for (int trk = 0; trk < (int)m_tbvs[tbv].m_trks.size(); ++trk)
      if (deltaR(m_tbvs[tbv].m_trks[trk], jet) <= m_dr)
	{m_nbvs.push_back(Svr(&m_tbvs[tbv])); break;}
  }
  //TODO Code below will remove 2-body vertices that are subsets of 3-body vertices
  //double n(0);
  //while (m_nbvs.size() != n) {
  //  n = m_nbvs.size();
  //  for (int nbv1 = 0; nbv1 < (int)m_nbvs.size(); ++nbv1) {
  //    for (int nbv2 = 0; nbv2 < (int)m_nbvs.size();) {
  //       if(nbv1==nbv2) ++nbv2;
  //       else if (m_nbvs[nbv1].removeSubset(m_nbvs[nbv2]))
  //        m_nbvs.erase(m_nbvs.begin() + nbv2);
  //      else ++nbv2;
  //    }
  //  }
  //}

  // Sort the n-body SVRs.
  for (int nbv = 0; nbv < (int)m_nbvs.size(); ++nbv) m_nbvs[nbv].calc(jet);
  if (m_nbvSort == "pt")
    sort(m_nbvs.begin(), m_nbvs.end(), compareSvrPt);
  else if (m_nbvSort == "bdt0")
    sort(m_nbvs.begin(), m_nbvs.end(), compareSvrBdt0);
  else if (m_nbvSort == "bdt1")
    sort(m_nbvs.begin(), m_nbvs.end(), compareSvrBdt1);
  return &m_nbvs;
}

//==========================================================================

// Internal SVR representation for the LoKi::BDTTag jet tagger class.

//--------------------------------------------------------------------------

// Two-body constructor.

LoKi::BDTTag::Svr::Svr(BDTTag *parent, const Particle *prt1,
		       const Particle *prt2,
		       const Particle *prt3) {
	//if(prt3) std::cout << *prt1 << "\n" << *prt2 << "\n" << *prt3 << "\n" << std::endl;
	//else std::cout << *prt1 << "\n" << *prt2 << "\n" << std::endl;
  m_parent = parent;
  m_trks.push_back(prt1);
  m_trks.push_back(prt2);
  if(prt3) m_trks.push_back(prt3);
  m_stored = false;
  //std::cout << "??";
  //for (int trk = 0; trk < (int)m_trks.size(); ++trk) {
  //  std::cout << "\t" << m_trks[trk]->key();
  //}
  //std::cout << std::endl;
}

//--------------------------------------------------------------------------

// N-body constructor.

LoKi::BDTTag::Svr::Svr(const Svr *svr) {
  m_parent = svr->m_parent;
  m_trks   = svr->m_trks;
  if (svr->m_tbvs.size() == 0) m_tbvs.push_back(svr);
  else m_tbvs = svr->m_tbvs;
  m_stored = false;
}

//--------------------------------------------------------------------------

// Check if the two-body SVR passes the quality requirements.

bool LoKi::BDTTag::Svr::pass() {
  if (!calc())                          return false;
  if (!m_fit)                           return false;
  if (!(m_docaMin < 0.2))               return false;
  if (!(m_vrt.chi2() < 10))             return false;
  if (!(m_prt.momentum().M() < 5000))   return false;
  if (!(m_mCor > 600))                  return false;
  if (!(m_fdChi2 > 25))                 return false;
  if (!(m_tz >= 0 && m_tz < 10))        return false;
  if (!(m_hits == 0))                   return false;
  //std::cout << "!!";
  //for (int trk = 0; trk < (int)m_trks.size(); ++trk) {
  //  std::cout << "\t" << m_trks[trk]->key();
  //}
  //std::cout << std::endl;
  return true;
}

//--------------------------------------------------------------------------

// Write the n-body SVR info, if it passes the quality requirements.

bool LoKi::BDTTag::Svr::info(int idx, const Particle *jet,
			     map<string, double> &props) {

  // Check if the SVR passes the n-body requirements.
  bool pass(true);
  if (!calc(jet)) return false;
  if (m_trks.size() == 2 && !(abs(m_prt.momentum().M() - 500.0) > 20))
    pass = false;
  else if (!(m_fdrMin < 15))                pass = false;
  else if (!(m_nTrkJet > 1))                pass = false;
  else if (!(m_vrt.position().Z() < 200))   pass = false;
  else if (!(m_mCor > 600))                 pass = false;
  else if (!(m_prt.momentum().Pt() > 2000)) pass = false;
  else if (!(m_tau < 1.5))                  pass = false;
  else if (!(m_fdChi2 > 32))                pass = false;
  else if (!(m_drSvrJet <= m_parent->m_dr)) pass = false;
  if (m_parent->m_nbvSelect && !pass) return false;

//  m_parent->error() << m_prt.momentum().Px() << "\t" << m_prt.momentum().Py() << "\t" << m_prt.momentum().Pz() << "\t" << m_prt.momentum().Pt() << "\t" << m_prt.momentum().E() << "\t" << m_prt.momentum().M() << endmsg;
//  std::cout << "===================================="<<std::endl;

  // Fill the BDT properties.
  stringstream pre; pre << "Tag" << idx << "_";
  props[pre.str() + "fdrMin"] 	   = m_fdrMin;
  props[pre.str() + "ptSvrJet"]    = m_ptSvrJet;
  props[pre.str() + "nTrk"]	   = m_trks.size();
  props[pre.str() + "nTrkJet"]     = m_nTrkJet;
  props[pre.str() + "drSvrJet"]    = m_drSvrJet;
  props[pre.str() + "absQSum"]     = m_absQSum;
  props[pre.str() + "m"]	   = m_prt.momentum().M();
  props[pre.str() + "mCor"]        = m_mCor;
  props[pre.str() + "mCorErr"]     = m_mCorErr;
  props[pre.str() + "mCorErrFull"] = m_mCorErrFull;
  props[pre.str() + "fdChi2"]      = m_fdChi2;
  props[pre.str() + "ipChi2Min"]   = m_ipChi2Min;
  props[pre.str() + "ipChi2Sum"]   = m_ipChi2Sum;
  props[pre.str() + "ipChi2MinTrk"]= m_ipChi2MinTrk;

  // Fill the additional properties.
  props[pre.str() + "bdt0"]      = m_bdt0;
  props[pre.str() + "bdt1"]      = m_bdt1;
  props[pre.str() + "pass"]      = pass;
  props[pre.str() + "tau"]       = m_tau;
  props[pre.str() + "z"]         = m_vrt.position().Z();
  props[pre.str() + "pt"]        = m_prt.momentum().Pt();
  props[pre.str() + "backwards"] = m_parent->m_backwards;

  // Properties needed for training.
  props[pre.str() + "x"]  = m_vrt.position().X();
  props[pre.str() + "y"]  = m_vrt.position().Y();
  props[pre.str() + "px"] = m_prt.momentum().Px();
  props[pre.str() + "py"] = m_prt.momentum().Py();
  props[pre.str() + "pz"] = m_prt.momentum().Pz();
  props[pre.str() + "idx_pvr"] = m_pvr->key();
  for (int trk = 0; trk < (int)m_trks.size(); ++trk) {
    stringstream pst; pst << "idx_trk" << trk;
    props[pre.str() + pst.str()] = m_trks[trk]->key();
  }

  props[pre.str() + "chi2"] = m_vrt.chi2();

  return true;
}

//--------------------------------------------------------------------------

// Calculate the stored variables. Returns true if successful.

bool LoKi::BDTTag::Svr::calc(const Particle *jet, bool force) {
  if (m_stored && !force) return m_stored;
  else m_stored = false;
  m_fv.SetXYZ(0, 0, 0);
  m_hits      = 0;  m_tz        = -1;
  m_fdrMin    = -1; m_absQSum   = 0;  m_ipMin     = -1; m_dsMax     = -1;
  m_ptSvrJet  = -1; m_mCor      = -1; m_ipChi2Min = -1;	m_tau       = -1;
  m_drSvrJet  = -1; m_fdChi2    = -1; m_docaMin   = -1;	m_bdt0      = -999;
  m_nTrkJet   = 0;  m_ipChi2Sum = 0;  m_drMax     = -1;	m_bdt1      = -999;
  if (!m_parent || !m_parent->m_fitTool || !m_parent->m_dstTool ||
      !m_parent->m_pvrs) return m_stored;
  else m_stored = true;
  m_fit = m_parent->m_fitTool->fit(m_trks, m_vrt, m_prt).isSuccess();

  // Determine position using two-body chi-squared weighted positions.
  TVector3 svrPos(0, 0, 0);
  if (m_tbvs.size() < 2)
    svrPos.SetXYZ(m_vrt.position().X(), m_vrt.position().Y(),
		  m_vrt.position().Z());
  else {
    double sumWeights(0), weight(0);
    for (int tbv1 = 0; tbv1 < (int)m_tbvs.size(); ++tbv1) {
      weight      = 1.0 / m_tbvs[tbv1]->m_vrt.chi2();
      sumWeights += weight;
      svrPos += TVector3(m_tbvs[tbv1]->m_vrt.position().X(),
			 m_tbvs[tbv1]->m_vrt.position().Y(),
			 m_tbvs[tbv1]->m_vrt.position().Z()) * weight;

      // Minimum two-body radial flight distance.
      if (m_fdrMin == -1 || m_tbvs[tbv1]->m_fdrMin < m_fdrMin)
	m_fdrMin = m_tbvs[tbv1]->m_fdrMin;

      // Minimum two-body separation and dR.
      for (int tbv2 = tbv1 + 1; tbv2 < (int)m_tbvs.size(); ++tbv2) {
	double ds = sqrt((m_tbvs[tbv1]->m_vrt.position() -
			  m_tbvs[tbv2]->m_vrt.position()).Mag2());
	double dr = deltaR(&m_tbvs[tbv1]->m_prt, &m_tbvs[tbv2]->m_prt);
	if (ds > m_dsMax) m_dsMax = ds;
	if (dr > m_drMax) m_drMax = dr;
      }
    }
    svrPos *= (1.0 / sumWeights);
    m_vrt.setPosition(Gaudi::XYZPoint(svrPos.X(), svrPos.Y(), svrPos.Z()));
    m_prt.setReferencePoint(Gaudi::XYZPoint(svrPos.X(), svrPos.Y(),
					    svrPos.Z()));
  }

  // Correct the momentum and check track hit positions.
  Gaudi::LorentzVector pv4(0, 0, 0, 0);
  for (int prt = 0; prt < (int)m_trks.size(); ++prt) {
    if (m_trks[prt]->proto() && !m_trks[prt]->proto()->track()) {
      const Track *trk = m_trks[prt]->proto()->track();
      for (vector<LHCbID>::const_iterator id = trk->lhcbIDs().begin(); 
	   id != trk->lhcbIDs().end(); ++id)
	if (id->isVelo() && m_parent->m_velStrc->
	    sensor(id->veloID())->z() < svrPos.z()) ++m_hits;
    }
    pv4 += m_trks[prt]->momentum();
  }
  m_prt.setMomentum(pv4);

  // Impact parameter and flight distance information.
  TVector3 pvrPos(0, 0, 0); double ip(-1), ipChi2(-1), fd(-1), fdr(-1);
  for (RecVertices::iterator pvr = m_parent->m_pvrs->begin();
       pvr != m_parent->m_pvrs->end(); ++pvr) {
    m_parent->m_dstTool->distance(&m_prt, *pvr, ip, ipChi2);
    if (ipChi2 >= 0 && (m_ipChi2Min == -1 || ipChi2 < m_ipChi2Min))
      m_ipChi2Min = ipChi2;
    if (ip >= 0 && (m_ipMin == -1 || ip < m_ipMin)) {
      m_pvr = *pvr;
      pvrPos.SetXYZ((*pvr)->position().X(), (*pvr)->position().Y(),
		    (*pvr)->position().Z());
      m_ipMin = ip;
      m_parent->m_dstTool->distance(&m_vrt, *pvr, fd, m_fdChi2);
    }
    fdr = sqrt((m_vrt.position() - (*pvr)->position()).perp2());
    if (m_tbvs.size() < 2 && (m_fdrMin == -1 || fdr < m_fdrMin)) m_fdrMin = fdr;
  }

  // Calculate the pseudo-lifetime.
  m_tz = (svrPos.Z() - pvrPos.Z())*m_prt.momentum().M()/m_prt.momentum().Pz()
    /(3e11)*(1e12);

  // DOCA information.
  double doca;
  for (int trk1 = 0; trk1 < (int)m_trks.size(); ++trk1)
    for (int trk2 = trk1 + 1; trk2 < (int)m_trks.size(); ++trk2)
      if (m_parent->m_dstTool->
	  distance(m_trks[trk1], m_trks[trk2], doca).isSuccess())
  	if (m_docaMin == -1 || doca < m_docaMin) m_docaMin = doca;

  double minall(-1);
  // IP chi-squared and absolute charge sums.
  for (int trk = 0; trk < (int)m_trks.size(); ++trk) {
    double min(-1), ip, chi2;
    m_absQSum += m_trks[trk]->charge();
    for (RecVertices::iterator pvr = m_parent->m_pvrs->begin();
	 pvr != m_parent->m_pvrs->end(); ++pvr) {
      m_parent->m_dstTool->distance(m_trks[trk], *pvr, ip, chi2);
      if (chi2 >= 0 && (min == -1 || chi2 < min)) min = chi2;
    }
    m_ipChi2Sum += min;
    if(minall == -1 || min < minall) minall = min;
  }
  m_ipChi2MinTrk = minall;
  m_absQSum = abs(m_absQSum);

  // Corrected mass.
  m_fv = svrPos - pvrPos;
  TVector3 pv3(pv4.Px(), pv4.Py(), pv4.Pz());
  double pt2 = pv3.Cross(m_fv.Unit()).Mag2();
  m_mCor = sqrt(pv4.M2() + pt2) + sqrt(pt2);
  m_tau  = m_fv.Mag() / pv3.Mag() * 1000.0;

  //Corrected mass error //TODO
  mCorErrors(svrPos, pvrPos, pv4, m_prt.covMatrix(), m_pvr->covMatrix());
  //TODO

  // Calculate jet variables.
  if (!jet) return m_stored;
  m_ptSvrJet = m_prt.momentum().Pt() / jet->momentum().Pt();
  m_drSvrJet = deltaR(this, jet, m_parent->m_backwards ? -1 : 1);
  for (int trk = 0; trk < (int)m_trks.size(); ++trk)
    if (deltaR(m_trks[trk], jet) <= m_parent->m_dr) ++m_nTrkJet;

  // Fill the BDT structure.
  m_parent->m_vars.m_fdrMin    = m_fdrMin;
  m_parent->m_vars.m_ptSvrJet  = m_ptSvrJet;
  m_parent->m_vars.m_nTrk      = m_trks.size();
  m_parent->m_vars.m_nTrkJet   = m_nTrkJet;
  m_parent->m_vars.m_drSvrJet  = m_drSvrJet;
  m_parent->m_vars.m_absQSum   = m_absQSum;
  m_parent->m_vars.m_m         = m_prt.momentum().M();
  m_parent->m_vars.m_mCor      = m_mCor;
  m_parent->m_vars.m_fdChi2    = m_fdChi2;
  m_parent->m_vars.m_ipChi2Sum = m_ipChi2Sum;

  // Calculate the BDT responses.
  double bdt0(m_parent->m_reader->EvaluateMVA("BDT0"));
  double bdt1(m_parent->m_reader->EvaluateMVA("BDT1"));
  if (bdt0 != -999)
    m_bdt0 = (1 - TMath::Prob(pow(abs(bdt0)/0.2, 2), 1)) * (bdt0 < 0 ? -1 : 1);
  if (bdt1 != -999)
    m_bdt1 = (1 - TMath::Prob(pow(abs(bdt1)/0.2, 2), 1)) * (bdt1 < 0 ? -1 : 1);
  return m_stored;
}

//--------------------------------------------------------------------------

// Link an n-body SVR with another n-body SVR, failure returns false.

bool LoKi::BDTTag::Svr::link(Svr &nbv) {
  if (m_tbvs.size() == 0 || nbv.m_tbvs.size() == 0) return false;
  set<const Particle*> cmb;
  for (int trk = 0; trk < (int)m_trks.size(); ++trk)
    cmb.insert(m_trks[trk]);
  for (int trk = 0; trk < (int)nbv.m_trks.size(); ++trk)
    cmb.insert(nbv.m_trks[trk]);
  if (cmb.size() >= m_trks.size() + nbv.m_trks.size()) return false;
  m_trks.clear();
  for (set<const Particle*>::iterator trk = cmb.begin(); trk != cmb.end();
       ++trk) m_trks.push_back(*trk);
  for (int tbv = 0; tbv < (int)nbv.m_tbvs.size(); ++tbv)
    m_tbvs.push_back(nbv.m_tbvs[tbv]);
  m_stored = false;
  return true;
}

//==========================================================================
bool LoKi::BDTTag::Svr::removeSubset(Svr &nbv) {
  if (m_tbvs.size() == 0 || nbv.m_tbvs.size() == 0) return false;

  if (m_trks.size() <= nbv.m_trks.size()) return false;

  //check if nbv is a subset of this Svr
  for (int trk1 = 0; trk1 < (int)nbv.m_trks.size(); ++trk1) {
  	bool found(false);
  	for (int trk2 = 0; trk2 < (int)m_trks.size(); ++trk2) {
  		if(nbv.m_trks[trk1] == m_trks[trk2]) {
  			found=true;
  			break;
  		}
  	}
  	if(!found) return false;
  }
  return true;
}

//==========================================================================
//// Include files
//
//// from Gaudi
//#include "GaudiKernel/ToolFactory.h"
//
//// local
//#include "TupleToolSLTools.h"
//
//#include <Kernel/GetIDVAlgorithm.h>
//#include <Kernel/IDVAlgorithm.h>
//#include <Kernel/IDistanceCalculator.h>
//#include <Kernel/IVertexFit.h>
//#include <Kernel/ILifetimeFitter.h>
//#include "Kernel/IPVReFitter.h"
//#include "GaudiAlg/Tuple.h"
//#include "GaudiAlg/TupleObj.h"
//
//#include "TROOT.h"
//#include "Event/Particle.h"
//#include "TLorentzVector.h"
//#include "TVector3.h"
//
//using namespace LHCb;
//
////-----------------------------------------------------------------------------
//// Implementation file for class : TupleToolSLTools
////
//// @author Mitesh Patel, Patrick Koppenburg
//// @date   2008-04-15
////-----------------------------------------------------------------------------
//
//// Declaration of the Tool Factory
//// actually acts as a using namespace TupleTool
//DECLARE_TOOL_FACTORY( TupleToolSLTools )
//
////=============================================================================
//// Standard constructor, initializes variables
////=============================================================================
//  TupleToolSLTools::TupleToolSLTools( const std::string& type,
//                                        const std::string& name,
//                                        const IInterface* parent )
//    : TupleToolBase ( type, name , parent )
//    , m_dva(0)
//    , m_dist(0)
//    , m_pVertexFit(0)
//{
//  declareInterface<IParticleTupleTool>(this);
//  declareProperty("VertexFitter", m_typeVertexFit = "LoKi::VertexFitter" );
//  declareProperty("Bmass", m_Bmass = 5620.2 );
//  declareProperty("VertexCov",m_vcov=false);
//  declareProperty("MomCov",m_momcov=false);
//}
//
////=============================================================================
//
//StatusCode TupleToolSLTools::initialize()
//{
//  const StatusCode sc = TupleToolBase::initialize();
//  if ( sc.isFailure() ) return sc;
//
//  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc(), this ) ;
//  if (!m_dva) return Error("Couldn't get parent DVAlgorithm",
//                           StatusCode::FAILURE);
//
//  m_dist = m_dva->distanceCalculator ();
//  if( !m_dist ){
//    Error("Unable to retrieve the IDistanceCalculator tool");
//    return StatusCode::FAILURE;
//  }
//
//  //m_pVertexFit= m_dva->vertexFitter();
//  m_pVertexFit  = tool<IVertexFit>(m_typeVertexFit,this);
//  //m_pVertexFit  = tool<IVertexFit>("OfflineVertexFitter",this);
//  if( !m_pVertexFit ){
//    Error("Unable to retrieve the IVertexFit tool");
//    return StatusCode::FAILURE;
//  }
//  m_ltfit = tool<ILifetimeFitter>( "PropertimeFitter", this );
//  if( !m_ltfit ){
//    Error("Unable to retrieve the ILifetimeFitter tool");
//    return StatusCode::FAILURE;
//  }
//  m_pvReFitter = tool<IPVReFitter>("AdaptivePVReFitter", this );
//
//  return sc;
//}
//
////=============================================================================
//
//StatusCode TupleToolSLTools::fill( const Particle* mother
//                                    , const Particle* P
//                                    , const std::string& head
//                                    , Tuples::Tuple& tuple )
//{
//
//  const std::string prefix=fullName(head);
//  Assert( P && mother && m_dist
//          , "This should not happen, you are inside TupleToolSLTools.cpp :(" );
//
//
//  bool test=true;
//  // --------------------------------------------------
//
//  // find the origin vertex. Either the primary or the origin in the
//  // decay
//  /*
//    const VertexBase* vtx = 0;
//    if( mother != P ) vtx = originVertex( mother, P );
//    if( !vtx ){
//    Error("Can't retrieve the origin vertex for " + prefix );
//    return StatusCode::FAILURE;
//    }
//  */
//  const LHCb::Particle* mu_part;
//  LHCb::Particle::ConstVector source;
//  LHCb::Particle::ConstVector target;
//
//  double mcorr, q2_one, q2_two, mcorrerr, mcorrfullerr;
//  std::vector<double> mcorr_errors;
//  std::vector<TLorentzVector> nu_slns;
//  const LHCb::Vertex* pmu_vert;
//  const LHCb::VertexBase* PV;
//
//
//  if (P->isBasicParticle()){
//    source.push_back(mother);
//  }
//  else{
//    source.push_back(P);
//  }
//
//  do {
//    target.clear();
//    for(LHCb::Particle::ConstVector::const_iterator isource = source.begin(); 
//        isource != source.end(); isource++){
//      
//      if(!((*isource)->daughters().empty())){
//        
//        LHCb::Particle::ConstVector tmp = (*isource)->daughtersVector();
//        
//        for( LHCb::Particle::ConstVector::const_iterator itmp = tmp.begin(); 
//             itmp!=tmp.end(); itmp++){
//          target.push_back(*itmp);
//          // Add the final states, i.e. particles with proto and ignoring gammas
//          if((*itmp)->proto() && 22 != (*itmp)->particleID().pid()){
//		if ((*itmp)->particleID().abspid() == 13)
//		{
//		mu_part = (*itmp);
//		}	
//           }
//        }
//      } // if endVertex
//    } // isource
//    source = target;
//  } while(target.size() > 0);
//
//
//  // vetex selected proton and muon 
//  
//
//    pmu_vert = P->endVertex();
//    PV = m_dva->bestVertex ( P );
//
//    Gaudi::LorentzVector LV_P,LV_mu;
//    LV_P = P->momentum();
//    TLorentzVector TLV_P(LV_P.px(),LV_P.py(),LV_P.pz(),LV_P.E());
//    // muon LV
//    LV_mu = mu_part->momentum();
//    TLorentzVector TLV_mu(LV_mu.px(),LV_mu.py(),LV_mu.pz(),LV_mu.E());
//
//
//    TVector3 TV3_PV(PV->position().X(),PV->position().Y(), PV->position().Z());
//    TVector3 TV3_SV(pmu_vert->position().X(),pmu_vert->position().Y(), pmu_vert->position().Z());
//    
//    TVector3 Mdirn = (TV3_SV - TV3_PV).Unit();
//
//    mcorr = Mcorr(TLV_P, Mdirn);
//    mcorr_errors = Mcorr_errors(TV3_SV, TV3_PV, TLV_P, P->covMatrix(), PV->covMatrix());
//    nu_slns = recoNu(TLV_P, Mdirn, m_Bmass); 
//
//    fillVertex(pmu_vert, prefix +"_SV", tuple);
//    fillVertex(PV, prefix + "_PV", tuple);
//
//    if(m_momcov)
//    {   
//    const Gaudi::SymMatrix4x4 & mom_cov = P->momCovMatrix ();
//    const Gaudi::Matrix4x3 & posmom_cov = P->posMomCovMatrix ();
//    test &= tuple->matrix( prefix +  "_MOM_COV_", mom_cov );
//    test &= tuple->matrix( prefix +  "_POSMOM_COV_", posmom_cov );
//}
//
//    test &= tuple->column(  prefix + "_MCORR", mcorr);
//    if(mcorr_errors.size() ==  2)
//    {
//    mcorrerr = mcorr_errors[0];
//    mcorrfullerr = mcorr_errors[1];
//    test &= tuple->column(  prefix + "_MCORRERR", mcorrerr);
//    test &= tuple->column(  prefix + "_MCORRFULLERR", mcorrfullerr);
//    }
//    else
//    {
//    test &= tuple->column(  prefix + "_MCORRERR", -1000.);
//    test &= tuple->column(  prefix + "_MCORRFULLERR", -1000.);
//    }
//
//    if(nu_slns.size() ==  2)
//    {
//    q2_one = (nu_slns[0] + TLV_mu).M2();
//    q2_two = (nu_slns[1] + TLV_mu).M2();
//
//    test &= tuple->column(  prefix + "_Q2SOL1", q2_one);
//    test &= tuple->column(  prefix + "_Q2SOL2", q2_two);
//    }
//    else
//    {
//    test &= tuple->column(  prefix + "_Q2SOL1", -1000.);
//    test &= tuple->column(  prefix + "_Q2SOL2", -1000.);
//    }
//    debug() << "PAST FILL " << endreq;
//    //}
//    // replaced by V.B. 20.Aug.2k+9: (parts2Vertex,vtxWithExtraTrack);
//    // Remove the added track from parts2Vertex
//
//    
//
//  return StatusCode(test);
//}
//
////=========================================================================
////
////=========================================================================
//StatusCode TupleToolSLTools::fillVertex( const LHCb::VertexBase* vtx,
//                                          const std::string& vtx_name,
//                                          Tuples::Tuple& tuple ) const
//{
//  bool test = true ;
//
//  // decay vertex information:
//  if ( !vtx )
//  {
//    Gaudi::XYZPoint pt(-999.,-999.,-999.) ; // arbitrary point
//    test &= tuple->column(  vtx_name+"_", pt );
//    test &= tuple->column(  vtx_name + "_XERR", -999. );
//    test &= tuple->column(  vtx_name + "_YERR", -999. );
//    test &= tuple->column(  vtx_name + "_ZERR", -999. );
//    test &= tuple->column(  vtx_name + "_CHI2", -999. );
//    test &= tuple->column(  vtx_name + "_NDOF", -1 );
//    test &= tuple->matrix(  vtx_name + "_COV_", Gaudi::SymMatrix3x3()  );
//  }
//  else
//  {
//    test &= tuple->column( vtx_name+"_", vtx->position() );
//    test &= tuple->column(  vtx_name + "_CHI2", vtx->chi2() );
//    test &= tuple->column(  vtx_name + "_NDOF", vtx->nDoF() );
//    if(m_vcov)
//    {
//    const Gaudi::SymMatrix3x3 & m = vtx->covMatrix ();
//    test &= tuple->matrix(  vtx_name + "_COV_", m );
//    }
////    test &= tuple->matrix(  vtx_name + "_COV_", m );
//  }
//  if (!test) Warning("Error in fillVertex "+vtx_name,1).ignore();
//  return StatusCode(test) ;
//}
//
//double TupleToolSLTools::Mcorr(TLorentzVector Y, TVector3 M_dirn)
//{
//     double P_T = (Y.Vect()).Perp(M_dirn);
//     return  sqrt(Y.M2() + P_T * P_T) + P_T;
//}
//
//
//
//std::vector<TLorentzVector> TupleToolSLTools::recoNu(TLorentzVector Y, TVector3 M_dirn, double mass)
//    {
//	std::vector<TLorentzVector> p_nu;
//	// Get 3 Vector from Y
//	TVector3 Y3  = (Y).Vect(); 
//	// Construct combined pmu/Kmu 4 vector in new rotated frame
//	TLorentzVector Ppmu(Y3.Dot(M_dirn),(Y).Perp(M_dirn),0.0,(Y).E() );
//	// Get component of Y3 perpendicular to *M_dirn
//	TVector3 Perp_dirn  = ( Y3  - (M_dirn) * Y3.Dot((M_dirn))).Unit();  
//	// Calculate neutrino energy in mother rest frame
//	double E_nurest = ( mass  * mass  - (Y).M2())/ (2 * mass ); 
//	// Calculate pmu " "
//	double E_pmurest = sqrt(E_nurest * E_nurest + (Y).M2());
//	double px_rest;
//
//	// Find magnitude of momentum along mother dirn in rest frame
//
//	if( E_nurest * E_nurest  - Ppmu.Py() * Ppmu.Py() >= 0)
//	{
//	    px_rest = sqrt( E_nurest * E_nurest  - Ppmu.Py() * Ppmu.Py() ); 
//	}
//	else
//	{  
//	    return p_nu;
//	}
//
//	//px_rest = sqrt( E_nurest * E_nurest  - Ppmu.Py() * Ppmu.Py() ); 
//	// quadratic coefficients
//	double A = (Y).E() * (Y).E()  - Ppmu.Px() * Ppmu.Px();
//	double B = - 2 * Ppmu.Px()  * (E_nurest * E_pmurest + px_rest * px_rest );
//	double C =  - (E_nurest * E_pmurest + px_rest * px_rest )*
//	    (E_nurest * E_pmurest + px_rest * px_rest ) +  (Y).E() * (Y).E() * Ppmu.Py() * Ppmu.Py();
//
//	if (B*B < 4*A*C) { 
//	    return  p_nu;
//	}
//	// Two neutrino E/p solutions in Lab Frame
//	double p1nu = (- B + sqrt(B*B - 4 * A * C) )/(2*A);
//	double p2nu = (- B - sqrt(B*B - 4 * A * C) )/(2*A);
//
//	// reconstruct neutrino 3 vectors and 4 vectors
//	TVector3 P_nu_recon_3V1 = (M_dirn) * p1nu  + Perp_dirn * -Ppmu.Py(); 
//	TVector3 P_nu_recon_3V2 = (M_dirn) * p2nu  + Perp_dirn * -Ppmu.Py(); 
//	p_nu.push_back(TLorentzVector(P_nu_recon_3V1, sqrt(p1nu*p1nu + Ppmu.Py() * Ppmu.Py())));
//	p_nu.push_back(TLorentzVector(P_nu_recon_3V2, sqrt(p2nu*p2nu + Ppmu.Py() * Ppmu.Py())));
//	return p_nu;
//    }
//
bool LoKi::BDTTag::Svr::mCorErrors(TVector3 v1, TVector3 v2, Gaudi::LorentzVector p,  Gaudi::SymMatrix7x7 cov1, Gaudi::SymMatrix3x3 cov2 )
{
    //std::cout << cov1[0][1] << "  "  << cov1[0][2] << "   "  << cov1[1][2] << std::endl;
   // std::cout << cov1[1][0] << "  "  << cov1[2][0] << "   "  << cov1[2][1] << std::endl;
    double x  = v1.Px();
    double y  = v1.Py();
    double z  = v1.Pz();
    double xp  = v2.Px();
    double yp  = v2.Py();
    double zp  = v2.Pz();
    double px = p.Px();
    double py = p.Py();
    double pz = p.Pz();
    double k = p.E();
//     double P_T = (p.Vect()).Perp((v1 - v2));
  //   double M_corr =sqrt(p.M2() + P_T * P_T) + P_T;

double dMcdpx = (2*(1 - pow(x - xp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
	             (2*(x - xp)*(y - yp)*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/
					               (pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		           (2*(x - xp)*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			             (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (-2*px + 2*(1 - pow(x - xp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
		         (2*(x - xp)*(y - yp)*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/
					                   (pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
			       (2*(x - xp)*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdpy = (2*(1 - pow(y - yp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
	             (2*(x - xp)*(y - yp)*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/
					               (pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		           (2*(y - yp)*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			             (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (-2*py + 2*(1 - pow(y - yp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
		         (2*(x - xp)*(y - yp)*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/
					                   (pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
			       (2*(y - yp)*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdpz = (2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (1 - pow(z - zp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
	             (2*(x - xp)*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		               (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		           (2*(y - yp)*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			             (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (-2*pz + 2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (1 - pow(z - zp,2)/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
		         (2*(x - xp)*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			           (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
			       (2*(y - yp)*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdE = k/sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
	     pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
	          pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		       pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2));

double dMcdx = (2*((2*pow(x - xp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
	             (px*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		              (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
	             2*((2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
			          (px*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
			          2*((2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
				               (px*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*((2*pow(x - xp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
	                (px*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
			         (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		         2*((2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
			              (px*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
				      2*((2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
					           (px*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdy = (2*((2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
	             (py*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
	             2*((2*pow(y - yp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
			          (py*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
				           (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
			          2*((2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
				               (py*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*((2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
	                (py*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		         2*((2*pow(y - yp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
			              (py*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
				               (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
				      2*((2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) - 
					           (py*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdz = (2*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (-((pz*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		         (2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
	             2*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            (-((pz*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
			              (2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
			          2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (-((px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
					           (pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
						            (2*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*pow(z - zp,2))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (-((pz*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		             (2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
		         2*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        (-((pz*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
				          (2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
				      2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             (-((px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) - 
					               (pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
						                (2*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*pow(z - zp,2))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdxp = (2*((-2*pow(x - xp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
	             (px*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
		              (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
	             2*((-2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
			          (px*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
			          2*((-2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
				               (px*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*((-2*pow(x - xp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
	                (px*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
			         (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		         2*((-2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
			              (px*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
				      2*((-2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
					           (px*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdyp = (2*((-2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
	             (py*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
	             2*((-2*pow(y - yp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
			          (py*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
				           (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
			          2*((-2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
				               (py*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*((-2*(x - xp)*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
	                (py*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           (px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
		         2*((-2*pow(y - yp,2)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
			              (py*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
				               (px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        (py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))) + 
				      2*((-2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2) + 
					           (py*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             (pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2))))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double dMcdzp = (2*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	       ((pz*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		         (2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
	             2*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
		            ((pz*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
			              (2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
			          2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				         ((px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
					           (pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
						            (2*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*pow(z - zp,2))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)))/
        (2.*sqrt(pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
		         pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
			         pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2))) + 
	   (2*(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
	           ((pz*(x - xp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
		             (2*(x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
		         2*(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
			        ((pz*(y - yp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
				          (2*(y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)) + 
				      2*(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)))*
				             ((px*(x - xp) + py*(y - yp) + pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) + 
					               (pz*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)) - 
						                (2*(px*(x - xp) + py*(y - yp) + pz*(z - zp))*pow(z - zp,2))/pow(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2),2)))/
	       (2.*sqrt(pow(k,2) - pow(px,2) - pow(py,2) - pow(pz,2) + 
			        pow(px - ((x - xp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
				        pow(py - ((y - yp)*(px*(x - xp) + py*(y - yp) + pz*(z - zp)))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2) + 
					        pow(pz - ((px*(x - xp) + py*(y - yp) + pz*(z - zp))*(z - zp))/(pow(x - xp,2) + pow(y - yp,2) + pow(z - zp,2)),2)));

double vertex_errsq = 
	cov1(0,0) * dMcdx * dMcdx + cov1(1,1) * dMcdy * dMcdy  + cov1(2,2) * dMcdz * dMcdz + 
	cov2(0,0) * dMcdxp * dMcdxp + cov2(1,1) * dMcdyp * dMcdyp  + cov2(2,2) * dMcdzp * dMcdzp  + 2.* dMcdx *dMcdy * cov1(0,1)  + 2. * dMcdx * dMcdz * cov1(0,2) + 2. * dMcdy * dMcdz * cov1(1,2)  + 2.* dMcdxp *dMcdyp * cov2(0,1)  + 2. * dMcdxp * dMcdzp * cov2(0,2) + 2. * dMcdyp * dMcdzp * cov2(1,2);


double momentum_errsq =
	  dMcdpx * dMcdpx * cov1(3,3) +   dMcdpy * dMcdpy * cov1(4,4) + dMcdpz * dMcdpz * cov1(5,5) + dMcdE * dMcdE * cov1(6,6) + 
        // mom v mom cross terms
	+  2. * dMcdpx * dMcdpy * cov1(3,4) +   2.* dMcdpx * dMcdpz * cov1(3,5) + 2. * dMcdpx * dMcdE * cov1(3,6) + 2. * dMcdpy * dMcdpz * cov1(4,5) + 2. * dMcdpy * dMcdE * cov1(4,6)         + 2. *  dMcdpz * dMcdE * cov1(5,6) +  
	// mom vs positon terms
	2 * dMcdx * dMcdpx *  cov1(0,3) + 2. * dMcdy * dMcdpx *  cov1(1,3)  + 2. * dMcdz * dMcdpx *  cov1(2,3) +
	 2 *dMcdx * dMcdpy *  cov1(0,4) + 2. * dMcdy * dMcdpy *  cov1(1,4)  + 2. * dMcdz * dMcdpy *  cov1(2,4) +
	2*  dMcdx * dMcdpz *  cov1(0,5) + 2. * dMcdy * dMcdpz *  cov1(1,5)  + 2. * dMcdz * dMcdpz *  cov1(2,5) +
	2*  dMcdx * dMcdE *  cov1(0,6) + 2. * dMcdy * dMcdE *  cov1(1,6)  + 2. * dMcdz * dMcdE *  cov1(2,6) ;

m_mCorErr = sqrt(vertex_errsq);
m_mCorErrFull = sqrt(vertex_errsq + momentum_errsq);

return true;
}
