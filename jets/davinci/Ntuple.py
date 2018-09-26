# Access to classes.
from collections import OrderedDict
import ROOT, array, GaudiPython
from GaudiPython.Bindings import gbl
STD  = gbl.std
LHCB = gbl.LHCb

# Simple ntuple class.
class Ntuple:
    """
    Class to store an ntuple.
    """
    def __init__(self, name, tes, toolSvc, detSvc, jetPath, recSvrPath, recMuPath, hasL0=True):
        """
        Initialize the ntuple with the needed tools.
        """
        self.hasL0 = hasL0
        ROOT.gInterpreter.ProcessLine(
            "#include \"/cvmfs/lhcb.cern.ch/lib/lhcb/HLT/"
            "HLT_v25r4/Hlt/HltDisplVertices/Kernel/"
            "IMatterVeto.h\"")
        self.pvrTool = toolSvc.create(
            'GenericParticle2PVRelator<_p2PVWithIPChi2, '
            'OfflineDistanceCalculatorName>/P2PVWithIPChi2',
            interface = 'IRelatedPVFinder')
        self.tagTool = toolSvc.create(
            'LoKi::BDTTag',
            interface = 'IJetTagTool')
        self.genTool = toolSvc.create(
            'DaVinciSmartAssociator',
            interface = 'IParticle2MCWeightedAssociator')
        self.mtrTool = toolSvc.create(
            'MatterVetoTool',
            interface = 'IMatterVeto')
        self.velTool = toolSvc.create(
            'VeloExpectation',
            interface = 'IVeloExpectation')
        self.dstTool = toolSvc.create(
            'LoKi::TrgDistanceCalculator',
            interface = 'IDistanceCalculator')
        self.ltTool = toolSvc.create(
            'LoKi::LifetimeFitter',
            interface = 'ILifetimeFitter')
        self.trkTool = toolSvc.create(
            'TrackMasterExtrapolator',
            interface = 'ITrackExtrapolator')
        self.pidTool = toolSvc.create(
            'ANNGlobalPID::ChargedProtoANNPIDTool',
            interface = 'ANNGlobalPID::IChargedProtoANNPIDTool')
        self.detTool = detSvc[
            '/dd/Structure/LHCb/BeforeMagnetRegion/Velo']
        self.l0Tool = toolSvc.create(
                 'L0TriggerTisTos',
                 interface = 'ITriggerTisTos')
        self.hlt1Tool = toolSvc.create(
                 'TriggerTisTos/Hlt1TriggerTisTos',
                 interface = 'ITriggerTisTos')
        self.hlt2Tool = toolSvc.create(
                 'TriggerTisTos/Hlt2TriggerTisTos',
                 interface = 'ITriggerTisTos')

        self.stable = [11,-11,13,-13,211,-211,321,-321,2212,-2212,2112,-2112,22,111,310,130,311,-311]
        self.Ds = [411,-411,421,-421,431,-431,4122,-4122]
        self.tes     = tes
        self.jetPath = jetPath
        self.recSvrPath = recSvrPath
        self.recMuPath = recMuPath
        self.saved   = {}
        self.ntuple  = OrderedDict()
        self.tfile   = ROOT.TFile('output.root', 'RECREATE')
        self.ttree   = ROOT.TTree('data', 'data')
        self.vrs     = {}
        mom = ['p','pt','px', 'py', 'pz', 'e']
        pos = ['x', 'y', 'z']
        cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
        l0trig = ['%s_%s' % (line,dec) for line in ['l0_hadron','l0_photon','l0_electron','l0_muon','l0_dimuon'] for dec in ['dec','tis','tos']]
        hlt1trig = ["%s_%s" % (line,dec) for line in ['hlt1_track','hlt1_ditrack'] for dec in ['dec','tis','tos']]
        self.init('gen', ['idx_pvr', 'idx_jet', 'idx_prnt', 'pid', 'q'] + mom + pos + ['prnt_pid', 'res_pid', 'from_sig'])
        self.init('pvr', pos + cov)
        self.init('svr', ['idx_pvr', 'idx_jet'] + [
                'idx_trk%i' % i for i in range(0, 10)] + 
                  mom + pos + ['dx', 'dy', 'dz', 'm', 'm_cor', 'm_cor_err', 'm_cor_err_full', 'fd_min', 'fd_chi2', 'chi2', 'ip_chi2', 'ip_chi2_sum', 'ip_chi2_min_trk', 'abs_q_sum', 'tau', 'ntrk', 'ntrk_jet', 'jet_dr', 'jet_pt', 'pass', 'bdt0', 'bdt1', 'in_mtr', 'backwards', 'nTBVs'] + l0trig + hlt1trig)
        self.init('jet', ['idx_pvr', 'ntrk', 'nneu'] + mom + l0trig + hlt1trig)
        self.init('trk', ['idx_gen', 'idx_pvr', 'idx_jet'] + mom +
                  ['pid', 'q', 'ip', 'ip_chi2', 'pnn_e', 'pnn_mu', 'pnn_pi',
                   'pnn_k', 'pnn_p', 'pnn_ghost', 'ecal', 'hcal', 'prb_ghost', 'type', 'is_mu',
                   'vid', 'x', 'y', 'z'] + l0trig + hlt1trig)# + ['vhit%i' % i for i in range(0, 61)])# + [
                   #'vz%i' % i for i in range(0, 61)])
        self.init('neu', ['idx_gen', 'idx_jet'] + mom + ['pid'] + l0trig + hlt1trig)
        self.init('z0', ['idx_pvr','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_dr'] + mom + pos + ['m', 'ip', 'ip_chi2','dr_jet'] + ['idx_trk%i' % i for i in range(0, 2)] + l0trig + hlt1trig)
        self.init('jpsi', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 2)] + l0trig + hlt1trig)
        self.init('d0', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 2)] + l0trig + hlt1trig)
        self.init('dp', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_trk2', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)] + l0trig + hlt1trig)
        self.init('ds', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_trk2', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)] + l0trig + hlt1trig)
        self.init('lc', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_trk2', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)] + l0trig + hlt1trig)
        self.init('k3pi', ['idx_pvr','idx_jet','idx_jet_trk0', 'idx_jet_trk1', 'idx_jet_trk2', 'idx_jet_trk3', 'idx_jet_dr','dr_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 4)] + l0trig + hlt1trig)
        self.init('evt', ['dec'] + ['%s_%s' % (k,l) for k in ['j1', 'j2'] for l in ['idx','dR','nsv','nmu','ntrk','nneu'] + mom ] )
        self.ntuple['evt_pvr_n'] = array.array('d', [-1])
        self.ntuple['evt_trk_n'] = array.array('d', [-1])
        for key, val in self.ntuple.iteritems():
            if type(val) is array.array: self.ttree.Branch(key, val, key + '/D')
            else: self.ttree.Branch(key, val)
    def init(self, pre, vrs):
        """
        Initialize a set of variables.
        """
        self.saved[pre] = {}
        self.vrs[pre]   = vrs
        for v in vrs: self.ntuple['%s_%s' % (pre, v)] = ROOT.vector('double')()
    def key(self, obj):
        """
        Generate the key for an object.
        """
        key = None
        try: 
            key = (obj.momentum().Px(), obj.momentum().Py(), 
                   obj.momentum().Pz())
            try:
                trk = obj.proto().track()
                key = (trk.momentum().X(), trk.momentum().Y(), 
                       trk.momentum().Z())
            except:
                try:
                    pos = obj.proto().calo()[0].position()
                    key = (pos.x(), pos.y(), pos.z(), pos.e())
                except: pass
        except:
            try: key = (obj.position().X(), obj.position().Y(),
                        obj.position().Z())
            except: pass
        return key
    def close(self):
        """
        Close the ntuple.
        """
        self.tfile.Write(); self.tfile.Close()
    def clear(self):
        """
        Clear the ntuple.
        """
        for key, val in self.saved.iteritems(): val.clear()
        for key, val in self.ntuple.iteritems():
            if type(val) is array.array: val[0] = -1
            else: val.clear()
    def lookupVeloStation(self, z):
        return {   4 :  0,   6 :  1,  19 :  2,  21 :  3,  34 :  4,  36 :  5,
                  49 :  6,  51 :  7,  64 :  8,  66 :  9,  79 : 10,  81 : 11,
                  93 : 12,  96 : 13, 109 : 14, 111 : 15, 123 : 16, 126 : 17,
                 139 : 18, 141 : 19, 153 : 20, 156 : 21, 169 : 22, 171 : 23,
                 184 : 24, 186 : 25, 198 : 26, 201 : 27, 214 : 28, 216 : 29,
                 229 : 30, 231 : 31, 244 : 32, 246 : 33, 258 : 34, 261 : 35,
                 273 : 36, 276 : 37, 289 : 38, 291 : 39, 434 : 40, 436 : 41,
                 449 : 42, 451 : 43, 584 : 44, 586 : 45, 599 : 46, 601 : 47,
                 634 : 48, 636 : 49, 648 : 50, 651 : 51, 684 : 52, 686 : 53,
                 699 : 54, 701 : 55, 734 : 56, 736 : 57, 748 : 58, 751 : 59
               }.get(int(floor(z)),-1)
    def getJetMult(self, summary, nNeu=0, nChrg=0):
        if summary.summarizedObjectCLID() == 801:
            if 801 in [ss.summarizedObjectCLID() for ss in summary.substructure()]:##do any of the substructures contain a particle class
                for ss in summary.substructure():
                    nNeu,nChrg = self.getJetMult(ss, nNeu, nChrg)
            else:
                if summary.numericalInfo()["0#Particle.particleID.pid"] in [22,111]:
                    nNeu+=1
                else:
                    nChrg+=1
        return nNeu,nChrg
    def processHltSummary(self, summary, indent=""):
        if summary.summarizedObjectCLID() == 801:
            print indent+str(summary.numericalInfo()["0#Particle.particleID.pid"])
        if not summary.substructure().empty():
            for ss in summary.substructure():
                self.processHltSummary(ss,indent+" ")

    def children(self, gen):
        parts = []
        for vtx in gen.endVertices():
            for part in vtx.products():
                parts += [part]
        return parts

 # Return the final charged children of a generated particle.
    def finalChildren(self, gen):
        inters = []
        finals = []
        if not gen: return []
        final = []
        decay = self.children(gen)
        while len(decay) > 0:
            tmp = self.children(decay[-1])
            pid = int(decay[-1].particleID().abspid())
            if (len(tmp) == 0 or pid in self.stable):
                if pid != 22: finals += [pid]
                final += [decay[-1]]
                decay.pop()
            else:
                if pid != 22: inters += [pid]
                decay.pop()
                decay += tmp
        if gen.particleID().abspid() == 15:
            while finals[-1] != 16:
                final.pop(); finals.pop()
        finals.sort()
        inters.sort()
        return (final, inters, finals)

    def fill(self, key = None, val = None, idx = None, vrs = None):
        """
        Fill the ntuple for either an event or an object.
        """
        if key == None and val == None and idx == None and vrs == None:
            self.tfile.Cd(''); self.ttree.Fill()
        elif vrs != None and key != None:
            pre = key
            for key in self.vrs[pre]:
                val = vrs[key] if key in vrs else -1
                self.ntuple[pre + '_' + key].push_back(val)
        elif key in self.ntuple: 
            if idx == None: self.ntuple[key].push_back(val)
            elif idx < len(self.ntuple[key]): self.ntuple[key][idx] = val
    def fillMom(self, obj, vrs):
        if not obj: return
        vrs['p'] = obj.P()
        vrs['pt'] = obj.Pt()
        vrs['px'] = obj.Px()
        vrs['py'] = obj.Py()
        vrs['pz'] = obj.Pz()
        vrs['e' ] = obj.E()
    def fillPid(self, obj, vrs):
        if not obj: return
        vrs['pid'] = obj.pid()
        vrs['q'  ] = float(obj.threeCharge())/3.0
    def fillPro(self, obj, vrs):
        if not obj: return
        if obj.muonPID(): vrs['is_mu'] = obj.muonPID().IsMuon()
        vrs['pnn_e' ] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(11),   "MC15TuneV1").value #obj.info(700, -100)
        vrs['pnn_mu'] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(13),   "MC15TuneV1").value #obj.info(701, -100)
        vrs['pnn_pi'] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(211),  "MC15TuneV1").value #obj.info(702, -100)
        vrs['pnn_k' ] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(321),  "MC15TuneV1").value #obj.info(703, -100)
        vrs['pnn_p' ] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(2212), "MC15TuneV1").value #obj.info(704, -100)
        vrs['pnn_ghost' ] = self.pidTool.annPID(obj, ROOT.LHCb.ParticleID(0),"MC15TuneV1").value
        vrs['ecal']   = obj.info(332, -100)
        vrs['hcal']   = obj.info(333, -100)
    def fillTrk(self, obj, pid, vrs):
        if not obj or not pid: return
        vrs['prb_ghost'] = obj.ghostProbability()
        vrs['type'] = obj.type()
        vrs['vid'] = 1.0
        vids = []
        for vid in obj.lhcbIDs():
            if vid.isVelo():
                vrs['vid'] *= float(vid.veloID().channelID())/1000000.0
                vids += [(self.detTool.sensor(vid.veloID()).z(), 
                          vid.veloID().channelID())]
                #station = self.lookupVeloStation(vids[-1][0])
                #vrs['vz%i' % (station)] = vids[-1][0]
                #vrs['vhit%i' % (station)] = vids[-1][1]
        vids.sort()
        sta = LHCB.StateVector()
        if len(vids) > 0: self.trkTool.propagate(obj, vids[0][0], sta, pid)
        vrs['x'] = sta.x()
        vrs['y'] = sta.y()
        vrs['z'] = sta.z()
    def fillDst(self, obj, pvr, dst, vrs):
        if not obj or not pvr: return
        val, valChi2 = ROOT.Double(-1), ROOT.Double(-1)
        self.dstTool.distance(obj, pvr, val, valChi2)
        vrs[dst] = val;
        vrs[dst + '_chi2'] = valChi2;
    def fillLifetime(self, obj, pvr, vrs):
        if not obj or not pvr: return
        val, valErr, valChi2 = ROOT.Double(-1), ROOT.Double(-1), ROOT.Double(-1)
        self.ltTool.fit(pvr, obj, val, valErr, valChi2)
        vrs['tau'] = val;
        vrs['tau_err'] = valErr;
        vrs['tau_chi2'] = valChi2;
    def fillVtx(self, obj, vrs):
        vrs['vtx_chi2'] = obj.chi2()
        vrs['vtx_ndof'] = obj.nDoF()
    def fillPos(self, obj, vrs):
        if not obj: return
        pos = obj.position()
        vrs['x'] = pos.X()
        vrs['y'] = pos.Y()
        vrs['z'] = pos.Z()
        vrs['in_mtr'] = self.mtrTool.isInMatter(pos)
    def fillCov(self, obj, vrs):
        try:
            cov = obj.covMatrix()
            vrs['ndof'] = obj.nDoF()
            vrs['chi2'] = obj.chi2()
            vrs['dx']   = cov[0][0]
            vrs['dy']   = cov[1][1]
            vrs['dz']   = cov[2][2]
        except: pass
    def fillPvr(self, obj, vrs):
        if not obj: return
        vrs['idx_pvr'] = self.addPvr(obj)
    def fillGen(self, obj, vrs):
        if not obj: return
        gen = None; wgt = 0; rels = self.genTool.relatedMCPs(obj)
        for rel in rels: gen = rel.to() if rel.weight() > wgt else gen
        if gen: vrs['idx_gen'] = self.addGen(gen) 
    def fillTrig(self, obj, vrs):
        if not obj: return
        try:
            self.hlt1Tool.setOfflineInput()
            self.hlt1Tool.addToOfflineInput(obj)
            dec = self.hlt1Tool.triggerTisTos(obj,"Hlt1TrackMVADecision")
            vrs['hlt1_track_dec'] = dec.decision()
            vrs['hlt1_track_tis'] = dec.tis()
            vrs['hlt1_track_tos'] = dec.tos()
            dec = self.hlt1Tool.triggerTisTos(obj,"Hlt1TwoTrackMVADecision")
            vrs['hlt1_ditrack_dec'] = dec.decision()
            vrs['hlt1_ditrack_tis'] = dec.tis()
            vrs['hlt1_ditrack_tos'] = dec.tos()
        except:
            pass
        if not self.hasL0: return
        try:
            self.l0Tool.setOfflineInput()
            self.l0Tool.addToOfflineInput(obj)
            dec = self.l0Tool.triggerTisTos("L0HadronDecision")
            vrs['l0_hadron_dec'] = dec.decision()
            vrs['l0_hadron_tis'] = dec.tis()
            vrs['l0_hadron_tos'] = dec.tos()
            dec = self.l0Tool.triggerTisTos("L0PhotonDecision")
            vrs['l0_photon_dec'] = dec.decision()
            vrs['l0_photon_tis'] = dec.tis()
            vrs['l0_photon_tos'] = dec.tos()
            dec = self.l0Tool.triggerTisTos("L0ElectronDecision")
            vrs['l0_electron_dec'] = dec.decision()
            vrs['l0_electron_tis'] = dec.tis()
            vrs['l0_electron_tos'] = dec.tos()
            dec = self.l0Tool.triggerTisTos("L0MuonDecision")
            vrs['l0_muon_dec'] = dec.decision()
            vrs['l0_muon_tis'] = dec.tis()
            vrs['l0_muon_tos'] = dec.tos()
            dec = self.l0Tool.triggerTisTos("L0DiMuonDecision")
            vrs['l0_dimuon_dec'] = dec.decision()
            vrs['l0_dimuon_tis'] = dec.tis()
            vrs['l0_dimuon_tos'] = dec.tos()
        except:
            pass

    def addZ(self, obj, pre="z0"):
        vrs = {}
        trks = []
        try:
            jets = self.tes[self.jetPath]
            trkJets = { 0:-1, 1:-1}
            bestJet = -1
            bestDr = 10.0
            for idx, jet in enumerate(jets):
                jetTrkKeys = [ self.key(dau) for dau in jet.daughters() ]
                nIn=0
                for idau, dau in enumerate(obj.daughters()):
                    key = self.key(dau)
                    if key in jetTrkKeys:
                        nIn+=1
                        trkJets[idau] = idx
                if nIn==0:
                    dr = ROOT.TLorentzVector(jet.momentum().Px(),jet.momentum().Py(),jet.momentum().Pz(),0.).DeltaR(ROOT.TLorentzVector(obj.momentum().Px(),obj.momentum().Py(),obj.momentum().Pz(),0.))
                    if bestDr > dr:
                        bestDr = dr
                        bestJet = idx
            
            for dau in obj.daughters():
                trks.append(self.addTrk(dau))
            pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
            self.fillPvr(pvr, vrs)
            self.fillDst(obj, pvr, 'ip', vrs)
            ##self.fillDst(obj.endVertex(), pvr, 'fd', vrs)
            ##self.fillLifetime(obj, pvr, vrs)
            self.fillMom(obj.momentum(), vrs)
            ##self.fillPos(obj.endVertex(), vrs)
            ##self.fillVtx(obj.vertex(), vrs)
            self.fillTrig(obj, vrs)
            vrs['m'] = obj.momentum().M()
            vrs['idx_jet_trk0'] = trkJets[0]
            vrs['idx_jet_trk1'] = trkJets[1]
            vrs['idx_jet_dr'] = bestJet
            vrs['dr_jet'] = bestDr
            for idx, trk in enumerate(trks):
                vrs['idx_trk%i' % idx] = trk

            self.fill(pre, vrs = vrs)

            if bestDr < 10.0:
                return True

        except:
            pass

        return False
    def addDHad(self, obj, pre="d0"):
        vrs = {}
        trks = []
        try:
            jets = self.tes[self.jetPath]
            nInBest = 0
            bestJet = -1
            bestJetDr = -1
            bestDr = 10.0
            for idx, jet in enumerate(jets):
                jetTrkKeys = [ self.key(dau) for dau in jet.daughters() ]
                nIn=0
                for dau in obj.daughters():
                    key = self.key(dau)
                    if key in jetTrkKeys:
                        nIn+=1
                if nIn>nInBest:
                    nInBest=nIn
                    bestJet=idx
                dr = ROOT.TLorentzVector(jet.momentum().Px(),jet.momentum().Py(),jet.momentum().Pz(),0.).DeltaR(ROOT.TLorentzVector(obj.momentum().Px(),obj.momentum().Py(),obj.momentum().Pz(),0.))
                if bestDr > dr:
                    bestDr = dr
                    bestJetDr = idx
            
            for dau in obj.daughters():
                trks.append(self.addTrk(dau))
            pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
            self.fillPvr(pvr, vrs)
            self.fillDst(obj, pvr, 'ip', vrs)
            self.fillDst(obj.endVertex(), pvr, 'fd', vrs)
            self.fillLifetime(obj, pvr, vrs)
            self.fillMom(obj.momentum(), vrs)
            self.fillPos(obj.endVertex(), vrs)
            self.fillVtx(obj.vertex(), vrs)
            self.fillTrig(obj, vrs)
            vrs['m'] = obj.momentum().M()
            vrs['idx_jet'] = bestJet
            vrs['ntrk_jet'] = nInBest
            vrs['idx_jet_dr'] = bestJetDr
            vrs['dr_jet'] = bestDr
            for idx, trk in enumerate(trks):
                vrs['idx_trk%i' % idx] = trk

            self.fill(pre, vrs = vrs)
        except:
            pass

    def addTrigger(self, pre = 'evt'):
        for idec, dec in enumerate(['DiJet','DiJetSV','DiJetSVSV','DiJetSVMu','DiJetMuMu','DiJetLowPt','DiJetSVLowPt','DiJetSVSVLowPt','DiJetSVMuLowPt','DiJetMuMuLowPt','DiJetHighPt','DiJetSVHighPt','DiJetSVSVHighPt','DiJetSVMuHighPt','DiJetMuMuHighPt']):
            vrs = {}
            summaries = self.hlt2Tool.hltObjectSummaries('Hlt2Jets'+dec+'Decision')
            if summaries.size()>0:
                vrs['dec'] = idec
            else:
                continue
            for summary in summaries:
                if summary.numericalInfo()["0#Particle.particleID.pid"] != 97: return
                jet1 = summary.substructure().at(0)
                jet2 = summary.substructure().at(1)
                if jet1.numericalInfo()["0#Particle.particleID.pid"] != 98 or jet2.numericalInfo()["0#Particle.particleID.pid"] != 98: return
                j1tx = jet1.numericalInfo()["5#Particle.slopes.x"]
                j1ty = jet1.numericalInfo()["6#Particle.slopes.y"]
                j1p  = 1./jet1.numericalInfo()["7#Particle.1/p"]
                j1m  = jet1.numericalInfo()["1#Particle.measuredMass"]
                j1Denom = (1. + j1tx**2. + j1ty**2.)**.5 
                j1px = j1p * j1tx / j1Denom
                j1py = j1p * j1ty / j1Denom
                j1pz = j1p  / j1Denom
                j1pt = j1p * (j1tx**2. + j1ty**2.)**.5 / j1Denom
                j1e  = (j1p**2. + j1m**2.)**.5
                j2tx = jet2.numericalInfo()["5#Particle.slopes.x"]
                j2ty = jet2.numericalInfo()["6#Particle.slopes.y"]
                j2p  = 1./jet2.numericalInfo()["7#Particle.1/p"]
                j2m  = jet2.numericalInfo()["1#Particle.measuredMass"]
                j2Denom = (1. + j2tx**2. + j2ty**2.)**.5 
                j2px = j2p * j2tx / j2Denom
                j2py = j2p * j2ty / j2Denom
                j2pz = j2p  / j2Denom
                j2pt = j2p * (j2tx**2. + j2ty**2.)**.5 / j2Denom
                j2e  = (j2p**2. + j2m**2.)**.5
                nNeu1=0
                nChrg1=0
                nNeu1,nChrg1 = self.getJetMult(jet1)
                nNeu2,nChrg2 = self.getJetMult(jet2)

                vrs['j1_ntrk'] = nChrg1
                vrs['j1_nneu'] = nNeu1
                vrs['j1_px'] = j1px
                vrs['j1_py'] = j1py
                vrs['j1_pz'] = j1pz
                vrs['j1_e']  = j1e

                vrs['j2_ntrk'] = nChrg2
                vrs['j2_nneu'] = nNeu2
                vrs['j2_px'] = j2px
                vrs['j2_py'] = j2py
                vrs['j2_pz'] = j2pz
                vrs['j2_e']  = j2e

                j1_idx=-1
                j2_idx=-1
                j1_nsv=0
                j2_nsv=0
                j1_nmu=0
                j2_nmu=0
                j1_dr=0.5
                j2_dr=0.5
                
                try:
                    jets = self.tes[self.jetPath]
                    for idx, jet in enumerate(jets):
                        dr = ROOT.TLorentzVector(jet.momentum().px(),jet.momentum().py(),jet.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j1px,j1py,j1pz,0.))
                        if j1_dr > dr:
                            j1_dr = dr
                            j1_idx=idx
                        dr = ROOT.TLorentzVector(jet.momentum().px(),jet.momentum().py(),jet.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j2px,j2py,j2pz,0.))
                        if j2_dr > dr:
                            j2_dr = dr
                            j2_idx=idx
                except:
                    pass
                try:
                    svs = self.tes[self.recSvrPath]
                    for sv in svs:
                        pv = self.pvrTool.relatedPV(sv, 'Rec/Vertex/Primary')
                        fd = ROOT.TLorentzVector(sv.vertex().position().X() - pv.position().X(), sv.vertex().position().Y() - pv.position().Y(), sv.vertex().position().Z() - pv.position().Z(), 0.)
                        if 0.5 > fd.DeltaR(ROOT.TLorentzVector(j1px,j1py,j1pz,0.)):
                            j1_nsv+=1
                        if 0.5 > fd.DeltaR(ROOT.TLorentzVector(j2px,j2py,j2pz,0.)):
                            j2_nsv+=1
                except:
                    pass
                try:
                    mus = self.tes[self.recMuPath]
                    for mu in mus:
                        if 0.5 > ROOT.TLorentzVector(mu.momentum().px(),mu.momentum().py(),mu.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j1px,j1py,j1pz,0.)):
                            j1_nmu+=1
                        if 0.5 > ROOT.TLorentzVector(mu.momentum().px(),mu.momentum().py(),mu.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j2px,j2py,j2pz,0.)):
                            j2_nmu+=1
                except:
                    pass

                if idec==2 and j1_idx>=0 and j2_idx>=0 and (j1_nsv==0 or j2_nsv==0):
                    pass

                vrs['j1_idx'] = j1_idx
                vrs['j2_idx'] = j2_idx
                vrs['j1_dR'] = j1_dr
                vrs['j2_dR'] = j2_dr
                vrs['j1_nsv'] = j1_nsv
                vrs['j2_nsv'] = j2_nsv
                vrs['j1_nmu'] = j1_nmu
                vrs['j2_nmu'] = j2_nmu

                self.fill(pre, vrs = vrs)
           
    def addGen(self, obj, jet = -1, pre = 'gen', par = None):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        parent = obj.mother()
        res = None
        if par and parent.particleID().pid()!=par.particleID().pid(): #if they don't match we have a resonance
            res = parent
            parent = par
        parKey = -1
        parIdx = -1
        fromSig=0
        if obj.fromSignal(): fromSig = 1
        vrs = {}
        idx = len(self.saved[pre])
        self.fillPid(obj.particleID(), vrs)
        self.fillMom(obj.momentum(), vrs)
        self.fillPos(obj.originVertex(), vrs)
        self.fillPvr(obj.primaryVertex(), vrs)
        vrs['idx_jet'] = jet
        vrs['from_sig'] = fromSig
        if res:
            vrs['res_pid'] = res.particleID().pid()
        if parent:
            parKey = self.key(parent)
            if parKey in self.saved[pre]: parIdx = self.saved[pre][parKey]
            vrs['idx_prnt'] = parIdx
            vrs['prnt_pid'] = parent.particleID().pid()
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
        #also print immediate decay products of heavy hadrons
        pid = obj.particleID()
        if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
            if pid.abspid() in self.Ds:
                for part in self.finalChildren(obj)[0]:
                    self.addGen(part,par=obj)
            else:
                for part in self.children(obj):
                    self.addGen(part)
        if pid.pid() == 23:
            for part in self.children(obj):
                self.addGen(part)
            #for vtx in obj.endVertices():
            #    for part in vtx.products():
            #        self.addGen(part)
        return idx
    def addPvr(self, obj, pre = 'pvr'):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        self.fillPos(obj, vrs)
        self.fillCov(obj, vrs)
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
        return idx
#    def addSvr(self, obj, pre = 'svr'):
#        vrs = {}
#        pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
#        self.fillMom(obj.momentum(), vrs)
#        self.fillPos(obj.endVertex(), vrs)
#        self.fillCov(obj.endVertex(), vrs)
#        self.fillPvr(pvr, vrs)
#        self.fillDst(obj.endVertex(), pvr, 'fd', vrs)
#        trks = []
#        for dtr in obj.daughters():
#            if not dtr.proto() or not dtr.proto().track(): continue
#            trks += [[self.addTrk(dtr), dtr.proto().track()]]
#            vrs['idx_trk%i' % (len(trks) - 1)] = trks[-1][0]
#        self.fill(pre, vrs = vrs)
    def addTags(self, obj, jet = -1, pre = 'svr'):
        tags = STD.map('string', 'double')()
        if not self.tagTool.calculateJetProperty(obj, tags): return
        ntag = int(tags['Tag'])
        for itag in range(0, ntag):
            vrs = {}
            vrs['idx_pvr'] = self.addPvr(self.tes['Rec/Vertex/Primary']
                                         [int(tags['Tag%i_idx_pvr' % itag])])
            vrs['idx_jet'] = jet
            ntrk = int(tags['Tag%i_nTrk' % itag])
            for itrk in range(0, ntrk): vrs['idx_trk%i' % itrk] = self.addTrk(
                self.tes['Phys/StdAllNoPIDsPions/Particles']
                [int(tags['Tag%i_idx_trk%i' % (itag, itrk)])]);
            for vr in ['x', 'y', 'z', 'dx', 'dy', 'dz', 'px', 'py', 'pz', 'e']:
                vrs[vr] = tags['Tag%i_%s' % (itag, vr)]
            vrs['m']  = tags['Tag%i_m' % itag]
            vrs['m_cor']  = tags['Tag%i_mCor' % itag]
            vrs['m_cor_err']  = tags['Tag%i_mCorErr' % itag]
            vrs['m_cor_err_full']  = tags['Tag%i_mCorErrFull' % itag]
            vrs['pt']  = tags['Tag%i_pt' % itag]
            vrs['fd_min'] = tags['Tag%i_fdrMin' % itag]
            vrs['fd_chi2']  = tags['Tag%i_fdChi2' % itag]
            vrs['chi2']  = tags['Tag%i_chi2' % itag]
            vrs['ip_chi2']  = tags['Tag%i_ipChi2Min' % itag]
            vrs['ip_chi2_sum']  = tags['Tag%i_ipChi2Sum' % itag]
            vrs['ip_chi2_min_trk']  = tags['Tag%i_ipChi2MinTrk' % itag]
            vrs['abs_q_sum']  = tags['Tag%i_absQSum' % itag]
            vrs['tau']  = tags['Tag%i_tau' % itag]
            vrs['ntrk']  = tags['Tag%i_nTrk' % itag]
            vrs['ntrk_jet']  = tags['Tag%i_nTrkJet' % itag]
            vrs['jet_dr']  = tags['Tag%i_drSvrJet' % itag]
            vrs['jet_pt']  = tags['Tag%i_ptSvrJet' % itag]
            vrs['pass']  = tags['Tag%i_pass' % itag]
            vrs['bdt0']  = tags['Tag%i_bdt0' % itag]
            vrs['bdt1']  = tags['Tag%i_bdt1' % itag]
            vrs['backwards']  = tags['Tag%i_backwards' % itag]
            vrs['nTBVs']  = tags['Tag%i_nTBVs' % itag]
            self.fillTrig(obj, vrs)
            self.fill(pre, vrs = vrs)
    def addJet(self, obj, pre = 'jet'):
        vrs = {}
        idx = self.ntuple['jet_idx_pvr'].size();
        pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
        self.fillMom(obj.momentum(), vrs)
        self.fillPvr(pvr, vrs)
        self.fillTrig(obj, vrs)
        trks = []
        nneu=0
        for dtr in obj.daughters():
            try:
                #TODO honestly no idea what's going on here
                # if we add recSVs and recMus to the main sequence then obj.daughters contains Particles
                # if we don't, it contains smart poniters to Particles
                # try to treat them as pointers and if that fails assume they're particles
                dtr = dtr.target()
            except:
                pass
            if not dtr.proto() or not dtr.proto().track():
                self.addNeu(dtr, idx)
                nneu+=1
            else:
                trks += [[self.addTrk(dtr, idx), dtr.proto().track()]]
        vrs['ntrk'] = len(trks)
        vrs['nneu'] = nneu
        self.addTags(obj, idx)
        self.fill(pre, vrs = vrs)
    def addNeu(self, obj, jet = -1, pre = 'neu'):
        vrs = {}
        self.fillMom(obj.momentum(), vrs)
        self.fillPid(obj.particleID(), vrs)
        self.fillGen(obj, vrs)
        self.fillTrig(obj, vrs)
        vrs['idx_jet'] = jet
        self.fill(pre, vrs = vrs)
    def addTrk(self, obj, jet = -1, pre = 'trk'):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
        self.fillMom(obj.momentum(), vrs)
        self.fillPid(obj.particleID(), vrs)
        self.fillPro(obj.proto(), vrs)
        self.fillDst(obj, pvr, 'ip', vrs)
        self.fillTrk(obj.proto().track(), obj.particleID(), vrs)
        self.fillPvr(pvr, vrs)
        self.fillGen(obj, vrs)
        self.fillTrig(obj, vrs)
        vrs['idx_jet'] = jet
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
        return idx

