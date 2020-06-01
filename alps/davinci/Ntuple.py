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
    def __init__(self, name, tes, toolSvc, detSvc, hasL0=True, primaries='Rec/Vertex/Primary'):
        """
        Initialize the ntuple with the needed tools.
        """
        self.hasL0 = hasL0
        self.primaries = primaries
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
        self.saved   = {}
        self.ntuple  = OrderedDict()
        self.tfile   = ROOT.TFile('output.root', 'RECREATE')
        self.ttree   = ROOT.TTree('data', 'data')
        self.vrs     = {}
        mom = ['p','pt','px', 'py', 'pz', 'e']
        pos = ['x', 'y', 'z']
        cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
        l0trig = ['%s_%s' % (line,dec) for line in ['l0_hadron','l0_photon','l0_electron','l0_muon','l0_dimuon'] for dec in ['dec','tis','tos']]
        hlt1trig = []#["%s_%s" % (line,dec) for line in ['hlt1_multidimu','hlt1_dimu'] for dec in ['dec','tis','tos']]
        hlt2trig = []#["%s_%s" % (line,dec) for line in ['hlt2_quadmu'] for dec in ['dec','tis','tos']]

        #self.init('gen', ['idx_pvr', 'idx_prnt', 'pid', 'q'] + mom + pos + ['prnt_pid', 'res_pid', 'from_sig'])
        self.init('pvr', pos + cov)
        self.init('trk', ['idx_gen', 'idx_pvr'] + mom +
                  ['pid', 'q'] + ['%s%s' % (i,j) for i in ['ip', 'ip_chi2'] for j in ['','_best', '_bvtx']] + ['pnn_e', 'pnn_mu', 'pnn_pi',
                   'pnn_k', 'pnn_p', 'pnn_ghost', 'ecal', 'hcal', 'prb_ghost', 'type', 'is_mu',
                   'vid', 'x', 'y', 'z'] + l0trig + hlt1trig + hlt2trig)
        self.init('neu', ['idx_gen'] + mom + ['pid'] + l0trig + hlt1trig + hlt2trig)
        self.init('b', ['idx_pvr'] + mom + pos + ['m', 'vtx_chi2', 'vtx_ndof', 'ip', 'ip_chi2', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2'] + ['idx_kst', 'idx_k', 'idx_pi', 'idx_x'] + l0trig + hlt1trig + hlt2trig + ['mode'])
        self.init('kst', ['idx_pvr'] + mom + pos + ['m', 'vtx_chi2', 'vtx_ndof'] + ['%s%s' % (i,j) for i in ['ip', 'ip_chi2', 'fd', 'fd_chi2'] for j in ['','_best', '_bvtx']] + ['tau', 'tau_err', 'tau_chi2'] + ['idx_k', 'idx_pi'] + l0trig + hlt1trig + hlt2trig)
        self.init('x', ['idx_pvr'] + mom + pos + ['m', 'vtx_chi2', 'vtx_ndof'] + ['%s%s' % (i,j) for i in ['ip', 'ip_chi2', 'fd', 'fd_chi2'] for j in ['','_best', '_bvtx']] + ['tau', 'tau_err', 'tau_chi2'] + ['idx_%s%i' % (t,i) for t in ['trk','neu'] for i in range(0, 6)] + l0trig + hlt1trig + hlt2trig)
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
    def fillDst(self, obj, pvr, dst, vrs, mod=""):
        if not obj or not pvr: return
        val, valChi2 = ROOT.Double(-1), ROOT.Double(-1)
        self.dstTool.distance(obj, pvr, val, valChi2)
        vrs[dst + mod] = val;
        vrs[dst + '_chi2' + mod] = valChi2;
    def fillLifetime(self, obj, pvr, vrs):
        if not obj or not pvr: return
        val, valErr, valChi2 = ROOT.Double(-1), ROOT.Double(-1), ROOT.Double(-1)
        self.ltTool.fit(pvr, obj, val, valErr, valChi2)
        vrs['tau'] = val;
        vrs['tau_err'] = valErr;
        vrs['tau_chi2'] = valChi2;
    def fillVtx(self, obj, vrs):
        if not obj: return
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
            self.hlt2Tool.setOfflineInput()
            self.hlt2Tool.addToOfflineInput(obj)
            dec = self.hlt2Tool.triggerTisTos(obj,"Hlt2ExoticaQuadMuonNoIPDecision")
            vrs['hlt2_quadmu_dec'] = dec.decision()
            vrs['hlt2_quadmu_tis'] = dec.tis()
            vrs['hlt2_quadmu_tos'] = dec.tos()
        except:
            pass
        try:
            self.hlt1Tool.setOfflineInput()
            self.hlt1Tool.addToOfflineInput(obj)
            dec = self.hlt1Tool.triggerTisTos(obj,"Hlt1MultiDiMuonNoIPDecision")
            vrs['hlt1_multidimu_dec'] = dec.decision()
            vrs['hlt1_multidimu_tis'] = dec.tis()
            vrs['hlt1_multidimu_tos'] = dec.tos()
            dec = self.hlt1Tool.triggerTisTos(obj,"Hlt1DiMuonNoIPDecision")
            vrs['hlt1_dimu_dec'] = dec.decision()
            vrs['hlt1_dimu_tis'] = dec.tis()
            vrs['hlt1_dimu_tos'] = dec.tos()
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

    def addX(self, obj, head, hpv, pre="x"):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        trks = []
        neus = []
        decay = []
        try:
            decay += obj.daughters()
            while len(decay) > 0:
                prt = decay[-1]
                try:
                    prt = prt.target()
                except:
                    pass
                dpid = prt.particleID().abspid()
                if dpid in [11,13,211,321,2212]:
                    trks.append(self.addTrk(prt, head, hpv))
                    decay.pop()
                elif dpid in [22,111,221]:
                    decay.pop()
                    neus.append(self.addNeu(prt))
                else:
                    decay.pop()
                    decay += prt.daughters()
            bestpv = self.pvrTool.relatedPV(obj, self.primaries)
            self.fillPvr(bestpv, vrs)
            self.fillDst(obj, hpv, 'ip', vrs)
            self.fillDst(obj, head.endVertex(), 'ip', vrs, "_bvtx")
            self.fillDst(obj, bestpv, 'ip', vrs, "_best")
            self.fillDst(obj.endVertex(), hpv, 'fd', vrs)
            self.fillDst(obj.endVertex(), head.endVertex(), 'fd', vrs, "_bvtx")
            self.fillDst(obj.endVertex(), bestpv, 'fd', vrs, "_best")
            self.fillLifetime(obj, head.endVertex(), vrs)
            self.fillMom(obj.momentum(), vrs)
            self.fillPos(obj.endVertex(), vrs)
            self.fillVtx(obj.endVertex(), vrs)
            self.fillTrig(obj, vrs)
            vrs['m'] = obj.momentum().M()
            for trkidx, trk in enumerate(trks):
                vrs['idx_trk%i' % trkidx] = trk
            for neuidx, neu in enumerate(neus):
                vrs['idx_neu%i' % neuidx] = neu

            self.saved[pre][key] = idx
            self.fill(pre, vrs = vrs)
            return idx
        except:
            print "except addX"
            #raise
            return -1

    def addKst(self, obj, head, hpv, pre="kst"):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        try:
            for dau in obj.daughters():
                try:
                    dau = dau.target()
                except:
                    pass
                dpid = dau.particleID().abspid()
                if dpid==321:
                    vrs['idx_k'] = self.addTrk(dau, head, hpv)
                if dpid==211:
                    vrs['idx_pi'] = self.addTrk(dau, head, hpv)
            bestpv = self.pvrTool.relatedPV(obj, self.primaries)
            self.fillPvr(bestpv, vrs)
            self.fillDst(obj, hpv, 'ip', vrs)
            self.fillDst(obj, head.endVertex(), 'ip', vrs, "_bvtx")
            self.fillDst(obj, bestpv, 'ip', vrs, "_best")
            self.fillDst(obj.endVertex(), hpv, 'fd', vrs)
            self.fillDst(obj.endVertex(), head.endVertex(), 'fd', vrs, "_bvtx")
            self.fillDst(obj.endVertex(), bestpv, 'fd', vrs, "_best")
            self.fillLifetime(obj, head.endVertex(), vrs)
            self.fillMom(obj.momentum(), vrs)
            self.fillPos(obj.endVertex(), vrs)
            self.fillVtx(obj.endVertex(), vrs)
            self.fillTrig(obj, vrs)
            vrs['m'] = obj.momentum().M()

            self.fill(pre, vrs = vrs)
            self.saved[pre][key] = idx
            return idx
        except:
            print "except addKst"
            return -1

    def addB(self, obj, pre="b", mode=0):
        vrs = {}
        try:
            pvr = self.pvrTool.relatedPV(obj, self.primaries)
            for dau in obj.daughters():
                try:
                    # Sometimes we get smartpointers instead of particles
                    dau = dau.target()
                except:
                    pass
                dpid = dau.particleID().abspid()
                if dpid==321:
                    vrs['idx_k'] = self.addTrk(dau, obj, pvr)
                elif dpid==211:
                    vrs['idx_pi'] = self.addTrk(dau, obj, pvr)
                elif dpid==313:
                    vrs['idx_kst'] = self.addKst(dau, obj, pvr)
                    for gdau in dau.daughters():
                        gdpid = gdau.particleID().abspid()
                        if gdpid==321:
                            vrs['idx_k'] = self.addTrk(gdau, obj, pvr)
                        if gdpid==211:
                            vrs['idx_pi'] = self.addTrk(gdau, obj, pvr)
                elif dpid==22:
                    pass
                else:
                    vrs['idx_x'] = self.addX(dau, obj, pvr)
            self.fillPvr(pvr, vrs)
            self.fillDst(obj, pvr, 'ip', vrs)
            self.fillDst(obj.endVertex(), pvr, 'fd', vrs)
            self.fillLifetime(obj, pvr, vrs)
            self.fillMom(obj.momentum(), vrs)
            self.fillPos(obj.endVertex(), vrs)
            self.fillVtx(obj.endVertex(), vrs)
            self.fillTrig(obj, vrs)
            vrs['m'] = obj.momentum().M()
            vrs['mode'] = mode

            self.fill(pre, vrs = vrs)
        except:
            print "except addB"
            #raise
            pass

    def addGen(self, obj, pre = 'gen', par = None):
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
    def addNeu(self, obj, pre = 'neu'):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        self.fillMom(obj.momentum(), vrs)
        self.fillPid(obj.particleID(), vrs)
        #self.fillGen(obj, vrs)
        self.fillTrig(obj, vrs)
        self.fill(pre, vrs = vrs)
        return idx
    def addTrk(self, obj, head, hpv, pre = 'trk'):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        bestpv = self.pvrTool.relatedPV(obj, self.primaries)
        self.fillMom(obj.momentum(), vrs)
        self.fillPid(obj.particleID(), vrs)
        self.fillPro(obj.proto(), vrs)
        self.fillDst(obj, hpv, 'ip', vrs)
        self.fillDst(obj, head.endVertex(), 'ip', vrs, "_bvtx")
        self.fillDst(obj, bestpv, 'ip', vrs, "_best")
        self.fillTrk(obj.proto().track(), obj.particleID(), vrs)
        self.fillPvr(bestpv, vrs)
        #self.fillGen(obj, vrs)
        self.fillTrig(obj, vrs)
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
        return idx

