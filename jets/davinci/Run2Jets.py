#!/bin/env python

#from PhysConf.Filters import LoKi_Filters
#fltrs = LoKi_Filters (
#    STRIP_Code = """
#   HLT_PASS_RE ( 'StrippingHltQEEJetsDiJet.*'    )
#   """
#    )

# Data type configuration.
from GaudiKernel import SystemOfUnits as Units
Type     = 'MC'
JetPtMin = 10 * Units.GeV

## Data.
#from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000001_1.ldst'], #/tmp/dcraik/00042952_00000002_1.ldst'], #/data/dst/MC15.MD.49000004.1.00.dst'],
#                            clear = True)
##Type = 'MC'

# Create the generated jets.
from Configurables import McParticleFlow, McJetBuilder
genPF = McParticleFlow('genPF')
genPF.Inputs = [
    ['PID',        'ban',       '12,-12,14,-14,16,-16'],
    ['PID',        'particle',  '321,211,130,3222,310,3122,3112,3312,3322,'
     '30221,9010221,-321,-211,-130,-3222,-310,-3122,-3112,-3312,-3322,-30221,'
     '-9010221,443,-443'],
    ['MCParticle', 'daughters', 'MC/Particles']
    ]
genPF.Output = 'Phys/PF/MCParticles'
genJB = McJetBuilder('genJB')
genJB.JetPtMin = JetPtMin
genJB.JetR = 0.5
genJB.ChrVrt = True
genJB.NeuVrt = True
genJB.Inputs = [genPF.Output]
genJB.Output = 'Phys/JB/MCParticles'

# Create the reconstructed jets.
from Configurables import HltParticleFlow, HltJetBuilder
from StandardParticles import (StdLooseKsDD, StdLooseKsLL, StdLooseKsLD,
                               StdLooseLambdaDD, StdLooseLambdaLL, 
                               StdLooseLambdaLD)
recPF = HltParticleFlow('recPF')
recPF.Inputs = [
    ['Particle',       'particle', StdLooseKsDD.outputLocation()],
    ['Particle',       'particle', StdLooseKsLL.outputLocation()],
    ['Particle',       'particle', StdLooseKsLD.outputLocation()],
    ['Particle',       'particle', StdLooseLambdaDD.outputLocation()],
    ['Particle',       'particle', StdLooseLambdaLL.outputLocation()],
    ['Particle',       'particle', StdLooseLambdaLD.outputLocation()],
    ['ProtoParticle',  'best',     'Rec/ProtoP/Charged'],
    ['ProtoParticle',  'gamma',    'Rec/ProtoP/Neutrals']
    ]
recPF.Output = 'Phys/PF/Particles'
recPF.ProBestNames = ['mu+', 'e+', 'p+', 'K+', 'pi+']
recPF.ProBestKeys  = [701,   700,  704,  703,  702]
recPF.ProBestMins  = [0.5,   0.5,  0.5,  0.5,  0.5]
recPF.EcalBest = True
recPF.SprRecover = False
recPF.TrkLnErrMax = 10
recPF.TrkUpErrMax = 10
recPF.TrkDnErrMax = 10
recJB = HltJetBuilder('recJB')
recJB.JetEcPath = ''
recJB.Inputs = [recPF.Output]
recJB.Output = 'Phys/JB/Particles'
recJB.JetPtMin = JetPtMin

from StandardParticles import StdAllNoPIDsKaons    as loosekaons
from StandardParticles import StdAllNoPIDsPions    as loosepions
from StandardParticles import StdAllNoPIDsProtons  as looseprotons
from StandardParticles import StdLoosePions        as longpions
from StandardParticles import StdNoPIDsDownPions   as downpions
from StandardParticles import StdLooseProtons      as longprotons
from StandardParticles import StdNoPIDsDownProtons as downprotons
from StandardParticles import StdLooseKsLL      as looseKSLL
from StandardParticles import StdLooseKsDD      as looseKSDD
from StandardParticles import StdLooseLambdaLL  as looseLambdaLL
from StandardParticles import StdVeryLooseLambdaLL  as vlooseLambdaLL
from StandardParticles import StdLooseLambdaDD  as looseLambdaDD
from StandardParticles import StdAllNoPIDsMuons as loosemuons
from PhysSelPython.Wrappers import SimpleSelection, MergedSelection
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles
from itertools import combinations_with_replacement

pids   = ['K+', 'K-', 'KS0', 'Lambda0', 'Lambda~0']
basic  = "(ABSID=='K+')"
combos = list(combinations_with_replacement(pids, 2))
decays = ['K*(892)0 -> ' + ' '.join(combo) for combo in combos]

kaons = SimpleSelection (
    'kaons'           ,
    FilterDesktop   ,
    [ loosekaons ]    ,
    DecayDescriptor = "[K+]cc",
    Code = ("(MIPCHI2DV(PRIMARY)>16.) "
            "& (PT>500*MeV)"
            "& (TRGHOSTPROB<0.2)")
    )

combKS = MergedSelection (
    'combKS',
    RequiredSelections =  [looseKSLL, looseKSDD]
    )

combLambda = MergedSelection (
    'combLambda',
    RequiredSelections =  [vlooseLambdaLL, looseLambdaDD]
    )

allSVs = []

for decay in decays:
    inputs = []
    if 'K+' in decay or 'K-' in decay:
        inputs.append(kaons)
    if 'KS0' in decay:
        inputs.append(combKS)
    if 'Lambda0' in decay or 'Lambda~0' in decay:
        inputs.append(combLambda)
    allSVs.append(
        SimpleSelection (
        'allSVs%d' % (len(allSVs)),
        CombineParticles,
        inputs,
        DecayDescriptor = decay,
        CombinationCut = ("(APT > 2000 * MeV) "
                          "& (ANUM((ID=='KS0')|(ABSID=='Lambda0')) < 2) "
                          "& (ACUTDOCACHI2(1000, '')) "
                          "& ((AALLSAMEBPV(-1, -1, -1) | (AMINCHILD(MIPCHI2DV(PRIMARY)) > 16)) "
                          "| (ANUM((ID == 'KS0') | (ABSID == 'Lambda0')) > 0)) "
                          "& (AM < 10000 * MeV) "
                          "& (ANUM(" + basic + " & (MIPCHI2DV(PRIMARY) < 16)) < 2) "),
        MotherCut      =  ("(HASVERTEX)"
                           "& (VFASPF(VCHI2) < 1000) "
                           "& (BPVVDCHI2 > 16) "
                           "& (in_range(2, BPVETA, 5))"),
        )
    )

combSVs = MergedSelection (
    'combSVs',
    RequiredSelections =  allSVs
    )

pid  = "((ABSID=='K+') | (ID=='KS0') | (ABSID=='Lambda0'))"

recSVs = SimpleSelection (
    'recSVs'           ,
    FilterDesktop   ,
    [ combSVs ]    ,
    DecayDescriptor = "K*(892)0",
    Code = ("(MINTREE(" + pid + ",PT) > 500*MeV) "
            "& (MINTREE(ISBASIC,TRGHOSTPROB) < 0.2) "
            "& (MINTREE((ABSID=='K+'),MIPCHI2DV(PRIMARY)) > 16) "
            "& (HASVERTEX) & (VFASPF(VCHI2PDOF) < 10) "
            "& (BPVVDCHI2 > 25)")
    )


recMus = SimpleSelection (
    'recMus'         ,
    FilterDesktop   ,
    [ loosemuons ]    ,
    DecayDescriptor = "[mu+]cc",
    Code = ("(PROBNNmu>0.5) "
            "& (PT>1000*MeV)"
            "& (TRGHOSTPROB<0.2)")
    )

from PhysSelPython.Wrappers import SelectionSequence
recSVs_seq = SelectionSequence('recSVs_Seq', TopSelection=recSVs)
recMus_seq = SelectionSequence('recMus_Seq', TopSelection=recMus)


###########################
###### D0 candidates ######

charmPions = SimpleSelection (
    'charmPions'         ,
    FilterDesktop   ,
    [ loosepions ]    ,
    DecayDescriptor = "[pi+]cc",
    Code = ( #"(PIDK - PIDpi < 3) & "
            "(PT>200*MeV)"
            "& (TRGHOSTPROB<0.3)"
            "& (MIPCHI2DV(PRIMARY) > 4)")
    )

charmKaons = SimpleSelection (
    'charmKaons'         ,
    FilterDesktop   ,
    [ loosekaons ]    ,
    DecayDescriptor = "[K+]cc",
    Code = (#"(PIDK - PIDpi > 5) & "
            "(PT>200*MeV)"
            "& (TRGHOSTPROB<0.3)"
            "& (MIPCHI2DV(PRIMARY) > 4)")
    )

dcD0 = { }
for child in ['pi+','K+'] :
    dcD0[child] = "(PT > 250*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 16)"
#                  "& (TRCHI2 < 3)" \

combcutsD0 = "in_range(1784*MeV,  AM, 1944*MeV)" \
             "& (AMINDOCA('') < 0.1*mm )"

parentcutsD0 = "(VFASPF(VCHI2PDOF) < 10)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2> 49 )" \
               "& (BPVDIRA > 0.99985 )"

recD0 = SimpleSelection (
    'recD0',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D0 -> K- pi+]cc",
    DaughtersCuts   = dcD0,
    CombinationCut = (combcutsD0),
    MotherCut      =  (parentcutsD0),
)

dcDp = { }
for child in ['pi+','K+'] :
    dcDp[child] = "(PT > 200*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"
#                  "& (TRCHI2 < 3)" \

combcutsDp = "in_range(1789*MeV,  AM, 1949*MeV)" \
             "& (ANUM(PT > 400*MeV) > 1)" \
             "& (ANUM(PT > 1000*MeV) > 0)" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 10) > 1)" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 50) > 0)"

parentcutsDp = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVLTIME() > 0.150*ps )" \
               "& (BPVDIRA > 0.9994 )"

recDp = SimpleSelection (
    'recDp',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D+ -> K- pi+ pi+]cc",
    DaughtersCuts   = dcDp,
    CombinationCut = (combcutsDp),
    MotherCut      =  (parentcutsDp),
)

dcDs = { }
for child in ['pi+','K+'] :
    dcDs[child] = "(PT > 200*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"
#                  "& (TRCHI2 < 3)" \

combcutsDs = "in_range(1889*MeV,  AM, 2049*MeV)" \
             "& (ANUM(PT > 400*MeV) > 1)" \
             "& (ANUM(PT > 1000*MeV) > 0)" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 10) > 1)" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 50) > 0)"

parentcutsDs = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVLTIME() > 0.150*ps )" \
               "& (BPVDIRA > 0.9994 )"

recDs = SimpleSelection (
    'recDs',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D_s+ -> K- K+ pi+]cc",
    DaughtersCuts   = dcDs,
    CombinationCut = (combcutsDs),
    MotherCut      =  (parentcutsDs),
)

dcD02K3pi = { }
for child in ['pi+','K+'] :
    dcD02K3pi[child] = "(PT > 200*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"
#                  "& (TRCHI2 < 3)" \

combcutsD02K3pi = "in_range(1784*MeV,  AM, 1944*MeV)" \
             "& (AMINDOCA('') < 0.1*mm )"

parentcutsD02K3pi = "(VFASPF(VCHI2PDOF) < 10)" \
               "& BPVVALID()" \
               "& (BPVLTIME() > 0.10*ps )" \
               "& (BPVDIRA > 0.99985 )"

recD02K3pi = SimpleSelection (
    'recD02K3pi',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D0 -> K- pi+ pi+ pi-]cc",
    DaughtersCuts   = dcD02K3pi,
    CombinationCut = (combcutsD02K3pi),
    MotherCut      =  (parentcutsD02K3pi),
)

D0_seq = SelectionSequence('D0_Seq', TopSelection=recD0)
Dp_seq = SelectionSequence('Dp_Seq', TopSelection=recDp)
Ds_seq = SelectionSequence('Ds_Seq', TopSelection=recDs)
D02K3pi_seq = SelectionSequence('D2K3pi0_Seq', TopSelection=recD02K3pi)

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = True
DaVinci().appendToMainSequence([genPF, genJB, recPF, recJB])
DaVinci().appendToMainSequence([recSVs_seq.sequence(), recMus_seq.sequence(), D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(), D02K3pi_seq.sequence()])
DaVinci().DataType = '2015'
#DaVinci().EventPreFilters = fltrs.filters ('Filters')

# Configure the BDT tagger.
from Configurables import LoKi__BDTTag
tagger = LoKi__BDTTag()
tagger.NbvSelect = False

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
    def __init__(self, name, tes, toolSvc, detSvc):
        """
        Initialize the ntuple with the needed tools.
        """
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
        self.tes     = tes
        self.saved   = {}
        self.ntuple  = OrderedDict()
        self.tfile   = ROOT.TFile('output.root', 'RECREATE')
        self.ttree   = ROOT.TTree('data', 'data')
        self.vrs     = {}
        mom = ['px', 'py', 'pz', 'e']
        pos = ['x', 'y', 'z']
        cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
        self.init('gen', ['idx_pvr', 'idx_jet', 'pid', 'q'] + mom + pos)
        self.init('pvr', pos + cov)
        self.init('svr', ['idx_pvr', 'idx_jet'] + [
                'idx_trk%i' % i for i in range(0, 10)] + 
                  mom + pos + ['m', 'm_cor', 'pt', 'fd_min', 'fd_chi2', 'chi2', 'ip_chi2_sum', 'abs_q_sum', 'tau', 'ntrk', 'ntrk_jet', 'jet_dr', 'jet_pt', 'pass', 'bdt0', 'bdt1'])
        self.init('jet', ['idx_pvr', 'ntrk', 'nneu'] + mom) #+ ['idx_trk%i' % i for i in range(0, 40)])
        self.init('trk', ['idx_gen', 'idx_pvr', 'idx_jet'] + mom +
                  ['pid', 'q', 'ip', 'ip_chi2', 'pnn_e', 'pnn_mu', 'pnn_pi',
                   'pnn_k', 'pnn_p', 'pnn_ghost', 'ecal', 'hcal', 'prb_ghost', 'type', 'is_mu',
                   'vid', 'x', 'y', 'z']) #+ ['vhit%i' % i for i in range(0, 61)])# + [
                   #'vz%i' % i for i in range(0, 61)])
        self.init('neu', ['idx_gen', 'idx_jet'] + mom + ['pid'])
        self.init('d0', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 2)]) 
        self.init('dp', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)]) 
        self.init('ds', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)]) 
        self.init('d02k3pi', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 4)]) 
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
                #vrs['vz%i' % (len(vids) - 1)] = vids[-1][0]
                #vrs['vhit%i' % (len(vids) - 1)] = vids[-1][1]
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

    def addDHad(self, obj, pre="d0"):
        vrs = {}
        trks = []
        try:
            jets = tes[recJB.Output]
            nInBest = 0
            bestJet = -1
            for idx, jet in enumerate(jets):
                jetTrkKeys = [ self.key(dau) for dau in jet.daughters() ]
                nIn=0
                for dau in obj.daughters():
                    key = self.key(dau)
                    #if key in self.saved['trk']:
                    #    trks.append(self.saved['trk'][key])
                    #else :
                    #    trks.append(self.addTrk(dau))
                    if key in jetTrkKeys:
                        nIn+=1
                if nIn>nInBest:
                    nInBest=nIn
                    bestJet=idx
            
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
            vrs['m'] = obj.momentum().M()
            vrs['idx_jet'] = bestJet
            vrs['ntrk_jet'] = nInBest
            for idx, trk in enumerate(trks):
                vrs['idx_trk%i' % idx] = trk

            self.fill(pre, vrs = vrs)
        except:
            pass

    def addGen(self, obj, jet = -1, pre = 'gen'):
        key = self.key(obj)
        if key in self.saved[pre]: return self.saved[pre][key]
        vrs = {}
        idx = len(self.saved[pre])
        self.fillPid(obj.particleID(), vrs)
        self.fillMom(obj.momentum(), vrs)
        self.fillPos(obj.originVertex(), vrs)
        self.fillPvr(obj.primaryVertex(), vrs)
        vrs['idx_jet'] = jet
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
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
    def addSvr(self, obj, pre = 'svr'):
        vrs = {}
        pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
        self.fillMom(obj.momentum(), vrs)
        self.fillPos(obj.endVertex(), vrs)
        self.fillCov(obj.endVertex(), vrs)
        self.fillPvr(pvr, vrs)
        self.fillDst(obj.endVertex(), pvr, 'fd', vrs)
        trks = []
        for dtr in obj.daughters():
            if not dtr.proto() or not dtr.proto().track(): continue
            trks += [[self.addTrk(dtr), dtr.proto().track()]]
            vrs['idx_trk%i' % (len(trks) - 1)] = trks[-1][0]
        self.fill(pre, vrs = vrs)
    def addTags(self, obj, jet = -1, pre = 'svr'):
#        stash.append(pre)
        tags = STD.map('string', 'double')()
        if not self.tagTool.calculateJetProperty(obj, tags): return
        ntag = int(tags['Tag'])
        for itag in range(0, ntag):
#            stash.append("tag")
            vrs = {}
            vrs['idx_pvr'] = self.addPvr(self.tes['Rec/Vertex/Primary']
                                         [int(tags['Tag%i_idx_pvr' % itag])])
            vrs['idx_jet'] = jet
            ntrk = int(tags['Tag%i_nTrk' % itag])
#            stash.append(ntrk)
            for itrk in range(0, ntrk): vrs['idx_trk%i' % itrk] = self.addTrk(
                self.tes['Phys/StdAllNoPIDsPions/Particles']
                [int(tags['Tag%i_idx_trk%i' % (itag, itrk)])]);
            for vr in ['x', 'y', 'z', 'px', 'py', 'pz', 'e']:
                vrs[vr] = tags['Tag%i_%s' % (itag, vr)]
            vrs['m']  = tags['Tag%i_m' % itag]
            vrs['m_cor']  = tags['Tag%i_mCor' % itag]
            vrs['pt']  = tags['Tag%i_pt' % itag]
            vrs['fd_min'] = tags['Tag%i_fdrMin' % itag]
            vrs['fd_chi2']  = tags['Tag%i_fdChi2' % itag]
            vrs['chi2']  = tags['Tag%i_chi2' % itag]
            vrs['ip_chi2_sum']  = tags['Tag%i_ipChi2Sum' % itag]
            vrs['abs_q_sum']  = tags['Tag%i_absQSum' % itag]
            vrs['tau']  = tags['Tag%i_tau' % itag]
            vrs['ntrk']  = tags['Tag%i_nTrk' % itag]
            vrs['ntrk_jet']  = tags['Tag%i_nTrkJet' % itag]
            vrs['jet_dr']  = tags['Tag%i_drSvrJet' % itag]
            vrs['jet_pt']  = tags['Tag%i_ptSvrJet' % itag]
            vrs['pass']  = tags['Tag%i_pass' % itag]
            vrs['bdt0']  = tags['Tag%i_bdt0' % itag]
            vrs['bdt1']  = tags['Tag%i_bdt1' % itag]
            self.fill(pre, vrs = vrs)
    def addJet(self, obj, pre = 'jet'):
        vrs = {}
        idx = self.ntuple['jet_idx_pvr'].size();
        pvr = self.pvrTool.relatedPV(obj, 'Rec/Vertex/Primary')
        self.fillMom(obj.momentum(), vrs)
        self.fillPvr(pvr, vrs)
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
                #vrs['idx_trk%i' % (len(trks) - 1)] = trks[-1][0]
        vrs['ntrk'] = len(trks)
        vrs['nneu'] = nneu
        self.addTags(obj, idx)
        self.fill(pre, vrs = vrs)
    def addNeu(self, obj, jet = -1, pre = 'neu'):
        vrs = {}
        self.fillMom(obj.momentum(), vrs)
        self.fillPid(obj.particleID(), vrs)
        self.fillGen(obj, vrs)
        vrs['idx_jet'] = jet
        self.fill(pre, vrs = vrs)
    def addTrk(self, obj, jet = -1, pre = 'trk'):
#        printStash=False
#        if stash[-2] == "trk": printStash=True
        key = self.key(obj)
#        if printStash: stash.append("addTrk")
#        if printStash: stash.append(key)
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
        vrs['idx_jet'] = jet
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
#        if printStash: stash.append(vrs)
#        if printStash: stash.append("endTrk")
        return idx

# GaudiPython configuration.
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

# Run.
import sys, ROOT
from math import floor
#stash = []
evtmax = -1
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc())
while evtmax < 0 or evtnum < evtmax:
#    stash.append("...")
    gaudi.run(1)
    if not bool(tes['/Event']): break
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    try: ntuple.ntuple['evt_pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
    except: continue
    try: ntuple.ntuple['evt_trk_n'][0] = len(tes['Phys/StdAllNoPIDsPions/Particles'])
    except: continue

    # Fill generator level info.
    fill = False;
    gens = tes['MC/Particles']
    try:
        ntuple.addGen(gens[0])
        ntuple.addGen(gens[1])
    except: pass
    try:
        for gen in gens:
            pid = gen.particleID()
            if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
                ntuple.addGen(gen)
    except: pass
    try:
        jets = tes[genJB.Output]
        for jet in jets: ntuple.addGen(jet); fill = True
    except: pass

    # Fill reconstructed.
    try:
        jets = tes[recJB.Output]
        for jet in jets:
            ntuple.addJet(jet); fill = True;
    except: pass

    # fill other tracks
    try:
        for trk in tes['Phys/StdAllNoPIDsPions/Particles']:
            ntuple.addTrk(trk);
    except: pass

    # fill D's
    try:
        d0s = tes[recD0.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"d0")
        dps = tes[recDp.algorithm().Output]
        for dp in dps:
            ntuple.addDHad(dp,"dp")
        dss = tes[recDs.algorithm().Output]
        for ds in dss:
            ntuple.addDHad(ds,"ds")
        d0s = tes[recD02K3pi.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"d02k3pi")
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
