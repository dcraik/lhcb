#!/bin/env python
import sys
sys.path.append('.')

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
####DiJet*
    STRIP_Code = """
   HLT_PASS_RE ( 'StrippingHltQEEJetsDiJet.*LineDecision'    )
   """
    )
####DiJetSV{,HighPt,LowPt}
#    STRIP_Code = """
#   HLT_PASS_RE ( 'StrippingHltQEEJetsDiJetSVLineDecision' )
#   | HLT_PASS_RE ( 'StrippingHltQEEJetsDiJetSVHighPtLineDecision' )
#   | HLT_PASS_RE ( 'StrippingHltQEEJetsDiJetSVLowPtLineDecision' )
#   """
#    )
####DiJetSV*
#    STRIP_Code = """
#   HLT_PASS_RE ( 'StrippingHltQEEJetsDiJetSV.*LineDecision' )
#   """
#    )
####DiJet
#    STRIP_Code = """
#   HLT_PASS_RE ( 'StrippingHltQEEJetsDiJetLineDecision'    )
#   """
#    )


# Data type configuration.
from GaudiKernel import SystemOfUnits as Units
##Type     = 'MC'
JetPtMin = 10 * Units.GeV


## Data.
#from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/data/dijets/00069603_00000642_1.bhadroncompleteevent.dst'],
#                            clear = True)
##Type = 'MC'

# Create the generated jets.
from Configurables import McParticleFlow, McJetBuilder
genPF = McParticleFlow('genPF')
genPF.Inputs = [
    ['PID',        'ban',       '12,-12,14,-14,16,-16'],
    ['PID',        'particle',  '321,211,130,3222,310,3122,3112,3312,3322,'
     '-321,-211,-130,-3222,-310,-3122,-3112,-3312,-3322'],
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

#from JetAccessories import HltJetConf
#hltPF = HltJetConf.HltParticleFlowConf('hltPF',['Photons', 'ResolvedPi0s', 'MergedPi0s',
#                                                'Ks', 'Lambda', 'ChargedProtos', 'NeutralProtos',
#                                                'EcalClusters', #'HcalClusters',
#                                                'EcalMatch', #'HcalMatch',
#                                                'PrimaryVtxs' ])
#hltJB = HltJetConf.HltJetBuilderConf('hltJB',hltPF.getOutputLocation(),JetInfo=False,JetEcPath="",Output='Phys/JB2/Particles')

from commonSelections import *

from PhysSelPython.Wrappers import SelectionSequence
recSVs_seq = SelectionSequence('recSVs_Seq', TopSelection=recSVs)
recMus_seq = SelectionSequence('recMus_Seq', TopSelection=recMus)

Jpsi_seq = SelectionSequence('Jpsi_Seq', TopSelection=recJpsi)
D0_seq = SelectionSequence('D0_Seq', TopSelection=recD0)
Dp_seq = SelectionSequence('Dp_Seq', TopSelection=recDp)
Ds_seq = SelectionSequence('Ds_Seq', TopSelection=recDs)
Lc_seq = SelectionSequence('Lc_Seq', TopSelection=recLc)
D02K3pi_seq = SelectionSequence('D2K3pi0_Seq', TopSelection=recD02K3pi)

##########################

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = False
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
#DaVinci().appendToMainSequence([genPF, genJB, recPF, recJB])
DaVinci().appendToMainSequence([recPF, recJB])
#DaVinci().appendToMainSequence([hltPF.getSeq(), hltJB.getSeq()])
DaVinci().appendToMainSequence([recSVs_seq.sequence(), recMus_seq.sequence()])
DaVinci().appendToMainSequence([Jpsi_seq.sequence(),D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(), Lc_seq.sequence(), D02K3pi_seq.sequence()])
##TODO adding recSVs and recMus changes the daughters of jet objects from smart poniters to Particles
DaVinci().DataType = '2016'
DaVinci().EventPreFilters = fltrs.filters ('Filters')

from Configurables import LumiIntegrateFSR, LumiIntegratorConf
LumiIntegrateFSR('IntegrateBeamCrossing').SubtractBXTypes = ['None']

# Configure the BDT tagger.
from Configurables import LoKi__BDTTag
tagger = LoKi__BDTTag()
tagger.NbvSelect = False
tagger = LoKi__BDTTag("Backwards")
tagger.NbvSelect = False
tagger.Backwards = True

from Configurables import ToolSvc, TriggerTisTos
for stage in ('Hlt1', 'Hlt2', 'Strip/Phys'):
    ToolSvc().addTool(TriggerTisTos, stage + "TriggerTisTos")
    tool = getattr(ToolSvc(), stage + "TriggerTisTos")
    tool.HltDecReportsLocation = '/Event/' + stage + '/DecReports'
    tool.HltSelReportsLocation = '/Event/' + stage + '/SelReports'

# Access to classes.
from collections import OrderedDict
import ROOT, array, GaudiPython
from GaudiPython.Bindings import gbl
STD  = gbl.std
LHCB = gbl.LHCb

from Ntuple import Ntuple

# GaudiPython configuration.
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

# Run.
import sys, ROOT
from math import floor
evtmax = -1#1000#-1
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc(), recJB.Output, recSVs.outputLocation(), recMus.algorithm().Output, addBackTags=True)
while evtmax < 0 or evtnum < evtmax:
    gaudi.run(1)
    if not bool(tes['/Event']): break
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    try: ntuple.ntuple['evt_pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
    except: continue
    try: ntuple.addEventInfo();
    except: continue

    # Fill generator level info.
    fill = False;
    #gens = tes['MC/Particles']
    #try:
    #    ntuple.addGen(gens[0])
    #    ntuple.addGen(gens[1])
    #except: pass
    #try:
    #    for gen in gens:
    #        pid = gen.particleID()
    #        if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
    #            ntuple.addGen(gen)
    #except: pass
    #try:
    #    jets = tes[genJB.Output]
    #    for jet in jets: ntuple.addGen(jet); fill = True
    #except: pass

    # Fill reconstructed.
    try:
        jets = tes[recJB.Output]
        ntuple.addTrigger()
        for jet in jets:
            ntuple.addJet(jet); fill = True;
    except: pass

#    # Fill reconstructed.
#    try:
#        jets = tes[hltJB.jb.Output]
#        print jets.size()
#        for jet in jets:
#            ntuple.addJet(jet,'hlt_jet'); fill = True;
#    except: pass

    # fill other tracks
    try:
        for trk in tes['Phys/StdAllNoPIDsPions/Particles']:
            ntuple.addTrk(trk);
    except: pass

    # fill D's
    try:
        jpsis = tes[recJpsi.algorithm().Output]
        p2pvTable = tes[recJpsi.algorithm().Output[:-1]+'2VertexRelations']
        for jpsi in jpsis:
            relPVs = p2pvTable.relations(jpsi)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(jpsi,"jpsi",bestVertex)
    except: pass
    try:
        d0s = tes[recD0.algorithm().Output]
        p2pvTable = tes[recD0.algorithm().Output[:-1]+'2VertexRelations']
        for d0 in d0s:
            relPVs = p2pvTable.relations(d0)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(d0,"d0",bestVertex)
    except: pass
    try:
        dps = tes[recDp.algorithm().Output]
        p2pvTable = tes[recDp.algorithm().Output[:-1]+'2VertexRelations']
        for dp in dps:
            relPVs = p2pvTable.relations(dp)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(dp,"dp",bestVertex)
    except: pass
    try:
        dss = tes[recDs.algorithm().Output]
        p2pvTable = tes[recDs.algorithm().Output[:-1]+'2VertexRelations']
        for ds in dss:
            relPVs = p2pvTable.relations(ds)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(ds,"ds",bestVertex)
    except: pass
    try:
        lcs = tes[recLc.algorithm().Output]
        p2pvTable = tes[recLc.algorithm().Output[:-1]+'2VertexRelations']
        for lc in lcs:
            relPVs = p2pvTable.relations(lc)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(lc,"lc",bestVertex)
    except: pass
    try:
        d0s = tes[recD02K3pi.algorithm().Output]
        p2pvTable = tes[recD02K3pi.algorithm().Output[:-1]+'2VertexRelations']
        for d0 in d0s:
            relPVs = p2pvTable.relations(d0)
            bestVertex = relPVs.back().to()
            if relPVs.size()>1: print evtnum, relPVs.size()
            ntuple.addDHad(d0,"k3pi",bestVertex)
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
