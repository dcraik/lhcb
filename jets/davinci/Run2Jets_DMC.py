#!/bin/env python
import sys
sys.path.append('.')

#from PhysConf.Filters import LoKi_Filters
#fltrs = LoKi_Filters (
#    HLT_Code = """
#    HLT_PASS_RE ( 'L0DiMuonDecision' )
#    & HLT_PASS_RE ( 'Hlt1DiMuonHighMassDecision' )
#    & HLT_PASS_RE ( 'Hlt2DiMuonB.*Decision' )
#   """,
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
##IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.DST/00057115/0000/00057115_00000003_3.AllStreams.dst'],#/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'],#/tmp/dcraik/00042952_00000002_1.ldst'], #/data/dst/MC15.MD.49000004.1.00.dst'],
##IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00041855/0000/00041855_00000011_1.ldst'],#/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.DST/00057115/0000/00057115_00000003_3.AllStreams.dst'],#/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'],#/tmp/dcraik/00042952_00000002_1.ldst'], #/data/dst/MC15.MD.49000004.1.00.dst'],
#IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/2015/DST/00046982/0000/00046982_00000003_2.dst'],#K3pi
##IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.DST/00061410/0000/00061410_00000002_7.AllStreams.dst'],#Lc
##IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/2015/DST/00046269/0000/00046269_00000001_2.dst'],
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

from commonSelections import *

from PhysSelPython.Wrappers import SelectionSequence
#recSVs_seq = SelectionSequence('recSVs_Seq', TopSelection=recSVs)
#recMus_seq = SelectionSequence('recMus_Seq', TopSelection=recMus)

Jpsi_seq = SelectionSequence('Jpsi_Seq', TopSelection=recJpsi)
D0_seq = SelectionSequence('D0_Seq', TopSelection=recD0)
Dp_seq = SelectionSequence('Dp_Seq', TopSelection=recDp)
Ds_seq = SelectionSequence('Ds_Seq', TopSelection=recDs)
Lc_seq = SelectionSequence('Lc_Seq', TopSelection=recLc)
D02K3pi_seq = SelectionSequence('D2K3pi0_Seq', TopSelection=recD02K3pi)

##########################

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = True
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
#DaVinci().appendToMainSequence([genPF, genJB, recPF, recJB])
DaVinci().appendToMainSequence([recPF, recJB])
#DaVinci().appendToMainSequence([recSVs_seq.sequence(), recMus_seq.sequence()])
DaVinci().appendToMainSequence([Jpsi_seq.sequence(), D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(),  Lc_seq.sequence(), D02K3pi_seq.sequence()])
##TODO adding recSVs and recMus changes the daughters of jet objects from smart poniters to Particles
DaVinci().DataType = '2016'
#DaVinci().EventPreFilters = fltrs.filters ('Filters')

from Configurables import LumiIntegrateFSR, LumiIntegratorConf
LumiIntegrateFSR('IntegrateBeamCrossing').SubtractBXTypes = ['None']

# Configure the BDT tagger.
from Configurables import LoKi__BDTTag
tagger = LoKi__BDTTag()
tagger.NbvSelect = False

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
evtmax = -1
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc(), recJB.Output, recSVs.outputLocation(), recMus.algorithm().Output)
while evtmax < 0 or evtnum < evtmax:
    gaudi.run(1)
    if not bool(tes['/Event']): break
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    try: ntuple.ntuple['evt_pvr_n'][0] = len(tes['Rec/Vertex/Primary'])
    except: continue

    # Fill generator level info.
    fill = False;
    gens = tes['MC/Particles']
#    try:
#        ntuple.addGen(gens[0])
#        ntuple.addGen(gens[1])
#    except: pass
    try:
        for gen in gens:
            pid = gen.particleID()
            if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
                ntuple.addGen(gen); fill = True
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
            ntuple.addTrk(trk); #fill = True;
    except: pass

    # fill D's
    try:
        jpsis = tes[recJpsi.algorithm().Output]
        for jpsi in jpsis:
            ntuple.addDHad(jpsi,"jpsi")
        d0s = tes[recD0.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"d0"); fill = True;
        dps = tes[recDp.algorithm().Output]
        for dp in dps:
            ntuple.addDHad(dp,"dp"); fill = True;
        dss = tes[recDs.algorithm().Output]
        for ds in dss:
            ntuple.addDHad(ds,"ds"); fill = True;
        lcs = tes[recLc.algorithm().Output]
        for lc in lcs:
            ntuple.addDHad(lc,"lc"); fill = True;
        d0s = tes[recD02K3pi.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"k3pi"); fill = True;
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
