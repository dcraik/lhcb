#!/bin/env python
import sys
sys.path.append('.')

from PhysConf.Filters import LoKi_Filters

fltrs = LoKi_Filters (
    STRIP_Code = """
    HLT_PASS_RE ( 'StrippingZ02MuMu.*' )
    """
    )

# Data type configuration.
from GaudiKernel import SystemOfUnits as Units
##Type     = 'MC'
JetPtMin = 10 * Units.GeV

m4_ifdef(`TEST', from GaudiConf import IOHelper)
m4_ifdef(`TEST', IOHelper('ROOT').inputFiles(['/data/zjet/00103400_00000032_1.ew.dst'],clear = True))

from StandardParticles import StdAllNoPIDsMuons as loosemuons
from PhysSelPython.Wrappers import SimpleSelection, MergedSelection, DataOnDemand, Selection
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles


from commonSelections import *

from PhysSelPython.Wrappers import SelectionSequence
Z_seq = SelectionSequence('Z_Seq', TopSelection=Zs)

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
    #['Particle',       'daughters', 'Phys/Z02MuMuLine/Particles'],
    ['Particle',       'daughters', Zs.outputLocation()],
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
DaVinci().Simulation = False
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
DaVinci().appendToMainSequence([Z_seq.sequence()])
DaVinci().appendToMainSequence([recPF, recJB])
#DaVinci().appendToMainSequence([recSVs_seq.sequence(), recMus_seq.sequence()])
DaVinci().appendToMainSequence([Jpsi_seq.sequence(), D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(),  Lc_seq.sequence(), D02K3pi_seq.sequence()])
##TODO adding recSVs and recMus changes the daughters of jet objects from smart poniters to Particles
DaVinci().DataType = 'M4_DATAYEAR'
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
evtmax = m4_ifdef(`TEST',10000,-1)
#'
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc(), recJB.Output, recSVs.outputLocation(), recMus.algorithm().Output)
#gaudi.run(39000)#TODO
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

    fill = False;

    # Fill reconstructed.
    try:
        jets = tes[recJB.Output]
        #ntuple.addTrigger()
        for jet in jets:
            ntuple.addJet(jet); #fill = True;
    except: pass

    # fill other tracks
    try:
        for trk in tes['Phys/StdAllNoPIDsPions/Particles']:
            ntuple.addTrk(trk);
    except: pass

    ## fill Z
    try:
        for z in tes[Zs.outputLocation()]:
            if ntuple.addZ(z):
                fill = True;
    except: pass

    # fill D's
    try:
        jpsis = tes[recJpsi.algorithm().Output]
        for jpsi in jpsis:
            ntuple.addDHad(jpsi,"jpsi")
    except: pass
    try:
        d0s = tes[recD0.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"d0")
    except: pass
    try:
        dps = tes[recDp.algorithm().Output]
        for dp in dps:
            ntuple.addDHad(dp,"dp")
    except: pass
    try:
        dss = tes[recDs.algorithm().Output]
        for ds in dss:
            ntuple.addDHad(ds,"ds")
    except: pass
    try:
        lcs = tes[recLc.algorithm().Output]
        for lc in lcs:
            ntuple.addDHad(lc,"lc")
    except: pass
    try:
        d0s = tes[recD02K3pi.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"k3pi")
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
