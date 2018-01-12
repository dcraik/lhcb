#!/bin/env python
import sys
sys.path.append('.')

# Data type configuration.
from GaudiKernel import SystemOfUnits as Units
Type     = 'MC'
JetPtMin = 10 * Units.GeV

### Data.
#from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles([
##'/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000001_1.ldst' #/tmp/dcraik/00042952_00000002_1.ldst' #/data/dst/MC15.MD.49000004.1.00.dst'
##    '/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00042982/0000/00042982_00000002_1.ldst'#light
#    '/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00042950/0000/00042950_00000001_1.ldst'#charm
##    '/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00042972/0000/00042972_00000003_1.ldst'#beauty
#    ],
#    clear = True)
##Type = 'MC'

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
     '-321,-211,-130,-3222,-310,-3122,-3112,-3312'],
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
    ['Particle',       'daughters', Zs.outputLocation()],
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
recSVs_seq = SelectionSequence('recSVs_Seq', TopSelection=recSVs)
recMus_seq = SelectionSequence('recMus_Seq', TopSelection=recMus)

Jpsi_seq = SelectionSequence('Jpsi_Seq', TopSelection=recJpsi)
D0_seq = SelectionSequence('D0_Seq', TopSelection=recD0)
Dp_seq = SelectionSequence('Dp_Seq', TopSelection=recDp)
Ds_seq = SelectionSequence('Ds_Seq', TopSelection=recDs)
Lc_seq = SelectionSequence('Lc_Seq', TopSelection=recLc)
D02K3pi_seq = SelectionSequence('D2K3pi0_Seq', TopSelection=recD02K3pi)

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = True
DaVinci().appendToMainSequence([Z_seq.sequence()])
DaVinci().appendToMainSequence([genPF, genJB, recPF, recJB])
#DaVinci().appendToMainSequence([recSVs_seq.sequence(), recMus_seq.sequence()])
DaVinci().appendToMainSequence([Jpsi_seq.sequence(), D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(),  Lc_seq.sequence(), D02K3pi_seq.sequence()])
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
    try: ntuple.ntuple['evt_trk_n'][0] = len(tes['Phys/StdAllNoPIDsPions/Particles'])
    except: continue

    # Fill generator level info.
    fill = False;
    gens = tes['MC/Particles']
    try:
        for gen in gens:
            pid = gen.particleID()
            if pid.isQuark() or pid.pid() == 21:
                ntuple.addGen(gen)
        #ntuple.addGen(gens[0])
        #ntuple.addGen(gens[1])
    except: pass
    try:
        for gen in gens:
            pid = gen.particleID()
            if pid.pid() == 23:
                ntuple.addGen(gen)
    except: pass
    try:
        for gen in gens:
            pid = gen.particleID()
            if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
                ntuple.addGen(gen);# fill = True
    except: pass
    try:
        jets = tes[genJB.Output]
        for jet in jets: ntuple.addGen(jet);# fill = True
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
        d0s = tes[recD0.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"d0")
        dps = tes[recDp.algorithm().Output]
        for dp in dps:
            ntuple.addDHad(dp,"dp")
        dss = tes[recDs.algorithm().Output]
        for ds in dss:
            ntuple.addDHad(ds,"ds")
        lcs = tes[recLc.algorithm().Output]
        for lc in lcs:
            ntuple.addDHad(lc,"lc")
        d0s = tes[recD02K3pi.algorithm().Output]
        for d0 in d0s:
            ntuple.addDHad(d0,"k3pi")
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
