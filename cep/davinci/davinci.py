#!/bin/env python
import sys
sys.path.append('.')

# Test data
#from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst'], clear = True)
#

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = """
   HLT_PASS_RE ( 'StrippingLowMultLMR2HHLineDecision' )
   """
    )


from PhysSelPython.Wrappers import SimpleSelection, MergedSelection, DataOnDemand, Selection
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles

##other selections
from StandardParticles import StdAllNoPIDsKaons    as loosekaons

kaons = SimpleSelection (
    'kaons'         ,
    FilterDesktop   ,
    [ loosekaons ]    ,
    DecayDescriptor = "[K+]cc",
    Code = (#"(PIDK - PIDpi > 5) & "
            "(PT>50*MeV)"
            "& (TRGHOSTPROB<0.3)")
            #"& (MIPCHI2DV(PRIMARY) > 4)")
    )

dcDiK = { }
dcDiK['K+'] = "(PT > 50.0) & (P > 500.0) & (TRGHOSTPROB < 0.3)"
combcutsDiK = "(AMAXDOCA('') < 0.5) & (in_range( 0, AM, 6000 ))"
parentcutsDiK = "(VFASPF(VCHI2PDOF) < 10)"

phis = SimpleSelection (
    'phis',
    CombineParticles,
    [kaons],
    DecayDescriptor = "phi(1020) -> K+ K-",
    DaughtersCuts   = dcDiK,
    CombinationCut = (combcutsDiK),
    MotherCut      =  (parentcutsDiK),
)


from PhysSelPython.Wrappers import SelectionSequence

phi_seq = SelectionSequence('phi_seq', TopSelection=phis)

##########################

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = False
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
DaVinci().appendToMainSequence([phi_seq.sequence()])
DaVinci().DataType = '2016'
DaVinci().EventPreFilters = fltrs.filters ('Filters')

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
evtmax = -1 #TODO
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc())
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

    # fill combinations
    try:
        dis = tes[phis.algorithm().Output]
        for di in dis:
            ntuple.addDi(di,"phi")
            fill = True
    except: pass
    # fill other tracks
    try:
        for trk in tes['Phys/StdAllNoPIDsPions/Particles']:
            ntuple.addTrk(trk);
            #fill = True
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
