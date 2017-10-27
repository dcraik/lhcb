#!/bin/env python

from PhysConf.Filters import LoKi_Filters
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
#IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/Dev/LDST/00041855/0000/00041855_00000011_1.ldst'],#/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.DST/00057115/0000/00057115_00000003_3.AllStreams.dst'],#/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'],#/tmp/dcraik/00042952_00000002_1.ldst'], #/data/dst/MC15.MD.49000004.1.00.dst'],
##IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/MC/2015/DST/00046269/0000/00046269_00000001_2.dst'],
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

#THESE MATCH STDPARTICLES
#KSLL = SimpleSelection (
#    'KSLL',
#    CombineParticles,
#    [longpions],
#    DecayDescriptor = "KS0 -> pi+ pi-",
#    DaughtersCuts = {
#        "pi+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 9.)"
#        },
#    CombinationCut = ("(ADAMASS('KS0') < 50.*MeV) & (ADOCACHI2CUT(25, ''))"),
#    MotherCut      =  ("(ADMASS('KS0') < 35.*MeV) & (VFASPF(VCHI2) < 25.)"),
#)
#
#KSDD = SimpleSelection (
#    'KSDD',
#    CombineParticles,
#    [downpions],
#    DecayDescriptor = "KS0 -> pi+ pi-",
#    DaughtersCuts = {
#        "pi+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 4.)"
#        },
#    CombinationCut = ("(ADAMASS('KS0') < 80.*MeV) & (ADOCACHI2CUT(25, ''))"),
#    MotherCut      =  ("(ADMASS('KS0') < 64.*MeV) & (VFASPF(VCHI2) < 25.)"),
#)
#
#LambdaLL = SimpleSelection (
#    'LambdaLL',
#    CombineParticles,
#    [longpions,longprotons],
#    DecayDescriptor = "[Lambda0 -> p+ pi-]cc",
#    DaughtersCuts = {
#        "p+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 9.)",
#        "pi+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 9.)"
#        },
#    CombinationCut = ("(ADAMASS('Lambda0') < 50.*MeV) & (ADOCACHI2CUT(30, ''))"),
#    MotherCut      =  ("(ADMASS('Lambda0') < 35.*MeV) & (VFASPF(VCHI2) < 30.)"),
#)
#
#LambdaDD = SimpleSelection (
#    'LambdaDD',
#    CombineParticles,
#    [downpions,downprotons],
#    DecayDescriptor = "[Lambda0 -> p+ pi-]cc",
#    DaughtersCuts = {
#        "p+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 4.)",
#        "pi+" : "(P > 2.*GeV) & (MIPCHI2DV(PRIMARY) > 4.)"
#        },
#    CombinationCut = ("(ADAMASS('Lambda0') < 80.*MeV) & (ADOCACHI2CUT(25, ''))"),
#    MotherCut      =  ("(ADMASS('Lambda0') < 64.*MeV) & (VFASPF(VCHI2) < 25.)"),
#)

##THESE MATCH Hlt/Hlt2SharedParticles/python/Hlt2SharedParticles
#KSLL = SimpleSelection (
#    'KSLL',
#    CombineParticles,
#    [longpions],
#    DecayDescriptor = "KS0 -> pi+ pi-",
#    DaughtersCuts = {
#        "pi+" : "(TRCHI2DOF<3.)& (MIPCHI2DV(PRIMARY)>36)"
#        },
#    CombinationCut = ("(ADAMASS('KS0') < 50.*MeV)"),
#    MotherCut      =  ("(ADMASS('KS0') < 35.*MeV) & (VFASPF(VCHI2PDOF)<30) & (BPVLTIME() > 2.0*ps)"),
#)
#
#KSDD = SimpleSelection (
#    'KSDD',
#    CombineParticles,
#    [downpions],
#    DecayDescriptor = "KS0 -> pi+ pi-",
#    DaughtersCuts = {
#        "pi+" : "(TRCHI2DOF<4) & (P>3000*MeV) & (PT > 175.*MeV)"
#        },
#    CombinationCut = ("(ADAMASS('KS0') < 80.*MeV)"),
#    MotherCut      =  ("(ADMASS('KS0') < 64.*MeV) & (VFASPF(VCHI2PDOF)<30)  & (BPVVDZ > 400.0*mm)"),
#)
#
#LambdaLL = SimpleSelection (
#    'LambdaLL',
#    CombineParticles,
#    [longpions,longprotons],
#    DecayDescriptor = "[Lambda0 -> p+ pi-]cc",
#    DaughtersCuts = {
#        "p+" : "BPVVALID() & (TRCHI2DOF<4) & (MIPCHI2DV(PRIMARY)>36)",
#        "pi+" : "BPVVALID() & (TRCHI2DOF<4) & (MIPCHI2DV(PRIMARY)>36)"
#        },
#    CombinationCut = ("(ADAMASS('Lambda0') < 50.*MeV)"),
#    MotherCut      =  ("(ADMASS('Lambda0') < 20.*MeV) & (CHI2VXNDOF<30) & (BPVLTIME() > 2.0*ps)"),
#)
#
#LambdaDD = SimpleSelection (
#    'LambdaDD',
#    CombineParticles,
#    [downpions,downprotons],
#    DecayDescriptor = "[Lambda0 -> p+ pi-]cc",
#    DaughtersCuts = {
#        "p+" : "(TRCHI2DOF<4)& (P>3000*MeV)  & (PT>175.*MeV)",
#        "pi+" : "(TRCHI2DOF<4)& (P>3000*MeV)  & (PT>175.*MeV)"
#        },
#    CombinationCut = ("(ADAMASS('Lambda0') < 80.*MeV)"),
#    MotherCut      =  ("in_range(1095*MeV, M, 1140*MeV) & (CHI2VXNDOF<30) & (BPVVDZ > 400.0*mm)"),
#)

combKS = MergedSelection (
    'combKS',
    #RequiredSelections =  [KSLL, KSDD]
    RequiredSelections =  [looseKSLL, looseKSDD]
    )

combLambda = MergedSelection (
    'combLambda',
    #RequiredSelections =  [LambdaLL, LambdaDD]
    RequiredSelections =  [vlooseLambdaLL, looseLambdaDD]
#    RequiredSelections =  [looseLambdaLL, vlooseLambdaLL, looseLambdaDD]
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


#recSVs = SimpleSelection (
#    'recSVs',
#    CombineParticles,
#    [ kaons ],
#    DecayDescriptor = "K*(892)0 -> K+ K-",
##    CombinationCut  = "ALL",
#    MotherCut       = "(VFASPF(VCHI2PDOF) < 10.) & (BPVVDCHI2 > 25.)",
#)


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

charmProtons = SimpleSelection (
    'charmProtons'         ,
    FilterDesktop   ,
    [ looseprotons ]    ,
    DecayDescriptor = "[p+]cc",
    Code = (#"(PIDp - PIDpi > 5) & "
            "(PT>200*MeV)"
            "& (TRGHOSTPROB<0.3)"
            "& (MIPCHI2DV(PRIMARY) > 4)")
    )

dcD0 = { }
for child in ['pi+','K+'] :
    dcD0[child] = "(PT > 250*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"#16)"
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
             "& (ANUM(PT > 1000*MeV) > 0)" #\
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 10) > 1)" \
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 50) > 0)"

parentcutsDp = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVDIRA > 0.9994 )"
#               "& (BPVLTIME() > 0.150*ps )" \

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
             "& (ANUM(PT > 1000*MeV) > 0)" #\
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 10) > 1)" \
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 50) > 0)"

parentcutsDs = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVDIRA > 0.9994 )"
#               "& (BPVLTIME() > 0.150*ps )" \

recDs = SimpleSelection (
    'recDs',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D_s+ -> K- K+ pi+]cc",
    DaughtersCuts   = dcDs,
    CombinationCut = (combcutsDs),
    MotherCut      =  (parentcutsDs),
)

dcLc = { }
for child in ['pi+','K+','p+'] :
    dcLc[child] = "(PT > 200*MeV)" \
                  "& (P > 2*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"
#                  "& (TRCHI2 < 3)" \

combcutsLc = "in_range(2256*MeV,  AM, 2316*MeV)" \
             "& (ANUM(PT > 400*MeV) > 1)" \
             "& (ANUM(PT > 1000*MeV) > 0)" #\
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 10) > 1)" \
#             "& (ANUM(MIPCHI2DV(PRIMARY) > 50) > 0)"

parentcutsLc = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVDIRA > 0.9994 )"

recLc = SimpleSelection (
    'recLc',
    CombineParticles,
    [charmProtons,charmKaons,charmPions],
    DecayDescriptor = "[Lambda_c+ -> p+ K- pi+]cc",
    DaughtersCuts   = dcLc,
    CombinationCut = (combcutsLc),
    MotherCut      =  (parentcutsLc),
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
Lc_seq = SelectionSequence('Lc_Seq', TopSelection=recLc)
D02K3pi_seq = SelectionSequence('D2K3pi0_Seq', TopSelection=recD02K3pi)

##########################

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = True
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
#DaVinci().appendToMainSequence([genPF, genJB, recPF, recJB])
#DaVinci().appendToMainSequence([recPF, recJB])
DaVinci().appendToMainSequence([recPF, recJB, recSVs_seq.sequence(), recMus_seq.sequence(), D0_seq.sequence(), Dp_seq.sequence(), Ds_seq.sequence(), Lc_seq.sequence(), D02K3pi_seq.sequence()])
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
        self.l0Tool = toolSvc.create(
                 'L0TriggerTisTos',
                 interface = 'ITriggerTisTos')
        self.hlt1Tool = toolSvc.create(
                 'TriggerTisTos/Hlt1TriggerTisTos',
                 interface = 'ITriggerTisTos')
        self.hlt2Tool = toolSvc.create(
                 'TriggerTisTos/Hlt2TriggerTisTos',
                 interface = 'ITriggerTisTos')
        #self.hlt2Tool.setProperty('HltDecReportsLocation', '/Event/Hlt2/DecReports')
        #self.hlt2Tool.setProperty('HltSelReportsLocation', '/Event/Hlt2/SelReports')
        #print self.hlt2Tool.HltDecReportsLocation, self.hlt2Tool.HltSelReportsLocation
        #from Configurables import ToolSvc, TriggerTisTos
        #for stage in ('Hlt1', 'Hlt2', 'Strip/Phys'):
        #    toolSvc.addTool(TriggerTisTos, stage + "TriggerTisTos")
        #    tool = getattr(toolSvc, stage + "TriggerTisTos")
        #    tool.HltDecReportsLocation = '/Event/' + stage + '/DecReports'
        #    tool.HltSelReportsLocation = '/Event/' + stage + '/SelReports'
        #from Configurables import ToolSvc, TriggerTisTos
        #ToolSvc().addTool(TriggerTisTos, "Hlt2TriggerTisTos")
        #tool = getattr(ToolSvc(), "Hlt2TriggerTisTos")
        #tool.HltDecReportsLocation = '/Event/Hlt2/DecReports'
        #tool.HltSelReportsLocation = '/Event/Hlt2/SelReports'
        #self.hlt2Tool = tool

        self.stable = [11,-11,13,-13,211,-211,321,-321,2212,-2212,2112,-2112,22,111,310,130,311,-311]
        self.Ds = [411,-411,421,-421,431,-431,4122,-4122]
        self.tes     = tes
        self.saved   = {}
        self.ntuple  = OrderedDict()
        self.tfile   = ROOT.TFile('output.root', 'RECREATE')
        self.ttree   = ROOT.TTree('data', 'data')
        self.vrs     = {}
        mom = ['px', 'py', 'pz', 'e']
        pos = ['x', 'y', 'z']
        cov = ['dx', 'dy', 'dz', 'chi2', 'ndof']
        self.init('gen', ['idx_pvr', 'idx_jet', 'idx_prnt', 'pid', 'q'] + mom + pos + ['prnt_pid', 'res_pid', 'from_sig'])#, 'prnt_key', 'key'])
        self.init('pvr', pos + cov)
        self.init('svr', ['idx_pvr', 'idx_jet'] + [
                'idx_trk%i' % i for i in range(0, 10)] + 
                  mom + pos + ['m', 'm_cor', 'pt', 'fd_min', 'fd_chi2', 'chi2', 'ip_chi2_sum', 'abs_q_sum', 'tau', 'ntrk', 'ntrk_jet', 'jet_dr', 'jet_pt', 'pass', 'bdt0', 'bdt1'])
        self.init('jet', ['idx_pvr', 'ntrk', 'nneu'] + mom) #+ ['idx_trk%i' % i for i in range(0, 40)] + [
        #'%s_%s' % (i,j) for i in ['DiJet','DiJetSV','DiJetSVSV','DiJetSVMu','DiJetMuMu'] for j in ['Dec','pRatio1','pRatio2','neuMultRatio1','neuMultRatio2','chrgMultRatio1','chrgMultRatio2','dR1','dR2'] 
        #])
        self.init('trk', ['idx_gen', 'idx_pvr', 'idx_jet'] + mom +
                  ['pid', 'q', 'ip', 'ip_chi2', 'pnn_e', 'pnn_mu', 'pnn_pi',
                   'pnn_k', 'pnn_p', 'pnn_ghost', 'ecal', 'hcal', 'prb_ghost', 'type', 'is_mu',
                   'vid', 'x', 'y', 'z'])# + ['vhit%i' % i for i in range(0, 61)])# + [
                   #'vz%i' % i for i in range(0, 61)])
        self.init('neu', ['idx_gen', 'idx_jet'] + mom + ['pid'])
        self.init('d0', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 2)]) 
        self.init('dp', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)]) 
        self.init('ds', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)]) 
        self.init('lc', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 3)]) 
        self.init('d02k3pi', ['idx_pvr','idx_jet'] + mom + pos + ['m', 'ip', 'ip_chi2', 'vtx_chi2', 'vtx_ndof', 'fd', 'fd_chi2', 'tau', 'tau_err', 'tau_chi2', 'ntrk_jet'] + ['idx_trk%i' % i for i in range(0, 4)]) 
        self.init('evt', ['dec'] + ['%s_%s' % (k,l) for k in ['j1', 'j2'] for l in ['idx','dR','nsv','nmu','ntrk','nneu'] + mom ] )
        self.ntuple['evt_pvr_n'] = array.array('d', [-1])
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
    #def getFirstState(self, obj):
    #    if not obj: return
    #    print 0
    #    vids = []
    #    for vid in obj.lhcbIDs():
    #        if vid.isVelo():
    #            vrs['vid'] *= float(vid.veloID().channelID())/1000000.0
    #            vids += [(self.detTool.sensor(vid.veloID()).z(), 
    #                      vid.veloID().channelID())]
    #    print 1
    #    vids.sort()
    #    print 2
    #    sta = LHCB.StateVector()
    #    print 3
    #    if len(vids) > 0: self.trkTool.propagate(obj, vids[0][0], sta, pid)
    #    print 4
    #    return sta
    def getJetMult(self, summary, nNeu=0, nChrg=0):
        if summary.summarizedObjectCLID() == 801:
            #if not ROOT.TMath.Abs(summary.numericalInfo()["0#Particle.particleID.pid"]) in [22,111,211,321,2212]:
            #    if not ROOT.TMath.Abs(summary.numericalInfo()["0#Particle.particleID.pid"]) in [98,111,310,3122]:
            #        print "A", summary.numericalInfo()["0#Particle.particleID.pid"]
            #        self.processHltSummary(summary)
            #    for ss in summary.substructure():
            #        nNeu,nChrg = self.getJetMult(ss, nNeu, nChrg)
            if 801 in [ss.summarizedObjectCLID() for ss in summary.substructure()]:##do any of the substructures contain a particle class
            #if 801 in [summary.substructure()[i].summarizedObjectCLID() for i in range(summary.substructure().size())]:
                for ss in summary.substructure():
                    #if not ROOT.TMath.Abs(summary.numericalInfo()["0#Particle.particleID.pid"]) in [98,111,310,3122]:
                    #    print "A", summary.numericalInfo()["0#Particle.particleID.pid"]
                    #    self.processHltSummary(summary)
                    nNeu,nChrg = self.getJetMult(ss, nNeu, nChrg)
            else:
                if summary.numericalInfo()["0#Particle.particleID.pid"] in [22,111]:
                    nNeu+=1
                else:
                    #if not ROOT.TMath.Abs(summary.numericalInfo()["0#Particle.particleID.pid"]) in [211,321,2212]:
                    #    print "B", summary.numericalInfo()["0#Particle.particleID.pid"]
                    nChrg+=1
        return nNeu,nChrg
    def processHltSummary(self, summary, indent=""):
        #print summary.summarizedObjectCLID()
        if summary.summarizedObjectCLID() == 801:
            print indent+str(summary.numericalInfo()["0#Particle.particleID.pid"])
        #else:
        #    print summary.summarizedObjectCLID()
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
#    def fillJetTISTOS(self, obj, vrs, dtrs, pre = 'jet'):
#        if not obj: return
#        #print "foo"
#        for dec in ['DiJet','DiJetSV','DiJetSVSV','DiJetSVMu','DiJetMuMu']:
#            summaries = self.hlt2Tool.hltObjectSummaries(obj,'Hlt2Jets'+dec+'Decision')
#            if summaries.size()>0:
#                #print 1
#                vrs[dec+'_Dec'] = 1
#            else:
#                #print 0
#                vrs[dec+'_Dec'] = 0
#                continue
#            #print "foo"
#            nNeu=0
#            nChrg=0
#            for dtr in dtrs:
#                try:
#                    dtr = dtr.target()
#                except:
#                    pass
#                if not dtr.proto() or not dtr.proto().track(): nNeu+=1
#                else: nChrg+=1
#            #print nNeu, nChrg
#            #print dec, vrs[dec]
#            for summary in summaries:
#                if summary.numericalInfo()["0#Particle.particleID.pid"] != 97: return
#                jet1 = summary.substructure().at(0)
#                jet2 = summary.substructure().at(1)
#                if jet1.numericalInfo()["0#Particle.particleID.pid"] != 98 or jet2.numericalInfo()["0#Particle.particleID.pid"] != 98: return
#                j1tx = jet1.numericalInfo()["5#Particle.slopes.x"]
#                j1ty = jet1.numericalInfo()["6#Particle.slopes.y"]
#                j1p  = 1./jet1.numericalInfo()["7#Particle.1/p"]
#                j1Denom = (1. + j1tx**2. + j1ty**2.)**.5 
#                j1px = j1p * j1tx / j1Denom
#                j1py = j1p * j1ty / j1Denom
#                j1pz = j1p  / j1Denom
#                j1pt = j1p * (j1tx**2. + j1ty**2.)**.5 / j1Denom
#                j2tx = jet2.numericalInfo()["5#Particle.slopes.x"]
#                j2ty = jet2.numericalInfo()["6#Particle.slopes.y"]
#                j2p  = 1./jet2.numericalInfo()["7#Particle.1/p"]
#                j2Denom = (1. + j2tx**2. + j2ty**2.)**.5 
#                j2px = j2p * j2tx / j2Denom
#                j2py = j2p * j2ty / j2Denom
#                j2pz = j2p  / j2Denom
#                j2pt = j2p * (j2tx**2. + j2ty**2.)**.5 / j2Denom
#                nNeu1=0
#                nChrg1=0
#                #self.processHltSummary(jet1)
#                nNeu1,nChrg1 = self.getJetMult(jet1)
#                nNeu2,nChrg2 = self.getJetMult(jet2)
#                #print nNeu1, nChrg1
#                #matchFact1=0.
#                #matchFact2=0.
#                #if obj.p() > j1p:
#                #    matchFact1 += j1p/obj.p()
#                #else :
#                #    matchFact1 += obj.p()/j1p
#                ##if obj.pt() > j1pt:
#                ##    matchFact1 += 0.5*j1pt/obj.pt()
#                ##else :
#                ##    matchFact1 += 0.5*obj.pt()/j1pt
#                #if obj.p() > j2p:
#                #    matchFact2 += j2p/obj.p()
#                #else :
#                #    matchFact2 += obj.p()/j2p
#                ##if obj.pt() > j2pt:
#                ##    matchFact2 += 0.5*j2pt/obj.pt()
#                ##else :
#                ##    matchFact2 += 0.5*obj.pt()/j2pt
#                vrs[dec+"_pRatio1"] = obj.p()/j1p #matchFact1
#                vrs[dec+"_pRatio2"] = obj.p()/j2p #matchFact2
#                vrs[dec+"_neuMultRatio1"] = nNeu/float(nNeu1)
#                vrs[dec+"_neuMultRatio2"] = nNeu/float(nNeu2)
#                vrs[dec+"_chrgMultRatio1"] = nChrg/float(nChrg1)
#                vrs[dec+"_chrgMultRatio2"] = nChrg/float(nChrg2)
#                vrs[dec+"_dR1"] = ROOT.TLorentzVector(obj.momentum().px(),obj.momentum().py(),obj.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j1px,j1py,j1pz,0.))
#                vrs[dec+"_dR2"] = ROOT.TLorentzVector(obj.momentum().px(),obj.momentum().py(),obj.momentum().pz(),0.).DeltaR(ROOT.TLorentzVector(j2px,j2py,j2pz,0.))
#                #print "bar"
#                #if obj.p() > j1p*0.9 and obj.p() < j1p*1.1 and obj.pt() > j1pt*0.9 and obj.pt() < j1pt*1.1:
#                #    vrs["match1"] = 1
#                #if obj.p() > j2p*0.9 and obj.p() < j2p*1.1 and obj.pt() > j2pt*0.9 and obj.pt() < j2pt*1.1:
#                #    vrs["match2"] = 1
#                #print obj.p(), j1p, j2p
#                #print obj.pt(), j1pt, j2pt
#                #print obj.p(), 1./summary.substructure().at(0).numericalInfo()["7#Particle.1/p"], 1./summary.substructure().at(1).numericalInfo()["7#Particle.1/p"]
#                #print obj.momentum().Px()/obj.momentum().Pz(), summary.substructure().at(0).numericalInfo()["5#Particle.slopes.x"], summary.substructure().at(0).numericalInfo()["5#Particle.slopes.x"]
#                #print obj.momentum().Py()/obj.momentum().Pz(), summary.substructure().at(0).numericalInfo()["6#Particle.slopes.y"], summary.substructure().at(0).numericalInfo()["6#Particle.slopes.y"]
#                #print summary.summarizedObjectCLID()
#                #self.processHltSummary(summary)
#                #print summary.summarizedObjectCLID()
#                #sumobj = summary.summarizedObject()
#                #stash.append(sumobj)
#                #print sumobj.tes()
#                #print sumobj.tesLocation()
#                #print dir(sumobj)
#                #print sumobj.momentum().Px() #, sumobj.momentum.Py(), sumobj.momentum().Pz(), sumobj.momentum.Pt()
#        #summaries = self.hlt2Tool.hltObjectSummaries(".*")
#        #print summaries

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
                    #key = self.key(dau)
                    #if key in self.saved['trk']:
                    #    trks.append(self.saved['trk'][key])
                    #else :
                    #    trks.append(self.addTrk(dau))
                    #trks.append(self.addTrk(dau))
                    key = self.key(dau)
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

    def addTrigger(self, pre = 'evt'):
        #self.init('evt', ['%s_%s' % (i,j) for i in ['DiJet','DiJetSV','DiJetSVSV','DiJetSVMu','DiJetMuMu'] for j in ['Dec', 'jet1', 'jet2'] + ['%s_%s' % (k,l) for k in ['j1', 'j2'] for l in ['nsv','nmu','ntrk','nneu'] + mom ] ] )
        for idec, dec in enumerate(['DiJet','DiJetSV','DiJetSVSV','DiJetSVMu','DiJetMuMu']):
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
                    jets = tes[recJB.Output]
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
                    svs = tes[recSVs.outputLocation()]
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
                    mus = tes[recMus.algorithm().Output]
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
            #if parent.particleID().pid() == 431:
            #    print parIdx, obj.particleID().abspid()##TODO
#            if parKey: vrs['prnt_key'] = int(parKey[2]*10)*1e12 + int(parKey[1]*10)*1e6 + int(parKey[0]*10)
#        vrs['key'] = int(key[2]*10)*1e12 + int(key[1]*10)*1e6 + int(key[0]*10)
        self.saved[pre][key] = idx
        self.fill(pre, vrs = vrs)
        #also print immediate decay products of heavy hadrons
        pid = obj.particleID()
        if pid.isHadron() and (pid.hasCharm() or pid.hasBottom()):
            if pid.abspid() in self.Ds:
                #if pid.abspid() == 431:
                #    print idx, self.finalChildren(obj)##TODO
                for part in self.finalChildren(obj)[0]:
                    self.addGen(part,par=obj)
            else:
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
                #stash.append(trks[-1])
                #vrs['idx_trk%i' % (len(trks) - 1)] = trks[-1][0]
        vrs['ntrk'] = len(trks)
        vrs['nneu'] = nneu
        #self.fillJetTISTOS(obj, vrs, obj.daughters())
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
        ntuple.addTrigger()
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
            ntuple.addDHad(d0,"d02k3pi"); fill = True;
    except: pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
