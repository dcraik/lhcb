# Test data
from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/eos/lhcb/grid/prod/lhcb/LHCb/Protonion16/IFT.DST/00076144/0000/00076144_00000038_1.ift.dst'], clear = True)
#IOHelper('ROOT').inputFiles(['/tmp/dcraik/00076144_00000038_1.ift.dst'], clear = True)

## import DaVinci
from Configurables import DaVinci
from Configurables import CondDB

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = """
    HLT_PASS_RE ( 'StrippingHeavyIonTopologyLowActivityLineDecision' )
    """
)

dv = DaVinci (
    DataType        = '2016'              ,
    InputType       = 'DST'              ,
    TupleFile       = 'Tuples.root' , ## IMPORTANT
    HistogramFile   = 'Histos.root' , ## IMPORTANT
    Simulation      = False,
    Lumi            = True,
    EventPreFilters = fltrs.filters ('Filters'),
    #CondDBtag       = 'cond-20180625-1',
    #DDDBtag         = 'dddb-20180726-3',
    #DQFLAGStag      = 'dq-20170627',
    EvtMax          = -1 ##TODO
    #EvtMax          = 1000 ##TODO
    )
db = CondDB( LatestGlobalTagByDataType = '2016' )
db.LocalTags["DQFLAGS"] = [ "herschel-20161018" ]


from StandardParticles import StdAllNoPIDsKaons    as loosekaons

from PhysSelPython.Wrappers import SimpleSelection, MergedSelection, DataOnDemand, Selection
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles
#from Configurables import FilterDesktop

kaons = SimpleSelection (
    'kaons'         ,
    FilterDesktop   ,
    [ loosekaons ]    ,
    DecayDescriptor = "[K+]cc",
    Code = ("(PIDK - PIDpi > 5)"
            "& (TRGHOSTPROB<0.3)")
    )

dcDiK = { }
dcDiK['K+'] = "(TRGHOSTPROB < 0.3)"
combcutsDiK = "(in_range( 0, AM, 1200 ))"
parentcutsDiK = "ALL"

phis = SimpleSelection (
    'phis',
    CombineParticles,
    [kaons],
    DecayDescriptor = "phi(1020) -> K+ K-",
    DaughtersCuts   = dcDiK,
    CombinationCut = (combcutsDiK),
    MotherCut      =  (parentcutsDiK),
)


#from PhysSelPython.Wrappers import SelectionSequence
#phi_seq = SelectionSequence('phi_seq', TopSelection=phis)

## use selection as input for DecayTreeTuple
from GaudiConfUtils.ConfigurableGenerators import DecayTreeTuple
tup = SimpleSelection (
    'Tuple'    ,
    DecayTreeTuple ,
    [ phis ]    ,
    Decay    = "phi(1020) -> ^K+ ^K-" ,
    Branches = {
    "phi"               : "phi(1020)->  K+  K-",
    "K+"                : "phi(1020)-> ^K+  K-",
    "K-"                : "phi(1020)->  K+ ^K-",
    }
    )

## extra configurataion of tuple algorithm (if needed)
tuples = [tup]

import DecayTreeTuple.Configuration
from Configurables import TupleToolMCTruth, MCTupleToolKinematic, MCTupleToolHierarchy, TupleToolMCBackgroundInfo, TupleToolTrackInfo, TupleToolRICHPid, TupleToolMuonPid, TupleToolPid, TupleToolGeometry

for tup in tuples:
    tupalg = tup.algorithm()
    #tupalg.OutputLevel = MSG.DEBUG

    tupalg.ToolList +=  [
                  "TupleToolGeometry"
                , "TupleToolKinematic"
                , "TupleToolPrimaries"
                , "TupleToolEventInfo"
                , "TupleToolTrackInfo"
    #            , "TupleToolAngles"
                , "TupleToolPid"
#                , "TupleToolPropertime"
    #            , "TupleToolTrigger"
    #            , "TupleToolRICHPid"
    #            , "TupleToolMuonPid"
                , "TupleToolProtoPData"
    #            , "TupleToolMCBackgroundInfo"
                ]
    #
    #tool=tupalg.addTupleTool('TupleToolMCTruth')
    #tool.ToolList =  [
    #        "MCTupleToolKinematic"
    #      , "MCTupleToolHierarchy"
    #]

    #LoKi_MOTHER = LoKi__Hybrid__MCTupleTool("LoKi_MOTHER")
    #LoKi_MOTHER.Variables = {
    #    }
    #
    #tool.addTupleTool(LoKi_MOTHER)

    from Configurables import TupleToolTISTOS
    tool = tupalg.addTupleTool(TupleToolTISTOS)
    tool.Verbose=True
    tool.TriggerList = [ "L0DiHadron,lowMultDecision",
                         "L0HRCDiHadron,lowMultDecision",
                         "Hlt1LowMultPassThroughDecision"
                         "Hlt1LowMultDecision",
                         "Hlt1LowMultHerschelDecision",
                         "Hlt1LowMultVeloCut_HadronsDecision",
                         "Hlt1LowMultVeloAndHerschel_HadronsDecision",
                         "Hlt2LowMultLMR2HHDecision",
                         "Hlt2LowMultLMR2HHWSDecision"
                         ]

    #from Configurables import TupleToolL0Calo
    #tool = tupalg.addTupleTool(TupleToolL0Calo,name="L0Hcal")
    #tool = tupalg.addTupleTool(TupleToolL0Calo,name="L0Ecal")
    #tool.WhichCalo="ECAL"

    # RecoStats to filling SpdMult, etc
    from Configurables import TupleToolRecoStats
    tool = tupalg.addTupleTool(TupleToolRecoStats)
    #tool.Verbose=True

    # event tuple
    from Configurables import LoKi__Hybrid__EvtTupleTool
    tool = tupalg.addTupleTool( LoKi__Hybrid__EvtTupleTool, name = 'LoKi_EvtTuple')
    tool.VOID_Variables = {
        "LoKi_nPVs"                : "CONTAINS('Rec/Vertex/Primary')",
        "LoKi_nSpdMult"            : "CONTAINS('Raw/Spd/Digits')",
        "LoKi_nVeloClusters"       : "CONTAINS('Raw/Velo/Clusters')",
        "LoKi_nVeloLiteClusters"   : "CONTAINS('Raw/Velo/LiteClusters')",
        "LoKi_nITClusters"         : "CONTAINS('Raw/IT/Clusters')",
        "LoKi_nTTClusters"         : "CONTAINS('Raw/TT/Clusters')",
        "LoKi_nOThits"             : "CONTAINS('Raw/OT/Times')",
        "LoKi_nNeutrals"           : "CONTAINS('Rec/ProtoP/Neutrals')",
        "LoKi_nCharged"            : "CONTAINS('Rec/ProtoP/Charged')"
        }

    #HRC tool
    tool=tupalg.addTupleTool("TupleToolHerschel")

## for debugging:
from SelPy.graph import graph
for tup in tuples:
  graph( tup, format = 'png' )

from PhysSelPython.Wrappers                import SelectionSequence
#seq  = SelectionSequence ( 'LbKSSEQ' , tupLbKS )
#dv.UserAlgorithms += [ seq.sequence() ]
#seq  = SelectionSequence ( 'LbKstSEQ' , tupLbKst )
#dv.UserAlgorithms += [ seq.sequence() ]
#seq  = SelectionSequence ( 'LbKSJpsipKSEQ' , tupLbKS_jpsipk )
#dv.UserAlgorithms += [ seq.sequence() ]
#seq  = SelectionSequence ( 'LbKstJpsipKSEQ' , tupLbKst_jpsipk )
#dv.UserAlgorithms += [ seq.sequence() ]

from Configurables import HCRawBankDecoder
decoder = HCRawBankDecoder()
DaVinci().UserAlgorithms += [ decoder ];

for tup in tuples:
    seq  = SelectionSequence ( tup.name()+'_SEQ' , tup )
    dv.UserAlgorithms += [ seq.sequence() ]
