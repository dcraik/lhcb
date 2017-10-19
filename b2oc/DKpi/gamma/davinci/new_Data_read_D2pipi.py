########################################################################
#
#  an options file to read data and create an ntuple.
#
########################################################################


from Gaudi.Configuration import *
from Configurables import DecayTreeTuple, TupleToolTrigger, TupleToolDecay, SubstitutePID
from Configurables import LoKi__Hybrid__TupleTool, TupleToolDecayTreeFitter, TupleToolTISTOS
from Configurables import LoKi__Hybrid__EvtTupleTool, TupleToolEventInfo, TupleToolTrackIsolation
from Configurables import TupleToolMCTruth, TupleToolTrackInfo, TupleToolRICHPid, TupleToolMuonPid, TupleToolPid, TupleToolGeometry
from Configurables import DaVinci
from Configurables import PrintDecayTree
from DecayTreeTuple.Configuration import *

from Configurables import EventNodeKiller
NodeKiller = EventNodeKiller('NodeKiller')
NodeKiller.Nodes = ['DAQ', 'pRec']

#teslocation = "/Event/AllStreams/" #change to this for MC

teslocation = "/Event/Bhadron/" #change to this for Data

#
#Momentum scale correction
#

from Configurables import TrackScaleState

scaler = TrackScaleState( 'Scaler' , RootInTES = teslocation ) 
#    ScaleSlope = 1 - 1.0e-3 ,
#    DeltaSlope = 0.9e-3 ,
# )
#scaler.DeltaScale = 0.3e-3

from Configurables import CondDB
CondDB().UseLatestTags = ["2012"] 

myTriggerList = [ "L0HadronDecision","L0MuonDecision","L0DiMuonDecision","L0ElectronDecision","L0PhotonDecision",
                  "L0GlobalPi0Decision", "L0LocalPi0Decision",
                  "Hlt1TrackAllL0Decision", "Hlt1GlobalDecision", "Hlt2GlobalDecision",
                  "Hlt2Topo2BodySimpleDecision", "Hlt2Topo3BodySimpleDecision", "Hlt2Topo4BodySimpleDecision",
                  "Hlt2Topo2BodyBBDTDecision",   "Hlt2Topo3BodyBBDTDecision",
                  "Hlt2Topo4BodyBBDTDecision", "Hlt1TrackMuonDecision", "Hlt1DiMuonLowMassDecision",
                  "Hlt1DiMuonHighMassDecision", "Hlt2DiMuonDetachedJPsiDecision", "Hlt2DiMuonJPsiDecision"
                  ]

#### -- B mass cut -- ####

from Configurables import FilterDesktop, CombineParticles
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence

BMassSel = AutomaticData(Location = '/Event/Bhadron/Phys/B02D0KPiD2HHBeauty2CharmLine/Particles')


_bmassFilter = FilterDesktop('bmassFilter', Code = '(M<6000.*MeV) & (BPVIPCHI2()<15) & (BPVDIRA>0.9999)')

BMassFilterSel = Selection(name = 'BMassFilterSel',
                           Algorithm = _bmassFilter,
                           OutputBranch = teslocation+'Phys',
                           RequiredSelections = [ BMassSel ])

BMassSeq = SelectionSequence('SeqBMass',
                            TopSelection = BMassFilterSel)

BMassseq = BMassSeq.sequence()

### Tuple part
tuple = DecayTreeTuple()

tuple.RootInTES = teslocation
tuple.Inputs = [ 'Phys/BMassFilterSel/Particles' ]
#tuple.Inputs = [ '/Event/Bhadron/Phys/B02D0KPiD2HHBeauty2CharmLine/Particles' ]

tuple.ToolList +=  [
              "TupleToolGeometry"
            , "TupleToolRecoStats"
            , "TupleToolKinematic"
            , "TupleToolPrimaries"
            , "TupleToolEventInfo"
            , "TupleToolTrackInfo"
            , "TupleToolAngles"
            , "TupleToolPid"
            , "TupleToolPropertime"
            , "TupleToolTrigger"
            , "TupleToolRICHPid"
            , "TupleToolMuonPid"
            , "TupleToolProtoPData"
            ]

tuple.addTool(TupleToolGeometry,name="TupleToolGeometry")
tuple.TupleToolGeometry.Verbose = True

tuple.addTool(TupleToolTrackInfo,name="TupleToolTrackInfo")
tuple.TupleToolTrackInfo.Verbose = True

tuple.addTool(TupleToolPid,name="TupleToolPid")
tuple.TupleToolPid.Verbose = True

#tuple.UseLabXSyntax = True

#tuple.Decay = " [B0 -> ^[ D0 -> ^K- ^pi+ ]CC ^[ K*(892)0 -> ^K+ ^pi- ]CC]CC "
#tuple.Decay = " [B0 -> ^(D0 -> ^K- ^pi+) ^(K*(892)0 -> ^K+ ^pi-) ]CC "
tuple.Decay = "[B0 -> ^( D0 -> ^pi+ ^pi- ) [^( K*(892)0 -> ^K+ ^pi-)]CC ]CC"
tuple.addBranches ({
    "B"    : " [B0 ->  ( D0 ->  pi+  pi- ) [ ( K*(892)0 ->  K+  pi-)]CC ]CC ",
    "D"    : " [B0 -> ^( D0 ->  pi+  pi- ) [ ( K*(892)0 ->  K+  pi-)]CC ]CC ",
    "D0p"  : " [B0 ->  ( D0 -> ^pi+  pi- ) [ ( K*(892)0 ->  K+  pi-)]CC ]CC ",
    "D0m"  : " [B0 ->  ( D0 ->  pi+ ^pi- ) [ ( K*(892)0 ->  K+  pi-)]CC ]CC ",      
    "K"    : " [B0 ->  ( D0 ->  pi+  pi- ) [ ( K*(892)0 -> ^K+  pi-)]CC ]CC ",       
    "pi"   : " [B0 ->  ( D0 ->  pi+  pi- ) [ ( K*(892)0 ->  K+ ^pi-)]CC ]CC ",
    "Kst"  : " [B0 ->  ( D0 ->  pi+  pi- ) [^( K*(892)0 ->  K+  pi-)]CC ]CC "
})

tuple.addTupleTool("TupleToolTISTOS/TupleToolTISTOS")
tuple.TupleToolTISTOS.VerboseL0 = True
tuple.TupleToolTISTOS.VerboseHlt1 = True
tuple.TupleToolTISTOS.VerboseHlt2 = True
tuple.TupleToolTISTOS.TriggerList = myTriggerList

#Add Loki variables
LoKi_DOCA = LoKi__Hybrid__TupleTool("LoKi_DOCA")
LoKi_DOCA.Variables = {
    "MAXDOCA" : "LoKi.Particles.PFunA(AMAXDOCA('LoKi::TrgDistanceCalculator'))",
    "MINDOCA" : "LoKi.Particles.PFunA(AMINDOCA('LoKi::TrgDistanceCalculator'))",
    "DOCA12"  : "LoKi.Particles.PFunA(ADOCA(1,2))"
    }

LoKi_LT = LoKi__Hybrid__TupleTool("LoKi_LT")
LoKi_LT.Variables = {
    "LT_FITCHI2" : "BPVLTFITCHI2()",
    "LT" : "BPVLTIME()",
    "LT_CHI2" : "BPVLTCHI2()"
    }

LoKi_Iso =  LoKi__Hybrid__TupleTool( 'LoKi_Iso')
LoKi_Iso.Variables = {
    "LOKI_CONE1ANGLE" : "INFO(9000, 0.)",
    "LOKI_CONE1MULT"  : "INFO(9000+1, 0.)",
    "LOKI_CONE1PASY"  : "INFO(9000+10, 0.)",
    "LOKI_CONE1PTASY" : "INFO(9000+11, 0.)",
    }

LoKi_constrain = LoKi__Hybrid__TupleTool("LoKi_constrain")
LoKi_constrain.Preambulo += ['from LoKiTracks.decorators import *']
LoKi_constrain.Preambulo = [
]

LoKi_constrain.Variables = {
    'D0_B_CM_MM' : "DTF_FUN(MM, True,strings(['D0']))",
    'D0_B_CM_M' : "DTF_FUN(M, True,strings(['D0']))",
    'D0_B_CM_PX' : "DTF_FUN(PX, True,strings(['D0']))",
    'D0_B_CM_PY' : "DTF_FUN(PY, True,strings(['D0']))",
    'D0_B_CM_PZ' : "DTF_FUN(PZ, True,strings(['D0']))",
    'D0_B_CM_E' : "DTF_FUN(E, True,strings(['D0']))",
    'D0_B_CM_MINIPCHI2' : "DTF_FUN(MIPCHI2DV(), True,strings(['D0']))",
    'D0_B_CM_DIRA_OWNPV' : "DTF_FUN(BPVDIRA, True,strings(['D0']))",
    'D0_B_CM_ENDVERTEX_CHI2' : "DTF_FUN(VFASPF(VCHI2), True,strings(['D0']))",
    'D0_B_CM_ENDVERTEX_NDOF' : "DTF_FUN(VFASPF(VDOF), True,strings(['D0']))",
    'D0_B_CM_FDCHI2_OWNPV' : "DTF_FUN(BPVVDCHI2, True,strings(['D0']))",
    'D0_CM_DTF_CHI2' : "DTF_CHI2NDOF( True, strings(['D0']))",

    'D0_D_CM_MM' : "DTF_FUN(CHILD(1,MM), True,strings(['D0']))",
    'D0_D_CM_M' : "DTF_FUN(CHILD(1,M), True,strings(['D0']))",
    'D0_D_CM_PX' : "DTF_FUN(CHILD(1,PX), True,strings(['D0']))",
    'D0_D_CM_PY' : "DTF_FUN(CHILD(1,PY), True,strings(['D0']))",
    'D0_D_CM_PZ' : "DTF_FUN(CHILD(1,PZ), True,strings(['D0']))",
    'D0_D_CM_E' : "DTF_FUN(CHILD(1,E), True,strings(['D0']))",
    'D0_D_CM_MINIPCHI2' : "DTF_FUN(CHILD(1,MIPCHI2DV()), True,strings(['D0']))",
    'D0_D_CM_DIRA_OWNPV' : "DTF_FUN(CHILD(1,BPVDIRA), True,strings(['D0']))",
    'D0_D_CM_ENDVERTEX_CHI2' : "DTF_FUN(CHILD(1,VFASPF(VCHI2)), True,strings(['D0']))",
    'D0_D_CM_ENDVERTEX_NDOF' : "DTF_FUN(CHILD(1,VFASPF(VDOF)), True,strings(['D0']))",
    'D0_D_CM_FDCHI2_OWNPV' : "DTF_FUN(CHILD(1,BPVVDCHI2), True,strings(['D0']))",

    'D0_D0p_CM_ID':   "DTF_FUN(CHILD(1,CHILD(1,ID)),True,strings(['D0']))",
    'D0_D0p_CM_PX':   "DTF_FUN(CHILD(1,CHILD(1,PX)),True,strings(['D0']))",
    'D0_D0p_CM_PY':   "DTF_FUN(CHILD(1,CHILD(1,PY)),True,strings(['D0']))",
    'D0_D0p_CM_PZ':   "DTF_FUN(CHILD(1,CHILD(1,PZ)),True,strings(['D0']))",
    'D0_D0p_CM_E':   "DTF_FUN(CHILD(1,CHILD(1,E)),True,strings(['D0']))",
    'D0_D0p_CM_MINIPCHI2' : "DTF_FUN(CHILD(1,CHILD(1,MIPCHI2DV())), True,strings(['D0']))",

    'D0_D0m_CM_ID':   "DTF_FUN(CHILD(1,CHILD(2,ID)),True,strings(['D0']))",
    'D0_D0m_CM_PX':   "DTF_FUN(CHILD(1,CHILD(2,PX)),True,strings(['D0']))",
    'D0_D0m_CM_PY':   "DTF_FUN(CHILD(1,CHILD(2,PY)),True,strings(['D0']))",
    'D0_D0m_CM_PZ':   "DTF_FUN(CHILD(1,CHILD(2,PZ)),True,strings(['D0']))",
    'D0_D0m_CM_E':   "DTF_FUN(CHILD(1,CHILD(2,E)),True,strings(['D0']))",
    'D0_D0m_CM_MINIPCHI2' : "DTF_FUN(CHILD(1,CHILD(2,MIPCHI2DV())), True,strings(['D0']))",

    'D0_p_CM_ID':   "DTF_FUN(CHILD(2,CHILD(1,ID)),True,strings(['D0']))",
    'D0_p_CM_PX':   "DTF_FUN(CHILD(2,CHILD(1,PX)),True,strings(['D0']))",
    'D0_p_CM_PY':   "DTF_FUN(CHILD(2,CHILD(1,PY)),True,strings(['D0']))",
    'D0_p_CM_PZ':   "DTF_FUN(CHILD(2,CHILD(1,PZ)),True,strings(['D0']))",    
    'D0_p_CM_E':   "DTF_FUN(CHILD(2,CHILD(1,E)),True,strings(['D0']))",    
    'D0_p_CM_MINIPCHI2' : "DTF_FUN(CHILD(2,CHILD(1,MIPCHI2DV())), True,strings(['D0']))",

    'D0_m_CM_ID':   "DTF_FUN(CHILD(2,CHILD(2,ID)),True,strings(['D0']))",
    'D0_m_CM_PX':   "DTF_FUN(CHILD(2,CHILD(2,PX)),True,strings(['D0']))",
    'D0_m_CM_PY':   "DTF_FUN(CHILD(2,CHILD(2,PY)),True,strings(['D0']))",
    'D0_m_CM_PZ':   "DTF_FUN(CHILD(2,CHILD(2,PZ)),True,strings(['D0']))",
    'D0_m_CM_E':   "DTF_FUN(CHILD(2,CHILD(2,E)),True,strings(['D0']))",
    'D0_m_CM_MINIPCHI2' : "DTF_FUN(CHILD(2,CHILD(2,MIPCHI2DV())), True,strings(['D0']))",

    }

tuple.addTool(TupleToolDecay, name="B")
tuple.B.addTupleTool(LoKi_DOCA)
tuple.B.addTupleTool(LoKi_Iso)
tuple.B.addTupleTool(LoKi_LT)
tuple.B.addTupleTool(LoKi_constrain)

# Fit with B0 mass constraint
B0Const = tuple.B.addTupleTool('TupleToolDecayTreeFitter/ConstB0Fit')
B0Const.UpdateDaughters = True
B0Const.daughtersToConstrain += [ 'D0','B0' ]
# Fit with B_s0 mass constraint
BsConst = tuple.B.addTupleTool('TupleToolDecayTreeFitter/ConstBsFit')
BsConst.UpdateDaughters = True
BsConst.Substitutions = { 'B0': 'B_s0' }
BsConst.daughtersToConstrain += [ 'D0','B_s0' ]

#configure differently depending on which data you use
#DaVinci().DataType="2012"
#DaVinci().DDDBtag="dddb-20120831"
#DaVinci().CondDBtag="cond-20121025"
######################################################
DaVinci().DataType="2012"
#DaVinci().DDDBtag="dddb-20130111"
#DaVinci().CondDBtag="cond-20130114"
######################################################

DaVinci().EvtMax = -1
DaVinci().PrintFreq = 100
DaVinci().Simulation = False #changed to false for data
#DaVinci().RedoMCLinks = True
DaVinci().TupleFile="B_Data_nTuples.root"
#DaVinci().UserAlgorithms = [ NodeKiller, scaler, BMassseq, tuple ]
DaVinci().appendToMainSequence( [ NodeKiller, scaler, BMassseq, tuple] )
DaVinci().InputType='MDST'


