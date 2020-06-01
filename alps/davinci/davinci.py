#!/bin/env python
import sys
sys.path.append('.')

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = """
   HLT_PASS_RE ( '.*DarkBoson.*' )
   """
    )

stream='Bhadron'
lines = ['B2KX2KKPi','B2KpiX2KKPi','B2KX2KKPiM','B2KpiX2KKPiM','B2KX2PiPiPi','B2KpiX2PiPiPi','B2KX2PiPiPiM','B2KpiX2PiPiPiM','B2KX2EtaPiPi23PI','B2KpiX2EtaPiPi23PI','B2KX2EtaPiPi2GG','B2KpiX2EtaPiPi2GG',
         'B2KX24Pi','B2KpiX24Pi','B2KX26Pi','B2KpiX26Pi','B2KX22K2Pi','B2KpiX22K2Pi','B2KX24K','B2KpiX24K','B2KpiX2MuMuDD','B2KpiX2MuMuDDSS','B2KX2MuMuDD','B2KX2MuMuDDSS','B2KpiX2MuMu','B2KpiX2MuMuSS',
         'B2KX2MuMu','B2KX2MuMuSS','B2KpiX2PiPi','B2KX2PiPi','B2KX2PiPiSS','B2KpiX2KK','B2KX2KK','B2KX2KKSS','B2KpiX2EE','B2KpiX2EESS','B2KX2EE','B2KX2EESS','B2JK','B2JKst','B2KstX2GammaGamma',
         'B2KstX2PiGamma','B2KstX2PiGammaM']


## Data.
#from GaudiConf import IOHelper
#IOHelper('ROOT').inputFiles(['/data/alps/00103400_00000004_1.bhadron.mdst'],#'/tmp/dcraik/00095518_00000372_1.bhadron.mdst'],
#                            clear = True)
##Type = 'MC'


#from StandardParticles import StdAllLooseMuons as loosemuons
#from StandardParticles import StdAllNoPIDsKaons as loosekaons
#from StandardParticles import StdAllNoPIDsPions as loosepions
#from StandardParticles import StdLooseAllPhotons as loosephotons
#from PhysSelPython.Wrappers import SimpleSelection, MergedSelection, DataOnDemand, Selection
#from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles
#
###basic selection
#muons = SimpleSelection (
#    'muons'         ,
#    FilterDesktop   ,
#    [ loosemuons ]    ,
#    DecayDescriptor = "[mu+]cc",
#    Code = ( "(PIDmu - PIDpi > -5) & "
#            "(PT>100*MeV)"
#            "& (TRGHOSTPROB<0.3)"
#            "& (MIPCHI2DV(PRIMARY) > 9)")
#    )
#
#pions = SimpleSelection (
#    'pions'         ,
#    FilterDesktop   ,
#    [ loosepions ]    ,
#    DecayDescriptor = "[pi+]cc",
#    Code = ( "(PROBNNpi > 0.2) & "
#            "(PT>250*MeV)"
#            "& (P>3000*MeV)"
#            "& (TRGHOSTPROB<0.3)"
#            "& (MIPCHI2DV(PRIMARY) > 36)")
#    )
#
#kaons = SimpleSelection (
#    'kaons'         ,
#    FilterDesktop   ,
#    [ loosekaons ]    ,
#    DecayDescriptor = "[K+]cc",
#    Code = ("(PROBNNK > 0.1) & "
#            "(PT>250*MeV)"
#            "& (P>3000*MeV)"
#            "& (TRGHOSTPROB<0.3)"
#            "& (MIPCHI2DV(PRIMARY) > 25)")
#    )
#
#photons = SimpleSelection (
#    'photons'         ,
#    FilterDesktop   ,
#    [ loosephotons ]    ,
#    DecayDescriptor = "gamma",
#    Code = (
#            "(PT>500*MeV)"
#            "& (P>1000*MeV)"
#            "& (CL>0.3)")
#    )
#
#loosepizM = DataOnDemand(Location="Phys/StdLooseMergedPi0/Particles")
#loosepizR = DataOnDemand(Location="Phys/StdLooseResolvedPi0/Particles")
#looseeta  = DataOnDemand(Location="Phys/StdLooseResolvedEta/Particles")
#
#pizM = SimpleSelection (
#    'pizM'         ,
#    FilterDesktop   ,
#    [ loosepizM ]    ,
#    DecayDescriptor = "pi0",
#    Code = (
#            "(PT>500*MeV)"
#            "& (P>3000*MeV)"
#            "& (CL>0.1)")
#    )
#
#pizR = SimpleSelection (
#    'pizR'         ,
#    FilterDesktop   ,
#    [ loosepizR ]    ,
#    DecayDescriptor = "pi0",
#    Code = (
#            "(PT>500*MeV)"
#            "& (CL>0.1)")
#    )
#
#eta = SimpleSelection (
#    'eta'         ,
#    FilterDesktop   ,
#    [ looseeta ]    ,
#    DecayDescriptor = "eta",
#    Code = (
#            "(PT>500*MeV)"
#            "& (P>2000*MeV)"
#            "& (MM>450*MeV)"
#            "& (MM<650*MeV)"
#            "& (CL>0.2)")
#    )
#
###X selections
#combcutsLL = "(APT > 250.0*MeV)"
#parentcutsLL = "(VFASPF(VCHI2PDOF) < 10)"
#combcutsPiPi = "(APT > 250.0*MeV)"
#parentcutsPiPi = "(VFASPF(VCHI2PDOF) < 5)"
#
#dimuons = SimpleSelection (
#    'dimuons',
#    CombineParticles,
#    [muons],
#    DecayDescriptor = "KS0 -> mu+ mu-",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#dimuonsSS = SimpleSelection (
#    'dimuonsSS',
#    CombineParticles,
#    [muons],
#    DecayDescriptor = "[KS0 -> mu+ mu+]cc",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#dikaons = SimpleSelection (
#    'dikaons',
#    CombineParticles,
#    [kaons],
#    DecayDescriptor = "KS0 -> K+ K-",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#quadkaons = SimpleSelection (
#    'quadkaons',
#    CombineParticles,
#    [kaons],
#    DecayDescriptor = "KS0 -> K+ K- K+ K-",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#dikaonsSS = SimpleSelection (
#    'dikaonsSS',
#    CombineParticles,
#    [kaons],
#    DecayDescriptor = "[KS0 -> K+ K+]cc",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#dipions = SimpleSelection (
#    'dipions',
#    CombineParticles,
#    [pions],
#    DecayDescriptor = "KS0 -> pi+ pi-",
#    CombinationCut = (combcutsPiPi),
#    MotherCut      =  (parentcutsPiPi),
#)
#
#quadpions = SimpleSelection (
#    'quadpions',
#    CombineParticles,
#    [pions],
#    DecayDescriptor = "KS0 -> pi+ pi- pi+ pi-",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#diKdipi = SimpleSelection (
#    'diKdipi',
#    CombineParticles,
#    [kaons,pions],
#    DecayDescriptor = "KS0 -> K+ K- pi+ pi-",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#dipionsSS = SimpleSelection (
#    'dipionsSS',
#    CombineParticles,
#    [pions],
#    DecayDescriptor = "[KS0 -> pi+ pi+]cc",
#    CombinationCut = (combcutsPiPi),
#    MotherCut      =  (parentcutsPiPi),
#)
#
#parentcutsDiG = "(M < 5000*MeV) & (PT > 2000*MeV)"
#
#diphotons = SimpleSelection (
#    'diphotons',
#    CombineParticles,
#    [photons],
#    DecayDescriptor   = "KS0 -> gamma gamma",
#    MotherCut         = parentcutsDiG,
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#pizgammaM = SimpleSelection (
#    'pizgammaM',
#    CombineParticles,
#    [pizM,photons],
#    DecayDescriptor   = "eta -> pi0 gamma",
#    MotherCut         = parentcutsDiG,
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#pizgammaR = SimpleSelection (
#    'pizgammaR',
#    CombineParticles,
#    [pizR,photons],
#    DecayDescriptor   = "eta -> pi0 gamma",
#    MotherCut         = parentcutsDiG,
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#combcuts3H = "(APT > 500.0*MeV)"
#parentcuts3H = "(VFASPF(VCHI2PDOF) < 10)"
#
#threepiM = SimpleSelection (
#    'threepiM',
#    CombineParticles,
#    [pions,pizM],
#    DecayDescriptor = "KS0 -> pi+ pi- pi0",
#    CombinationCut = (combcuts3H),
#    MotherCut      =  (parentcuts3H),
#)
#
#threepiR = SimpleSelection (
#    'threepiR',
#    CombineParticles,
#    [pions,pizR],
#    DecayDescriptor = "KS0 -> pi+ pi- pi0",
#    CombinationCut = (combcuts3H),
#    MotherCut      =  (parentcuts3H),
#)
#
#twoKpiM = SimpleSelection (
#    'twoKpiM',
#    CombineParticles,
#    [kaons,pizM],
#    DecayDescriptor = "KS0 -> K+ K- pi0",
#    CombinationCut = (combcuts3H),
#    MotherCut      =  (parentcuts3H),
#)
#
#twoKpiR = SimpleSelection (
#    'twoKpiR',
#    CombineParticles,
#    [kaons,pizR],
#    DecayDescriptor = "KS0 -> K+ K- pi0",
#    CombinationCut = (combcuts3H),
#    MotherCut      =  (parentcuts3H),
#)
#
#hexpions = SimpleSelection (
#    'hexpions',
#    CombineParticles,
#    [pions,pizR],
#    DecayDescriptor = "KS0 -> pi+ pi- pi0 pi+ pi- pi0",
#    CombinationCut = (combcutsLL),
#    MotherCut      =  (parentcutsLL),
#)
#
#combcutsetapipi = "(APT > 2000.0*MeV)"
#parentcutsetapipi = "(VFASPF(VCHI2PDOF) < 10)"
#
#etapipi = SimpleSelection (
#    'etapipi',
#    CombineParticles,
#    [pions,eta],
#    DecayDescriptor = "KS0 -> eta pi+ pi-",
#    CombinationCut = (combcutsetapipi),
#    MotherCut      =  (parentcutsetapipi),
#)
#
###B->K(*)X selections
#pionsB = SimpleSelection (
#    'pionsB'         ,
#    FilterDesktop   ,
#    [ loosepions ]    ,
#    DecayDescriptor = "[pi+]cc",
#    Code = ( #"(PIDK - PIDpi < 3) & "
#            "(PT>250*MeV)"
#            "& (P>2000*MeV)"
#            "& (TRGHOSTPROB<0.3)"
#            "& (MIPCHI2DV(PRIMARY) > 9)")
#    )
#
#kaonsB = SimpleSelection (
#    'kaonsB'         ,
#    FilterDesktop   ,
#    [ loosekaons ]    ,
#    DecayDescriptor = "[K+]cc",
#    Code = (#"(PIDK - PIDpi > 5) & "
#            "(PT>250*MeV)"
#            "& (P>2000*MeV)"
#            "& (TRGHOSTPROB<0.3)"
#            "& (MIPCHI2DV(PRIMARY) > 9)")
#    )
#
#dcKst = { }
#dcKst['K+'] = "(PT > 250.0) & (P > 2000.0) & (TRGHOSTPROB < 0.3) & (PROBNNK > 0.2)"
#dcKst['pi+'] = "(PT > 250.0) & (P > 2000.0) & (TRGHOSTPROB < 0.3) & (PROBNNpi > 0.2)"
#combcutsKst = "(AMAXDOCA('') < 0.5)"
#parentcutsKst = "(VFASPF(VCHI2PDOF) < 10) & (MIPCHI2DV(PRIMARY)> 16)"
#
#kst = SimpleSelection (
#    'kst',
#    CombineParticles,
#    [kaonsB,pionsB],
#    DecayDescriptor = "[K*(892)0 -> K+ pi-]cc",
#    DaughtersCuts   = dcKst,
#    CombinationCut = (combcutsKst),
#    MotherCut      =  (parentcutsKst),
#)
#
##B selections
#
#combcutsB = "(AM>4800*MeV) & (AM<6000*MeV) & (AMAXDOCA('') < 0.5)"
#parentcutsB = "(VFASPF(VCHI2PDOF) < 10)"
#
#kstmumu = SimpleSelection (
#    'kstmumu',
#    CombineParticles,
#    [kst,dimuons],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kmumu = SimpleSelection (
#    'kmumu',
#    CombineParticles,
#    [kaonsB,dimuons],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kstmumuSS = SimpleSelection (
#    'kstmumuSS',
#    CombineParticles,
#    [kst,dimuonsSS],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kmumuSS = SimpleSelection (
#    'kmumuSS',
#    CombineParticles,
#    [kaonsB,dimuonsSS],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst2k = SimpleSelection (
#    'kst2k',
#    CombineParticles,
#    [kst,dikaons],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2k = SimpleSelection (
#    'k2k',
#    CombineParticles,
#    [kaonsB,dikaons],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2kSS = SimpleSelection (
#    'k2kSS',
#    CombineParticles,
#    [kaonsB,dikaonsSS],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#
#kst2pi = SimpleSelection (
#    'kst2pi',
#    CombineParticles,
#    [kst,dipions],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2pi = SimpleSelection (
#    'k2pi',
#    CombineParticles,
#    [kaonsB,dipions],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2piSS = SimpleSelection (
#    'k2piSS',
#    CombineParticles,
#    [kaonsB,dipionsSS],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#combcutsB2neutrals = "(AM>4800*MeV) & (AM<6000*MeV)"
#parentcutsB2neutrals = "(PT>3000*MeV)"
#
#kstgg = SimpleSelection (
#    'kstgg',
#    CombineParticles,
#    [kst,diphotons],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB2neutrals),
#    MotherCut      =  (parentcutsB2neutrals),
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#kstpizgM = SimpleSelection (
#    'kstpizgM',
#    CombineParticles,
#    [kst,pizgammaM],
#    DecayDescriptor = "[B0 -> K*(892)0 eta]cc",
#    CombinationCut = (combcutsB2neutrals),
#    MotherCut      =  (parentcutsB2neutrals),
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#kstpizgR = SimpleSelection (
#    'kstpizgR',
#    CombineParticles,
#    [kst,pizgammaR],
#    DecayDescriptor = "[B0 -> K*(892)0 eta]cc",
#    CombinationCut = (combcutsB2neutrals),
#    MotherCut      =  (parentcutsB2neutrals),
#    ParticleCombiners = {"" : "MomentumCombiner:PUBLIC"},
#    ReFitPVs          = False,
#)
#
#kst3piM = SimpleSelection (
#    'kst3piM',
#    CombineParticles,
#    [kst,threepiM],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k3piM = SimpleSelection (
#    'k3piM',
#    CombineParticles,
#    [kaonsB,threepiM],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst3piR = SimpleSelection (
#    'kst3piR',
#    CombineParticles,
#    [kst,threepiR],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k3piR = SimpleSelection (
#    'k3piR',
#    CombineParticles,
#    [kaonsB,threepiR],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst2kpiM = SimpleSelection (
#    'kst2kpiM',
#    CombineParticles,
#    [kst,twoKpiM],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2kpiM = SimpleSelection (
#    'k2kpiM',
#    CombineParticles,
#    [kaonsB,twoKpiM],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst2kpiR = SimpleSelection (
#    'kst2kpiR',
#    CombineParticles,
#    [kst,twoKpiR],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2kpiR = SimpleSelection (
#    'k2kpiR',
#    CombineParticles,
#    [kaonsB,twoKpiR],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst4pi = SimpleSelection (
#    'kst4pi',
#    CombineParticles,
#    [kst,quadpions],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k4pi = SimpleSelection (
#    'k4pi',
#    CombineParticles,
#    [kaonsB,quadpions],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst4k = SimpleSelection (
#    'kst4k',
#    CombineParticles,
#    [kst,quadkaons],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k4k = SimpleSelection (
#    'k4k',
#    CombineParticles,
#    [kaonsB,quadkaons],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst2k2pi = SimpleSelection (
#    'kst2k2pi',
#    CombineParticles,
#    [kst,diKdipi],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k2k2pi = SimpleSelection (
#    'k2k2pi',
#    CombineParticles,
#    [kaonsB,diKdipi],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kst6pi = SimpleSelection (
#    'kst6pi',
#    CombineParticles,
#    [kst,hexpions],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#k6pi = SimpleSelection (
#    'k6pi',
#    CombineParticles,
#    [kaonsB,hexpions],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#kstetapipi = SimpleSelection (
#    'kstetapipi',
#    CombineParticles,
#    [kst,etapipi],
#    DecayDescriptor = "[B0 -> K*(892)0 KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#ketapipi = SimpleSelection (
#    'ketapipi',
#    CombineParticles,
#    [kaonsB,etapipi],
#    DecayDescriptor = "[B+ -> K+ KS0]cc",
#    CombinationCut = (combcutsB),
#    MotherCut      =  (parentcutsB),
#)
#
#
#from PhysSelPython.Wrappers import SelectionSequence
#
#kstmumu_seq   = SelectionSequence('kstmumu_seq',   TopSelection=kstmumu)
#kmumu_seq     = SelectionSequence('kmumu_seq',     TopSelection=kmumu)
#kstmumuSS_seq = SelectionSequence('kstmumuSS_seq', TopSelection=kstmumuSS)
#kmumuSS_seq   = SelectionSequence('kmumuSS_seq',   TopSelection=kmumuSS)
#kst2k_seq     = SelectionSequence('kst2k_seq',     TopSelection=kst2k)
#k2k_seq       = SelectionSequence('k2k_seq',       TopSelection=k2k)
#k2kSS_seq     = SelectionSequence('k2kSS_seq',     TopSelection=k2kSS)
#kst2pi_seq    = SelectionSequence('kst2pi_seq',    TopSelection=kst2pi)
#k2pi_seq      = SelectionSequence('k2pi_seq',      TopSelection=k2pi)
#k2piSS_seq    = SelectionSequence('k2piSS_seq',    TopSelection=k2piSS)
#
#kstgg_seq     = SelectionSequence('kstgg_seq',     TopSelection=kstgg)
#kstpizgM_seq  = SelectionSequence('kstpizgM_seq',  TopSelection=kstpizgM)
#kstpizgR_seq  = SelectionSequence('kstpizgR_seq',  TopSelection=kstpizgR)
#
#kst3piM_seq   = SelectionSequence('kst3piM_seq',   TopSelection=kst3piM)
#kst3piR_seq   = SelectionSequence('kst3piR_seq',   TopSelection=kst3piR)
#k3piM_seq     = SelectionSequence('k3piM_seq',     TopSelection=k3piM)
#k3piR_seq     = SelectionSequence('k3piR_seq',     TopSelection=k3piR)
#kst2kpiM_seq  = SelectionSequence('kst2kpiM_seq',  TopSelection=kst2kpiM)
#kst2kpiR_seq  = SelectionSequence('kst2kpiR_seq',  TopSelection=kst2kpiR)
#k2kpiM_seq    = SelectionSequence('k2kpiM_seq',    TopSelection=k2kpiM)
#k2kpiR_seq    = SelectionSequence('k2kpiR_seq',    TopSelection=k2kpiR)
#
#kst4pi_seq    = SelectionSequence('kst4pi_seq',    TopSelection=kst4pi)
#k4pi_seq      = SelectionSequence('k4pi_seq',      TopSelection=k4pi)
#kst4k_seq     = SelectionSequence('kst4k_seq',     TopSelection=kst4k)
#k4k_seq       = SelectionSequence('k4k_seq',       TopSelection=k4k)
#kst2k2pi_seq  = SelectionSequence('kst2k2pi_seq',  TopSelection=kst2k2pi)
#k2k2pi_seq    = SelectionSequence('k2k2pi_seq',    TopSelection=k2k2pi)
#kst6pi_seq    = SelectionSequence('kst6pi_seq',    TopSelection=kst6pi)
#k6pi_seq      = SelectionSequence('k6pi_seq',      TopSelection=k6pi)
#
#kstetapipi_seq= SelectionSequence('kstetapipi_seq',TopSelection=kstetapipi)
#ketapipi_seq  = SelectionSequence('ketapipi_seq',  TopSelection=ketapipi)

##########################

# Turbo/DaVinci configuration.
from Configurables import DstConf, TurboConf, DaVinci
DaVinci().Simulation = False
DaVinci().Lumi = True
DaVinci().TupleFile = "LumiTuple.root"
#DaVinci().appendToMainSequence([kstmumu_seq.sequence(), kstmumuSS_seq.sequence()])
#DaVinci().appendToMainSequence([kmumu_seq.sequence(), kmumuSS_seq.sequence()])
#DaVinci().appendToMainSequence([kst2k_seq.sequence()])
#DaVinci().appendToMainSequence([k2k_seq.sequence(), k2kSS_seq.sequence()])
#DaVinci().appendToMainSequence([kst2pi_seq.sequence()])
#DaVinci().appendToMainSequence([k2pi_seq.sequence(), k2piSS_seq.sequence()])
#DaVinci().appendToMainSequence([kstgg_seq.sequence()])
#DaVinci().appendToMainSequence([kstpizgM_seq.sequence(),kstpizgR_seq.sequence()])
#DaVinci().appendToMainSequence([kst3piM_seq.sequence(),kst3piR_seq.sequence()])
#DaVinci().appendToMainSequence([k3piM_seq.sequence(),k3piR_seq.sequence()])
#DaVinci().appendToMainSequence([kst2kpiM_seq.sequence(),kst2kpiR_seq.sequence()])
#DaVinci().appendToMainSequence([k2kpiM_seq.sequence(),k2kpiR_seq.sequence()])
#DaVinci().appendToMainSequence([k4pi_seq.sequence(), kst4pi_seq.sequence()])
#DaVinci().appendToMainSequence([k4k_seq.sequence(), kst4k_seq.sequence()])
#DaVinci().appendToMainSequence([k2k2pi_seq.sequence(), kst2k2pi_seq.sequence()])
#DaVinci().appendToMainSequence([k6pi_seq.sequence(), kst6pi_seq.sequence()])
#DaVinci().appendToMainSequence([ketapipi_seq.sequence(), kstetapipi_seq.sequence()])
DaVinci().DataType = '2016'
DaVinci().InputType='MDST'
DaVinci().RootInTES='/Event/'+stream
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
    tool.HltDecReportsLocation = '/Event/'+stream+'/' + stage + '/DecReports'
    tool.HltSelReportsLocation = '/Event/'+stream+'/' + stage + '/SelReports'

# Access to classes.
from collections import OrderedDict
import ROOT, array, GaudiPython
from GaudiPython.Bindings import gbl
STD  = gbl.std
LHCB = gbl.LHCb

from Ntuple import Ntuple
#from NtupleMinimal import Ntuple

# GaudiPython configuration.
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

# Run.
import sys, ROOT
from math import floor
#evtmax = 10000 #TODO
evtmax = -1 #TODO
#try: evtmax = int(sys.argv[1])
#except: evtmax = float('inf')
evtnum = 0
ntuple = Ntuple('output.root', tes, gaudi.toolsvc(), gaudi.detSvc(), primaries='/Event/'+stream+'/Rec/Vertex/Primary')
while evtmax < 0 or evtnum < evtmax:
    gaudi.run(1)
    if not bool(tes['/Event']): break
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    try: ntuple.ntuple['evt_pvr_n'][0] = len(tes['/Event/'+stream+'/Rec/Vertex/Primary'])
    except: continue

    fill = False;

    # fill combinations
    for idx, line in enumerate(lines):
        try:
            bs = tes['/Event/'+stream+'/Phys/'+line+'DarkBosonLine/Particles']
            if bs:
                for b in bs:
                    ntuple.addB(b,mode=idx)
                    fill = True
        except:
            #raise
            pass

    # Fill the ntuple.
    if fill: ntuple.fill();
ntuple.close()
