#StdAllNoPIDs particles - inputs for D's and Z's
from StandardParticles import StdAllNoPIDsKaons    as loosekaons
from StandardParticles import StdAllNoPIDsPions    as loosepions
from StandardParticles import StdAllNoPIDsProtons  as looseprotons
from StandardParticles import StdAllLooseMuons as looseZmuons

#inputs for HLT-like SVs and muons
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

from PhysSelPython.Wrappers import SimpleSelection, MergedSelection, DataOnDemand, Selection
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles

###########################
###### HLT-like SVs #######

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

###########################
###### HLT-like mus #######

recMus = SimpleSelection (
    'recMus'         ,
    FilterDesktop   ,
    [ loosemuons ]    ,
    DecayDescriptor = "[mu+]cc",
    Code = ("(PROBNNmu>0.5) "
            "& (PT>1000*MeV)"
            "& (TRGHOSTPROB<0.2)")
    )

###########################
###### Z candidates #######

ZMuons = SimpleSelection (
    'ZMuons'         ,
    FilterDesktop   ,
    [ looseZmuons ]    ,
    DecayDescriptor = "[mu+]cc",
    Code = (
            "(PT>3*GeV)"
           )
    )

dcZ = { }
for child in ['mu+'] :
    dcZ[child] = "(PT > 3000*MeV)"

combcutsZ = "in_range(40000*MeV,  AM, 200000*MeV)" \

parentcutsZ = "MM > 40000*MeV"

looseZs = SimpleSelection (
    'looseZ',
    CombineParticles,
    [ZMuons],
    DecayDescriptor = "Z0 -> mu+ mu-",
    DaughtersCuts   = dcZ,
    CombinationCut = (combcutsZ),
    MotherCut      =  (parentcutsZ),
)


#looseZs = DataOnDemand('/Event/EW/Phys/Z02MuMuLine/Particles')

Zs = SimpleSelection (
    'Zs'           ,
    FilterDesktop   ,
    [ looseZs ]    ,
    DecayDescriptor = 'Z0 -> mu+ mu-',
    Code = "(MINTREE('mu+'==ABSID,PT) > 10.*GeV)"
    )


###########################
###### D candidates #######

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
                  "& (MIPCHI2DV(PRIMARY) > 4)" #16)"
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

parentcutsDp = "(VFASPF(VCHI2PDOF) < 6)" \
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

parentcutsDs = "(VFASPF(VCHI2PDOF) < 6)" \
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

dcLc = { }
for child in ['pi+','K+','p+'] :
    dcLc[child] = "(PT > 200*MeV)" \
                  "& (P > 3*GeV)" \
                  "& (MIPCHI2DV(PRIMARY) > 4)"
#                  "& (TRCHI2 < 3)" \
#                  "& (P > 2*GeV)" \

combcutsLc = "in_range(2206*MeV,  AM, 2366*MeV)" \
             "& (ANUM(PT > 400*MeV) > 1)" \
             "& (ANUM(PT > 1000*MeV) > 0)" \
             "& (AMINDOCA('') < 0.5*mm )" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 4) > 1)" \
             "& (ANUM(MIPCHI2DV(PRIMARY) > 6) > 0)"
#             "& (AMINDOCA('') < 0.1*mm )" \

parentcutsLc = "(VFASPF(VCHI2PDOF) < 25)" \
               "& BPVVALID()" \
               "& ((BPVVDCHI2 > 4.0) | (BPVLTIME() > 0.075*ps))" \
               "& (BPVDIRA > 0.9994 )"
               #"& (BPVVDCHI2 > 16 )" \

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
             "& (AMINDOCA('') < 0.1*mm )" \
             "& (ANUM(PT > 400*MeV) > 1)" \
             "& (ANUM(PT > 1000*MeV) > 0)"

parentcutsD02K3pi = "(VFASPF(VCHI2PDOF) < 10)" \
               "& BPVVALID()" \
               "& (BPVVDCHI2 > 16 )" \
               "& (BPVDIRA > 0.99985 )"
               #"& (BPVLTIME() > 0.10*ps )" \

recD02K3pi = SimpleSelection (
    'recD02K3pi',
    CombineParticles,
    [charmKaons,charmPions],
    DecayDescriptor = "[D0 -> K- pi+ pi+ pi-]cc",
    DaughtersCuts   = dcD02K3pi,
    CombinationCut = (combcutsD02K3pi),
    MotherCut      =  (parentcutsD02K3pi),
)

