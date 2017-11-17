import sys

polarity = str(sys.argv[1]) #MagDown, MagUp
mode = str(sys.argv[2]) #490000XX

script=str('Run2Jets_JMC.py')

job_name = "RIIJ_JMC" + str(polarity) + str(mode)
print job_name
print script

DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV.options = [script]
DV.useGaudiRun = False

#path = '/afs/cern.ch/user/d/dcraik/jets/MC_Dev_'+mode+'_Beam6500GeVRunII'+polarity+'Nu1.625nsPythia8_Sim08f_Reco15DEV_LDST.py'

#DV = GaudiPython(project="DaVinci", version="v41r4p2", script=script)
#DV = GaudiPython(project="DaVinci", version="v41r2p1", script=script)
#files = []
#
#if mode == "49000002":
#    if polarity == "MagDown":
#        files = ['LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000001_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000002_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000003_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000004_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000005_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000006_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000007_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000008_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042952/0000/00042952_00000009_1.ldst']
#    elif polarity == "MagUp":
#        files = []
#elif mode == "49000042":
#    if polarity == "MagDown":
#        files = ['LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000001_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000002_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000003_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000004_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000005_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000006_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000007_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000008_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042978/0000/00042978_00000009_1.ldst']
#    elif polarity == "MagUp":
#        files = []
#elif mode == "49000052":
#    if polarity == "MagDown":
#        files = ['LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000001_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000002_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000003_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000004_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000005_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000006_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000007_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000008_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000009_1.ldst',
#        'LFN:/lhcb/MC/Dev/LDST/00042968/0000/00042968_00000010_1.ldst']
#    elif polarity == "MagUp":
#        files = []
#
#import sys
#if len(files) < 1:
#    sys.exit()

BK_locations = ['/MC/Dev/Beam6500GeV-RunII-'+polarity+'-Nu1.6-25ns-Pythia8/Sim08f/Reco15DEV/'+mode+'/LDST']
data = LHCbDataset()
bk = BKQuery()

for path in BK_locations:
    bk.path = path
    tmp = bk.getDataset()
    print path, len(tmp.files)
    if len(tmp.files) > 0:
        data.extend( tmp )

import sys
if len(data.files) < 1:
    sys.exit()

j = Job(
  name           = job_name,
  application    = DV,
  inputdata      = data,
  #do_auto_resubmit = True,
  inputfiles = [LocalFile('commonSelections.py'), LocalFile('Ntuple.py')]

  splitter       = SplitByFiles(filesPerJob = 1, maxFiles = -1, ignoremissing=True),
  backend        = Dirac(),
  outputfiles    = [DiracFile("*.root")],
  
#  splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 3),
#  backend        = Local(),
#  outputfiles     = [LocalFile("*.root")],
  )
#j.application.readInputData(path)
j.parallel_submit = True

queues.add(j.submit)

