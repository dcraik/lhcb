import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <42112022|42112052|42112053> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
mode = str(sys.argv[2]) #4211205{2,3}
note = str(sys.argv[3])

script=str('Run2Jets_ZMC.py')

job_name = "RIIJ_ZMC" + str(note) + str(polarity) + str(mode)
print job_name
print script

#DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV = GaudiExec(directory="/workspace/DaVinciDev_v44r9", platform="x86_64-slc6-gcc62-opt")
DV.options = [script]
DV.useGaudiRun = False

BK_locations = ['/MC/2015/Beam6500GeV-2015-'+polarity+'-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1NoPrescalingFlagged/'+mode+'/ALLSTREAMS.DST']
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
  )

if note=="test":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 3)
    j.backend        = Local()
    j.outputfiles     = [LocalFile("*.root")]
elif note=="testgrid":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 3)
    j.backend        = Dirac()
    j.outputfiles    = [LocalFile("*.root")]
else:
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = -1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [DiracFile("*.root")]


#j.application.readInputData(path)
j.parallel_submit = True

queues.add(j.submit)

