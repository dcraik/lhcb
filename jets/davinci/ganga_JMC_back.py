import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <490000XX> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
mode = str(sys.argv[2]) #490000XX
note = str(sys.argv[3])

script=str('Run2Jets_JMC_back.py')

job_name = "RIIJ_JMC_back" + str(note) + str(polarity) + str(mode)
print job_name
print script

DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV.options = [script]
DV.useGaudiRun = False

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
  )

if note=="test":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 3, ignoremissing=True)
    j.backend        = Local()
    j.outputfiles     = [LocalFile("*.root")]
elif note=="testgrid":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [LocalFile("*.root")]
else:
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = -1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [DiracFile("*.root")]
    

#j.application.readInputData(path)
j.parallel_submit = True

queues.add(j.submit)

