import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <mode> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
mode = str(sys.argv[2]) #2015, 2016
note = str(sys.argv[3])

script=str('davinci-MC.py')

job_name = "CEP phi MC" +str(note) + str(polarity) + str(mode)
print job_name
print script

DV = GaudiExec(directory="/workspace/DaVinciDev_v44r10p2")
DV.options = [script]
DV.useGaudiRun = False

BK_locations = []
if mode=='30000000':
    BK_locations = ['/MC/2015/Beam6500GeV-2015-'+polarity+'-Nu1.6-25ns-Pythia8/Sim09b/Trig0x411400a2/Reco15a/Turbo02/Stripping24NoPrescalingFlagged/30000000/ALLSTREAMS.DST']
else:
    sys.exit()

data = LHCbDataset()
bk = BKQuery()

if note not in ["test","testgrid"]:
    for path in BK_locations:
        bk.path = path
        tmp = bk.getDataset()
        print path, len(tmp.files)
        if len(tmp.files) > 0:
            data.extend( tmp )

    if len(data.files) < 1:
        sys.exit()

j = Job(
  name           = job_name,
  application    = DV,
  inputdata      = data,
  #do_auto_resubmit = True,
  inputfiles = [LocalFile('Ntuple.py')]
  )

if note=="test":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4)
    j.backend        = Local()
    j.outputfiles     = [LocalFile("*.root")]
elif note=="testgrid":
    j.splitter       = SplitByFiles(filesPerJob = 25, maxFiles = 10, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles     = [LocalFile("*.root")]
else:
    j.splitter       = SplitByFiles(filesPerJob = 25, maxFiles = -1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [DiracFile("*.root")]

if note in ["testgrid"]:
        j.inputdata = LHCbDataset(['LFN:/lhcb/MC/2015/ALLSTREAMS.DST/00057780/0000/00057780_00000005_3.AllStreams.dst'])
if note in ["test"]:
        j.inputdata = LHCbDataset(['PFN:/data/cep-phi/00057780_00000005_3.AllStreams.dst'])

j.parallel_submit = True

queues.add(j.submit)

