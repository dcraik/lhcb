import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <2016> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2015, 2016
note = str(sys.argv[3])

script=str('Run2Jets_dijet.py')

job_name = "RIIJ_dijet" +str(note) + str(polarity) + str(year)
print job_name
print script

DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV.options = [script]
DV.useGaudiRun = False
#BK_locations = []

if year=='2016':
    path = "datasets/LHCb_Collision16_90000000_Beam6500GeVVeloClosed"+polarity+"_Real Data_Reco16_Stripping28_BHADRONCOMPLETEEVENT.DST.py"

j = Job(
  name           = job_name,
  application    = DV,
  #inputdata      = data,
  #do_auto_resubmit = True,
  inputfiles = [LocalFile('commonSelections.py'), LocalFile('Ntuple.py')]
  )

if note=="test":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4)
    j.backend        = Local()
    j.outputfiles     = [LocalFile("*.root")]
elif note=="testgrid":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4)
    j.backend        = Dirac()
    j.outputfiles     = [DiracFile("*.root")]
else:
    j.splitter       = SplitByFiles(filesPerJob = 25, maxFiles = -1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [DiracFile("*.root")]

if note in ["test","testgrid"]:
    j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'])
else:
    j.application.readInputData(path)
j.parallel_submit = True

queues.add(j.submit)

