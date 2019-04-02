import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <2016|2017|2018> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2015, 2016
note = str(sys.argv[3])

script=str('~/git/lhcb/jets/davinci/Run2Jets_jpsijet.py')
if year=='2017':
    script=str('~/git/lhcb/jets/davinci/Run2Jets_jpsijet17.py')
if year=='2018':
    script=str('~/git/lhcb/jets/davinci/Run2Jets_jpsijet18.py')

job_name = "RIIJ_jpsijet" +str(note) + str(polarity) + str(year)
print job_name
print script

#DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV = GaudiExec(directory="~/DaVinciDev_v44r9")
DV.options = [script]
DV.useGaudiRun = False
BK_locations = []

if year=='2016':
    BK_locations = ['/LHCb/Collision16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16/Stripping28r1/90000000/DIMUON.DST']
elif year=='2017':
    BK_locations = ['/LHCb/Collision17/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco17/Stripping29r2/90000000/DIMUON.DST']
elif year=='2018':
    BK_locations = ['/LHCb/Collision18/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco18/Stripping34/90000000/DIMUON.DST']
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
  inputfiles = [LocalFile('commonSelections.py'), LocalFile('Ntuple.py')]
  )

if note=="test":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4)
    j.backend        = Local()
    j.outputfiles     = [LocalFile("*.root")]
elif note=="testgrid":
    j.splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4)
    j.backend        = Dirac()
    j.outputfiles     = [LocalFile("*.root")]
else:
    j.splitter       = SplitByFiles(filesPerJob = 25, maxFiles = -1, ignoremissing=True)
    j.backend        = Dirac()
    j.outputfiles    = [DiracFile("*.root")]

if note in ["test","testgrid"]:
    j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/DIMUON.DST/00069603/0000/00069603_00000533_1.dimuon.dst'])

queues.add(j.submit)

