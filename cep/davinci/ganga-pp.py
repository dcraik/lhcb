import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <2015|2016> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2015, 2016
note = str(sys.argv[3])

script=str('DV_pp16.py')
if year=='2015':
    script = str('DV_pp15.py')

job_name = "CEP phi comb " +str(note) + str(polarity) + str(year)
print job_name
print script

DV = GaudiExec(directory="/workspace/DaVinciDev_v44r10p2")
DV.options = [script]

if year not in ['2015','2016']:
    sys.exit()

BK_locations = ['/LHCb/Collision16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16/Stripping28r1/90000000/EW.DST']
if year=='2015':
    BK_locations = ['/LHCb/Collision15/Beam2510GeV-VeloClosed-'+polarity+'/Real Data/Reco15a/Stripping22b/90000000/ALL.DST']
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
    j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst'])
if note in ["test"]:
    j.inputdata = [LocalFile('/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst')]

j.parallel_submit = True

queues.add(j.submit)

