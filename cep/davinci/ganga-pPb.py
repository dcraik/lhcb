import sys

if len(sys.argv)<4:
    print "Usage: ", sys.argv[0], "<MagUp|MagDown> <2016-pPb|2016-Pbp> <note>"
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2015, 2016
note = str(sys.argv[3])

script=str('davinci-pPb.py')

job_name = "CEP phi" +str(note) + str(polarity) + str(year)
print job_name
print script

DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV.options = [script]
DV.useGaudiRun = False

if year not in ['2016-pPb','2016-Pbp']:##TODO only set up for 2016 data so far
    sys.exit()

BK_locations = []
if year=='2016-Pbp':
    BK_locations += ['/LHCb/Ionproton16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16pLead/Stripping30r2/90000000/IFT.DST']
elif year=='2016-pPb':
    BK_locations +=['/LHCb/Protonion16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16pLead/Stripping30r3/90000000/IFT.DST']
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
        j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst'])
if note in ["test"]:
        j.inputdata = LHCbDataset(['PFN:/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst'])

j.parallel_submit = True

queues.add(j.submit)

