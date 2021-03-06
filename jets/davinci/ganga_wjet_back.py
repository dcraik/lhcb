import sys

if len(sys.argv)<4:
    print("Usage: ", sys.argv[0], "<MagUp|MagDown> <2016|2017|2018> <note>")
    sys.exit()
polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2016, 2017, 2018
note = str(sys.argv[3])

script=str('Run2Jets_wjet_back.py')
if year=='2015':
    script=str('Run2Jets_wjet15_back.py')
if year=='2017':
    script=str('Run2Jets_wjet17_back.py')
if year=='2018':
    script=str('Run2Jets_wjet18_back.py')

job_name = "RIIJ_Wjet_back" +str(note) + str(polarity) + str(year)
print(job_name)
print(script)

#DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV = GaudiExec(directory="/workspace/DaVinciDev_v44r9", platform="x86_64-slc6-gcc62-opt")
DV.options = [script]
DV.useGaudiRun = False

BK_locations = []

if year=='2015':
    BK_locations = ['/LHCb/Collision15/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco15a/Stripping24r1/90000000/EW.DST']
elif year=='2016':
    BK_locations = ['/LHCb/Collision16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16/Stripping28r1/90000000/EW.DST']
elif year=='2017':
    BK_locations = ['/LHCb/Collision17/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco17/Stripping29r2/90000000/EW.DST']
elif year=='2018':
    BK_locations = ['/LHCb/Collision18/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco18/Stripping34/90000000/EW.DST']
else:
    sys.exit()
data = LHCbDataset()
bk = BKQuery()

if note not in ["test","testgrid"]:
    for path in BK_locations:
        bk.path = path
        tmp = bk.getDataset()
        print(path, len(tmp.files))
        if len(tmp.files) > 0:
            data.extend( tmp )

    if len(data.files) < 1:
        sys.exit()

#if year=='2016':
#    path = "datasets/LHCb_Collision16_Beam6500GeVVeloClosed"+polarity+"_Real Data_Reco16_Stripping28r1_90000000_EW.DST.py"

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

if note in ["testgrid"]:
    j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/EW.DST/00070440/0000/00070440_00000660_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000264_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000058_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000481_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000308_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000474_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000564_1.ew.dst',
                               'LFN:/lhcb/LHCb/Collision16/EW.DST/00069527/0000/00069527_00000270_1.ew.dst'
                               ])
if note in ["test"]:
    j.inputdata = LHCbDataset(['PFN:/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00005740_1.ew.dst'])
#else:
#    j.application.readInputData(path)

j.parallel_submit = True

queues.add(j.submit)

