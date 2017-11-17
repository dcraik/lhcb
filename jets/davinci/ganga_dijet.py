import sys

polarity = str(sys.argv[1]) #MagDown, MagUp
year = str(sys.argv[2]) #2015, 2016

script=str('Run2Jets_dijet.py')

job_name = "RIIJ_dijet" + str(polarity) + str(year)
print job_name
print script

DV = GaudiExec(directory="~/DaVinciDev_v42r6p1")
DV.options = [script]
DV.useGaudiRun = False
#BK_locations = []

if year=='2016':
    path = "LHCb_Collision16_90000000_Beam6500GeVVeloClosed"+polarity+"_Real Data_Reco16_Stripping28_BHADRONCOMPLETEEVENT.DST.py"
#    BK_locations = ['/LHCb/Collision16/Beam6500GeV-VeloClosed-'+polarity+'/Real Data/Reco16/Stripping28/90000000/BHADRONCOMPLETEEVENT.DST']

#data = LHCbDataset()
#bk = BKQuery()

#for path in BK_locations:
#    bk.path = path
#    tmp = bk.getDataset()
#    print path, len(tmp.files)
#    if len(tmp.files) > 0:
#        data.extend( tmp )

#data = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'])
#data = LHCbDataset(['PFN:/eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRONCOMPLETEEVENT.DST/00059907/0001/00059907_00010184_1.bhadroncompleteevent.dst'])

#import sys
#if len(data.files) < 1:
#    sys.exit()

j = Job(
  name           = job_name,
  application    = DV,
  #inputdata      = data,
  #do_auto_resubmit = True,
  inputfiles = [LocalFile('commonSelections.py'), LocalFile('Ntuple.py')]

  splitter       = SplitByFiles(filesPerJob = 25),#, ignoremissing=True),
  backend        = Dirac(),
  outputfiles    = [DiracFile("*.root")],
  
#  splitter       = SplitByFiles(filesPerJob = 1, maxFiles = 4),
#  backend        = Local(),
#  outputfiles     = [LocalFile("*.root")],
  )
j.application.readInputData(path)
j.parallel_submit = True

queues.add(j.submit)

