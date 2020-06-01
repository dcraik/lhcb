import sys

if len(sys.argv)<5:
    print("Usage: ", sys.argv[0], "<bhadron/dimuon/leptonic> <year> <MagDown/MagUp> <note>")
    sys.exit()
which = str(sys.argv[1]) #bhadron, dimuon, leptonic
year = str(sys.argv[2]) #2011, 2012, 2015, 2016, 2018
pol = str(sys.argv[3]) #MagDown, MagUp
note = str(sys.argv[4])

if which not in ['bhadron','dimuon','leptonic']:
    print("unknown set ", which)
    sys.exit()

#script=str('dv-'+str(which)+'.py')
script=str('davinci.py')

job_name = "ALP" + str(note) + str(year) + str(pol) + str(which)
print(job_name)
print(script)

if year in ["2011","2012"]:
    DV = GaudiExec(directory="~/DaVinciDev_v39r1p6", platform="x86_64-slc6-gcc62-opt")
elif year in ["2015","2016","2017","2018"]:
    DV = GaudiExec(directory="/workspace/DaVinciDev_v44r10p3", platform="x86_64-slc6-gcc62-opt")
else:
    print("unknown year ", year)
    sys.exit()

DV.options = [script]
DV.useGaudiRun = False

if pol not in ["MagDown","MagUp"]:
    print("unknown polarity ". pol)
    sys.exit()

if year == "2011":
    #basepath="/validation/Collision11/Beam3500GeV-VeloClosed-MagUp/Real Data/Reco14//Stripping21r1p2/90000000/"
    basepath="/LHCb/Collision11/Beam3500GeV-VeloClosed-"+pol+"/Real Data/Reco14/Stripping21r1p2/90000000/"
elif year == "2012":
    #basepath="/validation/Collision12/Beam4000GeV-VeloClosed-MagDown/Real Data/Reco14//Stripping21r0p2/90000000/"
    basepath="/LHCb/Collision12/Beam4000GeV-VeloClosed-"+pol+"/Real Data/Reco14/Stripping21r0p2/90000000/"
elif year == "2015":
    #basepath="/validation/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a//Stripping24r2/90000000/"
    basepath="/LHCb/Collision15/Beam6500GeV-VeloClosed-"+pol+"/Real Data/Reco15a/Stripping24r2/90000000/"
elif year == "2016":
    #basepath="/validation/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16//Stripping28r2/90000000/"
    basepath="/LHCb/Collision16/Beam6500GeV-VeloClosed-"+pol+"/Real Data/Reco16/Stripping28r2/90000000/"
elif year == "2017":
    basepath="/LHCb/Collision17/Beam6500GeV-VeloClosed-"+pol+"/Real Data/Reco17/Stripping29r2p1/90000000/"
elif year == "2018":
    #basepath="/validation/Collision18/Collision18/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco18/Stripping34r0p1/90000000/"
    basepath="/LHCb/Collision18/Beam6500GeV-VeloClosed-"+pol+"/Real Data/Reco18/Stripping34r0p1/90000000/"
else:
    print("unknown year ", year)
    sys.exit()

BK_locations = []
if which == "bhadron":
    BK_locations += [basepath+'BHADRON.MDST']
elif which == "dimuon":
    BK_locations += [basepath+'DIMUON.DST']
elif which == "leptonic":
    BK_locations += [basepath+'LEPTONIC.MDST']
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
        print("no data")
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
    if which is "bhadron":
        j.inputdata = LHCbDataset(['LFN:/lhcb/validation/Collision18/BHADRON.MDST/00086038/0000/00086038_00000027_1.bhadron.mdst'])
    elif which is "dimuon":
        j.inputdata = LHCbDataset(['LFN:/lhcb/validation/Collision18/DIMUON.DST/00086038/0000/00086038_00000031_1.dimuon.dst'])
    elif which is "leptonic":
        j.inputdata = LHCbDataset(['LFN:/lhcb/validation/Collision18/DIMUON.DST/00086038/0000/00086038_00000025_1.leptonic.mdst'])
if note in ["test"]:
    if which is "bhadron":
        j.inputdata = LHCbDataset(['PFN:/tmp/dcraik/00086038_00000027_1.bhadron.mdst'])
    elif which is "dimuon":
        j.inputdata = LHCbDataset(['PFN:/tmp/dcraik/00086038_00000031_1.dimuon.dst'])
    elif which is "leptonic":
        j.inputdata = LHCbDataset(['PFN:/tmp/dcraik/00086038_00000025_1.leptonic.mdst'])

j.parallel_submit = True

queues.add(j.submit)

