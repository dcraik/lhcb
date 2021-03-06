import sys

if len(sys.argv)<3:
    print("Usage: ", sys.argv[0], "<pPb|Pbp> <note>")
    sys.exit()
dset = str(sys.argv[1]) #pPb Pbp
note = str(sys.argv[2])

#script=str('DV_pPb16.py')
script=str('DV_pPb_new.py')

job_name = "CEP phi new " + str(dset) + " " + str(note)
print(job_name)
print(script)

DV = GaudiExec(directory="/workspace/DaVinciDev_v44r10p2")
DV.options = [script]

if dset not in ['pPb','Pbp']:
    print("bad dataset", dset)
    sys.exit()

BK_locations = []
if dset=='Pbp':
    BK_locations += ['/LHCb/Ionproton16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16pLead/Stripping30r2/90000000/IFT.DST']
elif dset=='pPb':
    BK_locations +=['/LHCb/Protonion16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16pLead/Stripping30r3/90000000/IFT.DST']
data = LHCbDataset()
bk = BKQuery(dqflag=['OK','UNCHECKED'])

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
        #j.inputdata = LHCbDataset(['LFN:/lhcb/LHCb/Collision16/EW.DST/00069603/0000/00069603_00001133_1.ew.dst'])
        print("no grid test LFN set")
        sys.exit()
if note in ["test"]:
        j.inputdata = LHCbDataset(['PFN:/data/cep-phi/00076144_00000038_1.ift.dst'])

j.parallel_submit = True

queues.add(j.submit)

