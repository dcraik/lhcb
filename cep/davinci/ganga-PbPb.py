import sys

if len(sys.argv)<2:
    print("Usage: ", sys.argv[0], "<note>")
    sys.exit()
note = str(sys.argv[1])

script=str('DV_PbPb_new.py')

job_name = "CEP phi new PbPb " + str(note)
print(job_name)
print(script)

DV = GaudiExec(directory="/workspace/DaVinciDev_v44r10p2")
DV.options = [script]

BK_locations = []
BK_locations += ["/LHCb/Lead18/Beam6370GeV-VeloClosed-MagDown/Real Data/Reco18Lead18/Stripping35a/90000000/IFT.DST"]
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
        sys.exit()
if note in ["test"]:
        j.inputdata = LHCbDataset(['PFN:/data/cep-phi/00089031_00009979_1.ift.dst'])

j.parallel_submit = True

queues.add(j.submit)

