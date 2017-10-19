#b = Bender(version="v22r3", platform = 'x86_64-slc5-gcc46-opt')
#dv = DaVinci(version="v33r4", platform = 'x86_64-slc5-gcc46-opt')
dv = DaVinci(version="v34r1", platform = 'x86_64-slc6-gcc48-opt')
#dv.module = File("/afs/cern.ch/user/d/dcraik/strip20/new_Data_read.py")

t = JobTemplate(name="new_Data_read", application = dv)
#t.application.optsfile = "/home/phrkbf/B2DhhAnalysis/DCP/new_Data_read_Dpipi_D2KK.py"
t.application.optsfile = "/home/phrkbf/B2DhhAnalysis/DCP/new_Data_read_Dpipi_D2pipi.py"

#ds = dv.readInputData("2012/LHCbCollision12Beam4000GeV-VeloClosed-MagDownRealDataReco14Stripping2090000000BHADRONMDST.py");
#ds2 = dv.readInputData("2012/LHCbCollision12Beam4000GeV-VeloClosed-MagUpRealDataReco14Stripping2090000000BHADRONMDST.py");
#ds3 = dv.readInputData("2012/LHCbCollision1190000000Beam3500GeV-VeloClosed-MagDownRealDataReco14Stripping20r1BHADRONMDST_new.py");
#ds4 = dv.readInputData("2012/LHCbCollision1190000000Beam3500GeV-VeloClosed-MagUpRealDataReco14Stripping20r1BHADRONMDST_new.py");
#ds3 = dv.readInputData("2012/MDremaining.py");
#ds4 = dv.readInputData("2012/MUremaining.py");
#ds5 = dv.readInputData("2012/MDpart.py");

#ds  = dv.readInputData("test.py");
#ds  = dv.readInputData("BK_LHCb_Collision11_Beam3500GeV-VeloClosed-MagDown_RealData_Reco14_Stripping20r1_90000000_BHADRON.MDST.py");
#ds2 = dv.readInputData("BK_LHCb_Collision11_Beam3500GeV-VeloClosed-MagUp_RealData_Reco14_Stripping20r1_90000000_BHADRON.MDST.py");
ds  = dv.readInputData("BK_LHCb_Collision12_Beam4000GeV-VeloClosed-MagDown_RealData_Reco14_Stripping20_90000000_BHADRON.MDST.py");
ds2 = dv.readInputData("BK_LHCb_Collision12_Beam4000GeV-VeloClosed-MagUp_RealData_Reco14_Stripping20_90000000_BHADRON.MDST.py");

#mypfn = PhysicalFile('/data/lhcb/phrkbf/2012DST/00020198_00000014_1.bhadron.mdst')
#ds6 = LHCbDataset()
#ds6.files = [mypfn]

j = Job(t, 
#  backend = Local(), 
  backend = Dirac(),
  inputdata = ds, 
  outputfiles = [LocalFile('B_Data_nTuples.root'),LocalFile('summary.xml')], 
  splitter = SplitByFiles( filesPerJob = 20 ), 
)
#j.backend.diracOpts = 'j.setBannedSites(["LCG.CERN.ch"])'
queues.add(j.submit)

k = Job(t,
#  backend = Local(),
  backend = Dirac(),
  inputdata = ds2,
  outputfiles = [LocalFile('B_Data_nTuples.root'),LocalFile('summary.xml')],
  splitter = SplitByFiles( filesPerJob = 20 ),
)
##k.backend.diracOpts = 'j.setBannedSites(["LCG.GRIDKA.de"])'
queues.add(k.submit)
#
#i = Job(t,
#  backend = Dirac(),
#  inputdata = ds5,
#  outputfiles = ['B_Data_nTuples.root','summary.xml'],
#  splitter = SplitByFiles( filesPerJob = 1 ),
#)
##h.backend.diracOpts = 'j.setBannedSites(["LCG.GRIDKA.de"])'
#i.submit()

#i = Job(t,
#  backend = Dirac(),
#  inputdata = ds3,
#  outputdata = ['DKpi_BdBs_Fit_2011.root'],
#  splitter = DiracSplitter( filesPerJob = 20 ),
#)
#i.backend.diracOpts = 'j.setBannedSites(["LCG.GRIDKA.de"])'
#i.submit()
#
#k = Job(t,
#  backend = Dirac(),
#  inputdata = ds4,
#  outputdata = ['DKpi_BdBs_Fit_2011.root'],
#  splitter = DiracSplitter( filesPerJob = 20 ),
#)
#k.backend.diracOpts = 'j.setBannedSites(["LCG.GRIDKA.de"])'
#k.submit()

