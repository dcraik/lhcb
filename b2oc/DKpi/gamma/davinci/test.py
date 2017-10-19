#-- GAUDI jobOptions generated on Thu Nov  6 00:26:01 2014
#-- Contains event types : 
#--   90000000 - 3030 files - 241785801 events - 2668.20 GBytes


#--  Extra information about the data processing phases:


#--  Processing Pass Step-124248 

#--  StepId : 124248 
#--  StepName : Stripping20r1-Merging-DV-v32r2p3 
#--  ApplicationName : DaVinci 
#--  ApplicationVersion : v32r2p3 
#--  OptionFiles : $APPCONFIGOPTS/Merging/DV-Stripping-Merging.py 
#--  DDDB : dddb-20130111 
#--  CONDDB : cond-20130114 
#--  ExtraPackages : AppConfig.v3r159 
#--  Visible : N 

from Gaudi.Configuration import * 
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles([
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0005/00022727_00052112_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000014_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000027_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000040_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000053_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000066_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000079_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000092_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000105_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000118_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0000/00022727_00000131_1.bhadron.mdst',
'LFN:/lhcb/LHCb/Collision11/BHADRON.MDST/00022727/0005/00022727_00052127_1.bhadron.mdst'
], clear=True)
