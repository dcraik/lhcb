import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat, correlated_values
import uncertainties as unc
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

###Functions for the jet tagging part###
def loadTagFile(whichDir, whichUnfold):
    data = open("dijets_"+whichDir+"/resultsOutputFull_"+whichUnfold+".log")

    sv = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:8]]),
                    np.array([float(i) for i in data.readline().split(" ")[0:8]]))

    d0 = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:8]]),
                    np.array([float(i) for i in data.readline().split(" ")[0:8]]))

    ignore = data.readline()
    ignore = data.readline()
    ignore = data.readline()
    ignore = data.readline()

    dp = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:8]]),
                    np.array([float(i) for i in data.readline().split(" ")[0:8]]))

    ignore = data.readline()
    ignore = data.readline()
    ignore = data.readline()
    ignore = data.readline()

    return (sv,d0,dp)

def loadTagFiles(svDir, svUnfold, d0Dir, d0Unfold, dpDir, dpUnfold):
    if d0Dir==svDir and d0Unfold==svUnfold:
        if dpDir==svDir and dpUnfold==svUnfold:
            (nsv,nd0,ndp) = loadTagFile(svDir,svUnfold)
        else:
            (nsv,nd0, _ ) = loadTagFile(svDir,svUnfold)
            ( _ , _ ,ndp) = loadTagFile(dpDir,dpUnfold)
    elif dpDir==svDir and dpUnfold==svUnfold:
        (nsv, _ ,ndp) = loadTagFile(svDir,svUnfold)
        ( _ ,nd0, _ ) = loadTagFile(d0Dir,d0Unfold)
    elif dpDir==d0Dir and dpUnfold==d0Unfold:
        (nsv, _ , _ ) = loadTagFile(svDir,svUnfold)
        ( _ ,nd0,ndp) = loadTagFile(d0Dir,d0Unfold)
    else:
        (nsv, _ , _ ) = loadTagFile(svDir,svUnfold)
        ( _ ,nd0, _ ) = loadTagFile(d0Dir,d0Unfold)
        ( _ , _ ,ndp) = loadTagFile(dpDir,dpUnfold)
    return (nsv,nd0,ndp)

def calcEffs(sv,d0,dp,bfffd0,bfffdp,wd0,wdp):
    nc = (wd0*d0/bfffd0 + wdp*dp/bfffdp)/(wd0+wdp)
    return sv/nc

def addTagSystsFromMaxDiff(sv,d0,dp,systDirs, systUnfolds, label=""):
    systSVs = []
    systD0s = []
    systDPs = []
    for (whichDir,whichUnfold) in zip(systDirs,systUnfolds):
        if isinstance(whichDir,list):
            svDir = whichDir[0]
            d0Dir = whichDir[1]
            dpDir = whichDir[2]
        else:
            svDir = whichDir
            d0Dir = whichDir
            dpDir = whichDir
        if isinstance(whichUnfold,list):
            svUnfold = whichUnfold[0]
            d0Unfold = whichUnfold[1]
            dpUnfold = whichUnfold[2]
        else:
            svUnfold = whichUnfold
            d0Unfold = whichUnfold
            dpUnfold = whichUnfold
        (systSVs[len(systSVs):],systD0s[len(systD0s):],systDPs[len(systDPs):]) = tuple(zip(loadTagFiles(svDir,svUnfold,d0Dir,d0Unfold,dpDir,dpUnfold)))

    #print(np.amax(np.stack([np.abs(syst-d0) for syst in systD0s]),0)/d0)
    #print(np.amax(np.stack([np.abs(syst-dp) for syst in systDPs]),0)/dp)
    systsv = unp.uarray(np.array([1.]*len(sv)),
                        np.array([i.n for i in np.amax(np.stack([np.abs(syst-sv) for syst in systSVs]),0)/sv]))
    systd0 = unp.uarray(np.array([1.]*len(d0)),
                        np.array([i.n for i in np.amax(np.stack([np.abs(syst-d0) for syst in systD0s]),0)/d0]))
    systdp = unp.uarray(np.array([1.]*len(dp)),
                        np.array([i.n for i in np.amax(np.stack([np.abs(syst-dp) for syst in systDPs]),0)/dp]))
    if isinstance(label,list):
        svlabel=label[0]
        d0label=label[1]
        dplabel=label[2]
    else:
        svlabel=label
        d0label=label
        dplabel=label
    for i in systsv:
        i.tag=svlabel
    for i in systd0:
        i.tag=d0label
    for i in systdp:
        i.tag=dplabel

    sv = sv*systsv
    d0 = d0*systd0
    dp = dp*systdp

    return (sv,d0,dp)

def addTagSystFromFile(data, filename, label):
    f = open(filename)
    syst = unp.uarray(np.array([1.]*8), np.array([float(i) for i in f.readline().split(" ")[0:8]]))

    for i in syst:
        i.tag=label

    return data*syst

###Functions for the Z+c part###
def loadZjFile(whichDir):
    data = open("zjets_"+whichDir+"/resultsOutput.log")

    sv = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:24]]),
                    np.array([float(i) for i in data.readline().split(" ")[0:24]]))

    nj = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:12]]),
                    np.array([float(i) for i in data.readline().split(" ")[0:12]]))

    eff = unp.uarray(np.array([float(i) for i in data.readline().split(" ")[0:8]]),
                    (np.array([float(i) for i in data.readline().split(" ")[0:8]])**2 +
                     np.array([float(i) for i in data.readline().split(" ")[0:8]])**2)**.5)

    #axes are: 0 - flavour, 1 - pT, 2 - y
    sv = sv.reshape((2,4,3))
    nj = nj.reshape((4,3))
    nj = np.repeat(nj[np.newaxis,:,:], 2, axis=0)
    eff = eff.reshape((2,4))
    eff = np.repeat(eff[:,:,np.newaxis], 3, axis=2)

    return (sv,nj,eff)

def calcRatios(nZsv,nZj,effsv):
    nZQ = nZsv/effsv
    
    nZQ_pt_int = nZQ[:,1:,:].sum(axis=1)
    nZj_pt_int = nZj[:,1:,:].sum(axis=1)
    
    nZQ_y_int = nZQ.sum(axis=2)
    nZj_y_int = nZj.sum(axis=2)
    
    nZQ_int_int = nZQ[:,1:,:].sum(axis=(1,2))
    nZj_int_int = nZj[:,1:,:].sum(axis=(1,2))

    return (nZQ/nZj, nZQ_pt_int/nZj_pt_int, nZQ_y_int/nZj_y_int, nZQ_int_int/nZj_int_int)

def addFlatZjSyst(dif_dif,pt_int,y_int,int_int,syst,tag):
    s = ufloat(1.,syst,tag=tag)
    return tuple(s*z for z in (dif_dif,pt_int,y_int,int_int))

def addZjSystFromMaxDiff(dif_dif,pt_int,y_int,int_int,systDirs,tag):
    dif_dif_systs = []
    pt_int_systs = []
    y_int_systs = []
    int_int_systs = []

    for zjDir in systDirs:
        (s_dd,s_pti,s_yi,s_ii) = calcRatios(*loadZjFile(zjDir))
        dif_dif_systs.append(np.array([i.n for i in ((s_dd - dif_dif)/dif_dif).flat]).reshape(2,4,3))
        pt_int_systs.append(np.array([i.n for i in ((s_pti - pt_int)/pt_int).flat]).reshape(2,3))
        y_int_systs.append(np.array([i.n for i in ((s_yi - y_int)/y_int).flat]).reshape(2,4))
        int_int_systs.append(np.array([i.n for i in ((s_ii - int_int)/int_int).flat]).reshape(2))

    syst_dd = unp.uarray(np.ones((2,4,3)),
                         np.amax(np.stack(np.abs(dif_dif_systs)),0))
    syst_pti = unp.uarray(np.ones((2,3)),
                         np.amax(np.stack(np.abs(pt_int_systs)),0))
    syst_yi = unp.uarray(np.ones((2,4)),
                         np.amax(np.stack(np.abs(y_int_systs)),0))
    syst_ii = unp.uarray(np.ones((2)),
                         np.amax(np.stack(np.abs(int_int_systs)),0))

    for i in syst_dd.flat:
        i.tag = tag
    for i in syst_pti.flat:
        i.tag = tag
    for i in syst_yi.flat:
        i.tag = tag
    for i in syst_ii.flat:
        i.tag = tag

    return (dif_dif*syst_dd,pt_int*syst_pti,y_int*syst_yi,int_int*syst_ii)
    
def addSystsFromMaxInputs(dif_dif,pt_int,y_int,int_int,tagDirs,tagUnfolds,zjDirs,bfffd0,bfffdp,wd0,wdp,tag):
    dif_dif_systs = []
    pt_int_systs = []
    y_int_systs = []
    int_int_systs = []

    for (tagDir,tagUnfold,zjDir) in zip(tagDirs,tagUnfolds,zjDirs):
        (sv,d0,dp) = loadTagFile(tagDir, tagUnfold)
        eff = 0.985*calcEffs(sv,d0,dp,bfffd0,bfffdp,wd0,wdp).reshape((2,4))
        eff = np.repeat(eff[:,:,np.newaxis], 3, axis=2)

        (nZsv,nZj,_) = loadZjFile(zjDir)
        (s_dd,s_pti,s_yi,s_ii) = calcRatios(nZsv,nZj,eff)
        dif_dif_systs.append(np.array([i.n for i in ((s_dd - dif_dif)/dif_dif).flat]).reshape(2,4,3))
        pt_int_systs.append(np.array([i.n for i in ((s_pti - pt_int)/pt_int).flat]).reshape(2,3))
        y_int_systs.append(np.array([i.n for i in ((s_yi - y_int)/y_int).flat]).reshape(2,4))
        int_int_systs.append(np.array([i.n for i in ((s_ii - int_int)/int_int).flat]).reshape(2))

    syst_dd = unp.uarray(np.ones((2,4,3)),
                         np.amax(np.stack(np.abs(dif_dif_systs)),0))
    syst_pti = unp.uarray(np.ones((2,3)),
                         np.amax(np.stack(np.abs(pt_int_systs)),0))
    syst_yi = unp.uarray(np.ones((2,4)),
                         np.amax(np.stack(np.abs(y_int_systs)),0))
    syst_ii = unp.uarray(np.ones((2)),
                         np.amax(np.stack(np.abs(int_int_systs)),0))

    for i in syst_dd.flat:
        i.tag = tag
    for i in syst_pti.flat:
        i.tag = tag
    for i in syst_yi.flat:
        i.tag = tag
    for i in syst_ii.flat:
        i.tag = tag

    return (dif_dif*syst_dd,pt_int*syst_pti,y_int*syst_yi,int_int*syst_ii)

def printZQjSystVertical(dif_dif,pt_int,y_int,int_int,tagDirs,tagUnfolds,zjDirs,bfffd0,bfffdp,wd0,wdp,titles,pt_strings,y_strings):
    dif_dif_alts = []
    pt_int_alts = []
    y_int_alts = []
    int_int_alts = []
    dif_dif_systs = []
    pt_int_systs = []
    y_int_systs = []
    int_int_systs = []

    for (tagDir,tagUnfold,zjDir) in zip(tagDirs,tagUnfolds,zjDirs):
        (sv,d0,dp) = loadTagFile(tagDir, tagUnfold)
        eff = 0.985*calcEffs(sv,d0,dp,bfffd0,bfffdp,wd0,wdp).reshape((2,4))
        eff = np.repeat(eff[:,:,np.newaxis], 3, axis=2)

        (nZsv,nZj,_) = loadZjFile(zjDir)
        (s_dd,s_pti,s_yi,s_ii) = calcRatios(nZsv,nZj,eff)
        dif_dif_alts.append(s_dd)
        pt_int_alts.append(s_pti)
        y_int_alts.append(s_yi)
        int_int_alts.append(s_ii)

        dif_dif_systs.append(np.array([i.n for i in ((s_dd - dif_dif)/dif_dif).flat]).reshape(2,4,3))
        pt_int_systs.append(np.array([i.n for i in ((s_pti - pt_int)/pt_int).flat]).reshape(2,3))
        y_int_systs.append(np.array([i.n for i in ((s_yi - y_int)/y_int).flat]).reshape(2,4))
        int_int_systs.append(np.array([i.n for i in ((s_ii - int_int)/int_int).flat]).reshape(2))

    syst_dd = np.amax(np.stack(np.abs(dif_dif_systs)),0)
    syst_pti = np.amax(np.stack(np.abs(pt_int_systs)),0)
    syst_yi = np.amax(np.stack(np.abs(y_int_systs)),0)
    syst_ii = np.amax(np.stack(np.abs(int_int_systs)),0)

    print(("& & baseline " + len(titles)*"& %s " + "\\\\") % tuple(titles))
    for i in range(len(pt_strings)):
      for j in range(len(y_strings)):
        print("%s & %s & $%5.2f\\pm%5.2f$ " % (pt_strings[i], y_strings[j], 100.*dif_dif[0,i,j].n, 100.*dif_dif[0,i,j].s), end="")
        for alt in dif_dif_alts:
          if abs(alt[0,i,j].n - dif_dif[0,i,j].n)/dif_dif[0,i,j].n >= 0.999*syst_dd[0,i,j]:
            print("& {\\bf %5.2f} " % (100.*alt[0,i,j].n), end="")
          else:
            print("&      %5.2f  " % (100.*alt[0,i,j].n), end="")
        print("& $%4.1f\\,\\%%$ \\\\" % (syst_dd[0,i,j]*100.))
      print("\\midrule")
    for j in range(len(y_strings)):
      print("%s & %s & $%5.2f\\pm%5.2f$ " % ("20--100", y_strings[j], 100.*pt_int[0,j].n, 100.*pt_int[0,j].s), end="")
      for alt in pt_int_alts:
        if abs(alt[0,j].n - pt_int[0,j].n)/pt_int[0,j].n >= 0.999*syst_pti[0,j]:
          print("& {\\bf %5.2f} " % (100.*alt[0,j].n), end="")
        else:
          print("&      %5.2f  " % (100.*alt[0,j].n), end="")
      print("& $%4.1f\\,\\%%$ \\\\" % (syst_pti[0,j]*100.))
    print("\\midrule")
    for i in range(len(pt_strings)):
      print("%s & %s & $%5.2f\\pm%5.2f$ " % (pt_strings[i], "2.00--4.50", 100.*y_int[0,i].n, 100.*y_int[0,i].s), end="")
      for alt in y_int_alts:
        if abs(alt[0,i].n - y_int[0,i].n)/y_int[0,i].n >= 0.999*syst_yi[0,i]:
          print("& {\\bf %5.2f} " % (100.*alt[0,i].n), end="")
        else:
          print("&      %5.2f  " % (100.*alt[0,i].n), end="")
      print("& $%4.1f\\,\\%%$ \\\\" % (syst_yi[0,i]*100.))
    print("\\midrule")
    print("%s & %s & $%5.2f\\pm%5.2f$ " % ("20--100", "2.00--4.50", 100.*int_int[0].n, 100.*int_int[0].s), end="")
    for alt in int_int_alts:
      if abs(alt[0].n - int_int[0].n)/int_int[0].n >= 0.999*syst_ii[0]:
        print("& {\\bf %5.2f} " % (100.*alt[0].n), end="")
      else:
        print("&      %5.2f  " % (100.*alt[0].n), end="")
    print("& $%4.1f\\,\\%%$ \\\\" % (syst_ii[0]*100.))

def printTagEffSystHorizontal(nsv,nd0,ndp,tagDirs,tagUnfolds,bfffd0,bfffdp,wd0,wdp,titles):
    eff_alts = []
    eff_systs = []

    eff = 0.985*calcEffs(nsv,nd0,ndp,bfffd0,bfffdp,wd0,wdp)

    for (tagDir,tagUnfold) in zip(tagDirs,tagUnfolds):
        (s_sv,s_d0,s_dp) = loadTagFile(tagDir, tagUnfold)
        s_eff = 0.985*calcEffs(s_sv,s_d0,s_dp,bfffd0,bfffdp,wd0,wdp)

        eff_alts.append(np.array([i.n for i in s_eff]))
        eff_systs.append(np.array([i.n for i in (s_eff - eff)/eff]))

    syst_eff = np.amax(np.stack(np.abs(eff_systs)),0)

    print(("nominal             " + "& %.3f "*len(eff) + "\\\\") % tuple(e.n for e in eff))
    print("\\midrule")
    for (alt,title) in zip(eff_alts,titles):
      print("%20s" % title, end="")
      for (a,e,s) in zip(alt,eff,syst_eff):
        if abs(a - e)/e >= 0.999*s:
          print("& {\\bf %.3f} " % a, end="")
        else:
          print("&      %.3f  " % a, end="")
      print("\\\\")
    print("\\midrule")
    print(("assigned syst.      " + "& %4.1f\\,\\%% "*len(eff) + "\\\\") % tuple([100.*s for s in syst_eff]))

def printNDSystHorizontal(d0,dp,dirs, unfolds, baseline_title, titles):
    systd0 = np.zeros(len(d0))
    systdp = np.zeros(len(dp))
    print("\midrule")
    print("\multicolumn{%d}{c}{\it %s}\\\\" % (4+1, baseline_title))
    print(("$N(\\cquark\\to\\Dz\\to\Km\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in d0[:4] for e in (a.n, sum(error**2 for (var, error) in a.error_components().items() if var.tag is "statd0")**.5)]))
    print(("$N(\\bquark\\to\\Dz\\to\Km\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in d0[4:] for e in (a.n, sum(error**2 for (var, error) in a.error_components().items() if var.tag is "statd0")**.5)]))
    print(("$N(\\cquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in dp[:4] for e in (a.n, sum(error**2 for (var, error) in a.error_components().items() if var.tag is "statdp")**.5)]))
    print(("$N(\\bquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in dp[4:] for e in (a.n, sum(error**2 for (var, error) in a.error_components().items() if var.tag is "statdp")**.5)]))
    for (svDir,svUnfold,title) in zip(dirs,unfolds,titles):
        (_,altd0,altdp) = loadTagFile(svDir,svUnfold)
        for i in range(len(d0)):
            if abs(altd0[i]-d0[i])/d0[i] > systd0[i]:
                systd0[i] = (abs(altd0[i]-d0[i])/d0[i]).n
            if abs(altdp[i]-dp[i])/dp[i] > systdp[i]:
                systdp[i] = (abs(altdp[i]-dp[i])/dp[i]).n

        print("\midrule")
        print("\multicolumn{%d}{c}{\it %s}\\\\" % (4+1, title))
        print(("$N(\\cquark\\to\\Dz\\to\Km\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in altd0[:4] for e in (a.n, a.s)]))
        print(("$N(\\bquark\\to\\Dz\\to\Km\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in altd0[4:] for e in (a.n, a.s)]))
        print(("$N(\\cquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in altdp[:4] for e in (a.n, a.s)]))
        print(("$N(\\bquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%5.0f \pm %4.0f$ "*4 + "\\\\") % tuple([e for a in altdp[4:] for e in (a.n, a.s)]))
    print("\midrule")
    print("\multicolumn{%d}{c}{\it %s}\\\\" % (4+1, "Difference (assigned as a sytematic uncertainty)"))
    print(("$N(\\cquark\\to\\Dz\\to\Km\\pip)$ " + "& $%4.1f\\,\\%%$ "*4 + "\\\\") % tuple([100.*a for a in systd0[:4]]))
    print(("$N(\\bquark\\to\\Dz\\to\Km\\pip)$ " + "& $%4.1f\\,\\%%$ "*4 + "\\\\") % tuple([100.*a for a in systd0[4:]]))
    print(("$N(\\cquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%4.1f\\,\\%%$ "*4 + "\\\\") % tuple([100.*a for a in systdp[:4]]))
    print(("$N(\\bquark\\to\\Dp\\to\Km\\pip\\pip)$ " + "& $%4.1f\\,\\%%$ "*4 + "\\\\") % tuple([100.*a for a in systdp[4:]]))

def printSystSummary(quants,systCats):
  for (cat,systs) in systCats:
    print("\\midrule")
    if cat!="":
      print("& \multicolumn{%d}{c}{\it %s}\\\\" % (len(quants.flat), cat))
      print("\\midrule")
    for (name,tags) in systs:
      print("%30s " % name, end="")
      for q in quants.flat:
        print(("& %4.1f ") % (100.*sum(error**2 for (var, error) in q.error_components().items() if var.tag in tags)**.5/q.n) , end="")
      print("\\\\")


###Begin main###

defaultDir="20210608"
defaultUnfold="bayes2"
svDir=None
svUnfold=None
d0Dir=None
d0Unfold=None
dpDir=None
dpUnfold=None

zjDir="20210609"

import sys

if len(sys.argv)>1:
    svDir=sys.argv[1]
    if len(sys.argv)>2:
        svUnfold=sys.argv[2]
        if len(sys.argv)>3:
            d0Dir=sys.argv[3]
            if len(sys.argv)>4:
                d0Unfold=sys.argv[4]
                if len(sys.argv)>5:
                    dpDir=sys.argv[5]
                    if len(sys.argv)>6:
                        dpUnfold=sys.argv[6]

if not svDir:
    svDir=defaultDir
if not svUnfold:
    svUnfold=defaultUnfold
if not d0Dir:
    d0Dir=svDir
if not d0Unfold:
    d0Unfold=svUnfold
if not dpDir:
    dpDir=svDir
if not dpUnfold:
    dpUnfold=svUnfold

mpl.use('Agg')

###First load baseline jet tagging results###
bfd0 = np.ones(8)*ufloat(0.03950,0.00031,tag="bfd0")
bfdp = np.ones(8)*ufloat(0.0938,0.0016,tag="bfdp")
(ffc2d0, ffc2dp) = correlated_values([0.6017,0.2409],[[0.0075*0.0075,       -0.66*0.0075*0.0067],
                                                      [-0.66*0.0075*0.0067, 0.0067*0.0067]], tags=["ffd0","ffdp"])
ffb2d0 = ufloat(.577,.021,tag="ffb2d0")
ffb2dp = ufloat(.218,.011,tag="ffb2dp")

ffd0 = np.concatenate((np.ones(4)*ffc2d0,np.ones(4)*ffb2d0))
ffdp = np.concatenate((np.ones(4)*ffc2dp,np.ones(4)*ffb2dp))

bfffd0 = bfd0*ffd0
bfffdp = bfdp*ffdp

(nsv,nd0,ndp) = loadTagFiles(svDir, svUnfold, d0Dir, d0Unfold, dpDir, dpUnfold)

for i in nd0:
    i.tag = "statd0"
for i in ndp:
    i.tag = "statdp"
for i in nsv:
    i.tag = "statSV"

nsv_nosyst = nsv
nd0_nosyst = nd0
ndp_nosyst = ndp

###Load systematic uncertainties on jet tagging###
systDataSimPIDCalib = ufloat(1.,0.020,tag="data-sim & PIDCalib")
nd0 = nd0*systDataSimPIDCalib
ndp = ndp*systDataSimPIDCalib
systVeloErrParam    = ufloat(1.,0.015,tag="VELO error param.")
nd0 = nd0*systVeloErrParam
ndp = ndp*systVeloErrParam

dFitSysts = ["_DSyst_%s"% s for s in ["noptbins","comb","mass","prompt","disp_mean_up","disp_mean_down","disp_width_up","disp_width_down"]]
dFitSystNames = ["$\\pt(D)$ bins","comb. shape","mass width","prompt","displ. mean $+$","displ. mean $-$","displ. width $+$","displ. width $-$"]

(_,nd0,ndp) = addTagSystsFromMaxDiff(nsv,nd0,ndp,[svDir+syst for syst in dFitSysts], len(dFitSysts)*["bayes2"], ["None","D0fit","Dpfit"])
(_,nd0,ndp) = addTagSystsFromMaxDiff(nsv,nd0,ndp,[svDir+"_evtByEvtWeights"], ["bayes2"], ["None","D0weights","Dpweights"])
(_,nd0,ndp) = addTagSystsFromMaxDiff(nsv,nd0,ndp,[svDir+"_3DPIDCalib"], ["bayes2"], "PIDCalibBinning")

nd0 = addTagSystFromFile(nd0, "dijets_"+svDir+"/systMCStats.log", "D0MCstats")
ndp = addTagSystFromFile(ndp, "dijets_"+svDir+"/systMCStats_Dp.log", "DpMCstats")

###Perform combination of tagging from different decay channels###
nc_d0 = nd0/(bfffd0)
nc_dp = ndp/(bfffdp)

#sd0sq = unp.uarray(np.array([x.s**2 for x in nc_d0]),np.zeros(8))
#sdpsq = unp.uarray(np.array([x.s**2 for x in nc_dp]),np.zeros(8))
sd0sq = np.array([x.s**2 for x in nc_d0])
sdpsq = np.array([x.s**2 for x in nc_dp])
##TODO#D0 only
#sd0sq = np.array([0. for x in nc_d0])
#sdpsq = np.array([1. for x in nc_dp])
##TODO#D+ only
#sd0sq = np.array([1. for x in nc_d0])
#sdpsq = np.array([0. for x in nc_dp])

nc = ((sdpsq*nc_d0) + (sd0sq*nc_dp))/(sdpsq + sd0sq)
chisq = (nc-nc_d0)**2/sd0sq + (nc-nc_dp)**2/sdpsq
print([c.n for c in chisq],'\n',chisq.sum().n,'\n',nsv/nc_d0,'\n',nsv/nc_dp,'\n',nsv/nc)

#print(nc_dp/nc_d0)
#for i in (nc_dp/nc_d0):
#  print(i)
#  for (var,val) in i.error_components().items():
#    #print("%s : %.2f" % (var.tag,val/(i.s)))
#    print("%s : %.2f" % (var.tag,val/(i.n)))
#
#for i in (nsv/nc):
#  print(i)
#  for (var,val) in i.error_components().items():
#    if var.tag!="None" and val!=0:
#      print("%s : %.2f" % (var.tag,val/(i.s)))
#
#print(nsv/nc_d0)
#print(nsv/nc_dp)
#
#print([x.n for x in (nc_dp/nc_d0)])

###Output tagging combination tables and plots###
bins = np.array([15.,20.,30.,50.,100.])
mids = (bins[1:] + bins[:-1]) / 2
xerrs = mids - bins[:-1]
offsets = [-2.,2.,0.]
#print(bins,mids,xerrs)
#statErrors = [sum(error**2 for (var, error) in item.error_components().items() if var.tag in ["statd0","statdp"])**.5 for item in (nc_dp/nc_d0)]
nonScaleErrorSources = ["statd0","statdp","D0weights","D0fit","Dpweights","Dpfit"]
nonScaleErrors = [sum(error**2 for (var, error) in item.error_components().items() if var.tag in nonScaleErrorSources)**.5 for item in (nc_dp/nc_d0)]

statErrors = [sum(error**2 for (var, error) in item.error_components().items() if var.tag in ["statd0","statdp"])**.5 for item in (nc_dp/nc_d0)]
nonScaleSystErrors = [sum(error**2 for (var, error) in item.error_components().items() if var.tag in nonScaleErrorSources and var.tag not in ["statd0","statdp"])**.5 for item in (nc_dp/nc_d0)]
scaleSystErrors = [sum(error**2 for (var, error) in item.error_components().items() if var.tag not in nonScaleErrorSources)**.5 for item in (nc_dp/nc_d0)]

print("D0")
for i in (nsv/nc_d0):
    print("& $%.3f \pm %.3f \pm %.3f \pm %.3f$" % (i.n, sum(error**2 for (var, error) in i.error_components().items() if var.tag in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag in nonScaleErrorSources and var.tag not in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag not in nonScaleErrorSources)**.5))
print("D+")
for i in (nsv/nc_dp):
    print("& $%.3f \pm %.3f \pm %.3f \pm %.3f$" % (i.n, sum(error**2 for (var, error) in i.error_components().items() if var.tag in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag in nonScaleErrorSources and var.tag not in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag not in nonScaleErrorSources)**.5))
print("comb")
for i in (nsv/nc):
    print("& $%.3f \pm %.3f \pm %.3f \pm %.3f$" % (i.n, sum(error**2 for (var, error) in i.error_components().items() if var.tag in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag in nonScaleErrorSources and var.tag not in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag not in nonScaleErrorSources)**.5))
print("D+/D0")
for i in (nc_dp/nc_d0):
    print("& $%.3f \pm %.3f \pm %.3f \pm %.3f$" % (i.n, sum(error**2 for (var, error) in i.error_components().items() if var.tag in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag in nonScaleErrorSources and var.tag not in ["statd0","statdp"])**.5,
                                                        sum(error**2 for (var, error) in i.error_components().items() if var.tag not in nonScaleErrorSources)**.5))

sourceGroups = [["statSV","statd0","statdp"],
                ["D0fit","Dpfit"],
                ["D0weights","Dpweights"],
                ["ffd0","ffdp","bfd0","bfdp"],
                ["D0MCstats","DpMCstats"],
                ["VELO error param."],
                ["data-sim & PIDCalib"],
                ["PIDCalibBinning"],
                ]
allSources = sum(sourceGroups, [])

for g in sourceGroups:
    print(g)
    for j in [nsv/nc_d0,nsv/nc_dp,nsv/nc]:
        print("& % 5.1f & % 5.1f & % 5.1f & % 5.1f" % tuple(100.*sum(error**2 for (var, error) in i.error_components().items() if var.tag in g)**.5/i.n for i in j[:4]) )
        #print(g,tuple(sum(error**2 for (var, error) in i.error_components().items() if var.tag in g)**.5 for i in j) )
print("total syst")
for j in [nsv/nc_d0,nsv/nc_dp,nsv/nc]:
    print("& % 5.1f & % 5.1f & % 5.1f & % 5.1f" % tuple(100.*sum(error**2 for (var, error) in i.error_components().items() if var.tag in allSources and var.tag not in ["statSV","statd0","statdp"])**.5/i.n for i in j[:4]) )
print("total")
for j in [nsv/nc_d0,nsv/nc_dp,nsv/nc]:
    print("& % 5.1f & % 5.1f & % 5.1f & % 5.1f" % tuple(100.*sum(error**2 for (var, error) in i.error_components().items() if var.tag in allSources)**.5/i.n for i in j[:4]) )
print("other (should be zero)")
for j in [nsv/nc_d0,nsv/nc_dp,nsv/nc]:
    print("& %.1f & %.1f & %.1f & %.1f" % tuple(100.*sum(error**2 for (var, error) in i.error_components().items() if var.tag not in allSources)**.5/i.n for i in j[:4]) )

plt.rc('text', usetex=True)
fig, axs = plt.subplots(2)
axs[0].set_ylabel("$\\epsilon_{\\rm tag}^{\\rm SV}(c)$")
axs[1].set_ylabel("$N_c(D^+)/N_c(D^0)$")
axs[1].set_xlabel("$p_{\\rm T}(j) [{\\rm GeV}/c]$")
axs[0].errorbar(mids+offsets[0],[y.n for y in (nsv/nc_d0)[:4]],[y.s for y in (nsv/nc_d0)[:4]],(xerrs+offsets[0],xerrs-offsets[0]),'b.', capsize=2)
axs[0].errorbar(mids+offsets[1],[y.n for y in (nsv/nc_dp)[:4]],[y.s for y in (nsv/nc_dp)[:4]],(xerrs+offsets[1],xerrs-offsets[1]),'r.', capsize=2)
axs[0].errorbar(mids+offsets[2],[y.n for y in (nsv/nc)[:4]],   [y.s for y in (nsv/nc)[:4]],   (xerrs+offsets[2],xerrs-offsets[2]),'m.', capsize=2)
axs[1].errorbar(mids,[y.n for y in (nc_dp/nc_d0)[:4]],nonScaleErrors[:4],xerrs,'k.', capsize=2)
axs[1].errorbar(mids,[y.n for y in (nc_dp/nc_d0)[:4]],[y.s for y in (nc_dp/nc_d0)[:4]],xerrs,'k.', capsize=2)
axs[1].axhline(y=1.,linestyle="--")
plt.savefig("jetTagD0DpCompWithComb.pdf")

fig, axs = plt.subplots(2)
axs[0].set_ylabel("$\\epsilon_{\\rm tag}^{\\rm SV}(b)$")
axs[1].set_ylabel("$N_b(D^+)/N_b(D^0)$")
axs[1].set_xlabel("$p_{\\rm T}(j) [{\\rm GeV}/c]$")
axs[0].set_ylim([0.3,1.0])
axs[1].set_ylim([0.0,2.0])
axs[0].errorbar(mids,[y.n for y in (nsv/nc_d0)[4:]],[y.s for y in (nsv/nc_d0)[4:]],xerrs,'b.', capsize=2)
axs[0].errorbar(mids,[y.n for y in (nsv/nc_dp)[4:]],[y.s for y in (nsv/nc_dp)[4:]],xerrs,'r.', capsize=2)
axs[0].errorbar(mids,[y.n for y in (nsv/nc)[4:]],   [y.s for y in (nsv/nc)[4:]],   xerrs,'m.', capsize=2)
axs[1].errorbar(mids,[y.n for y in (nc_dp/nc_d0)[4:]],nonScaleErrors[4:],xerrs,'k.', capsize=2)
axs[1].errorbar(mids,[y.n for y in (nc_dp/nc_d0)[4:]],[y.s for y in (nc_dp/nc_d0)[4:]],xerrs,'k.', capsize=2)
axs[1].axhline(y=1.,linestyle="--")
plt.savefig("jetTagD0DpCompWithComb_beauty.pdf")

###Write out tagging results for use by Zj fits###
write = open("jetTagComb.txt","w")
for i in (nsv/nc):
    write.write("%.8f %.8f %.8f\n" % (i.n,
                                      sum(error**2 for (var, error) in i.error_components().items() if var.tag in ["statsv","statd0","statdp"])**.5,
                                      sum(error**2 for (var, error) in i.error_components().items() if var.tag not in ["statsv","statd0","statdp"])**.5))
write.close()

###Start the Z+c part###

effsv_fullerr = 0.985*(nsv/nc).reshape((2,4))
effsv_fullerr = np.repeat(effsv_fullerr[:,:,np.newaxis], 3, axis=2)

###Here we must worry about correlations between pT bins so we need to reapply some of the systematics after the integrations###
nsv = nsv_nosyst
nd0 = nd0_nosyst
ndp = ndp_nosyst

#These systematics are common to the pT bins so no need to treat differently
nd0 = nd0*systDataSimPIDCalib
ndp = ndp*systDataSimPIDCalib
nd0 = nd0*systVeloErrParam
ndp = ndp*systVeloErrParam

#These systematics are uncorrelated between pT bins so no need to treat differently
nd0 = addTagSystFromFile(nd0, "dijets_"+svDir+"/systMCStats.log", "D0MCstats")
ndp = addTagSystFromFile(ndp, "dijets_"+svDir+"/systMCStats_Dp.log", "DpMCstats")

#D fitting and D weighting remain - we will apply these directly to nQ

#TODO#dFitSysts = ["20210608_DSyst_%s"% s for s in ["comb","prompt","disp_width_up","disp_width_down","disp_mean_up","disp_mean_down","mass","noptbins"]]

#(_,nd0,ndp) = addTagSystsFromMaxDiff(nsv,nd0,ndp,dFitSysts, len(dFitSysts)*["bayes2"], ["None","D0fit","Dpfit"])
#(_,nd0,ndp) = addTagSystsFromMaxDiff(nsv,nd0,ndp,["20210608_evtByEvtWeights"], ["bayes2"], ["None","D0weights","Dpweights"])

nc_d0 = nd0/(bfffd0)
nc_dp = ndp/(bfffdp)

nc = ((sdpsq*nc_d0) + (sd0sq*nc_dp))/(sdpsq + sd0sq)
#nc = nc_d0
#nc = nc_dp


###Load baseline results and shape to (flav,pT,y) ###
effsv = 0.985*(nsv/nc).reshape((2,4))
effsv = np.repeat(effsv[:,:,np.newaxis], 3, axis=2)
(nZsv,nZj,_) = loadZjFile(zjDir)

for i in nZsv.flat:
    i.tag = "ZSVstat"
for i in nZj.flat:
    i.tag = "Zjstat"

#Efficiency correct
#nZQ = nZsv/effsv
#
#nZQ_pt_int = nZQ[:,1:,:].sum(axis=1)
#nZj_pt_int = nZj[:,1:,:].sum(axis=1)
#
#nZQ_y_int = nZQ.sum(axis=2)
#nZj_y_int = nZj.sum(axis=2)
#
#nZQ_int_int = nZQ[:,1:,:].sum(axis=(1,2))
#nZj_int_int = nZj[:,1:,:].sum(axis=(1,2))
#
#zQjSysts         = np.ones((2,4,3))
#zQjSysts_pt_int  = np.ones((2,3))
#zQjSysts_y_int   = np.ones((2,4))
#zQjSysts_int_int = np.ones((2))

(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = calcRatios(nZsv,nZj,effsv)
#print((zQj,zQj_pt_int,zQj_y_int,zQj_int_int))

###Load systematic uncertainties on Z+j###
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addFlatZjSyst(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,0.01,"jetReco")

#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addZjSystFromMaxDiff(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,[zjDir+"_jetEnergyScaleUp",zjDir+"_jetEnergyScaleDown",zjDir+"_jetEnergySmearUp",zjDir+"_jetEnergySmearDown"],"jetEnergy")
#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addZjSystFromMaxDiff(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,[zjDir+"_svSyst_bkgrnd",zjDir+"_svSyst_sig"],"SV fit")
#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addZjSystFromMaxDiff(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,[zjDir+"_svSyst_bkgrnd"],"SV fit (b)")
#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addZjSystFromMaxDiff(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,[zjDir+"_svSyst_sig"],"SV fit (s)")
jeSysts=["_jetEnergyScaleUp","_jetEnergyScaleDown","_jetEnergySmearUp","_jetEnergySmearDown"]
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                                                               [svDir+syst for syst in jeSysts],len(jeSysts)*["bayes2"],
                                                               [zjDir+syst for syst in jeSysts],
                                                               bfffd0,bfffdp,sdpsq,sd0sq,"jetEnergy")#Note weights are sdpsq for D0 and sd0sq for D+
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                                                               [svDir+"_svSyst_bkgrnd",svDir+"_svSyst_sig"],2*["bayes2"],
                                                               [zjDir+"_svSyst_bkgrnd",zjDir+"_svSyst_sig"],
                                                               bfffd0,bfffdp,sdpsq,sd0sq,"SV fit")#Note weights are sdpsq for D0 and sd0sq for D+
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                                                               [svDir+syst for syst in dFitSysts],len(dFitSysts)*["bayes2"],
                                                               len(dFitSysts)*[zjDir],
                                                               bfffd0,bfffdp,sdpsq,sd0sq,"D fits")#Note weights are sdpsq for D0 and sd0sq for D+
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                                                               [svDir+"_evtByEvtWeights"],["bayes2"],
                                                               [zjDir],
                                                               bfffd0,bfffdp,sdpsq,sd0sq,"D weights")#Note weights are sdpsq for D0 and sd0sq for D+
(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                                                               [svDir+"_3DPIDCalib"],["bayes2"],
                                                               [zjDir],
                                                               bfffd0,bfffdp,sdpsq,sd0sq,"PIDCalibBinning")#Note weights are sdpsq for D0 and sd0sq for D+

#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
#                                                               [svDir+"_svSyst_sig"],["bayes2"],
#                                                               [zjDir+"_svSyst_sig"],
#                                                               bfffd0,bfffdp,sdpsq,sd0sq,"SV fit")#Note weights are sdpsq for D0 and sd0sq for D+
#(zQj,zQj_pt_int,zQj_y_int,zQj_int_int) = addSystsFromMaxInputs(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
#                                                               [svDir+"_svSyst_bkgrnd"],["bayes2"],
#                                                               [zjDir+"_svSyst_bkgrnd"],
#                                                               bfffd0,bfffdp,sdpsq,sd0sq,"SV fitb")#Note weights are sdpsq for D0 and sd0sq for D+
###Results###
#print((zQj,zQj_pt_int,zQj_y_int,zQj_int_int))

sCats = [
("stat", ["Zjstat","ZSVstat"]),
("tag stat", ["statSV","statd0","statdp"]),
("tag bfff", ["ffd0","ffdp","bfd0","bfdp"]),
("tag syst", ["D0MCstats","DpMCstats","D fits","D weights","data-sim & PIDCalib","VELO error param.","PIDCalibBinning"]),
("jet reco", ["jetReco"]),
("jet energy", ["jetEnergy"]),
("sv fit", ["SV fit"]),
]

print(zQj_pt_int)
for i in ((zQj_pt_int).flat):
  print(i)
  for (cat,tags) in sCats:
    print("%s: %5.2f%%" % (cat, (100.*(sum(error**2 for (var, error) in i.error_components().items() if var.tag in tags))**.5)/i.n), end="; ")
  #tags = set()
  #for (var, error) in i.error_components().items():
  #  tags.add(var.tag)
  #for tag in tags:
  #  print("%s: %5.2f%%" % (tag, (100.*(sum(error**2 for (var, error) in i.error_components().items() if var.tag == tag))**.5)/i.n), end="; ")
  print("%s: %5.2f%%" % ("total", (100.*(sum(error**2 for (var, error) in i.error_components().items()))**.5)/i.n))

output = open("results.txt","w")
for item in zQj_pt_int.flat:
    output.write("%f " % item.n)
output.write("\n")
for item in zQj_pt_int.flat:
    output.write("%f " % (sum(error**2 for (var, error) in item.error_components().items() if var.tag in ["Zjstat","ZSVstat"]))**.5)
output.write("\n")
for item in zQj_pt_int.flat:
    output.write("%f " % item.s)
output.write("\n")
output.close()

tagsStat = ["Zjstat","ZSVstat"]
tagsTagStat = ["statSV","statd0","statdp"]
tagsTagSyst = ["D0MCstats","DpMCstats","D fits","D weights","data-sim & PIDCalib","VELO error param.","PIDCalibBinning"]
tagsSyst = ["jetReco","jetEnergy","SV fit"]

ptlabel=["15-- 20","20-- 30","30-- 50","50--100"]
ylabel=["2.00--2.75","2.75--3.50","3.50--4.50"]
for flav in [0,1]:
  for pt in range(0,4):
    print("\\midrule")
    for y in range(0,3):
      print("%s & %s & $%6.0f\\pm%4.0f$ & $%4.0f\\pm%3.0f$ & $%.3f\\pm%.3f\\pm%.3f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\" % (ptlabel[pt],ylabel[y],
                                                                                       nZj[flav,pt,y].n,nZj[flav,pt,y].s,
                                                                                       nZsv[flav,pt,y].n,nZsv[flav,pt,y].s,
                                                                                       effsv_fullerr[flav,pt,y].n,
                                                                                       sum(error**2 for (var, error) in effsv_fullerr[flav,pt,y].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                       sum(error**2 for (var, error) in effsv_fullerr[flav,pt,y].error_components().items() if var.tag not in tagsTagStat)**.5,
                                                                                       100.*zQj[flav,pt,y].n,
                                                                                       100.*sum(error**2 for (var, error) in zQj[flav,pt,y].error_components().items() if var.tag in tagsStat)**.5,
                                                                                       100.*sum(error**2 for (var, error) in zQj[flav,pt,y].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                       100.*sum(error**2 for (var, error) in zQj[flav,pt,y].error_components().items() if var.tag in tagsTagSyst)**.5,
                                                                                       100.*sum(error**2 for (var, error) in zQj[flav,pt,y].error_components().items() if var.tag in tagsSyst)**.5,
                                                                                       ))
  print("\\midrule")
  for y in range(0,3):
    print("%s & %s & $%6.0f\\pm%4.0f$ & $%4.0f\\pm%3.0f$ & --- & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\" % ("20--100",ylabel[y],
                                                                                     nZj[:,1:,:].sum(axis=1)[flav,y].n,nZj[:,1:,:].sum(axis=1)[flav,y].s,
                                                                                     nZsv[:,1:,:].sum(axis=1)[flav,y].n,nZsv[:,1:,:].sum(axis=1)[flav,y].s,
                                                                                     100.*zQj_pt_int[flav,y].n,
                                                                                     100.*sum(error**2 for (var, error) in zQj_pt_int[flav,y].error_components().items() if var.tag in tagsStat)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_pt_int[flav,y].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_pt_int[flav,y].error_components().items() if var.tag in tagsTagSyst)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_pt_int[flav,y].error_components().items() if var.tag in tagsSyst)**.5,
                                                                                     ))
  print("\\midrule")
  for pt in range(0,4):
    print("%s & %s & $%6.0f\\pm%4.0f$ & $%4.0f\\pm%3.0f$ & $%.3f\\pm%.3f\\pm%.3f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\" % (ptlabel[pt],"2.00--4.50",
                                                                                     nZj[:,:,:].sum(axis=2)[flav,pt].n,nZj[:,:,:].sum(axis=2)[flav,pt].s,
                                                                                     nZsv[:,:,:].sum(axis=2)[flav,pt].n,nZsv[:,:,:].sum(axis=2)[flav,pt].s,
                                                                                     effsv_fullerr[flav,pt,0].n,
                                                                                     sum(error**2 for (var, error) in effsv_fullerr[flav,pt,0].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                     sum(error**2 for (var, error) in effsv_fullerr[flav,pt,0].error_components().items() if var.tag not in tagsTagStat)**.5,
                                                                                     100.*zQj_y_int[flav,pt].n,
                                                                                     100.*sum(error**2 for (var, error) in zQj_y_int[flav,pt].error_components().items() if var.tag in tagsStat)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_y_int[flav,pt].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_y_int[flav,pt].error_components().items() if var.tag in tagsTagSyst)**.5,
                                                                                     100.*sum(error**2 for (var, error) in zQj_y_int[flav,pt].error_components().items() if var.tag in tagsSyst)**.5,
                                                                                     ))
  print("\\midrule")
  print("%s & %s & $%6.0f\\pm%4.0f$ & $%4.0f\\pm%3.0f$ & --- & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\" % ("20--100","2.00--4.50",
                                                                                   nZj[:,1:,:].sum(axis=(1,2))[flav].n,nZj[:,1:,:].sum(axis=(1,2))[flav].s,
                                                                                   nZsv[:,1:,:].sum(axis=(1,2))[flav].n,nZsv[:,1:,:].sum(axis=(1,2))[flav].s,
                                                                                   100.*zQj_int_int[flav].n,
                                                                                   100.*sum(error**2 for (var, error) in zQj_int_int[flav].error_components().items() if var.tag in tagsStat)**.5,
                                                                                   100.*sum(error**2 for (var, error) in zQj_int_int[flav].error_components().items() if var.tag in tagsTagStat)**.5,
                                                                                   100.*sum(error**2 for (var, error) in zQj_int_int[flav].error_components().items() if var.tag in tagsTagSyst)**.5,
                                                                                   100.*sum(error**2 for (var, error) in zQj_int_int[flav].error_components().items() if var.tag in tagsSyst)**.5,
                                                                                   ))

ybins = np.array([2.00,2.75,3.50,4.50])
ymids = (ybins[1:] + ybins[:-1]) / 2
yerrs = ymids - ybins[:-1]
yoffset = [-.15,-.05,.05,.15]

fig, axs = plt.subplots(1)
axs.set_ylabel("$\\it{\\sigma}_{\\it{Z}+\\it{c}}/\\it{\\sigma}_{\\it{Z}+\\it{j}}$",fontsize=20)
axs.set_xlabel("$\\it{y}_{\\it{Z}}$",fontsize=20)
axs.set_ylim([0.,0.12])
axs.tick_params(axis="both", which="both", direction="in", bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="major", length=7, labelsize=12, bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="minor", length=3, bottom=True, top=True, left=True, right=True)
axs.xaxis.set_minor_locator(AutoMinorLocator())
axs.yaxis.set_minor_locator(AutoMinorLocator())
axs.errorbar(ymids+yoffset[0],[y.n for y in (zQj)[0,0,:]],[y.s for y in (zQj)[0,0,:]],(yerrs+yoffset[0],yerrs-yoffset[0]),'b.', capsize=2)
axs.errorbar(ymids+yoffset[0],[y.n for y in (zQj)[0,0,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[0,0,:]**.5],(yerrs+yoffset[0],yerrs-yoffset[0]),'b.', capsize=2)
axs.errorbar(ymids+yoffset[1],[y.n for y in (zQj)[0,1,:]],[y.s for y in (zQj)[0,1,:]],(yerrs+yoffset[1],yerrs-yoffset[1]),'g.', capsize=2)
axs.errorbar(ymids+yoffset[1],[y.n for y in (zQj)[0,1,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[0,1,:]**.5],(yerrs+yoffset[1],yerrs-yoffset[1]),'g.', capsize=2)
axs.errorbar(ymids+yoffset[2],[y.n for y in (zQj)[0,2,:]],[y.s for y in (zQj)[0,2,:]],(yerrs+yoffset[2],yerrs-yoffset[2]),'r.', capsize=2)
axs.errorbar(ymids+yoffset[2],[y.n for y in (zQj)[0,2,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[0,2,:]**.5],(yerrs+yoffset[2],yerrs-yoffset[2]),'r.', capsize=2)
axs.errorbar(ymids+yoffset[3],[y.n for y in (zQj)[0,3,:]],[y.s for y in (zQj)[0,3,:]],(yerrs+yoffset[3],yerrs-yoffset[3]),'m.', capsize=2)
axs.errorbar(ymids+yoffset[3],[y.n for y in (zQj)[0,3,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[0,3,:]**.5],(yerrs+yoffset[3],yerrs-yoffset[3]),'m.', capsize=2)
plt.savefig("charmFracZj_sep.pdf")

fig, axs = plt.subplots(1)
axs.set_ylabel("$\\it{\\sigma}_{\\it{Z}+\\it{c}}/\\it{\\sigma}_{\\it{Z}+\\it{j}}$",fontsize=20)
axs.set_xlabel("$\\it{y}_{\\it{Z}}$",fontsize=20)
axs.set_ylim([0.,0.12])
axs.tick_params(axis="both", which="both", direction="in", bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="major", length=7, labelsize=12, bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="minor", length=3, bottom=True, top=True, left=True, right=True)
axs.xaxis.set_minor_locator(AutoMinorLocator())
axs.yaxis.set_minor_locator(AutoMinorLocator())
axs.errorbar(ymids,[y.n for y in (zQj_pt_int)[0,:]],[y.s for y in (zQj_pt_int)[0,:]],yerrs,'b.', capsize=2)
axs.errorbar(ymids,[y.n for y in (zQj_pt_int)[0,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj_pt_int)[0,:]**.5],yerrs,'k.', capsize=2)
plt.savefig("charmFracZj_combined.pdf")

fig, axs = plt.subplots(1)
axs.set_ylabel("$\\it{\\sigma}_{\\it{Z}+\\it{b}}/\\it{\\sigma}_{\\it{Z}+\\it{j}}$",fontsize=20)
axs.set_xlabel("$\\it{y}_{\\it{Z}}$",fontsize=20)
axs.set_ylim([0.,0.12])
axs.tick_params(axis="both", which="both", direction="in", bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="major", length=7, labelsize=12, bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="minor", length=3, bottom=True, top=True, left=True, right=True)
axs.xaxis.set_minor_locator(AutoMinorLocator())
axs.yaxis.set_minor_locator(AutoMinorLocator())
axs.errorbar(ymids+yoffset[0],[y.n for y in (zQj)[1,0,:]],[y.s for y in (zQj)[1,0,:]],(yerrs+yoffset[0],yerrs-yoffset[0]),'b.', capsize=2)
axs.errorbar(ymids+yoffset[0],[y.n for y in (zQj)[1,0,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[1,0,:]**.5],(yerrs+yoffset[0],yerrs-yoffset[0]),'b.', capsize=2)
axs.errorbar(ymids+yoffset[1],[y.n for y in (zQj)[1,1,:]],[y.s for y in (zQj)[1,1,:]],(yerrs+yoffset[1],yerrs-yoffset[1]),'g.', capsize=2)
axs.errorbar(ymids+yoffset[1],[y.n for y in (zQj)[1,1,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[1,1,:]**.5],(yerrs+yoffset[1],yerrs-yoffset[1]),'g.', capsize=2)
axs.errorbar(ymids+yoffset[2],[y.n for y in (zQj)[1,2,:]],[y.s for y in (zQj)[1,2,:]],(yerrs+yoffset[2],yerrs-yoffset[2]),'r.', capsize=2)
axs.errorbar(ymids+yoffset[2],[y.n for y in (zQj)[1,2,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[1,2,:]**.5],(yerrs+yoffset[2],yerrs-yoffset[2]),'r.', capsize=2)
axs.errorbar(ymids+yoffset[3],[y.n for y in (zQj)[1,3,:]],[y.s for y in (zQj)[1,3,:]],(yerrs+yoffset[3],yerrs-yoffset[3]),'m.', capsize=2)
axs.errorbar(ymids+yoffset[3],[y.n for y in (zQj)[1,3,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj)[1,3,:]**.5],(yerrs+yoffset[3],yerrs-yoffset[3]),'m.', capsize=2)
plt.savefig("beautyFracZj_sep.pdf")

fig, axs = plt.subplots(1)
axs.set_ylabel("$\\it{\\sigma}_{\\it{Z}+\\it{b}}/\\it{\\sigma}_{\\it{Z}+\\it{j}}$",fontsize=20)
axs.set_xlabel("$\\it{y}_{\\it{Z}}$",fontsize=20)
axs.set_ylim([0.,0.12])
axs.tick_params(axis="both", which="both", direction="in", bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="major", length=7, labelsize=12, bottom=True, top=True, left=True, right=True)
axs.tick_params(axis="both", which="minor", length=3, bottom=True, top=True, left=True, right=True)
axs.xaxis.set_minor_locator(AutoMinorLocator())
axs.yaxis.set_minor_locator(AutoMinorLocator())
axs.errorbar(ymids,[y.n for y in (zQj_pt_int)[1,:]],[y.s for y in (zQj_pt_int)[1,:]],yerrs,'b.', capsize=2)
axs.errorbar(ymids,[y.n for y in (zQj_pt_int)[1,:]],[sum(error**2 for (var, error) in y.error_components().items() if var.tag in tagsStat)**.5 for y in (zQj_pt_int)[1,:]**.5],yerrs,'k.', capsize=2)
plt.savefig("beautyFracZj_combined.pdf")

#Syst tables
print("TABLE28")
printZQjSystVertical(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                     [svDir+syst for syst in jeSysts],len(jeSysts)*["bayes2"],
                     [zjDir+syst for syst in jeSysts],
                     bfffd0,bfffdp,sdpsq,sd0sq,["scale+","scale-","smear+","smear-"],ptlabel,ylabel)#Note weights are sdpsq for D0 and sd0sq for D+
print("TABLE29")
printNDSystHorizontal(nd0,ndp,[svDir+"_evtByEvtWeights"], ["bayes2"], "Average-efficiency weighted", ["Event-by-event weighted"])
print("NEW TABLE 31")
printZQjSystVertical(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                     [svDir+"_3DPIDCalib"],["bayes2"],
                     [zjDir],
                     bfffd0,bfffdp,sdpsq,sd0sq,["track mult. binned"],ptlabel,ylabel)#Note weights are sdpsq for D0 and sd0sq for D+
print("TABLE31")
printTagEffSystHorizontal(nsv,nd0,ndp,[svDir+syst for syst in dFitSysts],len(dFitSysts)*["bayes2"],bfffd0,bfffdp,sdpsq,sd0sq,dFitSystNames)
                                                               
print("TABLE33")
printZQjSystVertical(zQj,zQj_pt_int,zQj_y_int,zQj_int_int,
                     [svDir+"_svSyst_bkgrnd",svDir+"_svSyst_sig"],2*["bayes2"],
                     [zjDir+"_svSyst_bkgrnd",zjDir+"_svSyst_sig"],
                     bfffd0,bfffdp,sdpsq,sd0sq,["c/b shapes", "mis-tag shape"],ptlabel,ylabel)#Note weights are sdpsq for D0 and sd0sq for D+
print("TABLE34")
systCats = [
("", [
    ("stat", ["Zjstat","ZSVstat"]),
]),
("Jet Reconstruction Efficiency", [
    ("jet reco. eff.", ["jetReco"]),
]),
("Jet Energy Scale \\& Resolution", [
    ("jet energy scale \\& res.", ["jetEnergy"]),
]),
("Flavour Tagging Efficiency", [
    ("dijet calib. sample size", ["statSV","statd0","statdp"]),
    ("BF \\& FF", ["ffd0","ffdp","bfd0","bfdp"]),
    ("weight method", ["D weights"]),
    ("MC sample size", ["D0MCstats","DpMCstats"]),
    ("data--sim \\& PIDCalib", ["data-sim & PIDCalib","PIDCalibBinning"]),
    ("$D$ fits", ["D fits"]),
    ("VELO error param.", ["VELO error param."]),
]),
("SV Fits", [
    ("SV fit", ["SV fit"]),
]),
("Totals", [
("tagging stat.", ["statSV","statd0","statdp"]),
("tagging syst.", ["ffd0","ffdp","bfd0","bfdp","D0MCstats","DpMCstats","D fits","D weights","data-sim & PIDCalib","VELO error param.","PIDCalibBinning"]),
("other syst.", ["jetReco","jetEnergy","SV fit"]),
]),
("", [
("total syst.", ["statSV","statd0","statdp","ffd0","ffdp","bfd0","bfdp","D0MCstats","DpMCstats","D fits","D weights","data-sim & PIDCalib","VELO error param.","PIDCalibBinning","jetReco","jetEnergy","SV fit"]),
("total", ["Zjstat","ZSVstat","statSV","statd0","statdp","ffd0","ffdp","bfd0","bfdp","D0MCstats","DpMCstats","D fits","D weights","data-sim & PIDCalib","VELO error param.","PIDCalibBinning","jetReco","jetEnergy","SV fit"]),
])
]
printSystSummary(zQj_pt_int,systCats)
print("\n")
printSystSummary(zQj,systCats)
print("\n")
printSystSummary(zQj_y_int,systCats)
print("\n")
printSystSummary(zQj_int_int,systCats)
