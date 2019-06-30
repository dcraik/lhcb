
skip=[]#39,40,41,61]#181,195,202,203,204,205,217,226,242,243,245,247,249,286,287,288,291,292,319,369,370,371,374,375,376,378,465,466,467,502,503,507,508,509,510,511,513,518,520,521,522,523,524,525,526,527,552,556,557,580,581,582,592,593,594,595,604,605,606]
#for jobid in [89,90,91,92,93,94,101,102]:
#for jobid in range(108,138):
#for jobid in range(157,249):
#for jobid in range(248,298):
#for jobid in range(300,316):
#for jobid in range(356,384):
#for jobid in range(379,404):
#for jobid in range(418,423):
#for jobid in range(512,532):
#for jobid in range(545,572):
#for jobid in range(577,607):
#for jobid in range(608,610):
#for jobid in range(633,634):
#for jobid in [4,5,6,7,8,9,10,11,12]:
#for jobid in range(19,79):
for jobid in range(71,79):
#for jobid in range(0,12):
    if jobid in skip: continue

    f=open("lfns/lfns%d.log"%jobid, "w")
    f2=open("lumi/lumi%d.log"%jobid, "w")

    nfound=0
    nsubjobs = len(jobs(jobid).subjobs)

    for i in range(nsubjobs):
        if jobs(jobid).subjobs(i).status=="completed":
            nfound+=1
            for output in jobs(jobid).subjobs(i).outputfiles:
                if output.namePattern=="output.root":
                    f.write("%s\n"%(output.lfn))
                if output.namePattern=="LumiTuple.root":
                    f2.write("%s\n"%(output.lfn))
        else:
            f.write("NOFILE\n")
            f2.write("NOFILE\n")

    f.close()
    f2.close()

    print "job %d (%s): found %d of %d" % (jobid,jobs(jobid).name,nfound,nsubjobs)
