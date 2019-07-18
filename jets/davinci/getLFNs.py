
skip=[0,1,2,3,4,13,14,17,18,19,20,29]
for jobid in range(42):
    if jobid in skip: continue

    f=open("lfns/lfnsX%d.log"%jobid, "w")
    f2=open("lumi/lumiX%d.log"%jobid, "w")

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
