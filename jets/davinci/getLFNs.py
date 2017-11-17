
skip=[181,195,202,203,204,205,217,226,242,243,245,247,249,286,287,288,291,292,319]
#for jobid in [89,90,91,92,93,94,101,102]:
#for jobid in range(108,138):
#for jobid in range(157,249):
#for jobid in range(248,298):
#for jobid in range(300,316):
for jobid in range(320,326):
    if jobid in skip: continue

    f=open("lfns/lfns%d.log"%jobid, "w")
    
    nfound=0
    nsubjobs = len(jobs(jobid).subjobs)
    
    for i in range(nsubjobs):
        if jobs(jobid).subjobs(i).status=="completed":
            nfound+=1
            for output in jobs(jobid).subjobs(i).outputfiles:
                if output.namePattern=="output.root":
                    f.write("%s\n"%(output.lfn))
                    break
        else:
            f.write("NOFILE\n")
    
    f.close()
    
    print "job %d: found %d of %d" % (jobid,nfound,nsubjobs)
