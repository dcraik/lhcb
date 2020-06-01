skip = [0,1,2]
#for jobid in [0,1,2]:
#for jobid in [4,5,6]:
for jobid in range(4):
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

    print("job %d: found %d of %d" % (jobid,nfound,nsubjobs))
