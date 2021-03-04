
skip=[0,1,2,3,4,13,14,17,18,19,20,29,101]
#skip=[]
#for jobid in range(42):
#for jobid in range(35,41):
#for jobid in range(42,58):
#for jobid in range(60,92):
for jobid in range(68,84):
#for jobid in range(84,92):
#for jobid in range(95,118):
#for jobid in range(131,157):
#for jobid in range(2):#gangadir-jets
    if jobid in skip: continue

    f=open("lfns/lfns%d.log"%jobid, "w")#gangadir-dcraik
    f2=open("lumi/lumi%d.log"%jobid, "w")#gangadir-dcraik
    f3=open("local/local%d.log"%jobid, "w")#gangadir-dcraik
    #f=open("lfns/lfnsX%d.log"%jobid, "w")#gangadir-jets
    #f2=open("lumi/lumiX%d.log"%jobid, "w")#gangadir-jets

    nfound=0
    nsubjobs = len(jobs(jobid).subjobs)

    for i in range(nsubjobs):
        if jobs(jobid).subjobs(i).status=="completed":
            nfound+=1
            for output in jobs(jobid).subjobs(i).outputfiles:
                if output.namePattern=="output.root":
                    f.write("%s\n"%(output.lfn))
                    f3.write("%s/output.root\n"%(output.localDir))
                if output.namePattern=="LumiTuple.root":
                    f2.write("%s\n"%(output.lfn))
        else:
            f.write("NOFILE\n")
            f2.write("NOFILE\n")
            f3.write("NOFILE\n")

    f.close()
    f2.close()
    f3.close()

    print("job %d (%s): found %d of %d" % (jobid,jobs(jobid).name,nfound,nsubjobs))
