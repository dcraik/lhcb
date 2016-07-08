f = open("figs.txt")
lines = f.readlines()
f.close()

fbash = open("figs.sh","w")
fbash.write("#!/bin/bash\n")

iFig=0
iSub=0

fbash.write("vim -E -s main.tex <<-EOF\n")
for line in lines:
    if line.find("begin{figure}") is not -1:
        iFig+=1
        iSub=0
    else:
        fig=line.split("includegraphics")[1].split("{")[1].split("}")[0].split(".")[0].split("/")[-1]
        fbash.write("\t:g/includegraphic/s/{\(figs\/\)\=%s\(.pdf\)\=}/{fig%d%c}\n" % (fig, iFig, chr(iSub+ord('a'))))
        iSub+=1

fbash.write("\t:g/graphicspath{/d\n");
fbash.write("\t:wq\n")
fbash.write("EOF\n")

iFig=0
iSub=0

for line in lines:
    if line.find("begin{figure}") is not -1:
        iFig+=1
        iSub=0
    else:
        fig=line.split("includegraphics")[1].split("{")[1].split("}")[0].split(".")[0].split("/")[-1]
        fbash.write("mv figs/%s.pdf fig%d%c.pdf\n" % (fig, iFig, chr(iSub+ord('a'))))
        iSub+=1

f.close()
fbash.close()
