f = open("inputs.txt")
lines = f.readlines()
f.close()

fbash = open("inputs.sh","w")
fbash.write("#!/bin/bash\n")

fbash.write("vim -E -s main.tex <<-EOF\n")
for line in lines:
    line = line.split("{")[1].split("}")[0].replace("/","\\/").split(".",1)
    base = line[0]
    if len(line) > 1:
        ext = line[1]
    else:
        ext = "tex"
    fbash.write("\t:g/\\\\input{%s/r %s.%s\n" %(base,base,ext))
    fbash.write("\t:g/\\\\input{%s/s/^.*$//\n" %(base))

fbash.write("\t:wq\n")
fbash.write("EOF\n")
