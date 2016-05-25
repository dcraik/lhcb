f = open("inputs.txt")
lines = f.readlines()
f.close()

fbash = open("inputs.sh","w")
fbash.write("#!/bin/bash\n")

fbash.write("vim -E -s main.tex <<-EOF\n")
for line in lines:
    line = line.split("{")[1].split("}")[0]
    if line.endswith('.tex'):
        line = line[:-4]
    fbash.write("\t:g/\\\\input{%s/r %s.tex\n" %(line,line))
    fbash.write("\t:g/\\\\input{%s/s/^.*$//\n" %(line))

fbash.write("\t:wq\n")
fbash.write("EOF\n")
