import sys
import re

args = sys.argv

for arg in args[1:]:
  f1 = open("tables2012L0HadronTOS_S20realET/tables_S20/"+arg,"r")
  f2 = open("tables2012L0HadronTOS_S20realET/"+arg+".dat","w")
  for line in f1:
    a = re.findall("\-{0,1}[0-9]+\.{0,1}[0-9]*e*-{0,1}[0-9]*|-",line)
    if len(a) == 0:
      continue
    b = [a[0]]
    b += a[2:]
    f2.write("\t".join(b))
    f2.write("\n")
  f1.close()
  f2.close()
