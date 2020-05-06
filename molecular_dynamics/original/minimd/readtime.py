import os
import sys
import re

""" First Simulation and then Analysis will be printed with space """

argc = len(sys.argv)



if argc !=2 :
    print("Please Specify the file to read time")
    sys.exit(1)

path = sys.argv[1]
if os.path.exists(path) == False:
    print("File doen't exists")
    sys.exit(1)


fd = open(path,"r")
lines = fd.readlines()

simt = 0.0
anat = 0.0

for l in lines:
    x = re.search("^Analysis Time",l)
    if x:
        """ current line is simulation"""
        l = l.split()
        anat = float(l[2])
    else:
        """current line is analysis"""
        x = re.search("^Simulation Time",l)
        if x:
            l = l.split()
            simt = float(l[2])
            
print(str(simt))
