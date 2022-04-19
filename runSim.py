#!/usr/bin/python

import shutil
import os

logf=file("log.txt","w")
step0=5
for i in range(0,250,5):
    outfile="gps%dMeV.mac"%i
    simfile="gps%dMeV.root"%i
    shutil.copy("testgps.mac",outfile)
    os.system("sed -i \"s/250/%d/g\" %s"%(i,outfile))
    print "../g4Radem %s %s"%(simfile,outfile)
    logf.write("./g4Radem %s %s\n"%(simfile,outfile))
    os.system("./g4Radem %s %s"%(simfile,outfile))




