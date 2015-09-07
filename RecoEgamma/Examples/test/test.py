#!/usr/bin/env python
# -*- coding:Utf-8 -*-
import os, sys
import glob
import re
import subprocess
import shutil

if len(sys.argv) > 2:
#     stripScriptFromHtml( sys.argv[1], sys.argv[2] )
     print sys.argv[1], sys.argv[2]
else:
     if len(sys.argv) > 1:
#        stripscripts( sys.argv[1] )
        print sys.argv[1]
        ind = int(sys.argv[1])
     else:
        print "rien"
        
dirname = os.getcwd()
fileList = glob.glob("*.py")

print fileList
print

filename = 'toto_file'

for x in range(0, 10+ind):
    print x
    print "x = %003d" % x
    print filename + '-%003d'%x + '-old.txt'

print "fin"