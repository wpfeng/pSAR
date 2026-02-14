#!/usr/bin/env python
#
#
##########################################
import sys
import os
import glob
import xml.dom.minidom
import zipfile
import pS1
#########################################
if len(sys.argv) < 2:
   #
   helpstr = \
   '''
      %s <mode, 0 or 1> -dir [current folder in default]
      ++++++++++++++++++++++++++++++++++++++++++
      To move S1 TOPS zip files into separate folders based on their own track informtion
      
      Only 0 or 1 is required as input.
      0 means to make a soft link of data into its track folder... 
      1 means to move the data into its track folder

      Provided by Wanpeng Feng, @CCRS,NRCan, 2015
   '''
   print(helpstr % sys.argv[0])
   sys.exit(-1)
# 
#########################################
# kmlstr_s,kmlstr_e,polystr_s,polystr_e = pDATA.kml_poly()
#
model = int(sys.argv[1])
#
###########################################
#
datadir = os.getcwd()
for i,key in enumerate(sys.argv):
    if key == '-dir':
       datadir = sys.argv[i+1]
#
#
zfiles = glob.glob(datadir+'/S1*.zip')
#
#
#
strs = []
for czip in zfiles:
    #
  if zipfile.is_zipfile(czip):
     #
     cxmlinfo = pS1.s1zip2manifestString(czip)
  else:
     #
     print("%s is broken." % czip)
     #
     cxmlinfo = None
  xmlbasename= os.path.basename(czip).split('.')[0]
  if cxmlinfo is not None:
    DOMTree    = xml.dom.minidom.parseString(cxmlinfo)
    collection = DOMTree.documentElement
    #
    # time information
    timeinformation   = 'XXXXX'
    tags       = collection.getElementsByTagName("safe:acquisitionPeriod")
    if len(tags) > 0:
       conts = tags[0].getElementsByTagName("safe:startTime")
       if len(conts) > 0:
          timeinformation = ("%s" % conts[0].childNodes[0].data)
    # safe:platform
    platform   = 'Sentinel-1'
    tags       = collection.getElementsByTagName("safe:platform")
    if len(tags) > 0:
       conts = tags[0].getElementsByTagName("safe:number")
       if len(conts) > 0:
          platform = ("%s%s" % (platform,conts[0].childNodes[0].data))
    #
    # RelativeTracknumber
    relativetrack   = '000'
    tags       = collection.getElementsByTagName("safe:orbitReference")
    if len(tags) > 0:
       conts = tags[0].getElementsByTagName("safe:relativeOrbitNumber")
       if len(conts) > 0:
          relativetrack = ("%s" % conts[0].childNodes[0].data)
    #
    descri = ("%-25s : %s\n %-25s : %s\n %-25s : %s" % ("AcquisitionTime",timeinformation,\
                                                        "Platform",       platform,\
                                                        "RelativeOrbitNumber", relativetrack))
    #
    print(" %s : %s" % (czip,relativetrack))
       
    
    # 
    #
    tags       = collection.getElementsByTagName("safe:footPrint")
    #
    if len(tags) > 0:
       conts = tags[0].getElementsByTagName("gml:coordinates")[0]
    footPrintData = conts.childNodes[0].data.split()
    spoly = []
    counter=0 
    #
    if not os.path.exists('T%s' % relativetrack):
       os.makedirs('T%s' % relativetrack)
    #
    # czip = czip.split('.zip')[0]+'.zip'
    cfile = os.path.abspath(czip)
    #
    if model == 0:
      if os.path.exists(czip):
         #
         if not os.path.exists('T%s' % relativetrack+'/'+os.path.basename(cfile)):
            os.symlink(cfile,'T%s' % relativetrack+'/'+os.path.basename(cfile))
    else:
      if not os.path.exists('T%s' % relativetrack+'/'+os.path.basename(cfile)):
          os.rename(cfile,'T%s' % relativetrack+'/'+os.path.basename(cfile))
    #
