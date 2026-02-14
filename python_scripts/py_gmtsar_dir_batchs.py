#!/usr/bin/env python

#
#
import pSAR
import os
import sys
import numpy as np
import glob
#

if (len(sys.argv)<2):
    #
    helpstr = \
            '''
            %s <roi,e.g. minlon,maxlon,minlat,maxlat> <s1zipsdir> [-update]
                          [-rlook 8]
                          [-azloook 2]
                          [-season 0]
                          [-srtm]
            ++++++++++++++++++++++++++++++++++++++++++++++++++++
            To conduct a batch job to process all S1*.zip files for InSAR

            Note that <s1zips> can have S1 TOPS SAR data from different tracks.
            The script will sort data first into local zips/T*
            
            -update If given, "-update" will be sent to pSAR_gmtsar_s1.py as well.
            -rlook  look number in range, 8 in default
            -azlook look number in azimuth, 2 in default
            Written by wafeng, @SYSU, Guangzhou, 2025/09/10
            Updated by wafeng, @SYSU, Guangzhou, 2025/09/17

            '''
    print(helpstr % sys.argv[0])
    sys.exit(-1)
    #
#
pSAR.util.log(sys.argv)
#
# 
# added by wafeng, @SYSU, Guanghzou, 2025/09/17
#
rlook = 8
azlook = 2
season = 0
srtm=False
#
s_update = ""
for i, ckey in enumerate(sys.argv):
    #
    if ckey == '-srtm':
        srtm = True
    if ckey == '-season':
        season = int(sys.argv[i+1])
    if ckey == '-update':
        s_update = "-update"
    if ckey == '-rlook':
        rlook = int(sys.argv[i+1])
        #
    if ckey == '-azlook':
        azlook = int(sys.argv[i+1])
#     
if True:
    #
    rois = sys.argv[1]
    indir= sys.argv[2]
    #
    indir = os.path.abspath(indir)
     
    topdir = os.getcwd()
    #
    if not os.path.exists('zips'):
       os.mkdir('zips')
    #
    os.system('ln -s %s/S*.zip zips/. ' % indir)
    os.chdir('zips/')
    os.system('pSAR_s1sorting.py 0')
    os.chdir(topdir)
    #
    #
    indir= os.getcwd()+'/zips'
    #
    #
    dirs = glob.glob(indir+'/T*/')
    #
    topdir=os.getcwd()
    #
    #
    for cdir in dirs:
        #
        track = os.path.basename(os.path.dirname(cdir))
        #
        #
        if not os.path.exists(track):
            #
            os.makedirs(track)
            #
        os.chdir(track)
        #
        srtm_str = ''
        if srtm:
            srtm_str = '-dem_source srtm'
            #
        #
        goStr = 'pSAR_gmtsar_s1.py -roi %s -s1dir %s -roidir %s -grd %s -rlook %d -azlook %d -season %d %s' % (rois,cdir,track,s_update,rlook,azlook,season,srtm_str)
        print(goStr)
        os.system(goStr)
        #
        os.chdir(topdir)
        #

        
