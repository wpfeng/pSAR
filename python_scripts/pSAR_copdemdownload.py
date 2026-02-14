#!/usr/bin/env python
#
#
import os
import numpy 
import sys
import glob
import pSAR
#
#
if len(sys.argv)<2:
    helpstr = \
            '''
            %s lon0,lon1,lat0,lat1 -output <dem.grd>
               -copdem_dir [/home/wafeng/NAS80T_2/DEM/COPD, in default]
               -mode [30M or 90m, 30M in default]
               -download_only

            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            To download and merge Copernicus DEM data (30 meters or 90 meters) with ESA API...
            
            It seems that Copernicus DEM can be much better than SRTM dataset for InSAR applications...

            Wanpeng Feng, @SYSU, Guangzhou, 2023/11/20
            '''
    print(helpstr % sys.argv[0])
    #
    sys.exit(-1)
    #
#
#
pSAR.util.log(sys.argv)
#
download_only = False
#
roi = [float(croi) for croi in sys.argv[1].split(',')]
output = 'dem.grd'
#
copdem_dir = '/home/wafeng/NAS80T_2/DEM/COPD'
copdem_dir = '/home/wafeng/NAS80T_2/COPDEM/COP30M/'
mode = '30M'
fmt = 'GMT'
#
#
for i,key in enumerate(sys.argv):
    #
    if key == '-fmt':
        fmt = sys.argv[i+1]
    if key == '-output':
        output = sys.argv[i+1]
    if key == '-copdem_dir':
        copdem_dir = sys.argv[i+1]
    if key == "-mode":
        mode = sys.argv[i+1].upper()
        #
    if key == '-download_only':
        download_only = True
    #
#
#
if True:
    #
    #
    d_lon0 = int(roi[0]-1)
    d_lon1 = int(roi[1]+1)
    d_lat0 = int(roi[2]-1)
    d_lat1 = int(roi[3]+1)
    #
    topdir = os.getcwd()
    #
    cdem_dir = srtm_home = os.environ['COPDEM_HOME']
    if not os.path.exists(cdem_dir):
        os.makedirs(cdem_dir)
        #
    #
    os.chdir(cdem_dir)
    #
    #
    goStr = 'pSAR_COP30M_downloader.py -lon0 %d -lon1 %d -lat0 %d -lat1 %d' % (d_lon0,d_lon1,d_lat0,d_lat1)
    print(goStr)
    os.system(goStr)
    #
    os.chdir(topdir)
    #
    #
    if not download_only:
      # pSAR_COP30M_tar2grd.py <in_dir> -roi <lon0,lon1,lat0,lat1> -outgrd
      goStr = 'pSAR_COP30M_tar2grd.py $PWD -roi %s -outgrd %s' % ('%f,%f,%f,%f' % (roi[0],roi[1],roi[2],roi[3]),output)
      print(goStr)
      os.system(goStr)
      #
      #  
      os.chdir(topdir)
      #
      if os.path.exists(cdem_dir+'/%s' % output):
          #
          goStr = 'ln -s %s .' % os.path.abspath(cdem_dir+'/%s' % output)
          print(goStr)
          os.system(goStr)
    #
    #
    if fmt.upper() == 'ROI':
       # print(" FWP", os.getcwd())
       goStr = 'pSAR_imgformat.py %s %s -of EHdr' % (output,topdir+'/'+output.split('.')[0]+'.img')
       print(goStr)
       os.system(goStr)
    #
    if fmt.upper() == 'GMT' or fmt.upper() == 'NETCDF':
       goStr = 'pSAR_imgformat.py %s %s -of netCDF' % (output+'.tif',topdir+'/'+output+'.grd')
       print(goStr)
       os.system(goStr)
       #
    #
    os.chdir(topdir)
