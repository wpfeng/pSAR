#!/usr/bin/env python
#
#
import os
import sys
import glob
import pSAR
#
#
#
#
# https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-90-DGED__2022_1/Copernicus_DSM_30_N30_00_E060_00.tar
#
if len(sys.argv)<2:
    helpstr = \
            '''
            %s [run] -lon0 [80 in default] -lon1 [120, in default] 
                     -lat0 [20 in default] -lat1 [55, in default]
                     -mode [30m in default]
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            To download Copernicus DEM (30m) for use in SAR interferometry
            
            -mode  option for 30 meter or 90 meter
                   30 meter in default, as 30m
            URL:
            90 meter:
            https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-90-DGED__2022_1/Copernicus_DSM_30_N30_00_E060_00.tar 
            
            30 meter:
            https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-30-DGED__2022_1/Copernicus_DSM_10_N30_00_E060_00.tar

            Example:
               %s -lon0 80 -lon1 82 -lat0 23 -lat1 24

            '''
    print(helpstr % (sys.argv[0],os.path.basename(sys.argv[0])))
    #
    sys.exit(-1)
    #
#
pSAR.util.log(sys.argv)
#   
mode = '30m'
lon0 = 80
lon1 = 120
lat0 = 20
lat1 = 55
#
for i,key in enumerate(sys.argv):
    #
    if key == '-mode':
        mode = sys.argv[i+1]
    if key == '-lon0':
        lon0 = int(sys.argv[i+1])
    if key == '-lon1':
        lon1 = int(sys.argv[i+1])
    if key == '-lat0':
        lat0 = int(sys.argv[i+1])
    if key == '-lat1':
        lat1 = int(sys.argv[i+1])
#
url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-30-DGED__2022_1/'
dem_prefix = 'Copernicus_DSM_10_'
#
if mode.upper()=='90M':
    url = 'https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-90-DGED__2022_1/'
    #
    dem_prefix = 'Copernicus_DSM_30_'
#
#
if True:
    #
    s_lat = 'N'
    s_lon = 'E'
    #
    for clon in range(lon0,lon1,1):
        for clat in range(lat0,lat1,1):
            #
            s_lat = 'N'
            s_lon = 'E'
            #
            if clat <0:
                s_lat = 'S'
            if clon <0:
                s_lon = 'W'
            #
            outname = '%s%s%02d_00_%s%03d_00.tar' % (dem_prefix,s_lat,abs(clat),s_lon,abs(clon))
            #
            maxlon = clon+1
            maxlat = clat+1
            #
            print(clon,clat,outname)
            goStr = "wget -c %s/%s" % (url,outname) 
            #
            #os.stat(outname).st_size
            #
            if not os.path.exists(outname):
                #
                print(goStr)
                os.system(goStr)
                #
