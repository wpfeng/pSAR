#!/usr/bin/env python
#
#
#
############################################
import sys, numpy as np
import pSAR
import pDATA
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
############################################
def usage():
    log =\
    ''' INFO    : @NRCan, part of Automated processing system for RCM, by Feng, Wanpeng
        Usage   : %s <roi_data_file> 
                     -dtype [float in default] 
                     -poly [point or polygon, point in default]
                     -scale [10, in default]  [-nowrap] 
                     -outloc <pySAR_locs.list in default> -pos [none in default]
                     -is180

     +Updated by FWP, @SYSU, Guangzhou, 2020/06/10
        since this version, GRD format file is supported...

     (c) 2015-2018 Canada Center for Mapping and Earth Observations, NRCan, Ottawa Canada.\n
    '''
    print(log % sys.argv[0])
    #
    #
############################################
#
if len(sys.argv)<2:
   usage()
   sys.exit(-1)
#
roi_file = sys.argv[1]
#
dtype = 'float32'
mode  = "POINT"
outloc= "pySAR_locs.list"
pos   = 999999. 
minv  = -5
maxv  = 5
wrap  = True
scale = 50
is180 = False
###########################################
counter = 0
for ckey in sys.argv:
   counter += 1
   if ckey == '-is180':
       is180 = True
   if ckey == "-dtype":
      dtype = sys.argv[counter]
   if ckey == "-outloc":
      outloc = sys.argv[counter]
   if ckey == "-poly":
      mode = "POLYGON"
   if ckey == "-scale":
      scale = int(float(sys.argv[counter]))
   if ckey == "-pos":
      pos = np.float32(sys.argv[counter])
   if ckey == "-nowrap":
      wrap = False
###########################################
if 'grd' in roi_file:
    # 
    print(roi_file)
    info,ext,data = pDATA.grd_read(roi_file)
else:
    rsc_file  = roi_file + '.rsc'
    info,ext  = pSAR.roipac.rsc_read(rsc_file)
    data      = pSAR.roipac.roi_read(roi_file,dtype=dtype)
#
fig = plt.figure()
ax  = fig.add_subplot(111)
if wrap:
  ax.imshow(pSAR.roipac.wrap(data[::scale,::scale],minv=minv,maxv=maxv),cmap='jet',interpolation='none',extent=ext)
else:
  ax.imshow(data,interpolation='none',extent=ext)
#
ispoly=False
if mode.upper() == "POINT":
  # fig = plt.subplot(111)
  output = pSAR.ui.pickup_points(fig,outloc);
  plt.show()
  print(dir(output))
  #
else:
  
  output = pSAR.ui.draw_rect();
  plt.show()
  outp = np.array([output.poly_x,output.poly_y])
  #print(outp)
  #
  
  minlon = np.min(outp[0,:])
  maxlon = np.max(outp[0,:])
  minlat = np.min(outp[1,:])
  maxlat = np.max(outp[1,:])
  if pos != 999999.:
     nx = np.round((maxlon-minlon)/pos)
     maxlon = minlon + (nx-1) * pos
     ny = np.round((maxlat-minlat)/pos)
     minlat = maxlat - (ny-1) * pos
     #
  if is180:
      if minlon>180:
          minlon = minlon - 360
      if maxlon > 180:
          maxlon = maxlon - 360 
  #
  print(" Poly: %f,%f,%f,%f" % (minlon,maxlon,minlat,maxlat))
  print(" GMT_FORMAT: -R%f/%f/%f/%f" % (minlon,maxlon,minlat,maxlat))
  ispoly=True
#
poly = pSAR.roipac.ext2polygon([minlon,maxlon,minlat,maxlat])
np.savetxt(outloc,poly,fmt='%f %f')
#
plt.close(fig)
#
#############################################
# return mean value of target region region if ispoly
if ispoly:
   #
   #
   ext = [minlon,maxlon,minlat,maxlat]
   if 'grd' in roi_file:
       print(' Not supported yet to calculate mean value for a grd file... :(')
       sys.exit(-1)
   #
   sd,sx,sy = pSAR.roipac.roi_subgrid(roi_file,ext)
   #
   sd[np.isnan(sd)] = 0
   subdata = sd[sd !=0 ]
   meanv   = np.mean(sd[np.nonzero(sd)])
   sdlos = sd * float(info['WAVELENGTH']) * -1. / (4. * np.pi) * 100.
   sdlos = sdlos[sdlos != 0.]
   stdv  = np.std(sdlos)
   print(" Mean value in the region: %f and standard deviation: %f" %  (meanv,stdv))
#
