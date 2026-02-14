#!/usr/bin/env python
#
#
# by Feng, Wanpeng, @NRCan, 2016-03-09
# to visually see the offsets of coregistration
#
#########################################
#
from __future__ import print_function
import sys
import os
import pSAR
import glob
import numpy as np
import matplotlib.mlab as mlab 
import matplotlib.pyplot as plt
import pDATA
#
#######################################
if len(sys.argv) < 3:
   #
   helpstr = \
   '''
    %s <inphase> -mask [minlon,maxlon,minlat,maxlat] -iters [1 in default]
            -refroi [] -npara [3 in default] -output [none, in default] 
            -theilsen -plot  -dem [demdata] -scale [1,in default]
            -outdir [None in default]
            -zshift [False in default]
    ******************************************************************
     To correct orbital errors with a polynominal fitting
     V0.1, Wanpeng Feng, @NRCan, 2016-05-17
     ---
     -sim      only error SIM will be saved...
     -iters    An iteration was proposed since 2017-08-21, very useful
     -dem      a dem file to estimate dem-correlated signals
     -theilsen Applied since this version. If given, median will be used 
               to exclude outliers from the reference points.
     Updated by Wanpeng Feng, @NRCan, 2017-03-01
     ---
     -npara  3,4,6 can be allowed.
     added by FWP, 2017-04-07
     Updated by Wanpeng Feng,@CCRS/NRCan, 2017-08-21

     ---
     <in_phs> in grd can be supported since this version.
     Updated by FWP, @SYSU, Guangzhou, 2019/07/07

   '''
   print(helpstr % sys.argv[0])
   #
   sys.exit(-1)
#######################################
#
pSAR.util.log(sys.argv)
#
#
cphs   = sys.argv[1]
roi    = [0,0,0,0] 
scale  = 5
demfile = None
offs   = 0.3
plot   = False
model  = 0
maxfev = 1000
ramp   = False
output = None
npara  = 3
refroi = 'NULL'
theilsen = False
iters    = 1
outdir = None
isim = False
zshift = False
#
counter = 0
for ckeyword in sys.argv:
   counter += 1
   if ckeyword == '-zshift':
      zshift = True
   if ckeyword == '-sim':
      isim = True
   if ckeyword == '-outdir':
      outdir = sys.argv[counter]
      #
   if ckeyword == '-iters':
      iters = int(sys.argv[counter])
      #
   if ckeyword == "-theilsen":
      theilsen = True
      #
   if ckeyword == "-dem":
      demfile = sys.argv[counter]
   if ckeyword == "-refroi":
      refroi = sys.argv[counter]
   if ckeyword == '-npara':
      npara = float(sys.argv[counter])
   if ckeyword == '-output':
      output = sys.argv[counter]
   if ckeyword == '-ramp':
      ramp = True
   if ckeyword == '-mask':
      tmp = sys.argv[counter].split(',')
      minlon = float(tmp[0])
      maxlon = float(tmp[1])
      minlat = float(tmp[2])
      maxlat = float(tmp[3])
      roi    = [minlon,maxlon,minlat,maxlat]
   if ckeyword == '-scale':
      scale = int(sys.argv[counter])
   if ckeyword == '-maxfev':
      maxfev = int(sys.argv[counter])
   if ckeyword == '-off':
      offs = np.float32(sys.argv[counter])
   if ckeyword == '-plot':
      plot = True
   if ckeyword == '-model':
      model = int(sys.argv[counter])

#######################################
#
mphs = glob.glob(cphs)
for cphs in mphs:
  print(" Working %s now" % cphs)
  #
  fext = os.path.basename(cphs).split('.')[-1]
  #
  # info,ext  = pSAR.roipac.rsc_read(rsc)
  # lonm,latm = pSAR.roipac.rsc_to_llm(rsc)
  #
  toGRD = 0
  #
  if fext.upper() == 'GRD':
     #
     info,ext,allphs = pDATA.grd_read(cphs)
     _,_,_,coord =  pDATA.grd_read_xr(cphs,proj=True)
     #
     _,lons,lats = pDATA.grd_read_xr(cphs)
     #
     #
     lonm,latm = pSAR.roipac.info_to_llm(info)
     rsc = cphs.split('.grd')[0]+'.rsc'
     #
     #
     if not os.path.exists(rsc):
         #
         pSAR.roipac.info_to_rsc(info,rsc)
     #
     goGRD = 1
     #
  else:
     #
     rsc       = cphs + '.rsc'
     info,ext  = pSAR.roipac.rsc_read(rsc)
     lonm,latm = pSAR.roipac.rsc_to_llm(rsc)
     fmt    = pSAR.roipac.roi_to_fmt(cphs)
     allphs = pSAR.roipac.roi_read(cphs,dtype=fmt)
     #
  #
  # print(info)
  # sys.exit(-1)
  #
  if os.path.exists(refroi):
     rfext = os.path.basename(refroi).split('.')[-1]
     if rfext.upper() == "GRD":
        refdata,_,_ = pDATA.grd_read_xr(refroi)
     else:
        refdata = pSAR.roipac.roi_read(refroi,dtype=fmt)
  else:
     refdata = np.copy(allphs)
  #
  refdata[np.isnan(refdata)] = 0.
  #
  if zshift:
     #
     allphs[allphs!=0]   = allphs[allphs!=0]   - float(info['Z_SHIFT'])
     refdata[refdata!=0] = refdata[refdata!=0] - float(info['Z_SHIFT'])
  #
  # 
  maskdata = np.copy(refdata)
  #
  if (roi[0] != 0 and roi[1] !=0):
     #
     maskdata = pSAR.roipac.roi_maskdata(maskdata,rsc,roi)
  #
  if (demfile is not None and os.path.exists(demfile) ):
     demext = demfile.split('.')[-1]
     #
     if demext.upper() == 'GRD':
         _,_,demdata = pDATA.grd_read(demfile)
         #
     else: 
        demdata = pSAR.roipac.roi_read(demfile)
  else:
     demdata = "NULL"
  #
  #
  # Updatedy by Wanpeng Feng, @NRCan, 2017-03-01
  # sbas has been changed to ts since this version.
  # Anyway, sbas will still be kept for a while...
  #
  citer = 0
  cor = np.copy(allphs)
  #
  sim = np.copy(cor)*0
  sim[np.isnan(sim)] = 0
  #
  while citer < iters:
    print(" NO: %d" % citer)
    Dcor,ramp = pSAR.ts.orbitalCor(maskdata,lonm,latm,dem=demdata,\
                     istheilsen=theilsen,scale=scale,npara=npara,demexp=False)
    citer += 1
    sim = sim + ramp
    #
    maskdata[maskdata!=0] = maskdata[maskdata!=0] - ramp[maskdata!=0]
    cor = cor - ramp
  #
  # plt.imshow(sim)
  # plt.show()
  # sim = allphs - cor
  if isim is not True:
     sim[allphs==0] = 0
  #
  cor[allphs==0] = 0
  #
  #
  if output is not None:
     outrsc = output + '.rsc'
  else:
     if outdir is None:
        outdir = os.path.basename(cphs).split('.')[0]+'-ORB'
     output = outdir + '/' + os.path.basename(cphs)
     outrsc = output + '.rsc'
     if not os.path.exists(outdir):
         os.makedirs(outdir)
  #
  # make sure to edit siminfo only and no change occurrs for info
  #
  siminfo = info.copy()
  if isim is not True:
     if goGRD == 0:
       print(" pSAR: output a file to %s" % output)
       pSAR.roipac.info_to_rsc(info,outrsc)
       pSAR.roipac.roi_write(cor,output)
       #
       #
       siminfo['Z_SHIFT'] = 0
       pSAR.roipac.info_to_rsc(siminfo,output+'.sim.rsc')
       pSAR.roipac.roi_write(sim,output+'.sim')
     else:
       pDATA.grd_write_xr(cor,lons,lats,output,coord=coord)
     #
     ######################################
     #
     # fig1 = plt.figure()
     outfig = output+'.jpg'
     pSAR.ui.datashow(allphs,title='raw',r=3,l=1,no=1)
     pSAR.ui.datashow(sim,title='SIM',r=3,l=1,no=2)
     pSAR.ui.datashow(cor,title='cor',r=3,l=1,no=3)
     fig = plt.gcf()
     fig.tight_layout()
     plt.savefig(outfig)
     output = None
  else:
     #
     # plt.imshow(ramp)
     # plt.show()
     #
     #
     if toGRD == 0:
       print(" pSAR: %s is done!" % (output+'.sim'))
       pSAR.roipac.roi_write(sim,output+'.sim')
       siminfo['Z_SHIFT'] = 0
       pSAR.roipac.info_to_rsc(siminfo,output+'.sim.rsc')
     else:
       #
       #
       pDATA.grd_write_xr(sim,lons,lats,output,coord=coord)
       #
