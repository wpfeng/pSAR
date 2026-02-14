#!/usr/bin/env python
#
import pSAR
import os
import sys
#
#
if True:
    #
    ingrd = sys.argv[1]
    inrsc = sys.argv[2]#
    #
    # info_from_grd(ingrd,azi=-12,inc=32,wavelength=0.055246,inrsc=None)
    #
    #
    info = pSAR.roipac.info_from_grd(ingrd,inrsc=inrsc)
    #
    outrsc = ingrd.replace('.grd','.rsc')
    #

    pSAR.roipac.info_to_rsc(info,outrsc)
