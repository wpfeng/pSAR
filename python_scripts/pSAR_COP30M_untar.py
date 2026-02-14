#!/usr/bin/env python
#
#
import glob
import os
import sys
#
#
if __name__ == '__main__':
    #
    #
    tars = glob.glob('Coper*.tar')
    #
    #
    for ctar in tars:
        #
        outdir = os.path.basename(ctar).split('.tar')[0]
        print(outdir)
        if not os.path.exists(outdir):
            #
            goStr = 'tar xvf %s' % ctar
            os.system(goStr)

