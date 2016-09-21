#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np
import datetime as dt
import re

#import geodetic as gd
#import gpsTime as gpsT
import sys
sys.path.append("../rinex/")
sys.path.append("../util/")

import Observation as rnxO
import sinex as snx
import argparse

#================================================================================
def dt2hours(str):
    #print("DURATION:", duration) 

    return(str)
    
#================================================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='checkAusposSolutions',description='Compare the observation time span of the submitted RINEX files with the observation span in the SINEX file')

    parser.add_argument('-s',dest='snxfile',nargs='+',default='',help="SINEX file")
    parser.add_argument('-f', '--file', dest='rnxFiles',nargs='+',default='',help="Submitted RINEX files")
    parser.add_argument('-d',dest='duration',default=False,action='store_true',help="Check the duration of the observatiosn in the SINEX file")
    parser.add_argument('-t',dest='threshold',default=3,help="Check the duration of the observatiosn in the SINEX file")

    args = parser.parse_args()

#================================================================================


    if args.rnxFiles:
        rctr = 1
        for rnxfile in args.rnxFiles:
            obs = rnxO.parseRinexObsFile(rnxfile)
            start = obs['epochs'][0]['time']
            end   = obs['epochs'][-1]['time'] 
            duration = end - start
            #duration = dt2hours(duration)
            if (duration.seconds + duration.days*3600*24) < 3600 * args.threshold :
                print("ERROR ***** ",rnxfile,duration)
            else:
                print(rnxfile,duration)

    if args.snxfile:
        #sinex_data = []
        skipcova = True
        for sf in args.snxfile:
            sinex_data = (snx.readSINEX(sf,skipcova)[0])

            for stn in sinex_data['stations']:
                start = sinex_data[stn]['start_epochs'][0]
                end   = sinex_data[stn]['end_epochs'][0]
                duration = end - start
                if (duration.seconds + duration.days*3600*24) < 3600 * args.threshold :
                    print("ERROR ****** ",sf,stn,duration)
                else:
                    print(sf,stn,duration)


