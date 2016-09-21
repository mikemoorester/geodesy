#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import re
import numpy as np

import datetime as dt
import argparse

import gpsTime as gt

def readSINEX(sinex_file,skip_cova=False):
    """
    Read the SINEX file into a data structure:

    sinex,cova = readSINEX(filename)

    """
    sinex = {}
    cova = []
    sinex['stations'] = []
    last_index = 0
    #============================================
    header_RGX       = re.compile('^\*')
    siteID_START_RGX = re.compile('^\+SITE\/ID')
    siteID_END_RGX   = re.compile('^\-SITE\/ID')
    siteID_Flag      = 0
    #
    receiver_START_RGX = re.compile('^\+SITE\/RECEIVER')
    receiver_END_RGX   = re.compile('^\-SITE\/RECEIVER')
    receiver_Flag      = 0
    #
    antenna_START_RGX = re.compile('^\+SITE\/ANTENNA')
    antenna_END_RGX   = re.compile('^\-SITE\/ANTENNA')
    antenna_Flag      = 0
    #
    solution_epochs_START_RGX = re.compile('^\+SOLUTION\/EPOCHS')
    solution_epochs_END_RGX   = re.compile('^\-SOLUTION\/EPOCHS')
    solution_epochs_Flag      = 0
    #
    solution_estimate_START_RGX = re.compile('^\+SOLUTION\/ESTIMATE')
    solution_estimate_END_RGX   = re.compile('^\-SOLUTION\/ESTIMATE')
    solution_estimate_Flag      = 0
    #
    solution_cova_START_RGX = re.compile('^\+SOLUTION/MATRIX_ESTIMATE L COVA')
    solution_cova_END_RGX   = re.compile('^\-SOLUTION/MATRIX_ESTIMATE L COVA')
    solution_cova_Flag      = 0
    #============================================
    with open(sinex_file,'r') as snxFH:
        for line in snxFH:
            if header_RGX.search(line):
                continue
            elif siteID_START_RGX.search(line):
                siteID_Flag = 1
                continue
            elif siteID_END_RGX.search(line):
                siteID_Flag = 0
            elif siteID_Flag == 1:
                # check if this is a header
                code = line[1:5].upper()
                sinex['stations'].append(code)
                sinex[code] = {}
                sinex[code]['mean_epochs'] = []
                sinex[code]['start_epochs'] = []
                sinex[code]['end_epochs'] = []
                sinex[code]['domes'] = line[9:19]
                sinex[code]['description'] = line[21:44].lstrip()
                lon = float(line[44:48])
                mm  = float(line[48:50])/60.
                ss  = float(line[51:56])/3600.
                sinex[code]['lon'] = lon + mm + ss

                lat = float(line[56:60])
                mm  = float(line[60:63])/60.
                ss  = float(line[63:68])/3600.
                sinex[code]['lat'] = lat + mm + ss
                sinex[code]['h'] = float(line[68:76])
            elif receiver_START_RGX.search(line):
                receiver_Flag = 1
            elif receiver_END_RGX.search(line):
                receiver_Flag = 0
            elif receiver_Flag == 1 :
                code = line[1:5].upper()
                sinex[code]['syy']  = int(line[16:18])
                sinex[code]['sdoy'] = int(line[19:22])
                sinex[code]['ssec'] = int(line[23:28])
                sinex[code]['eyy']  = int(line[29:31])
                sinex[code]['edoy'] = int(line[32:35])
                sinex[code]['esec'] = int(line[36:41])
                sinex[code]['receiver'] = line[42:62].rstrip()
                #print(code,syy,sdoy,ssec,eyy,edoy,esec,rec)
            elif antenna_START_RGX.search(line):
                antenna_Flag = 1
            elif antenna_END_RGX.search(line):
                antenna_Flag = 0
            elif antenna_Flag == 1 :
                code = line[1:5].upper()
                sinex[code]['syy']  = int(line[16:18])
                sinex[code]['sdoy'] = int(line[19:22])
                sinex[code]['ssec'] = int(line[23:28])
                sinex[code]['eyy']  = int(line[29:31])
                sinex[code]['edoy'] = int(line[32:35])
                sinex[code]['esec'] = int(line[36:41])
                sinex[code]['antenna'] = line[42:62].rstrip()
                #print(code,sinex[code]['antenna'])
            elif solution_epochs_START_RGX.search(line):
                solution_epochs_Flag = 1
            elif solution_epochs_END_RGX.search(line):
                solution_epochs_Flag = 0
            elif solution_epochs_Flag == 1 :
                code = line[1:5].upper()
                sinex[code]['point'] = line[6:8].rstrip().lstrip()
                sinex[code]['solution'] = int(line[9:13])
                #sol = sinex[code]['solution']
                # get the start epoch for this solution
                yy  = int(line[16:18])
                doy = int(line[19:22])
                sec = int(line[23:28])
                #print("start:",code,sol,yy,doy,sec)
                dto = gt.yds2dt(yy,doy,sec)
                sinex[code]['start_epochs'].append(dto)

                # get the final epoch for this solution
                yy  = int(line[29:31])
                doy = int(line[32:35])
                sec = int(line[36:41])
                #print("end:",code,sol,yy,doy,sec)
                dto = gt.yds2dt(yy,doy,sec)
                sinex[code]['end_epochs'].append(dto)

                # get the mean epoch for this solution
                yy  = int(line[42:44])
                doy = int(line[45:48])
                sec = int(line[49:54])
                #print("mean:",code,sol,yy,doy,sec)
                dto = gt.yds2dt(yy,doy,sec)
                sinex[code]['mean_epochs'].append(dto)
            elif solution_estimate_START_RGX.search(line):
                solution_estimate_Flag = 1
            elif solution_estimate_END_RGX.search(line):
                solution_estimate_Flag = 0
            elif solution_estimate_Flag == 1 :
                code = line[14:18].upper()
                etype = line[7:13].rstrip().lstrip()
                index = int(line[1:6]) - 1
                #sol  = int(line[])
                if etype[0:3] == 'STA' or etype[0:3] == 'VEL':
                    #unit = line[40:44].rstrip()
                    estimate = float(line[47:68])
                    stdev    = float(line[69:80])
                    param    = etype[0:4].upper() 
                   
                    # check to see if this data structure needs to be created..
                    if param not in sinex[code]:
                        sinex[code][param] = []
                        sinex[code][param+'_index'] = []
                        sinex[code]['s'+param] = []
                        
                    sinex[code][param].append(estimate)
                    sinex[code]['s'+param].append(stdev)
                    sinex[code][param+'_index'].append(index)
                    last_index = index
            elif solution_cova_START_RGX.search(line):
                if skip_cova :
                    return sinex,cova
                solution_cova_Flag = 1
                #print("Will need a matrix which is ",last_index+1,"by",last_index+1)
                cova = np.zeros((last_index+1,last_index+1))
            elif solution_cova_END_RGX.search(line):
                solution_cova_Flag = 0
            elif solution_cova_Flag == 1 :
                index1 = int(line[1:6]) - 1 
                index2 = int(line[7:12]) - 1 
                value1 = float(line[13:34])
                #cova[index1-1,index2-1] = value1
                cova[index1,index2] = value1

                if len(line) > 55:
                    value2 = float(line[35:56])
                    #cova[index1-1,index2]   = value2
                    cova[index1,index2+1]   = value2
                if len(line) > 58:
                    value3 = float(line[57:78])
                    #cova[index1-1,index2+1] = value3
                    cova[index1,index2+2] = value3

    return sinex,cova 

def sinex2pos(sinex):
    """

    convert sinex data structure to pbo time series format

    pos = sinex2pos(SINEX) 

    """
    pos = {}
    for snx in sinex:
        for key in snx:
            if key not in pos:
                pos[key] = []
            pos[key].append(snx[key])


    key = 'ALIC'
    with open(key+".aus.final_frame.pos","w") as fpos :
        print("PBO Station Position Time Series. Reference Frame : IGb08",file=fpos)
        print("Format Version: 1.1.0",file=fpos)
        print("4-character ID: ",key,file=fpos)
        print("Station name  : ",key,file=fpos)
        print("First Epoch   : 19980101 115900",file=fpos)
        print("Last Epoch    : 20131228 115900",file=fpos)
        print("Release Data  : 20150225 110056",file=fpos)
        stax = pos[key][0]['STAX'][0]
        stay = pos[key][0]['STAY'][0]
        staz = pos[key][0]['STAZ'][0]
        print("XYZ Reference position : {:>15.5f} {:>14.5f} {:>14.5f} (IGb08)".format(stax,stay,staz),file=fpos)
        lat  = pos[key][0]['lat']
        lon  = pos[key][0]['lon']
        height = pos[key][0]['h']
        #print("NEU Reference position :   -29.0465539119  115.3469773647  241.27517 (IGb08/WGS84)",file=fpos)
        print("NEU Reference position :  {:>15.10f} {:>15.10f} {:>10.5f} (IGb08/WGS84)".format(lat,lon,height),file=fpos)
        print("Start Field Description",file=fpos)
        print("YYYYMMDD      Year, month, day for the given position epoch",file=fpos)
        print("HHMMSS        Hour, minute, second for the given position epoch",file=fpos)
        print("JJJJJ.JJJJJ   Modified Julian day for the given position epoch",file=fpos)
        print("X             X coordinate, Specified Reference Frame, meters",file=fpos)
        print("Y             Y coordinate, Specified Reference Frame, meters",file=fpos)
        print("Z             Z coordinate, Specified Reference Frame, meters",file=fpos)
        print("Sx            Standard deviation of the X position, meters",file=fpos)
        print("Sy            Standard deviation of the Y position, meters",file=fpos)
        print("Sz            Standard deviation of the Z position, meters",file=fpos)
        print("Rxy           Correlation of the X and Y position",file=fpos)
        print("Rxz           Correlation of the X and Z position",file=fpos)
        print("Ryz           Correlation of the Y and Z position",file=fpos)
        print("Nlat          North latitude, WGS-84 ellipsoid, decimal degrees",file=fpos)
        print("Elong         East longitude, WGS-84 ellipsoid, decimal degrees",file=fpos)
        print("Height (Up)   Height relative to WGS-84 ellipsoid, m",file=fpos)
        print("dN            Difference in North component from NEU reference position, meters",file=fpos)
        print("dE            Difference in East component from NEU reference position, meters",file=fpos)
        print("du            Difference in vertical component from NEU reference position, meters",file=fpos)
        print("Sn            Standard deviation of dN, meters",file=fpos)
        print("Se            Standard deviation of dE, meters",file=fpos)
        print("Su            Standard deviation of dU, meters",file=fpos)
        print("Rne           Correlation of dN and dE",file=fpos)
        print("Rnu           Correlation of dN and dU",file=fpos)
        print("Reu           Correlation of dEand dU",file=fpos)
        print("Soln          \"rapid\", \"final\", \"suppl/suppf\", \"campd\", or \"repro\" corresponding to products  generated with rapid or final orbit products, in supplemental processing, campaign data processing or reprocessing",file=fpos)
        print("End Field Description",file=fpos)
        print("*YYYYMMDD HHMMSS JJJJJ.JJJJ         X             Y             Z            Sx        Sy       Sz     Rxy   Rxz    Ryz            NLat         Elong         Height         dN        dE        dU         Sn       Se       Su      Rne    Rnu    Reu  Soln",file=fpos)

        print("POS:ALIC",pos['ALIC'][0])
        for i in range(0,np.size(pos['ALIC'])):
            print("I:",i)
            YYYYMMDD = pos[key][i]['mean_epochs'][0].strftime("%Y%m%d")
            HHMMSS   = pos[key][i]['mean_epochs'][0].strftime("%H%M%S")
            JJJ      = pos[key][i]['mean_epochs'][0].strftime("%j")+".50000"
            X        = pos[key][i]['STAX'][0]
            Y        = pos[key][i]['STAY'][0]
            Z        = pos[key][i]['STAZ'][0]
            sX       = pos[key][i]['sSTAX'][0]
            sY       = pos[key][i]['sSTAY'][0]
            sZ       = pos[key][i]['sSTAZ'][0]
            Rxy      = 0.0
            Rxz      = 0.0
            Ryz      = 0.0
            Nlat     = pos[key][i]['lat'] # should be dN
            Elong    = pos[key][i]['lon'] # should be dE
            h        = pos[key][i]['h']   # should be dh
            dN       = 0.0
            dE       = 0.0
            dH       = 0.0
            sN       = 0.0
            sE       = 0.0 
            sU       = 0.0
            Rne      = 0.0
            Rnu      = 0.0
            Reu      = 0.0
            Soln     = 'final'

            print(YYYYMMDD,HHMMSS,JJJ,"{:>14.5f} {:>14.5f} {:>14.5f} {:>8.5f} {:>8.5f} {:>8.5f} {:>6.3f} {:>6.3f} {:>6.3f} {:>18.10f} {:>18.10f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>6.3f} {:>6.3f} {:>6.3f} {:s}".format(X,Y,Z,sX,sY,sZ,Rxy,Rxz,Ryz,Nlat,Elong,h,dN,dE,dH,sN,sE,sU,Rne,Rnu,Reu,Soln),file=fpos)
            #print(YYYYMMDD,HHMMSS,JJJ,file=fpos)

    return pos

def removeSolution(remove_SOLN,snx,cova):
    '''
    removeSolution remove a selected list of stations from the SINEX file
    
    This will remove all solutions associated with a station. For example
    ALIC will remove the veolcity solution ALIC_0 and ALIC_1, etc..
    
    However if you specify ALIC_0,only this velocity will be removed
    
    Currently this will only remove the Covariance information, and does not
    rigorously remove the information from the SINEX file.
    
    '''
    del_ind = []
    
    for code in remove_SOLN:
                
        # find all of the codes which match
        for i in range(0,np.size(snx[code]['VELX_index'])):
            # work out what indicies need to be removed
            del_ind.append(snx[code]['STAX_index'][i])
            del_ind.append(snx[code]['STAY_index'][i])
            del_ind.append(snx[code]['STAZ_index'][i])
            del_ind.append(snx[code]['VELX_index'][i])
            del_ind.append(snx[code]['VELY_index'][i])
            del_ind.append(snx[code]['VELZ_index'][i])
          
    del_ind = np.sort(del_ind)
    del_ind[:] = del_ind[::-1]
    #print(del_ind)

    #print("COVA",np.shape(cova))    
    for d in del_ind:
        cova = np.delete(cova,d,0)
        cova = np.delete(cova,d,1)    
    #print("COVA",np.shape(cova))    

    return snx,cova
#================================================

if __name__ == "__main__":

    #=======
    parser = argparse.ArgumentParser(description='Calculate the plate model from a set of GPS velocities')

    parser.add_argument('-f',dest='snxfile',nargs='+',default='',help="SINEX file")
    parser.add_argument('--skip_cova',dest='skip_cova',default=False,action='store_true',help="Skip covariance")

    parser.add_argument('--pos',dest='snx2pos',default=False,action='store_true',help="Convert SINEX to POS format")
    #parser.add_argument('-y',dest='year',default=2010,type=int)
    parser.add_argument('-d',dest='duration',default=False,action='store_true',help="Check the duration of the observatiosn in the SINEX file")

    args = parser.parse_args()
    #======
    sinex = []
    if args.snxfile:
        ctr = 1
        for sf in args.snxfile:
            sinex.append(readSINEX(sf,args.skip_cova)[0])
            ctr = ctr + 1

        if args.snx2pos:
            sinex2pos(sinex)

    if args.duration:        
        for s in sinex:
            for stn in s['stations']:
                #print("stn:",stn
                start = s[stn]['start_epochs'][0]
                end   = s[stn]['end_epochs'][0]
                duration = end - start
                #print("For station",stn,start,end)
                print("For station",stn,duration)
                #val = str(dt.timedelta(seconds=int(duration)))
                #print("For station",stn," the duration is:",val)

    #print(np.size(sinex))

