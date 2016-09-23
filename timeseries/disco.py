from __future__ import print_function
import numpy as np
import re 

import sys
import os

import ftplib
import filecmp

from datetime import date
from shutil import copyfile

sys.path.insert(0,os.path.abspath('../util/'))

import gpsTime as gt

#==============================================================================
def backup_file(srcdir,backupdir,filename):
    d = date.today()
    src = srcdir+filename
    backup=backupdir+d.isoformat()+"_"+filename

    print("Backing up:",src,"to:",backup)
    copyfile(src,backup)

    return backup

#==============================================================================
def download_IGN_disco():
    '''

    '''

    # Download the latest copy from IGN
    ftp = ftplib.FTP('igs-rf.ign.fr','anonymous','acc@igs.org')
    ftp.cwd("/pub/IGS14P")
    fp = open("./downloads/soln_IGS14P.snx",'w')
    ftp.retrlines('RETR '+"soln_IGS14P.snx",lambda s, w = fp.write: w(s + '\n') )
    ftp.quit()
    fp.close()

    if filecmp.cmp("./downloads/soln_IGS14P.snx","./data/soln_IGS14P.snx"):
        #print("File compare is True files are the same - no need to do an update")
        # no need to do an update
        return 0
    else:
        #print("File compare is false")
        # first backup the orginal file
        backup_file("./data/","./data/backup/","soln_IGS14P.snx")

        # Now copy the file to the data archive 
        copyfile("./downloads/soln_IGS14P.snx","./data/soln_IGS14P.snx")

        return 1
    
#==============================================================================
def parse_disco_file(discofile) :
    '''
    disco = parse_disco_file(discofile) 
    
    read in discontinuity file in SINEX format 

    input:    
    path to disco file => soln_IGS14P.snx

    output:
    

    '''
    discoStartFlag = 0
    discoEndFlag   = 0

    discoStart_RGX = re.compile('^\+SOLUTION\/DISCONTINUITY')
    discoEnd_RGX   = re.compile('^\-SOLUTION\/DISCONTINUITY')

    disco = {}

    with open(discofile, 'r') as DFID:
        for line in DFID:
            if discoStart_RGX.search(line):
                discoStartFlag  = 1
                continue
            elif discoEnd_RGX.search(line) :
                discoEndFlag = 1
                continue
            elif (discoStartFlag == 1) & (discoEndFlag == 0) :
                station  = line[1:5]
                if station not in disco:
                    disco[station] = []

                d = {}
                d['marker']    = str(line[6:11]).rstrip().lstrip()
                d['dnum']     = int(line[11:13].rstrip().lstrip())
                d['flag']     = str(line[14])
                d['startyy']  = int(line[16:18])
                d['startdoy'] = int(line[19:22])
                d['startsec'] = int(line[23:28])
                d['endyy']    = int(line[29:31])
                d['enddoy']   = int(line[32:35])
                d['endsec']   = int(line[36:41])
                d['dtype']    = line[42]
                d['reason']   = line[44:].rstrip().lstrip()
                disco[station].append(d)
                #print("Station:",station,domes,dnum,startyy,startdoy,startsec,endyy,enddoy,endsec,dtype,reason)
                
    return disco 

#==============================================================================
def print_disco_file(disco,fname):
    """

    """
    WFID = open(fname, 'w')

    print('%=SNX 2.02 GA 16:187:29323 IGN 94:002:00000 15:046:00000 P     0 0 S',file=WFID)
    print('*-------------------------------------------------------------------------------',file=WFID)
    print('+SOLUTION/DISCONTINUITY',file=WFID)

    stations = disco.keys()
    stations.sort()
    for station in stations:
        for d in disco[station]:
            #print("station:<",station,">")
            #print("marker :<",d['marker'],">")
            #print("dnum   :<",d['dnum'],">")
            #print("flag   :<",d['flag'],">")
            line = ' {:4s}  {:<3s} {:2d} {:1s}'.format(station,d['marker'],d['dnum'],d['flag'])
            line = line+' {:02d}:{:03d}:{:05d}'.format(d['startyy'],d['startdoy'],d['startsec'])
            line = line+' {:02d}:{:03d}:{:05d}'.format(d['endyy'],d['enddoy'],d['endsec'])
            line = line+' {:s}'.format(d['dtype'])
            if 'reason' in d:
                line = line+' {:s}'.format(d['reason'])
            else:       
                line = line+' -'

            print(line, file=WFID)
    print('-SOLUTION/DISCONTINUITY',file=WFID)

    return 1

#==============================================================================
def str2num(val):
    """

    """
    dict = { 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 
             'S':19, 'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26 }
    return int(dict[val])

#==============================================================================
def num2str(val):
    """

    """
    
    if val < 10:
        return val

    dict = { '10':'J', '11':'K', '12':'L', '13':'M', '14':'N', '15':'O', '16':'P', '17':'Q', '18':'R', 
             '19':'S', '20':'T', '21':'U', '22':'V', '23':'W', '24':'X', '25':'Y', '26':'Z' }
    return dict[str(val)]


def print_tsview_file(disco,fname):
    """

    """

    WFID = open(fname, 'w')

    for station in disco:
        if len(disco[station]) < 3:
            continue

        for i in range(1,len(disco[station])):
            # check to see if this is a velocity discontinuity, if 
            # so move onto the next one
            d = disco[station][i]
            #e = disco[station][i+1]
            if d['dtype'] == 'V': # or e['dtype'] == 'V':
                continue

            line = ' rename {:4s}    {:4s}'.format(station,station)

            if d['dnum'] < 10:
                line = line+'_{:1d}PS '.format(d['dnum'])
            else:
                line = line+'_{:1s}PS '.format(num2str(d['dnum']))

            SYYYY = gt.yy2yyyy(d['startyy']) 

            yyyy,MM,DD =  gt.yyyydoy2cal(SYYYY,d['startdoy'])
            hh = 0
            mm = 0

            line = line+'{:4d} {:02d} {:02d} {:02d} {:02d} '.format(SYYYY,MM,DD,hh,mm)

            EYYYY = gt.yy2yyyy(d['endyy']) 
            if EYYYY < SYYYY:
                EYYYY = 2100

            yyyy,MM,DD =  gt.yyyydoy2cal(EYYYY,d['enddoy'])
            line = line+'{:4d} {:02d} {:02d} {:02d} {:02d}'.format(EYYYY,MM,DD,hh,mm)

            print(line,file=WFID)

    return 1

#==============================================================================

def read_renamefile(renamefile):
    """


    """

    disco = {}

    XPS_RGX = re.compile('_XPS ')

    #==========================================================================
    # Read in the file and parse it into a data structure 
    #==========================================================================
    with open(renamefile, 'r') as RFID:
        for line in RFID:
            if XPS_RGX.search(line):
                continue

            vals = line.split()

            station = vals[1][0:4]
            if station not in disco:
                disco[station] = []

            try:
                disco_num = int(vals[2][5])
            except ValueError:
                disco_num = str2num(vals[2][5])

            start_yyyy = vals[3]
            start_mon  = vals[4]
            start_dd   = vals[5]
            start_hh   = vals[6]
            start_mm   = vals[7]

            print("station:",station,disco_num,start_yyyy,start_mon,start_dd,start_hh,start_mm)
            
            d = {}
            d['marker']   = 'A' 
            d['dnum']     = disco_num
            d['flag']     = 'P'
                
            if disco_num == 1:
                d['startyy']  = 00
                d['startdoy'] = 000
                d['startsec'] = 0

                d['endyy']    = gt.yyyy2yy( int(start_yyyy) )
                d['enddoy']   = gt.ymd2yyyyddd(int(start_yyyy),int(start_mon),int(start_dd))[1]
            else: 
                print("disco_num:<",disco_num,">",np.size(disco[station]))
                print(disco[station])
                dp = disco[station][disco_num-2]
                d['startyy']  = dp['endyy']
                d['startdoy'] = dp['enddoy'] 
                d['startsec'] = 0

                d['endyy']  = gt.yyyy2yy( int(start_yyyy) )
                d['enddoy'] = gt.ymd2yyyyddd(int(start_yyyy),int(start_mon),int(start_dd))[1]
                
            d['endsec']   = 0
            d['dtype']    = 'P'

            disco[station].append(d)

    #==========================================================================
    # Add in the null velocity discontinuity needed by CATREF
    #==========================================================================
    for station in disco:
        d = {}
        d['marker']   = 'A' 
        d['dnum']     = 1
        d['flag']     = 'P'
        d['startyy']  = 00
        d['startdoy'] = 000
        d['startsec'] = 0
        d['endyy']    = 00
        d['enddoy']   = 000
        d['endsec']   = 0
        d['dtype']     = 'V'

        disco[station].append(d)

    return disco

#==============================================================================

def read_eqfile(eqfile):
    """
    read_eqfile(filename)

    read in a eqfile (GLOBK) formatted disconintuity file into the 'disco'
    data structure

    """

    disco = {}

    rename_RGX = re.compile('^ rename')

    with open(eqfile, 'r') as EQFID:
        for line in EQFID:
            if rename_RGX.search(line):
                station  = line[8:12]
                if station not in disco:
                    disco[station] = []
                #print("Station:",station)
                disco_num  = line[21]
                try:
                    disco_num = int(line[21])
                except ValueError:
                    disco_num = str2num(line[21])

                start_yyyy = line[25:29]
                start_mon  = line[30:32]
                start_dd   = line[33:35]
                start_hh   = line[36:38]
                start_mm   = line[39:41]
                #print("Start:",start_yyyy,start_mon,start_dd,start_hh,start_mm)

                end_yyyy   = line[42:46]
                end_mon    = line[47:49]
                end_dd     = line[50:52]     
                end_hh     = line[53:55]
                end_mm     = line[56:58]
                #print("End:",end_yyyy,end_mon,end_dd,end_hh,end_mm)
                
                if disco_num == 2:
                    d = {}
                    d['marker']   = 'A' 
                    d['dnum']     = 1
                    d['flag']     = 'P'
                    d['startyy']  = 00
                    d['startdoy'] = 000
                    d['startsec'] = 0
                    d['endyy']    = gt.yyyy2yy( int(start_yyyy) )
                    d['enddoy']   = gt.ymd2yyyyddd(int(start_yyyy),int(start_mon),int(start_dd))[1]
                    d['endsec']   = 0
                    d['dtype']     = 'P'
                    disco[station].append(d)
                
                    
                d = {}
                d['marker']   = 'A' 
                d['dnum']     = disco_num
                d['flag']     = 'P'
                d['startyy']  = gt.yyyy2yy( int(start_yyyy) )
                d['startdoy'] = gt.ymd2yyyyddd(int(start_yyyy),int(start_mon),int(start_dd))[1]
                d['startsec'] = 0

                if int(end_yyyy) == 2100:
                    d['endyy']  = 0
                    d['enddoy'] = 0
                    d['endsec'] = 0 
                    d['dtype']   = 'P'
                else:
                    d['endyy']  = gt.yyyy2yy( int(end_yyyy) )
                    d['enddoy'] = gt.ymd2yyyyddd(int(end_yyyy),int(end_mon),int(end_dd))[1]
                    d['endsec'] = 0 
                    d['dtype']   = 'P'
                
                disco[station].append(d)

    # Add in the null velocity discontinuity needed by CATREF
    for station in disco:
        d = {}
        d['marker']   = 'A' 
        d['dnum']     = 1
        d['flag']     = 'P'
        d['startyy']  = 00
        d['startdoy'] = 000
        d['startsec'] = 0
        d['endyy']    = 00
        d['enddoy']   = 000
        d['endsec']   = 0
        d['dtype']     = 'V'

        disco[station].append(d)

    return disco


#==============================================================================

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="read in discontinuity files in SINEX format, or the GLOBK .eq format, so that they can be merged or updated")

    parser.add_argument("-f", dest="filename", nargs='+',
                                help="input files (can be more than one) -f <file1> <file2>", metavar="FILE")

    parser.add_argument("-u", dest="update",default=False,action='store_true',
                                help="Update Discontinuities, -f <file1> <file2>, update <file1> based on unique discontinuties found in <file2>")

    parser.add_argument("-d", dest="downloadIGN",default=False,action='store_true',
                                help="Download IGN discontinuity file")

    parser.add_argument("-t", dest="convertToEq",default=False,action='store_true',
                                help="Convert to TSVIEW format")

    parser.add_argument("-e", dest="eqfile", nargs='+',
                                help="read in globk format discontinuity", metavar="FILE")

    parser.add_argument("-r", dest="renames", nargs='+',
                                help="read in tsview renames discontinuity file", metavar="FILE")

    parser.add_argument("-o", dest="outfile",
                                help="output file to print merged or updated disco file", metavar="FILE")
    
    parser.add_argument("-m", dest="merge",default=False,action='store_true',
                                help="Merge Discontinuities - only use this option when your creating a brand new file. A Merge will only add new stations into the original discontinuity file")

    args = parser.parse_args()

    #==========================================================================
    discos = []

    if args.downloadIGN:
        if download_IGN_disco():
            print("Will attempt to do an update:")
            args.update = True
            args.filename = []
            args.filename.append('./data/GA_discontinuity_file.snx')
            args.filename.append('./data/soln_IGS14P.snx')
            args.outfile = './data/GA_discontinuity_file.snx'
        else:
            print("No change in IGN discontinuity file")
            raise SystemExit

    if args.eqfile:
        for eqfile in args.eqfile:
            disco = read_eqfile(eqfile)
            print_disco_file(disco,eqfile+'.snx')

    #==========================================================================
    # Option to parse the tsview 'saved' discontinuities into
    # SINEX format
    #==========================================================================
    if args.renames:
        for renamefile in args.renames:
            disco = read_renamefile(renamefile)
            print_disco_file(disco,renamefile+'.snx')

    if args.filename:
        for discofile in args.filename:
            disco = parse_disco_file(discofile) 
            discos.append(disco)

    if args.merge:
        ndiscos = np.size(discos)
        if ndiscos <= 1:
            print("Need to have more than one disconinuity file if you want to merge them together")
            raise SystemExit

        print("Merging the following files:",args.filename[1:],"into:",args.filename[0])
        #======================================================
        # Take the first file as the orginal disconinuity file
        #======================================================
        odisco = discos[0]
        for n in range(1,ndiscos):
            mdisco = discos[n]
            #print("Looking at N:",n)
            for station in mdisco:
                if station in odisco:
                    md = np.size(mdisco[station])
                    od = np.size(mdisco[station])
                    if md > od :
                        print("******************************************************************************************")
                        print("             WARNING ")
                        print("This station already exists in the original file but the new file has more discontinuities")
                        print("Original:",odisco[station])
                        print("Merge   :",mdisco[station])
                        print(" Consider doing an update (-u) ...")
                        print("******************************************************************************************")

                else:
                    #print("Adding station",station)
                    odisco[station] = mdisco[station]

        # print the file out
        if args.outfile:
            fname = args.outfile
        else:
            fname = 'merged.snx'

        print_disco_file(odisco,fname)

    #==========================================================================
    #
    # UPDATE:
    #  - parse both files into a data structure
    #  - check the station disconitnuities of each file
    #  - if there are more discontinutities for a station in the second file:
    #    keep the new discontintuities, discard the original
    # 
    #  - copy the original file and backit up, dump the new updated 
    #    discontinuities into the original filename
    #  
    #==========================================================================
    if args.update:
        ndiscos = np.size(discos)
        if ndiscos <= 1:
            print("Need to have more than one disconinuity file if we want to update")
            raise SystemExit

        print("Will update the following file",args.filename[0],"with",args.filename[1:])
        # Take the first file as the orginal disconinuity file
        odisco = discos[0]
        for n in range(1,ndiscos):
            mdisco = discos[n]
            for station in mdisco:
                if station in odisco:
                    sm = np.size(mdisco[station])
                    so = np.size(odisco[station])

                    # If there are more disconintuites in the new file
                    # for a stations, store these into the official
                    # file so that they are updated
                    if sm > so :
                        print("Updating the discontinuties for station:",station)
                        # odisco[station] = mdisco[station]
                        # Try and preserver the 'reason' field
                        # of the original file
                        # .EQ files don't provide this, and update
                        # using this file will 'delete' previous 
                        # comments unless the following works:
                        for i in range(0,so):
                            md = mdisco[station][i]
                            #print("MD:",md)
                            od = odisco[station][i]
                            if ('reason' in md) and (md['reason'] != '-'):
                                continue
                            if 'reason' in od:
                                #print("Adding reason!",od['reason'])
                                md['reason'] = od['reason']  
                        odisco[station] = mdisco[station]

                else:
                    odisco[station] = mdisco[station]

        dirname  = os.path.split(args.filename[0])[0]
        filename = os.path.split(args.filename[0])[1]
        if dirname == '':
            srcdir = './'
        else:
            srcdir = dirname+'/'

        print("before applying update we will back up the file:",args.filename[0])
        backup_file(srcdir,srcdir,filename)

        #print_disco_file(odisco,filename)
        print_disco_file(odisco,args.filename[0])

    if args.convertToEq:
        # print the file out
        if args.outfile:
            fname = args.outfile
        else:
            fname = 'tsview.eq'

        print_tsview_file(disco,fname)
