#!/usr/bin/env python2

"""
Read AusCOM/CICE restart file and
reset date to required.
"""

__Author__ = 'Petteri Uotila <petteri.uotila@csiro.au>'
__Version__ = 0.1
__Date__ = '07/04/2011'

import sys
import struct
import getopt

class CICEdump:
    def __init__(self,nx=360,ny=300,ncat=5,\
                 nilyr=4,nslyr=1,ocnmixed=False):
        #AusCOM grid specifics:
        self.nx = nx # no points in x
        self.ny = ny # no points in y
        self.ncat = ncat # no ice categories
        self.nilyr = nilyr # no of layers in ice
        self.nslyr = nslyr # no of layers in snow
        self.ntilyr = ncat*nilyr
        self.ntslyr = ncat*nslyr
        self.bsize = nx*ny
        self.doublesize = 8  
        # if CICE standalone set this True:
        self.ocnmixed = ocnmixed
        nr = 24 if self.ocnmixed else 22
        self.nrecs = self.ncat*4+self.ntilyr+self.ntslyr+nr

    def read_buffer(self,fi):
        """
        read record and separator as strings without unpack
        """
        recsep = fi.read(self.doublesize)
        buf = fi.read(self.doublesize*self.bsize)
        return recsep, buf

    def read_buffers(self,fin='iced.sis2-cnyf2-41.00100101'):
        """
        read buffers and separators to lists
        """
        fi = open(fin,'rb')
        bufs = []; rseps = []
        # read timing information
        self.header = fi.read(24)
        for n in range(self.nrecs):
            recsep,buf = self.read_buffer(fi)
            rseps.append(recsep)
            bufs.append(buf)
        self.footer = fi.read(4)
        fi.close()
        self.rseps = rseps
        self.bufs = bufs

    def write_buffers(self,fon):
        """
         write buffers to fon
        """
        fp = open(fon,'wb')
        fp.write(self.header)
        for n in range(self.nrecs):
            fp.write(self.rseps[n])
            fp.write(self.bufs[n])
        fp.write(self.footer)
        fp.close()

    def print_header(self):
        bint, istep0, time, time_forc = \
             struct.unpack('>iidd',self.header)
        print "istep0=%d" % istep0
        print "time=%f" % time
        print "time_forc=%f" % time_forc

    def change_header(self,istep0=0,time=0.0,time_forc=0.0):
        """ Change binary header
        """
        self.oldheader = self.header
        bint, istep0old, timeold, time_forcold = \
             struct.unpack('>iidd',self.header)
        self.header = struct.pack('>iidd',bint,istep0,time,time_forc)

if __name__=='__main__':
    def usage():
        print """%s:
        Change time of the CICE dumpfile.
        -i input file
        -o output file [if not given shows timestamp of input file]
        -h this help
        -v verbose output
        --istep0 time step to be written [default 0]
        --time to be written [default 0.0]
        --time_forc to be written [default 0.0]
        """ % sys.argv[0]

    def change_time(ocnmixed=False):
        dump = CICEdump(ocnmixed=ocnmixed)
        dump.read_buffers(config["-i"])
        if config.has_key("-v"):
            print "Read %s" % config["-i"]
            print "Time information:"
            dump.print_header()
        dump.change_header(istep0=int(config["--istep0"]),\
                           time=float(config["--time"]),\
                           time_forc=float(config["--time_forc"]))
        if config.has_key("-v"):
            print "Time changed to:"
            dump.print_header()
        dump.write_buffers(config["-o"])

    # default values for command line args
    config = { "-i" : "/short/p66/pju565/OUTPUT/AusCOM1.0/sis3-cnyf2-23/restart/cice/iced.00110101",\
               "--istep0" : "0",\
               "--time" : "0.0",\
               "--time_forc" : "0.0"}
    try:
        if len(sys.argv)>1:
            options, operands = getopt.getopt( sys.argv[1:],'hvi:o:',\
                    ['istep0=','time=','time_forc='])
            config.update( dict(options) )
    except:
        usage(); sys.exit(0)
    if config.has_key("-h"):
        usage(); sys.exit(0)
    if not config.has_key("-o"):
        dump = CICEdump()
        dump.read_buffers(config["-i"])
        dump.print_header()
        sys.exit(0)
    try:
        change_time(ocnmixed=False)
    except:
        if config.has_key("-v"):
            print "This dump seems to have ocean mixed layer data."
        change_time(ocnmixed=True)
    print "Finished"
