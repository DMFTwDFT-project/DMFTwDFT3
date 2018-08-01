#!/usr/bin/env python
import sys,re,os
import optparse
from scipy import *

import numpy
if numpy.__version__ == '1.0.1':
    loadtxt = io.read_array
    def savetxt(filename, data):
        io.write_array(filename, data, precision=16)  

if __name__=='__main__':
    """ Takes several self-energy files and produces an average over these self-energy files
    """
    usage = """usage: %prog [ options ] argumens

    The script takes several self-energy files and produces an average self-energy

    arguments  -- all input self-energies
    option -o  -- output self-energy file
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--osig", dest="osig", default='sig.inpx', help="filename of the output self-energy file. Default: 'sig.inp'")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-n", "--nexecute", dest="nexecute", action="store_true", default=False,  help="Do not execute the comments")
    parser.add_option("-s", dest="stdev", action='store_true', default=False, help="Computes also the standard deviation - the error of self-energy")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    print 'files to average over:', args
    print 'output: ', options.osig
    ws_oo=[]
    wEdc=[]
    wdata=[]
    for f in args:
        dat = open(f,'r').readlines()
        s_oo = None
        Edc = None
        data=[]
        for line in dat:
            m = re.search('#(.*)',line)
            if m is not None and not options.nexecute:
                exec(m.group().strip())
            else:
                data.append( map(float, line.split() ) )

        if s_oo is not None: ws_oo.append(s_oo)
        if Edc is not None: wEdc.append(Edc)
        wdata.append(data)    
    
    fout = open(options.osig, 'w')
    
    if len(ws_oo):
        ws_oo = array(ws_oo)
        as_oo=[]
        for i in range(shape(ws_oo)[1]):  as_oo.append( sum(ws_oo[:,i])/len(ws_oo) )
        print 's_oo=', as_oo
        print >> fout, '# s_oo=', as_oo

    if len(wEdc):
        wEdc = array(wEdc)
        aEdc=[]
        for i in range(shape(wEdc)[1]):  aEdc.append( sum(wEdc[:,i])/len(wEdc) )
        print 'Edc=', aEdc
        print >> fout, '# Edc=', aEdc

    wdata = array(wdata)
    wres = zeros(shape(wdata)[1:], dtype=float)
    for i in range(len(wdata)): wres[:,:] += wdata[i,:,:]
    wres *= 1./len(wdata)
    
    if options.stdev:
        sw = shape(wdata)
        wstd = zeros((sw[1],sw[2]-1), dtype=float) # no frequency
        for i in range(len(wdata)): wstd[:,:] += wdata[i,:,1:]**2
        wstd *= 1./len(wdata)
        wstd[:,:] = sqrt(wstd[:,:] - wres[:,1:]**2)
        
        wres = hstack( (wres, wstd) )
        
    savetxt(fout,wres)
    
    
    
