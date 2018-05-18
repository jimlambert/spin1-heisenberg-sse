#-------------------------------------------------------------------------------
# Author: James Lambert
# Date: February 22, 2018
#
# Generate table of system commands for serial farming.
#-------------------------------------------------------------------------------

import numpy as np
import sys
import argparse

def main(argv):
    parser = argparse.ArgumentParser(description='Produce run table')
    parser.add_argument('--tablename', type=str, dest='table',
                        default='runtable.dat',
                        help='name of run table')
    parser.add_argument('--sizerng', nargs='+', type=int, dest='sizerng',
                        help='range of sizes')
    parser.add_argument('--equil', type=int, dest='es', default=10000,
                        help='equilibration steps')
    parser.add_argument('--simul', type=int, dest='ss', default=10000,
                        help='simulation steps')
    parser.add_argument('--temprng', nargs='+', type=float, dest='temprng',
                        help='range of temperatures')
    parser.add_argument('--dfldrng', nargs='+', type=float, dest='dfldrng',
                        help='range of D-fields')
    parser.add_argument('--epsilon', type=float, dest='eps', default=0.0,
                        help='epsilon value')
    parser.add_argument('--pbc', type=bool, dest='pbc', default=True,
                        help='boundary conditions')
    parser.add_argument('--odir', type=str, dest='odir', default="./",
                        help='default directory for data output')
    
    args = parser.parse_args()
    aSizes = np.arange(args.sizerng[0], args.sizerng[1]+args.sizerng[2], 
                       args.sizerng[2])
    aDflds = np.linspace(args.dfldrng[0], args.dfldrng[1],
                         args.dfldrng[2])
    aTemps = np.linspace(args.temprng[0], args.temprng[1],
                         args.temprng[2])
    rtable = open(args.table, 'w')
    for s in aSizes:
        for t in aTemps:
            for d in aDflds:
                cmdstr = "./HESSE_S1 "
                cmdstr += "--size "  + str(s) + " "
                cmdstr += "--equil " + str(args.es) + " "
                cmdstr += "--simul " + str(args.ss) + " "
                cmdstr += "--temp "  + str(t) + " "
                cmdstr += "--dfld "  + str(d) + " "
                cmdstr += "--eps "   + str(args.eps) + " "
                ofilename = args.odir + "ssedata"
                if(args.pbc):
                    cmdstr += "--pbc " + str(1) + " "
                    ofilename += "-pbc"
                else: 
                    cmdstr += "--pbc " + str(0) + " "
                    ofilename += "-obc"
                ofilename += "-" + str(s) + "-" + str(t)
                if d >= 0.0:
                    ofilename += "-" + str(0)
                else:
                    ofilename += "-" + str(1)
                ofilename += "-" + str(abs(d)) + ".dat"
                cmdstr += "--of " + ofilename
                rtable.write(cmdstr + '\n') 
    return 0

if __name__=='__main__':
    main(sys.argv[1:])
