#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:06:15 2019

@author: nvthaomy
"""
import os

def LoadFF(watermodel,pdbMix,libs,ffs):
    """write tleap input file to load forcefield and generate .parm7 and .crd"""
    top = pdbMix[:pdbMix.index('pdb')] + 'parm7'
    crd = pdbMix[:pdbMix.index('pdb')] + 'crd'
    s = ""
    for ff in ffs:
        s += 'source leaprc.{}\n'.format(ff)
    s += 'source leaprc.water.{}\n'.format(watermodel)
    for lib in libs:
        s += 'loadOFF {}\n'.format(lib)       
    s += """\nx=loadpdb {pdb}
addions x Na+ 0
addions x Cl- 0
setbox x vdw 1
saveamberparm x {top} {crd}
savepdb x {pdb}
quit""".format(pdb=pdbMix, top=top, crd=crd)
    file = open('loadFF.in','w')
    file.write(s)
#    os.system('tleap -s -f loadFF.in > loadFF.out')
    return 'loadFF.in',top,crd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", required=True, help="mixture pdb")
    parser.add_argument("-w","--watermodel", required = True,
                        help="Water model for simulation (opc,tip3p,spce), default = opc")
    parser.add_argument("-l", nargs='+', 
			help="paths to tleap libraries for all monomers in mixture")
    parser.add_argument("-ff", nargs='+', help = 'amber forcefields: gaff2, ff99')
    args = parser.parse_args() 
    pdbMix = args.pdb
    watermodel = args.watermodel
    libs = args.l
    ffs = args.ff
    LoadFF(watermodel, pdbMix, libs, ffs)
