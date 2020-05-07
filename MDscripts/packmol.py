#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 11:43:53 2019

@author: nvthaomy
"""
import sys
import numpy
import os 

def PackMol(nMols, pdbs, x, y, z, name): #(f,N,np,nw,watermodel,pdbList,x,y,z,ns,s, saltPdbDicts = {'Na+':'na.pdb', 'Cl-':'cl.pdb'}):
    """nMols: list of number of molecules for each species
    pdbs: list of pdbs of each species"""
    #convert to angstrom and subtrac 1 angstrom
    x,y,z = ((x-.2)*10,(y-.2)*10,(z-.2)*10)
    pdbMix = name + '.pdb'
    mixFile = name + '.inp'
    s = """tolerance 2.0
filetype pdb
output {}\n""".format(name+'.pdb')
    
    for i, nMol in enumerate(nMols):
        if nMol > 0:
            s += """structure {pdb}
\tnumber {n}
\tresnumbers 2
\tinside box 2 2 2 {x} {y} {z}
\tend structure\n""".format(n = nMol, pdb=pdbs[i], x=x, y=y, z=z)

    file = open(mixFile,'w')
    file.write(s)

    #writing job file to pack molecules
#    os.system('packmol < {}'.format(mixFile))      
    return pdbMix        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", nargs='+', type=int, required=True,
                        help="list of number of all molecule types")
    parser.add_argument("-pdb",nargs='+',required=True,
                        help="list of pdb files of all molecule types, in the same order as -n")
    parser.add_argument("-xyz",nargs ="+",type = float,required=True,
                        help="periodic box size in nm, dimensions used in packmol will be ~2A smaller")
    parser.add_argument("-o", default = 'mixture', help = 'name of mixture') 
    args = parser.parse_args()
    nMols = args.n
    pdbs = args.pdb
    x,y,z = (args.xyz[0], args.xyz[1], args.xyz[2])
    name = args.o
    PackMol(nMols, pdbs, x, y, z, name) 
    
    
    
