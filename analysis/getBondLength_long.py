#!/usr/bin/env python
import sys
sys.path.append('/home/mnguyen/bin/sim_git')
import pickleTraj
import numpy as np
def getMap(traj, top, nameMap):
    """map AA traj to CG traj
    assuming map one AA residues to one CG bead
    input: nameMap = {AAres:CGatomType}
    outputs: [CGatomTypes], [[AtomIds in bead 1],[AtomIds in bead 2], ...]"""    
    import mdtraj as md
    traj = md.load(traj, top = top)
    top = traj.topology
    Mass1List = []
    AAatomId = []
    AAres = []
    CGatomTypes = []
    for res in top.residues:
        AAres.append(res.name)
        CGatomTypes.append(nameMap[res.name])
        #get atom indices of all atoms in this residue
        Atoms1 = []
        Mass1 = []
        for atom in res.atoms:
            Atoms1.append(atom.index)
            Mass1.append(atom.element.mass)
        AAatomId.append(Atoms1)
        Mass1List.append(Mass1)
    Mol = []
    for bond in top.bonds:
        atom1,atom2 = bond[0],bond[1]
        res1,res2 = atom1.residue, atom2.residue
        if not res1 in Mol:
            Mol.append(res1)
        if res1 != res2 and not res2 in Mol:
            Mol.append(res2)
            
    return AAatomId, CGatomTypes, AAres, Mass1List 

def convertTraj(traj, top, lengthScale = 1., stride = 1, outTrajExt = '.lammpstrj.gz'):
    """convert trajectory to a specified format and scale box with lengthScale"""
    import mdtraj as md
    import os,shutil
    #to make traj in the current dir
    cwd = os.getcwd()
    outTraj = traj.split('/')[-1]
    outTraj = '.'.join(outTraj.split('.')[:-1]) + outTrajExt
    traj = md.load(traj, top = top, stride = stride)
    if lengthScale != 1.:
        print('Scaling positions and box size by 1/{}'.format(lengthScale))
    traj.xyz /= lengthScale
    traj.unitcell_lengths /= lengthScale
    traj.save(outTraj)
    shutil.move(outTraj,os.path.join(cwd,outTraj))
    print('moving traj to {}'.format(os.path.join(cwd,outTraj))) 
    outTraj = os.path.join(cwd,outTraj)
    return outTraj

def mapTraj(traj, top, nameMap, lengthScale, stride=1, outExt = '_mapped.lammpstrj.gz'):
    import sim
    """ CGatomTypes: 1 by n list of CG atom types
        AAatomID: n by x list of indices of AA atoms in CG beads
        lengthScale (nanometer): divide positions and box dimensions by this scale, for conversion between real and dimensionless units"""
    AAatomId, CGatomTypes, AAres, Mass1List = getMap(traj, top, nameMap)
    AtomTypes, counts = np.unique(CGatomTypes, return_counts = True)
    # ===== create mapped object =====
    print("\n ===== Creating Index Mapper =====")
    Map = sim.atommap.PosMap()
    for i, CGatomType in enumerate(CGatomTypes):
        Atoms1 = AAatomId[i]
        Atom2 = i
        Mass1 = Mass1List[i]
        this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2, Mass1=Mass1)
        Map += [this_Map]
    
    print("\n ===== Converting AA traj to lammpstrj format  =====")
    traj = convertTraj(traj, top, lengthScale = lengthScale, stride = stride, outTrajExt = '.lammpstrj')
    outTraj = traj.split('.lammpstrj')[0] + outExt    
    outPdb = traj.split('.lammpstrj')[0] + '_mapped.pdb'
    # ===== read AA traj =====
    print("\n ===== Reading scaled AA Traj =====")
    Trj = pickleTraj(traj)
    BoxL = Trj.FrameData['BoxL']
    
    # ===== write out new mapped traj =====
    print("\n ===== Mapping and Writing Trajectory =====")
    MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = CGatomTypes, BoxL = BoxL)
    #MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)
    print('mapped traj name {}'.format(outTraj)) 
    # ===== convert to lammps =====
    print("\n ===== Converting to LAMMPS =====")
    #sim.traj.base.Convert(MappedTrj, sim.traj.LammpsWrite, FileName = outTraj, Verbose = True)    
    sim.traj.base.Convert(MappedTrj, sim.traj.PdbWrite, FileName = outPdb, Verbose = True)
    print("\nCG Atom\tcount:")
    for i, atom in enumerate(AtomTypes):
        print('%s\t%i'%(atom,counts[i]))
    return  CGatomTypes, AAatomId, MappedTrj, BoxL, outPdb 
    
def GetBond(CGtrajPdbF, NP, DOP, a1, a2, nbins):
    import matplotlib.pyplot as plt
    import matplotlib
    showPlots = True
    try:
      os.environ["DISPLAY"] #Detects if display is available
    except KeyError:
      showPlots = False
      matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
    matplotlib.rc('font', size=7)
    matplotlib.rc('axes', titlesize=7)

    traj = md.load(CGtrajPdbF)
    top = traj.topology
    print('=== computing bond length ===')
    pairIds = []
    for i,atom in enumerate(top.atoms):
        if i+1 > NP*DOP:
            break
        if i == 0:
            an0 = atom.name
            ai0 = atom.index
        else:
            an1 = atom.name
            ai1 = atom.index
            if (i) % DOP > 0.: #this atom is not head 
                if (an0,an1) == (a1,a2) or (an0,an1) == (a2,a1): #check if this is the right bond
                    pairIds.append([ai0,ai1])
            an0 = an1
            ai0 = ai1
    #get bond lengths
    bs = md.compute_distances(traj,pairIds)
    bs = np.ravel(bs)
    hist,bins = np.histogram(bs, bins=nbins, density=True)
    binmid = 0.5*(bins[1:]+bins[0:-1])
    print('Max bond length {:3.2f} nm'.format(bs.max()))
    data = np.vstack([binmid,hist]).T
    np.savetxt('bond{}_{}.dat'.format(a1,a2), data, header= 'BondLength Probability')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(binmid, hist, marker=None,ls='-',lw=1, c = 'k')
    plt.xlim(0,1)
    plt.ylim(0)
    plt.xlabel('Bond Length')
    plt.ylabel('Histogram')
    title = 'bond {}_{} Histogram'.format(a1,a2)
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
    return 

if __name__ == "__main__":
    import argparse as ap
    import mdtraj as md
    import numpy as np
    import os, sys, re
    import argparse as ap

    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('np', type = float, help = 'number of chains')
    parser.add_argument('DOP', type = float, help = 'DOP')
    parser.add_argument('-atoms', required = True, type=str, nargs=2, help='CG atom names to calculate bond length')
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-nbins', type=float, default=1000)
    args = parser.parse_args()

    traj = sys.argv[1]
    top = sys.argv[2]
    NP = float(sys.argv[3])
    DOP = float(sys.argv[4])
    stride = args.stride
    a1 = args.atoms[0]
    a2 = args.atoms[1]
    nbins = args.nbins
    nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH',
               'ATP':'A', 'AHP':'A', 'AP': 'A', 'ATD': 'A-', 'AHD': 'A-', 'AD': 'A-',
               'NTP':'B+', 'NHP':'B+', 'NP': 'B+', 'NTD': 'B', 'NHD': 'B', 'ND': 'B'}
    lengthScale = 1.

    print('=== mapping trajectory ===')
    CGatomTypes, AAatomId, CGtraj, BoxL, CGtrajPdbF = mapTraj(traj,top,nameMap, lengthScale, stride = stride, outExt ='.lammpstrj.gz')
    GetBond(CGtrajPdbF, NP, DOP, a1, a2, nbins)
 
