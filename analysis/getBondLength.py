import os, sys, re
import numpy as np
import mdtraj as md
    
def GetBond(traj, top, stride, NP, DOP, r1, r2, nbins):
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

    traj = md.load(traj,top = top, stride=stride)
    top = traj.topology

    AAres = []
    AAatomId = []  
    for res in top.residues:
        if not res.name in ['HOH', 'Na+', 'Cl-', 'WAT']:
            AAres.append(res.name)
            #get atom indices of all atoms in this residue
            Atoms1 = []
            for atom in res.atoms:
                Atoms1.append(atom.index)
            AAatomId.append(Atoms1)
    COMs = []
    for atomId in AAatomId:
        COM = md.compute_center_of_mass(traj.atom_slice(atomId)) #time series of COM for this residue 
        COMs.append(COM)
    
    print('=== computing bond length ===')
    pairIds = []
    bs = []
    for i, res in enumerate(AAres):
        if i+1 > NP*DOP:
            break
        if i == 0:
            rn0 = res
            COM0 = COMs[i]
        else:
            rn1 = res
            COM1 = COMs[i]
            if (i) % DOP > 0.: #this atom is not head 
                if (rn0,rn1) == (r1,r2) or (rn0,rn1) == (r2,r1): #check if this is the right bond
                    #compute bond length between these 2 residues in all frames
                    b = (COM1-COM0)
                    b = np.sqrt(b[:,0]**2 +b[:,1]**2 + b[:,2]**2)
                    bs.append(b)
            rn0 = rn1
            COM0 = COM1

    bs = np.ravel(bs)
    hist,bins = np.histogram(bs, bins=nbins, density=True)
    binmid = 0.5*(bins[1:]+bins[0:-1])
    #normalize histogram
    #hist = hist/float(sum(hist))
    print('Max bond length {:3.2f} nm'.format(bs.max()))
    data = np.vstack([binmid,hist]).T
    np.savetxt('bond{}_{}.dat'.format(r1,r2), data, header= 'BondLength ProbabilityDensity')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(binmid, hist, marker=None,ls='-',lw=1, c = 'k')
    plt.xlim(0.9*bs.min(), 1.1*bs.max())
#    plt.ylim(0)
    plt.xlabel('Bond Length')
    plt.ylabel('Histogram')
    title = 'bond {}_{} Histogram'.format(r1,r2)
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
    return bs
       
if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('np', type = float, help = 'number of chains')
    parser.add_argument('DOP', type = float, help = 'DOP')
    parser.add_argument('-res', required = True, type=str, nargs=2, help='residues to calculate bond length')
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-nbins', type=float, default=100)
    args = parser.parse_args()

    traj = sys.argv[1]
    top = sys.argv[2]
    NP = float(sys.argv[3])
    DOP = float(sys.argv[4])
    stride = args.stride
    r1 = args.res[0]
    r2 = args.res[1]
    nbins = args.nbins

    GetBond(traj, top, stride, NP, DOP, r1, r2, nbins)
 
