import os, sys, re
import numpy as np
import mdtraj as md

resDict = {
'PAA': ['ATP','AHP','AP','ATD','AHD','AD'],
'PAH': ['NTP','NHP','NP','NTD','NHD','ND'], 
'Na+': ['Na+'],
'CB' : ['CB'],
'CC' : ['CC'],
'O': ['O']}
#id of first residue in 'for' loop
res0Id = 0

#def GetCOMs(traj,top,stride,MolName,DOP,Ls_in):
#    AtomsInMol = []
#    AtomsInRes = []
#    resCount = 0
#    for res in top.residues:
#        if res.name in resDict[MolName]:
#            resCount += 1
#            #get atom indices of all atoms in this residue
#            for atom in res.atoms:
#                AtomsInRes.append(atom.index)
#            if resCount%DOP == 0.:
#                AtomsInMol.append(AtomsInRes)
#                AtomsInRes = []
#    COMs = [] #rows = chains, columns = frames
#    for atomId in AtomsInMol: #loop through molecules in with this MolName
#        COM = md.compute_center_of_mass(traj.atom_slice(atomId)) #time series of COM for this chain 
#        COMs.append(COM)
#    COMs = np.array(COMs)
#    Ls = traj.unitcell_lengths
#    if len(Ls) == 0:
#        if Ls_in == None:
#            raise Exception('No box dimensions in trajectory, need to manually input box dimension with flag -L')
#        else:
#            Ls = traj.n_frames * [Ls_in]
#            Ls = np.array(Ls)
#    #modify Ls so that dimension is the same as COMs
#    Ls_extend = np.array(len(COMs)*[Ls])
#    # get minimum image
#    COMs = COMs - Ls_extend * (COMs/Ls_extend).astype(int) 
#    return traj, top, COMs, Ls

def GetCOMs(traj,top,stride,MolName,DOP,Ls_in):
    AtomsInMol = []
    AtomsInRes = []
    resCount = 0
    for atom in top.atoms:
        if atom.name in resDict[MolName]:
            resCount += 1
            AtomsInRes.append(atom.index)
            if resCount%DOP == 0.:
                AtomsInMol.append(AtomsInRes)
                AtomsInRes = []
    COMs = [] #rows = chains, columns = frames
    for atomId in AtomsInMol: #loop through molecules in with this MolName
        COM = md.compute_center_of_mass(traj.atom_slice(atomId)) #time series of COM for this chain 
        COMs.append(COM)

    Ls = traj.unitcell_lengths
    if len(Ls) == 0:
        if Ls_in == None:
            raise Exception('No box dimensions in trajectory, need to manually input box dimension with flag -L')
        else:
            Ls = traj.n_frames * [Ls_in]
            Ls = np.array(Ls)
    #modify Ls so that dimension is the same as COMs
    Ls_extend = np.array(len(COMs)*[Ls])
    COMs = COMs - Ls_extend * (COMs/Ls_extend).astype(int)
    return traj, top, COMs, Ls

def GetRDF(traj, top, COMs, MolName, atomName, nbins, Ls, rmax):
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
    
    atomIds = top.select("name '{}'".format(atomName))
    xyz0 = traj.atom_slice(atomIds).xyz #positions in the central image
    xyz0 = np.transpose(xyz0,(1,0,2))
    #get minimum image
    Ls_extend = np.array(len(xyz0)*[Ls])
    xyz0 = xyz0 - Ls_extend * (xyz0/Ls_extend).astype(int)
    xyzs = [] #array of atom positions for images of ix,iy,iz = -1,0,1 (total of 27 images) 

    #get nearest images, if rmax < half box dimension, won't have to account for overcounting?
    for ix in [-1.,0.,1.]:
        for iy in [-1.,0.,1.]:
            for iz in [-1.,0.,1.]:
                xyz_temp = xyz0 + np.array([ix,iy,iz]) * Ls
                xyzs.append(xyz_temp)
#    xyzs.extend(xyz0)
    xyzs = np.array(xyzs) #dimensions 27 x nMol x nFrames x 3
    nImg = len(xyzs)
    nAtom = len(xyzs[0])
    #modify Ls so that dimension is the same as xyzs
    Ls_extend = nImg * [nAtom*[Ls]]
    Ls = np.array(Ls_extend)

    print('... computing distance between {} COMs and {} ...'.format(MolName, atomName))
    ds = []
    for i,COM in enumerate(COMs):
#        print('chain {}'.format(i))
        COM_tmp = np.array(nImg * [COM])
        ds_tmp = []
        ds_tmp.extend(COM-xyzs)
        ds_tmp = np.array(ds_tmp) #dimensions 27 x nAtom x nFrames x 3
#        ds_tmp = ds_tmp - Ls* (ds_tmp/Ls).astype(int)
        ds_tmp = np.sqrt(ds_tmp[:,:,:,0]**2 +ds_tmp[:,:,:,1]**2 + ds_tmp[:,:,:,2]**2)
#        print('ds_tmp before minimum shape {}'.format(ds_tmp.shape))
        # get minimum distance from 27 images
        ds_tmp = np.amin(ds_tmp,axis=0)
#        print('ds_tmp after minimum shape {}'.format(ds_tmp.shape))
        ds_tmp = ds_tmp[np.where(ds_tmp<=rmax)] #dimensions nAtom x nFrames x 3
#        print('ds_tmp after filter shape {}'.format(ds_tmp.shape))
        ds.extend(ds_tmp)

    #histogram    
    ds = np.ravel(np.array(ds)) 
    hist,bins = np.histogram(ds, range=(0.,rmax), bins=nbins, density=False)
    rs = np.array(bins[0:-1])
    dr = rmax/nbins
    
    #get rdf: normalize histogram by number of atoms in the shell at distance r as if atoms are uncorrelated
    #normalization = (avg number density) * 4 pi r^2 dr * traj.n_frames
    rho = len(atomIds)/traj.unitcell_volumes.mean()
    print('mean vol {} nm3'.format(traj.unitcell_volumes.mean()))
    norms = rho * 4* np.pi * rs**2 * dr * (traj.n_frames * len(COMs)) # * len(xyzs)/len(atomIds)) #later 2 factors account for overcounting from multiple frames, multiple chains
    gs = hist/norms

    data = np.vstack([rs,gs]).T
    np.savetxt('rdf_{}COM_{}.dat'.format(MolName,atomName), data, header= 'Distance RDF', fmt='%.5e')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(rs, gs, marker=None,ls='-',lw=1, c = 'k')
    plt.xlim(0., rmax)
    plt.ylim(0.)
    plt.xlabel('Distance (nm)')
    plt.ylabel('rdf {}-{}'.format(MolName,atomName))
    title = 'rdf {}COM {}'.format(MolName,atomName)
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
    return rs,gs
       
if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('mol',  type = str, help = 'PAA or PAH')
    parser.add_argument('atom',  type = str, help = 'atom to get rdf wrt polymer COMs')
    parser.add_argument('-rmax', type = float, default = 1., help = 'max distance in nm')
    parser.add_argument('-dop', required=True, type = float, help = 'DOP')
    parser.add_argument('-L', nargs = 3, default = None, type=float, help='box dimensions in nm if tracjectory does not have box info')
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-nbins', type=float, default=500)
    args = parser.parse_args()

    traj = sys.argv[1]
    top = sys.argv[2]
    MolName = sys.argv[3]
    atomName = sys.argv[4]
    rmax = args.rmax
    DOP = args.dop
    stride = args.stride
    nbins = args.nbins
    Ls_in = args.L

    print('... loading trajectory ...')
    traj = md.load(traj,top = top, stride=stride)
    top = traj.topology
    print('... done loading ...')
    traj, top, COMs, Ls = GetCOMs(traj,top,stride,MolName,DOP,Ls_in)
    GetRDF(traj, top, COMs, MolName, atomName, nbins, Ls, rmax)
 
