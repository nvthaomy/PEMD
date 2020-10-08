import os, sys, re
import numpy as np
import mdtraj as md

resDict = {
'PAA': ['ATP','AHP','AP','ATD','AHD','AD'],
'PAH': ['NTP','NHP','NP','NTD','NHD','ND'], 
'Na+': ['Na+'],
'CB' : ['CB'],
'CC' : ['CC'],
'O': ['O'],
'A-': ['A-']}
#id of first residue in 'for' loop
res0Id = 0

def GetRDF(traj, top, atom1, atom2, nbins, Ls_in, rmax):
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

    Ls = traj.unitcell_lengths
    if len(Ls) == 0:
        if Ls_in == None:
            raise Exception('No box dimensions in trajectory, need to manually input box dimension with flag -L')
        else:
            Ls = traj.n_frames * [Ls_in]
            Ls = np.array(Ls)

    atom1Ids = top.select("name '{}'".format(atom1))
    atom2Ids = top.select("name '{}'".format(atom2))
    xyz1 = traj.atom_slice(atom1Ids).xyz #positions in the central image
    xyz1 = np.transpose(xyz1,(1,0,2))
    xyz2 = traj.atom_slice(atom2Ids).xyz #positions in the central image
    xyz2 = np.transpose(xyz2,(1,0,2))

    #get minimum image
    Ls_extend = np.array(len(xyz1)*[Ls])
    xyz1 = xyz1 - Ls_extend * (xyz1/Ls_extend).astype(int)
    xyzs = [] #array of atom positions for images of ix,iy,iz = -1,0,1 (total of 27 images) 

    #get nearest images, if rmax < half box dimension, won't have to account for overcounting?
    for ix in [-1.,0.,1.]:
        for iy in [-1.,0.,1.]:
            for iz in [-1.,0.,1.]:
                xyz_temp = xyz1 + np.array([ix,iy,iz]) * Ls
                xyzs.append(xyz_temp)
    xyzs = np.array(xyzs) #dimensions 27 x nMol x nFrames x 3
    nImg = len(xyzs)
    nAtom = len(xyzs[0])
    #modify Ls so that dimension is the same as xyzs
    Ls_extend = nImg * [nAtom*[Ls]]
    Ls = np.array(Ls_extend)

    print('... computing distance between {} and {} ...'.format(atom1, atom2))
    ds = []
    for i,r2 in enumerate(xyz2):
        ds_tmp = []
        ds_tmp.extend(r2-xyzs)
        ds_tmp = np.array(ds_tmp) #dimensions 27 x nAtom x nFrames x 3
        ds_tmp = np.sqrt(ds_tmp[:,:,:,0]**2 +ds_tmp[:,:,:,1]**2 + ds_tmp[:,:,:,2]**2)
        # get minimum distance from 27 images
        ds_tmp = np.amin(ds_tmp,axis=0)
        ds_tmp = ds_tmp[np.where(ds_tmp<=rmax)] #dimensions nAtom x nFrames x 3
        ds.extend(ds_tmp)

    #histogram    
    ds = np.ravel(np.array(ds)) 
    hist,bins = np.histogram(ds, range=(0.,rmax), bins=nbins, density=False)
    rs = np.array(bins[0:-1])
    dr = rmax/nbins
    
    #get rdf: normalize histogram by number of atoms in the shell at distance r as if atoms are uncorrelated
    #normalization = (avg number density) * 4 pi r^2 dr * traj.n_frames
    rho = len(atom1Ids)/traj.unitcell_volumes.mean()
    print('mean vol {} nm3'.format(traj.unitcell_volumes.mean()))
#    norms = (traj.n_frames * len(xyz2)) #normalize by number of frames and number of atom2
    norms = traj.n_frames #normalize by number of frames 
#    norms = rho * 4* np.pi * rs**2 * dr * (traj.n_frames * len(xyz2)) # * len(xyzs)/len(atomIds)) #later 2 factors account for overcounting from multiple frames, multiple chains
    gs = hist /norms

    data = np.vstack([rs,gs]).T
    np.savetxt('dist_{}_{}.dat'.format(atom1,atom2), data, header= 'Distance Count', fmt='%.5e')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(rs, gs, marker=None,ls='-',lw=1, c = 'b')
    plt.xlim(0., rmax)
    plt.ylim(0., 1.1*gs[5:].max())
    plt.xlabel('Distance (nm)')
    plt.ylabel('Count {}-{}'.format(atom1, atom2))
    title = 'dist {} {}'.format(atom1, atom2)
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
    return rs,gs
       
if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('atom1',  type = str, help = 'atom to get rdf')
    parser.add_argument('atom2',  type = str, help = 'atom to get rdf')
    parser.add_argument('-rmax', type = float, default = 1., help = 'max distance in nm')
    parser.add_argument('-L', nargs = 3, default = None, type=float, help='box dimensions in nm if tracjectory does not have box info')
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-nbins', type=float, default=500)
    args = parser.parse_args()

    traj = sys.argv[1]
    top = sys.argv[2]
    atom1 = sys.argv[3]
    atom2 = sys.argv[4]
    rmax = args.rmax
    stride = args.stride
    nbins = args.nbins
    Ls_in = args.L

    print('... loading trajectory ...')
    traj = md.load(traj,top = top, stride=stride)
    top = traj.topology
    print('... done loading ...')
    GetRDF(traj, top, atom1, atom2, nbins, Ls_in, rmax)
 
