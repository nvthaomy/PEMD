import os, sys, re
import numpy as np
import mdtraj as md

resDict = {
'PAA': ['ATP','AHP','AP','ATD','AHD','AD'],
'PAH': ['NTP','NHP','NP','NTD','NHD','ND']}    
#id of first residue in 'for' loop
res0Id = 0

def GetPEDist(traj, top, stride, MolNames, NPs, DOPs, nbins, Ls_in):
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

    AtomsInMolDict = {}
    for i,NP in enumerate(NPs):
        DOP = DOPs[i]
        MolName = MolNames[i]
        AtomsInMol = []
        AtomsInRes = []
        resCount = 0
        for res in top.residues:
            if res.name in resDict[MolName]:
                resCount += 1
                #get atom indices of all atoms in this residue
                for atom in res.atoms:
                    AtomsInRes.append(atom.index)
                if resCount%DOP == 0.:
                    AtomsInMol.append(AtomsInRes)
                    AtomsInRes = []
        AtomsInMolDict.update({MolName:AtomsInMol})

    COMDict = {}
    for MolName,AtomsInMol in AtomsInMolDict.items():
        COMs = []
        for atomId in AtomsInMol: #loop through molecules in with this MolName
            COM = md.compute_center_of_mass(traj.atom_slice(atomId)) #time series of COM for this residue 
            COMs.append(COM)
        COMDict.update({MolName:COMs})

    print('=== computing distance between PAA and PAH COMs ===')
    ds_t = [] #time series of distance
    
    Ls = traj.unitcell_lengths
    if len(Ls) == 0:
        if Ls_in == None:
            raise Exception('No box dimensions in trajectory, need to manually input box dimension with flag -L')
        else:
            Ls = traj.n_frames * [Ls_in]
            Ls = np.array(Ls)
    for i in range(len(COMDict[MolNames[0]])):
        PAACOM = COMDict[MolNames[0]][i]
        if MolNames[0] == MolNames[1]: #if calculating COM distance of same species of molecule, ignore the distance of the same molecule which is 0 and ignore pairs that were calculated
            j0 = i+1
        else: 
            j0 = 0
        for j in range(j0,len(COMDict[MolNames[1]])): 
            PAHCOM0 = COMDict[MolNames[1]][j]
            d_temp = []
            # get possible COM for nearest images
            # x,y,z = x0 + i_x*L_x, y0 + i_y*L_y, z0 + i_z*L_z
            for ix in [-1.,0.,1.]:
                for iy in [-1.,0.,1.]:
                    for iz in [-1.,0.,1.]:
                        PAHCOM = PAHCOM0 + [ix,iy,iz] * Ls           
                        d = (PAACOM-PAHCOM)
                        d = np.sqrt(d[:,0]**2 +d[:,1]**2 + d[:,2]**2)            
                        d_temp.append(d)
            #use the minimum distance
            ds_t.append(np.nanmin(np.array(d_temp),axis = 0))

    #get trajectory of distance between two PEs, averaging over PAA-PAH pairs
    ds_t_avg = np.sum(ds_t, axis = 0)/float(len(ds_t))
    ds_t_avg_err = np.std(ds_t, axis = 0)/np.sqrt(float(len(ds_t)))
    xs = range(len(ds_t_avg))
    np.savetxt('dist_PE.dat', np.vstack([xs, ds_t_avg, ds_t_avg_err]).T, header = 'frame AvgDistance Err', fmt='%.5e')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    step = int(len(ds_t_avg)/200.)
    axs.plot(xs[::step], ds_t_avg[::step], marker=None,ls='-',lw=1, c = 'g')
    axs.fill_between(xs[::step], ds_t_avg[::step] - ds_t_avg_err[::step],  ds_t_avg[::step] + ds_t_avg_err[::step], interpolate=True, color = 'g', alpha = 0.4, linewidth=0.)
    plt.xlim(0)
    plt.xlabel('Frame')
    plt.ylabel('Distance (nm)')
    title = 'Distance PE Trajectory'
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

    #histogram    
    ds = np.ravel(ds_t) 
    hist,bins = np.histogram(ds, bins=nbins, density=True)
    binmid = 0.5*(bins[1:]+bins[0:-1])
    #normalize histogram
#    hist = hist/float(sum(hist))
    print('Max separation {:3.2f} nm'.format(ds.max()))
    data = np.vstack([binmid,hist]).T
    np.savetxt('dist_PE_hist.dat', data, header= 'Distance ProbabilityDensity', fmt='%.5e')

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(binmid, hist, marker=None,ls='-',lw=1, c = 'k')
    plt.xlim(0.9*ds.min(), 1.1*ds.max())
#    plt.ylim(0)
    plt.xlabel('Distance (nm)')
    plt.ylabel('Histogram')
    title = 'Distance PE Histogram'
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
    return ds
       
if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('-np', nargs='+', type = float, help = 'number of chains')
    parser.add_argument('-dop', nargs='+', type = float, help = 'DOP')
    parser.add_argument('-mol', nargs='+', type=str, help='molecules to calculate COM separation')
    parser.add_argument('-L', nargs = 3, default = None, type=float, help='box dimensions in nm if tracjectory does not have box info')
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-nbins', type=float, default=100)
    args = parser.parse_args()

    traj = sys.argv[1]
    top = sys.argv[2]
    MolNames = args.mol
    NPs = args.np 
    DOPs = args.dop
    stride = args.stride
    nbins = args.nbins
    Ls_in = args.L
    GetPEDist(traj, top, stride, MolNames, NPs, DOPs, nbins, Ls_in)
 
