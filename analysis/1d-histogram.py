#========================#
# Calculate 1D histogram of specific atom types, given in -atoms flag 
# modified version of Kevin's
#========================#
# Example usage:
# python 1d-histogram.py trajectory.dcd topology.pdb -axis 2 -atoms <atom_names>

import mdtraj,os,re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colormap = 'coolwarm'
import argparse as ap
parser = ap.ArgumentParser(description="Get 1D histogram, assuming single species")

parser.add_argument('coordfile',type=str, help="trajectory file")
parser.add_argument('topfile',type=str, help="topology file")
parser.add_argument('-axis',type=int, default=0, choices=[0,1,2], help="axis to bin")
parser.add_argument('-atoms', type = str, nargs = '+', help='names of atoms to histogram')
parser.add_argument('-Lx',type=float, default=10, help="Lx, default 10")
parser.add_argument('-Ly',type=float, default=10, help="Ly, default 10")
parser.add_argument('-Lz',type=float, default=10, help="Lz, default 10")
parser.add_argument('-nbins',type=int, default=1000, help="Number of bins for 1d histogram")
parser.add_argument('-nbins2',type=int, default=100, help="Number of bins for probability vs density plot")
parser.add_argument('-rup', type=float, help="max coordinate to calculate probability of density")
parser.add_argument('-rlow', type=float, help="min coordinate to calculate probability of density")
parser.add_argument('-com',action='store_true', help="substract median value of position in each frame")
parser.add_argument('-stride',type=int, default=1, help="stride")
parser.add_argument('-intv', type = int, default=None, help="interval to calculate average density profile, only available with -atoms flag")
args = parser.parse_args()

ax = args.axis
coordfile = args.coordfile
topfile = args.topfile
anames = args.atoms
com = args.com
stride = args.stride
intv = args.intv

print("... Loading Trajectory ...")
traj = mdtraj.load(coordfile,top=topfile,stride=stride)
top = traj.topology
print("... Done Loading {} frames ...".format(traj.n_frames))
Lx,Ly,Lz = traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2] #assuming constant box shape
box = np.array([traj.unitcell_lengths[0,0], traj.unitcell_lengths[0,1], traj.unitcell_lengths[0,2]]) #assuming constant box shape

V   = Lx*Ly*Lz
Ntot = traj.xyz.shape[1]
nbins = args.nbins
n_frames = traj.n_frames

print("box: [{},{},{}]".format(Lx, Ly, Lz))
print("V: {}".format(V))
L = box[ax]
A = np.prod(box)/L

xs = traj.xyz[:,:,ax]
if com:
    xmedian = np.median(xs,1)
    xs = xs - xmedian[:,None]
xs = np.ravel(xs)
dx = L/nbins
xs = np.mod(xs,L) #wrap pbc
#histA,bins = np.histogram(xs, bins=100, density=False)
#histA = histA/traj.n_frames/A/dx
histAll,bins = np.histogram(xs, bins=nbins, density=False)
binRange = (xs.min(), xs.max())

#average over frames
histAll = np.array(histAll)
histAll = histAll/ n_frames

binmidAll = 0.5*(bins[1:]+bins[0:-1])
Vbin = dx * A

#get number density
densityAll = histAll/Vbin
data = np.vstack([binmidAll,histAll,densityAll]).T
np.savetxt('{}hist_all.dat'.format(['x','y','z'][ax]),data,header='bin-midpt\tCount\tNumberDensity')

#histogramming masked trajectory
if anames:
    print("Masking trajectory with {} atom names".format(anames))
    maskId = []
    for aname in anames:
        for atom in top.atoms_by_name(aname):
            maskId.append(atom.index)
    traj = traj.atom_slice(maskId)
    traj.save('masked_'+coordfile)  
    traj[0].save('masked_top_'+coordfile.split('.')[0]+'.pdb')
    n_frames = traj.n_frames
    if not intv:
        intv = n_frames
    NtotMasked = traj.xyz.shape[1]

    nSeries = n_frames // intv 
    if n_frames%intv != 0:
        nSeries += 1
    data = []
    header = ''
    Xdata = []
    XdataHeader = ''
    colors = [matplotlib.cm.get_cmap(colormap)(x/nSeries) for x in range(nSeries+1)]
    for i in range(0,nSeries): 
        if n_frames%intv != 0 and i == nSeries-1: #plot all frames
            xs = traj.xyz[:,:,ax]
            norm = traj.xyz[:,:,ax].shape[0]
        else:
            xs = traj.xyz[intv*i:,:,ax]
            norm = traj.xyz[intv*i:,:,ax].shape[0]
        if com:
            xmedian = np.median(xs,1)
            xs = xs - xmedian[:,None]
        xs = np.ravel(xs)
        xs = np.mod(xs,L) #wrap pbc

        histMasked,bins = np.histogram(xs, range = binRange,  bins=nbins, density=False)
        histMasked = np.array(histMasked)
        histMasked = histMasked/norm

        binmidMasked = 0.5*(bins[1:]+bins[0:-1])
        if all(abs(binmidMasked - binmidAll) < 1e-3):
           Exception('Bin values between of masked traj and original traj do not match')

        #get number density
        density = histMasked/Vbin

        # probability of the density
        if args.rlow and args.rup:
            xLow = args.rlow
            xUp = args.rup
        else:
            xLow = min(binmidMasked)
            xUp = max(binmidMasked)
        id1 = np.where(binmidMasked > xLow)[0]
        id2 = np.where(binmidMasked < xUp)[0]
        id = np.intersect1d(id2, id1)
        X = np.ravel(density[id])
        dX = (np.max(X)-np.min(X))/(float(args.nbins2))
        P,Xbins = np.histogram(X, bins=args.nbins2, density=True)
        Xbinmid = 0.5*(Xbins[1:]+Xbins[0:-1])
        Xdata.extend([Xbinmid,P])
        XdataHeader += 'NumberDensity(#/Vol) Probability(fr{}-{}) '.format(intv*i-1,n_frames)

        if i == 0:
            data.append(binmidMasked)
            #plt.set_cmap = matplotlib.cm.get_cmap(name='coolwarm')
            color = plt.cm.coolwarm(np.linspace(0.05,0.95,nSeries)) # This returns RGBA; convert:
#            hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
#               tuple(color[:,0:-1]))
            #matplotlib.rcParams['axes.color_cycle'] = hexcolor
            fig,(ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize=[3,5])

        data.extend([histMasked,density])
        header += 'Count(fr0-{a})\tNumberDensity(fr0-{a})\t'.format(a=intv*i-1)

        if i % 5 == 0:
            ax1.plot(binmidMasked, density, c = colors[i], marker=None,ls='-',lw=1.25, mfc="None",ms=2,label = 'fr{}-{}'.format(intv*i-1,n_frames))
            ax2.plot(Xbinmid, P, c = colors[i],marker = None, ls='-',lw=1.25, mfc="None",ms=2,label = 'fr{}-{}'.format(intv*i-1,n_frames))
        else:
            ax1.plot(binmidMasked, density, c = colors[i], marker=None,ls='-',lw=1.25, mfc="None",ms=2)
            ax2.plot(Xbinmid, P, c = colors[i], marker = None, ls='-',lw=1.25, mfc="None",ms=2)
    ax1.legend(loc='best',prop={'size': 5})
    ax1.set_xlabel('{}'.format(['x','y','z'][ax]))
    ax1.set_ylabel('Number Density')

    ax2.legend(loc='best',prop={'size': 5})
    ax2.set_ylabel('$P(\\rho)$')
    ax2.set_xlabel('Number Density')

    plt.set_cmap = matplotlib.cm.get_cmap(name='coolwarm')
    fig.tight_layout()
    title ='{}histogram_{}'.format(['x','y','z'][ax],' '.join(args.atoms))
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

    data = np.vstack(data).T
    np.savetxt('{}hist_{}.dat'.format(['x','y','z'][ax],'_'.join(anames)),data,header='Masked atoms: {}\n bin-midpt\t{}'.format(args.atoms,header))
    Xdata = np.vstack(Xdata).T
#    np.savetxt('prob_{}hist.dat'.format(['x','y','z'][ax], Xdata, header='Masked atoms: {}\n{}'.format(args.atoms, XdataHeader)))
else:
    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(binmidAll, densityAll, marker=None,ls='-',lw=0.75,mfc="None",ms=2)
    #axs.legend(loc='best',prop={'size': 5})
    plt.xlabel('{}'.format(['x','y','z'][ax]))
    plt.ylabel('Number Density')
    title ='{}histogram_All'.format(['x','y','z'][ax])
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
