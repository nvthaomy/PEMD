import numpy as np
import sys, re
sys.path.append('/home/mnguyen/bin/scripts')
import stats
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('log', help = 'log file from omm')
parser.add_argument('template', help = 'template file to write simulation script')
parser.add_argument('-pe', type=str, help = 'header of PE')
parser.add_argument('-ke', type=str, help = 'header of KE')
parser.add_argument('-e', type=str, help = 'header of total energy')
parser.add_argument('-v', type=str, help = 'header of volume')
parser.add_argument('-x',type=float,  help= 'box dimension in the x direction (nm)')
parser.add_argument('-y',type=float,  help= 'box dimension in the y direction (nm)')
parser.add_argument('-z',type=float,  help= 'box dimension in the z direction (nm)')
args = parser.parse_args()

log = args.log
template = args.template
Vname = args.v
PEname = args.pe
KEname = args.ke
Ename = args.e
x = args.x
y = args.y
z = args.z

#check, can only specify 2 dimensions at max
if not None in [x,y,z] and Vname:
    raise Exception('can only specify 2 dimensions at max')

#do stats
def Average(log, obs):
    file = open(log,'r')
    lines = file.readlines()
    for line in lines:
        if line.startswith('#'):
            vals = line.split()
            col = vals.index(obs) - 1
            break
    warmup,Data,nwarmup = stats.autoWarmupMSER(file, col)
    print ("Auto warmup detection with MSER-5 => ",nwarmup)
    (nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data, False ,False,'_{}'.format(obs))
    return mean

if Vname:
    Vmean = Average(log, Vname)
    dict = {'x': x, 'y': y, 'z': z}
    Ls = np.array([x,y,z])
    count = 0
    for d in Ls:
        if d is not None:
            count += 1
    print('Provided {} dimensions'.format(count))
    if count == 0:
        L = Vmean**(1./3.)
        Ls = np.array([L,L,L])
    elif count == 2:
        i = np.where(Ls==None)
        j = [ii for ii in range(len(Ls)) if not ii in i[0]]
        Ls[i] = Vmean/np.product(Ls[j])
    elif count == 1:
        i = np.where(Ls==None)
        j = [ii for ii in range(len(Ls)) if not ii in i[0]]
        L = np.sqrt(Vmean/Ls[j])
        Ls[i] = L
    print('Avg vol is {:4.3f} nm^3'.format(Vmean))

    x,y,z = Ls[0],Ls[1],Ls[2]
 
    with open(template,'r') as myfile:
        ini=myfile.read()
        ini=re.sub('__x__',str(x),ini)
        ini=re.sub('__y__',str(y),ini)
        ini=re.sub('__z__',str(z),ini)
        simfile = open("sim.py","w")
        simfile.write(ini)
        simfile.close()

    writeSim = True
        
if Ename: 
    Emean = Average(log, Ename)
    print('Avg total energy is {:4.3f} kJ/mol'.format(Emean))
    if writeSim:
        template = 'sim.py'
    with open(template,'r') as myfile:
        ini=myfile.read()
        ini=re.sub('__Etar__',str(Emean),ini)
        simfile = open("sim.py","w")
        simfile.write(ini)
        simfile.close()
