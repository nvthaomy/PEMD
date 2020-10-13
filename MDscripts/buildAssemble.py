#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:03:04 2020

@author: nvthaomy
"""

import sys, os, time
from subprocess import call
sys.path.append('/home/mnguyen/bin/PEMD/MDscripts/')
import singleChain
sys.path.append('/home/mnguyen/bin/scripts/')
import mdtraj_tac

# name map: tleap nomenclature to assemble
#sequenceInTleap sequenceInAssemble monomerpdbname
PAAmap={
'AHP':     ['h',   'AHP.pdb'],
'ATP':     ['t',   'ATP.pdb'],
'un':      ['u',   'APup.pdb'],
'dn':      ['d',   'APdown.pdb'],
'AHD':     ['H',   'AHD.pdb'],
'ATD':     ['T',   'ATD.pdb'],
'ui':      ['U',   'ADup.pdb'],
'di':      ['D',   'ADdown.pdb']}
PAHmap = {
'NHD':     ['h',   'NHD.pdb'],
'NTD':     ['t',   'NTD.pdb'],
'un':      ['u',   'NDup.pdb'],
'dn':      ['d',   'NDdown.pdb'],
'NHP':     ['H',   'NHP.pdb'],
'NTP':     ['T',   'NTP.pdb'],
'ui':      ['U',   'NPup.pdb'],
'di':      ['D',   'NPdown.pdb']
}

pdbLib = '/home/mnguyen/bin/PEMD/monomer4assemble/'
assemblePath = '/home/mnguyen/bin/assemble/Assemble.py'
f = 1.0
N = 24
Pm = 0.44
nChain = 30 # 1 if Pm = 1 (isotactic) or 0 (syndiotactic)
monName = 'AH'
nameExt = '_a' #extension to chain name
pattern = 'even' # 'even' or random' #random deprotonation

################
neutralMons = []
ionizedMons = []
libs =  []
ff = []
    
if monName == 'AH':
    map = PAHmap
    headName = 'NH'
    tailName = 'NT'
    database = 'AHdatabase.txt'
    input = 'PAHin.txt'
    molName = 'PAH'
elif monName == 'AA':
    map = PAAmap
    headName = 'AH'
    tailName = 'AT'
    database = 'AAdatabase.txt'
    input = 'PAAin.txt'
    molName = 'PAA'

sequences = []
print('===========================')
print('...Generating sequences...')
print('===========================')
for i in range(nChain):
    _,_,_, sequence = singleChain.BuildPoly(f,N,Pm,pattern, monName, neutralMons, ionizedMons, headName,tailName,libs,ff,runTleap=False)
    sequences.append(sequence)

# write input files for assemble
s = ''
for key,(mon,pdb) in map.items():
    s+= '{}  {}\n'.format(mon,os.path.join(pdbLib,pdb))
file = open(database,'w')
file.write(s)
file.close()

print('=======================')
print('...Running Assemble...')
print('=======================')

s = """
mode pdb
database {database}
output_folder .
system_name {molName}
box_grid_shape 5 5 5\n\n""".format(database=database, molName=molName)

chainNames = []
for i, seq in enumerate(sequences):
    if f == 1.0 or f == 0.0:
        DOIstr = 'f'+str(int(f))
    else:
        DOIstr = 'f'+str(round(f,1)) 
    chainName = monName + str(N) + DOIstr + '{}'.format(nameExt) + str(i)
    chainNames.append(chainName)
    s += 'chain {} '.format(chainName)
    for key in seq.split():
        mon = map[key][0]
        s += mon
    s += '\n'
s += 'molecule {}'.format(' '.join(chainNames))   
file = open(input,'w')
file.write(s)
file.close() 

print('python {} {}'.format(assemblePath, input))
call('python {} {}'.format(assemblePath, input), shell=True)

cwd = os.getcwd()
try:
    os.mkdir('minimizeE')
except:
    pass
os.chdir('minimizeE')

print('=================')
print('... Minimizing...')
print('=================')

# write packmol
nIon = 0
if molName == 'PAA':
    nIon= int(f*N)
    ion = 'na'
elif molName == 'PAH':
    nIon = int(f*N)
    ion = 'cl'
   
tactLog = 'tacticity.txt'
log = open(tactLog,'w')
log.write('# chainpdb mesofaction tacticity\n')
log.flush()

for chainName in chainNames:
    print('\n ... Get minimized structure for chain {} ...'.format(chainName))
    s="""
tolerance 2.0
filetype pdb
output init.pdb
structure {cwd}/{molName}/{chainName}.pdb
        number 1
        resnumbers 2
        inside box 1 1 1 60 60 60
        end structure
structure /home/mnguyen/bin/PEMD/{ion}.pdb
        number {nIon}
        resnumbers 2
        inside box 1 1 1 60 60 60
        end structure      
structure /home/mnguyen/bin/PEMD/opc.pdb
        number 500
        resnumbers 2
        inside box 1 1 1 60 60 60
        end structure""".format(cwd=cwd, molName = molName, chainName = chainName, ion=ion, nIon=nIon)
    file=open('mix.inp','w') 
    file.write(s)
    file.close()
    
    print('Packing molecules ...')
    os.system('packmol < {} > {}'.format('mix.inp', 'packmol.log'))
    finished = False
    while not finished:
        f = open('packmol.log','r')
        if 'Success' in f.read():
            finished = True
        time.sleep(1)
    
    # tleap to get parm7 and crd
    s = """
source leaprc.gaff2
source oldff/leaprc.ff99
source leaprc.water.opc
loadOFF /home/mnguyen/bin/PEMD/MDscripts/lib/PAAP9N9c5.lib
loadOFF /home/mnguyen/bin/PEMD/MDscripts/lib/PAAD9N9c4.lib
loadOFF /home/mnguyen/bin/PEMD/MDscripts/lib/PAHP1N11c5.lib
loadOFF /home/mnguyen/bin/PEMD/MDscripts/lib/PAHD1N11c6.lib
x=loadpdb init.pdb
addions x Na+ 0
addions x Cl- 0
setbox x vdw 1
saveamberparm x top.parm7 crd
savepdb x init.pdb
quit"""
    
    file = open('loadFF.in','w')
    file.write(s)
    file.close()
    print('Loading forcefield ...')
    os.system('tleap -s -f  loadFF.in > loadFF.out')

    s = """   
import sys
#sys.path.append('/home/mnguyen/miniconda3/')
from sys import stdout
import numpy as np
import mdtraj, simtk
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import parmed
from parmed import gromacs
print(gromacs.__file__)

# Input Files

top_file = 'top.parm7'
incoord = 'crd'
prmtop = AmberPrmtopFile(top_file)
inpcrd = AmberInpcrdFile(incoord)
TrajFile = 'trajectory.dcd'
ThermoLog = 'log.txt'

# System Configuration
nonbondedMethod = LJPME
nonbondedCutoff = 1.0*nanometers
tailCorrection = False
ewaldErrorTolerance = 0.0001
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
box_vectors = np.diag([6.,6.,6.]) * nanometer

#Integration Options
dtsim = 0.001*picoseconds
temperature = 298.15*kelvin
meltTemperature = 320.0*kelvin
pressure = 1.0*atmospheres
barostatInterval = 25

#Simulation Options
steps = 3000
platform = Platform.getPlatformByName('CPU')

dcdReporter = mdtraj.reporters.DCDReporter(TrajFile, 100)
dataReporter = StateDataReporter(ThermoLog, 100, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

# Prepare the Simulation
print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions
#pdb = PDBFile(initpositions)
#positions = pdb.positions

def MakeSimulation(temperature, dt = 0.002*picoseconds):
    system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
    system.addForce(MonteCarloBarostat(pressure,temperature,barostatInterval))
    friction = 1.0/(100.*dt)
    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    simulation = Simulation(topology, system, integrator, platform) #, platformProperties)
    simulation.context.setPositions(positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.context.setPeriodicBoxVectors(*box_vectors)
    forces = system.getForces()
    for force in forces:
        if isinstance(force,simtk.openmm.openmm.NonbondedForce):
            nonbondedforce = force
    nonbondedforce.setUseDispersionCorrection(tailCorrection)
    nonbondedforce.updateParametersInContext(simulation.context)
    forces = system.getForces()
    for force in forces:
        if isinstance(force,simtk.openmm.openmm.NonbondedForce):
            nonbondedforce = force
    print('getUseDispersionCorrection')
    print(nonbondedforce.getUseDispersionCorrection())
    return simulation

#Restart and Check point
#to load state: simulation.loadState('output.xml')

if pressure == None:
    print('===Running NVT simulation===')
else:
    print('===Running NPT simulation===')

simulation = MakeSimulation(temperature, dt=dtsim)
print('Minimizing...')
simulation.minimizeEnergy() #maxIterations=1000)
simulation.context.setVelocitiesToTemperature(temperature)
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)

#Production
print('Simulating...')
simulation.reporters = []
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(CheckpointReporter('checkpnt.chk',50))
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
PDBFile.writeModel(simulation.topology, positions.value_in_unit(nanometer), open('minimized.pdb','w'))

# get pdb of chain
t=mdtraj.load('trajectory.dcd',top=top_file)
t=t[-1]
top = t.topology
s = top.select("resname != SOL and resname != WAT and resname != HOH and resname != 'NA+' and resname != 'Na+' and resname != 'Cl-' and resname != 'CL-'")
t1=t.atom_slice(s)
t1.save('{chainName}.pdb')
print('Finish')
""".format(chainName=chainName)
     
    file = open('sim.py','w')
    file.write(s)
    file.close()
    call('python sim.py', shell=True)
    # call('python sim.py > sim.log', shell=True)
    # finished = False
    # while not finished:
    #     f = open('sim.log','r')
    #     if 'Finish' in f.read():
    #         finished = True
    #     time.sleep(1)
    
    # get actual tacticity
    coordfile = chainName+'.pdb'
    topfile = chainName+'.pdb'
    resrange=  [0,N-1]
    frames = [-1]
    topdat = None
    stride = 1
    warmup = 0
    
    while not os.path.exists(coordfile):
        time.sleep(1)
    
    # update resid so that the first resid is 1
    FixResid = False
    file = open(coordfile,'r')
    lines = file.readlines()
    linenum = 0
    newlines = []
    for line in lines:
        if 'ATOM' in line or 'HETATM' in line or 'TER' in line:
            resid = int(line[22:26])
            if linenum == 0:
                if resid == 0:
                    FixResid = True
                    newid = resid + 1
                    line = line[:22] + ' '*(4-len(str(newid))) + str(newid) + line[26:]
                    print('-- Updating residue id for {} --'.format(fileName))
                elif resid!= 0:
                    SkipToNextFile = True
            elif linenum !=0 and FixResid:
                newid = resid + 1
                line = line[:22] + ' '*(4-len(str(newid))) + str(newid) + line[26:]
            linenum += 1
        newlines.append(line)
    
        if SkipToNextFile:
            break
    
    if not SkipToNextFile:
        newfile = open(coordfile,'w')
        newfile.write(''.join(newlines))
        newfile.close()

    fm,tact = mdtraj_tac.GetTac(coordfile,topfile, resrange, frames, topdat, stride, warmup)
    fm = fm[0]
    tac = tact[0]
    
    log = open(tactLog,'a')
    log.write('{} {} {} \n'.format(chainName, round(fm,5), tact))
    log.flush()

    
    
    
