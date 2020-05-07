#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:12:58 2019

@author: nvthaomy
"""
import sys, glob
def WriteJobFile(simName, jobName):
    '''write job file'''
    name = 'run.sh'
    s = """#!/bin/bash
#SBATCH --ignore-pbs
#SBATCH --nodes=1 --partition=gpu --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --job-name={jobName}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=my@ucsb.edu
/bin/hostname
srun --gres=gpu:1 /usr/bin/nvidia-smi
cd $SLURM_SUBMIT_DIR

#PBS -q gpuq
#PBS -V
#PBS -j oe
#PBS -N {jobName}
#PBS -M my@ucsb.edu
#PBS -m abe
#cd $PBS_O_WORKDIR

export PATH="/home/mnguyen/miniconda3/bin:$PATH"
python {simName}""".format(jobName=jobName, simName=simName)
    file =  open(name,'w')
    file.write(s)
    return name

def WriteOpenMMInput(top, crd, steps, msteps, csteps, x, y, z, anisoP, stride=5000,chkStride=2000, temp=298.15, pressure=1.0, meltTemp=400.):
    name = 'sim.py'
    if pressure != None:
        pressure = '{}*atmospheres'.format(pressure)
        if anisoP:
            sBaro = 'system.addForce(MonteCarloAnisotropicBarostat(3*[pressure],temperature,False,False,True,barostatInterval))'
        else:
            sBaro =  "system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))"
    elif pressure == None:
         sBaro = ""
    print("""Writing simulation script ...
Melting steps {0} at {1} K
Cooling steps {2}
Production steps {3} at {4} K
Pressure {5} atm
Anisotropic barostat {6}
Box {7} {8} {9} nm""".format(msteps, meltTemp, csteps, steps, temp, pressure, anisoP, x,y,z ))


    s = """from sys import stdout
import numpy as np
import mdtraj, simtk
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

# Input Files
topName = '{top}'
TrajFile = 'trajectory.dcd'
ThermoLog = 'log.txt'
prmtop = AmberPrmtopFile(\'{top}\')
inpcrd = AmberInpcrdFile(\'{crd}\')

# System Configuration
nonbondedMethod = LJPME
nonbondedCutoff = 1.0*nanometers
tailCorrection = True
ewaldErrorTolerance = 0.0001
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
box_vectors = np.diag([{x},{y},{z}]) * nanometer


#Integration Options
dt = 0.002*picoseconds
temperature = {temp}*kelvin
meltTemperature = {meltTemp}*kelvin
friction = 1.0/(100.*dt)
pressure = {pressure}
barostatInterval = 25

#Simulation Options
meltingSteps = {msteps}
steps = {steps}
coolingSteps = {csteps}
platform = Platform.getPlatformByName(\'CUDA\')
platformProperties = {{\'Precision\': \'mixed\'}}

dcdReporter = mdtraj.reporters.DCDReporter(TrajFile, {stride})
#dcdReporter = DCDReporter(\'trajectory.dcd\', {stride})'
dataReporter = StateDataReporter(ThermoLog, {stride}, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=\'\\t\')

dcdReporterCooling = mdtraj.reporters.DCDReporter('trajectory_cool.dcd', {stride})
dataReporterCooling = StateDataReporter('log_cool.txt', {stride}, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

dcdReporterMelting = mdtraj.reporters.DCDReporter('trajectory_melt.dcd', {stride})
dataReporterMelting = StateDataReporter('log_melt.txt', {stride}, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

# Prepare the Simulation
print(\'Building system...\')
topology = prmtop.topology
positions = inpcrd.positions

def MakeSimulation(temperature, pressure):
    system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
    if pressure:
        {sBaro}
    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    simulation = Simulation(topology, system, integrator, platform, platformProperties)
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
#to load state: simulation.loadState(\'output.xml\')

if pressure == None:
    print('===Running NVT simulation===')
else:
    print('===Running NPT simulation===')
#Melt
if meltingSteps > 0:
    simulation = MakeSimulation(meltTemperature,pressure)
    simulation.reporters.append(CheckpointReporter('checkpntMelt.chk', {chkStr}))
    #simulation.loadState('output.xml')
    #simulation.loadCheckpoint('checkpnt.chk')

    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(meltTemperature)

    print('Melting at ',meltTemperature.value_in_unit(kelvin),' K')
    simulation.currentStep = 0
    simulation.reporters.append(dcdReporterMelting)
    simulation.reporters.append(dataReporterMelting)
    simulation.step(meltingSteps)
    simulation.saveState('output_melted.xml')
    positions = simulation.context.getState(getPositions=True).getPositions()
    #PDBFile.writeFile(simulation.topology, positions, open('melted.pdb','w'))


    #Cooling
    simulation = MakeSimulation(temperature,pressure)
    simulation.reporters.append(CheckpointReporter('checkpnt.chk',{chkStr}))
    simulation.loadCheckpoint('checkpntMelt.chk')
    print('Cool down...')
    simulation.reporters.append(dcdReporterCooling)
    simulation.reporters.append(dataReporterCooling)
    simulation.step(coolingSteps)
    simulation.saveState('output_cooled.xml')
    positions = simulation.context.getState(getPositions=True).getPositions()
    #PDBFile.writeFile(simulation.topology, positions, open('cooled.pdb','w'))

else:
    simulation = MakeSimulation(temperature)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)

#Production
print('Simulating...')
simulation.reporters = []
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(CheckpointReporter('checkpnt.chk',1000))
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')

#analyzing
#import sys
#sys.path.append('/home/mnguyen/bin/PAAMD/')
#import analysis

#NAtomsPerChain = None #None if analyzing AA
#NP = 
#top = topName
#DOP = 
#fi = 'openmm' #'lammps' or 'openmm' 
#analysis.getStats(TrajFile, top, NP, ThermoLog, DOP = DOP, NAtomsPerChain = NAtomsPerChain, StatsFName = 'AllStats.dat',
#            RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',
#             fi = fi, obs = [ 'Potential_Energy_(kJ/mole)',   'Kinetic_Energy_(kJ/mole)',   'Temperature_(K)',   'Box_Volume_(nm^3)',   'Density_(g/mL)'], cols = None,
#             res0Id = 0, stride = 1, autowarmup = True, warmup = 100, plot = False)
""".format(top=top, crd=crd, stride=stride, temp=temp, steps=steps, msteps=msteps, csteps=csteps, x=x, y=y, z=z, sBaro=sBaro, chkStr=chkStride, meltTemp=meltTemp, pressure=pressure)
    file = open(name,'w')
    file.write(s)
    file.close()
    return name

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("top", type=str, help='topology')
    parser.add_argument("crd", type=str, help='coordinate')
    parser.add_argument("-jn", type=str, default='simMix', help='job name')
    parser.add_argument("-s", type=int, required=True, help='production steps')
    parser.add_argument("-ms", type=int, required=True, help='melting steps')
    parser.add_argument("-cs", type=int, required=True, help='cooling steps')
    parser.add_argument("-xyz",nargs =3,type = float,required=True,
                        help="periodic box size in nm")
    parser.add_argument("-a", action = "store_true",
                        help="if use anisotropic MC barostat, acts on the z direction")
    parser.add_argument('-stride', type=int, default=5000,help='stride of thermo log and trajectory')
    parser.add_argument('-cstride', type=int, default=2000, help='stride of check point')
    parser.add_argument('-temp', type=float, default=298.15, help='temperature in Kelvin')
    parser.add_argument('-mtemp', type=float, default=400., help='melting temperature in Kelvin')
    parser.add_argument('-p', help='pressure in atm to run NPT, run NVT if not provide')

    args = parser.parse_args()
    if not args.p:
        args.p = None
    simName = WriteOpenMMInput(args.top, args.crd, args.s, args.ms, args.cs, args.xyz[0], args.xyz[1], args.xyz[2], args.a, 
                               stride=args.stride, chkStride=args.cstride, temp=args.temp, pressure=args.p, meltTemp=args.mtemp)
    WriteJobFile(simName, args.jn) 
