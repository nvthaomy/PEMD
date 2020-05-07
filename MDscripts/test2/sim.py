from sys import stdout
import numpy as np
import mdtraj, simtk
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

# Input Files
topName = 'mixture.parm7'
TrajFile = 'trajectory.dcd'
ThermoLog = 'log.txt'
prmtop = AmberPrmtopFile('mixture.parm7')
inpcrd = AmberInpcrdFile('mixture.crd')

# System Configuration
nonbondedMethod = LJPME
nonbondedCutoff = 1.0*nanometers
tailCorrection = True
ewaldErrorTolerance = 0.0001
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
box_vectors = np.diag([5.0,5.0,5.0]) * nanometer


#Integration Options
dt = 0.002*picoseconds
temperature = 298.15*kelvin
meltTemperature = 400.0*kelvin
friction = 1.0/(100.*dt)
pressure = None
barostatInterval = 25

#Simulation Options
meltingSteps = 100
steps = 100
coolingSteps = 100
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}

dcdReporter = mdtraj.reporters.DCDReporter(TrajFile, 10)
#dcdReporter = DCDReporter('trajectory.dcd', 10)'
dataReporter = StateDataReporter(ThermoLog, 10, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

dcdReporterCooling = mdtraj.reporters.DCDReporter('trajectory_cool.dcd', 10)
dataReporterCooling = StateDataReporter('log_cool.txt', 10, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='	')

dcdReporterMelting = mdtraj.reporters.DCDReporter('trajectory_melt.dcd', 10)
dataReporterMelting = StateDataReporter('log_melt.txt', 10, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='	')

# Prepare the Simulation
print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions

def MakeSimulation(temperature):
    system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
    
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
#to load state: simulation.loadState('output.xml')

if pressure == None:
    print('===Running NVT simulation===')
else:
    print('===Running NPT simulation===')
#Melt
simulation = MakeSimulation(meltTemperature)
simulation.reporters.append(CheckpointReporter('checkpntMelt.chk', 10))
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
simulation = MakeSimulation(temperature)
simulation.reporters.append(CheckpointReporter('checkpnt.chk',10))
simulation.loadCheckpoint('checkpntMelt.chk')
print('Cool down...')
simulation.reporters.append(dcdReporterCooling)
simulation.reporters.append(dataReporterCooling)
simulation.step(coolingSteps)
simulation.saveState('output_cooled.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
#PDBFile.writeFile(simulation.topology, positions, open('cooled.pdb','w'))

#Production
print('Simulating...')
simulation.reporters = []
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
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
