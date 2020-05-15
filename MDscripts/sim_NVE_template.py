from sys import stdout
import numpy as np
import mdtraj, simtk
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

# Input Files
topName = '__top__'
TrajFile = 'trajectory.dcd'
ThermoLog = 'log.txt'
prmtop = AmberPrmtopFile('__top__')
inpcrd = AmberInpcrdFile('__crd__')
state2load = __state__
check2load = __check__

# System Configuration
nonbondedMethod = LJPME
nonbondedCutoff = 1.0*nanometers
tailCorrection = True
ewaldErrorTolerance = 0.0001
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
box_vectors = np.diag([__x__,__y__,__z__]) * nanometer


#Integration Options
dt = 0.002*picoseconds
temperature = __temp__*kelvin
meltTemperature = 400.0*kelvin
friction = 1.0/(100.*dt)
pressure = None 
barostatInterval = 25
checkStride = __checkStride__
trajStride = __trajStride__

#Simulation Options
steps = __steps__
equilibrationSteps = __eqSteps__
Etar = __Etar__
platformName = '__platformName__'
platform = Platform.getPlatformByName(platformName)
platformProperties = {'Precision': 'mixed'}

dcdReporter = mdtraj.reporters.DCDReporter(TrajFile, trajStride)
#dcdReporter = DCDReporter('trajectory.dcd', 5000)'
dataReporter = StateDataReporter(ThermoLog, trajStride, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

dcdReporterEquil = mdtraj.reporters.DCDReporter('trajectory_equil.dcd', trajStride)
dataReporterEquil = StateDataReporter('log_equil.txt', trajStride, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='	')


# Prepare the Simulation
print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions

def MakeSimulationNVE(positions,velocities=None, temperature = temperature):
    system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)                         
    topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())                                
    integrator = VerletIntegrator(dt)                                                                    
    integrator.setConstraintTolerance(constraintTolerance)                                               
    if platformName in ['CUDA', 'OpenCL']:                                                               
        simulation = Simulation(topology, system, integrator, platform, platformProperties)              
    else:
        simulation = Simulation(topology, system, integrator, platform)                                  
    simulation.context.setPositions(positions)
    if velocities is not None: 
        simulation.context.setVelocities(velocities)
    else: 
        simulation.context.setVelocitiesToTemperature(temperature)                                                       
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

print('===Running NVE simulation===')
simulation = MakeSimulationNVE(positions)
if state2load:
    simulation.loadState(state2load)
elif check2load:
    simulation.loadCheckpoint(check2load)
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
PE = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoules/mole)
KE = simulation.context.getState(getEnergy=True).getKineticEnergy().value_in_unit(kilojoules/mole)

#rescaling velocities
KEtar = Etar - PE
print('target KE {} kJ/mol'.format(KEtar))
print('current KE {} kJ/mol'.format(KE))
lambd = np.sqrt(KEtar/KE)
velocities *= lambd
simulation = MakeSimulationNVE(positions,velocities)

#test E
PE = simulation.context.getState(getEnergy=True).getPotentialEnergy()
KE = simulation.context.getState(getEnergy=True).getKineticEnergy()
E = PE + KE
print('E after rescaling velocities is off from E target by %2.3f percent'%(np.abs((E.value_in_unit(kilojoules/mole)-Etar)/Etar)*100.))
#print('E after rescaling velocities : {}'.format(E))
#print('E target : {} kJ/mol'.format(Etar))
simulation.saveState('output_ini.xml')

#Equilibration
print('Equilibrating...')
simulation.reporters.append(dcdReporterEquil)
simulation.reporters.append(dataReporterEquil) 
simulation.reporters.append(CheckpointReporter('checkpnt.chk',checkStride))                                     
simulation.currentStep = 0
simulation.step(equilibrationSteps)
simulation.saveState('output_warmup.xml') 
#Production
print('Simulating...')
simulation.reporters = []
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(CheckpointReporter('checkpnt.chk',checkStride))
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output.xml')

