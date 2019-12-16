import numpy as np
from sys import stdout
import mdtraj
import simtk
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

N_av = 6.022140857*10**23 /mole
kB = 1.380649*10**(-23)*joules/kelvin* N_av #joules/kelvin/mol

#output name
fName = '2Mnacl_opc_Uext1_298K'

# Input Files

prmtop = AmberPrmtopFile('6432opc_245nacl.parm7')
inpcrd = AmberInpcrdFile('6432opc_245nacl.crd')

# System Configuration

nonbondedMethod = LJPME
nonbondedCutoff = 0.9*nanometers
tailCorrection = False
ewaldErrorTolerance = 0.0001
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
box_vectors = np.diag([4.65,4.65,9.29388]) * nanometer


# Integration Options

dt = 0.002*picoseconds
temperature = 298.15*kelvin
friction = 1.0/(100.*dt)
pressure = 1.0*atmospheres
barostatInterval = 25

#External potential

Uext = 1 *kB*temperature
Nperiod = 1
axis  = 2
planeLoc = 0*nanometer
atomsInExtField = ['O']
resInExtField = ['HOH','WAT']
# Simulation Options

steps = 2e7
equilibrationSteps = 500000
#platform = Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}
dcdReporter = mdtraj.reporters.DCDReporter('trajectory_{}.dcd'.format(fName), 1000)
dataReporter = StateDataReporter('log_{}.txt'.format(fName), 1000, totalSteps=steps, step=True, speed=True, progress=True, remainingTime=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

# Prepare the Simulation

print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)

#===========================
# Create external potential
#===========================
external={"U":Uext,"NPeriod":Nperiod,"axis":axis ,"planeLoc":planeLoc}
direction=['x','y','z']
ax = external["axis"]
#atomsInExtField = [elementMap[atomname]]
if external["U"] > 0.0 * kilojoules_per_mole:
        print('Creating sinusoidal external potential in the {} direction'.format(direction[axis]))
        energy_function = 'U*sin(2*pi*NPeriod*({axis}-r0)/L)'.format(axis=direction[ax])
        fExt = openmm.CustomExternalForce(energy_function)
        fExt.addGlobalParameter("U", external["U"])
        fExt.addGlobalParameter("NPeriod", external["NPeriod"])
        fExt.addGlobalParameter("pi",np.pi)
        fExt.addGlobalParameter("r0",external["planeLoc"])
        fExt.addGlobalParameter("L",box_vectors[ax][ax])

        for ia,atom in enumerate(topology.atoms()):
             if atom.name in atomsInExtField and atom.residue.name in resInExtField:
                fExt.addParticle( ia,[] )
        system.addForce(fExt)
        print('Number of atoms in fExt %i' %fExt.getNumParticles())

integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator) #, platform, platformProperties)
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
#print('getLJPMEParametersInContext')
#print(nonbondedforce.getLJPMEParametersInContext(simulation.context))

#Restart and Check point
#simulation.loadState('output298NVT_minimized.xml')
#simulation.loadCheckpoint('checkpnt_2Mnacl_opc_Uext1_298K.chk')
simulation.reporters.append(CheckpointReporter('checkpnt_{}.chk'.format(fName), 5000))

# Minimize and Equilibrate
state = simulation.context.getState(getPositions=True)
print ("Periodic box vector: {}".format(state.getPeriodicBoxVectors()))
simulation.saveState('output_warmup_{}.xml'.format(fName))

print('Performing energy minimization...')
simulation.minimizeEnergy()
state = simulation.context.getState(getPositions=True)
PDBFile.writeModel(simulation.topology, state.getPositions(), open('traj_minimized_{}.pdb'.format(fName),'w'))

print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)
simulation.saveState('output_warmup_{}.xml'.format(fName))
# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)
simulation.saveState('output_{}.xml'.format(fName))
