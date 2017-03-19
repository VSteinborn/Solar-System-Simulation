"""
Solar System Simulation
////////////////////////

By: Tomasz Andrzejewski and Victor Steinborn
Matriculation numbers:
Victor Steinborn: ***REMOVED***
Tomasz Andrzejewski: ***REMOVED***


JH Computer Modeling
*****************

Assignment:
    Project A: Solar System

*****************
Log:
    10.02.2017
-VS: Created file
-----------------
This file contains the main script that is used to simulate the Solar System, using Velocity Verlet time integration techniques.


This file contains the main program that simulates the evolution of the system consisting of an arbitrary number of particles using the velocity Verlet time integration
algorithm, as a function of time in an external potential,
with given initial conditions.

To run the program, names of 4 output files must be written to the command line.
These 4 output files hold respectively data for:
trajectory
energy
apoapsides and periapsides output extremaCheck() method.
orbital period lengths periodCalculations() method

Current Assumptions:
-----------------
-The distance between two bodies will never be less than or equal to the sum of their radii. (point particle treatment)
-The orbital period is the amount of time it takes for the particle to complete an orbit about the central particle.

"""

# Imports
# ---------------

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D as P3D
from Celestial import Celestial as CEL

# Global Constants
# ---------------

# Universal Gravitational Constant G (m^3 kg^-1 s^-2)
G_Experimental = 6.674e-11
# Time step size dt (days)
dt = 1.0
# Total time of simulation (days) (It is assumed that dt is a factor of tTotal)
tTotal = 6*4*(365.25)
# Conversions:
Meters_In_AU = 149597870700.0
Kilograms_In_EarthMass = 5.972e24
Seconds_In_Year = 31557600.0 # We will be using Julian Years
Seconds_In_Day = 86400.0
# The value of G that will be used in the simulation (Simulation Units)
G = G_Experimental * Kilograms_In_EarthMass * Seconds_In_Day**2  / (Meters_In_AU**3 )

# Conversions
def metersToAU(m): return (m / Meters_In_AU)
def kilogramsToEarthMass(kg): return (kg / Kilograms_In_EarthMass)
def simulationEnergyToGJ(energy): return  energy * Kilograms_In_EarthMass* Meters_In_AU**2 / (Seconds_In_Day**2)*1.0e-9
def daysToYears(days): return days * Seconds_In_Day / Seconds_In_Year

# User Input and Program Execution
# ---------------

# Read name of 4 output files from command line: trajectory, energy, extrema and periods
if len(sys.argv) != 4:
    print    "Wrong number of arguments."
    print    "Usage: " + sys.argv[
        0] + " <trajectory output file>" + " <extrema output file>" + " <periods output file>"
    quit()
else:
    trajectory_File_Name = sys.argv[1]
    periAndApo_File_Name = sys.argv[2]
    periods_File_Name = sys.argv[3]
# Open output file for writing
trajectory_File_Handle = open(trajectory_File_Name, "w")
periAndApo_File_Handle = open(periAndApo_File_Name, "w")
periods_File_Handle = open(periods_File_Name, "w")

# Reading Input Files and Initializing Celestial bodies
# ---------------

# Attach file handle to input file
file_handle = open("particles.txt", "r")
# Return number of lines in the input file which represents the number of particles in the simulation
particleNumber = len(file_handle.readlines())



# Augment file handle with a list of particles so the particles' masses are in units of Earth masses
augmentedParticle_File_Handle = open("augmentedParticles.txt", "w")
with open("particles.txt", "r") as file_handle:
    for line in file_handle:
        comp = line.split()
        augmentedParticle_File_Handle.write(comp[0] +" "+ str(kilogramsToEarthMass(float(comp[1])))+" ")
        for i in range(2,5):
            augmentedParticle_File_Handle.write(comp[i]+" ")
        for i in range(5,8):
            augmentedParticle_File_Handle.write(comp[i]+" ")
        augmentedParticle_File_Handle.write(str(comp[8]))
        augmentedParticle_File_Handle.write("\n")
augmentedParticle_File_Handle.close()



# Initialize list that stores Celestial objects
celestialList = []

# Read the particle data for the all the particles from file and save it in celestialList
with open("augmentedParticles.txt", "r") as file_handle:
    for line in file_handle:
        components = line.split()
        celestialList = celestialList + [ CEL(P3D.from_line(line), components[8]) ]

# Initialize all the forces acting on each body
CEL.globalForceUpdate_Initialization(G)

# Change the frame of reference of the system to the centre-of-mass frame by correcting for the motion of the centre of mass
CEL.globalCOM_Correction()

# Associated with each orbiting body is a central body, about which the orbiting body will revolve about.
# This method pairs the orbiting body with its corresponding central body, by giving each orbiting body a "partnerObj"
# The partnerObj is the Celestial object of the central body.
# This makes it easier to reference the body about which another orbits.
CEL.globalCentralBody_Initialization()

# Find all the initial vectors from the central bodies to the orbiting bodies.
CEL.globalOrbitPosVecUpdate_Initialization()

# Set up data lists
# ---------------

# The following lists will store the potential energy, kinetic energy, and total energy of the system at various times
potentialEnergyList = []
kineticEnergyList = []
energyList = []

# Main simulation loop
# ---------------

# Simulation runs from t=0 to t=tTotal (Including boundaries) in time steps of dt
timeArray=np.arange(0, tTotal + dt, dt) # returns an array of evenly spaced time steps: (0, dt,2*dt,...tTotal)
for t in timeArray:
    # Data
    # ---------------
    CEL.globalPositionPrint_XYZ(trajectory_File_Handle, t)
    potentialEnergyList.append( CEL.totalPotentialEnergy(G))
    kineticEnergyList.append( CEL.totalKineticEnergy())
    energyList.append( potentialEnergyList[-1] + kineticEnergyList[-1])
    CEL.globalOrbitPosVecUpdate()
    CEL.globalOrbitSeparationUpdate()

    # Dynamics
    # ---------------
    CEL.globalLeapPosition(dt)
    CEL.globalForceUpdate(G)
    CEL.globalLeapVelocity(dt)
    #CEL.globalAngle_Check_and_Update(t)

# Data presentation
# ---------------

# Convert Energies to SI
energyInGJList = [simulationEnergyToGJ(x) for x in energyList]
keInGJList = [simulationEnergyToGJ(x) for x in kineticEnergyList]
peInGJList = [simulationEnergyToGJ(x) for x in potentialEnergyList]
# Convert Times to (Julian) Years
timeInYears = [daysToYears(t) for t in np.arange(0, tTotal +dt , dt)]

# Plot energies (KE,PE and total) of the system (in SI units)
pyplot.plot(timeInYears, energyInGJList)
pyplot.plot(timeInYears, keInGJList)
pyplot.plot(timeInYears, peInGJList)

pyplot.title("The total energy of the Solar System as a function of time")
pyplot.xlabel("Time (Years)")
pyplot.ylabel("Total energy (GJ)")
pyplot.axhline(y=0, color='b', linestyle='dashed')
# Plot additional horizontal line line showing the initial value of the energy.
pyplot.axhline(y=energyList[0], color='r', linestyle='dashed')

pyplot.legend(['Total Energy', 'Kinetic Energy', 'Potential Energy', 'E=0', 'Initial Total Energy'], loc='upper right')
#               pyplot.savefig('SolarSystem_energy.png')
pyplot.figure()

# Plot total energy of the system (in SI units)

pyplot.plot(timeInYears, energyInGJList)
pyplot.title("The total energy of the Solar System as a function of time")
pyplot.xlabel("Time (Years)")
pyplot.ylabel("Total energy (GJ)")
# Plot additional horizontal line line showing the initial value of the energy.
pyplot.figure()

# Plot orbital separations
CEL.globalApoAndPeriapsesIndexSearch()
for obj in CEL.objReg:
    if obj.orbitingAround != 'NONE':
        pyplot.plot(timeInYears, obj.orbitSeparation)
        pyplot.title("Separation of "+ obj.P3D.label +" from the " + obj.orbitingAround + " as a function of time")
        pyplot.xlabel("Time (Years)")
        pyplot.ylabel("Distance (AU)")

        periAndApo_File_Handle.write("Periapsis Times (days): " + obj.P3D.label)
        periAndApo_File_Handle.write('\n')
        for i in obj.perhapsesIndex:
            periAndApo_File_Handle.write(str(timeArray[i]) + "\n")
        periAndApo_File_Handle.write('\n')
        
        periods_File_Handle.write("Periods (Earth days): " + obj.P3D.label)
        periods_File_Handle.write('\n')
        for i in obj.periodTimes:
            periods_File_Handle.write(str(timeArray[i]) + "\n")
        periods_File_Handle.write('\n')

        periAndApo_File_Handle.write("Apoapsis Times (days): " + obj.P3D.label)
        periAndApo_File_Handle.write('\n')
        for i in obj.apoapsisIndex:
            periAndApo_File_Handle.write(str(timeArray[i]) + "\n")
        periAndApo_File_Handle.write('\n')
        pyplot.figure()
pyplot.show()
