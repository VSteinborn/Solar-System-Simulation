# Solar-System-Simulation

Simulation of Solar System using velocity verlet time integration techniques and a classical Newtonian gravitational model.

Only gravitational interactions are considered between the simulated particles of interest.

By: Tomasz Andrzejewski and Victor Steinborn

Matriculation Numbers:
Victor Steinborn: ***REMOVED***
Tomasz Andrzejewski: ***REMOVED***

How to execute the program:
////////////////////
To execute the program, type in terminal:
python main_solarSystemSimulation <trajectory output file> <extrema output file> <periods output file> <energy output file>

where the output file names can be specified by the user in the command line.

For how long the simulation should run for, and with what time step can be determined by the user in simParameter.txt.
The particles of interest that should be simulated can be specified by the user in particle.txt

particle.txt:
///////////////////
The particle.txt file is formatted in the folowing way:

label1 mass1 x1 y1 z1 vx1 vy1 vz1 label_of_Particle_About_Which_Particle1_Orbits
label2 mass2 x2 y2 z2 vx2 vy2 vz2 label_of_Particle_About_Which_Particle2_Orbits
.
.
.

where for the ith particle:

labeli: label of the particle (ie. "Sun")
massi: mass of the particle
xi, yi, zi: x y and z position coordinates of the particle respectively
vxi, vyi, vzi: x y and z velocity components of the particle respectively
label_of_Particle_About_Which_Particlei_Orbits: string label of particle about which particle i orbits (ie. Earth for the Moon)

(!)If a particle does not orbit about another one, then write "NONE" for label_of_Particle_aboutWiich_Pariclei_Orbits

  Units:
  - for particle.txt:
      -all masses should be given in kg
      -all position coordinates should be in au (astronomical units)
      -all velocity coordinates should be in au/day
  - The simulation will run using Au, days and earth masses
  
It is assumed that the gravitational constant is (exactly):
G=6.674e-11 m^3 kg^-1 s^-2

simParameter.txt:
////////////////////
The simParameter.txt file is formatted in the folowing way:

tTotal dt

where:
tTotal: The total time the simulation should run for (in days)
dt: Time step of the simulation (in days)
