"""

Particle3D (Class)
/////////////////

By: Tomasz Andrzejewski and Victor Steinborn


JH Computer Modeling
*****************

Current assignment:
    Project A: Solar System

*****************

Adapted from:
    Exercise 3
-----------------

This file contains the Particle3D class.
The instances of this class are representations of classical point like particles.
Some methods of this class give expressions that can be used to  aproximate their interactions and dynamics,
using time integration. Other methods can be used to create particle instances from files, or give information
on some of their various physical properties.

Assumptions and Conventions:
-It is assumed all numerical values that are passed though this class are floats.
-Particle3D object instances may be referred to as "particles".
-A single Particle3D object instance may be referred to as a "particle".
-We will be using a cartesian representation to describe the vectors associated with the particles.

Exercise 3:
    26.10.2016

Edited for Solar System Simulation:
    8.02.2017
"""

# Imports
# ----------------
import math
import numpy as np


# Particle3D Class
# ----------------
class Particle3D(object):
    # The _init_() command initializes a Particle3D instance
    # The particle will have a:
    # string label, a float mass, a 1x3 position array (vector) and a 1x3 velocity array.
    def __init__(self, label, mass, position, velocity):
        self.label = label
        self.mass = mass
        self.position = position
        self.velocity = velocity

    # _str_() gives information on a particle's position.
    # Returns: a string containing the particle's label and its position array.
    # This output is formatted for XYZ files
    def __str__(self):
        return str(self.label) + " " + str(self.position[0]) + " " + str(self.position[1])+ " " + str(self.position[2])

    # Instance Methods
    # ************************

    # kineticEnergy
    # Returns: The classical kinetic energy of a Particle3D object, as a float.
    def kineticEnergy(self):
        return 0.5 * self.mass * (self.velocity[0] ** 2 + self.velocity[1] ** 2 + self.velocity[2] ** 2)

    # momentum
    # Returns: The classical momentum of a particle3D object as a float.
    def momentum(self):
        return float(self.mass) * self.velocity

    # Time Integration Methods
    # -------------------------

    # leapVelocity
    # Updates the velocity of the particle using the Velocity Verlet method
    # Arguments- Floats: dt; 1x3 Arrays: force, new_force
    def leapVelocity(self, dt, force, new_force):
        self.velocity = self.velocity + 0.5 * dt * ((force + new_force) / (self.mass))

    # leapPosition
    # Updates the position of the particle to second order using the Velocity Verlet method.
    # Arguments- Floats: dt; 1x3 Arrays: force
    def leapPosition(self, dt, force):
        self.position = self.position + dt * self.velocity + 0.5 * dt ** 2 * force / self.mass

    # Static Methods
    # **************************

    # vector_position_wrt
    # Returns the vector position of a particle (p1)  with respect to another particle (p2).
    # 'wrt' stands for 'with respect to'
    # Arguments- Particle3D objects: p1, p2
    @staticmethod
    def vector_position_wrt(p1, p2):
        return np.array(
            [p1.position[0] - p2.position[0], p1.position[1] - p2.position[1], p1.position[2] - p2.position[2]])

    # separation
    # Returns the magnitude of the vector separation between two particles (p1 and p2)
    # Arguments- Particle3D objects: p1, p2
    @staticmethod
    def separation(p1, p2):
        vecPosition = Particle3D.vector_position_wrt(p1, p2)
        return math.sqrt(vecPosition[0] ** 2 + vecPosition[1] ** 2 + vecPosition[2] ** 2)

    # potential_energy
    # Returns the gravitational potential energy arising from the gravitational interaction between
    # two particles (p1 and p2)
    # Caution: p1 and p2 must be at different spatial positions, or this method will divide by zero.
    # Arguments- Particle3D objects: p1, p2; floats: G
    @staticmethod
    def potential_energy(p1, p2, G):
        return -1.0 * G * p1.mass * p2.mass / (Particle3D.separation(p1, p2))

    # force_vector
    # Returns the gravitational force on a particle (p1) due to another particle (p2)
    # Caution: p1 and p2 must be at different spatial positions, or this method will divide by zero.
    # Arguments- Particle3D objects: p1, p2; floats: G
    @staticmethod
    def force_vector(p1, p2, G):
        return (G * p1.mass * p2.mass / Particle3D.separation(p1, p2) ** 3) * Particle3D.vector_position_wrt(p2, p1)

    # from_line
    # Static method to create a particle from a file entry.
    # It takes a line of a file handle as argument and returns a Particle 3D object
    # The file line should be formatted in the following way:
    # label mass x1 x2 x3 v1 v2 v3 (other information may follow)
    @staticmethod
    def from_line(file_line):
        # Split the line into its components
        components = file_line.split()
        # We will create vectors that will represent the particle's position and velocity.
        # We will also ensure that the numerical values of the Particle3D instance are floats.
        position = np.array([float(components[2]), float(components[3]), float(components[4])])
        velocity = np.array([float(components[5]), float(components[6]), float(components[7])])
        return Particle3D(components[0], float(components[1]), position, velocity)
