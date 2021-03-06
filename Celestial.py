"""

Celestial (Class)
/////////////////

By: Tomasz Andrzejewski and Victor Steinborn


JH Computer Modeling
*****************

Current assignment:
    Project A: Solar System

*****************

This file contains the Celestial class.
The instances of this class are representations of celestial bodies in the solar system.
This class is an extension of the Particle3D class.
While Particle3D focused on kinematics, Celestial will focus on astrophysical data
and the interactions between the Particle3D objects.

Assumptions and Conventions:
- All Celestial bodies are distinct (i.e. there are no duplicate particles)



"""

# Imports
# ----------------
import math
import numpy as np
from Particle3D import Particle3D as P3D
from scipy.signal import argrelextrema

# Celestial Class
# ----------------
class Celestial (object):

    # Initialise a class level list for all existing Celestial instances
    objReg=[]

    def __init__(self, P3D, orbitingAround):
        # Register the instance in the Celestial Registrar
        Celestial.objReg.append(self)

        # Assign the relevant values of the variables to the instance 
        # (some objects will not use all variables, and are there as placeholders
        
        # For the dynamics of the particles involved:
        self.P3D = P3D # Each instance has a P3D object associated with it
        self.force_History = np.array([0,0,0]) # The forces that act on the particle at each time step will be saved in this array
        # For determining which body orbits about each other:
        self.orbitingAround = orbitingAround # A string that specifies about which the body of interest orbits 
        self.centralObj = 0 # centralObj is the (programming level) object of the body about wich the body of interest orbits
        # For maxima and minima separations:
        self.orbitVec_History = np.array([0,0,0]) # Saves the position vectors of the central body, wrt to the central body
        self.orbitSeparation = np.array([]) # Saves the magnitude of the separations between central bodies and orbital bodies
        self.perhapsesIndex=[] # List that will contain on which (integer) time step the periapsides occour
        self.apoapsisIndex=[] # List that will contain on which (integer) time step the apoapsides occour
        # For periods:
        self.periodTimes = [] # List that saves the times when the particle of interest compleded a full period.
        self.periods = [] # List that stores the time of the periods for orbiting bodies
        self.angle = 0.0 # Angular displacement of orbiting body, around central body, wrt the orbiting body's initial position. 
        self.averagePeriod = 0.0 # A float that will hold the calculated mean period of the orbiting body.
        



    # Static Methods
    # ----------------

    # The totalKineticEnergy method returns the total kinetic energy of all physical existing Celestial bodies
    @staticmethod
    def totalKineticEnergy():
        totalKE = 0.0
        # Cycle through every Celestial object that has been called into existence in the program so far
        for obj in Celestial.objReg:
            # Add the kinetic energy of all the bodies to the total energy
            totalKE = totalKE + P3D.kineticEnergy(obj.P3D)
        return totalKE

    # The totalPotentialEnergy method returns the total potential energy of all physical existing Celestial bodies.
    @staticmethod
    def totalPotentialEnergy(G):
        totalPE = 0.0
        for obj in Celestial.objReg:
            # Now go through all the pair potentials and add their contribution to the total energy
            # Note the sum below is one of the type SIGMA_{j<i} ie. there is no double counting
            for obj2 in Celestial.objReg[:Celestial.objReg.index(obj)]:
                totalPE = totalPE + P3D.potential_energy(obj.P3D, obj2.P3D, G)
        return totalPE

    # The totalEnergy method returns the total energy of the system (KE and pair PE of all existing Celestial instances)
    # The gravitational constant has to be specified as an input parameter
    @staticmethod
    def totalEnergy(G):
        totalEnergy=0.0
        # Cycle through every Celestial object that has been called into existence in the program so far
        for obj in Celestial.objReg:
            # Add the kinetic energy of all the bodies to the total energy
            totalEnergy = totalEnergy + P3D.kineticEnergy(obj.P3D)

            # Now go through all the pair potentials and add their contribution to the total energy
            # Note the sum below is one of the type SIGMA_{j<i} ie. there is no double counting
            for obj2 in Celestial.objReg[:Celestial.objReg.index(obj)]:
                totalEnergy = totalEnergy + P3D.potential_energy(obj.P3D, obj2.P3D, G)
        return totalEnergy

    # The globalForceUpdate_Initialization method determines the initial net force acting on each body in existence,
    # due to all other existing bodies.
    # The gravitational constant has to be specified as an input parameter.
    @staticmethod
    def globalForceUpdate_Initialization(G):
        # Each body will be treated in turn
        for obj in Celestial.objReg:
            # Initialize a force array
            force = np.array([0.0,0.0,0.0])

            # Consider all other bodies that are different from the body of interest
            for obj2 in Celestial.objReg:
                if obj!=obj2:
                    # Do a vector sum of the individual forces acting on the body of interest due to all other bodies
                    force = force + P3D.force_vector(obj.P3D, obj2.P3D, G)

            # Update the net force acting on the body of interest
            obj.force_History = force

    # The globalForceUpdate method is similar to the globalForceUpdate_Initialization method, but this one
    # stores the updated net force acting on the body of interest in a new row of the Nx3 force_History matrix
    # This way the past, present and future forces on the bodies can be easily known.
    # Note that the convention is so that the LAST ROW of force_History is the force acting on the particle at time t+dt
    # the PENULTIMATE row of force_History is the force acting on the particle at time t.
    @staticmethod
    def globalForceUpdate(G):
        for obj in Celestial.objReg:

            force = np.array([0.0,0.0,0.0])

            for obj2 in Celestial.objReg:
                if obj!=obj2:
                    force = force + P3D.force_vector(obj.P3D, obj2.P3D, G)

            # updated force vector is stacked at the bottom row of the force_History array.
            obj.force_History = np.vstack((obj.force_History, force))

    # The global Leap Velocity method updates the velocities of all existing bodies, given a time step
    @staticmethod
    def globalLeapVelocity(dt):
        # Each particle will be treated individually
        for obj in Celestial.objReg:
            # Update the velocity of the object of interest, given its force_History and time step dt using V-Verlet
            P3D.leapVelocity(obj.P3D, dt, obj.force_History[-2], obj.force_History[-1])

    # The globalLeapPosition method updates the positions of all existing bodies, given a time step dt
    @staticmethod
    def globalLeapPosition(dt):
        # Each particle will be treated individually
        for obj in Celestial.objReg:
            # Update the position of the object of interest, given its force_History and time step dt using V-Verlet
            P3D.leapPosition(obj.P3D, dt, obj.force_History[-1])

    # The globalCOM_Correction method corrects for the velocity of the centre of mass of the system of existing bodies.
    @staticmethod
    def globalCOM_Correction():

        # Initialize the momentum of the centre of mass (COM) and the total mass of the system
        momentumCOM=np.array([0,0,0])
        totalMass=0

        # Go through every particle and add their contributions to the momentum of the COM and the total system mass
        for obj in Celestial.objReg:
            momentumCOM = momentumCOM + P3D.momentum(obj.P3D)
            totalMass = totalMass + obj.P3D.mass

        # Define the centre of mass velocity as the the momentum of the COM divided by the system mass
        velocityCOM=momentumCOM / totalMass

        # Go through every particle and subtract off the COM velocity contribution from their velocity.
        # This changes the frame of reference of our system to one who's origin coincides with the position of the COM.
        for obj in Celestial.objReg:
            obj.P3D.velocity = obj.P3D.velocity - velocityCOM

    # The globalPositionPrint_XYZ method writes positions of all the existing Celestial objects to the given file handle
    # The positions are printed in the XYZ file format, that are also used for chemical simulations.
    # The method writes the positions of all the bodies, at the time it is called.
    # It also notes the time the method was called in the comment line.
    @staticmethod
    def globalPositionPrint_XYZ(trajectory_File_Handle, t):
        # Note the number of objects that are being written to the file
        trajectory_File_Handle.write(str(len(Celestial.objReg)))
        trajectory_File_Handle.write("\n")
        # Note the time the method was called (the time is given as an input, when calling the function)
        trajectory_File_Handle.write("Time = " + str(t))
        trajectory_File_Handle.write("\n")
        # Go through every existing object and write its position to the input file handle
        for obj in Celestial.objReg:
            # Note that the str method for particle3D objects has been defined to return a string in the XYZ file format
            trajectory_File_Handle.write(str(obj.P3D))
            trajectory_File_Handle.write("\n")
            
    # Determines the object around which a specific object orbits  
    @staticmethod
    def globalCentralBody_Initialization():
        for obj in Celestial.objReg:
            for obj2 in Celestial.objReg:
                if obj != obj2 and obj.orbitingAround != 'NONE' and obj2.P3D.label == obj.orbitingAround :
                    obj.centralObj = obj2

    # Determines the initial vector position of each object wrt the central body it orbits.
    # If the body does not orbit any other body in the system its position is left unchanged.
    @staticmethod
    def globalOrbitPosVecUpdate_Initialization():
        for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE': # Select only orbiting bodies
                fromCentre = P3D.vector_position_wrt(obj.P3D, obj.centralObj.P3D)
                obj.orbitVec_History = fromCentre
                
    # The globalOrbitPosVecUpdate method is similar to the globalOrbitPosVecUpdate_Initialization method, but this one
    # stores the updated position vector of the body of interest wrt to the central body it orbits in a new row of the Nx3 force_History matrix
    # This way the past, present and future positions of the bodies can be easily known.
    # Note that the convention is so that the LAST ROW of orbitVec_History is the position of the particle at time t+dt
    # the PENULTIMATE row of orbitVec_History is the position of the particle at time t.
    @staticmethod
    def globalOrbitPosVecUpdate():
        for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE':
                fromCentre = P3D.vector_position_wrt(obj.P3D, obj.centralObj.P3D)
                # updated position vector is stacked at the bottom row of the orbitVec_History array
                obj.orbitVec_History = np.vstack((obj.orbitVec_History, fromCentre))
                
    # Determines the magnitude of the (most recent) vector position of each orbiting object wrt the central body it orbits.
    # In other words, it determines the separation of each orbiting object and the central body it orbits.
    @staticmethod
    def globalOrbitSeparationUpdate():
        for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE':
                # Appends the separation each time it calculates to a numpy array
                obj.orbitSeparation = np.append( obj.orbitSeparation, np.linalg.norm( obj.orbitVec_History[-1]))
                
    # Determines the times at which apsides (both periapsides and apoapsides) occur of each orbiting object
    @staticmethod
    def globalApsidesIndexSearch():
        for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE':
                # Periapsis (Closest)
                periapsisIndexTuple = argrelextrema(obj.orbitSeparation, np.less) # Find the (local) minima of the separations 
                obj.periapsisIndex = periapsisIndexTuple[0].tolist() # Convert the array of the minima to a list
                # Apoapsis (Furthest)
                apoapsisIndexTuple = argrelextrema(obj.orbitSeparation, np.greater) # determine the (local) maximum separations
                obj.apoapsisIndex = apoapsisIndexTuple[0].tolist() # Convert the array of the minima to a list

    # Determines times at which one full period is made and stores them in the list.
    @staticmethod   
    def globalAngle_Check_and_Update(t):
        # Return the unit vector of the vector
        def unit_vector(vector):
            norm = np.linalg.norm(vector)
            if norm==0:
                return vector
            else:
                return vector / norm
        # Return the angle in radians between vectors 'v1' and 'v2'
        def angle_between(v1, v2):
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        # For all orbiting objects compute the angle between two latest position vectors.
        # Then sum the computed angles until full angle is obtained at which point start again from a 0-angle. 
        for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE':
                deltaAngle =angle_between(obj.orbitVec_History[-1], obj.orbitVec_History[-2])
                obj.angle = obj.angle + deltaAngle
                if obj.angle >= 2*math.pi: # If an orbit is complete, save the time where the orbit is completed to periodTimes
                    obj.angle = 0.0
                    obj.periodTimes = obj.periodTimes + [t]
                    
    # Given the periodTimes (times determined by the previous method) globalPeriodCalculation calculates the difference 
    # between the consecutive times and thus returns the periods of each body. It also calculates the average period.
    @staticmethod
    def globalPeriodCalculation():
         for obj in Celestial.objReg:
            if obj.orbitingAround != 'NONE':
                if len(obj.periodTimes)!=0:
                    for i in range (0,len(obj.periodTimes)):
                        if i==0: # 1st entry to periodTimes list gives a value of orbital period
                            period = obj.periodTimes[i]
                            obj.periods = obj.periods + [period]
                        else: # For the other entries
                            period = obj.periodTimes[i] - obj.periodTimes[i-1] # Take difference between the times periods occour to obtain orbital period.
                            obj.periods = obj.periods + [period]
                        obj.averagePeriod = np.mean(obj.periods) # calculate the average period
              
