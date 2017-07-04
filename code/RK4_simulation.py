# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 15:55:16 2017

@author: Dan
"""

import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from astropy.io import ascii


class point:
    def __init__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z

class body:
    def __init__(self, location, mass, velocity, name = ""):
        self.location = location
        self.mass = mass
        self.velocity = velocity
        self.name = name
        
def partial_step(point1, point2, time_step):
    ret = point(0,0,0)
    ret.x = point1.x + point2.x * time_step
    ret.y = point1.y + point2.y * time_step
    ret.z = point1.z + point2.z * time_step
    return ret
 
def potential_z_force(x,y,z):
    
    Sun_mass = 1.99e30
    pc_to_m = 3.086e16
    G = 6.67408e-11 #m3 kg-1 s-2


    M_bulge = 1.40e10*Sun_mass
    a_bulge = 0.0
    b_bulge = 350*pc_to_m
    
    M_disk = 7.91e10*Sun_mass
    a_disk = 3500*pc_to_m
    b_disk = 250*pc_to_m
    
    M_halo = 6.98e11*Sun_mass
    a_halo = 0.0
    b_halo = 24000*pc_to_m
    
    K_z_bulge = (G * z * (M_bulge * (a_bulge + np.sqrt(np.power(b_bulge,2) + np.power(z,2)) ) )) / (np.power((np.square(x) + np.square(y) + np.power(a_bulge + np.sqrt( np.power(b_bulge,2) + np.power(z,2) ),2)),(3/2)) * np.sqrt( np.power(b_bulge,2) + np.power(z,2)))

    K_z_disk = (G * z * (M_disk * (a_disk + np.sqrt(np.power(b_disk,2) + np.power(z,2)) ) )) / (np.power((np.square(x) + np.square(y) + np.power(a_disk + np.sqrt( np.power(b_disk,2) + np.power(z,2) ),2)),(3/2)) * np.sqrt( np.power(b_disk,2) + np.power(z,2)))

    K_z_halo = (G * z * (M_halo * (a_halo + np.sqrt(np.power(b_halo,2) + np.power(z,2)) ) )) / (np.power((np.square(x) + np.square(y) + np.power(a_halo + np.sqrt( np.power(b_halo,2) + np.power(z,2) ),2)),(3/2)) * np.sqrt( np.power(b_halo,2) + np.power(z,2)))

    return (-(K_z_bulge + K_z_disk + K_z_halo))

def potential_x_force(x,y,z):

    Sun_mass = 1.99e30
    pc_to_m = 3.086e16
    G = 6.67408e-11 #m3 kg-1 s-2


    M_bulge = 1.40e10*Sun_mass
    b_bulge = 350*pc_to_m
    
    M_disk = 7.91e10*Sun_mass
    a_disk = 3500*pc_to_m
    b_disk = 250*pc_to_m
    
    M_halo = 6.98e11*Sun_mass
    b_halo = 24000*pc_to_m
    
    K_x_disk = (G*M_disk*x) / np.power(np.square(a_disk + np.sqrt(np.square(b_disk) + np.square(z))) + np.square(x) + np.square(y),(3/2))
    K_x_bulge = (G*M_bulge*x) / np.power(np.square(b_bulge) + np.square(x) + np.square(y) + np.square(z),(3/2))
    k_x_halo = (G*M_halo*x) / np.power(np.square(b_halo) + np.square(x) + np.square(y) + np.square(z),(3/2))

    return(-(K_x_disk + K_x_bulge + k_x_halo))


def potential_y_force(x,y,z):

    Sun_mass = 1.99e30
    pc_to_m = 3.086e16
    G = 6.67408e-11 #m3 kg-1 s-2


    M_bulge = 1.40e10*Sun_mass
    b_bulge = 350*pc_to_m
    
    M_disk = 7.91e10*Sun_mass
    a_disk = 3500*pc_to_m
    b_disk = 250*pc_to_m
    
    M_halo = 6.98e11*Sun_mass
    b_halo = 24000*pc_to_m
    
    K_y_disk = (G*M_disk*y) / np.power(np.square(a_disk + np.sqrt(np.square(b_disk) + np.square(z))) + np.square(x) + np.square(y),(3/2))
    K_y_bulge = (G*M_bulge*y) / np.power(np.square(b_bulge) + np.square(x) + np.square(y) + np.square(z),(3/2))
    k_y_halo = (G*M_halo*y) / np.power(np.square(b_halo) + np.square(x) + np.square(y) + np.square(z),(3/2))

    return(-(K_y_disk + K_y_bulge + k_y_halo))


def calculate_single_body_acceleration(bodies, body_index,time_step):
    acceleration = point(0,0,0)
    target_body = bodies[body_index]

    k1 = point (0,0,0)
    k2 = point (0,0,0)
    k3 = point (0,0,0)
    k4 = point (0,0,0)
    tmp_loc = point (0,0,0)
    tmp_vel = point (0,0,0)

    #k1 - regular Euler acceleration
    
    k1.x = potential_x_force(target_body.location.x,target_body.location.y,target_body.location.z)
    k1.y = potential_y_force(target_body.location.x,target_body.location.y,target_body.location.z)
    k1.z = potential_z_force(target_body.location.x,target_body.location.y,target_body.location.z)

    #k2 - acceleration 0.5 timesteps in the future based on k1 acceleration value
    
    tmp_vel = partial_step(target_body.velocity, k1, 0.5)
    tmp_loc = partial_step(target_body.location, tmp_vel, 0.5 * time_step)
    
    k2.x = potential_x_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k2.y = potential_y_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k2.z = potential_z_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    
    #k3 acceleration 0.5 timesteps in the future using k2 acceleration
    
    tmp_vel = partial_step(target_body.velocity, k2, 0.5)
    tmp_loc = partial_step(target_body.location, tmp_vel, 0.5 * time_step)
    
    k3.x = potential_x_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k3.y = potential_y_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k3.z = potential_z_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)

    #k4 - location 1 timestep in the future using k3 acceleration
    
    tmp_vel = partial_step(target_body.velocity, k3, 1)
    tmp_loc = partial_step(target_body.location, tmp_vel, time_step)
    
    k4.x = potential_x_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k4.y = potential_y_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)
    k4.z = potential_z_force(tmp_loc.x,tmp_loc.y,tmp_loc.z)


    acceleration.x += (k1.x + k2.x * 2 + k3.x * 2 + k4.x) / 6;
    acceleration.y += (k1.y + k2.y * 2 + k3.y * 2 + k4.y) / 6;
    acceleration.z += (k1.z + k2.z * 2 + k3.z * 2 + k4.z) / 6;
    
    return acceleration

def compute_velocity(bodies, time_step):
    for body_index, target_body in enumerate(bodies):
        
        acceleration = calculate_single_body_acceleration(bodies, body_index,time_step)
        target_body.velocity.x += acceleration.x * time_step
        target_body.velocity.y += acceleration.y * time_step
        target_body.velocity.z += acceleration.z * time_step 


def update_location(bodies, time_step):
    for target_body in bodies:
        target_body.location.x += target_body.velocity.x * time_step
        target_body.location.y += target_body.velocity.y * time_step
        target_body.location.z += target_body.velocity.z * time_step

def compute_gravity_step(bodies, time_step):
    compute_velocity(bodies, time_step = time_step)
    update_location(bodies, time_step = time_step)

def plot_output(bodies, outfile = None):
    fig = plot.figure()
    colours = ['r','b','g','y','m','c']
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_range = 0
    for current_body in bodies: 
        max_dim = max(max(current_body["x"]),max(current_body["y"]),max(current_body["z"]))
        if max_dim > max_range:
            max_range = max_dim
        ax.plot(current_body["x"], current_body["y"], current_body["z"], c = random.choice(colours), label = current_body["name"])        
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    #ax.set_xlim([-max_range,max_range])    
    #ax.set_ylim([-max_range,max_range])
    #ax.set_zlim([-max_range,max_range])
    ax.legend()        

    if outfile:
        plot.savefig(outfile)
    else:
        plot.show()

def run_simulation(bodies, names = None, time_step = 1, number_of_steps = 10000, report_freq = 100):
    #create output container for each body
    body_locations_hist = []
    for current_body in bodies:
        body_locations_hist.append({"x":[], "y":[], "z":[], "name":current_body.name})
        
    for i in range(0,number_of_steps):
        compute_gravity_step(bodies, time_step)    
        if i % report_freq == 0:
            print(number_of_steps - i,"steps to go")
            for index, body_location in enumerate(body_locations_hist):
                body_location["x"].append(bodies[index].location.x) #bodies[0] is the LSR
                body_location["y"].append(bodies[index].location.y)           
                body_location["z"].append(bodies[index].location.z)       
    return body_locations_hist        
              
#cluster data (location (m), mass (kg), velocity (m/s)
pc_to_m = 3.086e16

LSR = {"location":point(8000.0*pc_to_m,0.0*pc_to_m,0.0*pc_to_m),"mass":0.0, "velocity":point(0.0e3,229.76e3,0.0e3)}   # this rotational speed comes from the fact we know the Sun is moving around galaxy at 242 km/s (C. A. L. Bailer-Jones, 2014) and that the Sun is moving 12.24 km/s faster than LSR (Ralph Schonrich, James Binney and Walter Dehnen 2009)
Sun = {"location":point(8000.0*pc_to_m,0.0*pc_to_m,10.0*pc_to_m),"mass":0.0, "velocity":point(-11.1e3,242e3,7.25e3)}

USco = {"location":point(7866.9*pc_to_m,-20.4*pc_to_m,59.5*pc_to_m),"mass":0.0, "velocity":point(-4.4e3,226.0e3,-0.75e3)} # these vels are then caluclated based on their values realtive to Suns motion in V
UCL = {"location":point(7879.6*pc_to_m,-65.5*pc_to_m,48.3*pc_to_m),"mass":0.0, "velocity":point(-4.3e3,222.7e3,1.55e3)}   
LCC = {"location":point(7934.0*pc_to_m,-98.5*pc_to_m,28.7*pc_to_m),"mass":0.0, "velocity":point(-2.9e3,223.4e3,0.85e3)}   
IC_2602 = {"location":point(7950.8*pc_to_m,-139.0*pc_to_m,-2.5*pc_to_m),"mass":0.0, "velocity":point(-1.9e3,223.5e3,6.95e3)}              
IC_2391 = {"location":point(7998.6*pc_to_m,-171.2*pc_to_m,-10.8*pc_to_m),"mass":0.0,"velocity":point(12.9e3,229.5e3,2.15e3)}             
NGC_2451A = {"location":point(8056.5*pc_to_m,-175.9*pc_to_m,-13.5*pc_to_m),"mass":0.0,"velocity":point(14.2e3,230.2e3,-4.95e3)}            
Collinder_135 = {"location":point(8087.6*pc_to_m,-222.6*pc_to_m,-37.5*pc_to_m),"mass":0.0,"velocity":point(5.5e3,231.6e3,-5.15e3)}   

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    #build list of clusters in the simulation
    bodies = [
        body( location = LSR["location"], mass = LSR["mass"], velocity = LSR["velocity"], name = "LSR"),
        body( location = Sun["location"], mass = Sun["mass"], velocity = Sun["velocity"], name = "Sun"),
        body( location = USco["location"], mass = USco["mass"], velocity = USco["velocity"], name = "USco"),
        body( location = UCL["location"], mass = UCL["mass"], velocity = UCL["velocity"], name = "UCL"),
        body( location = LCC["location"], mass = LCC["mass"], velocity = LCC["velocity"], name = "LCC"),
        body( location = IC_2602["location"], mass = IC_2602["mass"], velocity = IC_2602["velocity"], name = "IC_2602"),
        body( location = IC_2391["location"], mass = IC_2391["mass"], velocity = IC_2391["velocity"], name = "IC_2391"),
        body( location = Collinder_135["location"], mass = Collinder_135["mass"], velocity = Collinder_135["velocity"], name = "Collinder_135"),
        body( location = NGC_2451A["location"], mass = NGC_2451A["mass"], velocity = NGC_2451A["velocity"], name = "NGC_2451A")

        ]
    # time step is the length of the step in seconds, number of steps is the number of times this step length is repeated
    
    motions = run_simulation(bodies, time_step = -63080000000.0, number_of_steps = 10000, report_freq = 1)
    
    for index, body in enumerate(motions):  #writes a CSV file containing xyz positions for each cluster separately
        ascii.write([motions[index].get("x"),motions[index].get("y"),motions[index].get("z")], 'simulation_output_%s.csv'%body["name"], names=['x','y','z'], format='csv', fast_writer=False)  
        
    plot_output(motions, outfile = 'orbits.png')