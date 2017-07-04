# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:30:24 2017

@author: Dan
"""

from astropy.table import Table
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy import units as u
import numpy as np
from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mp
import moviepy.editor as mpy

USco = Table.read('simulation_output_USco.csv',format='csv')    
UCL = Table.read('simulation_output_UCL.csv',format='csv')
LCC = Table.read('simulation_output_LCC.csv',format='csv')
IC_2602 = Table.read('simulation_output_IC_2602.csv',format='csv')
IC_2391 = Table.read('simulation_output_IC_2391.csv',format='csv')
NGC_2451A = Table.read('simulation_output_NGC_2451A.csv',format='csv')
Collinder_135 = Table.read('simulation_output_Collinder_135.csv',format='csv')

Sun = Table.read('simulation_output_Sun.csv',format='csv')
LSR = Table.read('simulation_output_LSR.csv',format='csv')

m_to_pc = (u.meter).to(u.parsec)
seconds_to_myr = ((u.second).to(u.year))/1e6

def plot_around_gal_center_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    plt.plot(USco['x']*m_to_pc,USco['y']*m_to_pc,label='USco')
    #plt.plot(UCL['x']*m_to_pc,UCL['y']*m_to_pc,label='UCL')
    #plt.plot(LCC['x']*m_to_pc,LCC['y']*m_to_pc,label='LCC')
    #plt.plot(IC_2602['x']*m_to_pc,IC_2602['y']*m_to_pc,label='IC 2602')
    #plt.plot(IC_2391['x']*m_to_pc,IC_2391['y']*m_to_pc,label='IC 2391')
    #plt.plot(NGC_2451A['x']*m_to_pc,NGC_2451A['y']*m_to_pc,label='NGC 2451A')
    #plt.plot(Collinder_135['x']*m_to_pc,Collinder_135['y']*m_to_pc,label='Collinder 135')

    #plt.plot(Sun['x']*m_to_pc,Sun['y']*m_to_pc,label='Sun', )
    plt.plot(LSR['x']*m_to_pc,LSR['y']*m_to_pc,label='LSR', )

    plt.scatter(USco['x'][0]*m_to_pc,USco['y'][0]*m_to_pc,s=40)
    #plt.scatter(UCL['x'][0]*m_to_pc,UCL['y'][0]*m_to_pc,s=40)
    #plt.scatter(LCC['x'][0]*m_to_pc,LCC['y'][0]*m_to_pc,s=40)
    #plt.scatter(IC_2602['x'][0]*m_to_pc,IC_2602['y'][0]*m_to_pc,s=40)
    #plt.scatter(IC_2391['x'][0]*m_to_pc,IC_2391['y'][0]*m_to_pc,s=40)
    #plt.scatter(NGC_2451A['x'][0]*m_to_pc,NGC_2451A['y'][0]*m_to_pc,s=40)
    #plt.scatter(Collinder_135['x'][0]*m_to_pc,Collinder_135['y'][0]*m_to_pc,s=40)  

    #plt.scatter(Sun['x'][0]*m_to_pc,Sun['y'][0]*m_to_pc,s=40)
    plt.scatter(LSR['x'][0]*m_to_pc,LSR['y'][0]*m_to_pc,s=40)
    
    
    plt.xlabel('X (pc)')
    plt.ylabel('Y (pc)')
    plt.legend()
    plt.savefig('around_gal_center.jpg',format='jpg', dpi=600)
    
    plt.show()


def plot_LSR_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    plt.plot(-(USco['x'] - LSR['x'])*m_to_pc,(USco['y'] - LSR['y'])*m_to_pc,label='USco')
    plt.plot(-(UCL['x'] - LSR['x'])*m_to_pc,(UCL['y'] - LSR['y'])*m_to_pc,label='UCL')     # - x value here to flip the x axis back to +ve being towards the gal center
    plt.plot(-(LCC['x'] - LSR['x'])*m_to_pc,(LCC['y'] - LSR['y'])*m_to_pc,label='LCC')
    plt.plot(-(IC_2602['x'] - LSR['x'])*m_to_pc,(IC_2602['y'] - LSR['y'])*m_to_pc,label='IC 2602')
    plt.plot(-(IC_2391['x'] - LSR['x'])*m_to_pc,(IC_2391['y'] - LSR['y'])*m_to_pc,label='IC 2391')
    plt.plot(-(NGC_2451A['x'] - LSR['x'])*m_to_pc,(NGC_2451A['y'] - LSR['y'])*m_to_pc,label='NGC 2451A')
    plt.plot(-(Collinder_135['x'] - LSR['x'])*m_to_pc,(Collinder_135['y'] - LSR['y'])*m_to_pc,label='Collinder 135')

    plt.plot(-(Sun['x'] - LSR['x'])*m_to_pc,(Sun['y'] - LSR['y'])*m_to_pc,label='Sun')

    plt.scatter(-(USco['x'][0] - LSR['x'][0])*m_to_pc,(USco['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(UCL['x'][0] - LSR['x'][0])*m_to_pc,(UCL['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(LCC['x'][0] - LSR['x'][0])*m_to_pc,(LCC['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(IC_2602['x'][0] - LSR['x'][0])*m_to_pc,(IC_2602['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(IC_2391['x'][0] - LSR['x'][0])*m_to_pc,(IC_2391['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(NGC_2451A['x'][0] - LSR['x'][0])*m_to_pc,(NGC_2451A['y'][0] - LSR['y'][0])*m_to_pc,s=40)
    plt.scatter(-(Collinder_135['x'][0] - LSR['x'][0])*m_to_pc,(Collinder_135['y'][0] - LSR['y'][0])*m_to_pc,s=40)  

    plt.scatter(-(Sun['x'][0] - LSR['x'][0])*m_to_pc,(Sun['y'][0] - LSR['y'][0])*m_to_pc,s=40)    
    
    plt.xlabel('X (pc)')
    plt.ylabel('Y (pc)')
    plt.legend()
    plt.savefig('LSR_XY.jpg',format='jpg', dpi=600)
    
    plt.show()
    
def plot_LSR_XZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    plt.plot(-(USco['x'] - LSR['x'])*m_to_pc,(USco['z'] - LSR['z'])*m_to_pc,label='USco')
    plt.plot(-(UCL['x'] - LSR['x'])*m_to_pc,(UCL['z'] - LSR['z'])*m_to_pc,label='UCL')     # - x value here to flip the x axis back to +ve being towards the gal center
    plt.plot(-(LCC['x'] - LSR['x'])*m_to_pc,(LCC['z'] - LSR['z'])*m_to_pc,label='LCC')
    plt.plot(-(IC_2602['x'] - LSR['x'])*m_to_pc,(IC_2602['z'] - LSR['z'])*m_to_pc,label='IC 2602')
    plt.plot(-(IC_2391['x'] - LSR['x'])*m_to_pc,(IC_2391['z'] - LSR['z'])*m_to_pc,label='IC 2391')
    plt.plot(-(NGC_2451A['x'] - LSR['x'])*m_to_pc,(NGC_2451A['z'] - LSR['z'])*m_to_pc,label='NGC 2451A')
    plt.plot(-(Collinder_135['x'] - LSR['x'])*m_to_pc,(Collinder_135['z'] - LSR['z'])*m_to_pc,label='Collinder 135')

    plt.plot(-(Sun['x'] - LSR['x'])*m_to_pc,(Sun['z'] - LSR['z'])*m_to_pc,label='Sun')

    plt.scatter(-(USco['x'][0] - LSR['x'][0])*m_to_pc,(USco['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(UCL['x'][0] - LSR['x'][0])*m_to_pc,(UCL['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(LCC['x'][0] - LSR['x'][0])*m_to_pc,(LCC['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(IC_2602['x'][0] - LSR['x'][0])*m_to_pc,(IC_2602['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(IC_2391['x'][0] - LSR['x'][0])*m_to_pc,(IC_2391['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(NGC_2451A['x'][0] - LSR['x'][0])*m_to_pc,(NGC_2451A['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter(-(Collinder_135['x'][0] - LSR['x'][0])*m_to_pc,(Collinder_135['z'][0] - LSR['z'][0])*m_to_pc,s=40)  

    plt.scatter(-(Sun['x'][0] - LSR['x'][0])*m_to_pc,(Sun['z'][0] - LSR['z'][0])*m_to_pc,s=40)    
    
    plt.xlabel('X (pc)')
    plt.ylabel('Z (pc)')
    plt.legend()
    plt.savefig('LSR_XZ.jpg',format='jpg', dpi=600)
    
    plt.show()
    
def plot_LSR_YZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    plt.plot((USco['y'] - LSR['y'])*m_to_pc,(USco['z'] - LSR['z'])*m_to_pc,label='USco')
    plt.plot((UCL['y'] - LSR['y'])*m_to_pc,(UCL['z'] - LSR['z'])*m_to_pc,label='UCL')     # - x value here to flip the x axis back to +ve being towards the gal center
    plt.plot((LCC['y'] - LSR['y'])*m_to_pc,(LCC['z'] - LSR['z'])*m_to_pc,label='LCC')
    plt.plot((IC_2602['y'] - LSR['y'])*m_to_pc,(IC_2602['z'] - LSR['z'])*m_to_pc,label='IC 2602')
    plt.plot((IC_2391['y'] - LSR['y'])*m_to_pc,(IC_2391['z'] - LSR['z'])*m_to_pc,label='IC 2391')
    plt.plot((NGC_2451A['y'] - LSR['y'])*m_to_pc,(NGC_2451A['z'] - LSR['z'])*m_to_pc,label='NGC 2451A')
    plt.plot((Collinder_135['y'] - LSR['y'])*m_to_pc,(Collinder_135['z'] - LSR['z'])*m_to_pc,label='Collinder 135')

    plt.plot((Sun['y'] - LSR['y'])*m_to_pc,(Sun['z'] - LSR['z'])*m_to_pc,label='Sun')

    plt.scatter((USco['y'][0] - LSR['y'][0])*m_to_pc,(USco['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((UCL['y'][0] - LSR['y'][0])*m_to_pc,(UCL['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((LCC['y'][0] - LSR['y'][0])*m_to_pc,(LCC['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((IC_2602['y'][0] - LSR['y'][0])*m_to_pc,(IC_2602['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((IC_2391['y'][0] - LSR['y'][0])*m_to_pc,(IC_2391['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((NGC_2451A['y'][0] - LSR['y'][0])*m_to_pc,(NGC_2451A['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    plt.scatter((Collinder_135['y'][0] - LSR['y'][0])*m_to_pc,(Collinder_135['z'][0] - LSR['z'][0])*m_to_pc,s=40)  

    plt.scatter((Sun['y'][0] - LSR['y'][0])*m_to_pc,(Sun['z'][0] - LSR['z'][0])*m_to_pc,s=40)    
    
    plt.xlabel('Y (pc)')
    plt.ylabel('Z (pc)')
    plt.legend()
    plt.savefig('LSR_YZ.jpg',format='jpg', dpi=600)
    
    plt.show()
    
def plot_3D_LSR(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    ax.plot(-(USco['x'] - LSR['x'])*m_to_pc,(USco['y'] - LSR['y'])*m_to_pc,(USco['z'] - LSR['z'])*m_to_pc, label='USco',linewidth = 3)
    ax.plot(-(UCL['x'] - LSR['x'])*m_to_pc,(UCL['y'] - LSR['y'])*m_to_pc,(UCL['z'] - LSR['z'])*m_to_pc, label='UCL',linewidth = 3)    
    ax.plot(-(LCC['x'] - LSR['x'])*m_to_pc,(LCC['y'] - LSR['y'])*m_to_pc,(LCC['z'] - LSR['z'])*m_to_pc, label='LCC',linewidth = 3)
    ax.plot(-(IC_2602['x'] - LSR['x'])*m_to_pc,(IC_2602['y'] - LSR['y'])*m_to_pc,(IC_2602['z'] - LSR['z'])*m_to_pc, label='IC_2602',linewidth = 3)
    ax.plot(-(IC_2391['x'] - LSR['x'])*m_to_pc,(IC_2391['y'] - LSR['y'])*m_to_pc,(IC_2391['z'] - LSR['z'])*m_to_pc, label='IC_2391',linewidth = 3)
    ax.plot(-(NGC_2451A['x'] - LSR['x'])*m_to_pc,(NGC_2451A['y'] - LSR['y'])*m_to_pc,(NGC_2451A['z'] - LSR['z'])*m_to_pc, label='NGC_2451A',linewidth = 3)
    ax.plot(-(Collinder_135['x'] - LSR['x'])*m_to_pc,(Collinder_135['y'] - LSR['y'])*m_to_pc,(Collinder_135['z'] - LSR['z'])*m_to_pc, label='Collinder_135',linewidth = 3)
    
    ax.plot(-(Sun['x'] - LSR['x'])*m_to_pc,(Sun['y'] - LSR['y'])*m_to_pc,(Sun['z'] - LSR['z'])*m_to_pc, label='Sun',linewidth = 3)
 
    ax.scatter(-(USco['x'][0] - LSR['x'][0])*m_to_pc,(USco['y'][0] - LSR['y'][0])*m_to_pc,(USco['z'][0] - LSR['z'][0])*m_to_pc,s=40)    
    ax.scatter(-(UCL['x'][0] - LSR['x'][0])*m_to_pc,(UCL['y'][0] - LSR['y'][0])*m_to_pc,(UCL['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    ax.scatter(-(LCC['x'][0] - LSR['x'][0])*m_to_pc,(LCC['y'][0] - LSR['y'][0])*m_to_pc,(LCC['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    ax.scatter(-(IC_2602['x'][0] - LSR['x'][0])*m_to_pc,(IC_2602['y'][0] - LSR['y'][0])*m_to_pc,(IC_2602['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    ax.scatter(-(IC_2391['x'][0] - LSR['x'][0])*m_to_pc,(IC_2391['y'][0] - LSR['y'][0])*m_to_pc,(IC_2391['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    ax.scatter(-(NGC_2451A['x'][0] - LSR['x'][0])*m_to_pc,(NGC_2451A['y'][0] - LSR['y'][0])*m_to_pc,(NGC_2451A['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    ax.scatter(-(Collinder_135['x'][0] - LSR['x'][0])*m_to_pc,(Collinder_135['y'][0] - LSR['y'][0])*m_to_pc,(Collinder_135['z'][0] - LSR['z'][0])*m_to_pc,s=40)
    
    ax.scatter(-(Sun['x'][0] - LSR['x'][0])*m_to_pc,(Sun['y'][0] - LSR['y'][0])*m_to_pc,(Sun['z'][0] - LSR['z'][0])*m_to_pc,s=40)

    
    ax.set_xlim(-100,300)
    ax.set_ylim(-250, 150)
    ax.set_zlim(-200, 200)
    ax.set_xlabel('X (pc)')
    ax.set_ylabel('Y (pc)')
    ax.set_zlabel('Z (pc)')
    plt.legend()
    plt.show()
    
def Fernandez_2007_coords(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun):
    
    LSR_radius_val = []
    LSR_nu = []
    set_initial_radius = 8.5  #to be entered in kpc
    
    for i in range(0,len(LSR)):
        
        angle = np.arctan2(LSR["y"][i],LSR["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(LSR["x"][i],2) + np.power(LSR["y"][i],2))
        
        LSR_radius_val.append(((radius - set_initial_radius*3.086e+19)/3.086e+16)) 
        LSR_nu.append((angle*set_initial_radius*3.086e+19)/3.086e+16)
    
    USco_radius_val = []
    USco_nu = []
    
    
    for i in range(0,len(USco)):
        
        angle = np.arctan2(USco["y"][i],USco["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(USco["x"][i],2) + np.power(USco["y"][i],2))
        
        USco_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        USco_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    
    plt.plot(USco_nu,USco_radius_val,label = 'USco',linewidth=3)
    plt.scatter(USco_nu[0],USco_radius_val[0])
    
    UCL_radius_val = []
    UCL_nu = []
    
    for i in range(0,len(UCL)):
        
        angle = np.arctan2(UCL["y"][i],UCL["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(UCL["x"][i],2) + np.power(UCL["y"][i],2))
        
        UCL_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        UCL_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    
    plt.plot(UCL_nu,UCL_radius_val,label='UCL',linewidth=3)
    plt.scatter(UCL_nu[0],UCL_radius_val[0])
    
    LCC_radius_val = []
    LCC_nu = []
    
    for i in range(0,len(LCC)):
        
        angle = np.arctan2(LCC["y"][i],LCC["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(LCC["x"][i],2) + np.power(LCC["y"][i],2))
        
        LCC_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        LCC_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    plt.plot(LCC_nu,LCC_radius_val,label='LCC',linewidth=3)
    plt.scatter(LCC_nu[0],LCC_radius_val[0])
    
    IC_2602_radius_val = []
    IC_2602_nu = []
    
    for i in range(0,len(IC_2602)):
        
        angle = np.arctan2(IC_2602["y"][i],IC_2602["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(IC_2602["x"][i],2) + np.power(IC_2602["y"][i],2))
        
        IC_2602_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        IC_2602_nu.append((((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i]))
        
    
    plt.plot(IC_2602_nu,IC_2602_radius_val,label='IC 2602',linewidth=3)
    plt.scatter(IC_2602_nu[0],IC_2602_radius_val[0])

    IC_2391_radius_val = []
    IC_2391_nu = []
    
    for i in range(0,len(IC_2391)):
        
        angle = np.arctan2(IC_2391["y"][i],IC_2391["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(IC_2391["x"][i],2) + np.power(IC_2391["y"][i],2))
        
        IC_2391_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        IC_2391_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    
    plt.plot(IC_2391_nu,IC_2391_radius_val,label='IC 2391',linewidth=3)
    plt.scatter(IC_2391_nu[0],IC_2391_radius_val[0])
    
    
    NGC_2451A_radius_val = []
    NGC_2451A_nu = []
    
    for i in range(0,len(NGC_2451A)):
        
        angle = np.arctan2(NGC_2451A["y"][i],NGC_2451A["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(NGC_2451A["x"][i],2) + np.power(NGC_2451A["y"][i],2))
        
        NGC_2451A_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        NGC_2451A_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    
    plt.plot(NGC_2451A_nu,NGC_2451A_radius_val,label='NGC 2451A',linewidth=3)
    plt.scatter(NGC_2451A_nu[0],NGC_2451A_radius_val[0])
    
    Collinder_135_radius_val = []
    Collinder_135_nu = []
    
    for i in range(0,len(Collinder_135)):
        
        angle = np.arctan2(Collinder_135["y"][i],Collinder_135["x"][i]) # = 0 initially 
        radius = np.sqrt(np.power(Collinder_135["x"][i],2) + np.power(Collinder_135["y"][i],2))
        
        Collinder_135_radius_val.append(((radius - set_initial_radius*3.086e+19 )/3.086e+16) - LSR_radius_val[i]) 
        Collinder_135_nu.append(((angle*set_initial_radius*3.086e+19)/3.086e+16) - LSR_nu[i])
        
    
    plt.plot(Collinder_135_nu,Collinder_135_radius_val,label='Collinder 135',linewidth=3)
    plt.scatter(Collinder_135_nu[0],Collinder_135_radius_val[0])    
    
    
    
    
    #plt.ylim(-250,0)
    #plt.xlim(-150,250)
    plt.ylabel("\u03BE' (pc)")
    plt.xlabel("\u03B7' (pc)")
    plt.title('Fernández et al 2007 Comparison')
    plt.savefig('Fernández et al 2007 Comparison.jpg',format='jpg', dpi=600)

    #plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
def around_gal_center_animation(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,time_step_length):
    
    duration = 10 # time in seconds the animation lasts for
    fig_mpl, ax = plt.subplots(1,figsize=(5,5), facecolor='white')
  
    plt.scatter(USco['x'][0]*m_to_pc,USco['y'][0]*m_to_pc,s=20)
    plt.scatter(UCL['x'][0]*m_to_pc,UCL['y'][0]*m_to_pc,s=20)
    plt.scatter(LCC['x'][0]*m_to_pc,LCC['y'][0]*m_to_pc,s=20)
    plt.scatter(IC_2602['x'][0]*m_to_pc,IC_2602['y'][0]*m_to_pc,s=20)
    plt.scatter(IC_2391['x'][0]*m_to_pc,IC_2391['y'][0]*m_to_pc,s=20)
    plt.scatter(NGC_2451A['x'][0]*m_to_pc,NGC_2451A['y'][0]*m_to_pc,s=20)
    plt.scatter(Collinder_135['x'][0]*m_to_pc,Collinder_135['y'][0]*m_to_pc,s=20)  

    plt.scatter(Sun['x'][0]*m_to_pc,Sun['y'][0]*m_to_pc,s=20)
    plt.scatter(LSR['x'][0]*m_to_pc,LSR['y'][0]*m_to_pc,s=20)
    
    
    line, = ax.plot(USco['x']*m_to_pc,USco['y']*m_to_pc,label='USco',lw=2)
    line_2, = ax.plot(UCL['x']*m_to_pc,UCL['y']*m_to_pc,label='UCL',lw=2)
    line_3, = ax.plot(LCC['x']*m_to_pc,LCC['y']*m_to_pc,label='LCC',lw=2)
    line_4, = ax.plot(IC_2602['x']*m_to_pc,IC_2602['y']*m_to_pc,label='IC 2602',lw=2)
    line_5, = ax.plot(IC_2391['x']*m_to_pc,IC_2391['y']*m_to_pc,label='IC 2391',lw=2)
    line_6, = ax.plot(NGC_2451A['x']*m_to_pc,NGC_2451A['y']*m_to_pc,label='NGC 2451A',lw=2)
    line_7, = ax.plot(Collinder_135['x']*m_to_pc,Collinder_135['y']*m_to_pc,label='Collinder 135',lw=2)
    line_8, = ax.plot(LSR['x']*m_to_pc,LSR['y']*m_to_pc,label='LSR',lw=2)
    line_9, = ax.plot(Sun['x']*m_to_pc,Sun['y']*m_to_pc,label='Sun',lw=2)

    #ax.set_xlim(4250,9750)
    #ax.set_ylim(-5000,500)
    
    def make_frame_mpl(t):
        i = np.int(len(USco)*(t/duration))
        time = np.round((time_step_length*i)*seconds_to_myr,0)
        
        line.set_xdata(USco['x'][:i]*m_to_pc)
        line.set_ydata(USco['y'][:i]*m_to_pc)
        
        line_2.set_xdata(UCL['x'][:i]*m_to_pc)
        line_2.set_ydata(UCL['y'][:i]*m_to_pc)
        
        line_3.set_xdata(LCC['x'][:i]*m_to_pc)
        line_3.set_ydata(LCC['y'][:i]*m_to_pc)
        
        line_4.set_xdata(IC_2602['x'][:i]*m_to_pc)
        line_4.set_ydata(IC_2602['y'][:i]*m_to_pc)
        
        line_5.set_xdata(IC_2391['x'][:i]*m_to_pc)
        line_5.set_ydata(IC_2391['y'][:i]*m_to_pc)
        
        line_6.set_xdata(NGC_2451A['x'][:i]*m_to_pc)
        line_6.set_ydata(NGC_2451A['y'][:i]*m_to_pc)
        
        line_7.set_xdata(Collinder_135['x'][:i]*m_to_pc)
        line_7.set_ydata(Collinder_135['y'][:i]*m_to_pc)
        
        line_8.set_xdata(LSR['x'][:i]*m_to_pc)
        line_8.set_ydata(LSR['y'][:i]*m_to_pc)
        
        line_9.set_xdata(Sun['x'][:i]*m_to_pc)
        line_9.set_ydata(Sun['y'][:i]*m_to_pc)

        ax.set_title('Galactic Frame, Time: %s Myr'%time)
        plt.tight_layout()
        plt.legend(loc=0)
        ax.set_xlabel('X (pc)')
        ax.set_ylabel('Y (pc)')
        return mplfig_to_npimage(fig_mpl)
   
    animation =mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_gif("around_gal_center.gif", fps=20)


def LSR_animation_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,time_step_length):
    
    duration = 10 # time in seconds the animation lasts for
    fig_mpl, ax = plt.subplots(1,figsize=(5,5), facecolor='white')
    
    plt.scatter(-(USco['x'][0] - LSR['x'][0])*m_to_pc,(USco['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(UCL['x'][0] - LSR['x'][0])*m_to_pc,(UCL['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(LCC['x'][0] - LSR['x'][0])*m_to_pc,(LCC['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(IC_2602['x'][0] - LSR['x'][0])*m_to_pc,(IC_2602['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(IC_2391['x'][0] - LSR['x'][0])*m_to_pc,(IC_2391['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(NGC_2451A['x'][0] - LSR['x'][0])*m_to_pc,(NGC_2451A['y'][0] - LSR['y'][0])*m_to_pc,s=20)
    plt.scatter(-(Collinder_135['x'][0] - LSR['x'][0])*m_to_pc,(Collinder_135['y'][0] - LSR['y'][0])*m_to_pc,s=20)  
  
    line, = ax.plot(-(USco['x'] - LSR['x'])*m_to_pc,(USco['y'] - LSR['y'])*m_to_pc,label='USco',lw=2)
    line_2, = ax.plot(-(UCL['x'] - LSR['x'])*m_to_pc,(UCL['y'] - LSR['y'])*m_to_pc,label='UCL',lw=2)
    line_3, = ax.plot(-(LCC['x'] - LSR['x'])*m_to_pc,(LCC['y'] - LSR['y'])*m_to_pc,label='LCC',lw=2)
    line_4, = ax.plot(-(IC_2602['x'] - LSR['x'])*m_to_pc,(IC_2602['y'] - LSR['y'])*m_to_pc,label='IC 2602',lw=2)
    line_5, = ax.plot(-(IC_2391['x'] - LSR['x'])*m_to_pc,(IC_2391['y'] - LSR['y'])*m_to_pc,label='IC 2391',lw=2)
    line_6, = ax.plot(-(NGC_2451A['x'] - LSR['x'])*m_to_pc,(NGC_2451A['y'] - LSR['y'])*m_to_pc,label='NGC 2451A',lw=2)
    line_7, = ax.plot(-(Collinder_135['x'] - LSR['x'])*m_to_pc,Collinder_135['y']*m_to_pc,label='Collinder 135',lw=2)
    #line_8, = ax.plot(-(Sun['x'] - LSR['x'])*m_to_pc,(Sun['y'] - LSR['y'])*m_to_pc,label='Sun',lw=2)
   
    ax.set_xlim(-100,300)
    ax.set_ylim(-250,150)
    
    
    def make_frame_mpl(t):
        
        i = np.int(len(USco)*(t/duration))
        time = np.round((time_step_length*i)*seconds_to_myr,0)
        
        line.set_xdata(-(USco['x'][:i] - LSR['x'][:i])*m_to_pc)
        line.set_ydata((USco['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_2.set_xdata(-(UCL['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_2.set_ydata((UCL['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_3.set_xdata(-(LCC['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_3.set_ydata((LCC['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_4.set_xdata(-(IC_2602['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_4.set_ydata((IC_2602['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_5.set_xdata(-(IC_2391['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_5.set_ydata((IC_2391['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_6.set_xdata(-(NGC_2451A['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_6.set_ydata((NGC_2451A['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        line_7.set_xdata(-(Collinder_135['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_7.set_ydata((Collinder_135['y'][:i] - LSR['y'][:i])*m_to_pc)
        
        #line_8.set_xdata(-(Sun['x'][:i] - LSR['x'][:i])*m_to_pc)
        #line_8.set_ydata((Sun['y'][:i] - LSR['y'][:i])*m_to_pc)
        plt.tight_layout()
        plt.legend(loc=1)
        ax.set_title('LSR Frame, Time: %s Myr'%time)
        ax.set_xlabel('X (pc)')
        ax.set_ylabel('Y (pc)')
        return mplfig_to_npimage(fig_mpl)
   
    animation =mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_gif("LSR_XY.gif", fps=20)

def LSR_animation_XZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,time_step_length):
    
    duration = 10 # time in seconds the animation lasts for
    fig_mpl, ax = plt.subplots(1,figsize=(5,5), facecolor='white')
    
    plt.scatter(-(USco['x'][0] - LSR['x'][0])*m_to_pc,(USco['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(UCL['x'][0] - LSR['x'][0])*m_to_pc,(UCL['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(LCC['x'][0] - LSR['x'][0])*m_to_pc,(LCC['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(IC_2602['x'][0] - LSR['x'][0])*m_to_pc,(IC_2602['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(IC_2391['x'][0] - LSR['x'][0])*m_to_pc,(IC_2391['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(NGC_2451A['x'][0] - LSR['x'][0])*m_to_pc,(NGC_2451A['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter(-(Collinder_135['x'][0] - LSR['x'][0])*m_to_pc,(Collinder_135['z'][0] - LSR['z'][0])*m_to_pc,s=20)  
  
    line, = ax.plot(-(USco['x'] - LSR['x'])*m_to_pc,(USco['z'] - LSR['z'])*m_to_pc,label='USco',lw=2)
    line_2, = ax.plot(-(UCL['x'] - LSR['x'])*m_to_pc,(UCL['z'] - LSR['z'])*m_to_pc,label='UCL',lw=2)
    line_3, = ax.plot(-(LCC['x'] - LSR['x'])*m_to_pc,(LCC['z'] - LSR['z'])*m_to_pc,label='LCC',lw=2)
    line_4, = ax.plot(-(IC_2602['x'] - LSR['x'])*m_to_pc,(IC_2602['z'] - LSR['z'])*m_to_pc,label='IC 2602',lw=2)
    line_5, = ax.plot(-(IC_2391['x'] - LSR['x'])*m_to_pc,(IC_2391['z'] - LSR['z'])*m_to_pc,label='IC 2391',lw=2)
    line_6, = ax.plot(-(NGC_2451A['x'] - LSR['x'])*m_to_pc,(NGC_2451A['z'] - LSR['z'])*m_to_pc,label='NGC 2451A',lw=2)
    line_7, = ax.plot(-(Collinder_135['x'] - LSR['x'])*m_to_pc,(Collinder_135['z'] - LSR['z'])*m_to_pc,label='Collinder 135',lw=2)
    #line_8, = ax.plot(-(Sun['x'] - LSR['x'])*m_to_pc,(Sun['z'] - LSR['z'])*m_to_pc,label='Sun',lw=2)
   
    ax.set_xlim(-100,300)
    ax.set_ylim(-200,200)
    
    
    def make_frame_mpl(t):
        
        i = np.int(len(USco)*(t/duration))
        time = np.round((time_step_length*i)*seconds_to_myr,0)
        
        line.set_xdata(-(USco['x'][:i] - LSR['x'][:i])*m_to_pc)
        line.set_ydata((USco['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_2.set_xdata(-(UCL['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_2.set_ydata((UCL['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_3.set_xdata(-(LCC['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_3.set_ydata((LCC['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_4.set_xdata(-(IC_2602['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_4.set_ydata((IC_2602['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_5.set_xdata(-(IC_2391['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_5.set_ydata((IC_2391['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_6.set_xdata(-(NGC_2451A['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_6.set_ydata((NGC_2451A['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_7.set_xdata(-(Collinder_135['x'][:i] - LSR['x'][:i])*m_to_pc)
        line_7.set_ydata((Collinder_135['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        #line_8.set_xdata(-(Sun['x'][:i] - LSR['x'][:i])*m_to_pc)
        #line_8.set_ydata((Sun['z'][:i] - LSR['z'][:i])*m_to_pc)
        plt.tight_layout()
        ax.set_title('LSR Frame, Time: %s Myr'%time)
        ax.set_xlabel('X (pc)')
        ax.set_ylabel('Z (pc)')
        return mplfig_to_npimage(fig_mpl)
   
    animation =mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_gif("LSR_XZ.gif", fps=20)

def LSR_animation_YZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,time_step_length):
    
    duration = 10 # time in seconds the animation lasts for
    fig_mpl, ax = plt.subplots(1,figsize=(5,5), facecolor='white')
    
    plt.scatter((USco['y'][0] - LSR['y'][0])*m_to_pc,(USco['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((UCL['y'][0] - LSR['y'][0])*m_to_pc,(UCL['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((LCC['y'][0] - LSR['y'][0])*m_to_pc,(LCC['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((IC_2602['y'][0] - LSR['y'][0])*m_to_pc,(IC_2602['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((IC_2391['y'][0] - LSR['y'][0])*m_to_pc,(IC_2391['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((NGC_2451A['y'][0] - LSR['y'][0])*m_to_pc,(NGC_2451A['z'][0] - LSR['z'][0])*m_to_pc,s=20)
    plt.scatter((Collinder_135['y'][0] - LSR['y'][0])*m_to_pc,(Collinder_135['z'][0] - LSR['z'][0])*m_to_pc,s=20)  
  
    line, = ax.plot((USco['y'] - LSR['y'])*m_to_pc,(USco['z'] - LSR['z'])*m_to_pc,label='USco',lw=2)
    line_2, = ax.plot((UCL['y'] - LSR['y'])*m_to_pc,(UCL['z'] - LSR['z'])*m_to_pc,label='UCL',lw=2)
    line_3, = ax.plot((LCC['y'] - LSR['y'])*m_to_pc,(LCC['z'] - LSR['z'])*m_to_pc,label='LCC',lw=2)
    line_4, = ax.plot((IC_2602['y'] - LSR['y'])*m_to_pc,(IC_2602['z'] - LSR['z'])*m_to_pc,label='IC 2602',lw=2)
    line_5, = ax.plot((IC_2391['y'] - LSR['y'])*m_to_pc,(IC_2391['z'] - LSR['z'])*m_to_pc,label='IC 2391',lw=2)
    line_6, = ax.plot((NGC_2451A['y'] - LSR['y'])*m_to_pc,(NGC_2451A['z'] - LSR['z'])*m_to_pc,label='NGC 2451A',lw=2)
    line_7, = ax.plot((Collinder_135['y'] - LSR['y'])*m_to_pc,(Collinder_135['z'] - LSR['z'])*m_to_pc,label='Collinder 135',lw=2)
    #line_8, = ax.plot((Sun['y'] - LSR['y'])*m_to_pc,(Sun['z'][0] - LSR['z'][0])*m_to_pc,label='Sun',lw=2)
   
    ax.set_xlim(-250,150)
    ax.set_ylim(-200,200)
    
    
    def make_frame_mpl(t):
        
        i = np.int(len(USco)*(t/duration))
        time = np.round((time_step_length*i)*seconds_to_myr,0)
        
        line.set_xdata((USco['y'][:i] - LSR['y'][:i])*m_to_pc)
        line.set_ydata((USco['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_2.set_xdata((UCL['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_2.set_ydata((UCL['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_3.set_xdata((LCC['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_3.set_ydata((LCC['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_4.set_xdata((IC_2602['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_4.set_ydata((IC_2602['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_5.set_xdata((IC_2391['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_5.set_ydata((IC_2391['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_6.set_xdata((NGC_2451A['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_6.set_ydata((NGC_2451A['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        line_7.set_xdata((Collinder_135['y'][:i] - LSR['y'][:i])*m_to_pc)
        line_7.set_ydata((Collinder_135['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        #line_8.set_xdata((Sun['y'][:i] - LSR['y'][:i])*m_to_pc)
        #line_8.set_ydata((Sun['z'][:i] - LSR['z'][:i])*m_to_pc)
        
        plt.tight_layout()
        ax.set_title('LSR Frame, Time: %s Myr'%time)
        ax.set_xlabel('Y (pc)')
        ax.set_ylabel('Z (pc)')
        return mplfig_to_npimage(fig_mpl)
   
    animation =mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_gif("LSR_YZ.gif", fps=20)

def combine_2_gifs(gif_name_1,gif_name_2):
    clip_1 = mpy.VideoFileClip(gif_name_1)
    clip_2 = mpy.VideoFileClip(gif_name_2).resize(height=clip_1.h)
    animation = mpy.clips_array([[clip_1, clip_2]])
    animation.write_gif("synced_2.gif", fps=20)

def combine_3_gifs(gif_name_1,gif_name_2,gif_name_3):
    clip_1 = mpy.VideoFileClip(gif_name_1)
    clip_2 = mpy.VideoFileClip(gif_name_2).resize(height=clip_1.h)
    clip_3 = mpy.VideoFileClip(gif_name_3).resize(height=clip_1.h)

    animation = mpy.clips_array([[clip_1, clip_2,clip_3]])
    animation.write_gif("synced_3.gif", fps=20)


#plot_around_gal_center_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
plot_LSR_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
#plot_LSR_XZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
#plot_LSR_YZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
#plot_3D_LSR(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
#Fernandez_2007_coords(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun)
#LSR_animation_XY(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,-63080000000.0)
#LSR_animation_XZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,-63080000000.0)
#LSR_animation_YZ(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,-63080000000.0)


#around_gal_center_animation(USco,UCL,LCC,IC_2602,IC_2391,NGC_2451A,Collinder_135,LSR,Sun,-63080000000.0)
#combine_2_gifs('around_gal_center.gif','LSR.gif')
#combine_3_gifs("LSR_XY.gif","LSR_XZ.gif","LSR_YZ.gif")