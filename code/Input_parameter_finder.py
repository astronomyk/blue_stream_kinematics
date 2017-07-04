# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:20:41 2017

@author: Dan
"""
import astropy
from astropy.table import Table

def integrator_step_length(time_to_integrate,number_of_steps):

    seconds_in_a_year = 3.154e+7
    
    time_to_integrate_sec = time_to_integrate*seconds_in_a_year
    time_per_step_sec = time_to_integrate_sec/number_of_steps

    print(time_per_step_sec, 'time per step')

    
    
def pos_and_vel_cal_old_LSR(x,y,z,u,v,w):
    
    x_LSR = 8500 - x
    y_LSR = y
    z_LSR = z + 10
    
    u_LSR = -(u + 10.0) 
    v_LSR = v + 220 + 5.25
    w_LSR = w + 7.17
    
    print(x_LSR,'X LSR')
    print(y_LSR,'Y LSR')
    print(z_LSR,'Z LSR')
    
    print(u_LSR,'U LSR')
    print(v_LSR,'V LSR')
    print(w_LSR,'W LSR')
    
def pos_and_vel_cal_new_LSR(x,y,z,u,v,w):
    
    x_LSR = 8000 - x
    y_LSR = y
    z_LSR = z + 10
    
    u_LSR = -(u + 11.1) 
    v_LSR = v + 229.76 + 12.24
    w_LSR = w + 7.25
    
    print(x_LSR,'x LSR')
    print(y_LSR,'y LSR')
    print(z_LSR,'z LSR')
    
    print(u_LSR,'U LSR')
    print(v_LSR,'V LSR')
    print(w_LSR,'W LSR')



data = Table.read('Sco_cen_blue_stream_data_Sco_cen_from_S03.csv',format='csv')

for i in range(0,len(data)):
    print('')
    print(data['AltName'][i])
    pos_and_vel_cal_new_LSR(data['X'][i],data['Y'][i],data['Z'][i],data['U'][i],data['V'][i],data['W'][i])
    print('')

integrator_step_length(20e6,10000)
