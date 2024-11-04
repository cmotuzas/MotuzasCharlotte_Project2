# Project 2 - Charlotte Motuzas 
# Due Nov. 8th, 2024 
# Physics 3926 

import numpy as np 
import scipy as scp 
import matplotlib.pyplot as plt 
import scipy.integrate as scint 

# Background 
''' - Distribution of ball speed off of the bat, gaussian of mean 100 mph with standard deviation 15 mph
- launch angle gaussian mean of 45 degrees with standard deviation of 10 degrees 
- determine how good of a home run hitter this RDH actually is 
- compute expected ratio of at bats to home runs, the AB/HR ratio 
- top hitters have a ratio of about 10
- model the ball as a projectile moving under combined forces gravity and air resistance 
- motion confined to xy plane, y being the height  above the ground 

dv/dt = -g(yhat) - v(C_d*rho*A*v)/2m 
dr/dt = v

v seems to be the velocity, r is the displacement of the ball from the origin 
yhat is the unit vector in the y direction

m = 0.145 km. 
d = 7.4 cm 
g = 9.81 m/s^2
rho = 1.2 kg/m^3 
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 m'''

# Constants 

m = 0.145  
d = 7.4 
g = 9.81 
rho = 1.2 
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 


# Part 1 

def BaseBallRobot (tau,v0,theta_0,y0,AirResistance=True): 
    # Euler Method



    # Euler-Cromer Method 



    # Midpoint Method 
    return 1



# Part 2 

mu_v0 = 100 
std_v0 = 15
size0 = (1,)
v0 = np.random.normal(mu_v0, std_v0,size0)/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = np.random.normal(mu_theta_0,std_theta_0,size0)
