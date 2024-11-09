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

def BaseBallRobot (tau,m,d,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=True): 
    r0 = np.array([0,y0]) # initial position vector 
    speed = np.array([v0*np.cos(theta_0*np.pi/180),v0*np.sin(theta_0*np.pi/180)]) # initial velocity vector 
    r = np.copy(r0)   # Set initial position 
    v = np.copy(speed)   # Set initial velocity
    maxstep = 1000
    xplot = np.empty(maxstep)
    yplot = np.empty(maxstep)

    if AirResistance == True: 
        rho = 0
    else: 
        rho = rho_input

    air_const = -0.5*C_d*rho*A/m

    if str == "Euler": 
        # Euler Method
        print('Euler Method')

        for istep in range(maxstep): 
            #* Record position (computed and theoretical) for plotting
            xplot[istep] = r[0]   # Record trajectory for plot
            yplot[istep] = r[1]
            t = istep*tau         # Current time

            accel = air_const*np.linalg.norm(v)*v
            accel[1] = accel[1] - g

            #* Calculate the new position and velocity using Euler method
            r = r + tau*v                    # Euler step
            v = v + tau*accel  
            
            #* If ball reaches ground (y<0), break out of the loop
            if r[1] < 0 : 
                laststep = istep+1
                xplot[laststep] = r[0]  # Record last values computed
                yplot[laststep] = r[1]
                break                   # Break out of the for loop

    elif str == "Euler-Cromer": 
        # Euler-Cromer Method 
        print('Euler-Cromer Method')

        # v(n+1) = vn + tau*an
        # r(n+1) = rn + tau*(v(n+1)+v(n))/2

        for istep in range(maxstep): 
            #* Record position (computed and theoretical) for plotting
            xplot[istep] = r[0]   # Record trajectory for plot
            yplot[istep] = r[1]
            t = istep*tau         # Current time

            accel = air_const*np.linalg.norm(v)*v
            accel[1] = accel[1] - g

            #* Calculate the new position and velocity using Euler method
            v_old = v
            v = v_old + tau*accel  
            r = r + tau*(v+v_old)/2  

            print(v+v_old)              

            
            #* If ball reaches ground (y<0), break out of the loop
            if r[1] < 0 : 
                laststep = istep+1
                xplot[laststep] = r[0]  # Record last values computed
                yplot[laststep] = r[1]
                break                   # Break out of the for loop     
            #
    elif str == "Midpoint":
        # Midpoint Method 
        print('Midpoint Method')
        
        # v(n+1) = vn + tau*an
        # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2

        for istep in range(maxstep): 
            #* Record position (computed and theoretical) for plotting
            xplot[istep] = r[0]   # Record trajectory for plot
            yplot[istep] = r[1]
            t = istep*tau         # Current time

            accel = air_const*np.linalg.norm(v)*v
            accel[1] = accel[1] - g

            #* Calculate the new position and velocity using Euler method
            v_old = v
            v = v_old + tau*accel  
            r = r + tau*v_old + (1/2)*accel*tau**2
            
            #* If ball reaches ground (y<0), break out of the loop
            if r[1] < 0 : 
                laststep = istep+1
                xplot[laststep] = r[0]  # Record last values computed
                yplot[laststep] = r[1]
                break                   # Break out of the for loop

        #
    else: 
        print("This is not a valid entry for computation method. Please enter one of the following character strings: ")
        print("'Euler'")
        print("'Euler-Cromer'")
        print("'Midpoint'")

    #* Print maximum range and time of flight
    print('Maximum range is', r[0], 'meters')
    print('Time of flight is', laststep*tau , ' seconds')

    return xplot, yplot



# Part 2 

mu_v0 = 100 
std_v0 = 15
size0 = (1,)
v0 = np.random.normal(mu_v0, std_v0,size0)/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = np.random.normal(mu_theta_0,std_theta_0,size0)
