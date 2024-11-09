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

m = 0.145 kg. 
d = 7.4 cm 
g = 9.81 m/s^2
rho = 1.2 kg/m^3 
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 m'''

# Part 1 

def BaseBallRobot (tau,m,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=1): 
    r0 = np.array([0,y0]) # initial position vector 
    v0 = np.array([v0*np.cos(theta_0*np.pi/180),v0*np.sin(theta_0*np.pi/180)]) # initial velocity vector 
    r = np.copy(r0)   # Set initial position 
    v = np.copy(v0)   # Set initial velocity
    maxstep = 1000
    xplot = np.empty(maxstep)
    yplot = np.empty(maxstep)
    xNoAir = np.empty(maxstep)
    yNoAir = np.empty(maxstep)

    if AirResistance == 1: 
        rho = rho_input
    else: 
        rho = 0

    air_const = -0.5*C_d*rho*A/m

    print(AirResistance)

    for istep in range(maxstep): 
 
        accel = air_const*np.sqrt(v[0]**2 + v[1]**2)*v
        accel[1] = accel[1] - g

        if str == "Euler": 

            print(v)

            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v                   # Euler step
            v, r = v_step,r_step

        


        elif str == "Euler-Cromer": 
            # Euler-Cromer Method 

            # v(n+1) = vn + tau*an
            # r(n+1) = rn + tau*(v(n+1)+v(n))/

            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*(v_step+v)/2                  # Euler step
            v, r = v_step,r_step

            print(air_const)         
            
        elif str == "Midpoint":
            # Midpoint Method 
        
            # v(n+1) = vn + tau*an
            # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2
            
            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v + (1/2)*accel*tau**2
            v, r = v_step,r_step

        else: 
            print("This is not a valid entry for computation method. Please enter one of the following character strings: ")
            print("'Euler'")
            print("'Euler-Cromer'")
            print("'Midpoint'")
            break
        
        #* Record position (computed and theoretical) for plotting
        xplot[istep] = r[0]   # Record trajectory for plot
        yplot[istep] = r[1]
        t = istep*tau         # Current time
        xNoAir[istep] = r0[0] + v0[0]*t
        yNoAir[istep] = r0[1] + v0[1]*t - 0.5*g*t**2
        
        #* If ball reaches ground (y<0), break out of the loop
        if r[1] < 0 : 
            laststep = istep+1
            xplot[laststep] = r[0]  # Record last values computed
            yplot[laststep] = r[1]
            break                   # Break out of the for loop

    #* Print maximum range and time of flight
    print('Maximum range is', r[0], 'meters')
    print('Time of flight is', laststep*tau , ' seconds')

    return xplot, yplot, laststep, xNoAir, yNoAir

# Constants 

m = 0.145 
d = 7.4/100 
g = 9.81 
rho = 1.2
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 
tau = 0.1
theta_0 = 45
v0 = 50
y0 = 1

str = "Euler"

xplot, yplot, laststep, xNoAir, yNoAir = BaseBallRobot (tau,m,g,v0,theta_0,y0,A,C_d,rho,str,AirResistance=1)


# Graph the trajectory of the baseball
# Mark the location of the ground by a straight line
xground = np.array([0., xplot[laststep-1]])
yground = np.array([0., 0.])
# Plot the computed trajectory and parabolic, no-air curve
plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], 'b+',
         xNoAir[0:laststep], yNoAir[0:laststep], 'r-',
         xground,yground,'k-')
plt.legend(['Euler method', 'Theory (No air)']);
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.show()








# Part 2 

mu_v0 = 100 
std_v0 = 15
size0 = (1,)
v0 = np.random.normal(mu_v0, std_v0,size0)/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = np.random.normal(mu_theta_0,std_theta_0,size0)
