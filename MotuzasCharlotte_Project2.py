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
A = np.pi*(d/2)**2 m^2
y0 = 1.0 m'''

# Part 1 

def RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=1,printon = 1): 
    '''This function is designed to obtain the trajectory of a baseball as it would be hit by a robotic designated hitter. 
    This function takes the following inputs: 
    tau: time step length 
    m: mass of the baseball 
    g: gravitational acceleration 
    v0: initial velocity (scalar)
    theta_0: initial launch angle 
    y0: initial height of the ball 
    A: cross-sectional area of the baseball 
    C_d: drag coefficient in air resistance case
    rho_input: air density in air resistance case 
    str: method chosen to complete numerical solution; options are 'Euler', 'Euler-Cromer' and 'Midpoint'
    AirResistance: equal to 1 automatically, if not equal to 1, air resistance will be neglected. 
    printon: automatically set to 1, if not equal to one then the message displaying the range and time of flight will not be printed in the terminal. 
    
    This function returns xplot, yplot, laststep, xNoAir, yNoAir, maxrange. 
    xplot and yplot are the x and y values of the solution for the chosen method, laststep is the index of the last time step calculated, 
    xNoAir and yNoAir are the x and y values for the analytical solution in the negligible air resistance case, and maxrange is the maximum range travelled by the ball.'''
    
    
    r0 = np.array([0,y0]) # initial position vector 
    v0 = np.array([v0*np.cos(theta_0*np.pi/180),v0*np.sin(theta_0*np.pi/180)]) # initial velocity vector 
    r = np.copy(r0)   # Set initial position 
    v = np.copy(v0)   # Set initial velocity
    maxstep = 1000 # max iteration number 
    # initiate arrays for plotting 
    xplot = np.empty(maxstep) 
    yplot = np.empty(maxstep)
    xNoAir = np.empty(maxstep)
    yNoAir = np.empty(maxstep)

    if AirResistance == 1: # turns on air resistance 
        rho = rho_input
    elif rho_input or C_d == None: 
        rho = 0
    else: 
        rho = 0

    air_const = -0.5*C_d*rho*A/m # set of constants as needed to determine acceleration from air resistance 

    for istep in range(maxstep): 
 
        accel = air_const*np.sqrt(v[0]**2 + v[1]**2)*v # acceleration according to ODE defined in assignment outline 
        accel[1] = accel[1] - g # adding gravitational acceleration in the y direction 
        v_step = v + tau*accel  # each method uses the same numerical solution for v
        r_old = r       # recording for later range calculations 

        if str == "Euler": 

            #* Calculate the new position and velocity using Euler method
            r_step = r + tau*v                   # Euler step
            v, r = v_step,r_step

        elif str == "Euler-Cromer": 
            # Euler-Cromer Method 

            # v(n+1) = vn + tau*an
            # r(n+1) = rn + tau*(v(n+1)+v(n))/
            
            #* Calculate the new position and velocity using Euler method
            r_step = r + tau*v_step              # Euler step
            v, r = v_step,r_step
            
        elif str == "Midpoint":
            # Midpoint Method 
        
            # v(n+1) = vn + tau*an
            # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2
            
            #* Calculate the new position and velocity using Euler method
            
            r_step = r + tau*(v_step+v)/2
            v, r = v_step,r_step

        else: 
            print("This is not a valid entry for computation method. Please enter one of the following character strings: ")
            print("'Euler'")
            print("'Euler-Cromer'")
            print("'Midpoint'")
            break
        
        # Record position (computed and theoretical) for plotting
        xplot[istep] = r[0]   # Record trajectory for plot
        yplot[istep] = r[1]
        t = istep*tau         # Current time
        xNoAir[istep] = r0[0] + v0[0]*t
        yNoAir[istep] = r0[1] + v0[1]*t - 0.5*g*t**2
        
        # If ball reaches ground (y<0), break out of the loop
        if r[1] < 0 : 
            laststep = istep+1
            xplot[laststep] = r[0]  # Record last values computed
            yplot[laststep] = r[1]

            maxrange = -r_old[1]*(r[0]-r_old[0])/(r[1]-r_old[1]) + r_old[0]
            tmax = (maxrange-r_old[0])*(tau*laststep-t)/(r[0]-r_old[0])+t

            break                   # Break out of the for loop

    if printon == 1: 
        # Print maximum range and time of flight
        print('Maximum range is', maxrange, 'meters')
        print('Time of flight is', tmax, ' seconds')

    return xplot, yplot, laststep, xNoAir, yNoAir, maxrange

# Constants 

m = 0.145 
d = 7.4/100 
g = 9.81 
rho = 1.2
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 0 
tau = 0.1
theta_0 = 45
v0 = 15

str = ['Euler', 'Euler-Cromer','Midpoint']

# loop which returns the max range for the given conditions for each of the three solution methods 

for i in range (3): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho,str[i],AirResistance=0,printon=1)
    if i == 0: 
        # plotting the ground for the first method used, not completed every time 
        xground = np.array([0., xplot[laststep-1]])
        yground = np.array([0., 0.])
    
        plt.plot(xNoAir[0:laststep], yNoAir[0:laststep], 'r-', # plotting the analytical case 
         xground,yground,'k-') # plotting the ground 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+') # plotting the numerical solution 
    else: 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+') # plotting the numerical solution 

plt.legend(['Theory (No air)','Ground','Euler Method','Euler-Cromer Method','Midpoint Method'], loc = 'best');
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.grid()
plt.show()

# Figure 2.3

v0 = 50 # intial velocity changed 
y0 = 1.0 # initial height changed 

for i in range (3): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho,str[i],AirResistance=1,printon=1)
    if i == 0:
        # plotting the ground for the first method used, not completed every time 
        xground = np.array([0., xNoAir[laststep-1]])
        yground = np.array([0., 0.])
    
        plt.plot(xNoAir[0:laststep], yNoAir[0:laststep], 'r-', # analytical case
         xground,yground,'k-') # ground 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+') # numerical solution 
    else: 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+') # numerical solution 

plt.legend(['Theory (No air)','Ground','Euler Method','Euler-Cromer Method','Midpoint Method']);
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.grid()
plt.show()

# Part 2 

rho = 1.2 # air density, kg/m^3
y0 = 1.0 # initial height 

N = 100 # number of iterations 

mu_v0 = 100 # mean for initial velocity distribution 
std_v0 = 15 # standard deviation for initial velocity distribution 
v0 = (mu_v0 + std_v0*np.random.randn(N))/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45 # mean for initial angle distribution 
std_theta_0 = 10 # standard deviation for initial angle distribution 
theta_0 = (mu_theta_0 + std_theta_0*np.random.randn(N))
str = "Midpoint" # midpoint was selected for accuracy 

count = 0 # initial count is zero 
for i in range(N): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0[i],theta_0[i],y0,A,C_d,rho,str,AirResistance=1,printon=0)
    if maxrange*3.281 >= 400: # converting max range to feet, comparing to 400 ft home-run designation 
        count = count + 1

ABHR_Ratio = N/count # determining at bat - home run ratio 
print(ABHR_Ratio)

# Part 3 

def RDHFence(tau,m,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=1,printon = 1): 
    '''This function is designed to obtain the trajectory of a baseball as it would be hit by a robotic designated hitter. 
    This function takes the following inputs: 
    tau: time step length 
    m: mass of the baseball 
    g: gravitational acceleration 
    v0: initial velocity (scalar)
    theta_0: initial launch angle 
    y0: initial height of the ball 
    A: cross-sectional area of the baseball 
    C_d: drag coefficient in air resistance case
    rho_input: air density in air resistance case 
    str: method chosen to complete numerical solution; options are 'Euler', 'Euler-Cromer' and 'Midpoint'
    AirResistance: equal to 1 automatically, if not equal to 1, air resistance will be neglected. 
    printon: automatically set to 1, if not equal to one then the message displaying the range and time of flight will not be printed in the terminal. 
    
    This function returns the value 'height', which is the height of the baseball relative to the reference point y=0. may be negative or positive.'''
    
    
    r0 = np.array([0,y0]) # initial position vector 
    v0 = np.array([v0*np.cos(theta_0*np.pi/180),v0*np.sin(theta_0*np.pi/180)]) # initial velocity vector 
    r = np.copy(r0)   # Set initial position 
    v = np.copy(v0)   # Set initial velocity
    maxstep = 1000 # max number of iterations 
    # initialize arrays for plotting 
    xplot = np.empty(maxstep)
    yplot = np.empty(maxstep)
    xNoAir = np.empty(maxstep)
    yNoAir = np.empty(maxstep)

    if AirResistance == 1: # air resistance case 
        rho = rho_input
    elif rho_input == None: # negligible air resistance case 
        rho = 0
    else: 
        rho = 0

    air_const = -0.5*C_d*rho*A/m # constants for acceleration due to air resistance 

    for istep in range(maxstep): 
 
        accel = air_const*np.sqrt(v[0]**2 + v[1]**2)*v # acceleration from ODE given in assignment outline 
        accel[1] = accel[1] - g # adding gravitational acceleration in y direction 
        v_step = v + tau*accel  # v is the same for each of the solution methods 
        r_old = r  # used for determining the height later on 


        if str == "Euler": 

            #* Calculate the new position and velocity using Euler method
            r_step = r + tau*v                   # Euler step
            v, r = v_step,r_step

        elif str == "Euler-Cromer": 
            # Euler-Cromer Method 

            # v(n+1) = vn + tau*an
            # r(n+1) = rn + tau*(v(n+1)+v(n))/

            #* Calculate the new position and velocity using Euler method
            r_step = r + tau*v_step                  # Euler step
            v, r = v_step,r_step
            
        elif str == "Midpoint":
            # Midpoint Method 
        
            # v(n+1) = vn + tau*an
            # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2
            
            #* Calculate the new position and velocity using Euler method
            
            r_step = r + tau*(v_step+v)/2 
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
        if r[0] > 400/3.281 : 
            laststep = istep+1
            xplot[laststep] = r[0]  # Record last values computed
            yplot[laststep] = r[1]

            slope = (r[1]-r_old[1])/(r[0]-r_old[0])
            b = r_old[1] - slope*r_old[0]

            height = slope*400/3.281 + b

            break                   # Break out of the for loop
        elif r[1] < -50: 
            height = -50 # still a height that would discount it from being a home run, but cuts off the loop 
            break 
    
    if printon == 1: 
        #* Print maximum range and time of flight
        print('height at 400 feet', height,'meters')

    return height


# same setup as part 2

rho = 1.2
y0 = 1.0 
N = 100
mu_v0 = 100 
std_v0 = 15
v0 = (mu_v0 + std_v0*np.random.randn(N))/2.237 # divided by 2.237 to obtain value in m/s rather than mph
mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = (mu_theta_0 + std_theta_0*np.random.randn(N))
str = "Midpoint"

fenceheight = np.linspace(0.5,15,30) # array of fence heights from 0.5 to 15 

ABHR_Ratio = np.empty(30) # initializing array for AB/HR value for each fence height 

lessthan10 = np.zeros(30) # array that will be used to identify which fence height has an AB/HR less than 10 
for j in range(30): 
    count = 0
    for i in range(N): 
        height = RDHFence(tau,m,g,v0[i],theta_0[i],y0,A,C_d,rho,str,AirResistance=1,printon = 0)
        if height >= fenceheight[j]:  # if the height is greater than the fence height at 400 ft, count + 1
            count = count + 1

    ABHR_Ratio[j] = N/count # determining AB/HR ratio for each fence height 
    if ABHR_Ratio[j] < 10: 
        lessthan10[j] = j

minfenceheight = fenceheight[np.where(lessthan10 == max(lessthan10))] # identify the last fence height that has an AB/HR ratio less than 10 
print('Minimum fence height should be ',minfenceheight[0]+0.5, 'm')