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

def RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=1,printon = 1): 
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

    for istep in range(maxstep): 
 
        accel = air_const*np.sqrt(v[0]**2 + v[1]**2)*v
        accel[1] = accel[1] - g

        if str == "Euler": 

            r_old = r


            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v                   # Euler step
            v, r = v_step,r_step

        


        elif str == "Euler-Cromer": 
            # Euler-Cromer Method 

            # v(n+1) = vn + tau*an
            # r(n+1) = rn + tau*(v(n+1)+v(n))/

            r_old = r

            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v_step              # Euler step
            v, r = v_step,r_step
            
        elif str == "Midpoint":
            # Midpoint Method 
        
            # v(n+1) = vn + tau*an
            # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2
            
            #* Calculate the new position and velocity using Euler method
            
            r_old = r
            
            v_step = v + tau*accel  
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
        if r[1] < 0 : 
            laststep = istep+1
            xplot[laststep] = r[0]  # Record last values computed
            yplot[laststep] = r[1]

            maxrange = -r_old[1]*(r[0]-r_old[0])/(r[1]-r_old[1]) + r_old[0]
            tmax = (maxrange-r_old[0])*(tau*laststep-t)/(r[0]-r_old[0])+t

            break                   # Break out of the for loop


    if printon == 1: 
        #* Print maximum range and time of flight
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

for i in range (3): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho,str[i],AirResistance=0,printon=1)
    if i == 0:
        xground = np.array([0., xplot[laststep-1]])
        yground = np.array([0., 0.])
    
        plt.plot(xNoAir[0:laststep], yNoAir[0:laststep], 'r-',
         xground,yground,'k-')
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+')
    else: 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+')

plt.legend(['Theory (No air)','Ground','Euler Method','Euler-Cromer Method','Midpoint Method']);
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.grid()
plt.show()

v0 = 50
y0 = 1.0

for i in range (3): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0,theta_0,y0,A,C_d,rho,str[i],AirResistance=1,printon=1)
    if i == 0:
        xground = np.array([0., xNoAir[laststep-1]])
        yground = np.array([0., 0.])
    
        plt.plot(xNoAir[0:laststep], yNoAir[0:laststep], 'r-',
         xground,yground,'k-')
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+')
    else: 
        plt.plot(xplot[0:laststep+1], yplot[0:laststep+1], '+')

plt.legend(['Theory (No air)','Ground','Euler Method','Euler-Cromer Method','Midpoint Method']);
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.grid()
plt.show()




# Part 2 

m = 0.145 
d = 7.4/100 
g = 9.81 
rho = 1.2
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 
tau = 0.1

N = 100

mu_v0 = 100 
std_v0 = 15
v0 = (mu_v0 + std_v0*np.random.randn(N))/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = (mu_theta_0 + std_theta_0*np.random.randn(N))
str = "Midpoint"

rangearray = np.empty(N)
count = 0 
for i in range(N): 
    xplot, yplot, laststep, xNoAir, yNoAir, maxrange = RoboticDesignatedHitter(tau,m,g,v0[i],theta_0[i],y0,A,C_d,rho,str,AirResistance=1,printon=0)
    rangearray[i] = maxrange*3.281 # in feet 
    if maxrange*3.281 >= 400: 
        count = count + 1

ABHR_Ratio = N/count
print(ABHR_Ratio)

# Part 3 

def RDHFence(tau,m,g,v0,theta_0,y0,A,C_d,rho_input,str,AirResistance=1,printon = 1): 
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

    for istep in range(maxstep): 
 
        accel = air_const*np.sqrt(v[0]**2 + v[1]**2)*v
        accel[1] = accel[1] - g

        if str == "Euler": 

            r_old = r


            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v                   # Euler step
            v, r = v_step,r_step

        elif str == "Euler-Cromer": 
            # Euler-Cromer Method 

            # v(n+1) = vn + tau*an
            # r(n+1) = rn + tau*(v(n+1)+v(n))/

            r_old = r

            #* Calculate the new position and velocity using Euler method
            v_step = v + tau*accel  
            r_step = r + tau*v_step                  # Euler step
            v, r = v_step,r_step
            
        elif str == "Midpoint":
            # Midpoint Method 
        
            # v(n+1) = vn + tau*an
            # --> r(n+1) = rn + tau*vn + (1/2)*an*tau^2
            
            #* Calculate the new position and velocity using Euler method
            
            r_old = r
            
            v_step = v + tau*accel  
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


m = 0.145 
d = 7.4/100 
g = 9.81 
rho = 1.2
C_d = 0.35
A = np.pi*(d/2)**2
y0 = 1.0 
tau = 0.1

N = 100
mu_v0 = 100 
std_v0 = 15
v0 = (mu_v0 + std_v0*np.random.randn(N))/2.237 # divided by 2.237 to obtain value in m/s rather than mph

mu_theta_0 = 45
std_theta_0 = 10 
theta_0 = (mu_theta_0 + std_theta_0*np.random.randn(N))
str = "Midpoint"

fenceheight = np.linspace(0.5,15,30)

ABHR_Ratio = np.empty(30)

lessthan10 = np.zeros(30)
for j in range(30): 
    count = 0
    for i in range(N): 
        height = RDHFence(tau,m,g,v0[i],theta_0[i],y0,A,C_d,rho,str,AirResistance=1,printon = 0)
        if height >= fenceheight[j]: 
            count = count + 1

    ABHR_Ratio[j] = N/count
    if ABHR_Ratio[j] < 10: 
        lessthan10[j] = j

minfenceheight = fenceheight[np.where(lessthan10 == max(lessthan10))]
print('Minimum fence height should be ',minfenceheight[0]+0.5, 'm')