# RJ Meade and Noah Thompson, 1/5/2023
# Cornell University and Louisiana Stupid University (LSU)
# Creates a movie showing the movement of a rodket with a constant thrust along the axis of the rod
# x and y are the horizontal and vertical position of the rod center of mass, theta is its angle to vertical
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as rand

# Classical 4th Order Runge-Kutta
# yi+1 = yi + 1/6(k1+2k2+2k3+k4)dt
# k1 = f(ti,yi)
# k2 = f(ti+1/2dt,yi+1/2k1dt)
# k3 = f(ti+1/2dt,yi+1/2k2dt)
# k4 = f(ti+dt,yi+k3dt)

# Solving diff. eqs. with 4th order Runge-Kutta, has placeholder functions in case diff. eq's get more complicated later

def DxDt(xdot):
    return xdot

def DxdotDt(mass,thrust,theta):
    return -thrust[0]/mass*np.sin(theta+thrust[1]) # horizontal thrust component

def DyDt(ydot):
    return ydot

def DydotDt(mass, thrust, theta):
    return thrust[0]/mass*np.cos(theta+thrust[1]) - 9.81 # vertical thrust component and gravity

def DthetaDt(thetadot):
    return thetadot

def DthetadotDt(length, mass, thrust): # would be torque/moment of inertia but there is no torque right now
    return -thrust[0]*np.sin(thrust[1])*(length/2)/(1/12*mass*length**2)

#def get_thrust(x,dx,y,dy,theta,dtheta,mass):
def get_thrust():
    thrust_mag = rand.randint(4820000,8450000) # wikipedia merlin 1d sea level
    thrust_dir = rand.uniform(-.175,.175)
    return thrust_mag, thrust_dir

# Input: Initial conditions (assumes initial t=0) and the constant force of thrust and mass and time step dt and number of time steps to run for (including first time step)
# Solves DxDt, DxdotDt, etc with given initial conditions using Runge-Kutta
# Returns arrays of x, y, and theta where each index is a time step of dt
def solve_eom(x0, dx0, y0, dy0, theta0, dtheta0, thrust0, length, mass, dt, totaldt):
    x = np.zeros(totaldt)
    y = np.zeros(totaldt)
    theta = np.zeros(totaldt)
    dx = np.zeros(totaldt)
    dy = np.zeros(totaldt)
    dtheta = np.zeros(totaldt)
    thrust = np.zeros(shape = (totaldt, 2))
    x[0] = x0
    y[0] = y0
    theta[0] = theta0
    dx[0] = dx0
    dy[0] = dy0
    dtheta[0] = dtheta0
    thrust[0][0], thrust[0][1] = thrust0
    
    for i in range(1,totaldt):
        thrust_i = get_thrust()
        thrust[i][0], thrust[i][1] = thrust_i
        # Slopes at beginning of interval
        k11 = DxDt(dx[i-1])
        k12 = DxdotDt(mass,thrust_i,theta[i-1])
        k13 = DyDt(dy[i-1])
        k14 = DydotDt(mass, thrust_i, theta[i-1])
        k15 = DthetaDt(dtheta[i-1])
        k16 = DthetadotDt(length, mass, thrust_i)
        # Estimate values at midpoint using initial slopes
        xmid = x[i-1] + k11*dt/2
        dxmid = dx[i-1] + k12*dt/2
        ymid = y[i-1] + k13*dt/2
        dymid = dy[i-1] + k14*dt/2
        thetamid = theta[i-1] + k15*dt/2
        dthetamid = dtheta[i-1] + k16*dt/2
        # Find slopes at estimated midpoint values
        k21 = DxDt(dxmid)
        k22 = DxdotDt(mass,thrust_i,thetamid)
        k23 = DyDt(dymid)
        k24 = DydotDt(mass, thrust_i, thetamid)
        k25 = DthetaDt(dthetamid)
        k26 = DthetadotDt(length, mass, thrust_i)
        # Use midpoint slope predictions to get better estimate of midpoint values from start of interval
        xmid = x[i-1] + k21*dt/2
        dxmid = dx[i-1] + k22*dt/2
        ymid = y[i-1] + k23*dt/2
        dymid = dy[i-1] + k24*dt/2
        thetamid = theta[i-1] + k25*dt/2
        dthetamid = dtheta[i-1] + k26*dt/2
        # Use new midpoint values to get new estimates for midpoint slopes
        k31 = DxDt(dxmid)
        k32 = DxdotDt(mass,thrust_i,thetamid)
        k33 = DyDt(dymid)
        k34 = DydotDt(mass, thrust_i, thetamid)
        k35 = DthetaDt(dthetamid)
        k36 = DthetadotDt(length, mass, thrust_i)
        # Use new midpoint slopes as estimate for entire interval to get estimated values at end of interval
        xend = x[i-1] + k31*dt
        dxend = dx[i-1] + k32*dt
        yend = y[i-1] + k33*dt
        dyend = dy[i-1] + k34*dt
        thetaend = theta[i-1] + k35*dt
        dthetaend = dtheta[i-1] + k36*dt
        # Use endpoint value estimates to get endpoint slope estimates
        k41 = DxDt(dxend)
        k42 = DxdotDt(mass,thrust_i, thetaend)
        k43 = DyDt(dyend)
        k44 = DydotDt(mass, thrust_i, thetaend)
        k45 = DthetaDt(dthetaend)
        k46 = DthetadotDt(length, mass, thrust_i)
        # Best estimate of slope across entire interval is a weighted average of all four slope estimates, midpoint slopes weighted higher
        x[i] = x[i-1] + 1/6*(k11+2*k21+2*k31+k41)*dt
        dx[i] = dx[i-1] + 1/6*(k12+2*k22+2*k32+k42)*dt
        y[i] = y[i-1] + 1/6*(k13+2*k23+2*k33+k43)*dt
        dy[i] = dy[i-1] + 1/6*(k14+2*k24+2*k34+k44)*dt
        theta[i] = theta[i-1] + 1/6*(k15+2*k25+2*k35+k45)*dt
        dtheta[i] = dtheta[i-1] + 1/6*(k16+2*k26+2*k36+k46)*dt
    return x,y,theta

def animate(i):
    length = 47.7
    x_axis_lim = 3000
    y_axis_lim = 15000
    x1 = x[i]-length*np.sin(theta[i])
    y1 = y[i]+length*np.cos(theta[i])
    x2 = x[i]+length*np.sin(theta[i])
    y2 = y[i]-length*np.cos(theta[i])
    print("x1: {}, x2: {}, y1: {}, y2: {}".format(x1, x2, y1, y2))
    ax.clear()
    plt.plot([x1, x2], [y1, y2], linewidth = 5)
    plt.plot([-x_axis_lim, x_axis_lim], [0, 0], Color = 'black', linewidth = 2)
    ax.scatter(x2, y2, s=5, Color = 'red', marker = 'o')
    plt.xlim([-x_axis_lim, x_axis_lim])
    plt.ylim([-y_axis_lim, y_axis_lim])
    plt.rcParams["figure.figsize"] = (10, 15)



if __name__ == "__main__":
    x0 = 0 # Initial x
    dx0 = 0 # Initial horizontal velocity
    y0 = 10000 # Initial y
    dy0 = -343*3 # Initial vertical velocity
    theta0 = 0 # Initial angle
    dtheta0 = 0 # Initial angular velocity
    thrust0 = [4000000, 0] # magnitude, direction [angle phi]
    length = 47.7
    mass = 25600+395700*0.1 # estimate https://www.spaceflightinsider.com/hangar/falcon-9/
    dt = 0.05 # Length of time step
    total = 400 # Total number of time steps
    x,y,theta = solve_eom(x0,dx0,y0,dy0,theta0,dtheta0,thrust0,length,mass,dt,total)

    fig, ax = plt.subplots(figsize = (10,13.5))
#    line, = ax.plot(x, y)
    ani = animation.FuncAnimation(fig, animate, frames = total, interval=2, blit=False, repeat = False)
    plt.show()
    # # Making movie of solution
    # for i in range(0,len(x)):
    #     length = 5
    #     x1 = x[i]-length*np.sin(theta[i])
    #     y1 = y[i]+length*np.cos(theta[i])
    #     x2 = x[i]+length*np.sin(theta[i])
    #     y2 = y[i]-length*np.cos(theta[i])
    #     print("x1: {}, x2: {}, y1: {}, y2: {}".format(x1, x2, y1, y2))
    #     plt.plot([x1, x2], [y1, y2], linewidth = 40)
    #     plt.xlim([-8,8])
    #     plt.ylim([-8,8])


    print('X:',x)
    print('Y:',y)
    print('Theta:',theta)