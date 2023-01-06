# RJ Meade and Noah Thompson, 1/5/2023
# Cornell University and Louisiana Stupid University (LSU)
# Creates a movie showing the movement of a rodket with a constant thrust along the axis of the rod
# x and y are the horizontal and vertical position of the rod center of mass, theta is its angle to vertical
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    return -thrust/mass*np.sin(theta) # horizontal thrust component

def DyDt(ydot):
    return ydot

def DydotDt(mass, thrust, theta):
    return thrust/mass*np.cos(theta) - 9.81 # vertical thrust component and gravity

def DthetaDt(thetadot):
    return thetadot

def DthetadotDt(): # would be torque/moment of inertia but there is no torque right now
    return 0

# Input: Initial conditions (assumes initial t=0) and the constant force of thrust and mass and time step dt and number of time steps to run for (including first time step)
# Solves DxDt, DxdotDt, etc with given initial conditions using Runge-Kutta
# Returns arrays of x, y, and theta where each index is a time step of dt
def solve_eom(x0, dx0, y0, dy0, theta0, dtheta0, thrust, mass, dt, totaldt):
    x = np.zeros(totaldt)
    y = np.zeros(totaldt)
    theta = np.zeros(totaldt)
    dx = np.zeros(totaldt)
    dy = np.zeros(totaldt)
    dtheta = np.zeros(totaldt)
    x[0] = x0
    y[0] = y0
    theta[0] = theta0
    dx[0] = dx0
    dy[0] = dy0
    dtheta[0] = dtheta0
    for i in range(1,totaldt):
        # Slopes at beginning of interval
        k11 = DxDt(dx[i-1])
        k12 = DxdotDt(mass,thrust,theta[i-1])
        k13 = DyDt(dy[i-1])
        k14 = DydotDt(mass, thrust, theta[i-1])
        k15 = DthetaDt(dtheta[i-1])
        k16 = DthetadotDt()
        # Estimate values at midpoint using initial slopes
        xmid = x[i-1] + k11*dt/2
        dxmid = dx[i-1] + k12*dt/2
        ymid = y[i-1] + k13*dt/2
        dymid = dy[i-1] + k14*dt/2
        thetamid = theta[i-1] + k15*dt/2
        dthetamid = dtheta[i-1] + k16*dt/2
        # Find slopes at estimated midpoint values
        k21 = DxDt(dxmid)
        k22 = DxdotDt(mass,thrust,thetamid)
        k23 = DyDt(dymid)
        k24 = DydotDt(mass, thrust, thetamid)
        k25 = DthetaDt(dthetamid)
        k26 = DthetadotDt()
        # Use midpoint slope predictions to get better estimate of midpoint values from start of interval
        xmid = x[i-1] + k21*dt/2
        dxmid = dx[i-1] + k22*dt/2
        ymid = y[i-1] + k23*dt/2
        dymid = dy[i-1] + k24*dt/2
        thetamid = theta[i-1] + k25*dt/2
        dthetamid = dtheta[i-1] + k26*dt/2
        # Use new midpoint values to get new estimates for midpoint slopes
        k31 = DxDt(dxmid)
        k32 = DxdotDt(mass,thrust,thetamid)
        k33 = DyDt(dymid)
        k34 = DydotDt(mass, thrust, thetamid)
        k35 = DthetaDt(dthetamid)
        k36 = DthetadotDt()
        # Use new midpoint slopes as estimate for entire interval to get estimated values at end of interval
        xend = x[i-1] + k31*dt
        dxend = dx[i-1] + k32*dt
        yend = y[i-1] + k33*dt
        dyend = dy[i-1] + k34*dt
        thetaend = theta[i-1] + k35*dt
        dthetaend = dtheta[i-1] + k36*dt
        # Use endpoint value estimates to get endpoint slope estimates
        k41 = DxDt(dxend)
        k42 = DxdotDt(mass,thrust, thetaend)
        k43 = DyDt(dyend)
        k44 = DydotDt(mass, thrust, thetaend)
        k45 = DthetaDt(dthetaend)
        k46 = DthetadotDt()
        # Best estimate of slope across entire interval is a weighted average of all four slope estimates, midpoint slopes weighted higher
        x[i] = x[i-1] + 1/6*(k11+2*k21+2*k31+k41)*dt
        dx[i] = dx[i-1] + 1/6*(k12+2*k22+2*k32+k42)*dt
        y[i] = y[i-1] + 1/6*(k13+2*k23+2*k33+k43)*dt
        dy[i] = dy[i-1] + 1/6*(k14+2*k24+2*k34+k44)*dt
        theta[i] = theta[i-1] + 1/6*(k15+2*k25+2*k35+k45)*dt
        dtheta[i] = dtheta[i-1] + 1/6*(k16+2*k26+2*k36+k46)*dt
    return x,y,theta

if __name__ == "__main__":
    x0 = 0 # Initial x
    dx0 = 0 # Initial horizontal velocity
    y0 = 0 # Initial y
    dy0 = 0 # Initial vertical velocity
    theta0 = np.sqrt(2)/2 # Initial angle
    dtheta0 = 0 # Initial angular velocity
    thrust = 9.81
    mass = 1
    dt = 1 # Length of time step
    total = 11 # Total number of time steps
    x,y,theta = solve_eom(x0,dx0,y0,dy0,theta0,dtheta0,thrust,mass,dt,total)

    # Making movie of solution

    print('X:',x)
    print('Y:',y)
    print('Theta:',theta)