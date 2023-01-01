import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def dVydt(m,thrust,theta):
    return -m*9.81+thrust*np.cos(theta)

def dVxdt(m,thrust,theta):
    return thrust*np.sin(theta)

def dThetadt(m,thrust,theta):
    return 0

def ode45(t_step, X0, m, thrust):
    k11 = dVydt(m,thrust,X0(1))
    k12 = dXydt()
    k21 = dVydt(m,thrust,X0(1)+k1*1/2*t_step)
    k31 = dVydt(m,thrust,X0(1)+k2*1/2*t_step)

def solve_eom():
    dVydt(thrust())

if __name__ == "__main__":
    print("hello")
    
    solve_eom()