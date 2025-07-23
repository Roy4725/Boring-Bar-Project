
import numpy as np
import matplotlib.pyplot as plt

# Problem parameters
k1=1281 # mass 1 spring constant (N/mm)
k2=115 # mass 2 spring constant (N/mm)
b1=10.65 # mass 1 viscous damping coefficient (Ns/mm)
b2=1.62 # mass 2 viscous damping coefficient (Ns/mm)
#m1=6.94 # mass 1 mass (kg)
m1 = 7.41
m2 = 2.35
#m2=2.28 # mass 2 mass (kg)
x10=0 # mass 1 initial position (mm)
x20=0 # mass 2 initial position (mm)
v10=0 # mass 1 initial velocity (mm/s)
v20=0 # mass 2 initial velocity (mm/s)
F=27 # external force (N)
omega=1470 # external force frequency (Hz)

# Set time step stuff
simTime=10 # simulation time (s)
tStep=0.001 # simulation time step
iterations=int(simTime/tStep) # total number of iterations
t = np.linspace(0, simTime, iterations) # time array

F_0 = F*np.sin(2*np.pi*omega*t) # external force function

# Pre-allocate variables for speed and add initial conditions
x1=np.zeros((iterations,1))
x1[0,:]=x10
x2=np.zeros((iterations,1))
x2[0,:]=x20
v1=np.zeros((iterations,1))
v1[0,:]=v10
v2=np.zeros((iterations,1))
v2[0,:]=v20
a1=np.zeros((iterations,1))
a1[0,:]=(F_0[0] + k2*x20 + b2*v20 -(k1+k2)*x10 -(b1+b2)*v10) /m1
a2=np.zeros((iterations,1))
a2[0,:]=(-b2*v20 - k2*x20 + b2*v10 + k2*x10) /m2

# Solve the ODE's with Euler's Method
for n in range(1,iterations):
  x1[n,:]=x1[n-1,:]+v1[n-1,:]*tStep # mass 1 position
  x2[n,:]=x2[n-1,:]+v2[n-1,:]*tStep # mass 2 position
  v1[n,:]=v1[n-1,:]+a1[n-1,:]*tStep # mass 1 velocity
  v2[n,:]=v2[n-1,:]+a2[n-1,:]*tStep # mass 2 velocity
  # Find mass accelerations
  a1[n,:]=(F_0[n-1] + k2*x2[n-1,:] + b2*v2[n-1,:] - (k1+k2)*x1[n-1,:] - (b1+b2)*v1[n-1,:])/m1
  a2[n,:]=(-b2*v2[n-1,:] - k2*x2[n-1,:] + b2*v1[n-1,:] + k2*x1[n-1,:])/m2

x30=0 # mass 1 initial position (mm)
v30=0 # mass 1 initial velocity (mm/s)
Om2 = np.sqrt((k1*1000)/m1) / (2*np.pi) # external force frequency (Hz)
F_1 = F*np.sin(Om2*t) # external force function

# Pre-allocate variables for speed and add initial conditions
x1_no_low=np.zeros((iterations,1))
x1_no_low[0,:]=x30
v1_no_low=np.zeros((iterations,1))
v1_no_low[0,:]=v30
a1_no_low=np.zeros((iterations,1))
a1_no_low[0,:]=(F_1[0] -k1*x30 -b1*v30) /m1

# Solve the ODE's with Euler's Method
for n in range(1,iterations):
  x1_no_low[n,:]=x1_no_low[n-1,:]+v1_no_low[n-1,:]*tStep # mass 1 position
  v1_no_low[n,:]=v1_no_low[n-1,:]+a1_no_low[n-1,:]*tStep # mass 1 velocity
  # Find mass accelerations
  a1_no_low[n,:]=(F_1[n-1] - k1*x1_no_low[n-1,:] - b1*v1_no_low[n-1,:])/m1

# Plot results
plt.rcParams["figure.figsize"] = (15,10) # resizes figures for viewing
plt.figure()
plt.plot(t,x1,'b',label='Bar')
#plt.plot(t,x2,'m',label='Capsul Mass')
#plt.plot(t,x1_no_low, 'g', label='Bar without Capsule')
plt.xlabel('Time (s)')
plt.ylabel('Position (mm)')
plt.title('Displacement')
#plt.xlim(3, 6)
plt.xlim(0, 10)
#plt.ylim(-0.000003, 0.000003)
plt.legend()
plt.show()






