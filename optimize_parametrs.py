
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Fixed parameters
x10, x20 = 0, 0  # initial positions (mm)
v10, v20 = 0, 0  # initial velocities (mm/s)
F = 27  # external force (N)
L = 546  # bar length (mm)
D_o = 50  # outer diameter of the bar (mm)
density1, density2 = 0.00000785, 0.0000193  # densities (kg/mm^3)
damping_ratio_2 = 0.05  # damping ratio of mass 2
E = 220000  # Young's modulus (N/mm^2)
cutting_speed = 942.5  # cutting speed (m/min)
fead_rate = 0.5  # feed rate (mm/rev)
depth_of_cut = 0.05  # depth of cut (mm)
omega1 = 17 - 0.566*cutting_speed + 3971*fead_rate + 155*depth_of_cut  # frequency (Hz)

# Layer 1: Inner diameter optimization
D_i_range = np.arange(5, 40, 5)  # inner diameter range (mm)
best_D_i, min_max_disp = 0, float('inf')

results = []

for D_i in D_i_range:
    I = np.pi * (D_o**4 - D_i**4) / 64  # moment of inertia (mm^4)
    m1 = (density1 * np.pi * L * ((D_o / 2)**2 - (D_i / 2)**2))  # mass 1 (kg)
    if m1 <= 0:  # Skip invalid configurations
        continue
    k1 = (3 * E * I) / (L**3)  # spring constant (N/mm)
    Y = (F * L**3) / (3 * E * I)  # deflection (mm)
    damping_ratio_1 = Y / np.sqrt((4 * np.pi**2) + Y**2)  # damping ratio of mass 1
    c1 = 2 * damping_ratio_1 * np.sqrt(k1 * m1)  # damping coefficient (Ns/mm)

    # Layer 2: Mass 2 parameters
    d_range = np.arange(30, 40, 1)  # diameter of mass 2 (mm)
    l_range = np.arange(50, 120, 10)  # length of mass 2 (mm)
    for d in d_range:
        for l in l_range:
            m2 = (density2 * np.pi * l * (d / 2)**2) # mass 2 (kg)
            if m2 <= 0:  # Skip invalid configurations
                continue
            k2_range = np.arange(115, 365, 25)  # spring constant range (N/mm)
            for k2 in k2_range:
                c2 = 2 * damping_ratio_2 * np.sqrt(k2 * m2)  # damping coefficient (Ns/mm)

                # Set time step and simulation range
                simTime=10 # simulation time (s)
                tStep=0.001 # simulation time step
                iterations=int(simTime/tStep) # total number of iterations
                t = np.linspace(0, simTime, iterations) # time array

                F_0 = F * np.sin(2*np.pi* omega1 * t)  # external force

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
                a1[0,:]=(F_0[0] + k2*x20 + c2*v20 -(k1+k2)*x10 -(c1+c2)*v10) /m1
                a2=np.zeros((iterations,1))
                a2[0,:]=(-c2*v20 - k2*x20 + c2*v10 + k2*x10) /m2

                # Solve the ODE's with Euler's Method
                for n in range(1,iterations):
                    x1[n,:]=x1[n-1,:]+v1[n-1,:]*tStep # mass 1 position
                    x2[n,:]=x2[n-1,:]+v2[n-1,:]*tStep # mass 2 position
                    v1[n,:]=v1[n-1,:]+a1[n-1,:]*tStep # mass 1 velocity
                    v2[n,:]=v2[n-1,:]+a2[n-1,:]*tStep # mass 2 velocity
                    # Find mass accelerations
                    a1[n,:]=(F_0[n-1] + k2*x2[n-1,:] + c2*v2[n-1,:] - (k1+k2)*x1[n-1,:] - (c1+c2)*v1[n-1,:])/m1
                    a2[n,:]=(-c2*v2[n-1,:] - k2*x2[n-1,:] + c2*v1[n-1,:] + k2*x1[n-1,:])/m2

                # Evaluate optimization: minimize max displacement
                t_start_idx = np.searchsorted(t, 6)  # מציאת האינדקס של t=6 שניות
                max_disp = np.max(np.abs(x1[t_start_idx:]))  # after initial transient 

                results.append({
                    "m1 (kg)": m1,
                    "m2 (kg)": m2,
                    "k1 (N/mm)": k1,
                    "k2 (N/mm)": k2,
                    "c1 (Ns/mm)": c1,
                    "c2 (Ns/mm)": c2,
                    "max_disp (mm)": max_disp
                }) 

                #max_disp = np.max(np.abs(x1))  # after initial transient
                if max_disp < min_max_disp:
                    min_max_disp = max_disp
                    best_D_i, best_d, best_l, best_k2, best_k1, best_m1, best_m2, best_c1, best_c2, best_I, best_ratio = D_i, d, l, k2, k1, m1, m2, c1, c2, I, damping_ratio_1

# Output optimized parameters
optimized_params = {
    "Best D_i (mm)": best_D_i,
    "Best d (mm)": best_d,
    "Best l (mm)": best_l,
    "Best k2 (N/mm)": best_k2,
    "Min Max Displacement (mm)": min_max_disp,
    'm1 (kg)': best_m1,
    'm2 (kg)': best_m2,
    'k1 (N/m)' : best_k1,
    'c1 (Ns/mm)': best_c1,
    'c2 (Ns/mm)': best_c2,
    'damper_ratio_1': best_ratio,
    'I': best_I
}
print("Optimized Parameters:", optimized_params)

# Save to CSV
pd.DataFrame([optimized_params]).to_csv("optimized_params.csv", index=False)

df_results = pd.DataFrame(results)
df_results.to_csv("all_max_disp_results.csv", index=False)













